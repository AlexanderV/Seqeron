using System.Text;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    private const string BiopythonUnknownDescription = "<unknown description>";

    /// <summary>The lines of a text as a Python text-mode stream gives them (universal newlines, the newline kept).</summary>
    private sealed class LineStream
    {
        private readonly List<string> _lines;
        private int _pos;

        public LineStream(string text) => _lines = [.. Lines(text.Replace("\r\n", "\n").Replace('\r', '\n'))];

        /// <summary>Python's <c>stream.readline()</c>: the next line or the empty string.</summary>
        public string ReadLine() => _pos < _lines.Count ? _lines[_pos++] : "";

        /// <summary>Python's <c>for line in stream</c> (one step); null at the end.</summary>
        public string? Next() => _pos < _lines.Count ? _lines[_pos++] : null;
    }

    private static string RStripNl(string s) => s.TrimEnd('\r', '\n');

    /// <summary>Python's <c>str.split(None, 1)</c> with a mandatory second part.</summary>
    private static (string, string) SplitWs1(string s)
    {
        string t = s.TrimStart();
        int i = 0;
        while (i < t.Length && !char.IsWhiteSpace(t[i]))
            i++;
        if (i == 0 || i >= t.Length)
            throw new ArgumentException("not enough values to unpack (expected 2, got 1)");
        return (t[..i], t[i..].TrimStart());
    }

    /// <summary><c>Alignment.length</c>: the number of printed columns of an alignment given by its coordinates.</summary>
    private static long AlignmentLength(long[][] coordinates)
    {
        int n = coordinates.Length;
        if (n == 0)
            return 0;
        int m = coordinates[0].Length - 1;
        var steps = new long[n][];
        for (int i = 0; i < n; i++)
        {
            steps[i] = new long[Math.Max(m, 0)];
            for (int c = 0; c < m; c++)
                steps[i][c] = coordinates[i][c + 1] - coordinates[i][c];
        }

        var aligned = new bool[Math.Max(m, 0)];
        for (int c = 0; c < m; c++)
            aligned[c] = Enumerable.Range(0, n).Count(i => steps[i][c] != 0) > 1;
        for (int i = 0; i < n; i++)
        {
            bool allNonNeg = true, allNonPos = true;
            for (int c = 0; c < m; c++)
            {
                if (!aligned[c])
                    continue;
                if (steps[i][c] < 0)
                    allNonNeg = false;
                if (steps[i][c] > 0)
                    allNonPos = false;
            }

            if (allNonNeg)
                continue;
            if (allNonPos)
            {
                for (int c = 0; c < m; c++)
                    steps[i][c] = -steps[i][c];
            }
            else
            {
                throw new ArgumentException($"Inconsistent steps in row {i}");
            }
        }

        long total = 0;
        for (int c = 0; c < m; c++)
        {
            long gap = long.MinValue;
            for (int i = 0; i < n; i++)
                gap = Math.Max(gap, steps[i][c]);
            for (int i = 0; i < n; i++)
            {
                if (!(steps[i][c] == gap || steps[i][c] <= 0))
                    throw new ArgumentException("Unequal step sizes in alignment");
            }

            total += gap;
        }

        return total;
    }

    /// <summary><c>Alignment(records, coordinates)</c> with <c>coordinates=None</c>: the sequences must have the same length (a record without a sequence has length 0).</summary>
    private static BioCoordAlignment AlignmentWithoutCoordinates(IReadOnlyList<BioCoordSequence> records, JsonObject attrs)
    {
        var lengths = records.Select(r => Math.Max(r.Length, 0)).Distinct().ToList();
        if (lengths.Count > 1)
            throw new ArgumentException("sequences must have the same length if coordinates is None");
        long[][] co = lengths.Count == 0 ? [] : [.. records.Select(_ => new long[] { 0, lengths[0] })];
        return new BioCoordAlignment { Sequences = records, Coordinates = co, Attributes = attrs };
    }

    private static bool PyIsDigits(string s) => s.Length > 0 && s.All(char.IsDigit);

    // ------------------------------------------------------------------ MSF

    /// <summary><c>Bio/Align/msf.py</c> <c>AlignmentIterator</c> (one alignment per file).</summary>
    private static CoordinateAlignmentFile ParseMsf(string text)
    {
        var stream = new LineStream(text.Replace("\r\n", "\n").Replace('\r', '\n'));
        var warnings = new List<string>();
        string line = stream.ReadLine();
        if (line.Length == 0)
            throw new ArgumentException("Empty file.");
        string[] knownHeaders = ["!!NA_MULTIPLE_ALIGNMENT", "!!AA_MULTIPLE_ALIGNMENT", "PileUp"];
        var first = Words(line);
        if (first.Length == 0)
            throw new ArgumentException("IndexError: list index out of range");
        if (!knownHeaders.Contains(first[0]))
            throw new ArgumentException($"{first[0]} is not a known GCG MSF header: {string.Join(", ", knownHeaders)}");
        bool found = false;
        string? l;
        while ((l = stream.Next()) is not null)
        {
            line = l.TrimEnd('\n');
            if (line.Contains("MSF: ", StringComparison.Ordinal) && line.EndsWith("..", StringComparison.Ordinal))
            {
                found = true;
                break;
            }
        }

        if (!found)
            throw new ArgumentException("Reached end of file without MSF/Type/Check header line");
        var parts = Words(line);
        int offset = Array.IndexOf(parts, "MSF:");
        if (offset < 0)
            throw new ArgumentException("'MSF:' is not in list");
        if (offset + 3 >= parts.Length || parts.Length < 3)
            throw new ArgumentException("IndexError");
        if (parts[offset + 2] != "Type:" || parts[^3] is not ("Check:" or "CompCheck:"))
            throw new ArgumentException("GCG MSF header line should be '<optional text> MSF: <int> Type: <letter> <optional date> Check: <int> ..'");
        long alnLength = PyInt(parts[offset + 1]);
        string seqType = parts[offset + 3];
        if (seqType is not ("P" or "N"))
            throw new ArgumentException($"GCG MSF header line should have 'Type: P' (protein) or 'Type: N' (nucleotide), not 'Type: {seqType}'");
        var names = new List<string>();
        var remaining = new List<long>();
        var weights = new List<double>();
        bool closed = false;
        while ((l = stream.Next()) is not null)
        {
            line = l.Trim();
            if (line == "//")
            {
                closed = true;
                break;
            }

            if (line.StartsWith("Name: ", StringComparison.Ordinal))
            {
                var words = Words(line);
                int iName = Array.IndexOf(words, "Name:"), iLen = Array.IndexOf(words, "Len:"), iWeight = Array.IndexOf(words, "Weight:"), iCheck = Array.IndexOf(words, "Check:");
                if (iName < 0 || iLen < 0 || iWeight < 0 || iCheck < 0)
                    throw new ArgumentException($"Malformed GCG MSF name line: '{line}'");
                string name = words[iName + 1];
                long length = PyInt(words[iLen + 1]);
                double weight = PyFloatParse(words[iWeight + 1]);
                _ = words[iCheck + 1];
                if (names.Contains(name))
                    throw new ArgumentException($"Duplicated ID of '{name}'");
                names.Add(name);
                remaining.Add(length);
                weights.Add(weight);
            }
        }

        if (!closed)
            throw new ArgumentException("End of file while looking for end of header // line.");
        line = stream.ReadLine();
        if (line.Length == 0)
            throw new ArgumentException("End of file after // line, expected sequences.");
        if (line.Trim().Length > 0)
            throw new ArgumentException("After // line, expected blank line before sequences.");
        var seqs = names.Select(_ => new StringBuilder()).ToList();
        bool done = false;
        while ((l = stream.Next()) is not null)
        {
            var words = Words(l);
            if (words.Length == 0)
                continue;
            string name = words[0];
            int index = names.IndexOf(name);
            if (index < 0)
            {
                if (words.All(PyIsDigits))
                    continue;
                throw new ArgumentException($"Unexpected line '{l}' in input");
            }

            string seq = string.Concat(words.Skip(1));
            long length = remaining[index] - (seq.Length - seq.Count(c => c == '-'));
            if (length < 0)
                throw new ArgumentException($"Received longer sequence than expected for {name}");
            seqs[index].Append(seq);
            remaining[index] = length;
            if (remaining.All(x => x == 0))
            {
                done = true;
                break;
            }
        }

        if (!done)
            throw new ArgumentException("End of file where expecting sequence data.");
        while ((l = stream.Next()) is not null)
        {
            if (l.Trim().Length > 0)
                throw new ArgumentException("AssertionError");
        }

        if (seqs.Count == 0)
            throw new ArgumentException("max() arg is an empty sequence");
        long maxLen = seqs.Max(s => s.Length);
        if (maxLen != alnLength)
        {
            warnings.Add($"GCG MSF headers said alignment length {alnLength}, but found {maxLen}");
            alnLength = maxLen;
        }

        var rows = new List<string>();
        foreach (var sb in seqs)
        {
            string seq = sb.ToString().Replace('~', '-').Replace('.', '-');
            if (seq.Length < alnLength)
                seq += new string('-', (int)(alnLength - seq.Length));
            rows.Add(seq);
        }

        var (ungapped, co) = ParsePrintedAlignment(rows);
        var records = new List<BioCoordSequence>();
        for (int i = 0; i < names.Count; i++)
        {
            records.Add(new BioCoordSequence
            {
                Id = names[i],
                Length = ungapped[i].Length,
                Parts = SeqParts([new KeyValuePair<long, string>(0, ungapped[i])], ungapped[i].Length),
                Description = names[i],
                Annotations = new JsonObject { ["weight"] = weights[i] },
            });
        }

        long columns = AlignmentLength(co);
        if (columns != alnLength)
            throw new ArgumentException($"GCG MSF headers said alignment length {alnLength}, but found {columns}");
        return new CoordinateAlignmentFile
        {
            Alignments = [new BioCoordAlignment { Sequences = records, Coordinates = co }],
            Metadata = [],
            Warnings = warnings,
        };
    }

    // ------------------------------------------------------------------ EMBOSS

    private const string EmbossHash = "########################################";
    private const string EmbossEq = "#=======================================";

    /// <summary><c>Bio/Align/emboss.py</c> <c>AlignmentIterator</c>.</summary>
    private static CoordinateAlignmentFile ParseEmboss(string text)
    {
        var stream = new LineStream(text);
        string line = stream.ReadLine();
        if (line.Length == 0)
            throw new ArgumentException("Empty file.");
        if (line.TrimEnd() != EmbossHash)
            throw new ArgumentException("TypeError: unsupported operand type(s) for %: 'ValueError' and 'str'");
        var metadata = new JsonObject { ["Align_format"] = "srspair" };
        StringBuilder? commandline = null;
        string? l;
        while ((l = stream.Next()) is not null)
        {
            line = l;
            if (line.TrimEnd() == EmbossHash)
                break;
            if (!line.StartsWith("# ", StringComparison.Ordinal))
                throw new ArgumentException("TypeError: unsupported operand type(s) for %: 'ValueError' and 'str'");
            if (commandline is not null)
            {
                if (line.StartsWith("#    ", StringComparison.Ordinal))
                {
                    commandline.Append(' ').Append(line[1..].Trim());
                    continue;
                }

                metadata["Command line"] = commandline.ToString();
                commandline = null;
            }

            int colon = line.IndexOf(':', 2);
            if (colon < 0)
                throw new ArgumentException("not enough values to unpack (expected 2, got 1)");
            string key = line[2..colon], value = line[(colon + 1)..];
            switch (key)
            {
                case "Program": metadata["Program"] = value.Trim(); break;
                case "Rundate": metadata["Rundate"] = value.Trim(); break;
                case "Report_file": metadata["Report_file"] = value.Trim(); break;
                case "Align_format": metadata["Align_format"] = value.Trim(); break;
                case "Commandline": commandline = new StringBuilder(value.Trim()); break;
            }
        }

        string alignFormat = metadata["Align_format"]!.GetValue<string>();
        var alignments = new List<BioCoordAlignment>();
        while (true)
        {
            var a = ReadEmbossAlignment(stream, alignFormat);
            if (a is null)
                break;
            alignments.Add(a);
        }

        return new CoordinateAlignmentFile { Alignments = alignments, Metadata = metadata };
    }

    private static BioCoordAlignment? ReadEmbossAlignment(LineStream stream, string alignFormat)
    {
        var annotations = new JsonObject();
        long? numberOfSequences = null;
        string? l;
        var identifiers = new List<string>();
        long? ncols = null;
        while ((l = stream.Next()) is not null)
        {
            string line = RStripNl(l);
            if (line.Length == 0)
                continue;
            if (line.StartsWith("#---------------------------------------", StringComparison.Ordinal))
                continue;
            if (line.StartsWith(EmbossEq, StringComparison.Ordinal))
                break;
            throw new ArgumentException($"Unexpected line: {line}");
        }

        bool closed = false;
        while ((l = stream.Next()) is not null)
        {
            string line = RStripNl(l);
            if (line == EmbossEq)
            {
                closed = true;
                break;
            }

            if (line.Trim() == "#")
                continue;
            if (!line.StartsWith("# ", StringComparison.Ordinal))
                throw new ArgumentException("TypeError: unsupported operand type(s) for %: 'ValueError' and 'str'");
            string key, value;
            int colon = line.IndexOf(':', 2);
            if (colon >= 0)
            {
                key = line[2..colon];
                value = line[(colon + 1)..];
            }
            else
            {
                int eq = line.IndexOf(" = ", 2, StringComparison.Ordinal);
                if (eq < 0)
                    throw new ArgumentException("not enough values to unpack (expected 2, got 1)");
                key = line[2..eq];
                value = line[(eq + 3)..];
            }

            switch (key)
            {
                case "Aligned_sequences":
                    numberOfSequences = PyInt(value.Trim());
                    if (identifiers.Count != 0)
                        throw new ArgumentException("AssertionError");
                    int i = 0;
                    string? il;
                    while ((il = stream.Next()) is not null)
                    {
                        if (!il.StartsWith("# ", StringComparison.Ordinal))
                            throw new ArgumentException("TypeError: unsupported operand type(s) for %: 'ValueError' and 'str'");
                        var pr = il[2..].Split(':');
                        if (pr.Length != 2)
                            throw new ArgumentException("ValueError: unpacking");
                        if (i + 1 != PyInt(pr[0]))
                            throw new ArgumentException("AssertionError");
                        identifiers.Add(pr[1].Trim());
                        i++;
                        if (identifiers.Count == numberOfSequences)
                            break;
                    }

                    break;
                case "Matrix": annotations[key] = value.Trim(); break;
                case "Gap_penalty":
                case "Extend_penalty":
                case "Score":
                    annotations[key] = PyFloatParse(value.Trim());
                    break;
                case "Length": ncols = PyInt(value.Trim()); break;
                case "Identity":
                case "Similarity":
                case "Gaps":
                    annotations[key] = PyInt(value.Trim().Split('/')[0]);
                    break;
                case "Longest_Identity":
                case "Longest_Similarity":
                case "Shortest_Identity":
                case "Shortest_Similarity":
                    annotations[key] = value.Trim();
                    break;
                default:
                    throw new ArgumentException($"Failed to parse line '{line}'");
            }
        }

        if (!closed)
            return null;
        if (identifiers.Count == 0)
            throw new ArgumentException("Number of sequences missing!");
        if (ncols is null)
            throw new ArgumentException("Length of alignment missing!");
        int n = (int)(numberOfSequences ?? throw new ArgumentException("TypeError"));
        var sequences = Enumerable.Range(0, n).Select(_ => new StringBuilder()).ToArray();
        var aligned = Enumerable.Range(0, n).Select(_ => new StringBuilder()).ToArray();
        var consensus = new StringBuilder();
        var starts = new long[n];
        var ends = new long[n];
        long column = 0;
        int index = 0;
        while ((l = stream.Next()) is not null)
        {
            string line = RStripNl(l);
            if (line.Length == 0)
            {
                if (index == n)
                {
                    if (column == ncols)
                        break;
                    index = 0;
                }

                continue;
            }

            string prefix = (line.Length > 21 ? line[..21] : line).Trim();
            if (prefix.Length == 0)
            {
                consensus.Append(line.Length > 21 ? line[21..Math.Min(71, line.Length)] : "");
                continue;
            }

            var (identifier, startText) = SplitWs1(prefix);
            if (index >= identifiers.Count)
                throw new ArgumentException("IndexError");
            if (!identifiers[index].StartsWith(identifier, StringComparison.Ordinal))
                throw new ArgumentException("AssertionError");
            var (alignedSequence, endText) = SplitWs1(line.Length > 21 ? line[21..] : "");
            long start = PyInt(startText), end = PyInt(endText);
            long length = sequences[index].Length;
            string sequence = alignedSequence.Replace("-", "");
            if (length == 0 && sequence.Length > 0)
            {
                if (start < end)
                {
                    start--;
                    if (end != start + sequence.Length)
                        throw new ArgumentException("AssertionError");
                }
                else
                {
                    end--;
                    if (end != start - sequence.Length)
                        throw new ArgumentException("AssertionError");
                }

                starts[index] = start;
            }
            else if (starts[index] <= ends[index])
            {
                if (alignFormat == "srspair" && sequence.Length == 0)
                {
                    if (start != ends[index] || end != start)
                        throw new ArgumentException("AssertionError");
                }
                else
                {
                    start--;
                    if (end != start + sequence.Length)
                        throw new ArgumentException("AssertionError");
                }
            }
            else if (alignFormat == "srspair" && sequence.Length == 0)
            {
                if (start - 1 != ends[index] || end != start)
                    throw new ArgumentException("AssertionError");
            }
            else
            {
                end--;
                if (end != start - sequence.Length)
                    throw new ArgumentException("AssertionError");
            }

            ends[index] = end;
            sequences[index].Append(sequence);
            aligned[index].Append(alignedSequence);
            if (index == 0)
                column += alignedSequence.Length;
            else if (column != aligned[index].Length)
                throw new ArgumentException("AssertionError");
            index++;
        }

        var (data, co) = ParsePrintedAlignment([.. aligned.Select(a => a.ToString())]);
        var records = new List<BioCoordSequence>();
        for (int i = 0; i < data.Count; i++)
        {
            long start = starts[i], end = ends[i];
            string d = data[i];
            BioCoordSequence rec;
            if (start == 0)
            {
                rec = new BioCoordSequence { Id = null, Length = d.Length, Parts = SeqParts([new KeyValuePair<long, string>(0, d)], d.Length) };
            }
            else if (start < end)
            {
                for (int c = 0; c < co[i].Length; c++)
                    co[i][c] += start;
                rec = new BioCoordSequence { Id = null, Length = end, Parts = SeqParts([new KeyValuePair<long, string>(start, d)], end) };
            }
            else
            {
                (start, end) = (end, start);
                for (int c = 0; c < co[i].Length; c++)
                    co[i][c] = end - co[i][c];
                rec = ReverseComplement(new BioCoordSequence { Id = null, Length = end, Parts = SeqParts([new KeyValuePair<long, string>(0, d)], end) });
            }

            if (i >= identifiers.Count)
                throw new ArgumentException("IndexError");
            records.Add(new BioCoordSequence { Id = identifiers[i], Length = rec.Length, Parts = rec.Parts, Description = BiopythonUnknownDescription });
        }

        var attrs = new JsonObject();
        if (annotations.Count > 0)
            attrs["annotations"] = annotations;
        if (consensus.Length > 0)
            attrs["column_annotations"] = new JsonObject { ["emboss_consensus"] = consensus.ToString() };
        return new BioCoordAlignment { Sequences = records, Coordinates = co, Attributes = attrs };
    }
}
