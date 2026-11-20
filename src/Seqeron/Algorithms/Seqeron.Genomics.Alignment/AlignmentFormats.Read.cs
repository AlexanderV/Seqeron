using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// The readers of Bio.Align (fasta.py, a2m.py, clustal.py, phylip.py, nexus.py + the Bio.Nexus matrix parser, stockholm.py): what Bio.Align.parse returns, as BioAlignment.
public static partial class AlignmentFormats
{
    /// <summary>Reads the alignments of a file as <c>Bio.Align.parse(handle, format)</c> does (FASTA, A2M, PHYLIP and NEXUS files hold one alignment; the Clustal reader reads to the end of the file as one).</summary>
    public static IReadOnlyList<BioAlignment> Parse(string text, AlignmentFileFormat format)
    {
        ArgumentNullException.ThrowIfNull(text);
        var lines = SplitLines(text);
        if (lines.Count == 0)
            throw new FormatException("Empty file.");
        return format switch
        {
            AlignmentFileFormat.Fasta => ParseFasta(lines),
            AlignmentFileFormat.A2m => ParseA2m(lines),
            AlignmentFileFormat.Clustal => ParseClustal(lines),
            AlignmentFileFormat.Phylip => ParsePhylip(lines),
            AlignmentFileFormat.Nexus => ParseNexus(lines),
            _ => ParseStockholm(lines),
        };
    }

    // the lines as a Python file iterator yields them: universal newlines, each keeping its trailing "\n" (the last one only if the file ends with it)
    private static List<string> SplitLines(string text)
    {
        text = text.Replace("\r\n", "\n", StringComparison.Ordinal).Replace('\r', '\n');
        var lines = new List<string>();
        int pos = 0;
        while (pos < text.Length)
        {
            int nl = text.IndexOf('\n', pos);
            if (nl < 0)
            {
                lines.Add(text[pos..]);
                break;
            }

            lines.Add(text.Substring(pos, nl - pos + 1));
            pos = nl + 1;
        }

        return lines;
    }

    // str.split(None, 1) of a stripped string
    private static string[] SplitWhitespace(string s, int maxSplit)
    {
        var parts = new List<string>();
        int i = 0;
        while (i < s.Length && char.IsWhiteSpace(s[i]))
            i++;
        while (i < s.Length)
        {
            if (parts.Count == maxSplit)
            {
                parts.Add(s[i..]);
                return [.. parts];
            }

            int j = i;
            while (j < s.Length && !char.IsWhiteSpace(s[j]))
                j++;
            parts.Add(s[i..j]);
            i = j;
            while (i < s.Length && char.IsWhiteSpace(s[i]))
                i++;
        }

        return [.. parts];
    }

    // Alignment.parse_printed_alignment: the rows must have one length (Biopython's parser reads uninitialised memory otherwise); columns that are gaps in every row do not exist in the alignment
    private static List<string> PrintedRows(IReadOnlyList<string> rows)
    {
        if (rows.Count == 0)
            return [];
        int len = rows[0].Length;
        if (rows.Any(r => r.Length != len))
            throw new FormatException("Sequences of different length in the alignment.");
        var keep = new List<int>();
        for (int c = 0; c < len; c++)
        {
            if (rows.Any(r => r[c] != '-'))
                keep.Add(c);
        }

        return [.. rows.Select(r => new string([.. keep.Select(c => r[c])]))];
    }

    private static BioAlignment MakeAlignment(IReadOnlyList<string> ids, IReadOnlyList<string> descriptions, IReadOnlyList<string> rows, JsonObject? column = null, JsonObject? annotations = null, Func<int, JsonObject>? recordAnnotations = null)
    {
        var records = new List<BioAlignmentRecord>();
        for (int i = 0; i < ids.Count; i++)
            records.Add(new BioAlignmentRecord { Id = ids[i], Description = descriptions[i], Row = rows[i], Annotations = recordAnnotations?.Invoke(i) ?? [] });
        return new BioAlignment { Records = records, ColumnAnnotations = column ?? [], Annotations = annotations ?? [] };
    }

    // ----------------------------------------------------------------- fasta ----

    private static (List<string> Names, List<string> Descriptions, List<string> Lines) ReadFastaRecords(IReadOnlyList<string> fileLines)
    {
        var names = new List<string>();
        var descriptions = new List<string>();
        var rows = new List<string>();
        foreach (string line in fileLines)
        {
            if (line.StartsWith('>'))
            {
                var parts = SplitWhitespace(line[1..].TrimEnd(), 1);
                names.Add(parts.Length > 0 ? parts[0] : "");
                descriptions.Add(parts.Length == 2 ? parts[1] : "");
                rows.Add("");
            }
            else
            {
                if (rows.Count == 0)
                    throw new FormatException("Sequence data before the first '>' line.");
                rows[^1] += line.Trim();
            }
        }

        return (names, descriptions, rows);
    }

    private static List<BioAlignment> ParseFasta(IReadOnlyList<string> fileLines)
    {
        var (names, descriptions, rows) = ReadFastaRecords(fileLines);
        if (rows.Count == 0)
            return [];
        return [MakeAlignment(names, descriptions, PrintedRows(rows))];
    }

    // ------------------------------------------------------------------- a2m ----

    private static List<BioAlignment> ParseA2m(IReadOnlyList<string> fileLines)
    {
        var (names, descriptions, lines) = ReadFastaRecords(fileLines);
        if (lines.Count == 0)
            return [];
        var state = new StringBuilder();
        foreach (char c in lines[0])
        {
            if (c == '-' || char.IsUpper(c))
                state.Append('D');
            else if (c == '.' || char.IsLower(c))
                state.Append('I');
            else
                throw new FormatException($"Unexpected letter '{c}' in alignment");
        }

        for (int i = 1; i < lines.Count; i++)
        {
            for (int k = 0; k < Math.Min(lines[i].Length, state.Length); k++)
            {
                char c = lines[i][k];
                if (state[k] == 'D' && !(c == '-' || char.IsUpper(c)))
                    throw new FormatException("Inconsistent match state in an A2M row.");
                if (state[k] == 'I' && !(c == '.' || char.IsLower(c)))
                    throw new FormatException("Inconsistent insertion state in an A2M row.");
            }
        }

        var rows = PrintedRows([.. lines.Select(l => l.ToUpperInvariant().Replace('.', '-'))]);
        return [MakeAlignment(names, descriptions, rows, new JsonObject { ["state"] = state.ToString() })];
    }

    // --------------------------------------------------------------- clustal ----

    private static readonly string[] ClustalPrograms = ["CLUSTAL", "PROBCONS", "MUSCLE", "MSAPROBS", "Kalign", "Biopython"];

    private static List<BioAlignment> ParseClustal(IReadOnlyList<string> fileLines)
    {
        var words = SplitWhitespace(fileLines[0], int.MaxValue);
        if (words.Length == 0 || !ClustalPrograms.Contains(words[0]))
            throw new FormatException($"{(words.Length > 0 ? words[0] : "")} is not known to generate CLUSTAL files: {string.Join(", ", ClustalPrograms)}");
        var ids = new List<string>();
        var aligned = new List<string>();
        string consensus = "";
        int index = -1;
        string lastAligned = "";
        int li = 1;
        bool broke = false;
        while (li < fileLines.Count)
        {
            string line = fileLines[li++];
            if (line.StartsWith(' '))
            {
                if (ids.Count == 0 || index < 0)
                    throw new FormatException("Consensus line before the sequences.");
                int length = lastAligned.Length;
                consensus = SafeSubstring(line, index, index + length);
                broke = true;
                break;
            }

            if (line.Trim().Length > 0)
            {
                var fields = SplitWhitespace(line, int.MaxValue);
                if (fields.Length < 2 || fields.Length > 3)
                    throw new FormatException($"Could not parse line:\n{line}");
                ids.Add(fields[0]);
                aligned.Add(fields[1]);
                lastAligned = fields[1];
                if (index < 0)
                    index = line.IndexOf(fields[1], fields[0].Length, StringComparison.Ordinal);
                CheckClustalCount(fields, fields[1], line);
            }
            else if (index > 0)
            {
                broke = true;
                break;
            }
        }

        if (!broke)
            return [];
        int length0 = aligned[0].Length;
        if (aligned.Any(s => s.Length != length0))
            throw new FormatException("Sequences of different length in the first block.");
        int n = aligned.Count;
        int i2 = 0;
        var consensusBuilder = new StringBuilder(consensus);
        for (; li < fileLines.Count; li++)
        {
            string line = fileLines[li];
            if (line.StartsWith(' '))
            {
                consensusBuilder.Append(SafeSubstring(line, index, index + lastAligned.Length));
            }
            else if (line.Trim().Length == 0)
            {
                continue;
            }
            else
            {
                var fields = SplitWhitespace(line, int.MaxValue);
                if (fields.Length < 2 || fields.Length > 3)
                    throw new FormatException($"Could not parse line:\n{line}");
                if (fields[0] != ids[i2])
                    throw new FormatException("Sequence identifiers differ between blocks.");
                lastAligned = fields[1];
                aligned[i2] += fields[1];
                CheckClustalCount(fields, aligned[i2], line);
                i2++;
                if (i2 == n)
                    i2 = 0;
            }
        }

        consensus = consensusBuilder.ToString();
        var rows = PrintedRows(aligned);
        JsonObject? column = null;
        if (consensus.Length > 0)
        {
            int columns = rows[0].Length;
            if (consensus.Length != columns)
                throw new FormatException($"Alignment has {columns} columns, consensus length is {consensus.Length}, '{consensus}'");
            column = new JsonObject { ["clustal_consensus"] = consensus };
        }

        return [MakeAlignment(ids, ids.Select(_ => "").ToList(), rows, column)];
    }

    private static void CheckClustalCount(string[] fields, string aligned, string line)
    {
        if (fields.Length != 3)
            return;
        if (!int.TryParse(fields[2], NumberStyles.Integer, CultureInfo.InvariantCulture, out int count))
            throw new FormatException($"Could not parse line, bad sequence count:\n{line}");
        if (aligned.Count(c => c != '-') != count)
            throw new FormatException($"Could not parse line, incorrect sequence count:\n{line}");
    }

    // ---------------------------------------------------------------- phylip ----

    private static List<BioAlignment> ParsePhylip(IReadOnlyList<string> fileLines)
    {
        var header = SplitWhitespace(fileLines[0], int.MaxValue);
        if (header.Length != 2 || !int.TryParse(header[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out int nSeqs) || !int.TryParse(header[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int seqLength))
            throw new FormatException($"Expected two integers in the first line, received '{fileLines[0]}'");
        int pos = 1;
        string ReadLine() => pos < fileLines.Count ? fileLines[pos++] : "";
        var names = new List<string>();
        var seqs = new List<List<string>>();
        var firstLines = new List<string>();
        for (int i = 0; i < nSeqs; i++)
            firstLines.Add(ReadLine());
        string next = ReadLine();
        if (next.TrimEnd().Length > 0)
        {
            firstLines.Add(next);
            int length = ParsePhylipSequential(firstLines, seqs, names, 0, seqLength);
            var rest = new List<string>();
            while (pos < fileLines.Count)
                rest.Add(fileLines[pos++]);
            ParsePhylipSequential(rest, seqs, names, length, seqLength);
        }
        else
        {
            foreach (string l in firstLines)
            {
                string line = l.TrimEnd();
                names.Add((line.Length > PhylipIdWidth ? line[..PhylipIdWidth] : line).Trim());
                seqs.Add([(line.Length > PhylipIdWidth ? line[PhylipIdWidth..] : "").Trim().Replace(" ", "", StringComparison.Ordinal)]);
            }

            int i = 0;
            while (pos < fileLines.Count)
            {
                string line = fileLines[pos++].TrimEnd();
                if (line.Length == 0)
                {
                    if (i != nSeqs)
                        throw new FormatException("Unexpected file format");
                    i = 0;
                }
                else
                {
                    seqs[i].Add(line.Replace(" ", "", StringComparison.Ordinal));
                    i++;
                }
            }

            if (i != 0 && i != nSeqs)
                throw new FormatException("Unexpected file format");
        }

        var joined = seqs.Select(s => string.Concat(s)).ToList();
        if (joined.Count != nSeqs)
            throw new FormatException($"Found {joined.Count} records in this alignment, told to expect {nSeqs}");
        foreach (string s in joined)
        {
            if (s.Length != seqLength)
                throw new FormatException($"Expected all sequences to have length {seqLength}; found {s.Length}");
            if (s.Contains('.', StringComparison.Ordinal))
                throw new FormatException("PHYLIP format no longer allows dots in sequence");
        }

        return [MakeAlignment(names, names.Select(_ => "").ToList(), PrintedRows(joined))];
    }

    private static int ParsePhylipSequential(IEnumerable<string> lines, List<List<string>> seqs, List<string> names, int length, int seqLength)
    {
        foreach (string l in lines)
        {
            string seq;
            if (length == 0)
            {
                string line = l.TrimEnd();
                names.Add((line.Length > PhylipIdWidth ? line[..PhylipIdWidth] : line).Trim());
                seq = (line.Length > PhylipIdWidth ? line[PhylipIdWidth..] : "").Trim();
                seqs.Add([]);
            }
            else
            {
                seq = l.Trim();
            }

            seq = seq.Replace(" ", "", StringComparison.Ordinal);
            seqs[^1].Add(seq);
            length += seq.Length;
            if (length == seqLength)
                length = 0;
        }

        return length;
    }

    // ----------------------------------------------------------------- nexus ----

    private const string NexusPunct = "()[]{}\\,;:=*\\'\"`+-<>";

    // Bio.Nexus.CharBuffer.next_word on one line: returns the word (quotes kept, doubled quotes collapsed) and the rest of the line
    private static (string? Word, string Remainder) NexusNextWord(string line)
    {
        int i = 0;
        while (i < line.Length && NexusWhitespace.Contains(line[i], StringComparison.Ordinal))
            i++;
        if (i >= line.Length)
            return (null, "");
        var word = new StringBuilder();
        char first = line[i++];
        word.Append(first);
        char quoted = '\0';
        if (first is '\'' or '"')
            quoted = first;
        else if (NexusPunct.Contains(first, StringComparison.Ordinal))
            return (first.ToString(), line[i..]);
        while (true)
        {
            char c = i < line.Length ? line[i] : '\0';
            if (c == quoted && quoted != '\0')
            {
                word.Append(c);
                i++;
                if (i < line.Length && line[i] == quoted)
                    i++;
                else
                    break;
            }
            else if (quoted != '\0')
            {
                if (i >= line.Length)
                    break;
                word.Append(line[i++]);
            }
            else if (c == '\0' || NexusPunct.Contains(c, StringComparison.Ordinal) || NexusWhitespace.Contains(c, StringComparison.Ordinal))
            {
                break;
            }
            else
            {
                word.Append(line[i++]);
            }
        }

        return (word.ToString(), line[i..]);
    }

    private static string? NexusQuoteStrip(string? word)
    {
        if (string.IsNullOrEmpty(word))
            return null;
        while ((word.StartsWith('\'') && word.EndsWith('\'') && word.Length >= 2) || (word.StartsWith('"') && word.EndsWith('"') && word.Length >= 2))
            word = word[1..^1];
        return word;
    }

    // the commands of a NEXUS file: text up to ';' outside quotes and [comments]
    private static List<string> NexusCommands(string text)
    {
        var commands = new List<string>();
        var cur = new StringBuilder();
        char quote = '\0';
        int comment = 0;
        foreach (char c in text)
        {
            if (comment > 0)
            {
                if (c == '[')
                    comment++;
                else if (c == ']')
                    comment--;
                continue;
            }

            if (quote != '\0')
            {
                cur.Append(c);
                if (c == quote)
                    quote = '\0';
                continue;
            }

            if (c == '[')
            {
                comment = 1;
            }
            else if (c is '\'' or '"')
            {
                quote = c;
                cur.Append(c);
            }
            else if (c == ';')
            {
                commands.Add(cur.ToString());
                cur.Clear();
            }
            else
            {
                cur.Append(c);
            }
        }

        if (cur.ToString().Trim().Length > 0)
            commands.Add(cur.ToString());
        return commands;
    }

    private static List<BioAlignment> ParseNexus(IReadOnlyList<string> fileLines)
    {
        if (fileLines[0].Trim() != "#NEXUS")
            throw new FormatException("File does not start with NEXUS header.");
        var commands = NexusCommands(string.Concat(fileLines.Skip(1)));
        bool inData = false;
        int ntax = 0;
        int nchar = 0;
        string datatype = "dna";
        bool interleave = false;
        string? matrixText = null;
        foreach (string raw in commands)
        {
            string cmd = raw.Trim();
            string lower = cmd.ToLowerInvariant();
            if (lower.StartsWith("begin", StringComparison.Ordinal))
            {
                inData = Regex.IsMatch(lower, @"^begin\s+(data|characters)\b");
            }
            else if (lower.StartsWith("end", StringComparison.Ordinal) || lower.StartsWith("endblock", StringComparison.Ordinal))
            {
                inData = false;
            }
            else if (inData && lower.StartsWith("dimensions", StringComparison.Ordinal))
            {
                var m1 = Regex.Match(cmd, @"ntax\s*=\s*(\d+)", RegexOptions.IgnoreCase);
                var m2 = Regex.Match(cmd, @"nchar\s*=\s*(\d+)", RegexOptions.IgnoreCase);
                ntax = m1.Success ? int.Parse(m1.Groups[1].Value, CultureInfo.InvariantCulture) : ntax;
                nchar = m2.Success ? int.Parse(m2.Groups[1].Value, CultureInfo.InvariantCulture) : nchar;
            }
            else if (inData && lower.StartsWith("format", StringComparison.Ordinal))
            {
                var m = Regex.Match(cmd, @"datatype\s*=\s*(\S+)", RegexOptions.IgnoreCase);
                if (m.Success)
                    datatype = m.Groups[1].Value.ToLowerInvariant();
                interleave = Regex.IsMatch(cmd, @"\binterleave\b", RegexOptions.IgnoreCase);
            }
            else if (inData && lower.StartsWith("matrix", StringComparison.Ordinal))
            {
                matrixText = cmd[6..];
            }
        }

        if (matrixText is null)
            return [];
        if (ntax == 0 || nchar == 0)
            throw new FormatException("Dimensions must be specified before matrix!");
        var lines = matrixText.Split('\n').Select(l => l.Trim()).Where(l => l.Length > 0).ToList();
        var labels = new List<string>();
        var matrix = new List<string>();
        int taxcount = 0;
        int li = 0;
        bool firstBlock = true;
        while (true)
        {
            if (li >= lines.Count)
            {
                if (taxcount < ntax)
                    throw new FormatException("Not enough taxa in matrix.");
                if (taxcount > ntax)
                    throw new FormatException("Too many taxa in matrix.");
                break;
            }

            string line = lines[li++];
            taxcount++;
            if (taxcount > ntax)
            {
                if (!interleave)
                    throw new FormatException("Too many taxa in matrix - should matrix be interleaved?");
                taxcount = 1;
                firstBlock = false;
            }

            var (word, remainder) = NexusNextWord(line);
            string id = NexusQuoteStrip(word) ?? "";
            line = remainder.Trim();
            var charsBuilder = new StringBuilder(string.Concat(SplitWhitespace(line, int.MaxValue)));
            if (interleave)
            {
                if (line.Length == 0)
                {
                    if (li >= lines.Count)
                        throw new FormatException("Not enough lines in matrix.");
                    charsBuilder.Append(string.Concat(SplitWhitespace(lines[li++], int.MaxValue)));
                }
            }
            else
            {
                while (charsBuilder.Length < nchar)
                {
                    if (li >= lines.Count)
                        throw new FormatException("Not enough taxa in matrix.");
                    charsBuilder.Append(string.Concat(SplitWhitespace(lines[li++], int.MaxValue)));
                }
            }

            string chars = charsBuilder.ToString();

            if (firstBlock)
            {
                labels.Add(id);
                matrix.Add(chars);
            }
            else
            {
                matrix[taxcount - 1] += chars;
            }
        }

        string? molecule = datatype switch { "dna" or "nucleotide" => "DNA", "rna" => "RNA", "protein" => "protein", _ => null };
        var rows = PrintedRows(matrix);
        // the NEXUS reader builds SeqRecord(seq, id=name, annotations=...) without a description: Biopython's default text
        return [MakeAlignment(labels, labels.Select(_ => "<unknown description>").ToList(), rows, recordAnnotations: _ => molecule is null ? [] : new JsonObject { ["molecule_type"] = molecule })];
    }
}
