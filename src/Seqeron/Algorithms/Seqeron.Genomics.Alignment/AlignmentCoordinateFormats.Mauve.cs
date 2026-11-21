using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    /// <summary><c>Bio/Align/mauve.py</c> <c>AlignmentWriter.write_header</c>.</summary>
    private static string MauveHeader(CoordinateWriteOptions options)
    {
        var metadata = options.Metadata ?? throw new ArgumentException("The Mauve writer needs the metadata of the file (alignments.metadata).");
        var identifiers = options.Identifiers ?? throw new ArgumentException("The Mauve writer needs the identifiers of the sequences (alignments.identifiers).");
        var sb = new StringBuilder();
        sb.Append("#FormatVersion ").Append(metadata.TryGetPropertyValue("FormatVersion", out var fv) && fv is not null ? MetaString(fv) : "Mauve1").Append('\n');
        string? filename = metadata.TryGetPropertyValue("File", out var fn) && fn is not null ? MetaString(fn) : null;
        if (filename is null)
        {
            for (int index = 0; index < identifiers.Count; index++)
            {
                int number = index + 1;
                sb.Append("#Sequence").Append(Str(number)).Append("File\t").Append(identifiers[index]).Append('\n');
                sb.Append("#Sequence").Append(Str(number)).Append("Format\tFastA\n");
            }
        }
        else
        {
            for (int index = 0; index < identifiers.Count; index++)
            {
                if (index != PyInt(identifiers[index]))
                    throw new ArgumentException("AssertionError");
                int number = index + 1;
                sb.Append("#Sequence").Append(Str(number)).Append("File\t").Append(filename).Append('\n');
                sb.Append("#Sequence").Append(Str(number)).Append("Entry\t").Append(Str(number)).Append('\n');
                sb.Append("#Sequence").Append(Str(number)).Append("Format\tFastA\n");
            }
        }

        if (metadata.TryGetPropertyValue("BackboneFile", out var bb) && bb is not null)
            sb.Append("#BackboneFile\t").Append(MetaString(bb)).Append('\n');
        return sb.ToString();
    }

    /// <summary><c>Bio/Align/mauve.py</c> <c>AlignmentWriter.format_alignment</c>.</summary>
    private static string FormatMauve(BioCoordAlignment alignment, CoordinateWriteOptions options)
    {
        var metadata = options.Metadata ?? throw new ArgumentException("The Mauve writer needs the metadata of the file.");
        int n = alignment.Coordinates.Length;
        int m = n == 0 ? 0 : alignment.Coordinates[0].Length;
        if (n == 0)
            throw new ArgumentException("Must have at least one sequence");
        if (m == 0)
            throw new ArgumentException("Non-empty sequences are required");
        string? filename = metadata.TryGetPropertyValue("File", out var fn) && fn is not null ? MetaString(fn) : null;
        var sb = new StringBuilder();
        for (int i = 0; i < n; i++)
        {
            string identifier = alignment.Sequences[i].Id ?? "None";
            long start = alignment.Coordinates[i][0], end = alignment.Coordinates[i][^1];
            string strand = "+";
            if (start > end)
            {
                strand = "-";
                (start, end) = (end, start);
            }

            if (start == end)
            {
                if (start != 0)
                    throw new ArgumentException("AssertionError");
            }
            else
            {
                start++;
            }

            string sequence = GetRow(alignment, i);
            if (filename is null)
            {
                int index = options.Identifiers?.ToList().IndexOf(identifier) ?? -1;
                if (index < 0)
                    throw new ArgumentException($"'{identifier}' is not in list");
                sb.Append("> ").Append(Str(index + 1)).Append(':').Append(Str(start)).Append('-').Append(Str(end)).Append(' ').Append(strand).Append(' ').Append(identifier).Append('\n');
            }
            else
            {
                sb.Append("> ").Append(Str(PyInt(identifier) + 1)).Append(':').Append(Str(start)).Append('-').Append(Str(end)).Append(' ').Append(strand).Append(' ').Append(filename).Append('\n');
            }

            sb.Append(sequence).Append('\n');
        }

        sb.Append("=\n");
        return sb.ToString();
    }

    /// <summary><c>Bio/Align/mauve.py</c> <c>AlignmentIterator</c>.</summary>
    private static List<BioCoordAlignment> ParseMauve(string text)
    {
        var all = Lines(text).ToList();
        int p = 0;
        string[] suffixes = ["File", "Entry", "Format"];
        var idInfo = suffixes.ToDictionary(x => x, _ => new List<string>());
        bool metadataAny = false;
        string? pending = null;
        while (p < all.Count)
        {
            string line = all[p++];
            if (!line.StartsWith('#'))
            {
                pending = line.Trim();
                break;
            }

            var kv = Words(line[1..]);
            if (kv.Length != 2)
                throw new ArgumentException("not enough values to unpack");
            string key = kv[0], value = kv[1];
            metadataAny = true;
            if (key.StartsWith("Sequence", StringComparison.Ordinal))
            {
                string? suffix = suffixes.FirstOrDefault(x => key.EndsWith(x, StringComparison.Ordinal)) ?? throw new ArgumentException($"Unexpected keyword '{key}'");
                if (suffix == "Entry")
                    value = Str(PyInt(value) - 1);
                long seqNum = PyInt(key["Sequence".Length..^suffix.Length]);
                idInfo[suffix].Add(value);
                if (seqNum != idInfo[suffix].Count)
                    throw new ArgumentException("AssertionError");
            }
        }

        if (pending is null && !metadataAny)
            throw new ArgumentException("Empty file.");
        var files = idInfo["File"];
        List<string> identifiers;
        if (files.Distinct().Count() == 1)
        {
            identifiers = idInfo["Entry"];
        }
        else
        {
            if (files.Distinct().Count() != files.Count)
                throw new ArgumentException("AssertionError");
            identifiers = files;
        }

        var result = new List<BioCoordAlignment>();
        var descriptions = new List<(string Id, long Start, long End, string Strand, string Comments)>();
        var seqs = new List<StringBuilder>();
        (string, long, long, string, string) Describe(string line)
        {
            var w = line[1..].Split((char[]?)null, 3, StringSplitOptions.RemoveEmptyEntries);
            if (w.Length != 3)
                throw new ArgumentException("not enough values to unpack");
            var locus = w[0].Split(':');
            if (locus.Length != 2)
                throw new ArgumentException("not enough values to unpack");
            long seqNum = PyInt(locus[0]) - 1;
            if (seqNum < 0 || seqNum >= identifiers.Count)
                throw new ArgumentException("list index out of range");
            if (w[1] is not ("+" or "-"))
                throw new ArgumentException("AssertionError");
            var se = locus[1].Split('-');
            if (se.Length != 2)
                throw new ArgumentException("not enough values to unpack");
            long start = PyInt(se[0]), end = PyInt(se[1]);
            if (start == 0)
            {
                if (end != 0)
                    throw new ArgumentException("AssertionError");
            }
            else
            {
                start--;
            }

            return (identifiers[(int)seqNum], start, end, w[1], w[2]);
        }

        if (pending is not null)
        {
            var d = Describe(pending);
            descriptions.Add(d);
            seqs.Add(new StringBuilder());
        }

        while (p < all.Count)
        {
            string line = all[p++].Trim();
            if (line.StartsWith('='))
            {
                var (rows, co) = ParsePrintedAlignment(seqs.Select(x => x.ToString()).ToList());
                var records = new List<BioCoordSequence>();
                for (int index = 0; index < descriptions.Count; index++)
                {
                    var (identifier, start, end, strand, comments) = descriptions[index];
                    string seq = rows[index];
                    if (seq.Length != end - start)
                        throw new ArgumentException("AssertionError");
                    if (strand == "-")
                    {
                        for (int c = 0; c < co[index].Length; c++)
                            co[index][c] = seq.Length - co[index][c];
                    }

                    for (int c = 0; c < co[index].Length; c++)
                        co[index][c] += start;
                    BioCoordSequence rec;
                    if (start == 0)
                    {
                        rec = new BioCoordSequence { Id = identifier, Length = seq.Length, Parts = new SortedDictionary<long, string> { [0] = seq }, Description = comments };
                        if (strand == "-")
                            rec = ReverseComplement(rec);
                    }
                    else if (strand == "+")
                    {
                        rec = new BioCoordSequence { Id = identifier, Length = end, Parts = new SortedDictionary<long, string> { [start] = seq }, Description = comments };
                    }
                    else
                    {
                        rec = ReverseComplement(new BioCoordSequence { Id = identifier, Length = end, Parts = new SortedDictionary<long, string> { [0] = seq }, Description = comments });
                    }

                    records.Add(new BioCoordSequence { Id = rec.Id, Length = rec.Length, Parts = rec.Parts, Description = comments });
                }

                result.Add(new BioCoordAlignment { Sequences = records, Coordinates = co, Attributes = [] });
                descriptions = [];
                seqs = [];
            }
            else if (line.StartsWith('>'))
            {
                descriptions.Add(Describe(line));
                seqs.Add(new StringBuilder());
            }
            else
            {
                if (seqs.Count == 0)
                    throw new ArgumentException("list index out of range");
                seqs[^1].Append(line);
            }
        }

        return result;
    }
}
