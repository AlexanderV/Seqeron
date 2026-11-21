using System.IO.Compression;
using System.Text;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    /// <summary>The autoSQL declaration of <c>Bio/Align/bigmaf.py</c> (<c>bedMaf</c>: bed3 plus the <c>lstring</c> MAF block).</summary>
    public static BigBedDeclaration BigMafDeclaration { get; } = new(
        "bedMaf",
        "Bed3 with MAF block",
        [
            new BigBedField("string", "chrom", "Reference sequence chromosome or scaffold"),
            new BigBedField("uint", "chromStart", "Start position in chromosome"),
            new BigBedField("uint", "chromEnd", "End position in chromosome"),
            new BigBedField("lstring", "mafBlock", "MAF block"),
        ]);

    /// <summary>
    /// Writes the alignments as <c>Bio.Align.write(alignments, handle, "bigmaf")</c> does: every alignment becomes a bed3 item of its reference row (<c>reference.chromosome</c> -&gt; <c>chromosome</c>) whose <c>mafBlock</c> is the MAF text of the whole
    /// alignment with <c>;</c> for the line ends. The block size and items per slot of the wrapper are not passed on by Biopython (the inner writer uses its defaults), so they are not options here either.
    /// </summary>
    /// <param name="alignments">The alignments (rows = species; the first row is the reference <c>reference.chromosome</c>).</param>
    /// <param name="targets">The chromosomes (<c>alignments.targets</c>); the first one is named <c>reference.chromosome</c>, the others by their chromosome name.</param>
    /// <param name="options">Only <see cref="BigBedWriteOptions.Compress"/> and <see cref="BigBedWriteOptions.Compressor"/> are used.</param>
    public static byte[] WriteBigMaf(IReadOnlyList<BioCoordAlignment> alignments, IReadOnlyList<BioCoordSequence> targets, BigBedWriteOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(alignments);
        ArgumentNullException.ThrowIfNull(targets);
        try
        {
            return WriteBigMafCore(alignments, targets, options ?? new BigBedWriteOptions());
        }
        catch (Exception e) when (e is IndexOutOfRangeException or ArgumentOutOfRangeException or KeyNotFoundException or InvalidCastException or NullReferenceException or FormatException or OverflowException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, e);
        }
    }

    private static string ChromosomeOf(string id)
    {
        int dot = id.IndexOf('.', StringComparison.Ordinal);
        if (dot < 0)
            throw new ArgumentException("ValueError: not enough values to unpack (expected 2, got 1)");
        return id[(dot + 1)..];
    }

    private static byte[] WriteBigMafCore(IReadOnlyList<BioCoordAlignment> alignments, IReadOnlyList<BioCoordSequence> targets, BigBedWriteOptions options)
    {
        var fixedAlignments = new List<BioCoordAlignment>();
        foreach (var alignment in alignments)
        {
            string mafBlock = FormatMaf(alignment)[..^1].Replace('\n', ';');
            if (alignment.Coordinates.Length == 0 || alignment.Coordinates[0].Length == 0)
                continue;
            string chromosome = ChromosomeOf(alignment.Sequences[0].Id ?? throw new ArgumentException("AttributeError: 'NoneType' object has no attribute 'split'"));
            var co = alignment.Coordinates;
            if (co[0][0] >= co[0][^1])
                throw new ArgumentException("AssertionError");
            var target = alignment.Sequences[0];
            var renamed = new BioCoordSequence { Id = chromosome, Length = target.Length, Parts = target.Parts, Description = target.Description, Annotations = target.Annotations };
            var sliced = new BioCoordAlignment
            {
                Sequences = [.. alignment.Sequences.Take(2).Select((s, i) => i == 0 ? renamed : s)],
                Coordinates = [.. co.Take(2).Select(r => (long[])r.Clone())],
                Attributes = new JsonObject { ["annotations"] = new JsonObject { ["mafBlock"] = mafBlock } },
            };
            fixedAlignments.Add(sliced);
        }

        var sorted = fixedAlignments.OrderBy(a => a.Sequences[0].Id, StringComparer.Ordinal).ThenBy(a => a.Coordinates[0][0]).ToList();
        var record = targets[0];
        string firstChromosome = ChromosomeOf(record.Id ?? throw new ArgumentException("AttributeError: 'NoneType' object has no attribute 'split'"));
        var fixedTargets = targets.ToList();
        fixedTargets[0] = new BioCoordSequence { Id = firstChromosome, Length = record.Length, Parts = record.Parts };
        var inner = new BigBedWriteOptions
        {
            BedN = 3,
            Compress = options.Compress,
            Compressor = options.Compressor,
            DeclarationFields = BigMafDeclaration.Fields,
            DeclarationName = BigMafDeclaration.Name,
            DeclarationComment = BigMafDeclaration.Comment,
        };
        return WriteBigBedCore(sorted, fixedTargets, inner);
    }

    /// <summary>Reads a bigMaf file as <c>Bio.Align.parse(handle, "bigmaf")</c>: the alignments (rows = the sequences of the MAF blocks), the chromosomes (the first one named <c>reference.chromosome</c>) and the declaration.</summary>
    public static CoordinateAlignmentFile ParseBigMaf(byte[] data)
    {
        ArgumentNullException.ThrowIfNull(data);
        try
        {
            return ParseBigBedCore(data, true);
        }
        catch (Exception e) when (e is IndexOutOfRangeException or ArgumentOutOfRangeException or KeyNotFoundException or InvalidCastException or NullReferenceException or FormatException or OverflowException or InvalidOperationException or InvalidDataException)
        {
            throw new ArgumentException(e.Message, e);
        }
    }

    /// <summary>The reference species of a bigMaf file: the part of the first <c>s</c> source name before the first dot (<c>AlignmentIterator._read_reference</c>).</summary>
    private static string BigMafReference(byte[] firstBlock)
    {
        int i = IndexOf(firstBlock, ";s"u8, 12);
        if (i < 0)
            throw new ArgumentException("ValueError: subsection not found");
        var words = PySplit(Encoding.Latin1.GetString(firstBlock, i + 1, firstBlock.Length - i - 1), -1);
        if (words.Count <= 2)
            throw new ArgumentException("ValueError: no sequence name found");
        string name = words[1];
        int dot = name.IndexOf('.', StringComparison.Ordinal);
        if (dot < 0)
            throw new ArgumentException("ValueError: not enough values to unpack (expected 2, got 1)");
        return Latin1ToUtf8(name[..dot]);
    }

    private static int IndexOf(byte[] data, ReadOnlySpan<byte> pattern, int from)
    {
        if (from > data.Length)
            return -1;
        int i = data.AsSpan(from).IndexOf(pattern);
        return i < 0 ? -1 : i + from;
    }

    private static string Latin1ToUtf8(string s) => Encoding.UTF8.GetString(Encoding.Latin1.GetBytes(s));

    private static bool IsPyAsciiSpace(char c) => c is ' ' or '\t' or '\n' or '\r' or '\v' or '\f';

    /// <summary>Python's <c>bytes.split(None, maxsplit)</c>: leading white space dropped, the remainder after <paramref name="maxSplit"/> splits keeps its trailing white space; negative = unlimited.</summary>
    private static List<string> PySplit(string s, int maxSplit)
    {
        var result = new List<string>();
        int pos = 0;
        while (true)
        {
            while (pos < s.Length && IsPyAsciiSpace(s[pos]))
                pos++;
            if (pos >= s.Length)
                break;
            if (maxSplit >= 0 && result.Count == maxSplit)
            {
                result.Add(s[pos..]);
                break;
            }

            int start = pos;
            while (pos < s.Length && !IsPyAsciiSpace(s[pos]))
                pos++;
            result.Add(s[start..pos]);
        }

        return result;
    }

    private static readonly Regex MafSourceLine = new(@"^s\s*\S*\s*\d*\s*\d*\s*[+-]\s*\d*\s*", RegexOptions.ECMAScript | RegexOptions.CultureInvariant);

    /// <summary><c>bigmaf.AlignmentIterator._create_alignment</c>: one MAF block (lines separated by <c>;</c>, ended by a NUL).</summary>
    private static BioCoordAlignment BbCreateMafAlignment(byte[] restBytes)
    {
        string buf = Encoding.Latin1.GetString(restBytes) + "\0";
        var records = new List<BioCoordSequence>();
        var strands = new List<string>();
        var rows = new List<string>();
        var annotations = new JsonObject();
        var empties = new JsonArray();
        double? score = null;
        string src = "";
        JsonObject? lastRecordAnnotations = null;
        int j = -1;
        while (true)
        {
            int i = j + 1;
            char prefix = i < buf.Length ? buf[i] : '￿';
            if (prefix == '\0')
                break;
            if (prefix == '#')
            {
                j = SemicolonEnd(buf, i);
            }
            else if (prefix == 'a')
            {
                j = SemicolonEnd(buf, i);
                foreach (string word in PySplit(buf[(i + 1)..j], -1))
                {
                    var kv = word.Split('=');
                    if (kv.Length != 2)
                        throw new ArgumentException("ValueError: not enough values to unpack");
                    if (kv[0] == "score")
                    {
                        score = ParseDouble(kv[1]);
                    }
                    else if (kv[0] == "pass")
                    {
                        long v = PyInt(kv[1]);
                        if (v <= 0)
                            throw new ArgumentException($"pass value must be positive (found {v})");
                        annotations["pass"] = v;
                    }
                    else
                    {
                        throw new ArgumentException($"Unknown annotation variable '{Latin1ToUtf8(kv[0])}'");
                    }
                }
            }
            else if (prefix == 's')
            {
                var m = MafSourceLine.Match(buf[i..]);
                if (!m.Success)
                    throw new ArgumentException("AttributeError: 'NoneType' object has no attribute 'span'");
                j = i + m.Length;
                var words = PySplit(buf[i..j], 5);
                if (words.Count != 6)
                    throw new ArgumentException("Error parsing alignment - 's' line must have 7 fields");
                src = Latin1ToUtf8(words[1]);
                long start = PyInt(words[2]), size = PyInt(words[3]), srcSize = PyInt(words[5]);
                string strand = words[4];
                i = j;
                int end = i;
                while (end < buf.Length && buf[end] is not (';' or '\0'))
                    end++;
                string text = buf[i..end];
                string sequence = text.Replace("-", "", StringComparison.Ordinal);
                if (sequence.Length != size)
                    throw new ArgumentException($"sequence size is incorrect (found {sequence.Length}, expected {size})");
                if (start < 0)
                    throw new ArgumentException("KeyError: 0");
                var parts = SeqParts([new KeyValuePair<long, string>(start, sequence)], srcSize);
                lastRecordAnnotations = [];
                records.Add(new BioCoordSequence { Id = src, Length = srcSize, Parts = parts, Annotations = lastRecordAnnotations });
                rows.Add(text);
                strands.Add(strand);
                j = end;
            }
            else if (prefix == 'i')
            {
                j = SemicolonEnd(buf, i);
                var words = PySplit(buf[i..j], 5);
                if (words.Count != 6 || Latin1ToUtf8(words[1]) != src || lastRecordAnnotations is null)
                    throw new ArgumentException("AssertionError");
                string[] statusChars = ["C", "I", "N", "n", "M", "T"];
                if (!statusChars.Contains(words[2]) || !statusChars.Contains(words[4]))
                    throw new ArgumentException("AssertionError");
                lastRecordAnnotations["leftStatus"] = words[2];
                lastRecordAnnotations["leftCount"] = PyInt(words[3]);
                lastRecordAnnotations["rightStatus"] = words[4];
                lastRecordAnnotations["rightCount"] = PyInt(words[5]);
            }
            else if (prefix == 'e')
            {
                j = SemicolonEnd(buf, i);
                var words = PySplit(buf[i..j], 6);
                if (words.Count != 7)
                    throw new ArgumentException("AssertionError");
                src = Latin1ToUtf8(words[1]);
                long start = PyInt(words[2]), size = PyInt(words[3]), srcSize = PyInt(words[5]);
                if (words[6] is not ("C" or "I" or "M" or "n"))
                    throw new ArgumentException("AssertionError");
                long end = start + size;
                long[] segment = words[4] == "+" ? [start, end] : [srcSize - start, srcSize - end];
                empties.Add(new JsonObject
                {
                    ["record"] = new JsonObject { ["id"] = src, ["length"] = srcSize, ["parts"] = null },
                    ["segment"] = new JsonArray(segment[0], segment[1]),
                    ["status"] = words[6],
                });
            }
            else if (prefix == 'q')
            {
                j = SemicolonEnd(buf, i);
                var words = PySplit(buf[i..j], 2);
                if (words.Count != 3 || Latin1ToUtf8(words[1]) != src || lastRecordAnnotations is null)
                    throw new ArgumentException("AssertionError");
                lastRecordAnnotations["quality"] = Latin1ToUtf8(words[2].Replace("-", "", StringComparison.Ordinal));
            }
            else
            {
                throw new ArgumentException("Error parsing alignment - unexpected line");
            }
        }

        if (empties.Count > 0)
            annotations["empty"] = empties;
        var coordinates = PrintedAlignmentShape(rows);
        for (int r = 0; r < records.Count; r++)
        {
            if (strands[r] == "+")
            {
                long first = records[r].Parts!.First().Key;
                for (int c = 0; c < coordinates[r].Length; c++)
                    coordinates[r][c] += first;
            }
            else
            {
                records[r] = ReverseComplement(records[r]);
                var p = records[r].Parts!.First();
                long last = p.Key + p.Value.Length;
                for (int c = 0; c < coordinates[r].Length; c++)
                    coordinates[r][c] = last - coordinates[r][c];
            }
        }

        var attrs = new JsonObject { ["annotations"] = annotations };
        if (score is not null)
            attrs["score"] = score.Value;
        return new BioCoordAlignment { Sequences = records, Coordinates = coordinates, Attributes = attrs };
    }

    private static int SemicolonEnd(string buf, int from)
    {
        int k = buf.IndexOf(';', from);
        return k < 0 ? buf.Length : k;
    }

    /// <summary>The coordinates <c>PrintedAlignmentParser.fill</c> produces for printed rows of equal width; unequal rows leave the C buffer uninitialised in Biopython, so they are refused.</summary>
    private static long[][] PrintedAlignmentShape(List<string> rows)
    {
        if (rows.Count == 0)
            return [];
        if (rows.Any(r => r.Length != rows[0].Length))
            throw new ArgumentException("the printed rows have different lengths (undefined in Biopython)");
        return ParsePrintedAlignment(rows).Coordinates;
    }
}
