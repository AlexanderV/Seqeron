using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The pairwise coordinate formats of Biopython 1.88 <c>Bio.Align</c> (<c>bed</c>, <c>psl</c>, <c>chain</c>): writers byte for byte as <c>Bio.Align.write</c> and readers that give what <c>Bio.Align.parse</c> gives
/// (sequences of the right length with unknown or block-wise known residues, coordinates, the attributes <c>score</c>, <c>matches</c>, ... of the alignment object). They are the genomic-coordinate counterpart of
/// <see cref="AlignmentFormats"/> (row-oriented text formats); the shared walk over the coordinate blocks is the one of Biopython's writers.
/// </summary>
public static partial class AlignmentCoordinateFormats
{
    /// <summary>Parses a Biopython format name (case-insensitive).</summary>
    public static bool TryParseFormat(string? name, out CoordinateAlignmentFormat format)
    {
        format = default;
        switch (name?.Trim().ToLowerInvariant())
        {
            case "bed": format = CoordinateAlignmentFormat.Bed; return true;
            case "psl": format = CoordinateAlignmentFormat.Psl; return true;
            case "chain": format = CoordinateAlignmentFormat.Chain; return true;
            case "maf": format = CoordinateAlignmentFormat.Maf; return true;
            default: return false;
        }
    }

    /// <summary>The lower-case Biopython name of the format.</summary>
    public static string FormatName(CoordinateAlignmentFormat format) => format.ToString().ToLowerInvariant();

    /// <summary>Writes the alignments as <c>Bio.Align.write(alignments, handle, format, ...)</c> does (the PSL header once, then one block / line per alignment).</summary>
    public static string Write(IReadOnlyList<BioCoordAlignment> alignments, CoordinateAlignmentFormat format, CoordinateWriteOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(alignments);
        options ??= new CoordinateWriteOptions();
        var sb = new StringBuilder();
        if (format == CoordinateAlignmentFormat.Psl && options.Header)
            sb.Append(PslHeader);
        if (format == CoordinateAlignmentFormat.Maf)
            sb.Append(MafHeader(options.Metadata));
        foreach (var a in alignments)
            sb.Append(Format(a, format, options));
        return sb.ToString();
    }

    /// <summary>Writes one alignment (<c>alignment.format(format)</c>: no file header).</summary>
    public static string Format(BioCoordAlignment alignment, CoordinateAlignmentFormat format, CoordinateWriteOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        options ??= new CoordinateWriteOptions();
        return format switch
        {
            CoordinateAlignmentFormat.Bed => FormatBed(alignment, options.BedN),
            CoordinateAlignmentFormat.Psl => FormatPsl(alignment, options),
            CoordinateAlignmentFormat.Maf => FormatMaf(alignment),
            _ => FormatChain(alignment),
        };
    }

    /// <summary>Parses the alignments of a text as <c>Bio.Align.parse(handle, format)</c> does.</summary>
    /// <exception cref="ArgumentException">Biopython raises <c>ValueError</c> for the text.</exception>
    public static IReadOnlyList<BioCoordAlignment> Parse(string text, CoordinateAlignmentFormat format)
    {
        ArgumentNullException.ThrowIfNull(text);
        return format switch
        {
            CoordinateAlignmentFormat.Bed => ParseBed(text),
            CoordinateAlignmentFormat.Psl => ParsePsl(text),
            CoordinateAlignmentFormat.Maf => ParseMaf(text),
            _ => ParseChain(text),
        };
    }

    // ------------------------------------------------------------------ shared helpers

    private static string Str(long v) => v.ToString(CultureInfo.InvariantCulture);

    private static long[][] Copy(long[][] c) => [.. c.Select(r => (long[])r.Clone())];

    private static void ReverseColumns(long[][] c)
    {
        foreach (var row in c)
            Array.Reverse(row);
    }

    private static bool TryAttr(BioCoordAlignment a, string key, out JsonNode node)
    {
        if (a.Attributes.TryGetPropertyValue(key, out var n) && n is not null)
        {
            node = n;
            return true;
        }

        node = null!;
        return false;
    }

    private static long Int(JsonNode n) => n.GetValueKind() == System.Text.Json.JsonValueKind.Number && n.AsValue().TryGetValue<long>(out var l) ? l : (long)n.AsValue().GetValue<double>();

    /// <summary>Python's <c>str(value)</c> of an attribute (numbers as Python prints them: integers plain, floats with <c>repr</c>).</summary>
    private static string PyStr(JsonNode n)
    {
        var v = n.AsValue();
        if (v.TryGetValue<string>(out var s))
            return s;
        if (v.TryGetValue<long>(out var l))
            return Str(l);
        double d = v.GetValue<double>();
        return PyFloat(d);
    }

    /// <summary>Python's <c>repr(float)</c>: the shortest round-trip digits, <c>.0</c> for integral values, exponents from 1e16.</summary>
    internal static string PyFloat(double d)
    {
        if (double.IsNaN(d))
            return "nan";
        if (double.IsInfinity(d))
            return d > 0 ? "inf" : "-inf";
        string r = d.ToString("R", CultureInfo.InvariantCulture);
        if (r.Contains('E'))
        {
            double abs = Math.Abs(d);
            if (abs >= 1e16 || abs < 1e-4)
            {
                int e = r.IndexOf('E');
                string mant = r[..e];
                string exp = r[(e + 1)..];
                int ex = int.Parse(exp, CultureInfo.InvariantCulture);
                return mant + "e" + (ex < 0 ? "-" : "+") + Math.Abs(ex).ToString("00", CultureInfo.InvariantCulture);
            }

            return d.ToString("0.0################", CultureInfo.InvariantCulture);
        }

        return r.Contains('.') ? r : r + ".0";
    }

    /// <summary>Python's <c>format(value, "g")</c>.</summary>
    private static string G(JsonNode n)
    {
        var v = n.AsValue();
        double d = v.TryGetValue<long>(out var l) ? l : v.GetValue<double>();
        return CFormat.General(d, 6);
    }

    private static string Id(BioCoordSequence s, string fallback) => s.Id ?? fallback;

    private static (BioCoordSequence Target, BioCoordSequence Query) Pair(BioCoordAlignment a)
    {
        if (a.Sequences.Count != 2)
            throw new ArgumentException("Expected an alignment of two sequences (target and query).", nameof(a));
        return (a.Sequences[0], a.Sequences[1]);
    }

    private static string ReverseComplement(string s)
    {
        bool rna = s.Contains('U') || s.Contains('u');
        var sb = new StringBuilder(s.Length);
        for (int i = s.Length - 1; i >= 0; i--)
        {
            char c = s[i];
            sb.Append(c switch
            {
                'A' => rna ? 'U' : 'T',
                'a' => rna ? 'u' : 't',
                'T' => 'A',
                't' => 'a',
                'U' => 'A',
                'u' => 'a',
                'C' => 'G',
                'c' => 'g',
                'G' => 'C',
                'g' => 'c',
                'R' => 'Y',
                'r' => 'y',
                'Y' => 'R',
                'y' => 'r',
                'K' => 'M',
                'k' => 'm',
                'M' => 'K',
                'm' => 'k',
                'B' => 'V',
                'b' => 'v',
                'V' => 'B',
                'v' => 'b',
                'D' => 'H',
                'd' => 'h',
                'H' => 'D',
                'h' => 'd',
                _ => c,
            });
        }

        return sb.ToString();
    }

    private static BioCoordSequence ReverseComplement(BioCoordSequence s)
    {
        if (s.Parts is null)
            return s;
        var parts = new SortedDictionary<long, string>();
        foreach (var (start, text) in s.Parts)
            parts[s.Length - (start + text.Length)] = ReverseComplement(text);
        return new BioCoordSequence { Id = s.Id, Length = s.Length, Parts = parts, Description = s.Description, Annotations = s.Annotations };
    }

    private static string[] Words(string line) => line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);

    private static IEnumerable<string> Lines(string text)
    {
        int pos = 0;
        while (pos < text.Length)
        {
            int nl = text.IndexOf('\n', pos);
            if (nl < 0)
            {
                yield return text[pos..];
                yield break;
            }

            yield return text[pos..(nl + 1)];
            pos = nl + 1;
        }
    }

    private static long PyInt(string s)
    {
        if (!long.TryParse(s.Trim(), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long v))
            throw new ArgumentException($"invalid literal for int() with base 10: '{s}'");
        return v;
    }

    private static List<long> IntList(string field)
    {
        var list = new List<long>();
        foreach (var p in field.TrimEnd(',').Split(','))
            list.Add(PyInt(p));
        return list;
    }

    private static BioCoordSequence Unknown(string? id, long length) => new() { Id = id, Length = length, Parts = null };
}
