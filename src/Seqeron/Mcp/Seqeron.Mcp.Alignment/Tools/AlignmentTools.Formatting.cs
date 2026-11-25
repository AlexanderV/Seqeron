using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-27b: the `format` parameter of the pairwise and multiple-alignment tools, mapped to the canonical alignment writers (AlignmentFormats = Biopython Bio.Align, EmbossAlignmentWriter = EMBOSS -aformat).
public partial class AlignmentTools
{
    public const string FormatParamDescription = "Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text.";

    private const string PairFormatParamDescription = FormatParamDescription + " This pairwise tool (it knows its sequences) also writes the alignment in the COORDINATE formats of Biopython 1.88 Alignment.format, byte for byte: bed (bed3 ... bed12), psl, chain, maf, sam (sam:md adds the MD tag), exonerate (exonerate:cigar), and the strings cigar (cigar:eqx = extended = / X; soft clips for the unaligned query), md (SAM MD tag) and cs (cs:long; minimap2 difference string). In those formats sequence1 is the target (reference) and sequence2 the query, except for the semi-global tools where sequence2 is the reference (target); the plain name sam is the Biopython SAM line (the EMBOSS one is 'emboss:sam'); a 'bio:' prefix forces the coordinate writer.";

    private const string TrimEndGapsParamDescription = "Coordinate formats (bed, psl, chain, maf, sam, exonerate, cigar, md, cs): drop the leading and trailing columns that have a gap in one row (the free end gaps of a semi-global result), so that the record starts and ends with an aligned pair (Biopython alignment[:, first:last]); false (default) keeps them as D / I blocks as Biopython's PairwiseAligner does (2D4M2D).";

    public const string NamesParamDescription = "Optional names of the rows for the formatted output (default seq1, seq2, ...).";

    // the Biopython Bio.Align writers (AlignmentFormats) for rows; a2m gets all-match states, nexus the molecule type guessed from the residues
    public static string WriteBiopythonFormat(AlignmentFileFormat bio, string[] rows, string[]? names, string[]? descriptions)
    {
        var records = new List<BioAlignmentRecord>();
        string all = string.Concat(rows).ToUpperInvariant();
        string molecule;
        if (all.All(c => "ACGTNRYSWKMBDHV-".Contains(c, StringComparison.Ordinal)))
            molecule = "DNA";
        else if (all.Contains('U', StringComparison.Ordinal) && all.All(c => "ACGUNRYSWKMBDHV-".Contains(c, StringComparison.Ordinal)))
            molecule = "RNA";
        else
            molecule = "protein";
        for (int i = 0; i < rows.Length; i++)
        {
            bool hasName = names is not null && i < names.Length && !string.IsNullOrWhiteSpace(names[i]);
            string name = hasName ? names![i] : "seq" + (i + 1).ToString(CultureInfo.InvariantCulture);
            string description = descriptions is not null && i < descriptions.Length ? descriptions[i] : "";
            var annotations = new System.Text.Json.Nodes.JsonObject();
            if (bio == AlignmentFileFormat.Nexus)
                annotations["molecule_type"] = molecule;
            records.Add(new BioAlignmentRecord
            {
                Id = name,
                Description = description ?? "",
                Row = rows[i] ?? throw new ArgumentException($"Row {i} is null.", nameof(rows)),
                Annotations = annotations,
            });
        }

        var column = new System.Text.Json.Nodes.JsonObject();
        if (bio == AlignmentFileFormat.A2m)
            column["state"] = new string('D', rows.Length > 0 ? rows[0].Length : 0);
        return AlignmentFormats.Write([new BioAlignment { Records = records, ColumnAnnotations = column }], bio);
    }

    /// <summary>
    /// Writer of the alignment formats this server does not know itself (the EMBOSS alignment formats pair, markx10, msf, ...): the EMBOSS servers install it at start-up
    /// (<c>EmbossSequenceFormatsTools.InstallAlignmentFormats</c>); without it the tools of this server write the Biopython formats only.
    /// </summary>
    public delegate string? ExternalAlignmentFormatWriter(string format, string[] rows, string[]? names, double score, int[]? offsets, int[]? offends, string program);

    /// <summary>The installed <see cref="ExternalAlignmentFormatWriter"/> (null = the Biopython formats only).</summary>
    public static ExternalAlignmentFormatWriter? ExternalFormatWriter { get; set; }

    public static string? FormatRowsForOutput(string? format, string[] rows, string[]? names, double score, int[]? offsets = null, int[]? offends = null, string program = "needle")
    {
        if (string.IsNullOrWhiteSpace(format))
            return null;
        string f = format.Trim();
        bool forceEmboss = false;
        if (f.StartsWith("emboss:", StringComparison.OrdinalIgnoreCase))
        {
            forceEmboss = true;
            f = f[7..];
        }

        if (!forceEmboss && AlignmentFormats.TryParseFormat(f, out var bio))
            return WriteBiopythonFormat(bio, rows, names, null);

        string? external = ExternalFormatWriter?.Invoke(f, rows, names, score, offsets, offends, program);
        if (external is not null)
            return external;
        throw new ArgumentException($"Unknown alignment format '{format}' (the EMBOSS alignment formats are written by the EMBOSS servers).", nameof(format));
    }

    // residues of the sequence before / after the aligned row: the start of the row clamped to what the row leaves of the sequence
    private static (int Offset, int Offend) RowFlanks(string row, int start, int length)
    {
        int consumed = 0;
        foreach (char c in row)
        {
            if (c is not ('-' or '.' or '~' or '?' or ' '))
                consumed++;
        }

        if (consumed >= length)
            return (0, 0);
        int offset = Math.Clamp(start, 0, length - consumed);
        return (offset, length - offset - consumed);
    }

    public static AlignmentResultDto FormatPair(AlignmentResultDto d, string? format, int length1, int length2, string[]? names, string? program = null, string? sequence1 = null, string? sequence2 = null, bool queryFirst = false, bool trimEndGaps = false)
    {
        if (string.IsNullOrWhiteSpace(format))
            return d;
        // R14-9: the coordinate formats (bed, psl, chain, maf, sam, exonerate, cigar, md, cs) of the pairwise tools that know their sequences: the bridge of the library, Biopython Alignment.format
        string f = format.Trim();
        string coordinate = f.StartsWith("bio:", StringComparison.OrdinalIgnoreCase) ? f[4..] : f;
        if (sequence1 is not null && sequence2 is not null && !f.StartsWith("emboss:", StringComparison.OrdinalIgnoreCase) && global::Seqeron.Genomics.Alignment.AlignmentResultCoordinates.IsFormatName(coordinate))
            return d with { Formatted = FormatCoordinates(d, coordinate, sequence1, sequence2, names, queryFirst, trimEndGaps) };
        if (f.StartsWith("bio:", StringComparison.OrdinalIgnoreCase))
            format = coordinate;
        double score = d.ExactScore != 0 ? d.ExactScore : d.Score;
        // the residues before / after each row come from the rows, not from the end coordinates
        var (offset1, offend1) = RowFlanks(d.AlignedSequence1, d.StartPosition1, length1);
        var (offset2, offend2) = RowFlanks(d.AlignedSequence2, d.StartPosition2, length2);
        int[] offsets = [offset1, offset2];
        int[] offends = [offend1, offend2];
        // the program whose report layout the alignment has: a local alignment is a water report (its SAM line is a local one, -aglobal N), every other one a needle report
        program ??= string.Equals(d.AlignmentType, "Local", StringComparison.OrdinalIgnoreCase) ? "water" : "needle";
        return d with { Formatted = FormatRowsForOutput(format, [d.AlignedSequence1, d.AlignedSequence2], names, score, offsets, offends, program) };
    }

    // the result as an AlignmentResult of the library and the text of one coordinate format; an alignment without a column (an empty local result) has no record to write
    private static string FormatCoordinates(AlignmentResultDto d, string format, string sequence1, string sequence2, string[]? names, bool queryFirst, bool trimEndGaps)
    {
        if (d.AlignedSequence1.Length == 0)
            return "";
        var type = Enum.TryParse<global::Seqeron.Genomics.Infrastructure.AlignmentType>(d.AlignmentType, true, out var t) ? t : global::Seqeron.Genomics.Infrastructure.AlignmentType.Global;
        var result = new global::Seqeron.Genomics.Infrastructure.AlignmentResult(d.AlignedSequence1, d.AlignedSequence2, d.Score, type, d.StartPosition1, d.StartPosition2, d.EndPosition1, d.EndPosition2) { ExactScore = d.ExactScore != 0 ? d.ExactScore : d.Score };
        string name1 = names is not null && names.Length > 0 && !string.IsNullOrWhiteSpace(names[0]) ? names[0] : "seq1";
        string name2 = names is not null && names.Length > 1 && !string.IsNullOrWhiteSpace(names[1]) ? names[1] : "seq2";
        var options = new global::Seqeron.Genomics.Alignment.AlignmentBridgeOptions { FirstIsQuery = queryFirst, TrimEndGaps = trimEndGaps };
        return global::Seqeron.Genomics.Alignment.AlignmentResultCoordinates.Format(result, format, name1, sequence1.ToUpperInvariant(), name2, sequence2.ToUpperInvariant(), options);
    }
}
