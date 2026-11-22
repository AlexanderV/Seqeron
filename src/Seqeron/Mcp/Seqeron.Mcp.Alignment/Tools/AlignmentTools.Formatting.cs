using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-27b: the `format` parameter of the pairwise and multiple-alignment tools, mapped to the canonical alignment writers (AlignmentFormats = Biopython Bio.Align, EmbossAlignmentWriter = EMBOSS -aformat).
public partial class AlignmentTools
{
    private const string FormatParamDescription = "Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text.";

    private const string NamesParamDescription = "Optional names of the rows for the formatted output (default seq1, seq2, ...).";

    // the Biopython Bio.Align writers (AlignmentFormats) for rows; a2m gets all-match states, nexus the molecule type guessed from the residues
    private static string WriteBiopythonFormat(AlignmentFileFormat bio, string[] rows, string[]? names, string[]? descriptions)
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

    private static string? FormatRowsForOutput(string? format, string[] rows, string[]? names, double score, int[]? offsets = null, int[]? offends = null, string program = "needle")
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

        return WriteAlignment(rows, f, names, null, offsets, offends, score, program: program).Text;
    }

    private static AlignmentResultDto FormatPair(AlignmentResultDto d, string? format, int length1, int length2, string[]? names, string? program = null)
    {
        if (string.IsNullOrWhiteSpace(format))
            return d;
        double score = d.ExactScore != 0 ? d.ExactScore : d.Score;
        // R14-1: the residues before / after each row come from the rows (EmbossAlignedSequence.Flanks), not from the end coordinates
        var (offset1, offend1) = EmbossAlignedSequence.Flanks(d.AlignedSequence1, d.StartPosition1, length1);
        var (offset2, offend2) = EmbossAlignedSequence.Flanks(d.AlignedSequence2, d.StartPosition2, length2);
        int[] offsets = [offset1, offset2];
        int[] offends = [offend1, offend2];
        // the program whose report layout the alignment has: a local alignment is a water report (its SAM line is a local one, -aglobal N), every other one a needle report
        program ??= string.Equals(d.AlignmentType, "Local", StringComparison.OrdinalIgnoreCase) ? "water" : "needle";
        return d with { Formatted = FormatRowsForOutput(format, [d.AlignedSequence1, d.AlignedSequence2], names, score, offsets, offends, program) };
    }

    private static MultipleAlignmentResultDto FormatMultiple(MultipleAlignmentResultDto d, string? format, string[]? names)
    {
        if (string.IsNullOrWhiteSpace(format))
            return d;
        double score = d.ExactTotalScore != 0 ? d.ExactTotalScore : d.TotalScore;
        return d with { Formatted = FormatRowsForOutput(format, d.AlignedSequences, names, score) };
    }
}
