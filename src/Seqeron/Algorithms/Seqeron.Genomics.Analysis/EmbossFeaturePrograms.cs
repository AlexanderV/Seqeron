using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>A run of <c>featcopy</c>.</summary>
/// <param name="Features">The value of <c>-features</c>: a file name with an optional <c>format::</c> prefix and <c>[begin:end:r]</c> suffix (the <c>-fformat</c> qualifier of the program has no effect).</param>
/// <param name="Files">The input files by name: a file that is not here does not exist.</param>
/// <param name="FeatureBegin"><c>-fbegin</c> (default 1: the first table is cut to the range, negative counts from the end).</param>
/// <param name="FeatureEnd"><c>-fend</c> (default: the length of the first table).</param>
/// <param name="OutFeat"><c>-outfeat</c>: the output file name, optionally with a <c>format::</c> prefix (default <c>unknown.gff</c>).</param>
/// <param name="OutFormat"><c>-offormat</c>: the output feature format (default <c>gff</c> = GFF3).</param>
/// <param name="OutOpenFile"><c>-ofopenfile</c>: the file the program opens instead of the one named by <c>-outfeat</c>.</param>
/// <param name="RunDate">The date of the run (GFF headers); null = now.</param>
public sealed record EmbossFeatcopyRequest(string Features, IReadOnlyDictionary<string, string> Files, int FeatureBegin = 1, int? FeatureEnd = null, string OutFeat = "unknown.gff", string? OutFormat = null, string? OutOpenFile = null, DateTime? RunDate = null);

/// <summary>The file <c>featcopy</c> wrote.</summary>
/// <param name="FileName">The name of the output file.</param>
/// <param name="Text">The tables in the output format (empty when nothing was written).</param>
/// <param name="Format">The output feature format as the program resolved it.</param>
/// <param name="Tables">The number of feature tables written.</param>
/// <param name="Messages">What the library printed on its error stream (<c>Warning: ...</c>, <c>Error: ...</c>), in order.</param>
public sealed record EmbossFeatcopyOutput(string FileName, string Text, string Format, int Tables, IReadOnlyList<string> Messages);

/// <summary>A run of <c>featmerge</c>.</summary>
/// <param name="AFeatures">The value of <c>-afeatures</c> (file with optional <c>format::</c> prefix and <c>[begin:end:r]</c> suffix).</param>
/// <param name="BFeatures">The value of <c>-bfeatures</c>.</param>
/// <param name="Files">The input files by name.</param>
/// <param name="AFeatureBegin"><c>-fbegin1</c>: the range of the first table (the merged table keeps it, the second table's range is not used).</param>
/// <param name="AFeatureEnd"><c>-fend1</c>.</param>
/// <param name="BFeatureBegin"><c>-fbegin2</c>.</param>
/// <param name="BFeatureEnd"><c>-fend2</c>.</param>
/// <param name="OutFeat"><c>-outfeat</c>: the output file name, optionally with a <c>format::</c> prefix (default <c>unknown.gff</c>).</param>
/// <param name="OutFormat"><c>-offormat</c>: the output feature format.</param>
/// <param name="OutOpenFile"><c>-ofopenfile</c>.</param>
/// <param name="RunDate">The date of the run (GFF headers); null = now.</param>
public sealed record EmbossFeatmergeRequest(string AFeatures, string BFeatures, IReadOnlyDictionary<string, string> Files, int AFeatureBegin = 1, int? AFeatureEnd = null, int BFeatureBegin = 1, int? BFeatureEnd = null,
    string OutFeat = "unknown.gff", string? OutFormat = null, string? OutOpenFile = null, DateTime? RunDate = null);

/// <summary>
/// The feature-table programs of EMBOSS 6.6.0 built on the feature library (<c>ajfeatread.c</c> readers, <c>ajfeatwrite.c</c> writers): <c>featcopy</c> (A5-8c3e, B09 F241): every feature table of a file in any input format (embl, genbank, refseq, genpept, refseqp,
/// gff2, gff3, the gff3 of EMBOSS before 6.4, swiss, pir, dasgff; with no format the library tries each in turn) is written in an output format (every output format of the library: gff / gff3, gff2, bed, embl, genbank, refseq, refseqp, pir, swiss, dasgff, draw, debug; see <see cref="EmbossReportWriter.FeatureTableFormats"/>).
/// </summary>
public static class EmbossFeaturePrograms
{
    /// <summary>Runs <c>featcopy</c>.</summary>
    /// <exception cref="EmbossFeatureDiedException">The input is no feature file ("Bad value for '-features' and no more retries").</exception>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops (a GFF2 group field its tag pattern cannot parse, a location that begins with a comma, a reader that returns tables for ever, a Swiss-Prot feature line the wrapper cannot break).</exception>
    public static EmbossFeatcopyOutput Featcopy(EmbossFeatcopyRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Features);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "featcopy" };
        var (index, file) = EmbossDiffseqPrograms.OpenFeatureOutput(r.OutFeat, new EmbossFeatureOutQualifiers(r.OutFormat, r.OutOpenFile));
        var input = new EmbossFeatureInput(r.Features, name => r.Files.TryGetValue(name, out string? t) ? t : null, r.FeatureBegin, r.FeatureEnd, context);
        string format = EmbossReportWriter.FeatureFormatName(index);
        var state = new EmbossFeatureWriteState();
        var sb = new System.Text.StringBuilder();
        int tables = 0;
        EmbossSeqFeatureTable? table;
        while ((table = input.Next()) is not null)
        {
            string? text = EmbossReportWriter.WriteFeatureTable(table, format, r.RunDate ?? DateTime.Now, "featcopy", context, state, tables);
            if (text is null)
                continue;
            sb.Append(text);
            tables++;
        }

        if (tables == 0)
            context.Warn($"No features written to output file '{file}'");
        return new EmbossFeatcopyOutput(file, sb.ToString(), format, tables, [.. context.Messages]);
    }

    /// <summary>
    /// Runs <c>featmerge</c> (<c>emboss/featmerge.c</c>): the first table of <c>-afeatures</c> and the first table of <c>-bfeatures</c> are read (each with its range), the features of the second are copied into the first
    /// (<c>ajFeattableMerge</c>: <c>ajFeatNewFeat</c> + <c>ajFeattableAdd</c>, the group numbers are kept) and the first is written in the output format (<c>ajFeattableWrite</c>).
    /// </summary>
    /// <exception cref="EmbossFeatureDiedException">An input is no feature file ("Bad value for '-afeatures' and no more retries").</exception>
    public static EmbossFeatcopyOutput Featmerge(EmbossFeatmergeRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.AFeatures);
        ArgumentNullException.ThrowIfNull(r.BFeatures);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "featmerge" };
        Func<string, string?> open = name => r.Files.TryGetValue(name, out string? t) ? t : null;
        var (index, file) = EmbossDiffseqPrograms.OpenFeatureOutput(r.OutFeat, new EmbossFeatureOutQualifiers(r.OutFormat, r.OutOpenFile));
        var a = new EmbossFeatureInput(r.AFeatures, open, r.AFeatureBegin, r.AFeatureEnd, context, "afeatures", 1).Next()!;
        var b = new EmbossFeatureInput(r.BFeatures, open, r.BFeatureBegin, r.BFeatureEnd, context, "bfeatures", 1).Next()!;
        a.Merge(b);
        string format = EmbossReportWriter.FeatureFormatName(index);
        string? text = EmbossReportWriter.WriteFeatureTable(a, format, r.RunDate ?? DateTime.Now, "featmerge", context, new EmbossFeatureWriteState(), 0);
        if (text is null)
            context.Warn($"No features written to output file '{file}'");
        return new EmbossFeatcopyOutput(file, text ?? "", format, text is null ? 0 : 1, [.. context.Messages]);
    }
}
