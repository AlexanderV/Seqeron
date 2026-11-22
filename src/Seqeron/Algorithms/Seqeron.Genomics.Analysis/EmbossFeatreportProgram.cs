using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>A run of <c>featreport</c> (A5-8c3h, B09 F262): the parameters of its ACD file.</summary>
/// <param name="Sequence">The value of <c>-sequence</c> (the USA).</param>
/// <param name="Features">The value of <c>-features</c>: a feature file with an optional <c>format::</c> prefix and range.</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers.</param>
/// <param name="Report">The report qualifiers (<c>-rformat</c>, <c>-raccshow</c> ...); null = the defaults (format gff).</param>
/// <param name="FeatureBegin"><c>-fbegin2</c>.</param>
/// <param name="FeatureEnd"><c>-fend2</c>.</param>
/// <param name="CommandLine">The entries of the command line the report header echoes (<see cref="EmbossCommandLineEcho"/>).</param>
/// <param name="RunDate">The date of the run; null = now.</param>
/// <param name="ReportFile">The name of the report file (<c>Report_file:</c>).</param>
public sealed record EmbossFeatreportRequest(string Sequence, string Features, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null, EmbossReportQualifiers? Report = null,
    int FeatureBegin = 1, int? FeatureEnd = null, IReadOnlyList<string>? CommandLine = null, DateTime? RunDate = null, string ReportFile = "stdout")
{
    /// <summary>The texts of the configuration files that define the databases of <c>db:id</c> USAs.</summary>
    public IReadOnlyList<string>? Definitions { get; init; }
}

/// <summary>The output of <c>featreport</c> and <c>twofeat</c>.</summary>
/// <param name="Report">The text of the report file.</param>
/// <param name="Messages">What the library printed on its error stream.</param>
public sealed record EmbossFeatureReportOutput(string Report, IReadOnlyList<string> Messages)
{
    /// <summary>The program died ("Died: ...").</summary>
    public bool Died { get; init; }
}

/// <summary>
/// <c>featreport</c> of EMBOSS 6.6.0 (<c>emboss/featreport.c</c>): the first feature table of a feature file reported for a sequence (<c>ajReportWrite</c>): the report formats that are feature formats (gff, gff3, embl, genbank, pir, swiss, dasgff) write the table
/// with the feature writers of the library, the others list the start, end, strand and score of its features.
/// </summary>
public static class EmbossFeatreportProgram
{
    /// <summary>Runs <c>featreport</c>.</summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    public static EmbossFeatureReportOutput Featreport(EmbossFeatreportRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Features);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "featreport" };
        var seqall = EmbossSeqall.Open("featreport", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, false, context, "any", true, r.Definitions);
        if (seqall is null)
        {
            var died = EmbossSeqretPrograms.Died(context);
            return new EmbossFeatureReportOutput("", died.Messages) { Died = true };
        }

        EmbossSeqFeatureTable? table;
        try
        {
            var input = new EmbossFeatureInput(r.Features, name => r.Files.TryGetValue(name, out string? t) ? t : null, r.FeatureBegin, r.FeatureEnd, context, "features", 1);
            table = input.Next();
        }
        catch (EmbossFeatureDiedException e)
        {
            return new EmbossFeatureReportOutput("", [.. e.Errors, e.Message]) { Died = true };
        }

        var seq = seqall.Peek();
        var options = EmbossFeatureReports.Options("featreport", "gff", r.Report, r.CommandLine, r.RunDate, r.ReportFile, context, false);
        if (EmbossFeatureReports.ValidateTags(options, context))
            return new EmbossFeatureReportOutput("", [.. context.Messages, "Died: Report option -outfile: Validation failed"]) { Died = true };
        var tables = new List<EmbossReport>();
        if (table is not null)
            tables.Add(EmbossFeatureReports.Build(seq, table, options.Tags));
        string text = EmbossFeatureReports.Write(options, tables, context, null);
        return new EmbossFeatureReportOutput(text, [.. context.Messages]);
    }
}

/// <summary>The reports of the programs that report the feature tables of the feature library (<c>ajReportWrite</c> with a table and a sequence).</summary>
internal static class EmbossFeatureReports
{
    /// <summary>The options of a report as <c>acdSetReport</c> sets them.</summary>
    public static EmbossReportOptions Options(string program, string defaultFormat, EmbossReportQualifiers? rq, IReadOnlyList<string>? commandLine, DateTime? runDate, string reportFile, EmbossFeatureReadContext context, bool footer)
    {
        rq ??= new EmbossReportQualifiers();
        string formatText = rq.Format ?? defaultFormat;
        if (!EmbossReportFormats.TryParse(formatText, out var format))
            throw new EmbossReportException($"Unknown report format '{formatText}'");
        return new EmbossReportOptions
        {
            Program = program,
            CommandLine = commandLine ?? [],
            RunDate = runDate ?? DateTime.Now,
            ReportFile = reportFile,
            Format = format,
            FormatText = formatText,
            Tags = EmbossReports.TagList(program),
            ShowScore = rq.ShowScore ?? true,
            ShowStrand = rq.ShowStrand ?? true,
            ShowAccession = rq.ShowAccession,
            ShowDescription = rq.ShowDescription,
            ShowUsa = rq.ShowUsa,
            MaxAll = rq.MaxAll,
            MaxSeq = rq.MaxSeq,
            Precision = 3,
            WriteFooter = footer,
            WriteEmptyTables = false,
            Context = context,
        };
    }

    /// <summary>The ACD validation of the report (<c>acdSetReport</c>): a format that needs more tags than the program has (diffseq) is an error and the program dies.</summary>
    public static bool ValidateTags(EmbossReportOptions options, EmbossFeatureReadContext context)
    {
        var flags = EmbossReportFormats.Flags(options.Format);
        if (flags.MinTags <= options.Tags.Count)
            return false;
        context.Error($"Report format '{options.FormatText ?? EmbossReportFormats.Name(options.Format)}' needs {flags.MinTags} tags, has only {options.Tags.Count}");
        return true;
    }

    /// <summary>The report table of a feature table: the full table for the feature formats, the flat list of its top-level features for the others.</summary>
    public static EmbossReport Build(EmbossSeq seq, EmbossSeqFeatureTable table, IReadOnlyList<EmbossReportTag> tags)
    {
        var features = table.Features.Select(f => new EmbossFeature(f.Type, f.Source, f.GetStart(), f.GetEnd(), f.Score, f.Strand, f.Frame, [.. f.Tags.Select(t => new EmbossFeatureTag(t.Name, t.Value))],
            f.Flags.HasFlag(EmbossLocationBits.BetweenSeq), f.Protein) { FromTable = true }).ToList();
        var sequence = new EmbossReportSequence(seq.Name, MakeUsa(seq), seq.Acc, seq.Desc, seq.Seq.ToString(), EmbossSeqRange.GetBegin(seq), EmbossSeqRange.GetEnd(seq), seq.Offset, seq.Reversed, seq.Type == "N");
        bool protein = table.Type == "P" || (table.Type.Length == 0 && seq.Type == "P");            // ajReportSetType: a table with no type takes the type of the sequence
        return new EmbossReport(sequence, features, protein, table.Type, unchecked((int)table.Start), unchecked((int)table.End), unchecked((int)table.Offset), table.Length, (int)table.Groups, table.FormatName ?? "<null>", tags)
        {
            FeatureTable = table,
            TypeUnknown = table.Type.Length == 0,
        };
    }

    /// <summary>Writes the report file; a format the table is not written in is an error of the library ("Report format ... not supported ...") and the file stays empty.</summary>
    public static string Write(EmbossReportOptions options, IReadOnlyList<EmbossReport> tables, EmbossFeatureReadContext context, Func<int, (long Sequences, long Length)>? totals)
    {
        try
        {
            return EmbossReportWriter.WriteFile(options, tables, totals);
        }
        catch (EmbossReportException e)
        {
            context.Error(e.Message);
            return "";
        }
    }

    // seqMakeUsa: the USA of a sequence (format::file:entry and the range)
    private static string MakeUsa(EmbossSeq t)
    {
        var sb = new System.Text.StringBuilder();
        if (t.Db.Length > 0)
            sb.Append(t.Db).Append("-id:").Append(t.Entryname);
        else if (t.Entryname.Length > 0)
            sb.Append(t.Formatstr).Append("::").Append(t.Filename).Append(':').Append(t.Entryname);
        else
            sb.Append(t.Formatstr).Append("::").Append(t.Filename);
        var range = new System.Text.StringBuilder("[");
        if (t.Reversed)
        {
            if (t.End != 0)
                range.Append(-t.End);
            range.Append(':');
            if (t.Begin != 0)
                range.Append(-t.Begin);
            range.Append(":r");
        }
        else
        {
            if (t.Begin != 0)
                range.Append(t.Begin);
            range.Append(':');
            if (t.End != 0)
                range.Append(t.End);
        }

        range.Append(']');
        if (range.Length > 3)
            sb.Append(range);
        return sb.ToString();
    }
}
