using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The adapter between the ACD command line of the stream runner (<see cref="EmbossAcdArgs"/>) and the report layer (<see cref="EmbossReportWriter"/>): <c>acdSetReport</c> of an ACD item <c>report</c>
/// (<c>-outfile</c>, its associated qualifiers <c>-rformat -rname -rextension -rdirectory -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq</c>, the tag list, the precision and the default
/// format of the program), the opening of the report file and the writing of the tables (B09 WP63, F282).
/// </summary>
internal static class EmbossReportArgs
{
    /// <summary>The associated qualifiers of the report item <paramref name="master"/> as typed.</summary>
    /// <exception cref="EmbossAcdDieException">A value the library refuses.</exception>
    public static EmbossReportQualifiers Qualifiers(EmbossAcdArgs a, string master) => new(
        Format: a.RawAssoc(master, "rformat"),
        ShowAccession: a.BoolAssoc(master, "raccshow") ?? false,
        ShowDescription: a.BoolAssoc(master, "rdesshow") ?? false,
        ShowScore: a.BoolAssoc(master, "rscoreshow"),
        ShowStrand: a.BoolAssoc(master, "rstrandshow"),
        ShowUsa: a.BoolAssoc(master, "rusashow") ?? false,
        MaxAll: a.IntAssoc(master, "rmaxall") ?? 0,
        MaxSeq: a.IntAssoc(master, "rmaxseq") ?? 0,
        Name: a.RawAssoc(master, "rname"),
        Extension: a.RawAssoc(master, "rextension"),
        Directory: a.RawAssoc(master, "rdirectory"));

    /// <summary>The directories of the work directory of the run: the directories of the input files (the only ones that exist).</summary>
    internal static IReadOnlyCollection<string> ExistingDirectories(IReadOnlyDictionary<string, string> files)
    {
        var dirs = new HashSet<string>(StringComparer.Ordinal);
        foreach (string name in files.Keys)
        {
            for (int slash = name.IndexOf('/', StringComparison.Ordinal); slash > 0; slash = name.IndexOf('/', slash + 1))
                dirs.Add(name[..slash]);
        }

        return dirs;
    }

    /// <summary>The report of an open ACD item: its options and the path of its file.</summary>
    internal sealed record Opened(EmbossReportOptions Options, string File);

    /// <summary>
    /// <c>acdSetReport</c>: the report qualifiers are read (a bad value dies), the format must exist and the tag list must have as many tags as the format needs, then the file is opened (<c>-outfile</c> as given, else
    /// <c>BASE.program</c> of <c>-rname</c> / the name of the first sequence; <c>-rdirectory</c> in front of it).
    /// </summary>
    /// <returns>The opened report; null when the program dies (the messages are in <paramref name="context"/>, the Died line is the caller's: <see cref="EmbossSeqretPrograms.Died(EmbossFeatureReadContext, string)"/>).</returns>
    public static Opened? Open(EmbossAcdArgs a, EmbossFeatureReadContext context, string master, string program, string defaultFormat, IReadOnlyList<EmbossReportTag> tags, int precision, string header, string? given, string baseName, DateTime runDate, IReadOnlyDictionary<string, string> files, bool multi = false)
    {
        var rq = Qualifiers(a, master);
        string formatText = rq.Format ?? defaultFormat;
        if (!EmbossReportFormats.TryParse(formatText, out var format))
        {
            context.Error($"Unknown report format '{formatText}'");
            Validation(context);
            return null;
        }

        var options = new EmbossReportOptions
        {
            Program = program,
            CommandLine = a.CommandLineEntries,
            RunDate = runDate,
            ReportFile = "stdout",
            Format = format,
            FormatText = formatText,
            Tags = tags,
            Header = header,
            ShowScore = rq.ShowScore ?? true,
            ShowStrand = rq.ShowStrand ?? true,
            ShowAccession = rq.ShowAccession,
            ShowDescription = rq.ShowDescription,
            ShowUsa = rq.ShowUsa,
            MaxAll = rq.MaxAll,
            MaxSeq = rq.MaxSeq,
            Precision = precision,
            WriteFooter = true,
            WriteEmptyTables = false,
            Context = context,
            Multi = multi,
        };
        if (EmbossFeatureReports.ValidateTags(options, context))
        {
            Validation(context);
            return null;
        }

        string? file = EmbossOutPath.OpenReport(context, program, given, rq, baseName, ExistingDirectories(files));
        if (file is null)
            return null;
        return new Opened(options with { ReportFile = file }, file);
    }

    // the ACD refuses the report: "Died: Report option -outfile: Validation failed" (no further death message)
    private static void Validation(EmbossFeatureReadContext context)
    {
        context.Messages.Add("Died: Report option -outfile: Validation failed");
        context.Fatal = true;
    }

    /// <summary>
    /// Writes the report file of the run (<c>ajReportWrite</c> of every table, <c>ajReportClose</c>); the totals of the footer are the sequences read (<c>ajReportSetSeqstats</c>); a report without a hit is the warning
    /// "No hits in output report file" of <c>ajReportClose</c>.
    /// </summary>
    public static string Write(Opened report, IReadOnlyList<EmbossReport> tables, EmbossFeatureReadContext context, long sequences, long length)
    {
        var options = report.Options;
        string text = EmbossFeatureReports.Write(options, tables, context, _ => (sequences, length), out bool failed);
        long reported = 0;
        foreach (var table in tables)
        {
            long n = table.Features.Count;
            if (options.MaxSeq != 0)
                n = Math.Min(n, options.MaxSeq);
            if (options.MaxAll != 0)
            {
                long remaining = options.MaxAll - reported;
                if (remaining <= 0)
                    break;
                n = Math.Min(n, remaining);
            }

            reported += n;
        }

        if (reported == 0 || failed)
            context.Warn($"No hits in output report file '{report.File}'");
        return text;
    }

    /// <summary>The result of a run that wrote a report file (<c>stdout</c> goes to the standard output) and perhaps more files.</summary>
    public static EmbossRunResult Result(EmbossFeatureReadContext context, string file, string text, IReadOnlyDictionary<string, string>? more = null)
    {
        var files = new Dictionary<string, string>();
        if (more is not null)
            foreach (var (k, v) in more)
                files[k] = v;
        string stdout = "";
        if (file == "stdout")
            stdout = text;
        else
            files[file] = text;
        return new EmbossRunResult(stdout, EmbossRunResult.Join(context.Messages), 0, files);
    }
}
