// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.FeatureTools, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.FeatureTools;

using Seqeron.Genomics.Analysis;

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

    /// <summary>The report file is not named (no <c>-outfile</c>): it is <c>BASE.featreport</c> with the base of <c>-rname</c> or of the name of the sequence and the extension of <c>-rextension</c> (<see cref="ReportFile"/> is not used).</summary>
    public bool DefaultReportFile { get; init; }

    /// <summary>The directories (without trailing slash) that exist where the program runs; a report file whose directory is not among them cannot be opened. Null: every directory exists.</summary>
    public IReadOnlyCollection<string>? ExistingDirectories { get; init; }
}

/// <summary>The output of <c>featreport</c> and <c>twofeat</c>.</summary>
/// <param name="Report">The text of the report file.</param>
/// <param name="Messages">What the library printed on its error stream.</param>
public sealed record EmbossFeatureReportOutput(string Report, IReadOnlyList<string> Messages)
{
    /// <summary>The program died ("Died: ...").</summary>
    public bool Died { get; init; }

    /// <summary>The path of the report file (<c>-outfile</c>, or the default name built from <c>-rname</c> / <c>-rextension</c> and put in <c>-rdirectory</c>); empty when the program died before the report was opened.</summary>
    public string ReportFile { get; init; } = "";
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
            var fq = r.Report?.FeatureInput;
            var input = new EmbossFeatureInput(r.Features, name => r.Files.TryGetValue(name, out string? t) ? t : null, fq?.Begin ?? r.FeatureBegin, fq?.End ?? r.FeatureEnd, context, "features", 1);
            table = input.Next();
        }
        catch (EmbossFeatureDiedException e)
        {
            return new EmbossFeatureReportOutput("", [.. e.Errors, e.Message]) { Died = true };
        }

        var seq = seqall.Peek();
        // acdSetReport: the report is validated (a format that needs more tags than the program has stops the program) and then its file is opened, after the inputs
        string defaultName = EmbossOutPath.DefaultBaseName(seq.Name, r.Sequence);
        var options = EmbossFeatureReports.Options("featreport", "gff", r.Report, r.CommandLine, r.RunDate, r.ReportFile, context, false);
        if (EmbossFeatureReports.ValidateTags(options, context))
            return new EmbossFeatureReportOutput("", [.. context.Messages, "Died: Report option -outfile: Validation failed"]) { Died = true };
        string? reportFile = EmbossOutPath.OpenReport(context, "featreport", r.DefaultReportFile ? null : r.ReportFile, r.Report, defaultName, r.ExistingDirectories);
        if (reportFile is null)
        {
            var died = EmbossSeqretPrograms.Died(context, "outfile");
            return new EmbossFeatureReportOutput("", died.Messages) { Died = true };
        }

        options = options with { ReportFile = reportFile };
        var tables = new List<EmbossReport>();
        if (table is not null)
            tables.Add(EmbossFeatureReports.Build(seq, table, options.Tags));
        string text = EmbossFeatureReports.Write(options, tables, context, null, out _);
        return new EmbossFeatureReportOutput(text, [.. context.Messages]) { ReportFile = reportFile };
    }
}
