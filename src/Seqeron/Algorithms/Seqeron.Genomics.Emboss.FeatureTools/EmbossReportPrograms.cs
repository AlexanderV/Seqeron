// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.FeatureTools, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.FeatureTools;

using Seqeron.Genomics.Analysis;

/// <summary>What the report of a program run needs besides the program's own parameters (audit A5-8b, B09 F184).</summary>
/// <param name="Fasta">The input file: one or more FASTA records (etandem and equicktandem read the first only).</param>
/// <param name="SequenceQualifiers">The sequence qualifiers as typed (<c>-sbegin 3 -send -2 -sreverse</c>, see <see cref="EmbossSequenceQualifiers.ParseInputsAndReport"/>), may also carry the report qualifiers.</param>
/// <param name="ReportQualifiers">The report qualifiers as typed (<c>-rformat table -raccshow -rmaxall 3</c>); null = the program's defaults.</param>
/// <param name="CommandLine">The command-line entries printed in the report header (<c>[-sequence] s.fa</c>, <c>-pattern ACGT</c>, ...: <see cref="EmbossCommandLineEcho"/>); null = none.</param>
/// <param name="RunDate">The run date of the header (<c>Rundate</c>, <c>#!Date</c>); null = now.</param>
/// <param name="FastaFile">The input file name (the USA of a sequence is <c>fasta::file:name</c>).</param>
/// <param name="ReportFile">The report file name (<c>Report_file:</c>).</param>
/// <param name="PatternName">The name of the pattern (<c>-pname</c>, default <c>pattern</c>).</param>
/// <param name="CancellationToken">Cancellation of the searches.</param>
public sealed record EmbossReportRequest(string Fasta, string? SequenceQualifiers = null, string? ReportQualifiers = null, IReadOnlyList<string>? CommandLine = null, DateTime? RunDate = null,
    string FastaFile = "s.fa", string ReportFile = "stdout", string PatternName = "pattern", CancellationToken CancellationToken = default);

/// <summary>The files a report program writes.</summary>
/// <param name="Report">The <c>-outfile</c> report text.</param>
/// <param name="OriginalFile">The <c>-origfile</c> text of etandem / equicktandem (their original output), else null.</param>
public sealed record EmbossProgramReport(string Report, string? OriginalFile);

/// <summary>The command-line echo of the report header (<c>ajacd.c</c> <c>acdArgsParse</c>): how the program writes the options it was given.</summary>
public static class EmbossCommandLineEcho
{
    /// <summary>The first parameter: <c>[-sequence] value</c> (quoted when it is not a word or holds one of <c>* ? [ ] { } | ! &amp; ^</c>).</summary>
    public static string Parameter(string name, string value) => "[-" + name + "] " + Quote(value);

    /// <summary>A qualifier as typed: a boolean given Y (or its <c>no</c> form given N) prints without a value.</summary>
    public static string Qualifier(string name, string? value, bool boolean = false)
    {
        if (value is null)
            return "-" + name;
        if (boolean && string.Equals(value, "Y", StringComparison.OrdinalIgnoreCase))
            return "-" + name;
        return "-" + name + " " + Quote(value);
    }

    public static string Quote(string v) => v.Length == 0 || v.Any(c => char.IsWhiteSpace(c) || "*?[]{}|!&^".Contains(c, StringComparison.Ordinal)) ? "\"" + v + "\"" : v;
}

/// <summary>
/// The report files of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem as the programs write them (<c>fuzznuc.c</c> ... <c>equicktandem.c</c>, EMBOSS-6.6.0): the search of the ported program, the feature tables of
/// <see cref="EmbossReports"/> and the one writer <see cref="EmbossReportWriter"/>, with every <c>-rformat</c>, the report options (<c>-raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq</c>), several sequences per file
/// and the extra <c>-origfile</c> of the tandem programs. Compared with the installed programs on TestData/Alignment/emboss_report_formats_oracle.jsonl.gz (audit A5-8b, B09 F184).
/// </summary>
public static class EmbossReportPrograms
{
    public static string Int(long v) => v.ToString(CultureInfo.InvariantCulture);

    public static string Bool(bool b) => b ? "Yes" : "No";

    public sealed record Setup(EmbossReportOptions Options, IReadOnlyList<EmbossSequenceQualifiers> Inputs);

    public static Setup Prepare(string program, EmbossReportRequest r, string defaultFormat, string header, bool footer, bool emptyTables)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Fasta);
        string text = (r.SequenceQualifiers ?? "") + " " + (r.ReportQualifiers ?? "");
        var (inputs, rq) = EmbossSequenceQualifiers.ParseInputsAndReport(program, text);
        string formatText = rq.Format ?? defaultFormat;
        if (!EmbossReportFormats.TryParse(formatText, out var format))
            throw new EmbossReportException($"Unknown report format '{formatText}'");
        var options = new EmbossReportOptions
        {
            Program = program,
            CommandLine = r.CommandLine ?? [],
            RunDate = r.RunDate ?? DateTime.Now,
            ReportFile = r.ReportFile,
            Format = format,
            FormatText = formatText,
            Tags = EmbossReportTagLists.TagList(program),
            Header = header,
            ShowScore = rq.ShowScore ?? true,
            ShowStrand = rq.ShowStrand ?? true,
            ShowAccession = rq.ShowAccession,
            ShowDescription = rq.ShowDescription,
            ShowUsa = rq.ShowUsa,
            MaxAll = rq.MaxAll,
            MaxSeq = rq.MaxSeq,
            Precision = 0,
            WriteFooter = footer,
            WriteEmptyTables = emptyTables,
        };
        return new Setup(options, inputs);
    }

    // the USA of a sequence read with a range: the qualifiers as typed
    public static string Usa(EmbossReportRequest r, string name, EmbossSequenceQualifiers q)
        => "fasta::" + r.FastaFile + ":" + name + UsaRange(q);

    /// <summary>The range part of the USA of a sequence read with <c>-sbegin</c> / <c>-send</c> / <c>-sreverse</c> (<c>[b:e:r]</c>), empty without them.</summary>
    public static string UsaRange(EmbossSequenceQualifiers q)
    {
        if (q.Begin == 0 && q.End == 0 && !q.Reverse)
            return "";
        string begin = q.Begin != 0 ? Int(q.Begin) : "";
        string end = q.End != 0 ? Int(q.End) : "";
        return "[" + begin + ":" + end + (q.Reverse ? ":r" : "") + "]";
    }

    // ajSeqGetBegin / ajSeqGetEnd: unset = whole, negative = from the end
    public static (int Begin, int End) Range(EmbossSequenceView v)
    {
        int len = v.Full.Length;
        return (Math.Max(1, Position(v.Begin, 1, len)), Math.Min(len, Position(v.End, len, len)));
    }

    public static int Position(int given, int unset, int len)
    {
        if (given == 0)
            return unset;
        return given < 0 ? len + given + 1 : given;
    }

    public static EmbossReportSequence ReportSequence(EmbossReportRequest r, EmbossFastaRecord rec, EmbossSequenceView view, EmbossSequenceQualifiers q, bool nucleotide)
    {
        var (b, e) = Range(view);
        return new EmbossReportSequence(rec.Name, Usa(r, rec.Name, q), rec.Accession, rec.Description, view.Full, b, e, 0, view.Reversed, nucleotide);
    }

    public static string Finish(Setup s, List<EmbossReport> tables, List<int> seqIndex, IReadOnlyList<long> lengths, int count)
    {
        return EmbossReportWriter.WriteFile(s.Options, tables, stoppedAt =>
        {
            int read = stoppedAt < 0 ? count : seqIndex[stoppedAt] + 1;
            return (read, lengths.Take(read).Sum());
        });
    }
}
