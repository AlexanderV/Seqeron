using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

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

    private static string Quote(string v) => v.Length == 0 || v.Any(c => char.IsWhiteSpace(c) || "*?[]{}|!&^".Contains(c, StringComparison.Ordinal)) ? "\"" + v + "\"" : v;
}

/// <summary>
/// The report files of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem as the programs write them (<c>fuzznuc.c</c> ... <c>equicktandem.c</c>, EMBOSS-6.6.0): the search of the ported program, the feature tables of
/// <see cref="EmbossReports"/> and the one writer <see cref="EmbossReportWriter"/>, with every <c>-rformat</c>, the report options (<c>-raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq</c>), several sequences per file
/// and the extra <c>-origfile</c> of the tandem programs. Compared with the installed programs on TestData/Alignment/emboss_report_formats_oracle.jsonl.gz (audit A5-8b, B09 F184).
/// </summary>
public static class EmbossReportPrograms
{
    private static string Int(long v) => v.ToString(CultureInfo.InvariantCulture);

    private static string Bool(bool b) => b ? "Yes" : "No";

    private sealed record Setup(EmbossReportOptions Options, IReadOnlyList<EmbossSequenceQualifiers> Inputs);

    private static Setup Prepare(string program, EmbossReportRequest r, string defaultFormat, string header, bool footer, bool emptyTables)
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
            Tags = EmbossReports.TagList(program),
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

    private static string PatternDoc(string patternName, int mismatch, string pattern)
        => string.Create(CultureInfo.InvariantCulture, $"{"Pattern_name",-12} {"Mismatch",8} {"Pattern"}\n{patternName,-12} {mismatch,8} {pattern}\n");

    // the USA of a sequence read with a range: the qualifiers as typed
    private static string Usa(EmbossReportRequest r, string name, EmbossSequenceQualifiers q)
    {
        string usa = "fasta::" + r.FastaFile + ":" + name;
        if (q.Begin != 0 || q.End != 0 || q.Reverse)
            usa += "[" + (q.Begin != 0 ? Int(q.Begin) : "") + ":" + (q.End != 0 ? Int(q.End) : "") + (q.Reverse ? ":r" : "") + "]";
        return usa;
    }

    // ajSeqGetBegin / ajSeqGetEnd: unset = whole, negative = from the end
    private static (int Begin, int End) Range(EmbossSequenceView v)
    {
        int len = v.Full.Length;
        return (Math.Max(1, Position(v.Begin, 1, len)), Math.Min(len, Position(v.End, len, len)));
    }

    private static int Position(int given, int unset, int len)
    {
        if (given == 0)
            return unset;
        return given < 0 ? len + given + 1 : given;
    }

    private static EmbossReportSequence ReportSequence(EmbossReportRequest r, EmbossFastaRecord rec, EmbossSequenceView view, EmbossSequenceQualifiers q, bool nucleotide)
    {
        var (b, e) = Range(view);
        return new EmbossReportSequence(rec.Name, Usa(r, rec.Name, q), rec.Accession, rec.Description, view.Full, b, e, 0, view.Reversed, nucleotide);
    }

    private static string Finish(Setup s, List<EmbossReport> tables, List<int> seqIndex, IReadOnlyList<long> lengths, int count)
    {
        return EmbossReportWriter.WriteFile(s.Options, tables, stoppedAt =>
        {
            int read = stoppedAt < 0 ? count : seqIndex[stoppedAt] + 1;
            return (read, lengths.Take(read).Sum());
        });
    }

    /// <summary>The report file of <c>fuzznuc</c> (<c>fuzznuc.c</c>): header text = pattern documentation and <c>Complement</c>; default format <c>seqtable</c>.</summary>
    /// <param name="r">The request.</param>
    /// <param name="pattern">The pattern (<c>-pattern</c>).</param>
    /// <param name="mismatch">Mismatches (<c>-pmismatch</c>).</param>
    /// <param name="complement">Also search the complement strand (<c>-complement</c>).</param>
    public static EmbossProgramReport Fuzznuc(EmbossReportRequest r, string pattern, int mismatch = 0, bool complement = false, IReadOnlyList<PatternHit>? hits = null)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        var s = Prepare("fuzznuc", r, "seqtable", PatternDoc(r.PatternName, mismatch, pattern) + "\nComplement: " + Bool(complement) + "\n", true, false);
        var records = EmbossFastaReader.Read(r.Fasta);
        var tables = new List<EmbossReport>();
        var index = new List<int>();
        var lengths = new List<long>();
        for (int i = 0; i < records.Count; i++)
        {
            var rec = records[i];
            var view = EmbossInputSequence.Read(rec.Text, s.Inputs[0], EmbossSequenceType.Dna);
            lengths.Add(view.Residues.Length);
            var found = i == 0 && hits is not null ? hits : GenomicAnalyzer.FindPatternHits(rec.Text, s.Inputs[0], pattern, mismatch, complement, r.CancellationToken);
            var table = EmbossReports.Fuzznuc(ReportSequence(r, rec, view, s.Inputs[0], view.IsNucleotide), pattern, found, r.PatternName, mismatch);
            if (table.Features.Count > 0)
            {
                tables.Add(table);
                index.Add(i);
            }
        }

        return new EmbossProgramReport(Finish(s, tables, index, lengths, records.Count), null);
    }

    /// <summary>The report file of <c>fuzzpro</c>; default format <c>seqtable</c>.</summary>
    public static EmbossProgramReport Fuzzpro(EmbossReportRequest r, string pattern, IReadOnlyList<ProteinPatternHit>? hits = null)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        var s = Prepare("fuzzpro", r, "seqtable", PatternDoc(r.PatternName, 0, pattern), true, false);
        var records = EmbossFastaReader.Read(r.Fasta);
        var tables = new List<EmbossReport>();
        var index = new List<int>();
        var lengths = new List<long>();
        for (int i = 0; i < records.Count; i++)
        {
            var rec = records[i];
            var view = EmbossInputSequence.Read(rec.Text, s.Inputs[0], EmbossSequenceType.Protein);
            lengths.Add(view.Residues.Length);
            var found = i == 0 && hits is not null ? hits : GenomicAnalyzer.FindProteinPatternHits(rec.Text, s.Inputs[0], pattern, r.CancellationToken);
            var table = EmbossReports.Fuzzpro(ReportSequence(r, rec, view, s.Inputs[0], false), pattern, found, r.PatternName, i == 0);
            if (table.Features.Count > 0)
            {
                tables.Add(table);
                index.Add(i);
            }
        }

        return new EmbossProgramReport(Finish(s, tables, index, lengths, records.Count), null);
    }

    /// <summary>The report file of <c>fuzztran</c>; default format <c>table</c>; the header names the genetic code table and the frames as typed.</summary>
    /// <param name="r">The request.</param>
    /// <param name="pattern">The protein pattern.</param>
    /// <param name="frame">The frames (<c>-frame</c>: 1 2 3 F -1 -2 -3 R 6).</param>
    /// <param name="table">The genetic code (<c>-table</c>).</param>
    public static EmbossProgramReport Fuzztran(EmbossReportRequest r, string pattern, string frame = "1", int table = 0, IReadOnlyList<TranslatedPatternHit>? hits = null)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(frame);
        var s = Prepare("fuzztran", r, "table", PatternDoc(r.PatternName, 0, pattern) + "TransTable: " + Int(table) + "\nFrames: " + frame + "\n", true, false);
        var records = EmbossFastaReader.Read(r.Fasta);
        var tables = new List<EmbossReport>();
        var index = new List<int>();
        var lengths = new List<long>();
        for (int i = 0; i < records.Count; i++)
        {
            var rec = records[i];
            var view = EmbossInputSequence.Read(rec.Text, s.Inputs[0], EmbossSequenceType.Dna);
            lengths.Add(view.Residues.Length);
            var found = i == 0 && hits is not null ? hits : GenomicAnalyzer.FindTranslatedPatternHits(rec.Text, s.Inputs[0], pattern, frame, table, r.CancellationToken);
            var t = EmbossReports.Fuzztran(ReportSequence(r, rec, view, s.Inputs[0], view.IsNucleotide), pattern, found, r.PatternName);
            if (t.Features.Count > 0)
            {
                tables.Add(t);
                index.Add(i);
            }
        }

        return new EmbossProgramReport(Finish(s, tables, index, lengths, records.Count), null);
    }

    private static EmbossFastaRecord FirstRecord(EmbossReportRequest r)
    {
        var records = EmbossFastaReader.Read(r.Fasta);
        if (records.Count == 0)
            throw new ArgumentException("The input has no FASTA record.", nameof(r));
        return records[0];
    }

    /// <summary>The report file of <c>etandem</c> (one sequence; also the <c>-origfile</c> text); default format <c>table</c>.</summary>
    public static EmbossProgramReport Etandem(EmbossReportRequest r, int minRepeat, int maxRepeat, int threshold = 20, bool mismatch = false, bool uniform = false, IReadOnlyList<EmbossTandemRepeat>? hits = null)
    {
        var s = Prepare("etandem", r, "table", $"Threshold: {Int(threshold)}\nMinrepeat: {Int(minRepeat)}\nMaxrepeat: {Int(maxRepeat)}\nMismatch: {Bool(mismatch)}\nUniform: {Bool(uniform)}\n", false, true);
        var rec = FirstRecord(r);
        var view = EmbossInputSequence.Read(rec.Text, s.Inputs[0], EmbossSequenceType.Dna);
        var found = hits ?? GenomicAnalyzer.FindTandemRepeatsEtandem(rec.Text, s.Inputs[0], minRepeat, maxRepeat, threshold, mismatch, uniform, r.CancellationToken);
        var table = EmbossReports.Etandem(ReportSequence(r, rec, view, s.Inputs[0], view.IsNucleotide), found);
        var orig = new StringBuilder();
        foreach (var h in table.Features)
        {
            // etandem_basicReport: "%6d %10d %10d %2d %3d %5.1f " and the consensus
            string size = h.Tags[1].Value[("*rpt_size ".Length)..], count = h.Tags[2].Value[("*rpt_count ".Length)..], ident = h.Tags[3].Value[("*identity ".Length)..], cons = h.Tags[4].Value[("*consensus ".Length)..];
            orig.Append(Int((long)h.Score).PadLeft(6)).Append(' ').Append(Int(h.Start).PadLeft(10)).Append(' ').Append(Int(h.End).PadLeft(10)).Append(' ').Append(size.PadLeft(2)).Append(' ').Append(count.PadLeft(3)).Append(' ').Append(ident.PadLeft(5)).Append(' ').Append(cons).Append('\n');
        }

        return Written(s, table, orig.ToString());
    }

    // the report of one table; the original file of the tandem programs is written before the report is refused
    private static EmbossProgramReport Written(Setup s, EmbossReport table, string orig)
    {
        try
        {
            return new EmbossProgramReport(EmbossReportWriter.WriteFile(s.Options, [table]), orig);
        }
        catch (EmbossReportException e)
        {
            throw new EmbossReportException(e.Message) { OriginalFile = orig };
        }
    }

    /// <summary>The report file of <c>equicktandem</c> (one sequence; also the <c>-origfile</c> text); default format <c>table</c>.</summary>
    public static EmbossProgramReport Equicktandem(EmbossReportRequest r, int maxRepeat, int threshold = 20, IReadOnlyList<EmbossTandemRepeat>? hits = null)
    {
        var s = Prepare("equicktandem", r, "table", $"Threshold: {Int(threshold)}\nMaxrepeat: {Int(maxRepeat)}\n", false, true);
        var rec = FirstRecord(r);
        var view = EmbossInputSequence.Read(rec.Text, s.Inputs[0], EmbossSequenceType.Dna);
        var found = hits ?? GenomicAnalyzer.FindTandemRepeatsEquicktandem(rec.Text, s.Inputs[0], maxRepeat, threshold, r.CancellationToken);
        var table = EmbossReports.Equicktandem(ReportSequence(r, rec, view, s.Inputs[0], view.IsNucleotide), found);
        var orig = new StringBuilder();
        foreach (var h in table.Features)
            orig.Append(Int((long)h.Score).PadLeft(6)).Append(' ').Append(Int(h.Start).PadLeft(10)).Append(' ').Append(Int(h.End).PadLeft(10)).Append(' ').Append(h.Tags[1].Value[("*rpt_size ".Length)..].PadLeft(2)).Append(' ').Append(h.Tags[2].Value[("*rpt_count ".Length)..].PadLeft(3)).Append('\n');
        return Written(s, table, orig.ToString());
    }
}
