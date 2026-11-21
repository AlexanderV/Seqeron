using System.ComponentModel;
using global::Seqeron.Genomics.Analysis;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Analysis.Tools;

// The report texts of the EMBOSS pattern / repeat programs (B09 A5-8 / A5-8b): every tool with a reportFormat delegates to EmbossReportPrograms (one writer, every -rformat, the sequence and report qualifiers).
public partial class AnalysisTools
{
    private static string Qualifier(string name, object value) => "-" + name + " " + string.Format(System.Globalization.CultureInfo.InvariantCulture, "{0}", value);

    // the qualifiers typed after the files (a name with its value is one entry of the report header)
    private static IEnumerable<string> HeaderEntries(string? text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return [];
        var tokens = text.Split([' ', '\t', '\r', '\n'], StringSplitOptions.RemoveEmptyEntries);
        var entries = new List<string>();
        int i = 0;
        while (i < tokens.Length)
        {
            string t = tokens[i];
            if (t.StartsWith('-') && i + 1 < tokens.Length && (!tokens[i + 1].StartsWith('-') || int.TryParse(tokens[i + 1], out _)))
            {
                entries.Add(t + " " + tokens[i + 1]);
                i += 2;
            }
            else
            {
                entries.Add(t);
                i++;
            }
        }

        return entries;
    }

    /// <summary>The report file text of a program run on one sequence (name + residues become a one-record FASTA file), in <paramref name="reportFormat"/> or with the report qualifiers of the text (<c>-rformat table -raccshow</c>).</summary>
    private static string EmbossReportText(string program, string sequenceName, string sequence, string? sequenceQualifiers, string reportFormat, IEnumerable<string> programEntries, Func<EmbossReportRequest, EmbossProgramReport> run, out string? original, CancellationToken cancellationToken = default)
    {
        string format = (reportFormat ?? "").Trim();
        bool qualifierText = format.StartsWith('-');
        var commandLine = new List<string> { global::Seqeron.Genomics.Analysis.EmbossCommandLineEcho.Parameter("sequence", "sequence.fasta") };
        commandLine.AddRange(programEntries);
        commandLine.AddRange(HeaderEntries(qualifierText ? format : "-rformat " + format));
        commandLine.AddRange(HeaderEntries(sequenceQualifiers));
        var request = new EmbossReportRequest(">" + sequenceName + "\n" + sequence + "\n", sequenceQualifiers, qualifierText ? format : "-rformat " + format, commandLine, DateTime.Now, "sequence.fasta", "stdout", CancellationToken: cancellationToken);
        try
        {
            var result = run(request);
            original = result.OriginalFile;
            return result.Report;
        }
        catch (global::Seqeron.Genomics.Alignment.EmbossReportException e)
        {
            throw new ArgumentException(e.Message, nameof(reportFormat));
        }
    }

    [McpServerTool(Name = "emboss_report_file", Title = "Genome — EMBOSS fuzznuc / fuzzpro / fuzztran / etandem / equicktandem report file", ReadOnly = true)]
    [Description("The whole report FILE of the EMBOSS 6.6.0 programs fuzznuc, fuzzpro, fuzztran, etandem and equicktandem on a FASTA text of one or several sequences (etandem / equicktandem read the first), byte for byte as the installed program writes it: every -rformat (simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq; diffseq is refused: it needs 7 report tags; embl / genbank are nucleotide only, pir / swiss protein only), the header with the command line and the pattern / parameter text, one block per sequence, the totals footer, the report options -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq, the sequence qualifiers -sbegin -send -sreverse and the extra -origfile text of etandem / equicktandem (originalFile). reportQualifiers is the text after the files, e.g. '-rformat table -raccshow -rmaxseq 2'; the default format is seqtable (fuzznuc, fuzzpro) or table (fuzztran, etandem, equicktandem). Delegates to EmbossReportPrograms and the one EmbossReportWriter.")]
    public static EmbossReportFileResult EmbossReportFile(
        [Description("fuzznuc | fuzzpro | fuzztran | etandem | equicktandem.")] string program,
        [Description("The input file text: FASTA records (name = first word of the header, an accession-like second word is the accession, the rest the description).")] string fasta,
        [Description("The pattern (fuzznuc / fuzzpro / fuzztran -pattern); give this or patternFile.")] string? pattern = null,
        [Description("fuzznuc -pmismatch.")] int maxMismatches = 0,
        [Description("The TEXT of a pattern file (-pattern @file; A5-8c): one pattern per line, or FASTA style (>name, optionally '<mismatch=n>' after the name, then the pattern); every pattern is a row of the header table and has its own features.")] string? patternFile = null,
        [Description("-pname: the pattern name (default 'pattern'; the lines of a simple file are numbered after it).")] string? patternName = null,
        [Description("-pformat: simple (one pattern per line) or fasta; default: decided by the first line of the file.")] string? patternFormat = null,
        [Description("fuzznuc -complement.")] bool complement = false,
        [Description("fuzztran -frame: 1 2 3 -1 -2 -3 F R 6 (default 1).")] string frame = "1",
        [Description("fuzztran -table (0 = standard).")] int table = 0,
        [Description("etandem -minrepeat (default 10).")] int minRepeat = 10,
        [Description("etandem / equicktandem -maxrepeat (default: etandem = minRepeat, equicktandem 600).")] int? maxRepeat = null,
        [Description("etandem / equicktandem -threshold (default 20).")] int threshold = 20,
        [Description("etandem -mismatch.")] bool mismatch = false,
        [Description("etandem -uniform.")] bool uniform = false,
        [Description("Sequence qualifiers as typed after the files: -sbegin 3 -send -2 -sreverse -slower -supper.")] string? sequenceQualifiers = null,
        [Description("Report qualifiers as typed after the files: -rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2 (negated booleans -norscoreshow).")] string? reportQualifiers = null,
        [Description("The command-line entries printed in the report header, one per line (e.g. '[-sequence] s.fa', '-pattern ACGT'); default: synthesized from the parameters.")] string[]? commandLine = null,
        [Description("The run date of the header (ISO 8601); default now.")] string? runDate = null,
        [Description("The input file name: sequence USAs read fasta::NAME:id (default sequence.fasta).")] string fastaFile = "sequence.fasta",
        [Description("The report file name printed in the header (default stdout).")] string reportFile = "stdout",
        CancellationToken cancellationToken = default)
    {
        string kind = (program ?? "").Trim().ToLowerInvariant();
        ArgumentException.ThrowIfNullOrEmpty(fasta);
        DateTime date = DateTime.Now;
        if (runDate is not null && !DateTime.TryParse(runDate, System.Globalization.CultureInfo.InvariantCulture, System.Globalization.DateTimeStyles.AdjustToUniversal | System.Globalization.DateTimeStyles.AssumeUniversal, out date))
            throw new ArgumentException("runDate is not a date", nameof(runDate));
        var entries = new List<string> { EmbossCommandLineEcho.Parameter("sequence", fastaFile) };
        switch (kind)
        {
            case "fuzznuc":
                entries.AddRange(PatternEntries(pattern, patternFile, patternName, patternFormat));
                entries.Add(Qualifier("pmismatch", maxMismatches));
                if (complement)
                    entries.Add("-complement");
                break;
            case "fuzzpro":
                entries.AddRange(PatternEntries(pattern, patternFile, patternName, patternFormat));
                break;
            case "fuzztran":
                entries.AddRange([.. PatternEntries(pattern, patternFile, patternName, patternFormat), Qualifier("frame", frame), Qualifier("table", table)]);
                break;
            case "etandem":
            case "equicktandem":
                entries.AddRange(TandemEntries(kind, minRepeat, maxRepeat ?? (kind == "equicktandem" ? 600 : minRepeat), threshold, mismatch, uniform));
                break;
            default:
                throw new ArgumentException("program must be fuzznuc, fuzzpro, fuzztran, etandem or equicktandem", nameof(program));
        }

        entries.AddRange(HeaderEntries(reportQualifiers));
        entries.AddRange(HeaderEntries(sequenceQualifiers));
        var request = new EmbossReportRequest(fasta, sequenceQualifiers, reportQualifiers, commandLine ?? [.. entries], date, fastaFile, reportFile, CancellationToken: cancellationToken);
        try
        {
            IReadOnlyList<EmbossProgramPattern> patterns = [];
            if (kind is "fuzznuc" or "fuzzpro" or "fuzztran")
            {
                string spec = patternFile is null ? pattern! : "@" + PatternFileName;
                patterns = EmbossPatternList.Read(spec, patternName, patternFormat, kind == "fuzznuc" ? maxMismatches : 0, _ => patternFile ?? "");
            }
            var result = kind switch
            {
                "fuzznuc" => EmbossReportPrograms.Fuzznuc(request, patterns, complement),
                "fuzzpro" => EmbossReportPrograms.Fuzzpro(request, patterns),
                "fuzztran" => EmbossReportPrograms.Fuzztran(request, patterns, frame, table),
                "etandem" => EmbossReportPrograms.Etandem(request, minRepeat, maxRepeat ?? minRepeat, threshold, mismatch, uniform),
                _ => EmbossReportPrograms.Equicktandem(request, maxRepeat ?? 600, threshold),
            };
            return new EmbossReportFileResult(result.Report, result.OriginalFile);
        }
        catch (global::Seqeron.Genomics.Alignment.EmbossReportException e)
        {
            throw new ArgumentException(e.Message, nameof(reportQualifiers));
        }
    }

    private const string PatternFileName = "patterns.txt";

    // -pattern (or -pattern @file) with its associated qualifiers as typed
    private static List<string> PatternEntries(string? pattern, string? patternFile, string? patternName, string? patternFormat)
    {
        if (string.IsNullOrWhiteSpace(patternFile) == string.IsNullOrWhiteSpace(pattern))
            throw new ArgumentException("Give exactly one of pattern and patternFile", nameof(pattern));
        var entries = new List<string> { EmbossCommandLineEcho.Qualifier("pattern", string.IsNullOrWhiteSpace(patternFile) ? pattern! : "@" + PatternFileName) };
        if (!string.IsNullOrEmpty(patternName))
            entries.Add(EmbossCommandLineEcho.Qualifier("pname", patternName));
        if (!string.IsNullOrEmpty(patternFormat))
            entries.Add(EmbossCommandLineEcho.Qualifier("pformat", patternFormat));
        return entries;
    }
}
