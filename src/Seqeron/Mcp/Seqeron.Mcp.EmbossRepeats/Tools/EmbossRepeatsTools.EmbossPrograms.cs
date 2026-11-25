// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossRepeats, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.EmbossRepeats.Tools;

public partial class EmbossRepeatsTools
{
    [McpServerTool(Name = "find_tandem_repeats_emboss", Title = "Genome — EMBOSS etandem / equicktandem", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 tandem-repeat programs ported statement by statement (diagonal-score models, not the exact-run finder find_tandem_repeats). program etandem (default): for every unit size minRepeat..maxRepeat a running per-column consensus score (+1 when the new base is the column's most frequent base, -1 otherwise; N counts as a fifth base, or as a plain -1 with mismatch) extends candidates in steps of the unit size; candidates need a best score >= threshold and, unless uniform, a consensus with two different bases; overlapping candidates lose to the highest score (ties: smaller unit). Reports score, size, count, identity = 100 * (score + n * (copies + 1)) / (2 * n * copies) and the lower-case consensus (phased as printed). program equicktandem: the quick pre-screen - the sequence compared with itself shifted by every gap 1..maxRepeat with a +1 / -1 diagonal score; reported stretches are masked for larger gaps; no consensus / identity. Hits are in the program's report order (etandem: highest score first). Non-ACGT symbols count as N (etandem) / one 'other' class (equicktandem). position is 0-based (EMBOSS Start - 1). Verified against etandem and equicktandem 6.6.0 on 1700 + 1300 random and planted sequences. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program). etandem with threshold = -2147483648 (A4-4): the program wraps threshold - 1 to 2147483647, reports its empty root and dies with SIGFPE (\"Floating point exception\") as soon as one repeat qualifies - refused here with InvalidOperationException; a sequence without a qualifying repeat returns no repeats as the program prints nothing; -2147483647 and above are exact.")]
    public static FindEmbossTandemRepeatsResult FindTandemRepeatsEmboss(
        [Description("Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive).")] string sequence,
        [Description("etandem (default) | equicktandem (case-insensitive).")] string program = "etandem",
        [Description("etandem -minrepeat: smallest unit size (>= 2, default 10 as in the program). Ignored by equicktandem.")] int minRepeat = 10,
        [Description("Largest unit size (-maxrepeat). Default: etandem = minRepeat (as in the program), equicktandem = 600.")] int? maxRepeat = null,
        [Description("Minimum score (-threshold, default 20; equicktandem needs >= 1).")] int threshold = 20,
        [Description("etandem -mismatch: count N as a mismatch instead of a fifth base (default false).")] bool mismatch = false,
        [Description("etandem -uniform: also report single-base consensus repeats (default false).")] bool uniform = false,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "etandem ignores -sbegin / -send (the whole sequence is searched) and numbers a -sreverse sequence on the reverse strand; equicktandem searches the range and numbers from its first position." + " " + EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null,
        [Description("EMBOSS report format (-rformat): any of simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq (diffseq needs 7 report tags: refused as the program does; embl / genbank are nucleotide only, pir / swiss protein only), or the report qualifiers as typed after the files, e.g. '-rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2'. The report FILE of the program (header with the command line, one block per sequence, footer; #Rundate / #!Date is now; sequence name from sequenceName) is returned in the result field report. Omit for no report text.")] string? reportFormat = null,
        [Description("Sequence name printed in the report (default sequence).")] string sequenceName = "sequence")
    {
        if (reportFormat is not null)
        {
            string kindName = (program ?? "etandem").Trim().ToLowerInvariant();
            if (kindName is not ("etandem" or "equicktandem"))
                throw new ArgumentException("program must be etandem or equicktandem", nameof(program));
            int max = maxRepeat ?? (kindName == "equicktandem" ? 600 : minRepeat);
            string[] cmd = TandemEntries(kindName, minRepeat, max, threshold, mismatch, uniform);
            if (sequenceQualifiers is null)
            {
                // one search: the report is built from the hits of the search the tool result is built from (A6-10)
                var (which, found) = RunEmbossTandem(sequence, program, minRepeat, maxRepeat, threshold, mismatch, uniform, cancellationToken);
                string text = EmbossReportText(which, sequenceName, sequence, null, reportFormat, cmd, r => RunTandemReport(which, r, minRepeat, max, threshold, mismatch, uniform, found), out var originalFile, cancellationToken);
                return TandemResult(which, found) with { Report = text, OriginalFile = originalFile };
            }

            var qualifiedResult = FindTandemRepeatsEmboss(sequence, program ?? "etandem", minRepeat, maxRepeat, threshold, mismatch, uniform, cancellationToken, sequenceQualifiers, null, sequenceName);
            return qualifiedResult with { Report = EmbossReportText(kindName, sequenceName, sequence, sequenceQualifiers, reportFormat, cmd, r => RunTandemReport(kindName, r, minRepeat, max, threshold, mismatch, uniform, null), out var original, cancellationToken), OriginalFile = original };
        }

        if (sequenceQualifiers is not null)
        {
            string kind = (program ?? "etandem").Trim().ToLowerInvariant();
            if (kind is not ("etandem" or "equicktandem"))
                throw new ArgumentException("program must be etandem or equicktandem", nameof(program));
            var q = EmbossSequenceQualifiers.ParseInputs(kind, sequenceQualifiers)[0];
            IReadOnlyList<EmbossTandemRepeat> qualified = kind switch
            {
                "etandem" => EmbossRepeatPrograms.FindTandemRepeatsEtandem(sequence, q, minRepeat, maxRepeat ?? minRepeat, threshold, mismatch, uniform, cancellationToken),
                "equicktandem" => EmbossRepeatPrograms.FindTandemRepeatsEquicktandem(sequence, q, maxRepeat ?? 600, threshold, cancellationToken),
                _ => throw new ArgumentException("program must be etandem or equicktandem", nameof(program)),
            };
            return new FindEmbossTandemRepeatsResult(kind, [.. qualified.Select(h => new EmbossTandemItem(h.Position, h.Length, h.Score, h.RepeatSize, h.RepeatCount, kind == "etandem" ? Math.Round((double)h.Identity, 1) : null, h.Consensus))]);
        }

        var (plainWhich, plainHits) = RunEmbossTandem(sequence, program, minRepeat, maxRepeat, threshold, mismatch, uniform, cancellationToken);
        return TandemResult(plainWhich, plainHits);
    }
    [McpServerTool(Name = "seqmatchall", Title = "Genome — EMBOSS seqmatchall (all-against-all word matches)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 seqmatchall: the all-against-all comparison of a set of sequences by exact word matches. For every pair i < j with length(i) > wordSize and length(j) >= wordSize every maximal exact match of at least wordSize symbols is reported (forward strand only; overlapping matches and repeated positions are all listed). Nucleotide or protein sequences (any letters, compared literally and case-insensitively); the program's word table leaves out words around its unknown symbol (N for a nucleotide sequence, X otherwise) with a scan that can keep a word that holds an N after a skipped one - reproduced exactly. Positions are 0-based (the program prints 1-based starts). Matches are listed per pair, by position (the program lists them in hash order: compare as sets). Verified against seqmatchall 6.6.0 on 330 random sets (2-5 sequences, nucleotide / N / lower case / protein, word sizes 2-12, 700+ pairs); on A/C/G/T input it equals find_common_regions mode MaximalExactMatches (MUMmer -maxmatch). Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program). More than 4 194 304 maximal word matches are refused with an ArgumentException (SequenceAligner.WordMatchMaxMatches).")]
    public static SeqMatchAllResult SeqMatchAll(
        [Description("The sequences (at least one; nucleotide or protein letters).")] string[] sequences,
        [Description("Word size (seqmatchall -wordsize, default 4, >= 1).")] int wordSize = 4,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "seqmatchall ignores -sbegin / -send; -sreverse (taken by the whole set when its first sequence is nucleotide) only labels the sequences reversed: the words are matched on the sequences as given and numbered from the end. The result then lists the program's start / end positions (reversed: start > end) in 'reportedMatches'." + " " + EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required", nameof(sequences));
        if (sequences.Any(string.IsNullOrEmpty))
            throw new ArgumentException("Sequences cannot be null or empty", nameof(sequences));
        if (sequenceQualifiers is not null)
        {
            var reported = EmbossRepeatPrograms.SeqMatchAll(sequences, EmbossSequenceQualifiers.ParseInputs("seqmatchall", sequenceQualifiers)[0], wordSize, cancellationToken)
                .Select(m => new ReportedWordMatchItem(m.Sequence1, m.Sequence2, m.Start1, m.End1, m.Start2, m.End2, m.Length, m.Reversed)).ToArray();
            return new SeqMatchAllResult([], reported);
        }
        var items = EmbossRepeatPrograms.SeqMatchAll(sequences, wordSize, cancellationToken)
            .Select(m => new WordMatchItem(m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length))
            .ToArray();
        return new SeqMatchAllResult(items);
    }
    [McpServerTool(Name = "wordmatch", Title = "Alignment — EMBOSS wordmatch (exact matches of a pattern set, Rabin-Karp)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 wordmatch: every exact match of at least wordSize symbols between a set of pattern sequences (a FASTA text) and large sequences (a FASTA text), found by Rabin-Karp (260 nucleotide / protein runs: alignments in match, pair, simple, fasta, srspair, markx3, score and srs format, the statistics log and both GFF3 feature files byte for byte; the feature files in every output feature format of the program (-offormat4 / -offormat5 via featureQualifiers, the format:: prefix of the file name): 500 further runs equal; the sequence qualifiers of both inputs -sbegin1 / -send1 -sreverse1 -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1 and the same with 2 via featureQualifiers: 330 further runs equal, the first input is a set whose ranges cut the residues, the second a stream whose -sbegin2 / -send2 the program ignores, a reversed sequence is reported on the reverse strand). The table holds every word of the pattern sequences (a pattern sequence is identified by its NAME; words around N, or X for a protein, are skipped as in seqmatchall); a window found in the table is verified and the match is extended to the right over the pattern sequence; a pattern sequence is not examined again before the end of its previous match in that scanned sequence. Returns the matches (0-based starts, per scanned sequence and pattern sequence in the program's order), the number of distinct words, the number and total length of the matches and the files of the program: the alignments (-outfile, alignFormat), the log (-logfile) and the feature files (-aoutfeat / -boutfeat). Distinct words of five or more symbols with equal hashes are not distinguished (the program picks an arbitrary one).")]
    public static WordmatchToolResult Wordmatch(
        [Description("The pattern sequences as FASTA text (-asequence); every record is a pattern sequence.")] string patternFasta,
        [Description("The sequences to scan as FASTA text (-bsequence); every record is scanned.")] string sequenceFasta,
        [Description("Word size (wordmatch -wordsize, default 4, at least 2).")] int wordSize = 4,
        [Description("Dump the matches as alignments (-dumpalign, default true).")] bool dumpAlign = true,
        [Description("Dump the matches as GFF3 feature files (-dumpfeat, default true).")] bool dumpFeatures = true,
        [Description("The alignment format (-aformat, default match).")] string alignFormat = "match",
        CancellationToken cancellationToken = default,
        [Description("The feature output qualifiers as typed on the command line: -offormat4 / -offormat5 (the format of -aoutfeat / -boutfeat: gff, gff2, gff3, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug), -ofopenfile4 / -ofopenfile5, -ofname, -ofextension, -ofdirectory, -ofsingle. Default: GFF3 files. The sequence qualifiers of the inputs can be given here too: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1 and -sbegin2 ... for the scanned sequences.")] string? featureQualifiers = null,
        [Description("The value of -aoutfeat (the file name, optionally with a 'format::' prefix; default unknown.gff).")] string? featureFile1 = null,
        [Description("The value of -boutfeat (default unknown.gff).")] string? featureFile2 = null)
    {
        if (string.IsNullOrWhiteSpace(patternFasta))
            throw new ArgumentException("The pattern set cannot be empty", nameof(patternFasta));
        if (string.IsNullOrWhiteSpace(sequenceFasta))
            throw new ArgumentException("The sequences cannot be empty", nameof(sequenceFasta));
        var set = EmbossFastaReader.Read(patternFasta);
        var scanned = EmbossFastaReader.Read(sequenceFasta);
        if (set.Count == 0 || scanned.Count == 0)
            throw new ArgumentException("An input has no FASTA record.", nameof(patternFasta));
        var files = EmbossWordmatchPrograms.Run(new EmbossWordmatchRequest(patternFasta, sequenceFasta, wordSize, dumpAlign, dumpFeatures, [EmbossCommandLineEcho.Parameter("asequence", "a.fasta"), EmbossCommandLineEcho.Parameter("bsequence", "b.fasta"), Qualifier("wordsize", wordSize)],
            DateTime.Now, "a.fasta", "b.fasta", "stdout", alignFormat, cancellationToken, featureQualifiers, featureFile1, featureFile2));
        var result = files.Result!;
        var matches = new List<WordmatchHitItem>();
        for (int q = 0; q < result.Matches.Count; q++)
            foreach (var perPattern in result.Matches[q])
                matches.AddRange(perPattern.Select(h => new WordmatchHitItem(set[h.Pattern].Name, h.PatternStart, scanned[q].Name, h.QueryStart, h.Length)));
        return new WordmatchToolResult([.. matches], result.PatternCount, result.AllMatches, result.SumLength, files.Alignments, files.Log, files.FeatureFile1, files.FeatureFile2, files.FeatureFileName1, files.FeatureFileName2);
    }
}
