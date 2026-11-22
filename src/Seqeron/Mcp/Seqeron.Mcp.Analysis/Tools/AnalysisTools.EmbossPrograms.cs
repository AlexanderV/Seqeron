using System.ComponentModel;
using ModelContextProtocol.Server;
using global::Seqeron.Genomics.Analysis;

namespace Seqeron.Mcp.Analysis.Tools;

// EMBOSS program ports on the Analysis server (B09 A3-8 / A3-12 / A3-13); every tool delegates to GenomicAnalyzer.
public partial class AnalysisTools
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
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "etandem ignores -sbegin / -send (the whole sequence is searched) and numbers a -sreverse sequence on the reverse strand; equicktandem searches the range and numbers from its first position." + " " + global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null,
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
            var q = global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.ParseInputs(kind, sequenceQualifiers)[0];
            IReadOnlyList<EmbossTandemRepeat> qualified = kind switch
            {
                "etandem" => GenomicAnalyzer.FindTandemRepeatsEtandem(sequence, q, minRepeat, maxRepeat ?? minRepeat, threshold, mismatch, uniform, cancellationToken),
                "equicktandem" => GenomicAnalyzer.FindTandemRepeatsEquicktandem(sequence, q, maxRepeat ?? 600, threshold, cancellationToken),
                _ => throw new ArgumentException("program must be etandem or equicktandem", nameof(program)),
            };
            return new FindEmbossTandemRepeatsResult(kind, [.. qualified.Select(h => new EmbossTandemItem(h.Position, h.Length, h.Score, h.RepeatSize, h.RepeatCount, kind == "etandem" ? Math.Round((double)h.Identity, 1) : null, h.Consensus))]);
        }

        var (plainWhich, plainHits) = RunEmbossTandem(sequence, program, minRepeat, maxRepeat, threshold, mismatch, uniform, cancellationToken);
        return TandemResult(plainWhich, plainHits);
    }

    private static string[] TandemEntries(string kind, int minRepeat, int max, int threshold, bool mismatch, bool uniform)
    {
        var entries = new List<string>();
        if (kind != "equicktandem")
            entries.Add(Qualifier("minrepeat", minRepeat));
        entries.Add(Qualifier("maxrepeat", max));
        entries.Add(Qualifier("threshold", threshold));
        if (kind != "equicktandem" && mismatch)
            entries.Add("-mismatch");
        if (kind != "equicktandem" && uniform)
            entries.Add("-uniform");
        return [.. entries];
    }

    private static EmbossProgramReport RunTandemReport(string kind, EmbossReportRequest r, int minRepeat, int max, int threshold, bool mismatch, bool uniform, IReadOnlyList<EmbossTandemRepeat>? hits)
    {
        if (kind == "equicktandem")
            return EmbossReportPrograms.Equicktandem(r, max, threshold, hits);
        return EmbossReportPrograms.Etandem(r, minRepeat, max, threshold, mismatch, uniform, hits);
    }

    private static (string Which, IReadOnlyList<EmbossTandemRepeat> Hits) RunEmbossTandem(string sequence, string? program, int minRepeat, int? maxRepeat, int threshold, bool mismatch, bool uniform, CancellationToken cancellationToken)
    {
        var nucleotides = RequireEmbossIupacDna(sequence, nameof(sequence));
        string which = (program ?? "etandem").Trim().ToLowerInvariant();
        return which switch
        {
            "etandem" => (which, GenomicAnalyzer.FindTandemRepeatsEtandem(nucleotides, minRepeat, maxRepeat ?? minRepeat, threshold, mismatch, uniform, cancellationToken)),
            "equicktandem" => (which, GenomicAnalyzer.FindTandemRepeatsEquicktandem(nucleotides, maxRepeat ?? 600, threshold, cancellationToken)),
            _ => throw new ArgumentException("program must be etandem or equicktandem", nameof(program)),
        };
    }

    private static FindEmbossTandemRepeatsResult TandemResult(string which, IReadOnlyList<EmbossTandemRepeat> found)
    {
        bool detailed = which == "etandem";
        return new FindEmbossTandemRepeatsResult(which, [.. found.Select(h => new EmbossTandemItem(h.Position, h.Length, h.Score, h.RepeatSize, h.RepeatCount, detailed ? Math.Round((double)h.Identity, 1) : null, h.Consensus))]);
    }

    [McpServerTool(Name = "seqmatchall", Title = "Genome — EMBOSS seqmatchall (all-against-all word matches)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 seqmatchall: the all-against-all comparison of a set of sequences by exact word matches. For every pair i < j with length(i) > wordSize and length(j) >= wordSize every maximal exact match of at least wordSize symbols is reported (forward strand only; overlapping matches and repeated positions are all listed). Nucleotide or protein sequences (any letters, compared literally and case-insensitively); the program's word table leaves out words around its unknown symbol (N for a nucleotide sequence, X otherwise) with a scan that can keep a word that holds an N after a skipped one - reproduced exactly. Positions are 0-based (the program prints 1-based starts). Matches are listed per pair, by position (the program lists them in hash order: compare as sets). Verified against seqmatchall 6.6.0 on 330 random sets (2-5 sequences, nucleotide / N / lower case / protein, word sizes 2-12, 700+ pairs); on A/C/G/T input it equals find_common_regions mode MaximalExactMatches (MUMmer -maxmatch). Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program). More than 4 194 304 maximal word matches are refused with an ArgumentException (SequenceAligner.WordMatchMaxMatches).")]
    public static SeqMatchAllResult SeqMatchAll(
        [Description("The sequences (at least one; nucleotide or protein letters).")] string[] sequences,
        [Description("Word size (seqmatchall -wordsize, default 4, >= 1).")] int wordSize = 4,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "seqmatchall ignores -sbegin / -send; -sreverse (taken by the whole set when its first sequence is nucleotide) only labels the sequences reversed: the words are matched on the sequences as given and numbered from the end. The result then lists the program's start / end positions (reversed: start > end) in 'reportedMatches'." + " " + global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required", nameof(sequences));
        if (sequences.Any(string.IsNullOrEmpty))
            throw new ArgumentException("Sequences cannot be null or empty", nameof(sequences));
        if (sequenceQualifiers is not null)
        {
            var reported = GenomicAnalyzer.SeqMatchAll(sequences, global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.ParseInputs("seqmatchall", sequenceQualifiers)[0], wordSize, cancellationToken)
                .Select(m => new ReportedWordMatchItem(m.Sequence1, m.Sequence2, m.Start1, m.End1, m.Start2, m.End2, m.Length, m.Reversed)).ToArray();
            return new SeqMatchAllResult([], reported);
        }
        var items = GenomicAnalyzer.SeqMatchAll(sequences, wordSize, cancellationToken)
            .Select(m => new WordMatchItem(m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length))
            .ToArray();
        return new SeqMatchAllResult(items);
    }

    [McpServerTool(Name = "fuzzpro", Title = "Protein — EMBOSS fuzzpro pattern search", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 fuzzpro: every window of a protein sequence that matches a PROSITE-like pattern: amino-acid letters (case-insensitive; B = D or N, Z = E or Q, J = I or L; x = any residue), classes [ILV], exclusions {P} (plain letters only: B, Z, J, X inside [] or {} are rejected as in the binary), repeats x(3), length ranges x(2,4) (every matching window length is reported; (0,n) makes the element optional), '-' separators, the anchors '<' (N terminus) and '>' (C terminus), a terminating '.'; e.g. [DE](2)HS{P}x(2)Px(2,4)C. The search is exact unless mismatch > 0 (fuzzpro -pmismatch; verified on 750 random patterns with 1 - 3 mismatches). A start-anchored pattern with a range reports only its longest match (the binary's regular-expression search class); a single trailing stop '*' is not part of the sequence. Positions are 0-based; hits are ordered by position, then length (with mismatch > 0 and a range: the program's order, a window listed once per matching unrolling). Verified against fuzzpro 6.6.0 on 900 random patterns (letters, classes, exclusions, repeats, ranges, anchors; sequences with X, B, Z, U, J, O, * and lower case). For nucleotide patterns see find_pattern_hits (fuzznuc), for a translated search fuzztran. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static FuzzproResult Fuzzpro(
        [Description("Protein sequence (letters; a trailing * is dropped).")] string sequence,
        [Description("The pattern (fuzzpro -pattern).")] string pattern,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "fuzzpro searches the range and counts positions from the start of the original sequence; a protein is never reversed." + " " + global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null,
        [Description("EMBOSS report format (-rformat): any of simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq (diffseq needs 7 report tags: refused as the program does; embl / genbank are nucleotide only, pir / swiss protein only), or the report qualifiers as typed after the files, e.g. '-rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2'. The report FILE of the program (header with the command line, one block per sequence, footer; #Rundate / #!Date is now; sequence name from sequenceName) is returned in the result field report. Omit for no report text.")] string? reportFormat = null,
        [Description("Sequence name printed in the report (default sequence).")] string sequenceName = "sequence",
        [Description("Maximum substituted residues per hit (fuzzpro -pmismatch, default 0 = exact; must be less than the real pattern length for a pattern with a class or an exclusion and no range: the program stops otherwise). Embpat's search engines: Baeza-Yates-Perleberg, Tarhio-Ukkonen-Bleasby (its skip table knows the letters A-Z only) or the recursive brute force, which with a range lets a mismatch fall only on the mandatory part and lists a window once per matching unrolling; each hit then carries its mismatches (the report: Score = length - mismatches and a *mismatch note).")] int mismatch = 0)
    {
        if (reportFormat is not null)
        {
            var plainResult = Fuzzpro(sequence, pattern, cancellationToken, sequenceQualifiers, null, sequenceName, mismatch);          // one search; the report reuses its hits
            var reportHits = plainResult.Hits.Select(h => new ProteinPatternHit(h.Position, h.Length, h.Mismatches)).ToList();
            return plainResult with
            {
                Report = EmbossReportText("fuzzpro", sequenceName, sequence, sequenceQualifiers, reportFormat, [EmbossCommandLineEcho.Qualifier("pattern", pattern), .. (mismatch > 0 ? new[] { Qualifier("pmismatch", mismatch) } : [])], r => EmbossReportPrograms.Fuzzpro(r, pattern, reportHits, mismatch), out _, cancellationToken),
            };
        }

        if (sequenceQualifiers is not null)
        {
            RequireFuzzproInput(sequence, pattern);
            return new FuzzproResult(ProteinPattern.Parse(pattern).ToString(), [.. GenomicAnalyzer.FindProteinPatternHits(sequence, global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.ParseInputs("fuzzpro", sequenceQualifiers)[0], pattern, cancellationToken, mismatch).Select(h => new ProteinPatternHitItem(h.Position, h.Length, h.Mismatches))]);
        }

        var (parsed, hits) = RunFuzzpro(sequence, pattern, cancellationToken, mismatch);
        return new FuzzproResult(parsed.ToString(), [.. hits.Select(h => new ProteinPatternHitItem(h.Position, h.Length, h.Mismatches))]);
    }

    private static void RequireFuzzproInput(string sequence, string pattern)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (string.IsNullOrWhiteSpace(pattern))
            throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
    }

    private static (ProteinPattern Parsed, IReadOnlyList<ProteinPatternHit> Hits) RunFuzzpro(string sequence, string pattern, CancellationToken cancellationToken, int mismatch = 0)
    {
        RequireFuzzproInput(sequence, pattern);
        var parsed = ProteinPattern.Parse(pattern);
        return (parsed, GenomicAnalyzer.FindProteinPatternHits(sequence, parsed, cancellationToken, mismatch));
    }

    [McpServerTool(Name = "fuzztran", Title = "Genome — EMBOSS fuzztran (protein pattern search after translation)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 fuzztran: translates a nucleotide sequence in the chosen frames (frame 1, 2, 3, -1, -2, -3, F = forward three, R = reverse three, 6 = all; transeq conventions: reverse frame -k is the reverse complement read from its k-th base, dangling bases give X) with an EMBOSS genetic-code table (0 = standard, else the NCBI table number) and searches every translation with a fuzzpro protein pattern (see fuzzpro: letters, x, [ ], { }, repeats, ranges, anchors; exact). Each hit gives the frame, the 1-based protein start / end, the matched residues, the forward-strand nucleotide location (position 0-based, length = 3 per residue) and embossStart / embossEnd, the two numbers the program prints (exact for the forward frames; EMBOSS 6.6.0 prints negative numbers for the reverse frames - a defect of its source-feature mapping, reproduced; use position / length). Hits per frame in the order 1, 2, 3, -1, -2, -3, by protein position, then length. Verified against fuzztran 6.6.0 on 760 random DNA / pattern / frame / table cases (all frame selections, tables 0-11, N / IUPAC codes and lower case). Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static FuzztranResult Fuzztran(
        [Description("Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive).")] string sequence,
        [Description("The protein pattern (fuzztran -pattern).")] string pattern,
        [Description("Frames to translate (fuzztran -frame): 1 (default), 2, 3, -1, -2, -3, F, R or 6.")] string frame = "1",
        [Description("EMBOSS genetic code (fuzztran -table): 0 = standard (default), otherwise an NCBI table number.")] int table = 0,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "fuzztran cuts the (reverse-complemented) text at the begin / end AS GIVEN, numbers the hits from the first position of that text and prints negative numbers for the reverse frames (a defect of the program, reproduced)." + " " + global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null,
        [Description("EMBOSS report format (-rformat): any of simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq (diffseq needs 7 report tags: refused as the program does; embl / genbank are nucleotide only, pir / swiss protein only), or the report qualifiers as typed after the files, e.g. '-rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2'. The report FILE of the program (header with the command line, one block per sequence, footer; #Rundate / #!Date is now; sequence name from sequenceName) is returned in the result field report. Omit for no report text.")] string? reportFormat = null,
        [Description("Sequence name printed in the report (default sequence).")] string sequenceName = "sequence",
        [Description("Maximum substituted residues per hit (fuzztran -pmismatch, default 0 = exact; see fuzzpro).")] int mismatch = 0)
    {
        if (reportFormat is not null)
        {
            var plainResult = Fuzztran(sequence, pattern, frame, table, cancellationToken, sequenceQualifiers, null, sequenceName, mismatch);          // one search; the report reuses its hits
            var reportHits = plainResult.Hits.Select(h => new TranslatedPatternHit(h.Frame, h.ProteinStart, h.ProteinEnd, h.Translation, h.Position, h.Length, h.EmbossStart, h.EmbossEnd, h.Mismatches)).ToList();
            return plainResult with
            {
                Report = EmbossReportText("fuzztran", sequenceName, sequence, sequenceQualifiers, reportFormat, [EmbossCommandLineEcho.Qualifier("pattern", pattern), Qualifier("frame", frame), Qualifier("table", table), .. (mismatch > 0 ? new[] { Qualifier("pmismatch", mismatch) } : [])],
                    r => EmbossReportPrograms.Fuzztran(r, pattern, frame, table, reportHits, mismatch), out _, cancellationToken),
            };
        }

        if (sequenceQualifiers is not null)
        {
            if (string.IsNullOrWhiteSpace(pattern))
                throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
            var qualifiedParsed = ProteinPattern.Parse(pattern);
            var qualifiedHits = GenomicAnalyzer.FindTranslatedPatternHits(sequence, global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.ParseInputs("fuzztran", sequenceQualifiers)[0], pattern, frame, table, cancellationToken, mismatch)
                .Select(TranslatedItem).ToArray();
            return new FuzztranResult(qualifiedParsed.ToString(), qualifiedHits);
        }

        var (parsed, hits) = RunFuzztran(sequence, pattern, frame, table, cancellationToken, mismatch);
        return new FuzztranResult(parsed.ToString(), [.. hits.Select(TranslatedItem)]);
    }

    private static TranslatedPatternHitItem TranslatedItem(TranslatedPatternHit h) => new(h.Frame, h.ProteinStart, h.ProteinEnd, h.Translation, h.Position, h.Length, h.EmbossStart, h.EmbossEnd, h.Mismatches);

    private static (ProteinPattern Parsed, IReadOnlyList<TranslatedPatternHit> Hits) RunFuzztran(string sequence, string pattern, string frame, int table, CancellationToken cancellationToken, int mismatch = 0)
    {
        var nucleotides = RequireEmbossIupacDna(sequence, nameof(sequence));
        if (string.IsNullOrWhiteSpace(pattern))
            throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
        var parsed = ProteinPattern.Parse(pattern);
        return (parsed, GenomicAnalyzer.FindTranslatedPatternHits(nucleotides, pattern, frame, table, cancellationToken, mismatch));
    }

    [McpServerTool(Name = "showalign", Title = "Alignment — EMBOSS showalign display", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 showalign: the pretty text display of a multiple alignment, identical to the program (600 random alignments with random options equal). Every row except the reference is converted against the reference (the consensus unless refseq names a row or a 1-based number): show N (default) non-identities, I identities, S similarities (matrix score > 0), D dissimilarities, A all; '.' marks suppressed residues; similarCase (default true) encodes similarity in the case (non-identity / similarity / all as in the program); order I (input), A (names), S (decreasing similarity to the reference). Per block of width columns (default 60): the position numbers (every 10th, relative to the reference residues when a row is the reference, with -n / +n flanks), the tick ruler, the reference first when it is a row, the rows, the consensus (consensus = false hides it; bottom = false puts it first when it is the reference). The consensus is the EMBOSS cons algorithm (plurality in percent, setCase, identity in percent, gaps: a gap can be the consensus) on the alignment as given (case kept); a lower-case letter of the input is kept lower case. margin = width of the name column (0 = no names, -1 = longest name + 1); uppercase = 1-based column ranges put in upper case, e.g. '3-10,20-25'. The set type is that of the first row (nucleotide: A/C/G/T/U + IUPAC letters); matrix defaults to EDNAFULL / EBLOSUM62; shorter rows are padded with trailing gaps. html = true wraps the text in <pre> ... </pre> and highlight (EMBOSS range syntax, a colour after each pair) colours the ranges of every sequence and consensus line with <font color=...> (A3-13c2: 300 random texts equal). Rows are read as EMBOSS reads an alignment (A4-3): letters and '*' are kept, '.' '~' '-' are gaps, '?' '#' '+' are the unknown residue, every other character (digits, blanks, '_', punctuation) is dropped, and rows that become unequal are padded with trailing gaps. A display text above 134 217 728 characters (blocks x lines x (margin + width + 1): a huge margin, width 1 over a very long alignment) is refused with an ArgumentException (MultipleAlignmentDisplay.ShowAlignMaxOutputCharacters); a width beyond the alignment shows one block. sbegin / send (A4-9, 400 random ranges equal the program) cut the displayed columns: the consensus is that of the whole alignment and the numbers and ticks keep the positions of the whole alignment.")]
    public static ShowAlignResult ShowAlign(
        [Description("Sequence names (one per row; the first word is used).")] string[] names,
        [Description("The aligned sequences (DNA or protein; gaps '-', '.', '~'; case is kept).")] string[] sequences,
        [Description("Reference: '0' = the consensus (default), a 1-based row number or a sequence name.")] string refseq = "0",
        [Description("What to show: N non-identities (default), I identities, S similarities, D dissimilarities, A all.")] string show = "N",
        [Description("Output order of the rows: I input (default), A alphabetical names, S similarity to the reference.")] string order = "I",
        [Description("Scoring matrix name (EBLOSUM62, EDNAFULL, EPAM250, ...); default by sequence type.")] string? matrix = null,
        [Description("Case encodes similarity (default true).")] bool similarCase = true,
        [Description("Show the consensus line (default true).")] bool consensus = true,
        [Description("Consensus at the bottom (default true) instead of the top.")] bool bottom = true,
        [Description("Position number line (default true).")] bool number = true,
        [Description("Tick ruler line (default true).")] bool ruler = true,
        [Description("Columns per block (default 60).")] int width = 60,
        [Description("Name column width: -1 = longest name + 1 (default), 0 = no names.")] int margin = -1,
        [Description("1-based inclusive column ranges to put in upper case, e.g. '3-10,20-25' (EMBOSS range syntax; '@file' = one range per line).")] string? uppercase = null,
        [Description("Consensus plurality in percent of the total weight (default 50).")] double plurality = 50,
        [Description("Consensus: weight at or below which the letter is lower case (default half the number of sequences).")] double? setCase = null,
        [Description("Consensus: percent of identities required at a position (default 0).")] double identity = 0,
        [Description("A gap character can be the consensus (default true).")] bool gaps = true,
        [Description("-sbegin: first column shown, 1-based (0 = the first; negative = counted from the end; beyond the alignment = the last).")] int sbegin = 0,
        [Description("-send: last column shown, 1-based (0 = the last; negative = counted from the end; before the first shown column = that column).")] int send = 0,
        [Description("-html (default false): wrap the display in <pre> ... </pre> and colour the highlight ranges with <font color=...> tags.")] bool html = false,
        [Description("-highlight: column ranges to colour with html, each followed by an HTML colour (an empty colour is red), e.g. '24-45 blue 56-78 orange'; '@file' = one range per line. Used only with html.")] string? highlight = null,
        CancellationToken cancellationToken = default)
    {
        if (names is null || sequences is null || sequences.Length == 0)
            throw new ArgumentException("Names and sequences are required", nameof(sequences));
        if (names.Length != sequences.Length)
            throw new ArgumentException("One name per sequence is required", nameof(names));
        // the ACD range type of the program (ajRangeNewString): pairs of positions with any separators, e.g. '3-10, 20-25' or '1:45 67=99;765..888'; '@file' = one range per line
        var ranges = EmbossRangeSpec.Parse(uppercase).Select(r => (r.Start, r.End)).ToList();
        var highlights = EmbossRangeSpec.Parse(highlight);

        if (show is not { Length: 1 } || "AINSD".IndexOf(char.ToUpperInvariant(show[0])) < 0)
            throw new ArgumentException("show must be A, I, N, S or D", nameof(show));
        if (order is not { Length: 1 } || "IAS".IndexOf(char.ToUpperInvariant(order[0])) < 0)
            throw new ArgumentException("order must be I, A or S", nameof(order));
        Seqeron.Genomics.Infrastructure.SubstitutionMatrix? sub = null;
        if (!string.IsNullOrWhiteSpace(matrix) && !Seqeron.Genomics.Infrastructure.SubstitutionMatrix.TryGet(matrix, out sub))
            throw new ArgumentException($"Unknown matrix '{matrix}'", nameof(matrix));
        var options = new ShowAlignOptions
        {
            Matrix = sub,
            RefSeq = refseq,
            Show = (ShowAlignContent)char.ToUpperInvariant(show[0]),
            Order = (ShowAlignOrder)char.ToUpperInvariant(order[0]),
            SimilarCase = similarCase,
            Consensus = consensus,
            Bottom = bottom,
            Number = number,
            Ruler = ruler,
            Width = width,
            Margin = margin,
            Uppercase = ranges,
            Plurality = (float)plurality,
            SetCase = (float?)setCase,
            Identity = (float)identity,
            Gaps = gaps,
            Begin = sbegin,
            End = send,
            Html = html,
            Highlight = highlights,
            CancellationToken = cancellationToken,
        };
        return new ShowAlignResult(MultipleAlignmentDisplay.ShowAlign(names, sequences, options));
    }

    [McpServerTool(Name = "megamerger", Title = "Alignment — EMBOSS megamerger", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 megamerger: merges two large overlapping DNA sequences by exact word matches (600 random related / unrelated pairs equal: merged sequence and report text). All exact matches of wordSize symbols (the maximal matches of seqmatchall) are ordered longest first and reduced to a minimal non-overlapping set: the longest is kept, matches inside its dead zone are dropped, those crossing it are truncated (dropped when shorter than the word), until none is left; the kept matches are ordered by position. The merged sequence = the part before the first match (from the sequence whose overlap starts at its own start; when neither does the longer prefix, upper case), the matches (first sequence, lower case), between matches the region of the sequence whose mismatch is farther from its ends (preferFirst: always the first sequence), upper case, and the part after the last match (the longer remainder in upper case when neither overlap is at an end). The report lists every overlap start / end, mismatch region and choice as the program writes it. A second sequence shorter than the word or without a match gives an empty merge; a first sequence shorter than the word (or without a word free of N) is an error (the program stops with 'No match found'). Note: this merges by anchors; the overlap merge of two short sequences by dynamic programming is the merger tool. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program). More than 4 194 304 maximal word matches of the pair are refused with an ArgumentException (SequenceAligner.WordMatchMaxMatches).")]
    public static MegaMergerResult Megamerger(
        [Description("First DNA sequence (letters, any case).")] string sequence1,
        [Description("Second DNA sequence.")] string sequence2,
        [Description("Word size (megamerger -wordsize, default 20, at least 2).")] int wordSize = 20,
        [Description("Use the first sequence over every mismatch region (megamerger -prefer, default false).")] bool preferFirst = false,
        [Description("Name of the first sequence in the report.")] string name1 = "seq1",
        [Description("Name of the second sequence in the report.")] string name2 = "seq2",
        [Description("EMBOSS sequence-input qualifiers with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 and the same with 2 for the second sequence. Measured on EMBOSS 6.6.0: megamerger merges the TRIMMED sequences (positions of the report are those of the trimmed texts; -sreverse1 reverse-complements before the range is cut); -sprotein is refused (the sequences are DNA). Omit for the whole sequences." + " " + global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null,
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence2));
        var q = sequenceQualifiers is null ? null : global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.ParseInputs("megamerger", sequenceQualifiers);
        var r = q is null
            ? GenomicAnalyzer.MegaMerge(sequence1, sequence2, wordSize, preferFirst, name1, name2, cancellationToken)
            : GenomicAnalyzer.MegaMerge(
                new global::Seqeron.Genomics.Infrastructure.EmbossSequence(name1, sequence1, "", q[0]),
                new global::Seqeron.Genomics.Infrastructure.EmbossSequence(name2, sequence2, "", q[1]),
                wordSize, preferFirst, cancellationToken);
        return new MegaMergerResult(r.MergedSequence, r.Report,
            r.Matches.Select(m => new WordMatchItem(m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length)).ToArray());
    }

    // plain letters become a FASTA record; the text of a FASTA, GenBank or EMBL file is the file
    private static string DiffseqInput(string name, string text)
    {
        string t = text.TrimStart();
        if (t.StartsWith('>') || t.StartsWith("LOCUS ", StringComparison.Ordinal) || t.StartsWith("ID   ", StringComparison.Ordinal))
            return text;
        return ">" + name + "\n" + text + "\n";
    }

    [McpServerTool(Name = "featcopy", Title = "Alignment — EMBOSS featcopy (read and write feature tables)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 featcopy: every feature table of a file in ANY input format is read by the feature readers of the EMBOSS library and written in an output format. Input formats: embl, genbank, refseq (nucleotide, feature lines of five blanks), genpept, refseqp (protein), gff2 (versions 1 / 2, ##Type, ##sequence-region, group fields, featflags), gff3 (ID / Parent grouping, several lines of one ID, percent-encoded values, several tables), the gff3 EMBOSS wrote before 6.4, swiss, pir, dasgff (XML); with no format the library tries embl, refseq, refseqp, genpept, gff3, gff2, swiss, pir and dasgff in turn on the start of the file and keeps the first that finds features, and a 'format::' prefix of features forces one (an unknown name is an error message and the formats are tried); the value may also be a LIST: '@listfile' or 'list::listfile' (files given in the files argument; a word of a line is a file, 'format::file', 'file[begin:end]' or a list again, '#' starts a comment, a line of several words is 'Bad record', a missing file or list is an error message; the queries of a nested list are read before the pending ones; 500 runs equal the installed program, the lists that name themselves, on which the program never stops, are refused), and 'asis::' is refused for features as the program refuses it. Output formats written for a feature table: gff3 (gff, unknown input = refused as the program does), gff2, embl (em), genbank (gb, ddbj, refseq), bed (track line, score 0, blocks of the sub-features), refseqp, pir (nbrf), swiss (sw, swissprot, uniprot: FT lines wrapped at 75 with the /FTId= rule), dasgff (XML, type category from the ontology table), draw (cirdna / lindna group), debug (the full internal trace: type, source, score, strand, positions, flags, group, exon, tags, sub-features); a format the program refuses for the kind of sequence (embl for a protein table, pir for a nucleotide one) is an error message and nothing is written. The first table is cut to the range featureBegin / featureEnd (a table with a range is trimmed and its features shifted: a feature cut at the end is flagged). The reader and writer pairs equal the installed program on 3 800 input files (embl, genbank, refseq, genpept, refseqp, gff2, gff3, gff3emboss, swiss, pir, dasgff: 236 or more each, damaged and mixed files, the format:: prefix, ranges, the files on which the program never stops are refused with an error) and 1 800 output files in gff3 / gff2 / embl / genbank and their aliases and 2 750 output files in bed, refseq, refseqp, pir / nbrf, swiss / sw / swissprot / uniprot, dasgff and draw (a Swiss-Prot feature line the library's wrapper cannot break - a word as long as the line - never ends in the program and is refused with an error); the messages (Warning: bad /protein_id value ..., Error: Unknown input format ...) are returned as the program prints them.")]
    public static FeatcopyToolResult Featcopy(
        [Description("The text of the feature file (a feature table of any input format, or an EMBL / GenBank entry; several tables may follow one another).")] string features,
        [Description("The value of -features: 'format::' prefix and [begin:end] suffix are read as the program reads them; default the file itself.")] string? featuresValue = null,
        [Description("The output feature format (-offormat): gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug.")] string outputFormat = "gff3",
        [Description("-fbegin: the first table is cut to this range start (default 1; negative counts from the end).")] int featureBegin = 1,
        [Description("-fend: the first table is cut to this range end (default: its length).")] int? featureEnd = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that featuresValue may name: a list file ('@list' or 'list::list': one query per word of a line, '#' starts a comment, the queries are files with an optional 'format::' prefix and [begin:end] suffix, or lists again), the files the list names; a file that is not here does not exist. 'asis::' is refused for features as the program refuses it.")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(features))
            throw new ArgumentException("The feature file cannot be empty", nameof(features));
        cancellationToken.ThrowIfCancellationRequested();
        string name = "features.dat";
        string value = string.IsNullOrWhiteSpace(featuresValue) ? name : featuresValue;
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [name] = features };
        try
        {
            var o = EmbossFeaturePrograms.Featcopy(new EmbossFeatcopyRequest(value, all, featureBegin, featureEnd, "unknown.gff", outputFormat, null, DateTime.Now));
            return new FeatcopyToolResult(o.FileName, o.Text, o.Format, o.Tables, [.. o.Messages]);
        }
        catch (global::Seqeron.Genomics.Infrastructure.EmbossFeatureDiedException e)
        {
            throw new ArgumentException(e.Message + (e.Errors.Count > 0 ? " (" + string.Join("; ", e.Errors.Where(m => m.Length > 0)) + ")" : ""), nameof(features));
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(outputFormat));
        }
    }

    [McpServerTool(Name = "featmerge", Title = "Alignment — EMBOSS featmerge (combine two feature tables)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 featmerge: the first feature table of featuresA and the first of featuresB (each in any input format of the EMBOSS feature readers: embl, genbank, refseq, genpept, refseqp, gff2, gff3, swiss, pir, dasgff; with no format every format is tried; 'format::' prefix and [begin:end] suffix in the *Value arguments) are read, a copy of every feature of B (every tag checked again, the group numbers kept) is added to A, and A is written in an output format: gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug. The range of A (-fbegin1 / -fend1) cuts and shifts the merged table as it does in the program; the range of B is read and not used. Equals the installed program on 450 runs (every pair of input formats, ranges of both inputs, every output format, missing files, the inputs on which the program never stops refused with an error); the messages (Warning: ..., Error: Failed to open filename ...) are returned as the program prints them.")]
    public static FeatcopyToolResult Featmerge(
        [Description("The text of the first feature file (-afeatures).")] string featuresA,
        [Description("The text of the second feature file (-bfeatures).")] string featuresB,
        [Description("The value of -afeatures: 'format::' prefix and [begin:end] suffix; default the file itself.")] string? featuresAValue = null,
        [Description("The value of -bfeatures.")] string? featuresBValue = null,
        [Description("The output feature format (-offormat): gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug.")] string outputFormat = "gff3",
        [Description("-fbegin1: the range start of the first table (default 1).")] int featureBegin = 1,
        [Description("-fend1: the range end of the first table (default: its length).")] int? featureEnd = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that featuresAValue / featuresBValue may name, list files included ('@list').")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(featuresA))
            throw new ArgumentException("The feature file cannot be empty", nameof(featuresA));
        if (string.IsNullOrWhiteSpace(featuresB))
            throw new ArgumentException("The feature file cannot be empty", nameof(featuresB));
        cancellationToken.ThrowIfCancellationRequested();
        try
        {
            var o = EmbossFeaturePrograms.Featmerge(new EmbossFeatmergeRequest(string.IsNullOrWhiteSpace(featuresAValue) ? "a.dat" : featuresAValue, string.IsNullOrWhiteSpace(featuresBValue) ? "b.dat" : featuresBValue,
                new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { ["a.dat"] = featuresA, ["b.dat"] = featuresB }, featureBegin, featureEnd, 1, null, "unknown.gff", outputFormat, null, DateTime.Now));
            return new FeatcopyToolResult(o.FileName, o.Text, o.Format, o.Tables, [.. o.Messages]);
        }
        catch (global::Seqeron.Genomics.Infrastructure.EmbossFeatureDiedException e)
        {
            throw new ArgumentException(e.Message + (e.Errors.Count > 0 ? " (" + string.Join("; ", e.Errors.Where(m => m.Length > 0)) + ")" : ""), nameof(featuresA));
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(outputFormat));
        }
    }

    [McpServerTool(Name = "feattext", Title = "Alignment — EMBOSS feattext (original text of a feature table)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 feattext: every feature table of a file (any input format of the EMBOSS feature readers; with no format every format is tried; 'format::' prefix and [begin:end] suffix in featuresValue) is read with the input storing the text it reads, and the text of each entry is returned as the library stored it: the EMBL / GenBank reader and the GFF3 reader keep the lines they consumed (a GFF3 table stops at the next ##gff-version line, which is given back with the library's off-by-the-line-ending cut of the stored text), the readers of the other formats keep none (nothing is returned for their tables). Equals the installed program on 400 runs (every input format, damaged files, the format:: prefix, the runs on which the program never stops refused with an error); the messages are returned as the program prints them.")]
    public static FeattextToolResult Feattext(
        [Description("The text of the feature file (a feature table of any input format, or an EMBL / GenBank entry; several tables may follow one another).")] string features,
        [Description("The value of -features: 'format::' prefix and [begin:end] suffix are read as the program reads them; default the file itself.")] string? featuresValue = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that featuresValue may name, list files included ('@list').")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(features))
            throw new ArgumentException("The feature file cannot be empty", nameof(features));
        cancellationToken.ThrowIfCancellationRequested();
        string name = "features.dat";
        try
        {
            var o = EmbossFeaturePrograms.Feattext(new EmbossFeattextRequest(string.IsNullOrWhiteSpace(featuresValue) ? name : featuresValue, new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [name] = features }));
            return new FeattextToolResult(o.Text, o.Tables, [.. o.Messages]);
        }
        catch (global::Seqeron.Genomics.Infrastructure.EmbossFeatureDiedException e)
        {
            throw new ArgumentException(e.Message + (e.Errors.Count > 0 ? " (" + string.Join("; ", e.Errors.Where(m => m.Length > 0)) + ")" : ""), nameof(features));
        }
        catch (InvalidOperationException e)
        {
            throw new ArgumentException(e.Message, nameof(features));
        }
    }

    [McpServerTool(Name = "seqret", Title = "Alignment — EMBOSS seqret (read and write sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 seqret: the sequences of a file are read by the sequence readers of the EMBOSS library (ajseqread.c) and written in FASTA (with the GFF3 or another feature format of the features the file carries). Input formats read byte for byte as the installed program reads them (200 inputs of each, with the library's warnings and errors): fasta (ncbi, pearson, gifasta, dbid), gcg (gcg8), embl (em), swiss (sw, swissprot, uniprot; the DE description of RecName / AltName / SubName / Includes / Contains / Flags is assembled as the library does), nbrf (pir), genbank (gb, ddbj, refseq), refseqp, genpept, codata, strider, ig, igstrict, iguspto, staden, textonly (plain, asis, unknown), fastq, fastq-sanger, fastq-illumina, fastq-solexa, sam, gff2, gff3 (gff; the sequence after ##FASTA / the ##DNA lines), raw. With no format every format the program tries is tried in turn (a text no format reads is read as plain text). The features of embl, genbank, refseq, refseqp, swiss, nbrf, gff2 and gff3 files are returned with featureFormat when feature is true. Not read yet (alignment and binary formats: clustal, msf, phylip, nexus, mega, selex, stockholm, mase, fitch, treecon, hennig86, jackknifer, gde, acedb, ace, pdb, abi, bam, nibble, scf, das, ensembl): the tool refuses an input of such a format. A5-8c3g, B09 F248.")]
    public static SeqretToolResult Seqret(
        [Description("The text of the sequence file (any input format above), or a sequence as letters. A value 'asis::ACGT' is the sequence itself.")] string sequence,
        [Description("The input format (-sformat): the names above, any case, or an EDAM term; default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Read and write the features of the sequences (-feature).")] bool feature = false,
        [Description("Read one sequence and stop (-firstonly).")] bool firstOnly = false,
        [Description("The feature output format (-offormat2): gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug.")] string featureFormat = "gff3",
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -slower1 -supper1 -snucleotide1 -sprotein1 (also -sformat1 F); the range qualifiers are not applied by this tool yet.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequence may name with a 'format::' prefix.")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        string name = "sequence.dat";
        var q = global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(qualifiers) with { Format = sequenceFormat };
        string spec = sequence.TrimStart().StartsWith("asis::", StringComparison.OrdinalIgnoreCase) ? sequence.Trim() : name;
        try
        {
            var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [name] = sequence };
            var o = EmbossSeqretPrograms.Seqret(new EmbossSeqretRequest(spec, all, q, feature, firstOnly, "stdout", "fasta", "outfeat.gff", featureFormat, DateTime.Now));
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new SeqretToolResult(o.SequenceText, o.FeatureText, o.Sequences, [.. o.Messages]);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }

    [McpServerTool(Name = "maskfeat", Title = "Alignment — EMBOSS maskfeat (mask the regions of features)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 maskfeat: the regions of the features of the wanted types are masked in the sequences of a file (any input format of seqret that carries features: embl, genbank, refseq, refseqp, swiss, nbrf, gff2, gff3; the feature types are those of the vocabulary, wildcards * and ? allowed, several types separated by blanks, commas, semicolons or bars; the default repeat* masks repeat_region and repeat_unit): the residues become the mask character (default N for nucleotide and X for protein sequences) or lower case (toLower, or a blank or no mask character), and the sequences are returned in FASTA (the sequence type of the program is any: every character that is no letter and no * is dropped). Equals the installed program on 1 200 runs (every input format with features, damaged files, the runs on which the program refuses an input or dies by a fault). A5-8c3g, B09 F249.")]
    public static MaskfeatToolResult Maskfeat(
        [Description("The text of the sequence file with features (see the formats above), or a sequence as letters (nothing to mask).")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("The feature types to mask (-type): names separated by blanks, commas, semicolons or bars; * and ? wildcards; default repeat*.")] string type = "repeat*",
        [Description("Make the masked regions lower case (-tolower).")] bool toLower = false,
        [Description("The mask character (-maskchar, one character; none or a blank = lower case); default N (nucleotide) / X (protein).")] string? maskChar = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -slower1 -supper1 -snucleotide1 -sprotein1.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequence may name with a 'format::' prefix.")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        string name = "sequence.dat";
        var q = global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(qualifiers) with { Format = sequenceFormat };
        string spec = sequence.TrimStart().StartsWith("asis::", StringComparison.OrdinalIgnoreCase) ? sequence.Trim() : name;
        try
        {
            var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [name] = sequence };
            var o = EmbossMaskfeatPrograms.Maskfeat(new EmbossMaskfeatRequest(spec, all, q, type, toLower, maskChar));
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new MaskfeatToolResult(o.SequenceText, o.Sequences, [.. o.Messages]);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }

    [McpServerTool(Name = "descseq", Title = "Alignment — EMBOSS descseq (alter the name or description of a sequence)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 descseq: the name and / or the description of the FIRST sequence of a file (any input format of seqret) is replaced by name / description, or appended to the existing ones (append), and the sequence is returned in FASTA (the header is the name, the sequence version or accession, and the description). Equals the installed program on 1 000 runs (19 input formats, names and descriptions with blanks, tabs and punctuation, damaged files). A5-8c3g, B09 F250.")]
    public static DescseqToolResult Descseq(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("The new name (-name); empty = the name is kept.")] string name = "",
        [Description("The new description (-description); empty = the description is kept.")] string description = "",
        [Description("Append the name and the description to the existing ones (-append).")] bool append = false,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -slower1 -supper1 -snucleotide1 -sprotein1.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequence may name with a 'format::' prefix.")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        string file = "sequence.dat";
        var q = global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(qualifiers) with { Format = sequenceFormat };
        string spec = sequence.TrimStart().StartsWith("asis::", StringComparison.OrdinalIgnoreCase) ? sequence.Trim() : file;
        try
        {
            var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [file] = sequence };
            var o = EmbossDescseqPrograms.Descseq(new EmbossDescseqRequest(spec, all, q, name ?? "", description ?? "", append));
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new DescseqToolResult(o.SequenceText, o.Sequences, [.. o.Messages]);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }

    [McpServerTool(Name = "diffseq", Title = "Alignment — EMBOSS diffseq (differences between two sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 diffseq: the regions of difference between two nearly identical sequences (336 nucleotide / protein pairs, every -rformat: report and both GFF3 feature files byte for byte). The exact word matches of wordSize symbols (the engine and the minimal non-overlapping set of megamerger) are the regions of identity; the stretches between them are the differences (a difference with length 0 in one sequence is an insertion in the other); with globalDifferences the differences before the first and after the last match are reported too (a start or end that matches by itself is not counted). A nucleotide pair also gets the SNP count and its transitions / transversions. Returns the differences and the program's files: the report in the default format 'diffseq' (7 report tags; reportFormat takes any -rformat or the report qualifiers as typed, e.g. '-rformat table -raccshow'), and the GFF3 files -aoutfeat / -boutfeat. The sequence qualifiers of both inputs (-sbegin1 -send1 -sreverse1 -slower1 -supper1 -snucleotide1 ... and the same with 2, given in reportFormat as typed after the files, e.g. '-sbegin1 20 -send1 200 -sreverse2 -rformat table') are applied as the program applies them: the sequence is cut to the range (positions, lengths and the report numbers refer to the cut sequence, the header shows the original range), reversed sequences are reported on the reverse strand; 260 / 260 runs equal the installed program. The inputs may be the text of a GenBank or EMBL file: its features (type, location with complement / join / order, qualifiers) are read as the program reads them, follow -sbegin / -send / -sreverse, and the features that overlap a difference are listed with it (first_feature / second_feature, the 'Feature:' lines of the diffseq report); the feature files -aoutfeat / -boutfeat are written in any output feature format of the program (-offormat4 / -offormat5: gff, gff2, gff3, bed, embl, genbank, ddbj, refseq, refseqp, pir, swiss, dasgff, draw, debug and their aliases, or the 'format::' prefix of the file name; a format the program refuses for the kind of sequence - embl for protein, pir for nucleotide - gives an empty file; 1 100 runs equal the installed program); 870 runs equal the installed program (420 with the usual qualifiers, 450 with the tags the program checks with a special validator - bad values become notes, correctable ones are corrected), including its loss of the features of a GenBank entry whose SOURCE block is followed directly by FEATURES. The features of a sequence can come from a SEPARATE file: -ufo1 / -ufo2 ([format::]file), -fformat1, -fopenfile1 (and the same with 2) in reportFormat name files given in featureFiles (file name -> text), read in every input format of the EMBOSS feature readers (embl, genbank, refseq, genpept, refseqp, gff2, gff3, the gff3 of EMBOSS before 6.4, swiss, pir, dasgff; with no format every format is tried; featcopy / seqret traces of 7 800 files equal, 450 diffseq runs equal); a UFO switches the features of the sequence file off. Not read: SwissProt / GFF / GenPept files as the SEQUENCE input.")]
    public static DiffseqToolResult Diffseq(
        [Description("First sequence: nucleotide or protein letters, or the text of a FASTA, GenBank (LOCUS) or EMBL (ID) file, whose features are read.")] string sequence1,
        [Description("Second sequence (the same kind, or any of the file texts).")] string sequence2,
        [Description("Word size (diffseq -wordsize, default 10, at least 2).")] int wordSize = 10,
        [Description("Report the differences at the start and the end (diffseq -globaldifferences).")] bool globalDifferences = false,
        [Description("Name of the first sequence.")] string name1 = "seqa",
        [Description("Name of the second sequence.")] string name2 = "seqb",
        [Description("The report format (-rformat name) or the qualifiers as typed after the files (e.g. '-rformat table -raccshow'; also the feature output qualifiers '-offormat4 embl -offormat5 gff2' / '-ofopenfile4 name', and the feature input qualifiers '-ufo1 f1.gff -fformat2 swiss -fopenfile2 f2.dat'); default format diffseq, feature files GFF3.")] string? reportFormat = null,
        CancellationToken cancellationToken = default,
        [Description("The feature files the inputs name with -ufo / -fformat / -fopenfile, by file name (for example {\"f1.gff\": \"##gff-version 3\\n...\"}); a file that is not here does not exist.")] IReadOnlyDictionary<string, string>? featureFiles = null)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence2));
        string? qualifiers = null;
        if (!string.IsNullOrWhiteSpace(reportFormat))
            qualifiers = reportFormat.TrimStart().StartsWith('-') ? reportFormat : "-rformat " + reportFormat;
        var entries = new List<string> { EmbossCommandLineEcho.Parameter("asequence", "a.fasta"), EmbossCommandLineEcho.Parameter("bsequence", "b.fasta"), Qualifier("wordsize", wordSize) };
        if (globalDifferences)
            entries.Add("-globaldifferences");
        entries.AddRange(HeaderEntries(qualifiers));
        try
        {
            var files = EmbossDiffseqPrograms.Run(new EmbossDiffseqRequest(DiffseqInput(name1, sequence1), DiffseqInput(name2, sequence2), wordSize, globalDifferences, qualifiers, entries, DateTime.Now, "a.fasta", "b.fasta", "stdout", name1 + ".diffgff", name2 + ".diffgff", cancellationToken, featureFiles));
            var result = files.Result ?? new DiffseqResult([], -1, 0, 0, 0, 0, 0, 0);
            return new DiffseqToolResult([.. result.Differences.Select(d => new DiffseqDifferenceItem(d.Start1, d.End1, d.Len1, d.Start2, d.End2, d.Len2))], result.OverlapStart1, result.OverlapEnd1, result.OverlapStart2, result.OverlapEnd2,
                result.Snps, result.Transitions, result.Transversions, files.Report, files.FeatureFile1, files.FeatureFile2, files.FeatureFileName1, files.FeatureFileName2);
        }
        catch (global::Seqeron.Genomics.Alignment.EmbossReportException e)
        {
            throw new ArgumentException(e.Message, nameof(reportFormat));
        }
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
        var set = global::Seqeron.Genomics.Infrastructure.EmbossFastaReader.Read(patternFasta);
        var scanned = global::Seqeron.Genomics.Infrastructure.EmbossFastaReader.Read(sequenceFasta);
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

    [McpServerTool(Name = "nucmer_clusters", Title = "Genome — MUMmer nucmer Anchor Clusters", ReadOnly = true)]
    [Description("The anchor-clustering stage of MUMmer 3.23 nucmer = `mummer -mum|-mumreference|-maxmatch -b -l L -n | mgaps -l c -s g -d D -f F` (700 related pairs, both strands, all three match modes equal; 758 further mgaps inputs equal): the exact matches (mode 'mum' = maximal unique matches in both sequences, 'mumreference' = unique in the reference, 'maxmatch' = all maximal matches; minMatchLength = -l; only A/C/G/T match) are clustered: matches are sorted by their start in the query, matches inside a repeat are removed and overlapping matches of one diagonal merged, neighbours whose diagonals differ by at most max(diagonalDifference, diagonalFraction * separation) within maxSeparation are joined, and each component is cut into chains of maximum score (match length - overlap - diagonal difference); a chain is reported when its length (sum of the match lengths, or the extent with useExtent) reaches minClusterLength. Coordinates are 1-based as printed by mummer / mgaps; on the reverse strand the query start is a position in the reverse-complemented query; a match's overlap is the number of bases cut off its start because it overlaps the previous match. This is the stage that decides which anchors belong together; the banded extension (postnuc) and delta-filter are not part of this tool. Every option value mgaps accepts is accepted and gives the program's clusters (A4-4, measured on the installed mgaps 3.23: negative or zero separation / cluster length / diagonal difference, NaN, infinite or huge diagonal fraction - the C conversion of the product to long is LONG_MIN, so no fraction is added). One cluster of more than 100 000 matches is refused (SequenceAligner.NucmerMaxClusterMatches: the chain dynamic program of mgaps is quadratic in the cluster).")]
    public static NucmerClustersResult NucmerClusters(
        [Description("Reference DNA sequence.")] string reference,
        [Description("Query DNA sequence.")] string query,
        [Description("Minimum exact match length (nucmer -l, default 20).")] int minMatchLength = 20,
        [Description("Match mode: 'mum' (default), 'mumreference' or 'maxmatch'.")] string mode = "mum",
        [Description("Also cluster the reverse-complement strand of the query (nucmer default true).")] bool bothStrands = true,
        [Description("Largest separation between matches of a cluster (nucmer -g, default 90).")] long maxSeparation = 90,
        [Description("Smallest reported cluster length (nucmer -c, default 65).")] long minClusterLength = 65,
        [Description("Fixed diagonal difference that joins matches (nucmer -D, default 5).")] int diagonalDifference = 5,
        [Description("Fraction of the separation added to the allowed diagonal difference (nucmer -d, default 0.12).")] double diagonalFraction = 0.12,
        [Description("Cluster length = extent in the reference instead of the sum of the match lengths (mgaps -e, default false).")] bool useExtent = false,
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(query));
        var matchMode = (mode ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "mum" => CommonRegionMode.MaximalUnique,
            "mumreference" => CommonRegionMode.MaximalUniqueReference,
            "maxmatch" => CommonRegionMode.MaximalExactMatches,
            _ => throw new ArgumentException("mode must be 'mum', 'mumreference' or 'maxmatch'", nameof(mode)),
        };
        var options = new Seqeron.Genomics.Alignment.NucmerClusterOptions
        {
            MaxSeparation = maxSeparation, MinClusterLength = minClusterLength, DiagonalDifference = diagonalDifference, DiagonalFraction = diagonalFraction, UseExtent = useExtent,
        };
        var clusters = GenomicAnalyzer.FindNucmerClusters(reference, query, minMatchLength, matchMode, options, bothStrands ? StrandSelection.Both : StrandSelection.Forward, cancellationToken)
            .Select(c => new NucmerClusterItem(c.IsReverseComplement ? "reverse" : "forward", c.Cluster.Matches.Select(m => new NucmerClusterMatchItem(m.Start1, m.Start2, m.Length, m.Overlap)).ToArray()))
            .ToArray();
        return new NucmerClustersResult(clusters);
    }
}
