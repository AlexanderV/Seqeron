using System.ComponentModel;
using ModelContextProtocol.Server;

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
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "etandem ignores -sbegin / -send (the whole sequence is searched) and numbers a -sreverse sequence on the reverse strand; equicktandem searches the range and numbers from its first position.")] string? sequenceQualifiers = null,
        [Description("EMBOSS report format (-rformat): gff or excel. The report text of the program (the whole file; #!Date is today; sequence name from sequenceName) is returned in the result field report; only for a search without sequenceQualifiers. The other formats of the program (seqtable, table, embl, genbank, motif, ...) are not written yet. Omit for no report text.")] string? reportFormat = null,
        [Description("Sequence name printed in the report (default sequence).")] string sequenceName = "sequence")
    {
        if (reportFormat is not null)
        {
            var format = ParseEmbossReportFormat(reportFormat);
            if (sequenceQualifiers is not null)
                throw new ArgumentException("reportFormat needs no sequenceQualifiers", nameof(reportFormat));
            var plain = FindTandemRepeatsEmboss(sequence, program, minRepeat, maxRepeat, threshold, mismatch, uniform, cancellationToken, null, null);
            string reportProgram = (program ?? "etandem").Trim().ToLowerInvariant();
            string clean = global::Seqeron.Genomics.Infrastructure.EmbossSequenceReader.Clean(sequence);
            var hits = reportProgram == "etandem"
                ? GenomicAnalyzer.FindTandemRepeatsEtandem(RequireEmbossIupacDna(sequence, nameof(sequence)), minRepeat, maxRepeat ?? minRepeat, threshold, mismatch, uniform, cancellationToken)
                : GenomicAnalyzer.FindTandemRepeatsEquicktandem(RequireEmbossIupacDna(sequence, nameof(sequence)), maxRepeat ?? 600, threshold, cancellationToken);
            var report = reportProgram == "etandem" ? global::Seqeron.Genomics.Analysis.EmbossReports.Etandem(sequenceName, clean, hits) : global::Seqeron.Genomics.Analysis.EmbossReports.Equicktandem(sequenceName, clean, hits);
            return plain with { Report = global::Seqeron.Genomics.Infrastructure.EmbossReportWriter.Write(report, format, DateTime.Now) };
        }

        if (sequenceQualifiers is not null)
        {
            var q = global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(sequenceQualifiers);
            string kind = (program ?? "etandem").Trim().ToLowerInvariant();
            IReadOnlyList<EmbossTandemRepeat> qualified = kind switch
            {
                "etandem" => GenomicAnalyzer.FindTandemRepeatsEtandem(sequence, q, minRepeat, maxRepeat ?? minRepeat, threshold, mismatch, uniform, cancellationToken),
                "equicktandem" => GenomicAnalyzer.FindTandemRepeatsEquicktandem(sequence, q, maxRepeat ?? 600, threshold, cancellationToken),
                _ => throw new ArgumentException("program must be etandem or equicktandem", nameof(program)),
            };
            return new FindEmbossTandemRepeatsResult(kind, [.. qualified.Select(h => new EmbossTandemItem(h.Position, h.Length, h.Score, h.RepeatSize, h.RepeatCount, kind == "etandem" ? Math.Round((double)h.Identity, 1) : null, h.Consensus))]);
        }

        var nucleotides = RequireEmbossIupacDna(sequence, nameof(sequence));
        string which = (program ?? "etandem").Trim().ToLowerInvariant();
        IReadOnlyList<EmbossTandemRepeat> found;
        switch (which)
        {
            case "etandem":
                found = GenomicAnalyzer.FindTandemRepeatsEtandem(nucleotides, minRepeat, maxRepeat ?? minRepeat, threshold, mismatch, uniform, cancellationToken);
                break;
            case "equicktandem":
                found = GenomicAnalyzer.FindTandemRepeatsEquicktandem(nucleotides, maxRepeat ?? 600, threshold, cancellationToken);
                break;
            default:
                throw new ArgumentException("program must be etandem or equicktandem", nameof(program));
        }

        bool detailed = which == "etandem";
        var items = found
            .Select(h => new EmbossTandemItem(h.Position, h.Length, h.Score, h.RepeatSize, h.RepeatCount,
                detailed ? Math.Round((double)h.Identity, 1) : null, h.Consensus))
            .ToArray();
        return new FindEmbossTandemRepeatsResult(which, items);
    }

    [McpServerTool(Name = "seqmatchall", Title = "Genome — EMBOSS seqmatchall (all-against-all word matches)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 seqmatchall: the all-against-all comparison of a set of sequences by exact word matches. For every pair i < j with length(i) > wordSize and length(j) >= wordSize every maximal exact match of at least wordSize symbols is reported (forward strand only; overlapping matches and repeated positions are all listed). Nucleotide or protein sequences (any letters, compared literally and case-insensitively); the program's word table leaves out words around its unknown symbol (N for a nucleotide sequence, X otherwise) with a scan that can keep a word that holds an N after a skipped one - reproduced exactly. Positions are 0-based (the program prints 1-based starts). Matches are listed per pair, by position (the program lists them in hash order: compare as sets). Verified against seqmatchall 6.6.0 on 330 random sets (2-5 sequences, nucleotide / N / lower case / protein, word sizes 2-12, 700+ pairs); on A/C/G/T input it equals find_common_regions mode MaximalExactMatches (MUMmer -maxmatch). Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program). More than 4 194 304 maximal word matches are refused with an ArgumentException (SequenceAligner.WordMatchMaxMatches).")]
    public static SeqMatchAllResult SeqMatchAll(
        [Description("The sequences (at least one; nucleotide or protein letters).")] string[] sequences,
        [Description("Word size (seqmatchall -wordsize, default 4, >= 1).")] int wordSize = 4,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "seqmatchall ignores -sbegin / -send; -sreverse (taken by the whole set when its first sequence is nucleotide) only labels the sequences reversed: the words are matched on the sequences as given and numbered from the end. The result then lists the program's start / end positions (reversed: start > end) in 'reportedMatches'.")] string? sequenceQualifiers = null)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required", nameof(sequences));
        if (sequences.Any(string.IsNullOrEmpty))
            throw new ArgumentException("Sequences cannot be null or empty", nameof(sequences));
        if (sequenceQualifiers is not null)
        {
            var reported = GenomicAnalyzer.SeqMatchAll(sequences, global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(sequenceQualifiers), wordSize, cancellationToken)
                .Select(m => new ReportedWordMatchItem(m.Sequence1, m.Sequence2, m.Start1, m.End1, m.Start2, m.End2, m.Length, m.Reversed)).ToArray();
            return new SeqMatchAllResult([], reported);
        }
        var items = GenomicAnalyzer.SeqMatchAll(sequences, wordSize, cancellationToken)
            .Select(m => new WordMatchItem(m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length))
            .ToArray();
        return new SeqMatchAllResult(items);
    }

    [McpServerTool(Name = "fuzzpro", Title = "Protein — EMBOSS fuzzpro pattern search", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 fuzzpro: every window of a protein sequence that matches a PROSITE-like pattern: amino-acid letters (case-insensitive; B = D or N, Z = E or Q, J = I or L; x = any residue), classes [ILV], exclusions {P} (plain letters only: B, Z, J, X inside [] or {} are rejected as in the binary), repeats x(3), length ranges x(2,4) (every matching window length is reported; (0,n) makes the element optional), '-' separators, the anchors '<' (N terminus) and '>' (C terminus), a terminating '.'; e.g. [DE](2)HS{P}x(2)Px(2,4)C. The search is exact (fuzzpro has no mismatch option). A start-anchored pattern with a range reports only its longest match (the binary's regular-expression search class); a single trailing stop '*' is not part of the sequence. Positions are 0-based; hits are ordered by position, then length. Verified against fuzzpro 6.6.0 on 900 random patterns (letters, classes, exclusions, repeats, ranges, anchors; sequences with X, B, Z, U, J, O, * and lower case). For nucleotide patterns see find_pattern_hits (fuzznuc), for a translated search fuzztran. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static FuzzproResult Fuzzpro(
        [Description("Protein sequence (letters; a trailing * is dropped).")] string sequence,
        [Description("The pattern (fuzzpro -pattern).")] string pattern,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "fuzzpro searches the range and counts positions from the start of the original sequence; a protein is never reversed.")] string? sequenceQualifiers = null,
        [Description("EMBOSS report format (-rformat): gff or excel. The report text of the program (the whole file; #!Date is today; sequence name from sequenceName) is returned in the result field report; only for a search without sequenceQualifiers. The other formats of the program (seqtable, table, embl, genbank, motif, ...) are not written yet. Omit for no report text.")] string? reportFormat = null,
        [Description("Sequence name printed in the report (default sequence).")] string sequenceName = "sequence")
    {
        if (reportFormat is not null)
        {
            var format = ParseEmbossReportFormat(reportFormat);
            if (sequenceQualifiers is not null)
                throw new ArgumentException("reportFormat needs no sequenceQualifiers", nameof(reportFormat));
            var plain = Fuzzpro(sequence, pattern, cancellationToken, null, null);
            var hitList = GenomicAnalyzer.FindProteinPatternHits(sequence, ProteinPattern.Parse(pattern), cancellationToken);
            return plain with { Report = global::Seqeron.Genomics.Infrastructure.EmbossReportWriter.Write(global::Seqeron.Genomics.Analysis.EmbossReports.Fuzzpro(sequenceName, sequence, pattern, hitList), format, DateTime.Now) };
        }

        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (string.IsNullOrWhiteSpace(pattern))
            throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
        var parsed = ProteinPattern.Parse(pattern);
        if (sequenceQualifiers is not null)
            return new FuzzproResult(parsed.ToString(), [.. GenomicAnalyzer.FindProteinPatternHits(sequence, global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(sequenceQualifiers), pattern, cancellationToken).Select(h => new ProteinPatternHitItem(h.Position, h.Length))]);
        var hits = GenomicAnalyzer.FindProteinPatternHits(sequence, parsed, cancellationToken)
            .Select(h => new ProteinPatternHitItem(h.Position, h.Length)).ToArray();
        return new FuzzproResult(parsed.ToString(), hits);
    }

    [McpServerTool(Name = "fuzztran", Title = "Genome — EMBOSS fuzztran (protein pattern search after translation)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 fuzztran: translates a nucleotide sequence in the chosen frames (frame 1, 2, 3, -1, -2, -3, F = forward three, R = reverse three, 6 = all; transeq conventions: reverse frame -k is the reverse complement read from its k-th base, dangling bases give X) with an EMBOSS genetic-code table (0 = standard, else the NCBI table number) and searches every translation with a fuzzpro protein pattern (see fuzzpro: letters, x, [ ], { }, repeats, ranges, anchors; exact). Each hit gives the frame, the 1-based protein start / end, the matched residues, the forward-strand nucleotide location (position 0-based, length = 3 per residue) and embossStart / embossEnd, the two numbers the program prints (exact for the forward frames; EMBOSS 6.6.0 prints negative numbers for the reverse frames - a defect of its source-feature mapping, reproduced; use position / length). Hits per frame in the order 1, 2, 3, -1, -2, -3, by protein position, then length. Verified against fuzztran 6.6.0 on 760 random DNA / pattern / frame / table cases (all frame selections, tables 0-11, N / IUPAC codes and lower case). Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static FuzztranResult Fuzztran(
        [Description("Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive).")] string sequence,
        [Description("The protein pattern (fuzztran -pattern).")] string pattern,
        [Description("Frames to translate (fuzztran -frame): 1 (default), 2, 3, -1, -2, -3, F, R or 6.")] string frame = "1",
        [Description("EMBOSS genetic code (fuzztran -table): 0 = standard (default), otherwise an NCBI table number.")] int table = 0,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "fuzztran cuts the (reverse-complemented) text at the begin / end AS GIVEN, numbers the hits from the first position of that text and prints negative numbers for the reverse frames (a defect of the program, reproduced).")] string? sequenceQualifiers = null,
        [Description("EMBOSS report format (-rformat): gff or excel. The report text of the program (the whole file; #!Date is today; sequence name from sequenceName) is returned in the result field report; only for a search without sequenceQualifiers. The other formats of the program (seqtable, table, embl, genbank, motif, ...) are not written yet. Omit for no report text.")] string? reportFormat = null,
        [Description("Sequence name printed in the report (default sequence).")] string sequenceName = "sequence")
    {
        if (reportFormat is not null)
        {
            var format = ParseEmbossReportFormat(reportFormat);
            if (sequenceQualifiers is not null)
                throw new ArgumentException("reportFormat needs no sequenceQualifiers", nameof(reportFormat));
            var plain = Fuzztran(sequence, pattern, frame, table, cancellationToken, null, null);
            var hitList = GenomicAnalyzer.FindTranslatedPatternHits(RequireEmbossIupacDna(sequence, nameof(sequence)), pattern, frame, table, cancellationToken);
            return plain with { Report = global::Seqeron.Genomics.Infrastructure.EmbossReportWriter.Write(global::Seqeron.Genomics.Analysis.EmbossReports.Fuzztran(sequenceName, pattern, hitList), format, DateTime.Now) };
        }

        if (sequenceQualifiers is not null)
        {
            if (string.IsNullOrWhiteSpace(pattern))
                throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
            var qualifiedParsed = ProteinPattern.Parse(pattern);
            var qualifiedHits = GenomicAnalyzer.FindTranslatedPatternHits(sequence, global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(sequenceQualifiers), pattern, frame, table, cancellationToken)
                .Select(h => new TranslatedPatternHitItem(h.Frame, h.ProteinStart, h.ProteinEnd, h.Translation, h.Position, h.Length, h.EmbossStart, h.EmbossEnd)).ToArray();
            return new FuzztranResult(qualifiedParsed.ToString(), qualifiedHits);
        }

        var nucleotides = RequireEmbossIupacDna(sequence, nameof(sequence));
        if (string.IsNullOrWhiteSpace(pattern))
            throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
        var parsed = ProteinPattern.Parse(pattern);
        var hits = GenomicAnalyzer.FindTranslatedPatternHits(nucleotides, pattern, frame, table, cancellationToken)
            .Select(h => new TranslatedPatternHitItem(h.Frame, h.ProteinStart, h.ProteinEnd, h.Translation, h.Position, h.Length, h.EmbossStart, h.EmbossEnd)).ToArray();
        return new FuzztranResult(parsed.ToString(), hits);
    }

    [McpServerTool(Name = "showalign", Title = "Alignment — EMBOSS showalign display", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 showalign: the pretty text display of a multiple alignment, identical to the program (600 random alignments with random options equal). Every row except the reference is converted against the reference (the consensus unless refseq names a row or a 1-based number): show N (default) non-identities, I identities, S similarities (matrix score > 0), D dissimilarities, A all; '.' marks suppressed residues; similarCase (default true) encodes similarity in the case (non-identity / similarity / all as in the program); order I (input), A (names), S (decreasing similarity to the reference). Per block of width columns (default 60): the position numbers (every 10th, relative to the reference residues when a row is the reference, with -n / +n flanks), the tick ruler, the reference first when it is a row, the rows, the consensus (consensus = false hides it; bottom = false puts it first when it is the reference). The consensus is the EMBOSS cons algorithm (plurality in percent, setCase, identity in percent, gaps: a gap can be the consensus) on the alignment as given (case kept); a lower-case letter of the input is kept lower case. margin = width of the name column (0 = no names, -1 = longest name + 1); uppercase = 1-based column ranges put in upper case, e.g. '3-10,20-25'. The set type is that of the first row (nucleotide: A/C/G/T/U + IUPAC letters); matrix defaults to EDNAFULL / EBLOSUM62; shorter rows are padded with trailing gaps. HTML output and highlight are not provided. Rows are read as EMBOSS reads an alignment (A4-3): letters and '*' are kept, '.' '~' '-' are gaps, '?' '#' '+' are the unknown residue, every other character (digits, blanks, '_', punctuation) is dropped, and rows that become unequal are padded with trailing gaps. A display text above 134 217 728 characters (blocks x lines x (margin + width + 1): a huge margin, width 1 over a very long alignment) is refused with an ArgumentException (MultipleAlignmentDisplay.ShowAlignMaxOutputCharacters); a width beyond the alignment shows one block. sbegin / send (A4-9, 400 random ranges equal the program) cut the displayed columns: the consensus is that of the whole alignment and the numbers and ticks keep the positions of the whole alignment.")]
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
        [Description("1-based inclusive column ranges to put in upper case, e.g. '3-10,20-25'.")] string? uppercase = null,
        [Description("Consensus plurality in percent of the total weight (default 50).")] double plurality = 50,
        [Description("Consensus: weight at or below which the letter is lower case (default half the number of sequences).")] double? setCase = null,
        [Description("Consensus: percent of identities required at a position (default 0).")] double identity = 0,
        [Description("A gap character can be the consensus (default true).")] bool gaps = true,
        [Description("-sbegin: first column shown, 1-based (0 = the first; negative = counted from the end; beyond the alignment = the last).")] int sbegin = 0,
        [Description("-send: last column shown, 1-based (0 = the last; negative = counted from the end; before the first shown column = that column).")] int send = 0,
        CancellationToken cancellationToken = default)
    {
        if (names is null || sequences is null || sequences.Length == 0)
            throw new ArgumentException("Names and sequences are required", nameof(sequences));
        if (names.Length != sequences.Length)
            throw new ArgumentException("One name per sequence is required", nameof(names));
        var ranges = new List<(int, int)>();
        foreach (var part in (uppercase ?? string.Empty).Split([',', ';', ' '], StringSplitOptions.RemoveEmptyEntries))
        {
            var ends = part.Split(['-', ':', '='], StringSplitOptions.RemoveEmptyEntries);
            if (ends.Length is < 1 or > 2 || !int.TryParse(ends[0], out int a) || !int.TryParse(ends[^1], out int b) || a < 1 || b < a)
                throw new ArgumentException($"Invalid uppercase range '{part}' (use 3-10)", nameof(uppercase));
            ranges.Add((a, b));
        }

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
        [Description("EMBOSS sequence-input qualifiers with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 and the same with 2 for the second sequence. Measured on EMBOSS 6.6.0: megamerger merges the TRIMMED sequences (positions of the report are those of the trimmed texts; -sreverse1 reverse-complements before the range is cut); -sprotein is refused (the sequences are DNA). Omit for the whole sequences.")] string? sequenceQualifiers = null,
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence2));
        var r = sequenceQualifiers is null
            ? GenomicAnalyzer.MegaMerge(sequence1, sequence2, wordSize, preferFirst, name1, name2, cancellationToken)
            : GenomicAnalyzer.MegaMerge(
                new global::Seqeron.Genomics.Infrastructure.EmbossSequence(name1, sequence1, "", global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(sequenceQualifiers, "1")),
                new global::Seqeron.Genomics.Infrastructure.EmbossSequence(name2, sequence2, "", global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(sequenceQualifiers, "2")),
                wordSize, preferFirst, cancellationToken);
        return new MegaMergerResult(r.MergedSequence, r.Report,
            r.Matches.Select(m => new WordMatchItem(m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length)).ToArray());
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
