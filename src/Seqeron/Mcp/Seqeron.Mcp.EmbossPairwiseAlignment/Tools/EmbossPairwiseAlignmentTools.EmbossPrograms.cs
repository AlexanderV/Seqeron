// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Alignment.Tools;
using static Seqeron.Mcp.Alignment.Tools.AlignmentTools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

public partial class EmbossPairwiseAlignmentTools
{
    [McpServerTool(Name = "supermatcher_align", Title = "Pairwise — Seeded Local Alignment of Large Sequences (EMBOSS supermatcher)", ReadOnly = true)]
    [Description("EMBOSS 6.6 supermatcher: approximate local alignment of large sequences. The words of wordLength residues of the TARGET sequences (-bsequence) are searched in every QUERY sequence (-asequence) with the Rabin-Karp multi-pattern search of embword.c; the longest word match (the first of equal length) fixes the diagonal, and a banded Smith-Waterman of `width` cells (embAlignPathCalcSWFast, single precision) along that diagonal, from the sequence starts to the end of the shorter overhang, gives the alignment, reported only when its score exceeds minScore. A pair without a shared word is reported in errorText ('No wordmatch start points for A vs B. No alignment'). The result is the program's heuristic alignment (the best one inside the band), not the Smith-Waterman optimum; use local_align_affine / local_align_top_n for that. Alignments are returned in the program's order (for each query, every target). VERIFICATION: output files (every -aformat) and error files identical to the installed supermatcher 6.6.0 on 260 / 260 committed random DNA / IUPAC / protein sets (1-3 queries x 1-4 targets, 4-600 residues, every option: matrix, gapOpen, gapExtend, wordLength 3-15, width 1-100, minScore). Use write_alignment (program 'supermatcher') for the EMBOSS output files. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static EmbossProgramAlignResult SupermatcherAlign(
        [Description("The -asequence sequences (first member of each pair).")] string[] querySequences,
        [Description("The -bsequence sequence set (second member of each pair).")] string[] targetSequences,
        [Description("Optional names of the query sequences (default query1, query2, ...).")] string[]? queryNames = null,
        [Description("Optional names of the target sequences (default target1, ...).")] string[]? targetNames = null,
        [Description("Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62.")] string? matrix = null,
        [Description("Gap opening penalty (-gapopen, 0 - 100): the cost of a gap of length 1.")] double gapOpen = 10,
        [Description("Gap extension penalty (-gapextend, 0 - 10).")] double gapExtend = 0.5,
        [Description("Word length of the seed matching (-wordlen, at least 3).")] int wordLength = 6,
        [Description("Alignment band width (-width, at least 1).")] int width = 16,
        [Description("Minimum alignment score to report an alignment (-minscore).")] double minScore = 0,
        CancellationToken cancellationToken = default,
        [Description(SequenceQualifiersDescription + " Here 1 = the query sequences, 2 = the target sequences (every member takes the same qualifiers).")] string? sequenceQualifiers = null)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var queries = WithQualifiers(ToEmbossSequences(querySequences, queryNames, "query", nameof(querySequences)), "supermatcher", sequenceQualifiers, 0);
        var targets = WithQualifiers(ToEmbossSequences(targetSequences, targetNames, "target", nameof(targetSequences)), "supermatcher", sequenceQualifiers, 1);
        var r = EmbossPairwisePrograms.Supermatcher(queries, targets, ResolveMatrix(matrix), gapOpen, gapExtend, wordLength, width, minScore);
        return ToProgramResult(r);
    }
    [McpServerTool(Name = "needleall_align", Title = "Pairwise — Many-to-Many Global Alignment (EMBOSS needleall)", ReadOnly = true)]
    [Description("EMBOSS 6.6 needleall: the global alignment (Needleman-Wunsch with affine gaps, the program's own scoring: gap of length k costs gapOpen + (k-1)*gapExtend, end gaps free unless endWeight) of every sequence of the set A against every sequence of B, for each B sequence every A sequence, exactly as needle aligns each pair (same rows, ties and single-precision scores). A pair whose score is not above minScore is not reported; the program's error-file line ('Alignment score (S) is less than minimum score(M) for sequences A vs B') is returned in errorText. brief = false adds the Longest_/Shortest_ Identity / Similarity percentages to the sub-header of the formatted output. VERIFICATION: output files (every -aformat) and error files identical to the installed needleall 6.6.0 on 260 / 260 committed random DNA / IUPAC / protein sets (1-3 x 1-4 sequences, every option: matrix, gapOpen, gapExtend, endWeight / endOpen / endExtend, minScore, brief). Use write_alignment (program 'needle') for the EMBOSS output files. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static EmbossProgramAlignResult NeedleallAlign(
        [Description("The -asequence set (first member of each pair).")] string[] setA,
        [Description("The -bsequence sequences (second member of each pair).")] string[] setB,
        [Description("Optional names of the A sequences (default a1, a2, ...).")] string[]? namesA = null,
        [Description("Optional names of the B sequences (default b1, b2, ...).")] string[]? namesB = null,
        [Description("Optional EMBOSS data matrix (-datafile). Omit for EDNAFULL when every A sequence is nucleotide, EBLOSUM62 otherwise.")] string? matrix = null,
        [Description("Gap opening penalty (-gapopen, 1 - 100).")] double gapOpen = 10,
        [Description("Gap extension penalty (-gapextend, 0 - 10).")] double gapExtend = 0.5,
        [Description("-endweight: apply the end-gap penalties.")] bool endWeight = false,
        [Description("-endopen (1 - 100).")] double endOpen = 10,
        [Description("-endextend (0 - 10).")] double endExtend = 0.5,
        [Description("-minscore (-10 - 100): minimum score to report an alignment.")] double minScore = 0,
        [Description("-brief: false adds the Longest_/Shortest_ Identity / Similarity lines to the sub-header.")] bool brief = true,
        CancellationToken cancellationToken = default,
        [Description(SequenceQualifiersDescription + " Here 1 = the A set, 2 = the B sequences (every member takes the same qualifiers).")] string? sequenceQualifiers = null)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var a = WithQualifiers(ToEmbossSequences(setA, namesA, "a", nameof(setA)), "needleall", sequenceQualifiers, 0);
        var b = WithQualifiers(ToEmbossSequences(setB, namesB, "b", nameof(setB)), "needleall", sequenceQualifiers, 1);
        var r = EmbossPairwisePrograms.NeedleAll(a, b, ResolveMatrix(matrix), gapOpen, gapExtend, endWeight, endOpen, endExtend, minScore, brief);
        return ToProgramResult(r);
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
        [Description("EMBOSS sequence-input qualifiers with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 and the same with 2 for the second sequence. Measured on EMBOSS 6.6.0: megamerger merges the TRIMMED sequences (positions of the report are those of the trimmed texts; -sreverse1 reverse-complements before the range is cut); -sprotein is refused (the sequences are DNA). Omit for the whole sequences." + " " + EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null,
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence2));
        var q = sequenceQualifiers is null ? null : EmbossSequenceQualifiers.ParseInputs("megamerger", sequenceQualifiers);
        var r = q is null
            ? EmbossPairwisePrograms.MegaMerge(sequence1, sequence2, wordSize, preferFirst, name1, name2, cancellationToken)
            : EmbossPairwisePrograms.MegaMerge(
                new EmbossSequence(name1, sequence1, "", q[0]),
                new EmbossSequence(name2, sequence2, "", q[1]),
                wordSize, preferFirst, cancellationToken);
        return new MegaMergerResult(r.MergedSequence, r.Report,
            r.Matches.Select(m => new WordMatchItem(m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length)).ToArray());
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
        catch (EmbossReportException e)
        {
            throw new ArgumentException(e.Message, nameof(reportFormat));
        }
    }
}
