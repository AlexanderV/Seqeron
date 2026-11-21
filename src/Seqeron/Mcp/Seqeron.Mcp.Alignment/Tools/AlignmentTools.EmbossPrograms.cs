using System.ComponentModel;
using System.Globalization;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-28c / A2-29b: EMBOSS supermatcher and needleall, and the ksw2 extension functions; delegate to SequenceAligner.Supermatcher / NeedleAll and Ksw2Aligner.
public partial class AlignmentTools
{
    private static global::Seqeron.Genomics.Infrastructure.EmbossSequence[] WithQualifiers(global::Seqeron.Genomics.Infrastructure.EmbossSequence[] sequences, string program, string? sequenceQualifiers, int input)
    {
        if (sequenceQualifiers is null)
            return sequences;
        var q = global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.ParseInputs(program, sequenceQualifiers)[input];
        return [.. sequences.Select(s => s with { Qualifiers = q })];
    }

    private static global::Seqeron.Genomics.Infrastructure.EmbossSequence[] ToEmbossSequences(string[] sequences, string[]? names, string prefix, string paramName)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        if (sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required.", paramName);
        var list = new global::Seqeron.Genomics.Infrastructure.EmbossSequence[sequences.Length];
        for (int i = 0; i < list.Length; i++)
        {
            RequireSequence(sequences[i], paramName);
            string name = names is not null && i < names.Length && !string.IsNullOrWhiteSpace(names[i]) ? names[i] : prefix + (i + 1).ToString(CultureInfo.InvariantCulture);
            list[i] = new global::Seqeron.Genomics.Infrastructure.EmbossSequence(name, sequences[i]);
        }

        return list;
    }

    private static EmbossProgramAlignResult ToProgramResult(global::Seqeron.Genomics.Infrastructure.EmbossProgramResult r, string? formatted = null)
    {
        var items = r.Alignments.Select(a => new PairAlignmentItem(a.First.Name, a.Second.Name, a.First.Row, a.Second.Row, a.First.Offset, a.Second.Offset, a.First.Offend, a.Second.Offend, a.Score)).ToArray();
        return new EmbossProgramAlignResult(items, items.Length, r.ErrorText, formatted);
    }

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
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.Supermatcher(queries, targets, ResolveMatrix(matrix), gapOpen, gapExtend, wordLength, width, minScore);
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
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.NeedleAll(a, b, ResolveMatrix(matrix), gapOpen, gapExtend, endWeight, endOpen, endExtend, minScore, brief);
        return ToProgramResult(r);
    }

    [McpServerTool(Name = "ksw2_extend", Title = "Pairwise — ksw2 Extension Alignment (z-drop, end bonus, band)", ReadOnly = true)]
    [Description("The ksw2 / minimap2 extension-alignment functions: global or extension alignment of a query against a target with affine (q + k*e) or dual-affine (min(q + k*e, q2 + k*e2)) gap cost, a diagonal band, the Z-drop rule, the end bonus, the ambiguous-base (N) penalty and CIGAR output. variant: 'extz' (ksw_extz, scalar affine), 'extz2' (ksw_extz2_sse, minimap2 affine), 'extd' (ksw_extd, scalar dual affine), 'extd2' (ksw_extd2_sse, minimap2 -O/-E dual affine; the default). Bases are encoded A,C,G,T = 0..3 and anything else = N (4); match / mismatch / ambiguous give the matrix as minimap2 builds it (-A, -B as a penalty, -N penalty of the wildcard). Penalties are POSITIVE (minimap2 -O q,q2 -E e,e2; a gap of length k costs q + k*e). w = band half-width (-1: none; the SIMD variants need w >= 1), zdrop = off-diagonal drop-off to stop an extension (-1: none; minimap2 -z), endBonus (minimap2 -s... as the extension end bonus; SIMD variants, extensionOnly). extensionOnly = the path ends at the best cell (or at the query end when mqe + endBonus > max); otherwise it is global (corner to corner). Returned: max, zdropped, maxQ/maxT, mqe (best score reaching the query end) with mqeT, mte with mteQ, score (corner score; -1073741824 = not reached), reachEnd and the CIGAR (M, I = query only, D = target only). VERIFICATION: exact cases (no Z-drop, no approximate mode) equal an independent dynamic program on 330 cases per variant (max, mqe, mte, score; every CIGAR re-scored along its path); the Z-drop and approximate cases have not yet been compared with the real C program (it could not be run in the development environment).")]
    public static Ksw2ExtendResult Ksw2Extend(
        [Description("Query sequence (ACGT, other letters = N).")] string query,
        [Description("Target sequence.")] string target,
        [Description("'extd2' (default), 'extz2', 'extd' or 'extz'.")] string variant = "extd2",
        [Description("Score of a match (minimap2 -A).")] int match = 2,
        [Description("Penalty of a mismatch, positive (minimap2 -B).")] int mismatch = 4,
        [Description("Penalty of the ambiguous base N, positive (minimap2 -N / sc_ambi); 0 gives -e.")] int ambiguous = 1,
        [Description("Gap open penalty q (minimap2 -O first value), >= 0.")] int gapOpen = 4,
        [Description("Gap extension penalty e (minimap2 -E first value), >= 1.")] int gapExtend = 2,
        [Description("Second gap open penalty q2 (dual variants; minimap2 -O second value).")] int gapOpen2 = 24,
        [Description("Second gap extension penalty e2 (dual variants; minimap2 -E second value).")] int gapExtend2 = 1,
        [Description("Band half-width (-1 = none).")] int band = -1,
        [Description("Z-drop (-1 = none; minimap2 -z).")] int zdrop = -1,
        [Description("End bonus (SIMD variants, extensionOnly).")] int endBonus = 0,
        [Description("Extension alignment: the path ends at the best cell (or the query end with the end bonus).")] bool extensionOnly = false,
        [Description("Right-align gaps (KSW_EZ_RIGHT).")] bool rightAlignGaps = false,
        [Description("Score only, no CIGAR.")] bool scoreOnly = false,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(query, nameof(query));
        RequireSequence(target, nameof(target));
        if (gapOpen < 0 || gapExtend < 1 || gapOpen > 127 || gapExtend > 127)
            throw new ArgumentOutOfRangeException(nameof(gapOpen), "Gap open must be 0..127 and gap extension 1..127 (they are bytes in ksw2).");
        var flags = global::Seqeron.Genomics.Alignment.Ksw2Option.None;
        if (extensionOnly)
            flags |= global::Seqeron.Genomics.Alignment.Ksw2Option.ExtensionOnly;
        if (rightAlignGaps)
            flags |= global::Seqeron.Genomics.Alignment.Ksw2Option.Right;
        if (scoreOnly)
            flags |= global::Seqeron.Genomics.Alignment.Ksw2Option.ScoreOnly;
        byte[] q = global::Seqeron.Genomics.Alignment.Ksw2Aligner.EncodeNucleotides(query);
        byte[] t = global::Seqeron.Genomics.Alignment.Ksw2Aligner.EncodeNucleotides(target);
        var mat = global::Seqeron.Genomics.Alignment.Ksw2Aligner.SimpleMatrix(5, match, -mismatch, ambiguous);
        var r = (variant ?? "extd2").Trim().ToLowerInvariant() switch
        {
            "extz" => global::Seqeron.Genomics.Alignment.Ksw2Aligner.ExtendZ(q, t, 5, mat, gapOpen, gapExtend, band, zdrop, flags),
            "extz2" => global::Seqeron.Genomics.Alignment.Ksw2Aligner.ExtendZ2(q, t, 5, mat, gapOpen, gapExtend, band, zdrop, endBonus, flags),
            "extd" => global::Seqeron.Genomics.Alignment.Ksw2Aligner.ExtendD(q, t, 5, mat, gapOpen, gapExtend, gapOpen2, gapExtend2, band, zdrop, flags),
            "extd2" => global::Seqeron.Genomics.Alignment.Ksw2Aligner.ExtendD2(q, t, 5, mat, gapOpen, gapExtend, gapOpen2, gapExtend2, band, zdrop, endBonus, flags),
            _ => throw new ArgumentException($"variant must be extz, extz2, extd or extd2, got '{variant}'.", nameof(variant)),
        };
        return new Ksw2ExtendResult(r.Max, r.ZDropped, r.MaxQ, r.MaxT, r.Mqe, r.MqeT, r.Mte, r.MteQ, r.Score, r.ReachEnd, r.CigarString);
    }
}
