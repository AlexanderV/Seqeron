using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Alignment area (global Needleman–Wunsch, local
/// Smith–Waterman).
///
/// Algebraic testing pins the laws the alignment score must obey: symmetry of
/// the score function, the perfect-self-alignment identity, the all-gap identity
/// against the empty sequence, and the zero floor of a local alignment with no
/// shared content.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 35, 36.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Alignment")]
public class AlignmentAlgebraicTests
{
    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ALIGN-GLOBAL-001 — Global (Needleman–Wunsch) alignment (Alignment)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 35.
    //
    // Model: Needleman–Wunsch optimal end-to-end alignment under a linear gap
    //        penalty d = GapExtend (this implementation uses the linear model, so
    //        d is the per-residue gap cost). SimpleDna scoring = (+1, −1, gap −1).
    //   — docs/algorithms/Alignment; SequenceAligner.GlobalAlign.
    //
    // Laws under test (checklist row 35):
    //   • COMM — score(a, b) = score(b, a): the NW recurrence is symmetric.
    //   • ID   — align(x, x) = perfect score = |x| × Match (a fully matched
    //            diagonal, no gaps or mismatches).
    //   • ID   — align(x, "") = d × |x|: the only alignment to the empty sequence
    //            is all-gap, costing the linear gap penalty per residue. (Tested
    //            through the DnaSequence overload, which evaluates the true NW
    //            recurrence; the string overload deliberately short-circuits an
    //            empty argument to AlignmentResult.Empty.)
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// COMM: global alignment score is symmetric — score(a, b) = score(b, a).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GlobalAlign_Commutative_ScoreIsSymmetric()
    {
        return Prop.ForAll(DnaArbitrary(1), DnaArbitrary(1), (a, b) =>
        {
            int ab = SequenceAligner.GlobalAlign(a, b).Score;
            int ba = SequenceAligner.GlobalAlign(b, a).Score;
            return (ab == ba).Label($"score(a,b)={ab} != score(b,a)={ba}");
        });
    }

    /// <summary>
    /// ID: aligning a sequence to itself yields the perfect score |x| × Match
    /// (no mismatch, no gap), under multiple scoring schemes.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GlobalAlign_Identity_SelfAlignmentIsPerfectScore()
    {
        return Prop.ForAll(DnaArbitrary(1), x =>
        {
            int simple = SequenceAligner.GlobalAlign(x, x, SequenceAligner.SimpleDna).Score;
            int blast = SequenceAligner.GlobalAlign(x, x, SequenceAligner.BlastDna).Score;
            return (simple == x.Length * SequenceAligner.SimpleDna.Match
                    && blast == x.Length * SequenceAligner.BlastDna.Match)
                .Label($"self-align not perfect: simple={simple}, blast={blast}, len={x.Length}");
        });
    }

    /// <summary>
    /// ID: aligning a sequence to the empty sequence is the all-gap alignment,
    /// scoring d × |x| with d = GapExtend (linear gap penalty).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GlobalAlign_Identity_AlignToEmptyIsGapPenaltyTimesLength()
    {
        return Prop.ForAll(DnaArbitrary(1), x =>
        {
            var scoring = SequenceAligner.SimpleDna;
            int score = SequenceAligner.GlobalAlign(new DnaSequence(x), new DnaSequence(""), scoring).Score;
            return (score == x.Length * scoring.GapExtend)
                .Label($"align(x,\"\")={score} != d*len={x.Length * scoring.GapExtend}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ALIGN-LOCAL-001 — Local (Smith–Waterman) alignment (Alignment)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 36.
    //
    // Model: Smith–Waterman optimal local alignment; cell values are floored at 0,
    //        so the optimal local score is ≥ 0 and is 0 exactly when no positively
    //        scoring local alignment exists.
    //   — docs/algorithms/Alignment; SequenceAligner.LocalAlign.
    //
    // Laws under test (checklist row 36):
    //   • ID   — no shared content → score 0: two sequences over disjoint
    //            alphabets ({A,T} vs {C,G}) have no matching column, so every cell
    //            is ≤ 0 and the local score floors at 0.
    //   • COMM — score(a, b) = score(b, a): the Smith–Waterman score is symmetric
    //            (even though the reported aligned strings need not be).
    // ═══════════════════════════════════════════════════════════════════════

    private static Arbitrary<string> AtAlphabet(int minLen) =>
        Gen.Elements('A', 'T').ArrayOf().Where(a => a.Length >= minLen).Select(a => new string(a)).ToArbitrary();

    private static Arbitrary<string> CgAlphabet(int minLen) =>
        Gen.Elements('C', 'G').ArrayOf().Where(a => a.Length >= minLen).Select(a => new string(a)).ToArbitrary();

    /// <summary>
    /// ID: disjoint-alphabet sequences share no base, so the local alignment score
    /// is exactly 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LocalAlign_Identity_NoSharedContentIsZero()
    {
        return Prop.ForAll(AtAlphabet(1), CgAlphabet(1), (at, cg) =>
        {
            int score = SequenceAligner.LocalAlign(at, cg).Score;
            return (score == 0).Label($"local score {score} != 0 for disjoint \"{at}\" / \"{cg}\"");
        });
    }

    /// <summary>
    /// COMM: local alignment score is symmetric — score(a, b) = score(b, a).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LocalAlign_Commutative_ScoreIsSymmetric()
    {
        return Prop.ForAll(DnaArbitrary(1), DnaArbitrary(1), (a, b) =>
        {
            int ab = SequenceAligner.LocalAlign(a, b).Score;
            int ba = SequenceAligner.LocalAlign(b, a).Score;
            return (ab == ba && ab >= 0).Label($"local score(a,b)={ab} != score(b,a)={ba}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ALIGN-STATS-001 — Alignment statistics (Alignment), row 226.
    //
    // Model: alignment identity = identical columns / alignment length × 100. A
    //        self-alignment is 100% identical, and swapping the two sequences
    //        transposes the alignment columns, leaving the identity unchanged.
    //   — docs/algorithms/Alignment; SequenceAligner.CalculateStatistics.
    //
    // Laws (row 226): ID — identity(align(x, x)) = 100.
    //                 COMM — symmetric: the optimal alignment the statistics
    //                 summarise has a symmetric score (the underlying commutative
    //                 invariant). NOTE: the per-alignment Identity percentage itself
    //                 is NOT a symmetric function, because align(a,b) and align(b,a)
    //                 may select different co-optimal tracebacks with different gap
    //                 placements — only the optimal score is guaranteed symmetric.
    // ═══════════════════════════════════════════════════════════════════════

    [FsCheck.NUnit.Property]
    public Property AlignStats_Identity_SelfAlignmentIsFullIdentity()
    {
        return Prop.ForAll(DnaArbitrary(1), x =>
        {
            var stats = SequenceAligner.CalculateStatistics(SequenceAligner.GlobalAlign(x, x));
            return (System.Math.Abs(stats.Identity - 100.0) < 1e-9).Label($"identity(x,x)={stats.Identity}");
        });
    }

    /// <summary>
    /// COMM: the optimal alignment score that the statistics summarise is symmetric,
    /// score(a,b) = score(b,a). (Per-alignment identity is traceback-dependent and
    /// therefore not asserted symmetric — see the unit note above.)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlignStats_Commutative_OptimalScoreSymmetric()
    {
        return Prop.ForAll(DnaArbitrary(1), DnaArbitrary(1), (a, b) =>
        {
            var ab = SequenceAligner.GlobalAlign(a, b);
            var ba = SequenceAligner.GlobalAlign(b, a);
            return (ab.Score == ba.Score).Label($"score(a,b)={ab.Score} != score(b,a)={ba.Score}");
        });
    }

    /// <summary>
    /// COMM (statistics level): a self-alignment is symmetric and maximal — its
    /// identity is 100 and it has zero mismatches, independent of argument order.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlignStats_SelfAlignmentHasNoMismatches()
    {
        return Prop.ForAll(DnaArbitrary(1), x =>
        {
            var stats = SequenceAligner.CalculateStatistics(SequenceAligner.GlobalAlign(x, x));
            return (stats.Mismatches == 0 && stats.Gaps == 0).Label($"self-align stats: {stats.Mismatches} mm, {stats.Gaps} gaps");
        });
    }
}
