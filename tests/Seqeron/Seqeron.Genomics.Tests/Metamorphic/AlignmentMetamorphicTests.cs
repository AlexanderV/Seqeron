using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Alignment area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ALIGN-GLOBAL-001 — global alignment / Needleman–Wunsch (Alignment).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 35.
///
/// API under test (SequenceAligner.GlobalAlign):
///   GlobalAlign(s1, s2, scoring) computes the optimal end-to-end Needleman–Wunsch
///   alignment score F(m,n) over the recurrence
///       F(i,j) = max( F(i−1,j−1) + (match? Match : Mismatch),
///                     F(i−1,j) + GapExtend,
///                     F(i,j−1) + GapExtend )
///   with linear gap initialisation F(i,0) = i·GapExtend, F(0,j) = j·GapExtend.
///   NOTE: this implementation uses a LINEAR gap model — every gap column costs GapExtend
///   and the affine GapOpen term is NOT applied in the DP. The relations below are stated
///   against that actual cost model. Inputs are upper-cased; an empty input yields Score 0.
///
/// Relations (derived from the recurrence, NOT from output):
///   • SYM (score symmetry): the DP matrix for (b,a) is the transpose of the matrix for
///          (a,b) — Match/Mismatch and the single gap cost are symmetric — so the optimal
///          score is unchanged: score(a,b) = score(b,a).
///   • COMP (identity ⇒ maximum score): aligning s to itself scores exactly |s|·Match
///          (the all-match diagonal), and this is the global maximum over every equal-length
///          partner, since each of the |s| alignable positions contributes at most Match.
///   • MON (more matches ⇒ higher score): among equal-length partners, introducing one more
///          mismatch (one fewer match) strictly lowers the optimal score, because a matched
///          column (+Match) is replaced by the best non-match option, which is worth strictly
///          less than Match.
///   • INV (gap-only insert ⇒ score change = gap penalty): inserting a block of length g into
///          one sequence forces exactly g extra gap columns in the optimum, changing the score
///          by precisely g·GapExtend. Proof: for X (len L) and Y (len L+g), any alignment has
///          gaps_Y = gaps_X + g and ≤ L alignable columns, so
///          score ≤ L·Match + g·GapExtend + gaps_X·(2·GapExtend − Match); since GapExtend < 0
///          and Match > 0 the bound is maximised at gaps_X = 0 and is achieved by matching all
///          of X and spending the block as g gaps — independent of the inserted content.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class AlignmentMetamorphicTests
{
    #region Helpers

    private static readonly Random Rng = new(20260619);

    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>A base guaranteed to differ from <paramref name="c"/>, so flipping yields a mismatch.</summary>
    private static char Flip(char c) => c == 'A' ? 'C' : 'A';

    /// <summary>Returns a copy of <paramref name="seq"/> with exactly <paramref name="count"/> spread-out positions mutated to a different base.</summary>
    private static string MutatePositions(string seq, int count)
    {
        var chars = seq.ToCharArray();
        for (int n = 0; n < count; n++)
        {
            int pos = (int)((long)n * seq.Length / count);   // evenly spaced, distinct positions
            chars[pos] = Flip(chars[pos]);
        }
        return new string(chars);
    }

    private static IEnumerable<string> AlignBodies()
    {
        yield return "ACGTACGTAC";
        yield return "GATTACAGAT";
        yield return "AAAAAAAAAA";
        yield return RandomDna(24);
        yield return RandomDna(40);
    }

    /// <summary>A random core over {A,C} only — so 'G'/'T' flanks can never match anywhere inside it.</summary>
    private static string RandomCoreAC(int length)
    {
        const string bases = "AC";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>The reference characters aligned opposite a non-gap query character (the fitted core window).</summary>
    private static string AlignedCore(AlignmentResult r)
    {
        var sb = new System.Text.StringBuilder();
        for (int i = 0; i < r.AlignedSequence1.Length; i++)
            if (r.AlignedSequence1[i] != '-')
                sb.Append(r.AlignedSequence2[i]);
        return sb.ToString();
    }

    private static readonly ScoringMatrix[] ScoringMatrices =
    {
        SequenceAligner.SimpleDna,
        SequenceAligner.BlastDna,
        SequenceAligner.HighIdentityDna,
    };

    #endregion

    #region SYM — score(a,b) = score(b,a)

    [Test]
    [Description("SYM: the Needleman–Wunsch score is symmetric in its arguments, because swapping the sequences transposes the DP matrix without changing any cell's optimum.")]
    public void GlobalAlign_Score_IsSymmetric()
    {
        var bodies = AlignBodies().ToList();

        foreach (var scoring in ScoringMatrices)
        {
            for (int i = 0; i < bodies.Count; i++)
            {
                for (int j = i; j < bodies.Count; j++)
                {
                    int forward = SequenceAligner.GlobalAlign(bodies[i], bodies[j], scoring).Score;
                    int swapped = SequenceAligner.GlobalAlign(bodies[j], bodies[i], scoring).Score;

                    swapped.Should().Be(forward,
                        because: $"the DP matrix of ({bodies[j]},{bodies[i]}) is the transpose of ({bodies[i]},{bodies[j]}), so the corner optimum is identical");
                }
            }
        }
    }

    #endregion

    #region COMP — aligning a sequence to itself attains the maximum score |s|·Match

    [Test]
    [Description("COMP: self-alignment scores exactly |s|·Match and is ≥ the score against any equal-length partner, since each alignable position contributes at most Match.")]
    public void GlobalAlign_Identity_AttainsMaximumScore()
    {
        foreach (var scoring in ScoringMatrices)
        {
            foreach (var body in AlignBodies())
            {
                int selfScore = SequenceAligner.GlobalAlign(body, body, scoring).Score;

                selfScore.Should().Be(body.Length * scoring.Match,
                    because: "the all-match diagonal aligns every position as a match, contributing |s|·Match with no gaps");

                for (int d = 1; d <= Math.Min(3, body.Length); d++)
                {
                    string variant = MutatePositions(body, d);
                    int variantScore = SequenceAligner.GlobalAlign(body, variant, scoring).Score;

                    variantScore.Should().BeLessThanOrEqualTo(selfScore,
                        because: "no equal-length partner can beat the identity, whose every column already yields the maximal per-column score Match");
                }
            }
        }
    }

    #endregion

    #region MON — one fewer match (one more mismatch) strictly lowers the score

    [Test]
    [Description("MON: among equal-length partners, each additional mismatch strictly decreases the optimal global score, because a +Match column is replaced by a strictly cheaper option.")]
    public void GlobalAlign_MoreMismatches_StrictlyLowersScore()
    {
        foreach (var scoring in ScoringMatrices)
        {
            foreach (var body in AlignBodies())
            {
                int maxMut = Math.Min(5, body.Length);
                int previous = int.MaxValue;

                for (int d = 0; d <= maxMut; d++)
                {
                    string variant = MutatePositions(body, d);
                    int score = SequenceAligner.GlobalAlign(body, variant, scoring).Score;

                    if (d > 0)
                        score.Should().BeLessThan(previous,
                            because: $"the {d}-mismatch variant has one fewer match than the {d - 1}-mismatch variant, and the lost +Match column is replaced by a strictly cheaper alignment column");
                    previous = score;
                }
            }
        }
    }

    #endregion

    #region INV — inserting a block of length g changes the score by exactly g·GapExtend

    [Test]
    [Description("INV: inserting a length-g block (at the end or in the middle) forces exactly g gap columns in the optimum, so the score changes by precisely g·GapExtend — independent of the inserted content.")]
    public void GlobalAlign_GapOnlyInsert_ChangesScoreByGapPenalty()
    {
        foreach (var scoring in ScoringMatrices)
        {
            foreach (var body in AlignBodies())
            {
                int baseScore = SequenceAligner.GlobalAlign(body, body, scoring).Score;

                foreach (int g in new[] { 1, 2, 4 })
                {
                    string block = new string('A', g);

                    // Append the block.
                    string appended = body + block;
                    int appendedScore = SequenceAligner.GlobalAlign(body, appended, scoring).Score;
                    (appendedScore - baseScore).Should().Be(g * scoring.GapExtend,
                        because: $"the {g} trailing inserted bases can only be spent as {g} gap columns at GapExtend each — the {body.Length} real positions still match");

                    // Insert the block in the middle.
                    int mid = body.Length / 2;
                    string inserted = body[..mid] + block + body[mid..];
                    int insertedScore = SequenceAligner.GlobalAlign(body, inserted, scoring).Score;
                    (insertedScore - baseScore).Should().Be(g * scoring.GapExtend,
                        because: "a mid-sequence insertion still forces exactly g extra gap columns: matches cap at |s|, so the optimum matches all of s and spends the block as gaps");
                }
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ALIGN-LOCAL-001 — local alignment / Smith–Waterman (Alignment).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 36.
    //
    // API under test (SequenceAligner.LocalAlign):
    //   Smith–Waterman with the same linear gap cost as GlobalAlign, but with a 0 floor on
    //   every cell, F(i,j) = max(0, diag±, up+gap, left+gap); the reported Score is the
    //   maximum cell, and traceback runs back to the first 0 cell. An empty/best-empty
    //   alignment scores 0.
    //
    // Relations (derived from the recurrence, NOT from output):
    //   • SUB (non-negativity): the 0 floor means every reported local score is ≥ 0 — the
    //          empty alignment is always available as a fallback.
    //   • COMP (identity ⇒ equals global): for two identical sequences the full diagonal is
    //          all matches and never dips below 0, so the local optimum spans the whole
    //          sequence and equals both |s|·Match and the global score.
    //   • MON (extend the matching region ⇒ ≥ score): lengthening a shared identical core
    //          (flanked by symbols that occur in neither partner's counterpart) adds matched
    //          columns, so the local optimum is exactly |core|·Match and strictly increases
    //          with the core length.
    //   • INV (distant non-matching flank ⇒ same local alignment): appending/prepending flanks
    //          built from a symbol absent from the other sequence ('G' on one side, 'T' on the
    //          other) can match nothing, so they neither extend nor outscore the core hit — the
    //          local Score and the aligned substrings are byte-for-byte identical to the
    //          unflanked alignment.
    // ───────────────────────────────────────────────────────────────────────────

    #region SUB — every local score is ≥ 0

    [Test]
    [Description("SUB: the Smith–Waterman 0 floor guarantees a non-negative local score for every input pair, including unrelated sequences.")]
    public void LocalAlign_Score_IsNonNegative()
    {
        var bodies = AlignBodies().ToList();

        foreach (var scoring in ScoringMatrices)
        {
            foreach (var a in bodies)
            {
                foreach (var b in bodies)
                {
                    SequenceAligner.LocalAlign(a, b, scoring).Score.Should().BeGreaterThanOrEqualTo(0,
                        because: "the cell recurrence floors at 0, so the empty local alignment (score 0) is always a valid fallback");
                }
            }
        }
    }

    #endregion

    #region COMP — local alignment of identical sequences equals the global alignment

    [Test]
    [Description("COMP: for identical sequences the local optimum spans the whole sequence, equalling |s|·Match and the global score.")]
    public void LocalAlign_Identity_EqualsGlobalAndMaximum()
    {
        foreach (var scoring in ScoringMatrices)
        {
            foreach (var body in AlignBodies())
            {
                int local = SequenceAligner.LocalAlign(body, body, scoring).Score;
                int global = SequenceAligner.GlobalAlign(body, body, scoring).Score;

                local.Should().Be(body.Length * scoring.Match,
                    because: "the all-match diagonal stays positive throughout, so Smith–Waterman keeps the entire sequence");
                local.Should().Be(global,
                    because: "when the whole sequence is the best local region, the local and global optima coincide");
            }
        }
    }

    #endregion

    #region MON — extending the shared core increases the local score

    [Test]
    [Description("MON: lengthening an identical shared core (with non-matching flanks) yields a local score of exactly |core|·Match that strictly increases with the core length.")]
    public void LocalAlign_ExtendingMatchingRegion_IncreasesScore()
    {
        foreach (var scoring in ScoringMatrices)
        {
            int previous = int.MinValue;

            foreach (int coreLen in new[] { 2, 4, 6, 8, 12 })
            {
                string core = RandomCoreAC(coreLen);
                string seq1 = "G" + core + "G";   // 'G' absent from seq2
                string seq2 = "T" + core + "T";   // 'T' absent from seq1

                int score = SequenceAligner.LocalAlign(seq1, seq2, scoring).Score;

                score.Should().Be(coreLen * scoring.Match,
                    because: "the only matchable region is the shared core, so the local optimum is exactly |core|·Match");
                score.Should().BeGreaterThan(previous,
                    because: "a longer shared core adds matched columns, so the local optimum strictly increases");
                previous = score;
            }
        }
    }

    #endregion

    #region INV — a distant non-matching flank leaves the local alignment unchanged

    [Test]
    [Description("INV: prepending/appending flanks made of a symbol absent from the other sequence cannot match anything, so the local Score and aligned substrings are identical to the unflanked alignment.")]
    public void LocalAlign_DistantNonMatchingFlank_PreservesLocalAlignment()
    {
        foreach (var scoring in ScoringMatrices)
        {
            foreach (int coreLen in new[] { 6, 10, 16 })
            {
                string core = RandomCoreAC(coreLen);

                var bare = SequenceAligner.LocalAlign(core, core, scoring);

                foreach (int flank in new[] { 1, 3, 5 })
                {
                    string g = new string('G', flank);
                    string t = new string('T', flank);
                    string seq1 = g + core + g;   // 'G' never appears in seq2
                    string seq2 = t + core + t;   // 'T' never appears in seq1

                    var flanked = SequenceAligner.LocalAlign(seq1, seq2, scoring);

                    flanked.Score.Should().Be(bare.Score,
                        because: "the flanks match nothing, so they cannot extend or outscore the core's local hit");
                    flanked.AlignedSequence1.Should().Be(bare.AlignedSequence1,
                        because: "the optimal local region is still exactly the shared core, unaffected by unmatchable flanks");
                    flanked.AlignedSequence2.Should().Be(bare.AlignedSequence2,
                        because: "the optimal local region is still exactly the shared core, unaffected by unmatchable flanks");
                }
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ALIGN-SEMI-001 — semi-global / fitting alignment (Alignment).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 37.
    //
    // API under test (SequenceAligner.SemiGlobalAlign):
    //   A fitting alignment — the query seq1 is aligned end-to-end (no free gaps), while the
    //   reference seq2 has FREE end gaps: its prefix (row 0 initialised to 0) and its suffix
    //   (the optimum is the maximum of the LAST row) are skipped at no cost. Linear gap cost.
    //
    // Relations (derived from that fitting definition, NOT from output):
    //   • MON (more matching overlap ⇒ higher score): when the query occurs verbatim inside
    //          the reference, the whole query fits to that window, scoring exactly |query|·Match
    //          (the maximum, since the fully-aligned query offers at most |query| matched
    //          columns). A longer embedded query has more matched columns, so the score
    //          strictly increases.
    //   • INV (extend the non-overlapping part ⇒ same core alignment): the reference prefix
    //          and suffix are free-gap regions, so lengthening them — with symbols ('G') absent
    //          from the {A,C} query, ruling out any rival placement — leaves both the score and
    //          the fitted core window (reference bases opposite the query) unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region MON — a longer embedded query (more overlap) raises the fitting score

    [Test]
    [Description("MON: embedding a longer query verbatim in the reference yields a fitting score of exactly |query|·Match that strictly increases with the query length.")]
    public void SemiGlobalAlign_MoreMatchingOverlap_IncreasesScore()
    {
        foreach (var scoring in ScoringMatrices)
        {
            int previous = int.MinValue;

            foreach (int queryLen in new[] { 3, 5, 8, 12, 16 })
            {
                string query = RandomCoreAC(queryLen);
                string reference = "GGGG" + query + "GGGG";   // 'G' absent from the {A,C} query

                int score = SequenceAligner.SemiGlobalAlign(new DnaSequence(query), new DnaSequence(reference), scoring).Score;

                score.Should().Be(queryLen * scoring.Match,
                    because: "the fully-aligned query fits to its verbatim occurrence, matching every one of its |query| columns");
                score.Should().BeGreaterThan(previous,
                    because: "a longer embedded query contributes more matched columns, so the fitting score strictly increases");
                previous = score;
            }
        }
    }

    #endregion

    #region INV — extending the free-gap reference flanks preserves score and fitted core

    [Test]
    [Description("INV: lengthening the reference prefix/suffix (free-gap regions, built from a symbol absent from the query) changes neither the fitting score nor the reference window fitted to the query.")]
    public void SemiGlobalAlign_ExtendNonOverlappingFlank_PreservesCoreAlignment()
    {
        foreach (var scoring in ScoringMatrices)
        {
            foreach (int queryLen in new[] { 6, 10, 16 })
            {
                string query = RandomCoreAC(queryLen);

                int? refScore = null;
                string? refCore = null;

                foreach (int flank in new[] { 2, 4, 8, 16 })
                {
                    string g = new string('G', flank);
                    string reference = g + query + g;

                    var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(query), new DnaSequence(reference), scoring);

                    result.Score.Should().Be(queryLen * scoring.Match,
                        because: "the flanks are free-gap regions made of an unmatchable symbol, so they add nothing to the score");
                    AlignedCore(result).Should().Be(query,
                        because: "the query still fits exactly to its verbatim occurrence — the fitted reference window is the query itself");

                    refScore ??= result.Score;
                    refCore ??= AlignedCore(result);
                    result.Score.Should().Be(refScore!.Value,
                        because: "extending free-gap flanks must not change the optimal fitting score");
                    AlignedCore(result).Should().Be(refCore,
                        because: "extending free-gap flanks must not move or alter the fitted core window");
                }
            }
        }
    }

    #endregion
}
