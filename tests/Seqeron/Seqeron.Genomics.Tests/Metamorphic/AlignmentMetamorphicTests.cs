using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;
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
}
