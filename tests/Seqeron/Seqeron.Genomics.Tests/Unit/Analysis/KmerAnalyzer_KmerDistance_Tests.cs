// KMER-DIST-001 — K-mer Euclidean Distance
// Evidence: docs/Evidence/KMER-DIST-001-Evidence.md
// TestSpec: tests/TestSpecs/KMER-DIST-001.md
// Source: Zielezinski A, Vinga S, Almeida J, Karlowski WM (2017). Genome Biology 18:186 (Fig. 1).
//         Lau AK et al. (2022). NAR Genom Bioinform (k-mer frequency = count / (L - k + 1)).
//         Boden M et al. (2014). Bioinformatics 30(14) (Euclidean over relative-frequency vectors).

using System;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Canonical test class for KMER-DIST-001: K-mer Euclidean Distance.
/// Verifies <see cref="KmerAnalyzer.KmerDistance(string, string, int)"/> — the Euclidean
/// distance between normalized k-mer frequency vectors — against the Zielezinski et al. (2017)
/// Figure 1 worked example and derivations from the Lau et al. (2022) frequency definition.
/// </summary>
[TestFixture]
public class KmerAnalyzer_KmerDistance_Tests
{
    // Tolerance for floating-point Euclidean comparisons (exact values are derived in Evidence).
    private const double Tolerance = 1e-10;

    // Expected value for the Zielezinski 2017 Fig.1 example (x=ATGTGTG, y=CATGTG, k=3):
    // freq diffs (-0.05,-0.25,0.15,0.15) => sqrt(0.11).
    private static readonly double Fig1Distance = Math.Sqrt(0.11);

    #region KmerDistance — MUST

    // M1 — Zielezinski et al. (2017) Fig.1: c_X=(1,0,2,2), c_Y=(1,1,1,1) over {ATG,CAT,GTG,TGT};
    // frequency vectors f_X=(0.2,0,0.4,0.4), f_Y=(0.25,0.25,0.25,0.25) => sqrt(0.11).
    [Test]
    public void KmerDistance_ZielezinskiFig1Example_ReturnsSqrt011()
    {
        double distance = KmerAnalyzer.KmerDistance("ATGTGTG", "CATGTG", 3);

        Assert.That(distance, Is.EqualTo(Fig1Distance).Within(Tolerance),
            "Fig.1 frequency vectors give sum of squared diffs 0.11, so distance must be sqrt(0.11) ≈ 0.3316624790.");
    }

    // M2 — Identical sequences yield distance 0 (Zielezinski 2017 Fig.1: "identical sequences yield a distance of 0").
    [Test]
    public void KmerDistance_IdenticalSequences_ReturnsZero()
    {
        double distance = KmerAnalyzer.KmerDistance("ATGTGTG", "ATGTGTG", 3);

        Assert.That(distance, Is.EqualTo(0.0).Within(Tolerance),
            "Equal frequency vectors give zero sum of squares; identical sequences must give distance 0 (INV-01).");
    }

    // M3 — Single substitution, k=1. AAAA -> f(A)=1; AAAT -> f(A)=0.75, f(T)=0.25.
    // distance = sqrt((1-0.75)^2 + (0-0.25)^2) = sqrt(0.125).
    [Test]
    public void KmerDistance_SingleSubstitutionK1_ReturnsSqrt0125()
    {
        double distance = KmerAnalyzer.KmerDistance("AAAA", "AAAT", 1);

        Assert.That(distance, Is.EqualTo(Math.Sqrt(0.125)).Within(Tolerance),
            "k=1 frequencies (A:1) vs (A:0.75,T:0.25) give sum of squares 0.125, so distance must be sqrt(0.125) ≈ 0.3535533906.");
    }

    // M4 — Symmetry: d(x,y) == d(y,x) (Euclidean is a metric; INV-02).
    [Test]
    public void KmerDistance_SwappedArguments_ReturnsSameValue()
    {
        double forward = KmerAnalyzer.KmerDistance("ATGTGTG", "CATGTG", 3);
        double reverse = KmerAnalyzer.KmerDistance("CATGTG", "ATGTGTG", 3);

        Assert.Multiple(() =>
        {
            Assert.That(reverse, Is.EqualTo(forward).Within(Tolerance),
                "Euclidean distance is symmetric: swapping arguments must not change the value (INV-02).");
            Assert.That(forward, Is.EqualTo(Fig1Distance).Within(Tolerance),
                "Both orderings must equal the Fig.1 distance sqrt(0.11).");
        });
    }

    // M5 — Disjoint single-k-mer sequences: AAAA (k=2)->AA freq 1; TTTT->TT freq 1.
    // Vectors (1,0) and (0,1) => sqrt(2). Also demonstrates non-negativity (INV-03/INV-04).
    [Test]
    public void KmerDistance_DisjointSingleKmers_ReturnsSqrt2()
    {
        double distance = KmerAnalyzer.KmerDistance("AAAA", "TTTT", 2);

        Assert.Multiple(() =>
        {
            Assert.That(distance, Is.EqualTo(Math.Sqrt(2.0)).Within(Tolerance),
                "Disjoint single-k-mer frequency vectors (1,0) and (0,1) give distance sqrt(2) ≈ 1.4142135624 (INV-04).");
            Assert.That(distance, Is.GreaterThanOrEqualTo(0.0),
                "Euclidean distance is always non-negative (INV-03).");
        });
    }

    #endregion

    #region KmerDistance — SHOULD

    // S1 — Disjoint single-k-mer sequences at k=2 (same closed form as M5, distinct inputs).
    [Test]
    public void KmerDistance_AllSameOppositeBases_ReturnsSqrt2()
    {
        double distance = KmerAnalyzer.KmerDistance("CCCCC", "GGGGG", 2);

        Assert.That(distance, Is.EqualTo(Math.Sqrt(2.0)).Within(Tolerance),
            "CCCCC has only CC (freq 1) and GGGGG only GG (freq 1); disjoint vectors give sqrt(2) (INV-04).");
    }

    // S2 — One sequence too short for k => empty (zero) vector; other has one k-mer freq 1.
    // distance = sqrt(1^2) = 1.0 (ASSUMPTION A2).
    [Test]
    public void KmerDistance_OneSequenceShorterThanK_ReturnsNormOfOther()
    {
        double distance = KmerAnalyzer.KmerDistance("ACGT", "AAAAAA", 5);

        Assert.That(distance, Is.EqualTo(1.0).Within(Tolerance),
            "ACGT (len 4 < k=5) is the zero vector; AAAAAA gives AAAAA freq 1; distance = sqrt(1) = 1.0 (A2).");
    }

    // S3 — Case-insensitivity: lower-case input is upper-cased before counting (ASSUMPTION A1).
    [Test]
    public void KmerDistance_LowerCaseInput_EqualsUpperCaseResult()
    {
        double distance = KmerAnalyzer.KmerDistance("atgtgtg", "CATGTG", 3);

        Assert.That(distance, Is.EqualTo(Fig1Distance).Within(Tolerance),
            "Inputs are upper-cased before counting, so lower-case x must reproduce the Fig.1 distance sqrt(0.11) (A1).");
    }

    #endregion

    #region KmerDistance — COULD (validation / boundary)

    // C1 — k <= 0 is invalid and must throw.
    [Test]
    public void KmerDistance_NonPositiveK_Throws()
    {
        Assert.That(() => KmerAnalyzer.KmerDistance("ACGT", "ACGT", 0),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(),
            "k must be positive; k=0 must throw ArgumentOutOfRangeException.");
    }

    // C2 — Both sequences empty => both zero vectors => distance 0 (ASSUMPTION A2).
    [Test]
    public void KmerDistance_BothEmpty_ReturnsZero()
    {
        double distance = KmerAnalyzer.KmerDistance("", "", 3);

        Assert.That(distance, Is.EqualTo(0.0).Within(Tolerance),
            "Two empty sequences are both the zero frequency vector, so the distance is 0 (A2).");
    }

    // C3 — Null input is documented as allowed (empty/zero vector, ASSUMPTION A2). A null
    // sequence vs a populated one must equal the L2 norm of the populated one's frequency
    // vector: AAAAA (k=3) -> AAA freq 1 => distance = sqrt(1) = 1.0.
    [Test]
    public void KmerDistance_NullInput_TreatedAsZeroVector()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.KmerDistance(null!, "AAAAA", 3), Is.EqualTo(1.0).Within(Tolerance),
                "Null sequence is the zero vector; AAAAA gives AAA freq 1, so distance = sqrt(1) = 1.0 (A2).");
            Assert.That(KmerAnalyzer.KmerDistance(null!, null!, 3), Is.EqualTo(0.0).Within(Tolerance),
                "Two null sequences are both the zero frequency vector, so the distance is 0 (A2).");
        });
    }

    #endregion
}
