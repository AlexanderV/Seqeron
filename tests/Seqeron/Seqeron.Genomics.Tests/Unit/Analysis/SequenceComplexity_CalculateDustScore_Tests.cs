// SEQ-COMPLEX-DUST-001 — DUST Score + symmetric DUST (SDUST) masking
// Evidence: docs/Evidence/SEQ-COMPLEX-DUST-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-COMPLEX-DUST-001.md
// Source: Morgulis A, Gertz EM, Schäffer AA, Agarwala R (2006). J Comput Biol 13(5):1028–1040,
//         doi:10.1089/cmb.2006.13.1028; NCBI dustmasker (symdust.cpp); lh3/sdust (sdust.c).
//
// Spec: for ℓ overlapping triplets, score(x) = (Σ_t c_t·(c_t−1)/2) / (ℓ − 1), ℓ = L − 2.
// The ℓ − 1 normaliser is the one thresholded by both reference implementations
// (symdust: 10·r > level·(ℓ−1); sdust find_perfect: new_r·10 > T·new_l, new_l = ℓ − 1)
// and is confirmed numerically: sdust -t 20 masks a 7-base A-run (score 10/4 = 2.5 > 2)
// but not a 6-base run (6/3 = 2.0), and masks (AC)×6 (20/9 > 2) but not ACACACACACA (16/8 = 2.0).
// A divisor of ℓ (= L − 2) would leave both 7-A and (AC)×6 unmasked (2.0, 2.0) — these tests fail it.
// Masking expectations are the verbatim output of the compiled lh3/sdust reference binary
// (github.com/lh3/sdust master, gcc -O2), NOT the implementation's output.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceComplexity_CalculateDustScore_Tests
{
    #region CalculateDustScore(DnaSequence, k) — canonical exact values

    // M1 — AAAAAA (L=6, ℓ=4): AAA=4 ⇒ Σ=6, /(ℓ−1)=3 ⇒ 2.0 (= (L−2)/2).
    [Test]
    public void CalculateDustScore_Homopolymer6_Returns2()
    {
        double score = SequenceComplexity.CalculateDustScore(new DnaSequence("AAAAAA"));

        Assert.That(score, Is.EqualTo(2.0).Within(1e-10),
            "AAAAAA: AAA×4 ⇒ Σ c(c−1)/2 = 6, divided by ℓ−1 = 3 ⇒ 2.0.");
    }

    // M2 — ACGTACGT (ℓ=6): ACG=2,CGT=2,GTA=1,TAC=1 ⇒ Σ=2, /5 ⇒ 0.4.
    [Test]
    public void CalculateDustScore_RepeatedTetramer_Returns0Point4()
    {
        double score = SequenceComplexity.CalculateDustScore(new DnaSequence("ACGTACGT"));

        Assert.That(score, Is.EqualTo(0.4).Within(1e-10),
            "ACGTACGT: ACG=2,CGT=2 contribute 1 each ⇒ Σ=2, divided by ℓ−1 = 5 ⇒ 0.4.");
    }

    // M3 — ATGC (ℓ=2): all distinct ⇒ Σ=0 ⇒ 0 (INV-2).
    [Test]
    public void CalculateDustScore_AllDistinctTriplets_ReturnsZero()
    {
        double score = SequenceComplexity.CalculateDustScore(new DnaSequence("ATGC"));

        Assert.That(score, Is.EqualTo(0.0).Within(1e-10),
            "ATGC has two distinct triplets ⇒ Σ c(c−1)/2 = 0 ⇒ score 0.");
    }

    // M4 — ACACACAC (ℓ=6): ACA=3,CAC=3 ⇒ Σ=6, /5 ⇒ 1.2.
    [Test]
    public void CalculateDustScore_DinucleotideRepeat_Returns1Point2()
    {
        double score = SequenceComplexity.CalculateDustScore(new DnaSequence("ACACACAC"));

        Assert.That(score, Is.EqualTo(1.2).Within(1e-10),
            "ACACACAC: ACA=3,CAC=3 ⇒ Σ=6, divided by ℓ−1 = 5 ⇒ 1.2.");
    }

    // M5 — A×10 (ℓ=8): Σ=28, /7 ⇒ 4.0 = (L−2)/2 (INV-5).
    [Test]
    public void CalculateDustScore_Homopolymer10_Returns4()
    {
        double score = SequenceComplexity.CalculateDustScore(new DnaSequence("AAAAAAAAAA"));

        Assert.That(score, Is.EqualTo(4.0).Within(1e-10),
            "A×10: AAA×8 ⇒ Σ=28, divided by ℓ−1 = 7 ⇒ 4.0 = (L−2)/2.");
    }

    // M7 — AATAATAA (ℓ=6): AAT=2,ATA=2,TAA=2 ⇒ Σ=3, /5 ⇒ 0.6; INV-1 ≥ 0.
    [Test]
    public void CalculateDustScore_MixedSequence_ReturnsExactNonNegative()
    {
        double score = SequenceComplexity.CalculateDustScore(new DnaSequence("AATAATAA"));

        Assert.Multiple(() =>
        {
            Assert.That(score, Is.EqualTo(0.6).Within(1e-10),
                "AATAATAA: AAT=ATA=TAA=2 ⇒ Σ=3, divided by ℓ−1 = 5 ⇒ 0.6.");
            Assert.That(score, Is.GreaterThanOrEqualTo(0.0), "INV-1");
        });
    }

    // M8 — reference-discriminating values at the default level-20 threshold (2.0):
    // A×7 ⇒ 10/4 = 2.5 (sdust masks it); A×6 ⇒ 2.0 (not masked); (AC)×6 ⇒ 20/9 (masked);
    // ACACACACACA ⇒ 16/8 = 2.0 (not masked).
    [TestCase("AAAAAAA", 2.5)]
    [TestCase("AAAAAA", 2.0)]
    [TestCase("ACACACACACAC", 20.0 / 9.0)]
    [TestCase("ACACACACACA", 2.0)]
    public void CalculateDustScore_SdustThresholdBoundary_ExactValue(string s, double expected)
    {
        Assert.That(SequenceComplexity.CalculateDustScore(new DnaSequence(s)), Is.EqualTo(expected).Within(1e-10));
    }

    // M9 — exactly one triplet (ℓ = 1): ℓ − 1 = 0 ⇒ no pair exists ⇒ 0 (no division by zero).
    [Test]
    public void CalculateDustScore_SingleTriplet_ReturnsZero()
    {
        Assert.That(SequenceComplexity.CalculateDustScore(new DnaSequence("AAA")), Is.EqualTo(0.0));
    }

    #endregion

    #region CalculateDustScore(string) overload + normalization

    // M6 — DnaSequence and string overloads agree.
    [Test]
    public void CalculateDustScore_OverloadsAgree_SameScore()
    {
        double fromDna = SequenceComplexity.CalculateDustScore(new DnaSequence("AAAAAA"));
        double fromString = SequenceComplexity.CalculateDustScore("AAAAAA");

        Assert.That(fromString, Is.EqualTo(fromDna).Within(1e-10), "INV-4 (both 2.0).");
    }

    // S3 — string overload upper-cases input.
    [Test]
    public void CalculateDustScore_LowercaseString_MatchesUppercase()
    {
        Assert.That(SequenceComplexity.CalculateDustScore("aaaaaa"), Is.EqualTo(2.0).Within(1e-10));
    }

    #endregion

    #region Edge cases and validation

    // S1 — null DnaSequence throws.
    [Test]
    public void CalculateDustScore_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceComplexity.CalculateDustScore((DnaSequence)null!));
    }

    // S2 — null/empty string ⇒ 0.
    [Test]
    public void CalculateDustScore_NullOrEmptyString_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateDustScore((string)null!), Is.EqualTo(0.0));
            Assert.That(SequenceComplexity.CalculateDustScore(""), Is.EqualTo(0.0));
        });
    }

    // C1 — input shorter than wordSize ⇒ 0.
    [Test]
    public void CalculateDustScore_ShorterThanWordSize_ReturnsZero()
    {
        Assert.That(SequenceComplexity.CalculateDustScore("AT", wordSize: 3), Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateDustScore_InvalidWordSize_Throws()
    {
        var seq = new DnaSequence("AAAAAA");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateDustScore(seq, wordSize: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateDustScore("AAAAAA", wordSize: 0));
        });
    }

    #endregion

    #region MaskLowComplexity — symmetric DUST (SDUST), cross-checked against lh3/sdust

    // Each row: input, window W, threshold (= sdust -t T / 10), expected masked string produced
    // by `sdust -w W -t T` (intervals converted to 'N'). Includes non-default W and T.
    private static readonly object[] SdustReferenceCases =
    {
        new object[] { "TAGTATAATAGTATAATAGTATAATAGTATAATAGTATAATAGTATAATAGTATAATAGTAAAAATTGATTGCTTGATTGAT", 64, 2.0,
                       "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAATTGATTGCTTGATTGAT" },
        new object[] { "TGCCATGCCACGCCATGCCATGCCATGCCAAGTCACGTCGTCACGTCGTCACGTCGTCACGTCATCACGTCGTCACGTGGTCACGTCGT", 64, 2.0,
                       "TGCCATGCCACGCCATGCCATGCCATGCCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" },
        new object[] { "ACGAGCTTCTTCTTCTTCTTCTTCTTCTTCTTGAGCACGCGCCATGATCC", 64, 2.0,
                       "ACGAGNNNNNNNNNNNNNNNNNNNNNNNNNNNGAGCACGCGCCATGATCC" },
        new object[] { "AAAAAAAATGTGTGTGTGTGTGTGTGTGGTGTGTCTGTGTGTGTGTGTCCGTCGCACGCACGTACGCACGC", 30, 2.0,
                       "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCGTCGCACGCACGTACGCACGC" },
        new object[] { "CCGGAACCGGAACCGGAACCGGAAACGAGGGTAACGAGCTCTCTCTCTCTCTTGTTCTAAGAAAAAAAAAAAAAAAATGAATGAATGA", 16, 1.5,
                       "CCGGAACCGGAACCGGAACCGGAAACGAGGGTAACGAGNNNNNNNNNNNNNNTGTTCTAAGNNNNNNNNNNNNNNNNTGAATGAATGA" },
        new object[] { "TTTTTTTTATATATATATATATATATATATATATATATAAGTCCACAGGGCGTT", 8, 2.5,
                       "NNNNNNNNATATATATATATATATATATATATATATATAAGTCCACAGGGCGTT" },
        new object[] { "AGAAAAAAAAACTCTCTCTCTCTCTCTCTCTCTCTCTGTCCAAAGAGGTTAGGTTCATACTGTCAAATTAAACCCGCACCCGCACCCGCACCCGCCCCCG", 64, 3.0,
                       "AGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTCCAAAGAGGTTAGGTTCATACTGTCAAATTAAACCCGCACCCGCACCCGCACCCGCCCCCG" },
        new object[] { "AAATCAAATCGTGCGTGCCTGCGTGCGTGCGTGCGTGCGGGGGGGGGGGGGGTGGGGGGGTGGGGGCCCCACAATTGGTGAATTGGTGAATTGGTGAATAG", 5, 1.0,
                       "AAATCAAATCGTGCGTGCCTGCGTGCGTGCGTGCGTGCNNNNNNNNNNNNNNTNNNNNNNTNNNNNCCCCACAATTGGTGAATTGGTGAATTGGTGAATAG" },
        new object[] { "CATTTTGAGGTTTGGTAGATATATGAGTTCCACAGTGTTCCAAGCAATTCTAGACGGATTTCGCG", 64, 2.0,
                       "CATTTTGAGGTTTGGTAGATATATGAGTTCCACAGTGTTCCAAGCAATTCTAGACGGATTTCGCG" },
        // Threshold boundary (default W=64, T=20): run of 7 A (flank A + A×6) masked [19,26); run of 6 A unmasked.
        new object[] { "ACGTTGCAGTCATGCGATCAAAAAAATGCATCGGATCCTAGGCTAA", 64, 2.0,
                       "ACGTTGCAGTCATGCGATCNNNNNNNTGCATCGGATCCTAGGCTAA" },
        new object[] { "ACGTTGCAGTCATGCGATCAAAAAATGCATCGGATCCTAGGCTAA", 64, 2.0,
                       "ACGTTGCAGTCATGCGATCAAAAAATGCATCGGATCCTAGGCTAA" },
        // (AC)×6 masked [12,24); ACACACACACA (score exactly 2.0) not masked.
        new object[] { "GGTTCATGGCTTACACACACACACGGTTCATGGCTT", 64, 2.0,
                       "GGTTCATGGCTTNNNNNNNNNNNNGGTTCATGGCTT" },
        new object[] { "GGTTCATGGCTTACACACACACAGGTTCATGGCTT", 64, 2.0,
                       "GGTTCATGGCTTACACACACACAGGTTCATGGCTT" },
        // Sequence shorter than the window is still scanned (sdust: A×7 ⇒ [0,7)).
        new object[] { "AAAAAAA", 64, 2.0, "NNNNNNN" },
        new object[] { "AAAAAA", 64, 2.0, "AAAAAA" },
    };

    [TestCaseSource(nameof(SdustReferenceCases))]
    public void MaskLowComplexity_MatchesSdustReference(string input, int window, double threshold, string expected)
    {
        string masked = SequenceComplexity.MaskLowComplexity(new DnaSequence(input), window, threshold);

        Assert.That(masked, Is.EqualTo(expected));
    }

    // Symmetry (Morgulis 2006: the SDUST rule is symmetric w.r.t. sequence reversal):
    // masking the reverse equals the reverse of the mask.
    [Test]
    public void MaskLowComplexity_ReversalSymmetric()
    {
        const string s = "CCGGAACCGGAACCGGAACCGGAAACGAGGGTAACGAGCTCTCTCTCTCTCTTGTTCTAAGAAAAAAAAAAAAAAAATGAATGAATGA";
        string rev = new string(s.Reverse().ToArray());

        string maskedRev = SequenceComplexity.MaskLowComplexity(new DnaSequence(rev));
        string revMasked = new string(SequenceComplexity.MaskLowComplexity(new DnaSequence(s)).Reverse().ToArray());

        Assert.That(maskedRev, Is.EqualTo(revMasked));
    }

    [Test]
    public void MaskLowComplexity_InvalidParameters_Throw()
    {
        var seq = new DnaSequence("ACGTACGTAAAAAAAA");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.MaskLowComplexity(seq, windowSize: 2));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.MaskLowComplexity(seq, threshold: -1.0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.MaskLowComplexity(seq, threshold: double.NaN));
            Assert.Throws<ArgumentNullException>(() => SequenceComplexity.MaskLowComplexity((DnaSequence)null!));
        });
    }

    #endregion
}
