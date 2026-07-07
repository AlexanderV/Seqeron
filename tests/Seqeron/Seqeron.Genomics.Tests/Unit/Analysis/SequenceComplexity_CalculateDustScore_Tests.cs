// SEQ-COMPLEX-DUST-001 — DUST Score
// Evidence: docs/Evidence/SEQ-COMPLEX-DUST-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-COMPLEX-DUST-001.md
// Source: Morgulis A, Gertz EM, Schäffer AA, Agarwala R (2006). J Comput Biol 13(5):1028–1040,
//         doi:10.1089/cmb.2006.13.1028; Li H (2025) longdust arXiv:2509.07357; lh3/sdust.
//
// Spec: S(x) = (Σ_t c_t·(c_t−1)/2) / (L−2) over overlapping triplets (k=3). Higher ⇒ lower complexity.
// Expected values are derived independently from the formula by hand (see Evidence §Test Datasets),
// NOT from the implementation output. A wrong divisor (e.g. L−3) would FAIL these tests.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceComplexity_CalculateDustScore_Tests
{
    #region CalculateDustScore(DnaSequence, k) — canonical exact values

    // M1 — AAAAAA (L=6): AAA=4 ⇒ Σ=4·3/2=6, /(L−2)=4 ⇒ 1.5. A divisor of L−3=3 would give 2.0 (wrong).
    [Test]
    public void CalculateDustScore_Homopolymer6_Returns1Point5()
    {
        var seq = new DnaSequence("AAAAAA");

        double score = SequenceComplexity.CalculateDustScore(seq);

        Assert.That(score, Is.EqualTo(1.5).Within(1e-10),
            "AAAAAA has one triplet AAA with count 4: Σ c(c−1)/2 = 6, divided by L−2 = 4 ⇒ 1.5.");
    }

    // M2 — ACGTACGT (L=8): ACG=2,CGT=2,GTA=1,TAC=1 ⇒ Σ=1+1=2, /(L−2)=6 ⇒ 1/3.
    [Test]
    public void CalculateDustScore_RepeatedTetramer_ReturnsOneThird()
    {
        var seq = new DnaSequence("ACGTACGT");

        double score = SequenceComplexity.CalculateDustScore(seq);

        Assert.That(score, Is.EqualTo(1.0 / 3.0).Within(1e-10),
            "ACGTACGT: ACG=2,CGT=2 each contribute 1, GTA=TAC=1 contribute 0 ⇒ Σ=2, /6 = 0.3333…");
    }

    // M3 — ATGC (L=4): ATG=1,TGC=1 ⇒ Σ=0 ⇒ score 0 (all-distinct triplets, INV-2).
    [Test]
    public void CalculateDustScore_AllDistinctTriplets_ReturnsZero()
    {
        var seq = new DnaSequence("ATGC");

        double score = SequenceComplexity.CalculateDustScore(seq);

        Assert.That(score, Is.EqualTo(0.0).Within(1e-10),
            "ATGC has two distinct triplets each count 1; c(c−1)/2 = 0 for both ⇒ score 0 (max complexity).");
    }

    // M4 — ACACACAC (L=8): ACA=3,CAC=3 ⇒ Σ=3+3=6, /(L−2)=6 ⇒ 1.0.
    [Test]
    public void CalculateDustScore_DinucleotideRepeat_ReturnsOne()
    {
        var seq = new DnaSequence("ACACACAC");

        double score = SequenceComplexity.CalculateDustScore(seq);

        Assert.That(score, Is.EqualTo(1.0).Within(1e-10),
            "ACACACAC: ACA=3,CAC=3 ⇒ each 3·2/2=3, Σ=6, divided by L−2=6 ⇒ 1.0.");
    }

    // M5 — AAAAAAAAAA (L=10): AAA=8 ⇒ Σ=8·7/2=28, /(L−2)=8 ⇒ 3.5 (homopolymer (L−3)/2, INV-5).
    [Test]
    public void CalculateDustScore_Homopolymer10_Returns3Point5()
    {
        var seq = new DnaSequence("AAAAAAAAAA");

        double score = SequenceComplexity.CalculateDustScore(seq);

        Assert.That(score, Is.EqualTo(3.5).Within(1e-10),
            "AAAAAAAAAA: AAA count 8 ⇒ Σ=28, divided by L−2=8 ⇒ 3.5 = (L−3)/2.");
    }

    // M7 — INV-1: score ≥ 0, with an exact known value for a mixed sequence.
    // AATAATAA (L=8): AAT=2,ATA=2,TAA=2 ⇒ Σ=1+1+1=3, /(L−2)=6 ⇒ 0.5.
    [Test]
    public void CalculateDustScore_MixedSequence_ReturnsExactNonNegative()
    {
        var seq = new DnaSequence("AATAATAA");

        double score = SequenceComplexity.CalculateDustScore(seq);

        Assert.Multiple(() =>
        {
            Assert.That(score, Is.EqualTo(0.5).Within(1e-10),
                "AATAATAA: AAT=2,ATA=2,TAA=2 ⇒ Σ=3, divided by L−2=6 ⇒ 0.5.");
            Assert.That(score, Is.GreaterThanOrEqualTo(0.0),
                "INV-1: DUST score is a sum of non-negative terms over a positive divisor, so ≥ 0.");
        });
    }

    #endregion

    #region CalculateDustScore(string) overload + normalization

    // M6 — DnaSequence and string overloads agree on the same sequence.
    [Test]
    public void CalculateDustScore_OverloadsAgree_SameScore()
    {
        var seq = new DnaSequence("AAAAAA");

        double fromDna = SequenceComplexity.CalculateDustScore(seq);
        double fromString = SequenceComplexity.CalculateDustScore("AAAAAA");

        Assert.That(fromString, Is.EqualTo(fromDna).Within(1e-10),
            "INV-4: both overloads wrap one core; identical input ⇒ identical score (1.5).");
    }

    // S3 — string overload upper-cases input, so lowercase yields the same score as upper.
    [Test]
    public void CalculateDustScore_LowercaseString_MatchesUppercase()
    {
        double lower = SequenceComplexity.CalculateDustScore("aaaaaa");

        Assert.That(lower, Is.EqualTo(1.5).Within(1e-10),
            "string overload upper-cases input; 'aaaaaa' ⇒ same triplet AAA count 4 ⇒ 1.5.");
    }

    #endregion

    #region Edge cases and validation

    // S1 — null DnaSequence throws.
    [Test]
    public void CalculateDustScore_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => SequenceComplexity.CalculateDustScore((DnaSequence)null!),
            "null DnaSequence must raise ArgumentNullException per the validation contract.");
    }

    // S2 — null/empty string ⇒ 0.
    [Test]
    public void CalculateDustScore_NullOrEmptyString_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateDustScore((string)null!), Is.EqualTo(0.0).Within(1e-10),
                "null string yields 0: no words to score.");
            Assert.That(SequenceComplexity.CalculateDustScore(""), Is.EqualTo(0.0).Within(1e-10),
                "empty string yields 0: no words to score.");
        });
    }

    // C1 — input shorter than wordSize ⇒ 0 (no triplet exists).
    [Test]
    public void CalculateDustScore_ShorterThanWordSize_ReturnsZero()
    {
        double score = SequenceComplexity.CalculateDustScore("AT", wordSize: 3);

        Assert.That(score, Is.EqualTo(0.0).Within(1e-10),
            "'AT' is shorter than the triplet word; no word exists ⇒ score 0 (defined convention).");
    }

    // Validation — wordSize < 1 throws on both overloads.
    [Test]
    public void CalculateDustScore_InvalidWordSize_Throws()
    {
        var seq = new DnaSequence("AAAAAA");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => SequenceComplexity.CalculateDustScore(seq, wordSize: 0),
                "wordSize 0 is invalid (no word) ⇒ ArgumentOutOfRangeException.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => SequenceComplexity.CalculateDustScore("AAAAAA", wordSize: 0),
                "wordSize 0 is invalid on the string overload ⇒ ArgumentOutOfRangeException.");
        });
    }

    #endregion
}
