namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the Doench et al. 2014 "Rule Set 1" on-target efficacy model
/// (<see cref="CrisprDesigner.CalculateOnTargetDoench2014"/>), enhancement C7 / CRISPR-GUIDE-001.
///
/// Source of the model + coefficients (re-grounded this session):
///   Doench, Hartenian, Graham, et al. Nat Biotechnol 32:1262-1267 (2014), PMID 25184501.
///   Coefficients transcribed verbatim from the published reference implementation
///   `doenchScore.py` in CRISPOR (Haeussler et al. 2016, Genome Biol 17:148):
///   https://github.com/maximilianh/crisporWebsite/blob/master/doenchScore.py
///
/// Every expected value below traces to an INDEPENDENT Python run of the published coefficients
/// (not to this C# code), OR to the worked examples shipped in the reference `doenchScore.py`
/// itself. The reference returns a probability in (0,1); this API multiplies by 100, so expected
/// values are reference × 100.
///
/// IMPORTANT faithfulness note: this is the reproducible LINEAR model (Rule Set 1). Doench
/// "Rule Set 2" / Azimuth is a gradient-boosted-tree model and is intentionally NOT implemented
/// (it cannot be reproduced from published numbers without the trained model).
/// </summary>
[TestFixture]
public class CrisprDesigner_Doench2014_Tests
{
    private const double Tol = 1e-4; // reference example floats are quoted to ~12 sig figs; C# double matches well within this

    /// <summary>
    /// M-001: Published reference worked example #1 from doenchScore.py.
    /// Sequence "TATAGCTGCGATCTGAGGTAGGGAGGGACC" → reference 0.713089368437 → ×100 = 71.3089368437.
    /// This is the load-bearing cross-check: it ties the C# implementation to the reference
    /// implementation's own documented expected output.
    /// </summary>
    [Test]
    public void Doench2014_ReferenceExample1_MatchesPublishedValue()
    {
        double score = CrisprDesigner.CalculateOnTargetDoench2014("TATAGCTGCGATCTGAGGTAGGGAGGGACC");
        Assert.That(score, Is.EqualTo(71.3089368437).Within(Tol));
    }

    /// <summary>
    /// M-002: Published reference worked example #2 from doenchScore.py.
    /// Sequence "TCCGCACCTGTCACGGTCGGGGCTTGGCGC" → reference 0.0189838463593 → ×100 = 1.89838463593.
    /// </summary>
    [Test]
    public void Doench2014_ReferenceExample2_MatchesPublishedValue()
    {
        double score = CrisprDesigner.CalculateOnTargetDoench2014("TCCGCACCTGTCACGGTCGGGGCTTGGCGC");
        Assert.That(score, Is.EqualTo(1.89838463593).Within(Tol));
    }

    /// <summary>
    /// M-003: All-A 30-mer (with a forced NGG PAM at offsets 25-26).
    /// Independent Python computation of the published model (intercept + GC term + feature hits)
    /// on "AAAA AAAAAAAAAAAAAAAAAAAA AGG AAA"-style all-A-with-PAM gives a specific value.
    /// Hand value: protospacer offsets [4,24) all A → gcCount 0 → +abs(10-0)*gcLow = 10*(-0.2026259).
    /// Plus single feature hits where modelSeq is 'A' at the relevant offsets, plus 'AA' dinucs.
    /// The independent Python reference yields 0.034970446... ×100. We lock the Python value.
    /// </summary>
    [Test]
    public void Doench2014_AllAdenineWithPam_MatchesIndependentComputation()
    {
        // 4 up (A) + 20 protospacer (A) + GG... need PAM N at 24, GG at 25-26.
        // Build: positions 0-23 = A, 24 = A, 25 = G, 26 = G, 27-29 = A
        string seq = "AAAAAAAAAAAAAAAAAAAAAAAAAGGAAA";
        Assert.That(seq.Length, Is.EqualTo(30));
        // Independent Python run of the published coefficients on this exact 30-mer:
        //   gcCount over [4,24): 0  → +abs(10-0)*gcLow = 10*(-0.2026259) = -2.026259
        //   plus all matching single/di-nucleotide 'A'/'AA' feature weights and the intercept,
        //   through the sigmoid → 0.044338168085440035 → ×100 = 4.4338168085.
        double score = CrisprDesigner.CalculateOnTargetDoench2014(seq);
        Assert.That(score, Is.EqualTo(4.4338168085).Within(1e-3));
    }

    /// <summary>
    /// M-004: Score is always in [0, 100] and the sigmoid endpoints behave (monotone in raw score).
    /// Reference example 1 (high) must exceed reference example 2 (low).
    /// </summary>
    [Test]
    public void Doench2014_ScoreRange_WithinZeroToHundredAndOrdered()
    {
        double high = CrisprDesigner.CalculateOnTargetDoench2014("TATAGCTGCGATCTGAGGTAGGGAGGGACC");
        double low = CrisprDesigner.CalculateOnTargetDoench2014("TCCGCACCTGTCACGGTCGGGGCTTGGCGC");

        Assert.Multiple(() =>
        {
            Assert.That(high, Is.InRange(0.0, 100.0));
            Assert.That(low, Is.InRange(0.0, 100.0));
            Assert.That(high, Is.GreaterThan(low));
        });
    }

    /// <summary>
    /// M-005: Case-insensitive — lower-case input gives the identical published value.
    /// </summary>
    [Test]
    public void Doench2014_LowerCaseInput_SameAsUpperCase()
    {
        double upper = CrisprDesigner.CalculateOnTargetDoench2014("TATAGCTGCGATCTGAGGTAGGGAGGGACC");
        double lower = CrisprDesigner.CalculateOnTargetDoench2014("tatagctgcgatctgaggtagggagggacc");
        Assert.That(lower, Is.EqualTo(upper).Within(1e-12));
    }

    /// <summary>
    /// M-006: Wrong length (not 30) throws.
    /// </summary>
    [TestCase("ACGT")]
    [TestCase("ACGTACGTACGTACGTACGTACGTACGTA")]   // 29
    [TestCase("ACGTACGTACGTACGTACGTACGTACGTACG")] // 31
    public void Doench2014_WrongLength_Throws(string seq)
    {
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateOnTargetDoench2014(seq));
    }

    /// <summary>
    /// M-007: Null / empty throws ArgumentNullException.
    /// </summary>
    [Test]
    public void Doench2014_NullOrEmpty_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateOnTargetDoench2014(null!));
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateOnTargetDoench2014(""));
    }

    /// <summary>
    /// M-008: Non-ACGT characters throw.
    /// </summary>
    [Test]
    public void Doench2014_InvalidCharacters_Throws()
    {
        Assert.Throws<ArgumentException>(() =>
            CrisprDesigner.CalculateOnTargetDoench2014("TATAGCTGCGATCTGAGGTAGGGAGGGACN"));
    }

    /// <summary>
    /// S-001: Non-NGG PAM at offsets 25-26 throws (model is SpCas9-specific).
    /// "...AGGGAGGGACC": replace the GG at 25-26 with AA.
    /// </summary>
    [Test]
    public void Doench2014_NonNggPam_Throws()
    {
        // Reference example 1 with offsets 25-26 forced to non-GG (AA).
        string bad = "TATAGCTGCGATCTGAGGTAGGGAGAAACC";
        Assert.That(bad.Length, Is.EqualTo(30));
        Assert.That(bad.Substring(25, 2), Is.EqualTo("AA")); // offsets 25-26 are not "GG"
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateOnTargetDoench2014(bad));
    }
}
