// SEQ-ATSKEW-001 — AT Skew
// Evidence: docs/Evidence/SEQ-ATSKEW-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-ATSKEW-001.md
// Source: Charneski CA et al. (2011) PLoS Genet 7(9):e1002283; Lobry JR (1996) Mol Biol Evol 13(5):660-665.

using System;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class GcSkewCalculator_CalculateAtSkew_Tests
{
    #region CalculateAtSkew(string)

    // M1 — Bounds: T=0 => +1 (Wikipedia/Lobry range). A wrong sign or off-by-one would not give exactly 1.
    [Test]
    public void CalculateAtSkew_PureAdenine_ReturnsPlusOne()
    {
        double skew = GcSkewCalculator.CalculateAtSkew("AAAA");

        Assert.That(skew, Is.EqualTo(1.0).Within(1e-10),
            "Pure-A sequence has T=0, so (A-T)/(A+T) = 4/4 = +1 (upper bound)");
    }

    // M2 — Bounds: A=0 => -1 (Wikipedia/Lobry range).
    [Test]
    public void CalculateAtSkew_PureThymine_ReturnsMinusOne()
    {
        double skew = GcSkewCalculator.CalculateAtSkew("TTTT");

        Assert.That(skew, Is.EqualTo(-1.0).Within(1e-10),
            "Pure-T sequence has A=0, so (A-T)/(A+T) = -4/4 = -1 (lower bound)");
    }

    // M3 — Balanced A=T => 0 (formula numerator A-T=0; Charneski 2011).
    [Test]
    public void CalculateAtSkew_BalancedAt_ReturnsZero()
    {
        double skew = GcSkewCalculator.CalculateAtSkew("ATAT");

        Assert.That(skew, Is.EqualTo(0.0).Within(1e-10),
            "A=2,T=2 => (2-2)/4 = 0; balanced strand has zero AT skew");
    }

    // M4 — Asymmetric positive: A=3,T=1 => (3-1)/4 = 0.5 (Charneski 2011 formula).
    [Test]
    public void CalculateAtSkew_ExcessAdenine_ReturnsExactFraction()
    {
        double skew = GcSkewCalculator.CalculateAtSkew("AAAT");

        Assert.That(skew, Is.EqualTo(0.5).Within(1e-10),
            "A=3,T=1 => (3-1)/(3+1) = 0.5; exact fraction from (A-T)/(A+T)");
    }

    // M5 — Asymmetric negative: A=1,T=3 => (1-3)/4 = -0.5 (Charneski 2011 formula).
    [Test]
    public void CalculateAtSkew_ExcessThymine_ReturnsNegativeFraction()
    {
        double skew = GcSkewCalculator.CalculateAtSkew("ATTT");

        Assert.That(skew, Is.EqualTo(-0.5).Within(1e-10),
            "A=1,T=3 => (1-3)/(1+3) = -0.5; exact negative fraction");
    }

    // M6 — Zero denominator: no A and no T => 0 (Biopython ZeroDivisionError -> 0.0).
    [Test]
    public void CalculateAtSkew_NoAdenineOrThymine_ReturnsZero()
    {
        double skew = GcSkewCalculator.CalculateAtSkew("GGCC");

        Assert.That(skew, Is.EqualTo(0.0).Within(1e-10),
            "A+T = 0 => result is 0 (no exception/NaN), per Biopython convention");
    }

    // M7 — Only A/T counted; G/C and others ignored (Biopython "does NOT look at ambiguous nucleotides").
    [Test]
    public void CalculateAtSkew_IgnoresGcAndOtherSymbols()
    {
        // A=3, T=1; the six G/C bases must NOT enter the denominator.
        double skew = GcSkewCalculator.CalculateAtSkew("AAATGGGCCC");

        Assert.That(skew, Is.EqualTo(0.5).Within(1e-10),
            "G/C ignored: (3-1)/(3+1) = 0.5; counting them would give a smaller magnitude");
    }

    // S1 — Case-insensitive counting (Biopython; repo ToUpperInvariant).
    [Test]
    public void CalculateAtSkew_LowercaseInput_MatchesUppercase()
    {
        double skew = GcSkewCalculator.CalculateAtSkew("aaat");

        Assert.That(skew, Is.EqualTo(0.5).Within(1e-10),
            "Lowercase 'aaat' counts the same as 'AAAT' => 0.5 (case-insensitive)");
    }

    // S2 — Null string => 0 (documented validation).
    [Test]
    public void CalculateAtSkew_NullString_ReturnsZero()
    {
        double skew = GcSkewCalculator.CalculateAtSkew((string?)null!);

        Assert.That(skew, Is.EqualTo(0.0).Within(1e-10),
            "Null string is treated as no A/T => 0");
    }

    // S3 — Empty string => 0 (documented validation).
    [Test]
    public void CalculateAtSkew_EmptyString_ReturnsZero()
    {
        double skew = GcSkewCalculator.CalculateAtSkew("");

        Assert.That(skew, Is.EqualTo(0.0).Within(1e-10),
            "Empty string has A+T=0 => 0");
    }

    // S5 / INV-1 — mixed sequence: exact sourced value AND within [-1, 1].
    // "AAATTGCGCAATA" => A=6, T=3, G/C ignored => (6-3)/(6+3) = 1/3 (formula (A-T)/(A+T),
    // Charneski 2011); the result is also within the [-1,+1] range bound (Wikipedia/Lobry).
    [Test]
    public void CalculateAtSkew_MixedSequence_ReturnsExactValueWithinBounds()
    {
        double skew = GcSkewCalculator.CalculateAtSkew("AAATTGCGCAATA");

        Assert.Multiple(() =>
        {
            Assert.That(skew, Is.EqualTo(1.0 / 3.0).Within(1e-10),
                "A=6, T=3 => (6-3)/(6+3) = 1/3; G/C ignored (exact value from (A-T)/(A+T))");
            Assert.That(skew, Is.GreaterThanOrEqualTo(-1.0),
                "AT skew is bounded below by -1 (INV-01)");
            Assert.That(skew, Is.LessThanOrEqualTo(1.0),
                "AT skew is bounded above by +1 (INV-01)");
        });
    }

    #endregion

    #region CalculateAtSkew(DnaSequence)

    // S4 — Null DnaSequence => ArgumentNullException (documented validation).
    [Test]
    public void CalculateAtSkew_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => GcSkewCalculator.CalculateAtSkew((DnaSequence)null!),
            "Null DnaSequence must throw ArgumentNullException");
    }

    // C1 — Delegate equivalence: DnaSequence overload == string overload on the same data.
    [Test]
    public void CalculateAtSkew_DnaSequenceOverload_MatchesStringOverload()
    {
        const string seq = "AAATGGGCCC";
        double viaString = GcSkewCalculator.CalculateAtSkew(seq);
        double viaDna = GcSkewCalculator.CalculateAtSkew(new DnaSequence(seq));

        Assert.Multiple(() =>
        {
            // Exact sourced value: A=3, T=1, G/C ignored => (3-1)/(3+1) = 0.5.
            Assert.That(viaDna, Is.EqualTo(0.5).Within(1e-10),
                "DnaSequence overload computes (A-T)/(A+T) = (3-1)/4 = 0.5");
            Assert.That(viaDna, Is.EqualTo(viaString).Within(1e-10),
                "DnaSequence overload delegates to the same core => identical to string overload");
        });
    }

    #endregion
}
