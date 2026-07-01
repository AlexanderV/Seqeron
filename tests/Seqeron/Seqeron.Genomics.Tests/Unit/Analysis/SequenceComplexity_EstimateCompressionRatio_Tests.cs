// SEQ-COMPLEX-COMPRESS-001 — Compression Ratio (Lempel–Ziv complexity)
// Evidence: docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-COMPLEX-COMPRESS-001.md
// Source: Lempel A, Ziv J (1976). On the Complexity of Finite Sequences.
//         IEEE Trans. Inf. Theory 22(1):75–81, doi:10.1109/TIT.1976.1055501.
//         Reference impl. Naereen/Lempel-Ziv_Complexity (lempel_ziv_complexity.py);
//         normalization per entropy/antropy lziv_complexity + Zhang et al. (2009).
//
// Spec: c(S) = number of components from an exhaustive-history left-to-right parse.
// Worked values are the reference implementation's doctests (8, 7, 9, 10), derived
// independently of this code. A wrong parsing convention (e.g. the entropy variant
// giving 6 for the first string) would FAIL these tests.

using NUnit.Framework;
using System;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceComplexity_EstimateCompressionRatio_Tests
{
    #region CalculateLempelZivComplexity — canonical exact values

    // M1 — '1001111011000010' → 1/0/01/11/10/110/00/010 = 8 components (Naereen doctest).
    [Test]
    public void CalculateLempelZivComplexity_Doctest1_Returns8()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("1001111011000010");

        Assert.That(c, Is.EqualTo(8),
            "Exhaustive-history parse yields components 1/0/01/11/10/110/00/010 = 8 (Lempel-Ziv 1976; Naereen doctest). The entropy-convention value 6 would be wrong.");
    }

    // M2 — '1010101010101010' → 1/0/10/101/01/010/1010 = 7 components (Naereen doctest).
    [Test]
    public void CalculateLempelZivComplexity_Doctest2_Returns7()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("1010101010101010");

        Assert.That(c, Is.EqualTo(7),
            "Components 1/0/10/101/01/010/1010 = 7 (Naereen doctest).");
    }

    // M3 — '1001111011000010000010' → adds 000 → 9 components (Naereen doctest).
    [Test]
    public void CalculateLempelZivComplexity_Doctest3_Returns9()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("1001111011000010000010");

        Assert.That(c, Is.EqualTo(9),
            "Components 1/0/01/11/10/110/00/010/000 = 9 (Naereen doctest).");
    }

    // M4 — '100111101100001000001010' → adds 0101 → 10 components (Naereen doctest).
    [Test]
    public void CalculateLempelZivComplexity_Doctest4_Returns10()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("100111101100001000001010");

        Assert.That(c, Is.EqualTo(10),
            "Components …/000/0101 = 10 (Naereen doctest).");
    }

    // M5 — homopolymer '0'×16 → 0/00/000/0000/00000 = 5 (NOT 1; traced parser).
    [Test]
    public void CalculateLempelZivComplexity_Homopolymer16_Returns5()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("0000000000000000");

        Assert.That(c, Is.EqualTo(5),
            "A length-16 homopolymer parses as 0/00/000/0000/00000 = 5 components (INV-04 low complexity); a trivial 'returns 1' impl would fail.");
    }

    // M6 — 'ACGT' all-distinct symbols → each its own component = 4.
    [Test]
    public void CalculateLempelZivComplexity_AllDistinct_Returns4()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("ACGT");

        Assert.That(c, Is.EqualTo(4),
            "Four distinct symbols each form a new component A/C/G/T = 4 (max complexity for n=4).");
    }

    #endregion

    #region CalculateNormalizedLempelZivComplexity — exact values

    // M7 — '1001111011000010': n=16,b=2,c=8 → 8/(16/log2(16))=8/4=2.0.
    [Test]
    public void CalculateNormalizedLempelZivComplexity_Doctest1_Returns2()
    {
        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity("1001111011000010");

        Assert.That(norm, Is.EqualTo(2.0).Within(1e-10),
            "Normalized = c/(n/log_b n) = 8/(16/log2 16) = 8/(16/4) = 2.0 (entropy/antropy formula).");
    }

    // M8 — single-symbol input: the entropy/antropy reference clamps the log base to 2
    // (`base = 2 if base < 2 else base`), so normalized = c/(n/log2 n) = 5/(16/log2 16)
    // = 5/(16/4) = 1.25. (NOT the raw count; verified against antropy entropy.py.)
    [Test]
    public void CalculateNormalizedLempelZivComplexity_SingleSymbol_ClampsBaseToTwo()
    {
        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity("0000000000000000");

        Assert.That(norm, Is.EqualTo(1.25).Within(1e-10),
            "With one distinct symbol the reference clamps base to 2: 5/(16/log2 16) = 5/(16/4) = 1.25 (antropy lziv_complexity).");
    }

    #endregion

    #region EstimateCompressionRatio — delegation (smoke)

    // M9 — EstimateCompressionRatio delegates to normalized LZ (2.0 for doctest 1).
    [Test]
    public void EstimateCompressionRatio_Doctest1_EqualsNormalized()
    {
        double ratio = SequenceComplexity.EstimateCompressionRatio("1001111011000010");
        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity("1001111011000010");

        Assert.Multiple(() =>
        {
            Assert.That(ratio, Is.EqualTo(2.0).Within(1e-10),
                "EstimateCompressionRatio returns the normalized LZ complexity = 2.0 (INV-05).");
            Assert.That(ratio, Is.EqualTo(norm).Within(1e-10),
                "EstimateCompressionRatio must equal CalculateNormalizedLempelZivComplexity (delegation).");
        });
    }

    #endregion

    #region Edge cases

    // S1 — empty string → 0 components (INV-01).
    [Test]
    public void CalculateLempelZivComplexity_EmptyString_ReturnsZero()
    {
        Assert.That(SequenceComplexity.CalculateLempelZivComplexity(""), Is.EqualTo(0),
            "Empty input produces no components (INV-01).");
    }

    // S2 — null DnaSequence → ArgumentNullException.
    [Test]
    public void CalculateLempelZivComplexity_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => SequenceComplexity.CalculateLempelZivComplexity((DnaSequence)null!),
            "Null DnaSequence must throw ArgumentNullException (sibling-method convention).");
    }

    // S3 — single base → 1 component (INV-02).
    [Test]
    public void CalculateLempelZivComplexity_SingleBase_ReturnsOne()
    {
        Assert.That(SequenceComplexity.CalculateLempelZivComplexity("A"), Is.EqualTo(1),
            "A single symbol is one component (INV-02: c≥1 for non-empty).");
    }

    // S4 — DnaSequence overload parity: 'ACGT' via DnaSequence = 4.
    [Test]
    public void CalculateLempelZivComplexity_DnaSequenceOverload_MatchesString()
    {
        int viaSeq = SequenceComplexity.CalculateLempelZivComplexity(new DnaSequence("ACGT"));

        Assert.That(viaSeq, Is.EqualTo(4),
            "DnaSequence overload must agree with the string overload for ACGT (=4).");
    }

    #endregion

    #region Invariants / property

    // C1 — INV-04: homopolymer is strictly less complex than all-distinct of same length.
    [Test]
    public void CalculateLempelZivComplexity_HomopolymerLessThanAllDistinct()
    {
        int homo = SequenceComplexity.CalculateLempelZivComplexity("AAAA");
        int distinct = SequenceComplexity.CalculateLempelZivComplexity("ACGT");

        Assert.Multiple(() =>
        {
            Assert.That(homo, Is.EqualTo(2),
                "AAAA parses as A/AA = 2 components.");
            Assert.That(homo, Is.LessThan(distinct),
                "INV-04: a homopolymer (2) is strictly less complex than all-distinct ACGT (4) of the same length.");
        });
    }

    // C2 — DNA alphabet (b=4) normalization uses log base 4.
    // 'ACGTACGTACGTACGT': n=16,b=4,c=9 → 9/(16/log4 16)=9/(16/2)=9/8=1.125.
    [Test]
    public void CalculateNormalizedLempelZivComplexity_DnaAlphabet_UsesBase4()
    {
        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity("ACGTACGTACGTACGT");

        Assert.That(norm, Is.EqualTo(1.125).Within(1e-10),
            "n=16,b=4,c=9: normalized = 9/(16/log4 16) = 9/(16/2) = 9/8 = 1.125 (Zhang 2009 bio-sequence application).");
    }

    // Normalized DnaSequence overload parity: same input via DnaSequence must equal the
    // string overload (b=4 → 1.125). Exercises CalculateNormalizedLempelZivComplexity(DnaSequence).
    [Test]
    public void CalculateNormalizedLempelZivComplexity_DnaSequenceOverload_MatchesString()
    {
        double viaSeq = SequenceComplexity.CalculateNormalizedLempelZivComplexity(new DnaSequence("ACGTACGTACGTACGT"));
        double viaStr = SequenceComplexity.CalculateNormalizedLempelZivComplexity("ACGTACGTACGTACGT");

        Assert.Multiple(() =>
        {
            Assert.That(viaSeq, Is.EqualTo(1.125).Within(1e-10),
                "DnaSequence overload must give 1.125 for ACGT×4 (b=4).");
            Assert.That(viaSeq, Is.EqualTo(viaStr).Within(1e-10),
                "DnaSequence and string normalized overloads must agree.");
        });
    }

    // Normalized null DnaSequence → ArgumentNullException (overload guard parity).
    [Test]
    public void CalculateNormalizedLempelZivComplexity_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => SequenceComplexity.CalculateNormalizedLempelZivComplexity((DnaSequence)null!),
            "Null DnaSequence must throw ArgumentNullException.");
    }

    // Normalized empty string → 0 (length guard).
    [Test]
    public void CalculateNormalizedLempelZivComplexity_EmptyString_ReturnsZero()
    {
        Assert.That(SequenceComplexity.CalculateNormalizedLempelZivComplexity(""), Is.EqualTo(0.0),
            "Empty input → 0 (no components, length guard).");
    }

    // Degenerate single-character input: n=1 ⇒ log_b(1)=0 ⇒ raw count returned (=1),
    // avoiding division by zero. antropy would yield lz/(1/0)=0; we adopt the raw-count guard.
    [Test]
    public void CalculateNormalizedLempelZivComplexity_SingleChar_ReturnsRawCount()
    {
        Assert.That(SequenceComplexity.CalculateNormalizedLempelZivComplexity("A"), Is.EqualTo(1.0).Within(1e-10),
            "n=1: log_b(1)=0; the degenerate guard returns the raw count (1).");
    }

    // EstimateCompressionRatio DnaSequence overload delegates to the normalized value.
    [Test]
    public void EstimateCompressionRatio_DnaSequenceOverload_EqualsNormalized()
    {
        var seq = new DnaSequence("ACGTACGTACGTACGT");
        double ratio = SequenceComplexity.EstimateCompressionRatio(seq);
        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity(seq);

        Assert.Multiple(() =>
        {
            Assert.That(ratio, Is.EqualTo(1.125).Within(1e-10),
                "EstimateCompressionRatio(DnaSequence) returns normalized LZ = 1.125 for ACGT×4.");
            Assert.That(ratio, Is.EqualTo(norm).Within(1e-10),
                "EstimateCompressionRatio(DnaSequence) must equal CalculateNormalizedLempelZivComplexity (delegation).");
        });
    }

    #endregion
}
