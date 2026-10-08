// SEQ-COMPLEX-COMPRESS-001 — Compression Ratio (Lempel–Ziv complexity)
// Evidence: docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-COMPLEX-COMPRESS-001.md
// Source: Lempel A, Ziv J (1976). On the Complexity of Finite Sequences.
//         IEEE Trans. Inf. Theory 22(1):75–81, doi:10.1109/TIT.1976.1055501.
//         Kaspar F, Schuster HG (1987) Phys. Rev. A 36:842 (linear scan);
//         reference impl. antropy.lziv_complexity (doctests + cross-check);
//         normalization per Zhang et al. (2009).
//
// Spec: c(S) = number of components of the LZ76 exhaustive history of S.
// Expected values are the antropy doctests / published worked examples, derived
// independently of this code. The LZ78 incremental-parse values (8 for
// '1001111011000010', 5 for '0'×16) are WRONG for LZ76 and would FAIL these tests.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceComplexity_EstimateCompressionRatio_Tests
{
    #region CalculateLempelZivComplexity — canonical exact values

    // M1 — '1001111011000010' → 1/0/01/1110/1100/0010 = 6 (antropy doctest; Wikipedia LZ76 example).
    [Test]
    public void CalculateLempelZivComplexity_Doctest1_Returns6()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("1001111011000010");

        Assert.That(c, Is.EqualTo(6),
            "LZ76 exhaustive history 1/0/01/1110/1100/0010 = 6 (antropy doctest). The LZ78 incremental-parse value 8 would be wrong.");
    }

    // M2 — Lempel & Ziv (1976) paper example: 0001101001000101 = 0·001·10·100·1000·101 → 6.
    [Test]
    public void CalculateLempelZivComplexity_LempelZiv1976Example_Returns6()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("0001101001000101");

        Assert.That(c, Is.EqualTo(6),
            "Exhaustive history 0·001·10·100·1000·101 = 6 (Lempel & Ziv 1976 worked example).");
    }

    // M3 — Estévez-Rams et al. (arXiv:1311.0546) example: 010011101101100 = 0.1.00.11.101.101100 → 6.
    [Test]
    public void CalculateLempelZivComplexity_EstevezRamsExample_Returns6()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("010011101101100");

        Assert.That(c, Is.EqualTo(6),
            "Exhaustive history 0.1.00.11.101.101100 = 6; the trailing reproducible factor 101100 is counted.");
    }

    // M4 — period-2 string: 1 / 0 / 10101010101010 (self-overlapping copy) → 3.
    [Test]
    public void CalculateLempelZivComplexity_Period2_Returns3()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("1010101010101010");

        Assert.That(c, Is.EqualTo(3),
            "1/0/10101010101010 = 3 (antropy doctest pattern '1 / 0 / 10'; overlap copy allowed).");
    }

    // M5 — homopolymer '0'×16 → 0 / 000000000000000 (self-overlapping copy) = 2.
    [Test]
    public void CalculateLempelZivComplexity_Homopolymer16_Returns2()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("0000000000000000");

        Assert.That(c, Is.EqualTo(2),
            "A homopolymer's exhaustive history is 0 / 0…0 = 2 components (LZ76; antropy). LZ78 would give 5.");
    }

    // M5b — non-DNA text via the lenient string surface (antropy doctests).
    [TestCase("HELLO WORLD! HELLO WORLD! HELLO WORLD! HELLO WORLD!", 11, 0.38596001132145313)]
    [TestCase("ABCDEFGHIJKLMNOPQRSTUVWXYZ", 26, 1.0)]
    public void LempelZiv_AntropyTextDoctests(string s, int expectedRaw, double expectedNorm)
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateLempelZivComplexity(s), Is.EqualTo(expectedRaw),
                "antropy lziv_complexity doctest raw value");
            Assert.That(SequenceComplexity.CalculateNormalizedLempelZivComplexity(s), Is.EqualTo(expectedNorm).Within(1e-12),
                "antropy lziv_complexity(normalize=True) doctest value");
        });
    }

    // M6 — 'ACGT' all-distinct symbols → each its own component = 4.
    [Test]
    public void CalculateLempelZivComplexity_AllDistinct_Returns4()
    {
        int c = SequenceComplexity.CalculateLempelZivComplexity("ACGT");

        Assert.That(c, Is.EqualTo(4),
            "Four distinct symbols each form a new component A/C/G/T = 4 (max complexity for n=4).");
    }

    // R1 — reference dataset: raw + normalized values computed with antropy 0.2.2
    // (entropy._lz_complexity, Kaspar–Schuster LZ76; lziv_complexity(normalize=True)) on
    // random DNA (python random.seed(2026)) and on a CAG-repeat (low-complexity) sequence.
    // The same values were confirmed by a brute-force exhaustive-history parse.
    [TestCase("AGACTTTCAAAGATATGCTG", 9, 0.9724338213496565)]
    [TestCase("GGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGTGTTTCGGAA", 16, 0.9030169903639559)]
    [TestCase("CTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAATCTACCCCC", 30, 0.9965784284662087)]
    [TestCase("TGTTATGCGCGTTTGTCGTTAGACCAATGTCAGCGCAGCGGCAGATCAAGCAGGAGGCGGAATGTAAACAGAAGGTATGCTTAGGTGGATAGGGAGTGAGCAACAAACGGATCGTTTCTCCCATGCCAAGTTGGCACAGGGAACTACCTGCGGCGGTTTGCCTCTAGTACAGGGCAACGATTCAACTGGGACCGGGGCTC", 55, 1.0510302260940245)]
    [TestCase("ATGCGATCGATTAGCCGATAGCTAGGCTAACGTTAGCATG", 16, 1.0643856189774725)]
    [TestCase("CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG", 4, 0.2484555351907228)]
    public void LempelZiv_DnaReferenceDataset_MatchesAntropy(string s, int expectedRaw, double expectedNorm)
    {
        var dna = new DnaSequence(s);
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateLempelZivComplexity(dna), Is.EqualTo(expectedRaw),
                "raw LZ76 complexity must equal antropy _lz_complexity");
            Assert.That(SequenceComplexity.EstimateCompressionRatio(dna), Is.EqualTo(expectedNorm).Within(1e-12),
                "normalized LZ must equal antropy lziv_complexity(normalize=True)");
        });
    }

    // R2 — long-sequence reference: 20 000 bases from the LCG x ← (1103515245·x + 12345) mod 2³¹
    // (seed 12345), base = "ACGT"[(x >> 16) & 3]. antropy 0.2.2 gives c = 2756 and
    // lziv_complexity(normalize=True) = 0.984423382950957 for this exact string.
    [Test]
    public void LempelZiv_LongLcgDna_MatchesAntropy()
    {
        var chars = new char[20000];
        long x = 12345;
        for (int i = 0; i < chars.Length; i++)
        {
            x = (x * 1103515245 + 12345) % (1L << 31);
            chars[i] = "ACGT"[(int)((x >> 16) & 3)];
        }
        string s = new(chars);

        Assert.Multiple(() =>
        {
            Assert.That(s[..20], Is.EqualTo("AACGTCCGGCATGTTACACA"), "generator must reproduce the Python string");
            Assert.That(SequenceComplexity.CalculateLempelZivComplexity(s), Is.EqualTo(2756),
                "raw LZ76 complexity = antropy _lz_complexity");
            Assert.That(SequenceComplexity.CalculateNormalizedLempelZivComplexity(s), Is.EqualTo(0.984423382950957).Within(1e-12),
                "normalized LZ = antropy lziv_complexity(normalize=True)");
        });
    }

    #endregion

    #region CalculateNormalizedLempelZivComplexity — exact values

    // M7 — '1001111011000010': n=16,b=2,c=6 → 6/(16/log2(16))=6/4=1.5 (antropy doctest).
    [Test]
    public void CalculateNormalizedLempelZivComplexity_Doctest1_Returns1_5()
    {
        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity("1001111011000010");

        Assert.That(norm, Is.EqualTo(1.5).Within(1e-10),
            "Normalized = c/(n/log_b n) = 6/(16/log2 16) = 6/4 = 1.5 (antropy doctest).");
    }

    // M8 — single-symbol input: the entropy/antropy reference clamps the log base to 2
    // (`base = 2 if base < 2 else base`), so normalized = c/(n/log2 n) = 2/(16/log2 16)
    // = 2/(16/4) = 0.5. (NOT the raw count; verified against antropy entropy.py.)
    [Test]
    public void CalculateNormalizedLempelZivComplexity_SingleSymbol_ClampsBaseToTwo()
    {
        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity("0000000000000000");

        Assert.That(norm, Is.EqualTo(0.5).Within(1e-10),
            "With one distinct symbol the reference clamps base to 2: 2/(16/log2 16) = 2/(16/4) = 0.5 (antropy lziv_complexity).");
    }

    #endregion

    #region EstimateCompressionRatio — delegation (smoke)

    // M9 — EstimateCompressionRatio delegates to normalized LZ (1.5 for doctest 1).
    [Test]
    public void EstimateCompressionRatio_Doctest1_EqualsNormalized()
    {
        double ratio = SequenceComplexity.EstimateCompressionRatio("1001111011000010");
        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity("1001111011000010");

        Assert.Multiple(() =>
        {
            Assert.That(ratio, Is.EqualTo(1.5).Within(1e-10),
                "EstimateCompressionRatio returns the normalized LZ complexity = 1.5 (INV-05).");
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
                "AAAA parses as A/AAA = 2 components (LZ76).");
            Assert.That(homo, Is.LessThan(distinct),
                "INV-04: a homopolymer (2) is strictly less complex than all-distinct ACGT (4) of the same length.");
        });
    }

    // C2 — DNA alphabet (b=4) normalization uses log base 4.
    // 'ACGTACGTACGTACGT' = A/C/G/T/ACGTACGTACGT: n=16,b=4,c=5 → 5/(16/log4 16)=5/8=0.625.
    [Test]
    public void CalculateNormalizedLempelZivComplexity_DnaAlphabet_UsesBase4()
    {
        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity("ACGTACGTACGTACGT");

        Assert.That(norm, Is.EqualTo(0.625).Within(1e-10),
            "n=16,b=4,c=5: normalized = 5/(16/log4 16) = 5/(16/2) = 0.625 (antropy cross-check).");
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
            Assert.That(viaSeq, Is.EqualTo(0.625).Within(1e-10),
                "DnaSequence overload must give 0.625 for ACGT×4 (b=4).");
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
    // avoiding division by zero (antropy raises ZeroDivisionError here); we adopt the raw-count guard.
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
            Assert.That(ratio, Is.EqualTo(0.625).Within(1e-10),
                "EstimateCompressionRatio(DnaSequence) returns normalized LZ = 0.625 for ACGT×4.");
            Assert.That(ratio, Is.EqualTo(norm).Within(1e-10),
                "EstimateCompressionRatio(DnaSequence) must equal CalculateNormalizedLempelZivComplexity (delegation).");
        });
    }

    #endregion
}
