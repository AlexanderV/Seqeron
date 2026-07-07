// SEQ-COMPLEX-KMER-001 — K-mer Entropy
// Evidence: docs/Evidence/SEQ-COMPLEX-KMER-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-COMPLEX-KMER-001.md
// Source: Li H (2025) longdust, arXiv:2509.07357; Shannon CE (1948) A Mathematical Theory of Communication.
//
// Spec: H = -Σ p_i·log₂(p_i) over the N = L-k+1 overlapping k-mers, p_i = count_i / N (bits).
// Expected values are derived independently from the formula, NOT from the implementation.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceComplexity_CalculateKmerEntropy_Tests
{
    // Exact constants derived by hand from H = -Σ p log₂ p (see Evidence §"Test Datasets").
    private const double Log2Of3 = 1.5849625007211562;          // ACGT, k=2: 3 distinct dimers, uniform
    private const double BinaryEntropyOf06 = 0.9709505944546686; // ATATAT, k=2: AT=3,TA=2 (p=0.6/0.4)
    private const double MixedCountsEntropy = 1.9219280948873623; // AAACGT, k=2: 2,1,1,1 = log₂5 - 0.4

    #region CalculateKmerEntropy(DnaSequence, k) — canonical exact values

    // M1 — ACGT,k=1: 4 distinct monomers, uniform ⇒ H = log₂(4) = 2.0 (Çakır 2025 max; Shannon uniform).
    [Test]
    public void CalculateKmerEntropy_UniformMonomers_ReturnsLog2Of4()
    {
        var seq = new DnaSequence("ACGT");

        double entropy = SequenceComplexity.CalculateKmerEntropy(seq, k: 1);

        Assert.That(entropy, Is.EqualTo(2.0).Within(1e-10),
            "ACGT,k=1 has 4 equiprobable monomers; uniform Shannon entropy = log₂(4) = 2.0 bits.");
    }

    // M2 — ACGT,k=2: 3 distinct dimers AC,CG,GT each once ⇒ H = log₂(3) (all-distinct/uniform).
    [Test]
    public void CalculateKmerEntropy_AllDistinctDimers_ReturnsLog2OfN()
    {
        var seq = new DnaSequence("ACGT");

        double entropy = SequenceComplexity.CalculateKmerEntropy(seq, k: 2);

        Assert.That(entropy, Is.EqualTo(Log2Of3).Within(1e-10),
            "ACGT,k=2 yields N=3 distinct dimers, each p=1/3; uniform entropy = log₂(3).");
    }

    // M3 — ATATAT,k=2: AT=3,TA=2 over N=5 ⇒ binary entropy of 0.6 = 0.97095... (Li 2025 formula).
    [Test]
    public void CalculateKmerEntropy_NonUniformDimers_ReturnsExact()
    {
        var seq = new DnaSequence("ATATAT");

        double entropy = SequenceComplexity.CalculateKmerEntropy(seq, k: 2);

        Assert.That(entropy, Is.EqualTo(BinaryEntropyOf06).Within(1e-10),
            "ATATAT,k=2: AT=3,TA=2,N=5 ⇒ H=-0.6·log₂0.6-0.4·log₂0.4=0.9709505944546686.");
    }

    // M4 — AAAA,k=2: only AA (deterministic) ⇒ H = 0 (Shannon: certainty ⇒ entropy 0).
    [Test]
    public void CalculateKmerEntropy_SingleRepeatedDimer_ReturnsZero()
    {
        var seq = new DnaSequence("AAAA");

        double entropy = SequenceComplexity.CalculateKmerEntropy(seq, k: 2);

        Assert.That(entropy, Is.EqualTo(0.0).Within(1e-10),
            "AAAA,k=2 has a single dimer AA (p=1); a deterministic distribution has entropy 0.");
    }

    // M5 — AAACGT,k=2: AA=2,AC=1,CG=1,GT=1 over N=5 ⇒ 1.92192... = log₂5 - 0.4 (Li 2025 formula).
    [Test]
    public void CalculateKmerEntropy_MixedCounts_ReturnsExact()
    {
        var seq = new DnaSequence("AAACGT");

        double entropy = SequenceComplexity.CalculateKmerEntropy(seq, k: 2);

        Assert.That(entropy, Is.EqualTo(MixedCountsEntropy).Within(1e-10),
            "AAACGT,k=2: AA=2,AC=1,CG=1,GT=1,N=5 ⇒ H=-(0.4·log₂0.4+3·0.2·log₂0.2)=1.9219280948873623.");
    }

    // M6 — sequence shorter than k ⇒ no k-mers ⇒ 0 (documented boundary).
    [Test]
    public void CalculateKmerEntropy_SequenceShorterThanK_ReturnsZero()
    {
        var seq = new DnaSequence("AC");

        double entropy = SequenceComplexity.CalculateKmerEntropy(seq, k: 5);

        Assert.That(entropy, Is.EqualTo(0.0).Within(1e-10),
            "L=2 < k=5: no k-mers exist, so the entropy of the empty distribution is 0.");
    }

    // M7 — invalid k (<1) ⇒ ArgumentOutOfRangeException (contract).
    [Test]
    public void CalculateKmerEntropy_InvalidK_Throws()
    {
        var seq = new DnaSequence("ACGT");

        Assert.Throws<ArgumentOutOfRangeException>(
            () => SequenceComplexity.CalculateKmerEntropy(seq, k: 0),
            "k must be ≥ 1; k=0 is an invalid window length.");
    }

    // M8 — null DnaSequence ⇒ ArgumentNullException (contract).
    [Test]
    public void CalculateKmerEntropy_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => SequenceComplexity.CalculateKmerEntropy((DnaSequence)null!, k: 2),
            "A null DnaSequence must be rejected with ArgumentNullException.");
    }

    #endregion

    #region CalculateKmerEntropy(string, k) — delegate overload

    // S1 — string overload agrees with DnaSequence overload (delegation + INV-04).
    [Test]
    public void CalculateKmerEntropy_StringOverload_AgreesWithDnaSequence()
    {
        const string s = "ATATAT";

        double fromString = SequenceComplexity.CalculateKmerEntropy(s, k: 2);
        double fromDna = SequenceComplexity.CalculateKmerEntropy(new DnaSequence(s), k: 2);

        Assert.Multiple(() =>
        {
            Assert.That(fromString, Is.EqualTo(BinaryEntropyOf06).Within(1e-10),
                "string overload computes the same evidence value 0.9709505944546686.");
            Assert.That(fromString, Is.EqualTo(fromDna).Within(1e-10),
                "string overload must delegate to the same core as the DnaSequence overload.");
        });
    }

    // S2 — case-insensitivity: lowercase equals uppercase (INV-04; input upper-cased).
    [Test]
    public void CalculateKmerEntropy_LowercaseString_EqualsUppercase()
    {
        double lower = SequenceComplexity.CalculateKmerEntropy("atatat", k: 2);
        double upper = SequenceComplexity.CalculateKmerEntropy("ATATAT", k: 2);

        Assert.That(lower, Is.EqualTo(upper).Within(1e-10),
            "Input is normalised to upper-case, so entropy is case-insensitive.");
    }

    // S3 — null/empty string ⇒ 0 (string overload contract).
    [Test]
    public void CalculateKmerEntropy_NullOrEmptyString_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateKmerEntropy((string)null!, k: 2),
                Is.EqualTo(0.0).Within(1e-10), "null string ⇒ 0 per overload contract.");
            Assert.That(SequenceComplexity.CalculateKmerEntropy("", k: 2),
                Is.EqualTo(0.0).Within(1e-10), "empty string ⇒ 0 (no k-mers).");
        });
    }

    // S4 — string overload has its OWN k<1 guard (distinct code path from the DnaSequence
    // overload tested in M7); k=0 must throw ArgumentOutOfRangeException.
    [Test]
    public void CalculateKmerEntropy_StringOverload_InvalidK_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => SequenceComplexity.CalculateKmerEntropy("ACGT", k: 0),
            "string overload must reject k<1 with ArgumentOutOfRangeException.");
    }

    // S5 — string overload, non-empty sequence shorter than k ⇒ 0 (L<k core branch via string path).
    // "AC", k=5: L=2 < k=5 ⇒ no k-mers ⇒ 0 (independently confirmed: empty multiset entropy = 0).
    [Test]
    public void CalculateKmerEntropy_StringOverload_SequenceShorterThanK_ReturnsZero()
    {
        Assert.That(SequenceComplexity.CalculateKmerEntropy("AC", k: 5),
            Is.EqualTo(0.0).Within(1e-10),
            "string overload, L=2 < k=5: no k-mers exist ⇒ entropy 0.");
    }

    #endregion

    #region Invariant (C1)

    // C1 — INV-01: 0 ≤ H ≤ log₂(N), N = L-k+1, for any valid input (Shannon bounds).
    [TestCase("ACGTACGTAA", 2)]
    [TestCase("AAAAAAAA", 3)]
    [TestCase("ACGTACGTACGT", 4)]
    [TestCase("GCGCGCGCGCAT", 1)]
    public void CalculateKmerEntropy_BoundsInvariant_WithinRange(string sequence, int k)
    {
        int n = sequence.Length - k + 1;
        double maxEntropy = Math.Log2(n);

        double entropy = SequenceComplexity.CalculateKmerEntropy(new DnaSequence(sequence), k);

        Assert.Multiple(() =>
        {
            Assert.That(entropy, Is.GreaterThanOrEqualTo(0.0),
                "Shannon entropy is non-negative.");
            Assert.That(entropy, Is.LessThanOrEqualTo(maxEntropy + 1e-10),
                $"Entropy cannot exceed log₂(N)=log₂({n}) for N={n} k-mers.");
        });
    }

    #endregion
}
