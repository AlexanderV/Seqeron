// SEQ-COMPLEX-KMER-001 — K-mer Entropy
// Evidence: docs/Evidence/SEQ-COMPLEX-KMER-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-COMPLEX-KMER-001.md
// Source: Shannon CE (1948) A Mathematical Theory of Communication; block (k-word) entropy of DNA:
// Herzel, Ebeling & Schmitt (1994) Phys Rev E 50:5061; reference implementation: BBMap/BBDuk
// EntropyTracker (pk = count / (window − k + 1)). Cross-check: scipy.stats.entropy(counts, base=2).
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

    // M1 — ACGT,k=1: 4 distinct monomers, uniform ⇒ H = log₂(4) = 2.0 (Pastore et al. 2025 max; Shannon uniform).
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

    // M3 — ATATAT,k=2: AT=3,TA=2 over N=5 ⇒ binary entropy of 0.6 = 0.97095... (Shannon formula over overlapping k-mers).
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

    // M5 — AAACGT,k=2: AA=2,AC=1,CG=1,GT=1 over N=5 ⇒ 1.92192... = log₂5 - 0.4 (Shannon formula over overlapping k-mers).
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
    #region Reference cross-check (scipy.stats.entropy) and canonical-counter consistency

    // R1 — scipy.stats.entropy(Counter(overlapping k-mers).values(), base=2) (scipy 1.x), computed
    // independently in Python over the same overlapping k-mer multiset. The 200-bp sequence was drawn
    // with random.seed(20260928); random.choice("ACGT").
    private const string Random200 =
        "GACATGCGATAGTAACGACTGGCCCCACCGGTAAGACCATTTAAGGCCAAGGAACAGCATACGACACGACGGGGCCACTGATACCTATTGAG" +
        "GACTTTTTATATCTGAGTGCAAGGAATCTGGCCGATATTCTGGATCACTCTATGCCAGTTTGCCTTATTGCCCGCATACCGGATAAAGTGAAAT" +
        "TAGGCCACTGTTAT";

    [TestCase("ATGCATGCAT", 2, 1.974937501201927)]
    [TestCase("ATGCGATCGATCG", 2, 2.4591479170272446)]
    [TestCase("ATGCGATCGATCG", 3, 2.7321588913645702)]
    [TestCase(Random200, 1, 1.9964735194730474)]
    [TestCase(Random200, 2, 3.937571048725419)]
    [TestCase(Random200, 3, 5.6801547658649625)]
    [TestCase(Random200, 5, 7.382517114501296)]
    [TestCase(Random200, 8, 7.582094342967564)]
    public void CalculateKmerEntropy_MatchesScipyReference(string sequence, int k, double expected)
    {
        Assert.That(Random200, Has.Length.EqualTo(200));

        double entropy = SequenceComplexity.CalculateKmerEntropy(new DnaSequence(sequence), k);

        Assert.That(entropy, Is.EqualTo(expected).Within(1e-12),
            $"scipy.stats.entropy over the overlapping {k}-mer counts of the sequence gives {expected}.");
    }

    // R2 — no duplicated counting: H computed from the canonical KmerAnalyzer.CountKmers table
    // (KMER-COUNT-001) with p_i = n_i / (L − k + 1) equals the method's result.
    [TestCase(Random200, 3)]
    [TestCase("AAACGT", 2)]
    public void CalculateKmerEntropy_AgreesWithCanonicalKmerCounts(string sequence, int k)
    {
        var counts = KmerAnalyzer.CountKmers(sequence, k);
        int n = sequence.Length - k + 1;
        double expected = 0;
        foreach (int c in counts.Values)
        {
            double p = (double)c / n;
            expected -= p * Math.Log2(p);
        }

        Assert.Multiple(() =>
        {
            Assert.That(counts.Values.Sum(), Is.EqualTo(n), "Σ n_i must equal N = L − k + 1 (overlapping).");
            Assert.That(SequenceComplexity.CalculateKmerEntropy(sequence, k), Is.EqualTo(expected).Within(1e-12));
        });
    }

    // R3 — for k = 1 on an A/C/G/T-only sequence, k-mer entropy is the per-base Shannon entropy.
    [Test]
    public void CalculateKmerEntropy_K1_EqualsPerBaseShannonEntropy()
    {
        Assert.That(SequenceComplexity.CalculateKmerEntropy(Random200, 1),
            Is.EqualTo(SequenceComplexity.CalculateShannonEntropy(Random200)).Within(1e-12));
    }

    #endregion

    #region Bias corrections and normalisation (WP13, F54)

    // Reference values: Miller–Madow = R package entropy 1.3.2 entropy.MillerMadow(y, unit = "log2");
    // Grassberger = Grassberger (2003, arXiv:physics/0307138) H = ln N − (1/N) Σ n_i G(n_i),
    // G(n) = ψ(n) + ½(−1)ⁿ[ψ((n+1)/2) − ψ(n/2)], evaluated with mpmath at 40 digits (= the G_1 = −γ − ln 2,
    // G_{2m+1} = G_{2m}, G_{2m+2} = G_{2m} + 2/(2m+1) recurrence); normalised = value / log₂ N (BBTools
    // EntropyTracker.calcEntropy multiplier 1/ln N). Not derived from the implementation.
    [TestCase("ATATAT", 2, 0.97095059445466864, 1.115220098543565, 1.2692841903863027)]
    [TestCase("ACGTACGTAAAAAAAAACGTACGT", 3, 2.5318692569751747, 2.7286003989145788, 2.8297096977806522)]
    [TestCase("ATGCATGCAT", 2, 1.974937501201927, 2.2153866746834209, 2.1172810969412527)]
    public void CalculateKmerEntropy_Corrections_MatchReferenceEstimators(
        string sequence, int k, double plugin, double millerMadow, double grassberger)
    {
        var dna = new DnaSequence(sequence);
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateKmerEntropy(dna, k, KmerEntropyCorrection.None), Is.EqualTo(plugin).Within(1e-12));
            Assert.That(SequenceComplexity.CalculateKmerEntropy(dna, k, KmerEntropyCorrection.MillerMadow), Is.EqualTo(millerMadow).Within(1e-12));
            Assert.That(SequenceComplexity.CalculateKmerEntropy(dna, k, KmerEntropyCorrection.Grassberger), Is.EqualTo(grassberger).Within(1e-12));
            Assert.That(SequenceComplexity.CalculateKmerEntropy(sequence.ToLowerInvariant(), k, KmerEntropyCorrection.Grassberger),
                Is.EqualTo(grassberger).Within(1e-12), "string overload, case-insensitive");
            Assert.That(SequenceComplexity.CalculateKmerEntropy(dna, k, KmerEntropyCorrection.None),
                Is.EqualTo(SequenceComplexity.CalculateKmerEntropy(dna, k)), "None = legacy plug-in exactly");
        });
    }

    [TestCase("ATATAT", 2, 0.4181656600790516, 0.48029915353501278, 0.54665094633254617)]
    [TestCase("ACGTACGTAAAAAAAAACGTACGT", 3, 0.5677560446030244, 0.61187178821420702, 0.6345449240558931)]
    public void CalculateKmerEntropy_Normalized_DividesByLog2N(
        string sequence, int k, double plugin, double millerMadow, double grassberger)
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateKmerEntropy(sequence, k, KmerEntropyCorrection.None, normalize: true), Is.EqualTo(plugin).Within(1e-12));
            Assert.That(SequenceComplexity.CalculateKmerEntropy(sequence, k, KmerEntropyCorrection.MillerMadow, normalize: true), Is.EqualTo(millerMadow).Within(1e-12));
            Assert.That(SequenceComplexity.CalculateKmerEntropy(sequence, k, KmerEntropyCorrection.Grassberger, normalize: true), Is.EqualTo(grassberger).Within(1e-12));
        });
    }

    // BBTools 40.02 tracker.EntropyTracker.calcEntropy("ATATAT".getBytes(), null, 2) = 0.41816565f (float).
    [Test]
    public void CalculateKmerEntropy_Normalized_MatchesBbtoolsEntropyTracker()
    {
        double e = SequenceComplexity.CalculateKmerEntropy("ATATAT", 2, KmerEntropyCorrection.None, normalize: true);
        Assert.That((float)e, Is.EqualTo(0.41816565f));
    }

    // N = 1: Miller–Madow adds (1 − 1)/2 = 0; Grassberger gives ln 1 − G(1) = γ + ln 2 nats = 1.8327461772768672 bits
    // (mpmath); normalisation by log₂ 1 = 0 is defined as 0. Homopolymer A×10, k = 1: Grassberger = −0.0023880009817158878
    // (the estimator can be slightly negative; mpmath).
    [Test]
    public void CalculateKmerEntropy_Corrections_EdgeCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateKmerEntropy("ACG", 3, KmerEntropyCorrection.MillerMadow), Is.EqualTo(0.0));
            Assert.That(SequenceComplexity.CalculateKmerEntropy("ACG", 3, KmerEntropyCorrection.Grassberger), Is.EqualTo(1.8327461772768672).Within(1e-12));
            Assert.That(SequenceComplexity.CalculateKmerEntropy("ACG", 3, KmerEntropyCorrection.Grassberger, normalize: true), Is.EqualTo(0.0));
            Assert.That(SequenceComplexity.CalculateKmerEntropy("AAAAAAAAAA", 1, KmerEntropyCorrection.Grassberger), Is.EqualTo(-0.0023880009817158878).Within(1e-12));
            Assert.That(SequenceComplexity.CalculateKmerEntropy("AC", 3, KmerEntropyCorrection.Grassberger), Is.EqualTo(0.0), "L < k");
            Assert.That(SequenceComplexity.CalculateKmerEntropy((string)null!, 2, KmerEntropyCorrection.MillerMadow), Is.EqualTo(0.0));
        });
    }

    [Test]
    public void CalculateKmerEntropy_Corrections_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateKmerEntropy("ACGT", 0, KmerEntropyCorrection.MillerMadow));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateKmerEntropy("ACGT", 2, (KmerEntropyCorrection)7));
            Assert.Throws<ArgumentNullException>(() => SequenceComplexity.CalculateKmerEntropy((DnaSequence)null!, 2, KmerEntropyCorrection.None));
            Assert.Throws<ArgumentException>(() => SequenceComplexity.ParseKmerEntropyCorrection("chao-shen"));
            Assert.That(SequenceComplexity.ParseKmerEntropyCorrection("Miller-Madow"), Is.EqualTo(KmerEntropyCorrection.MillerMadow));
            Assert.That(SequenceComplexity.ParseKmerEntropyCorrection("grassberger"), Is.EqualTo(KmerEntropyCorrection.Grassberger));
            Assert.That(SequenceComplexity.ParseKmerEntropyCorrection(null), Is.EqualTo(KmerEntropyCorrection.None));
        });
    }

    // Miller–Madow = plug-in + (D − 1)/(2N ln 2) exactly; random 200-mer, k = 1..8.
    [Test]
    public void CalculateKmerEntropy_MillerMadow_IsPluginPlusMillerTerm([Range(1, 8)] int k)
    {
        int d = KmerAnalyzer.CountKmers(Random200, k).Count;
        int n = Random200.Length - k + 1;
        double expected = SequenceComplexity.CalculateKmerEntropy(Random200, k) + (d - 1) / (2.0 * n) / Math.Log(2);
        Assert.That(SequenceComplexity.CalculateKmerEntropy(Random200, k, KmerEntropyCorrection.MillerMadow), Is.EqualTo(expected).Within(1e-13));
    }

    #endregion
}
