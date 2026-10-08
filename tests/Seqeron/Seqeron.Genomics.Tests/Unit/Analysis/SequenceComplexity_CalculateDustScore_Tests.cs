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

    // F34: DUST is defined for triplets only (Morgulis 2006; NCBI symdust triplet_type; lh3/sdust
    // SD_WLEN = 3; longdust README "SDUST … hardcodes k = 3"). Other word sizes were an unsourced
    // extrapolation and are rejected; the sourced k-mer generalisation is longdust.
    [TestCase(1)]
    [TestCase(2)]
    [TestCase(4)]
    [TestCase(7)]
    public void CalculateDustScore_WordSizeNotThree_Throws(int wordSize)
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateDustScore(new DnaSequence("AAAAAAAA"), wordSize));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateDustScore("AAAAAAAA", wordSize));
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

    #region F35 — dustmasker linker + soft masking (NCBI dustmasker 2.12.0 binary output)

    // 96 bp: A-run, (CA)-run, A-run. dustmasker (-window 64 -level 20) -outfmt interval reports the
    // closed intervals 10–25, 43–58, 81–93 for linker 1..17, 10–58 + 81–93 for linker 18 and 19,
    // and 10–93 for linker 32 (symdust save_masked_regions: prev.last + linker >= next.first).
    private const string LinkerSeq =
        "ACGTGCATGCAAAAAAAAAAAAAAAAGCTAGCATCGACTGCAGCACACACACACACACAGATCGATCGTACGGTGCATGACAAAAAAAAAAAAACT";

    [TestCase(1, new[] { 10, 26, 43, 59, 81, 94 })]
    [TestCase(17, new[] { 10, 26, 43, 59, 81, 94 })]
    [TestCase(18, new[] { 10, 59, 81, 94 })]
    [TestCase(19, new[] { 10, 59, 81, 94 })]
    [TestCase(32, new[] { 10, 94 })]
    public void FindLowComplexityIntervals_Linker_MatchesDustmasker(int linker, int[] flat)
    {
        var expected = Enumerable.Range(0, flat.Length / 2).Select(i => (flat[2 * i], flat[2 * i + 1])).ToList();

        var got = SequenceComplexity.FindLowComplexityIntervals(LinkerSeq, linker: linker);

        Assert.That(got, Is.EqualTo(expected));
    }

    [Test]
    public void FindLowComplexityIntervals_DefaultLinker_EqualsSdustAndMask()
    {
        // lh3/sdust (-w 64 -t 20): x 10 26 / x 43 59 / x 81 94; linker 1 is the sdust merge rule.
        var intervals = SequenceComplexity.FindLowComplexityIntervals(new DnaSequence(LinkerSeq));
        string masked = SequenceComplexity.MaskLowComplexity(new DnaSequence(LinkerSeq));

        Assert.Multiple(() =>
        {
            Assert.That(intervals, Is.EqualTo(new[] { (10, 26), (43, 59), (81, 94) }));
            for (int i = 0; i < LinkerSeq.Length; i++)
                Assert.That(masked[i] == 'N', Is.EqualTo(intervals.Any(p => i >= p.Start && i < p.End)), $"position {i}");
        });
    }

    [Test]
    public void MaskLowComplexity_SoftMask_MatchesDustmaskerFasta()
    {
        // dustmasker -outfmt fasta (soft masking), linker 1 and 32; lower-case input is upper-cased.
        const string Linker1 = "ACGTGCATGCaaaaaaaaaaaaaaaaGCTAGCATCGACTGCAGcacacacacacacacaGATCGATCGTACGGTGCATGACaaaaaaaaaaaaaCT";
        const string Linker32 = "ACGTGCATGCaaaaaaaaaaaaaaaagctagcatcgactgcagcacacacacacacacagatcgatcgtacggtgcatgacaaaaaaaaaaaaaCT";

        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.MaskLowComplexity(LinkerSeq, softMask: true), Is.EqualTo(Linker1));
            Assert.That(SequenceComplexity.MaskLowComplexity(LinkerSeq.ToLowerInvariant(), softMask: true), Is.EqualTo(Linker1));
            Assert.That(SequenceComplexity.MaskLowComplexity(new DnaSequence(LinkerSeq), 64, 2.0, 'N', linker: 32, softMask: true),
                Is.EqualTo(Linker32));
            Assert.That(SequenceComplexity.MaskLowComplexity(LinkerSeq, maskChar: 'X', linker: 18),
                Is.EqualTo(LinkerSeq[..10] + new string('X', 49) + LinkerSeq[59..81] + new string('X', 13) + LinkerSeq[94..]));
        });
    }

    [Test]
    public void SdustLinkerAndString_InvalidParameters_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals(LinkerSeq, linker: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.MaskLowComplexity(LinkerSeq, linker: -1));
            // symdust accepts 1–32 only (dustmasker silently falls back to 1 outside it): rejected here.
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals(LinkerSeq, linker: 33));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.MaskLowComplexity(new DnaSequence(LinkerSeq), 64, 2.0, 'N', linker: 33));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.MaskLowComplexity(LinkerSeq, windowSize: 2));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals(LinkerSeq, threshold: double.PositiveInfinity));
            Assert.Throws<ArgumentNullException>(() => SequenceComplexity.MaskLowComplexity((string)null!));
            Assert.Throws<ArgumentNullException>(() => SequenceComplexity.FindLowComplexityIntervals((string)null!));
            Assert.That(SequenceComplexity.MaskLowComplexity(string.Empty), Is.Empty);
        });
    }

    #endregion

    #region F36 — SDUST on strings with N / IUPAC (each ACGT run = independent sdust input)

    // Expected = compiled lh3/sdust run on every maximal ACGT piece separately, shifted to input
    // coordinates (sdust's documented contract "N effectively breaks input into pieces of independent
    // sequences"). Upstream sdust_core on the whole string leaks its window across the N and reports
    // shifted intervals instead — for the first row (8,20),(21,34),(35,72) on a 53-bp input.
    private static readonly object[] SdustNCases =
    {
        new object[] { "ACGTNNAAAAAAAAAAAANACGTACACACACACACACANNGGGCCCTAGGTCA", 64, 2.0, new[] { 6, 18, 23, 38 } },
        new object[] { "ACGTNNAAAAAAAAAAAANACGTACACACACACACACANNGGGCCCTAGGTCA", 16, 1.5, new[] { 6, 18, 23, 38 } },
        new object[] { "ACGTTGCAGTCATGCGATCAAAAAAANTGCATCGGATCCAAAAAAATAGG", 64, 2.0, new[] { 19, 26, 39, 46 } },
        new object[] { "acgtgcatgcAAAAAAAAAAAAAAAAgctaRcatcgacACACACACACACACACAGATnGATCG", 64, 2.0, new[] { 10, 26, 36, 55 } },
        new object[] { "NNNN", 64, 2.0, Array.Empty<int>() },
    };

    [TestCaseSource(nameof(SdustNCases))]
    public void FindLowComplexityIntervals_NonAcgt_EqualsPerPieceSdust(string input, int window, double threshold, int[] flat)
    {
        var expected = Enumerable.Range(0, flat.Length / 2).Select(i => (flat[2 * i], flat[2 * i + 1])).ToList();

        var got = SequenceComplexity.FindLowComplexityIntervals(input, window, threshold);
        string masked = SequenceComplexity.MaskLowComplexity(input, window, threshold);

        Assert.Multiple(() =>
        {
            Assert.That(got, Is.EqualTo(expected));
            Assert.That(masked, Has.Length.EqualTo(input.Length));
            for (int i = 0; i < input.Length; i++)
            {
                bool inInterval = expected.Any(p => i >= p.Item1 && i < p.Item2);
                Assert.That(masked[i], Is.EqualTo(inInterval ? 'N' : char.ToUpperInvariant(input[i])), $"position {i}");
            }
        });
    }

    [Test]
    public void MaskLowComplexity_String_AcgtOnly_EqualsDnaSequencePath()
    {
        foreach (var row in SdustReferenceCases.Cast<object[]>())
        {
            var (input, window, threshold) = ((string)row[0], (int)row[1], (double)row[2]);
            Assert.That(SequenceComplexity.MaskLowComplexity(input.ToLowerInvariant(), window, threshold), Is.EqualTo((string)row[3]));
        }
    }

    [Test]
    public void MaskLowComplexity_SoftMask_KeepsNonAcgtAndLowercasesIntervals()
    {
        // pieces: A×12 at [6,18) and the (AC) run at [23,38); N outside intervals stays upper-case.
        Assert.That(SequenceComplexity.MaskLowComplexity("ACGTNNAAAAAAAAAAAANACGTACACACACACACACANNGGGCCCTAGGTCA", softMask: true),
            Is.EqualTo("ACGTNNaaaaaaaaaaaaNACGTacacacacacacacaNNGGGCCCTAGGTCA"));
    }

    #endregion

    #region F53 — DustEngine.Dustmasker: NCBI dustmasker 2.12.0 parity (symdust + GetDustMasks_SkipNs)

    // Expected values are the verbatim output of the NCBI dustmasker 2.12.0 binary (Debian ncbi-blast+,
    // `-outfmt interval` closed "a - b" → half-open [a, b+1); `-outfmt fasta` joined lines), NOT the
    // implementation's output. Sources: ncbi-cxx-toolkit-public algo/dustmask/symdust.cpp/.hpp,
    // app/dustmasker/dust_mask_app.cpp, util/random_gen.cpp/.hpp (public domain).
    private const string DmNSeq = "NNACGTTGCAAAAAAAAAAAACGTRRRRRRRRRRRRTGCANNNNNNNNNNACGTTGCAANNNN";
    private const string DmShortNSeq = "ACGTTGCAAGCTTCGATGCAAAAAAAAAAAAAAANNNNNAAAAAAAAAAAAAAAACGTTGCAAGCTTCGATGC";
    private const string DmPolyGSeq = "TACTGTCCGGTGATTGGTGTCTCTGTAACATTACTAACTTTTCCACGCTTGTTCTGTTCGGGGGGGGGGGGGGGGGGCAAGGGCGGAACGGGTCC";

    private static (int, int)[] Pairs(int[] flat) =>
        Enumerable.Range(0, flat.Length / 2).Select(i => (flat[2 * i], flat[2 * i + 1])).ToArray();

    // dustmasker -window 8 -level 20 (-linker L): leading/trailing N runs and the interior N run of 10 > W
    // are masked and cut the scan; R codes are scanned as A (R×12 masked like poly-A); s_InsertMerge joins an
    // N run only when prev.end + linker == N.start exactly (linker 5 bridges 35→40, linker 6 does not).
    [TestCase(1, new[] { 0, 2, 9, 21, 24, 36, 40, 50, 59, 63 })]
    [TestCase(3, new[] { 0, 2, 9, 21, 24, 36, 40, 50, 59, 63 })]
    [TestCase(4, new[] { 0, 2, 9, 36, 40, 50, 59, 63 })]
    [TestCase(5, new[] { 0, 2, 9, 50, 59, 63 })]
    [TestCase(6, new[] { 0, 2, 9, 36, 40, 50, 59, 63 })]
    [TestCase(10, new[] { 0, 2, 9, 36, 40, 50, 59, 63 })]
    public void FindLowComplexityIntervals_Dustmasker_NRunsAndIupac_MatchDustmasker(int linker, int[] flat)
    {
        var got = SequenceComplexity.FindLowComplexityIntervals(DmNSeq, 8, 2.0, linker, DustEngine.Dustmasker);
        Assert.That(got, Is.EqualTo(Pairs(flat)));
    }

    [Test]
    public void FindLowComplexityIntervals_Dustmasker_DiffersFromSdustOnNonAcgt()
    {
        Assert.Multiple(() =>
        {
            // sdust engine: R and N split the scan, so only the A run is masked.
            Assert.That(SequenceComplexity.FindLowComplexityIntervals(DmNSeq, 8, 2.0),
                Is.EqualTo(new[] { (9, 21) }));
            // dustmasker (default W 64): the short N run (5 ≤ W) is scanned as CRandom bases inside the poly-A
            // and the whole tract is one interval; sdust (per ACGT run) gives two.
            Assert.That(SequenceComplexity.FindLowComplexityIntervals(DmShortNSeq, engine: DustEngine.Dustmasker),
                Is.EqualTo(new[] { (19, 55) }));
            Assert.That(SequenceComplexity.FindLowComplexityIntervals(DmShortNSeq),
                Is.EqualTo(new[] { (19, 34), (39, 55) }));
            // dustmasker -window 8: the CRandom codes of NNNNN (positions 34–38) leave only N@34 inside a
            // perfect interval — locks the port of CRandom (LFG 33/13, Reset() table, GetRand() >> 1).
            Assert.That(SequenceComplexity.FindLowComplexityIntervals(DmShortNSeq, 8, 2.0, engine: DustEngine.Dustmasker),
                Is.EqualTo(new[] { (19, 35), (39, 55) }));
        });
    }

    // A window holding one triplet value only (homopolymer ≥ W) is masked by symdust's num_diff ≤ 1 shortcut
    // (shift_window / shift_high insert a perfect interval with no score test). With W = 8 the best possible
    // DUST score is 15/5 = 3.0, so at level 30 sdust masks nothing while dustmasker masks the G×18 run.
    // This — not the size of thresholds_ (W − 2 entries) — is the W = 8 dustmasker/sdust difference of F35.
    [Test]
    public void FindLowComplexityIntervals_Dustmasker_SingleTripletWindowShortcut()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.FindLowComplexityIntervals(DmPolyGSeq, 8, 3.0, engine: DustEngine.Dustmasker),
                Is.EqualTo(new[] { (59, 77) }));
            Assert.That(SequenceComplexity.FindLowComplexityIntervals(new DnaSequence(DmPolyGSeq), 8, 3.0, engine: DustEngine.Dustmasker),
                Is.EqualTo(new[] { (59, 77) }));
            Assert.That(SequenceComplexity.FindLowComplexityIntervals(DmPolyGSeq, 8, 3.0), Is.Empty);
        });
    }

    [Test]
    public void MaskLowComplexity_Dustmasker_SoftMask_MatchesDustmaskerFasta()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.MaskLowComplexity(DmNSeq, 8, 2.0, softMask: true, engine: DustEngine.Dustmasker),
                Is.EqualTo("nnACGTTGCaaaaaaaaaaaaCGTrrrrrrrrrrrrTGCAnnnnnnnnnnACGTTGCAAnnnn"));
            // Lower-case input is upper-cased, as dustmasker's FASTA reader does.
            Assert.That(SequenceComplexity.MaskLowComplexity(DmNSeq.ToLowerInvariant(), softMask: true, engine: DustEngine.Dustmasker),
                Is.EqualTo("nnACGTTGCaaaaaaaaaaaacgtrrrrrrrrrrrrTGCANNNNNNNNNNACGTTGCAAnnnn"));
            Assert.That(SequenceComplexity.MaskLowComplexity(DmShortNSeq, softMask: true, engine: DustEngine.Dustmasker),
                Is.EqualTo("ACGTTGCAAGCTTCGATGCaaaaaaaaaaaaaaannnnnaaaaaaaaaaaaaaaaCGTTGCAAGCTTCGATGC"));
            Assert.That(SequenceComplexity.MaskLowComplexity(new DnaSequence(DmPolyGSeq), 8, 3.0, 'X', 1, engine: DustEngine.Dustmasker),
                Is.EqualTo(DmPolyGSeq[..59] + new string('X', 18) + DmPolyGSeq[77..]));
        });
    }

    [Test]
    public void Dustmasker_AcgtDefaults_EqualSdustEngine()
    {
        // At W = 64 / level 20 the single-triplet shortcut never fires below the score test, so on ACGT input the
        // two engines agree (dustmasker 2.12.0: 10-25 43-58 81-93 and, -linker 32, 10-93).
        const string seq = "ACGTGCATGCAAAAAAAAAAAAAAAAGCTAGCATCGACTGCAGCACACACACACACACAGATCGATCGTACGGTGCATGACAAAAAAAAAAAAACT";
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.FindLowComplexityIntervals(seq, engine: DustEngine.Dustmasker),
                Is.EqualTo(new[] { (10, 26), (43, 59), (81, 94) }));
            Assert.That(SequenceComplexity.FindLowComplexityIntervals(seq, linker: 32, engine: DustEngine.Dustmasker),
                Is.EqualTo(new[] { (10, 94) }));
        });
    }

    [Test]
    public void Dustmasker_InvalidParameters_Throw()
    {
        Assert.Multiple(() =>
        {
            // symdust accepts window 8–64 and level 2–64 (dustmasker silently substitutes 64 / 20 otherwise).
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 7, 2.0, engine: DustEngine.Dustmasker));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 65, 2.0, engine: DustEngine.Dustmasker));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 64, 0.1, engine: DustEngine.Dustmasker));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 64, 6.5, engine: DustEngine.Dustmasker));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 64, 2.05, engine: DustEngine.Dustmasker));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 64, 2.0, 33, DustEngine.Dustmasker));
            // dustmasker's reader turns U into T, '-' into N and drops X/'*' (coordinates shift): rejected here.
            Assert.Throws<ArgumentException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGU", engine: DustEngine.Dustmasker));
            Assert.Throws<ArgumentException>(() => SequenceComplexity.MaskLowComplexity("ACG-T", engine: DustEngine.Dustmasker));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", engine: (DustEngine)7));
            Assert.That(SequenceComplexity.FindLowComplexityIntervals("", engine: DustEngine.Dustmasker), Is.Empty);
            Assert.That(SequenceComplexity.ParseDustEngine("DustMasker"), Is.EqualTo(DustEngine.Dustmasker));
            Assert.That(SequenceComplexity.ParseDustEngine(null), Is.EqualTo(DustEngine.Sdust));
            Assert.Throws<ArgumentException>(() => SequenceComplexity.ParseDustEngine("seg"));
        });
    }

    #endregion

    #region F34 — longdust (Li & Li 2025) k-mer generalisation, compiled lh3/longdust 1.4-r97 output

    private const string LdVntr =
        "ATCAGTCATTAAACTATAAACCACTTGAACCACAACGATGTCGTTTATAGCGCGCGGGGACGGCAGCTGCGATACCCCCTCGAATCCCCGGCGGCTCTCACCTGCAGGGTGGACGTTTG" +
        "GGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCA" +
        "GGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCA" +
        "ACCGAGCCTCAACGGAAAGGCGGCATTGGGCGTAGATCATTGTAAGAATTGAGAGGACTGAGGGATAGGGAAAGGTACGGGCCCCGATTTCCCATGCAGGCATCTCCAAGTGTAAGCACG";

    private const string LdStrN =
        "CGGACACACCCTCAACCAAGCGCGTTCCGCCGCGGCTGTACCACAGGCCTTTATGTCAGCAGAAAGAGGGCATACAGCGGCACACACACACACACACACACACACACACACACACACACACA" +
        "CACACACACACACACANNNNGCAATCAGACCGCTCTTAGCCCATACGCAATTTTCGGCGGAAGGCTTGCCTCCCGAGGACTTTTTTTTTTTTTTTTTTTTTTTTTGTCAGAGGGGTATCTG" +
        "TACTGAGTGCGTCGATATCATGTTTTCAGCAATAACAGGTTAGTCTCCTGCTTATCGTTGCCTG";

    [Test]
    public void FindLongdustRegions_MatchesCompiledLongdust()
    {
        Assert.Multiple(() =>
        {
            Assert.That(LdVntr, Has.Length.EqualTo(408));
            Assert.That(LdStrN, Has.Length.EqualTo(307));
            // longdust (defaults -k7 -w5000 -t0.6 -e50 -b3, both strands)
            Assert.That(SequenceComplexity.FindLongdustRegions(LdVntr), Is.EqualTo(new[] { (120, 288) }));
            Assert.That(SequenceComplexity.FindLongdustRegions(LdStrN), Is.EqualTo(new[] { (80, 138), (202, 227) }));
            // longdust -k5 -w100
            Assert.That(SequenceComplexity.FindLongdustRegions(LdVntr, k: 5, windowSize: 100), Is.EqualTo(new[] { (117, 291) }));
            // longdust -f (forward strand only)
            Assert.That(SequenceComplexity.FindLongdustRegions(LdStrN, forwardOnly: true), Is.EqualTo(new[] { (80, 138), (202, 227) }));
            // longdust -k4 -w200 -t1.0
            Assert.That(SequenceComplexity.FindLongdustRegions(LdVntr, k: 4, windowSize: 200, threshold: 1.0), Is.EqualTo(new[] { (120, 288) }));
            Assert.That(SequenceComplexity.FindLongdustRegions(string.Empty), Is.Empty);
            Assert.That(SequenceComplexity.FindLongdustRegions("NNNNNNNNNN"), Is.Empty);
        });
    }

    // S_L(x) = Σ_t log c(t)! − f(ℓ/4^k): values from longdust's own f() table (ld_cal_f / ld_cal_f2 with -g)
    // and the same in-order Σ log c accumulation, compiled from lh3/longdust 1.4-r97.
    [TestCase("AAAAAAAAAA", 3, null, 10.263913947178825)]
    [TestCase("AAAAAAAAAA", 3, 0.41, 10.230970927904716)]
    [TestCase("ACACACACACACACACACACACACACACAC", 7, null, 39.962247232206003)]
    [TestCase("ACGT", 2, null, -0.1900272033309148)]
    [TestCase("ACGTNACGTACGTACGT", 4, null, 4.9941414253104837)]
    [TestCase("ACGTNACGTACGTACGT", 4, 0.41, 4.9589302235838071)]
    public void CalculateLongdustScore_MatchesLongdustF(string seq, int k, double? gc, double expected)
    {
        Assert.That(SequenceComplexity.CalculateLongdustScore(seq, k, gc), Is.EqualTo(expected).Within(1e-12));
    }

    [Test]
    public void CalculateLongdustScore_VntrExceedsThresholdTimesLength()
    {
        // longdust calls x low-complexity when S_L(x) − T·ℓ(x) > 0 (T = 0.6): true for the VNTR region
        // [120, 288) that longdust reports, false for a random flank.
        string vntr = LdVntr[120..288];
        string flank = LdVntr[..120];
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateLongdustScore(vntr) - 0.6 * (vntr.Length - 6), Is.GreaterThan(0));
            Assert.That(SequenceComplexity.CalculateLongdustScore(flank) - 0.6 * (flank.Length - 6), Is.LessThan(0));
            Assert.That(SequenceComplexity.CalculateLongdustScore("ACGTAC"), Is.EqualTo(0.0)); // ℓ = 0
        });
    }

    [Test]
    public void Longdust_InvalidParameters_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => SequenceComplexity.FindLongdustRegions(null!));
            Assert.Throws<ArgumentNullException>(() => SequenceComplexity.CalculateLongdustScore(null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", k: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", k: 15));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", windowSize: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", windowSize: 65535));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", threshold: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", threshold: double.NaN));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", xdropLength: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", minStartCount: 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", gcContent: 1.0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateLongdustScore("ACGT", k: 15));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateLongdustScore("ACGT", gcContent: 0.0));
        });
    }

    #endregion
}
