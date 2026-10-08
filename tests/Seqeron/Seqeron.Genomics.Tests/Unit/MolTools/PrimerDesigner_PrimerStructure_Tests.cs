namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Canonical tests for PRIMER-STRUCT-001: Primer Structure Analysis.
/// Tests for secondary structure detection in PCR primers.
/// 
/// Methods under test:
/// - PrimerDesigner.FindLongestHomopolymer(seq)
/// - PrimerDesigner.FindLongestDinucleotideRepeat(seq)
/// - PrimerDesigner.HasHairpinPotential(seq, minStemLength)
/// - PrimerDesigner.HasPrimerDimer(primer1, primer2, minComplementarity)
/// - PrimerDesigner.Calculate3PrimeStability(seq)
/// 
/// Evidence sources:
/// - Wikipedia (Primer, Primer dimer, Stem-loop, Nucleic acid thermodynamics)
/// - Primer3 Manual (primer3.org)
/// - SantaLucia (1998) PNAS 95:1460-65
/// </summary>
[TestFixture]
public class PrimerDesigner_PrimerStructure_Tests
{
    #region FindLongestHomopolymer Tests

    /// <summary>
    /// Empty sequence should return 0.
    /// Source: Standard null/empty handling.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_EmptySequence_ReturnsZero()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("");
        Assert.That(result, Is.EqualTo(0));
    }

    /// <summary>
    /// Null sequence should return 0.
    /// Source: Standard null handling.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_NullSequence_ReturnsZero()
    {
        int result = PrimerDesigner.FindLongestHomopolymer(null!);
        Assert.That(result, Is.EqualTo(0));
    }

    /// <summary>
    /// Sequence with no runs (all different bases) should return 1.
    /// Source: Primer3 PRIMER_MAX_POLY_X behavior.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_NoRun_ReturnsOne()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("ACGT");
        Assert.That(result, Is.EqualTo(1));
    }

    /// <summary>
    /// Sequence with internal homopolymer run returns run length.
    /// Source: Primer3 PRIMER_MAX_POLY_X.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_InternalRun_ReturnsRunLength()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("ACAAAAGT");
        Assert.That(result, Is.EqualTo(4)); // AAAA
    }

    /// <summary>
    /// All same nucleotide returns full length.
    /// Source: Primer3 PRIMER_MAX_POLY_X.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_AllSame_ReturnsFullLength()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("AAAAAA");
        Assert.That(result, Is.EqualTo(6));
    }

    /// <summary>
    /// Case-insensitive matching for homopolymer detection.
    /// Source: Universal DNA sequence handling convention.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_MixedCase_IsCaseInsensitive()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("AaAaAa");
        Assert.That(result, Is.EqualTo(6));
    }

    /// <summary>
    /// Homopolymer at end of sequence is detected.
    /// Source: Edge case verification.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_RunAtEnd_ReturnsRunLength()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("ACGTTTTT");
        Assert.That(result, Is.EqualTo(5)); // TTTTT at end
    }

    /// <summary>
    /// Multiple runs returns the longest.
    /// Source: Algorithm correctness.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_MultipleRuns_ReturnsLongest()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("AAACCCCCGG");
        Assert.That(result, Is.EqualTo(5)); // CCCCC is longest
    }

    /// <summary>
    /// N is a worst-case wildcard, exactly as Primer3's _pr_violates_poly_x (libprimer3.cc; its
    /// header comment lists NNG 3, ANA 3, CNN 3, TNNG 3, GNGNG 5, ANGNG 4). The 18–21-mers were
    /// confirmed with primer3-py 2.3.1 check_primers (PRIMER_MAX_NS_ACCEPTED 10): the smallest
    /// PRIMER_MAX_POLY_X that does not reject them equals the value below.
    /// </summary>
    [TestCase("NNG", 3)]
    [TestCase("ANA", 3)]
    [TestCase("CNN", 3)]
    [TestCase("TNNG", 3)]
    [TestCase("GNGNG", 5)]
    [TestCase("ANGNG", 4)]
    [TestCase("CAGTCAGTCANGNGTCAGTC", 4)]
    [TestCase("CAGTCAGTCAGNGNGTCAGTC", 5)]
    [TestCase("CAGTCAGTCAANAACAGTC", 5)]
    [TestCase("NNGTCAGTCAGTCAGTCAG", 3)]
    [TestCase("CAGTCAGTCAGTCAGTCNN", 3)]
    public void FindLongestHomopolymer_NIsWorstCaseWildcard_MatchesPrimer3PolyX(string sequence, int expected)
    {
        Assert.That(PrimerDesigner.FindLongestHomopolymer(sequence), Is.EqualTo(expected));
    }

    #endregion

    #region FindLongestDinucleotideRepeat Tests

    /// <summary>
    /// Null, empty, or short (&lt;4 bp) sequence returns 0.
    /// Source: Implementation bounds (need at least 2 dinucleotide units = 4 bp).
    /// </summary>
    [TestCase(null, 0)]
    [TestCase("", 0)]
    [TestCase("ACG", 0)]
    public void FindLongestDinucleotideRepeat_InvalidInput_ReturnsZero(string? sequence, int expected)
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat(sequence!);
        Assert.That(result, Is.EqualTo(expected));
    }

    /// <summary>
    /// Sequence with no dinucleotide repeats returns 1.
    /// Source: Primer3 behavior — single dinucleotide = 1 repeat unit.
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_NoRepeat_ReturnsOne()
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("ACGT");
        Assert.That(result, Is.EqualTo(1));
    }

    /// <summary>
    /// ACACACAC contains 4 AC repeats.
    /// Source: Primer3 behavior for dinucleotide repeats.
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_AcRepeat_ReturnsCount()
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("ACACACACG");
        Assert.That(result, Is.EqualTo(4)); // ACACACAC = 4 x AC
    }

    /// <summary>
    /// AT repeat pattern is detected.
    /// Source: Common microsatellite pattern.
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_AtRepeat_ReturnsCount()
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("ATATATAT");
        Assert.That(result, Is.EqualTo(4)); // ATATATAT = 4 x AT
    }

    /// <summary>
    /// Returns longest dinucleotide repeat when multiple exist.
    /// Source: Algorithm correctness.
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_MultipleRepeats_ReturnsLongest()
    {
        // "ACACGCGCGCGC" = ACAC (2 AC) + GCGCGCGC (4 GC)
        // Implementation counts the number of times the 2-base pattern repeats
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("ACACGCGCGCGC");
        Assert.That(result, Is.EqualTo(4)); // GCGCGCGC = 4 x GC (8 bases / 2)
    }

    #endregion

    #region HasHairpinPotential Tests

    /// <summary>
    /// Null, empty, or too-short sequences cannot form hairpin.
    /// Source: Wikipedia Stem-loop (minimum structure = 2×stem + loop).
    /// </summary>
    [TestCase(null)]
    [TestCase("")]
    [TestCase("ACGT")]
    [TestCase("ACGTACGTAC")] // 10 bp < 2×4+3=11 for default minStemLength=4
    public void HasHairpinPotential_InvalidOrTooShort_ReturnsFalse(string? sequence)
    {
        bool result = PrimerDesigner.HasHairpinPotential(sequence!);
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Non-self-complementary sequence cannot form hairpin.
    /// Source: Wikipedia Stem-loop (requires complementary regions).
    /// </summary>
    [Test]
    public void HasHairpinPotential_NonSelfComplementary_ReturnsFalse()
    {
        // All A's cannot form complementary stems
        bool result = PrimerDesigner.HasHairpinPotential("AAAACCCCAAAA");
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Self-complementary sequence can form hairpin.
    /// Source: Wikipedia Stem-loop.
    /// </summary>
    [Test]
    public void HasHairpinPotential_SelfComplementary_ReturnsTrue()
    {
        // ACGT reversed = TGCA, which is complementary to ACGT
        bool result = PrimerDesigner.HasHairpinPotential("ACGTACGTACGT");
        Assert.That(result, Is.True);
    }

    /// <summary>
    /// Custom minStemLength is respected.
    /// Source: API contract.
    /// </summary>
    [Test]
    public void HasHairpinPotential_CustomMinStem_RespectsParameter()
    {
        // With minStemLength=6, needs 2×6+3=15 bases minimum
        // 12 bases should return false
        bool result = PrimerDesigner.HasHairpinPotential("ACGTACGTACGT", minStemLength: 6);
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Custom minLoopLength is respected.
    /// Source: Wikipedia Stem-loop — minimum 3 nt loop is steric constraint.
    /// </summary>
    [Test]
    public void HasHairpinPotential_CustomMinLoopLength_RespectsParameter()
    {
        // GCGCTTTTGCGC: stem=4 (GCGC), loop=4 (TTTT)
        // With minLoopLength=5, the 4-nt loop is insufficient (and 12 < 4*2+5=13)
        bool withLoop5 = PrimerDesigner.HasHairpinPotential("GCGCTTTTGCGC", minStemLength: 4, minLoopLength: 5);
        bool withLoop3 = PrimerDesigner.HasHairpinPotential("GCGCTTTTGCGC", minStemLength: 4, minLoopLength: 3);

        Assert.Multiple(() =>
        {
            Assert.That(withLoop5, Is.False, "Loop=4 < minLoopLength=5 → no hairpin");
            Assert.That(withLoop3, Is.True, "Loop=4 ≥ minLoopLength=3 → hairpin detected");
        });
    }

    /// <summary>
    /// M10b — <c>HasHairpinPotential</c> is a sequence-only library screen (default stem 4 unsourced, loop 3 =
    /// Primer3 thal.c min_hrpn_loop), not Primer3's hairpin screen: it flags both sequences although Primer3's
    /// ntthal hairpin (PRIMER_HAIRPIN_TH, limit 47 °C) accepts them.
    /// Source: primer3-py 2.3.1 <c>calc_hairpin</c> at 50 mM Na⁺ / 1.5 mM Mg²⁺ / 0.6 mM dNTP / 50 nM:
    /// AAAACCCTTTT → structure_found False (Tm 0); CAGTAAAACCCTTTTGCAGC → Tm 37.65 °C.
    /// </summary>
    [TestCase("AAAACCCTTTT", 0.0)]
    [TestCase("CAGTAAAACCCTTTTGCAGC", 37.65)]
    public void HasHairpinPotential_LibraryScreen_DiffersFromPrimer3NtthalHairpin(string sequence, double primer3HairpinTm)
    {
        var primer3 = PrimerDesigner.CalculatePrimer3OligoStructure(sequence);
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.HasHairpinPotential(sequence), Is.True, "library screen flags the exact 4-bp stem");
            Assert.That(primer3, Is.Not.Null);
            Assert.That(primer3!.Value.HairpinTh, Is.EqualTo(primer3HairpinTm).Within(0.005), "Primer3 ntthal hairpin Tm");
            Assert.That(primer3.Value.HairpinTh, Is.LessThanOrEqualTo(PrimerDesigner.Primer3MaxStructureTm),
                "Primer3 PRIMER_MAX_HAIRPIN_TH accepts it");
        });
    }

    /// <summary>
    /// Long sequence (>100bp) uses suffix tree optimization.
    /// Source: Performance optimization test.
    /// </summary>
    [Test]
    public void HasHairpinPotential_LongSequence_UsesSuffixTreeOptimization()
    {
        // Create 150bp sequence with hairpin potential
        // ACGT...ACGT pattern at start and end (complementary when reversed)
        var sb = new System.Text.StringBuilder();
        sb.Append("ACGTACGTACGT"); // 12bp stem region
        sb.Append(new string('A', 126)); // spacer (loop + filler)
        sb.Append("ACGTACGTACGT"); // 12bp complementary region
        string longSeq = sb.ToString(); // 150bp total

        // Should detect hairpin using suffix tree (>100bp threshold)
        bool result = PrimerDesigner.HasHairpinPotential(longSeq);
        Assert.That(result, Is.True);
    }

    /// <summary>
    /// Long sequence without hairpin returns false.
    /// Source: Performance optimization test.
    /// </summary>
    /// <summary>
    /// A5-2: only A·T / G·C pairs form a stem at every length. The ≥ 100-nt suffix-tree path used the IUPAC reverse
    /// complement (N→N), so two NNNN blocks counted as a stem there but not in the &lt; 100-nt scan.
    /// </summary>
    [TestCase(97)]
    [TestCase(113)]
    public void HasHairpinPotential_NonAcgtBlocks_NeverPair_OnBothPaths(int length)
    {
        string flank = string.Concat(Enumerable.Repeat("ACCA", (length - 13) / 8));
        string seq = flank + "NNNNACCACNNNN" + string.Concat(Enumerable.Repeat("CAAC", (length - 13 - flank.Length) / 4));
        Assert.Multiple(() =>
        {
            Assert.That(seq.Length, Is.EqualTo(length));
            Assert.That(PrimerDesigner.HasHairpinPotential(seq), Is.False, "A/C/N only: no Watson–Crick stem");
            Assert.That(PrimerDesigner.HasHairpinPotential(seq.Replace("NNNNACCACNNNN", "GGTTACCACAACC")), Is.True,
                "control: GGTT·AACC stem with a 5-nt loop");
        });
    }

    /// <summary>
    /// A6-1: the canonical reverse complement maps U→A, so before the fix a U stem window became "AAAA", passed the
    /// pattern-only non-ACGT guard and matched an A run on the ≥ 100-nt suffix-tree path, while the &lt; 100-nt scan
    /// (strict A·T / G·C) never pairs U. Lowercase input is upper-cased before either path, so it pairs identically.
    /// </summary>
    [TestCase("AAAAACCCCUUUUU", false, TestName = "U_never_pairs_with_A")]
    [TestCase("AAAAACCCCuuuuu", false, TestName = "lowercase_u_never_pairs_with_A")]
    [TestCase("AAAAACCCCTTTTT", true, TestName = "control_T_pairs_with_A")]
    [TestCase("aaaaaccccttttt", true, TestName = "lowercase_pairs_like_uppercase")]
    [TestCase("AAAAACCCCCCCCC", false, TestName = "control_no_stem")]
    public void HasHairpinPotential_UracilAndLowercase_SameOnBothPaths(string core, bool expected)
    {
        string longSeq = core + string.Concat(Enumerable.Repeat("AC", 60)); // 134 nt → suffix-tree path
        Assert.Multiple(() =>
        {
            Assert.That(core.Length, Is.LessThan(100));
            Assert.That(longSeq.Length, Is.GreaterThanOrEqualTo(100));
            Assert.That(PrimerDesigner.HasHairpinPotential(core), Is.EqualTo(expected), "< 100-nt scan");
            Assert.That(PrimerDesigner.HasHairpinPotential(longSeq), Is.EqualTo(expected), "≥ 100-nt suffix tree");
        });
    }

    /// <summary>
    /// A6-1: randomized equivalence of the two scan paths. Padding a &lt; 100-nt sequence with N (which never pairs on
    /// either path) up to ≥ 100 nt switches it to the suffix-tree path without adding or removing any A·T / G·C stem,
    /// so the result must be unchanged — over an alphabet that includes U, N, an IUPAC code and lowercase.
    /// </summary>
    [Test]
    public void HasHairpinPotential_ShortScanAndSuffixTree_AgreeOnRandomSequences()
    {
        const string alphabet = "ACGTACGTUNSacgtu";
        var rng = new Random(20261008);
        int positives = 0;
        for (int trial = 0; trial < 2000; trial++)
        {
            int len = rng.Next(11, 99);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            string s = new(chars);
            int stem = rng.Next(3, 6), loop = rng.Next(3, 5);
            bool shortResult = PrimerDesigner.HasHairpinPotential(s, stem, loop);
            bool longResult = PrimerDesigner.HasHairpinPotential(s + new string('N', 100), stem, loop);
            Assert.That(longResult, Is.EqualTo(shortResult), $"seq={s} stem={stem} loop={loop}");
            if (shortResult) positives++;
        }
        Assert.That(positives, Is.InRange(200, 1800), "both outcomes are exercised");
    }

    [Test]
    public void HasHairpinPotential_LongSequenceNoHairpin_ReturnsFalse()
    {
        // All A's cannot form hairpin (A is complementary to T, not A)
        string longSeq = new string('A', 150);
        bool result = PrimerDesigner.HasHairpinPotential(longSeq);
        Assert.That(result, Is.False);
    }

    #endregion

    #region HasPrimerDimer Tests

    /// <summary>
    /// Null or empty primer returns false.
    /// Source: Standard null guard.
    /// </summary>
    [TestCase(null, "ACGT")]
    [TestCase("", "ACGT")]
    [TestCase("ACGT", null)]
    [TestCase("ACGT", "")]
    public void HasPrimerDimer_NullOrEmptyPrimer_ReturnsFalse(string? p1, string? p2)
    {
        bool result = PrimerDesigner.HasPrimerDimer(p1!, p2!);
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Primers whose 3' ends cannot pair do not form dimers.
    /// Source: Primer3 alignment-mode PRIMER_PAIR_COMPL_END (libprimer3.cc characterize_pair, dpal.c
    /// DPAL_GLOBAL_END). primer3-py 2.3.1 check_primers (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0):
    /// AAAAAAAA + AAAACCCC → PRIMER_PAIR_0_COMPL_END = 0.0.
    /// </summary>
    [Test]
    public void HasPrimerDimer_NonComplementary3Ends_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimerDimerEndComplementarity("AAAAAAAA", "AAAACCCC"), Is.EqualTo(0.0));
            Assert.That(PrimerDesigner.HasPrimerDimer("AAAAAAAA", "AAAACCCC"), Is.False);
        });
    }

    /// <summary>
    /// Regression (DUP_MAP §10): identical poly-A primers cannot pair (A·A is not a base pair), so
    /// they are not a primer-dimer. The former implementation compared primer1's 3' window with
    /// revcomp(primer2) by complementarity (parallel pairing) and flagged A₈/A₈.
    /// Source: Primer3 compl_end = align(A₈, revcomp(A₈)) = 0 (dpal.c compiled from the primer3
    /// source; PRIMER_LEFT_0_SELF_END of AAAAAAAA = 0.0 in primer3-py check_primers);
    /// primer3-py calc_heterodimer(A₂₀, A₂₀) → structure_found = False.
    /// </summary>
    [Test]
    public void HasPrimerDimer_IdenticalPolyA_IsNotADimer()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimerDimerEndComplementarity("AAAAAAAA", "AAAAAAAA"), Is.EqualTo(0.0));
            Assert.That(PrimerDesigner.HasPrimerDimer("AAAAAAAA", "AAAAAAAA"), Is.False);
            Assert.That(PrimerDesigner.CalculateDimerThermodynamicsNtthal(
                new string('A', 20), new string('A', 20)), Is.Null);
        });
    }

    /// <summary>
    /// Primers with complementary 3' ends form dimers: A₈ and T₈ pair over all 8 bases.
    /// Source: Primer3 compl_end (dpal.c compiled from source) = 8.
    /// </summary>
    [Test]
    public void HasPrimerDimer_Complementary3Ends_ReturnsTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimerDimerEndComplementarity("AAAAAAAA", "TTTTTTTT"), Is.EqualTo(8.0));
            Assert.That(PrimerDesigner.HasPrimerDimer("AAAAAAAA", "TTTTTTTT"), Is.True);
        });
    }

    /// <summary>
    /// Regression (DUP_MAP §10): 3' ends ...GGCC / ...GGCC pair through a 4-base offset overlap
    /// (GGCC is its own reverse complement), which a fixed full-window comparison misses.
    /// Source: primer3-py check_primers (alignment mode): TTCAGTCAGTCAGTGGCC + ACTGACTGACTGAGGCC →
    /// PRIMER_PAIR_0_COMPL_END = 4.0 (> PRIMER_PAIR_MAX_COMPL_END 3.00: "high end compl");
    /// ntthal END1 Tm 10.94 °C (calc_end_stability) confirms the 3'-anchored duplex.
    /// </summary>
    [Test]
    public void HasPrimerDimer_SelfComplementaryGgccEnds_Detected()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimerDimerEndComplementarity("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC"),
                Is.EqualTo(4.0));
            Assert.That(PrimerDesigner.HasPrimerDimer("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC"), Is.True);
        });
    }

    /// <summary>
    /// Primer3 compl_end values for primer pairs (primer3-py 2.3.1 check_primers,
    /// PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0 → PRIMER_PAIR_0_COMPL_END).
    /// </summary>
    [TestCase("AACCGGTTAACCATCGATCG", "AACCGGTTAAGCTAGCTA", 1.0)]
    [TestCase("AACCGGTTAACCATCGATCG", "AACCGGTTAACGATCGAT", 8.0)]
    [TestCase("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC", 4.0)]
    [TestCase("AAAAAAAA", "AAAACCCC", 0.0)]
    public void CalculatePrimerDimerEndComplementarity_MatchesPrimer3ComplEnd(string p1, string p2, double expected)
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimerDimerEndComplementarity(p1, p2), Is.EqualTo(expected));
            // Primer3 takes the max over both orientations, so the score is symmetric.
            Assert.That(PrimerDesigner.CalculatePrimerDimerEndComplementarity(p2, p1), Is.EqualTo(expected));
        });
    }

    /// <summary>
    /// Primer3 self_end values (primer3-py check_primers, alignment mode, PRIMER_LEFT_0_SELF_END /
    /// PRIMER_RIGHT_0_SELF_END).
    /// </summary>
    [TestCase("AACCGGTTAACCATCGATCG", 6.0)]
    [TestCase("TTCAGTCAGTCAGTGGCC", 4.0)]
    [TestCase("AAAAAAAA", 0.0)]
    [TestCase("AAAACCCC", 0.0)]
    public void CalculatePrimerSelfEndComplementarity_MatchesPrimer3SelfEnd(string primer, double expected)
    {
        Assert.That(PrimerDesigner.CalculatePrimerSelfEndComplementarity(primer), Is.EqualTo(expected));
    }

    /// <summary>
    /// Custom minComplementarity is respected: the self-complementary palindrome ACGTACGT pairs
    /// with itself over all 8 bases (Primer3 compl_end = 8, dpal.c compiled from source).
    /// </summary>
    [Test]
    public void HasPrimerDimer_CustomMinComplementarity_RespectsParameter()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.HasPrimerDimer("ACGTACGT", "ACGTACGT", minComplementarity: 8), Is.True);
            Assert.That(PrimerDesigner.HasPrimerDimer("ACGTACGT", "ACGTACGT", minComplementarity: 9), Is.False);
        });
    }

    #endregion

    #region Calculate3PrimeStability Tests

    /// <summary>
    /// Null or empty sequence returns 0.
    /// </summary>
    [TestCase(null, 0.0)]
    [TestCase("", 0.0)]
    public void Calculate3PrimeStability_InvalidInput_ReturnsZero(string? sequence, double expected)
    {
        double result = PrimerDesigner.Calculate3PrimeStability(sequence!);
        Assert.That(result, Is.EqualTo(expected));
    }

    /// <summary>
    /// Primer3 end_oligodg(seq, 5, santalucia) (oligotm.c, compiled from the primer3 source and
    /// run on these inputs; Primer3 reports −ΔG, so the expected ΔG is the negated output). Covers
    /// primers shorter than 5 (the whole primer is scored, including the +0.43 symmetry term for
    /// the self-complementary GC), N (Primer3's N parameters, no A·T penalty), and a 20-mer
    /// (end_oligodg 3.25).
    /// </summary>
    [TestCase("ACGT", -2.56)]
    [TestCase("GC", 0.15)]
    [TestCase("A", 2.06)]
    [TestCase("AT", 1.61)]
    [TestCase("ACGTN", -3.62)]
    [TestCase("NNNNN", -0.36)]
    [TestCase("AAANA", -1.40)]
    [TestCase("GATCGAGGACTGCCTTGGTA", -3.25)]
    public void Calculate3PrimeStability_MatchesPrimer3EndOligoDg(string sequence, double expected)
    {
        Assert.That(PrimerDesigner.Calculate3PrimeStability(sequence), Is.EqualTo(expected).Within(1e-9));
    }

    /// <summary>
    /// A character other than A/C/G/T/N in the 3' window is Primer3's OLIGOTM_ERROR (NaN here);
    /// characters before the window are irrelevant.
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_InvalidCharacterInWindow_ReturnsNaN()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.Calculate3PrimeStability("ACGU"), Is.NaN);
            Assert.That(PrimerDesigner.Calculate3PrimeStability("XXXXXGCGCG"), Is.EqualTo(-6.86).Within(1e-9));
        });
    }

    /// <summary>
    /// Exact 5-base input produces correct ΔG.
    /// Source: SantaLucia (1998) + Primer3 Manual.
    /// TACGT = TA(-0.58) + AC(-1.44) + CG(-2.17) + GT(-1.44) + Init(A·T)(+1.03) + Init(A·T)(+1.03) = -3.57
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_Exact5Bases_ProducesCorrectDeltaG()
    {
        double result = PrimerDesigner.Calculate3PrimeStability("TACGT");
        Assert.That(result, Is.EqualTo(-3.57).Within(0.01));
    }

    /// <summary>
    /// GC-rich 3' end is more stable (more negative ΔG) than AT-rich — with exact values.
    /// Source: SantaLucia (1998) + Primer3 Manual.
    /// GCGCG = -6.86, TATAT = -0.86.
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_GcRich_MoreNegativeThanAtRich()
    {
        double gcRich = PrimerDesigner.Calculate3PrimeStability("ACGTGCGCG"); // ends with GCGCG
        double atRich = PrimerDesigner.Calculate3PrimeStability("ACGTATATAT"); // ends with TATAT

        Assert.Multiple(() =>
        {
            Assert.That(gcRich, Is.EqualTo(-6.86).Within(0.01));
            Assert.That(atRich, Is.EqualTo(-0.86).Within(0.01));
            Assert.That(gcRich, Is.LessThan(atRich));
        });
    }

    /// <summary>
    /// Case insensitive: upper and lower case return identical ΔG.
    /// Source: Universal DNA convention. Verified with exact GCGCG = -6.86.
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_MixedCase_ReturnsSameExactValue()
    {
        double upper = PrimerDesigner.Calculate3PrimeStability("ACGTGCGCG");
        double lower = PrimerDesigner.Calculate3PrimeStability("acgtgcgcg");

        Assert.Multiple(() =>
        {
            Assert.That(upper, Is.EqualTo(-6.86).Within(0.01));
            Assert.That(lower, Is.EqualTo(upper));
        });
    }

    /// <summary>
    /// GCGCG (most stable 5mer) produces exact Primer3 reference value.
    /// Source: Primer3 Manual PRIMER_MAX_END_STABILITY with SantaLucia (1998) parameters.
    /// GCGCG = GC(-2.24) + CG(-2.17) + GC(-2.24) + CG(-2.17) + Init(G·C)(+0.98) + Init(G·C)(+0.98) = -6.86
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_MostStable5mer_MatchesPrimer3()
    {
        double result = PrimerDesigner.Calculate3PrimeStability("AAAAAGCGCG");

        // Primer3 Manual: "most stable 5mer duplex = 6.86 kcal/mol (GCGCG)"
        Assert.That(result, Is.EqualTo(-6.86).Within(0.01));
    }

    /// <summary>
    /// TATAT (least stable 5mer) produces exact Primer3 reference value.
    /// Source: Primer3 Manual PRIMER_MAX_END_STABILITY with SantaLucia (1998) parameters.
    /// TATAT = TA(-0.58) + AT(-0.88) + TA(-0.58) + AT(-0.88) + Init(A·T)(+1.03) + Init(A·T)(+1.03) = -0.86
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_LeastStable5mer_MatchesPrimer3()
    {
        double result = PrimerDesigner.Calculate3PrimeStability("GGGGGTATAT");

        // Primer3 Manual: "most labile 5mer duplex = 0.86 kcal/mol (TATAT)"
        Assert.That(result, Is.EqualTo(-0.86).Within(0.01));
    }

    #endregion

    #region Integration Tests

    /// <summary>
    /// Well-designed primer produces expected exact metrics.
    /// Source: Primer3 primer evaluation workflow.
    /// ACGTACGTACGTACGTACGT: homopolymer=1, dinuc repeat=1 (no repeats),
    /// hairpin=true (ACGT pattern has reverse-complement matches),
    /// last 5 = TACGT → ΔG = -3.57 kcal/mol.
    /// </summary>
    [Test]
    public void PrimerStructureAnalysis_WellDesignedPrimer_ExactMetrics()
    {
        const string primer = "ACGTACGTACGTACGTACGT"; // 20 bp

        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindLongestHomopolymer(primer), Is.EqualTo(1));
            Assert.That(PrimerDesigner.FindLongestDinucleotideRepeat(primer), Is.EqualTo(1));
            Assert.That(PrimerDesigner.Calculate3PrimeStability(primer), Is.EqualTo(-3.57).Within(0.01));
        });
    }

    /// <summary>
    /// Problematic primer (20 G's) produces exact metrics showing all issues.
    /// Source: Primer3 PRIMER_MAX_POLY_X, SantaLucia (1998).
    /// GGGGG = GG(-1.84)×4 + Init(G·C)(+0.98)×2 = -5.40 kcal/mol.
    /// </summary>
    [Test]
    public void PrimerStructureAnalysis_ProblematicPrimer_ExactMetrics()
    {
        const string badPrimer = "GGGGGGGGGGGGGGGGGGGG"; // 20 G's

        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindLongestHomopolymer(badPrimer), Is.EqualTo(20));
            Assert.That(PrimerDesigner.Calculate3PrimeStability(badPrimer), Is.EqualTo(-5.40).Within(0.01));
        });
    }

    #endregion

    #region Primer3 thermodynamic structure screen (ntthal)

    /// <summary>
    /// Primer3's default (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1) per-primer and pair structure values.
    /// Source: primer3-py 2.3.1 design_primers, PRIMER_TASK=check_primers, default conditions
    /// (50 mM monovalent, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM): PRIMER_{LEFT,RIGHT}_0_SELF_ANY_TH,
    /// _SELF_END_TH, _HAIRPIN_TH, PRIMER_PAIR_0_COMPL_ANY_TH, PRIMER_PAIR_0_COMPL_END_TH.
    /// </summary>
    private static readonly object[] Primer3StructureCases =
    {
        new object[] { "GGGGAAAACCCCATATGCAG", "CTGCATATGGGGTTTTCCCA",
            new[] { 18.30105741642683, 0.0, 58.09337111624046, 6.200192585916625, 6.200192585916625, 66.14453842604473, 57.51407534362289, 57.51407534362289 } },
        new object[] { "TTCAGTCAGTCAGTAAAAGGGG", "ACTGACTGACTGACCCCTTTT",
            new[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 42.22316686789662, 43.12474243211591 } },
        new object[] { "ACGTACGTACGTACGTACGT", "TGCAGCATGCATGCAATGCA",
            new[] { 60.49048111491186, 60.49048111491186, 71.76211301142945, 41.76336394643283, 41.76336394643283, 76.36745306969078, 0.0, 0.0 } },
        new object[] { "TGTCGGAAAACAGCCAGGCTAACG", "CGGTTGATCGTGAGTCAAATTCAGC",
            new[] { 15.644106552541643, 12.3535987107486, 36.400151579667295, 0.0, 0.35225774955779343, 37.336259822139596, 0.0, 0.0 } },
    };

    [TestCaseSource(nameof(Primer3StructureCases))]
    public void Primer3ThermodynamicStructure_MatchesPrimer3CheckPrimers(string left, string right, double[] expected)
    {
        var l = PrimerDesigner.CalculatePrimer3OligoStructure(left)!.Value;
        var r = PrimerDesigner.CalculatePrimer3OligoStructure(right)!.Value;
        var pair = PrimerDesigner.CalculatePrimer3PairComplementarity(left, right)!.Value;
        double[] actual = { l.SelfAnyTh, l.SelfEndTh, l.HairpinTh, r.SelfAnyTh, r.SelfEndTh, r.HairpinTh, pair.ComplAnyTh, pair.ComplEndTh };
        Assert.That(actual, Is.EqualTo(expected).Within(1e-9));
    }

    /// <summary>
    /// ntthal END1/END2 alignment types and the divalent/dNTP salt term.
    /// Source: primer3-py 2.3.1 calc_end_stability(a, b) (END1; END2(a, b) = END1(b, a)) and
    /// calc_heterodimer at mv 50, dv 1.5, dntp 0.6, dna 50 nM.
    /// </summary>
    [TestCase("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC", PrimerDesigner.NtthalAlignmentMode.End1, 10.942030550739162, -5.051857560300428)]
    [TestCase("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC", PrimerDesigner.NtthalAlignmentMode.End2, 13.481280192208146, -5.668962560300446)]
    [TestCase("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC", PrimerDesigner.NtthalAlignmentMode.Any, 41.167106798577606, -12.54145474705899)]
    [TestCase("TTCAGTCAGTCAGTAAAAGGGG", "ACTGACTGACTGACCCCTTTT", PrimerDesigner.NtthalAlignmentMode.End1, 22.803405505733053, -8.035343602454173)]
    public void CalculateDimerThermodynamicsNtthal_AlignmentModes_MatchPrimer3Py(
        string a, string b, PrimerDesigner.NtthalAlignmentMode mode, double tm, double dg)
    {
        var d = PrimerDesigner.CalculateDimerThermodynamicsNtthal(a, b, mode, 0.050, 0.0015, 0.0006, 50e-9)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(d.TmCelsius, Is.EqualTo(tm).Within(1e-9));
            Assert.That(d.DeltaG37, Is.EqualTo(dg).Within(1e-9));
        });
    }

    /// <summary>
    /// Divalent cations enter ntthal hairpins only through saltCorrectS.
    /// Source: primer3-py 2.3.1 calc_hairpin('GGGGAAAACCCCATATGCAG') = 58.09337111624046 °C at
    /// dv 1.5/dntp 0.6 and 54.28830868675732 °C at dv = dntp = 0.
    /// </summary>
    [Test]
    public void CalculateHairpinThermodynamicsNtthal_Divalent_MatchesPrimer3Py()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGGAAAACCCCATATGCAG", 0.050, 0.0015, 0.0006)!.Value.TmCelsius,
                Is.EqualTo(58.09337111624046).Within(1e-9));
            Assert.That(PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGGAAAACCCCATATGCAG", 0.050, 0.0, 0.0)!.Value.TmCelsius,
                Is.EqualTo(54.28830868675732).Within(1e-9));
        });
    }

    /// <summary>
    /// Poly-A cannot form any structure: every Primer3 structure value is 0 (align_thermod).
    /// Non-ACGT input has no ntthal structure (null).
    /// </summary>
    [Test]
    public void Primer3ThermodynamicStructure_NoStructureAndInvalidInput()
    {
        var a = PrimerDesigner.CalculatePrimer3OligoStructure(new string('A', 20))!.Value;
        var pair = PrimerDesigner.CalculatePrimer3PairComplementarity(new string('A', 20), new string('A', 20))!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(a, Is.EqualTo(new PrimerDesigner.Primer3OligoStructure(0, 0, 0)));
            Assert.That(pair, Is.EqualTo(new PrimerDesigner.Primer3PairComplementarity(0, 0)));
            Assert.That(PrimerDesigner.CalculatePrimer3OligoStructure("ACGTNACGT"), Is.Null);
            Assert.That(PrimerDesigner.CalculatePrimer3PairComplementarity("ACGT", ""), Is.Null);
        });
    }

    /// <summary>
    /// EvaluatePrimer applies Primer3's default limits (47 °C): ACGTACGTACGTACGTACGT has
    /// self-any 60.49, self-end 60.49 and hairpin 71.76 °C (primer3-py check_primers), so it is
    /// rejected for all three; the heuristic screen reports only the stem-loop and no Tm values.
    /// </summary>
    [Test]
    public void EvaluatePrimer_ThermodynamicScreen_ReportsPrimer3Values()
    {
        var param = PrimerDesigner.DefaultParameters with { MinGcContent = 0, MaxGcContent = 100, MinTm = 0, MaxTm = 100 };
        var thermo = PrimerDesigner.EvaluatePrimer("ACGTACGTACGTACGTACGT", 0, true, param);
        var heur = PrimerDesigner.EvaluatePrimer("ACGTACGTACGTACGTACGT", 0, true,
            param with { StructureScreen = PrimerStructureScreen.Heuristic });
        var relaxed = PrimerDesigner.EvaluatePrimer("ACGTACGTACGTACGTACGT", 0, true, param with { MaxStructureTm = 100 });
        Assert.Multiple(() =>
        {
            Assert.That(thermo.SelfAnyTh, Is.EqualTo(60.49048111491186).Within(1e-9));
            Assert.That(thermo.SelfEndTh, Is.EqualTo(60.49048111491186).Within(1e-9));
            Assert.That(thermo.HairpinTh, Is.EqualTo(71.76211301142945).Within(1e-9));
            Assert.That(thermo.HasHairpin, Is.True);
            Assert.That(thermo.IsValid, Is.False);
            Assert.That(thermo.Issues.Count(i => i.Contains("Primer3 PRIMER_MAX_")), Is.EqualTo(3));
            Assert.That(heur.SelfAnyTh, Is.Null);
            Assert.That(heur.HasHairpin, Is.True);
            Assert.That(heur.Issues, Has.Some.EqualTo("Potential hairpin structure detected"));
            Assert.That(relaxed.HasHairpin, Is.False);
            Assert.That(relaxed.Issues.Any(i => i.Contains("Primer3 PRIMER_MAX_")), Is.False);
        });
    }

    #endregion
}
