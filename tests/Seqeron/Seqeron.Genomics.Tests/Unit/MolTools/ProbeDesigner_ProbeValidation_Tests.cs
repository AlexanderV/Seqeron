namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for ThermodynamicProbeDesigner.ValidateProbe and CheckSpecificity.
/// Test Unit: PROBE-VALID-001
/// 
/// Evidence Sources:
/// - Wikipedia: Hybridization probe (cross-hybridization, stringency)
/// - Wikipedia: DNA microarray (probe specificity)
/// - Wikipedia: Off-target genome editing (mismatch tolerance 1-5 bp)
/// - Wikipedia: BLAST (approximate matching algorithms)
/// </summary>
[TestFixture]
public class ProbeDesigner_ProbeValidation_Tests
{

    // Standard probe for validation tests
    private const string StandardProbe = "ACGTACGTACGTACGTACGT";

    // Self-complementary (palindromic) probe
    private const string PalindromicProbe = "GCGCGCGCGCGCGCGCGCGC";

    // Unique probe that appears once in reference
    private const string UniqueProbe = "ATCGATCGATCGATCGATCG";

    // Reference containing the unique probe once
    private static readonly string[] SingleMatchReference = new[]
    {
        "NNNNNATCGATCGATCGATCGATCGNNNN"
    };

    // Reference with repeated sequence (multiple matches)
    private static readonly string[] MultipleMatchReference = new[]
    {
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    };



    [Test]
    public void CheckSpecificity_UniqueSequence_ReturnsOne()
    {
        // M10: Unique match returns 1.0
        string genome = "NNNNNATCGATCGATCGATCGATCGNNNN";
        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        double specificity = ProbeDesigner.CheckSpecificity(UniqueProbe, genomeIndex);

        Assert.That(specificity, Is.EqualTo(1.0),
            "Unique probe should have specificity 1.0");
    }

    [Test]
    public void CheckSpecificity_MultipleOccurrences_ReturnsOneOverCount()
    {
        // M11: Multiple matches returns 1.0/count
        string repeatedSequence = "ACGT";
        string genome = "ACGTACGTACGTACGT"; // Contains ACGT 4 times
        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        double specificity = ProbeDesigner.CheckSpecificity(repeatedSequence, genomeIndex);
        var positions = genomeIndex.FindAllOccurrences(repeatedSequence);

        Assert.Multiple(() =>
        {
            Assert.That(positions.Count, Is.GreaterThan(1),
                "Should find multiple occurrences");
            Assert.That(specificity, Is.EqualTo(1.0 / positions.Count).Within(0.001),
                "Specificity should be 1.0 / count");
        });
    }

    [Test]
    public void CheckSpecificity_NoMatch_ReturnsZero()
    {
        // M9 variant: No match returns 0.0
        string genome = "AAAAAAAAAAAAAAAAAAAAAA";
        string probe = "GCGCGCGCGC"; // Won't match in all-A genome
        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        double specificity = ProbeDesigner.CheckSpecificity(probe, genomeIndex);

        Assert.That(specificity, Is.EqualTo(0.0),
            "Non-matching probe should have specificity 0.0");
    }



    // PROBE-VALID-001 — gapped off-target scan + on/off-target separation
    // Evidence: docs/Evidence/PROBE-VALID-001-Evidence.md
    // Sources:
    //   Smith TF, Waterman MS (1981) J Mol Biol 147(1):195-197 (local-alignment recurrence)
    //   Altschul SF et al. (1990) J Mol Biol 215(3):403-410 (gapped local alignment finds indels the ungapped scan misses)
    //   Kane MD et al. (2000) Nucleic Acids Res 28(22):4552-4557 (>75% identity over the probe → off-target)

    // Probe used across the gapped-scan tests.
    private const string GappedProbe = "ACGTACGTACGT"; // 12 nt

    // Reference containing the EXACT on-target at start 5 and an indel off-target
    // ("ACGTACTGTACGT" = probe with a 'T' inserted after position 6) at start 27.
    private static readonly string[] IndelOffTargetReference = new[]
    {
        "NNNNN" + GappedProbe + "NNNNNNNNNN" + "ACGTACTGTACGT" + "NNNNN"
    };

    [Test]
    public void ScanOffTargetsGapped_OnTargetExactMatch_NotCountedAsOffTarget()
    {
        // MG2: The intended on-target (perfect, ungapped, full-coverage exact match) is reported
        // separately and excluded from the off-target count — fixing the OffTargetHits pooling.
        var gapped = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, IndelOffTargetReference);

        var onTarget = gapped.OnTargetHits[0];

        Assert.Multiple(() =>
        {
            Assert.That(onTarget.Start, Is.EqualTo(5),
                "On-target exact match begins at reference position 5");
            Assert.That(onTarget.End, Is.EqualTo(16),
                "On-target exact 12-mer match ends (inclusive) at reference position 16");
            Assert.That(onTarget.Identity, Is.EqualTo(1.0).Within(1e-10),
                "On-target is a 12/12 exact match → identity 1.0");
            Assert.That(onTarget.Coverage, Is.EqualTo(1.0).Within(1e-10),
                "On-target covers the full probe → coverage 1.0");
            Assert.That(onTarget.HasGaps, Is.False,
                "On-target exact match has no gaps");
            Assert.That(gapped.OffTargetHits, Has.None.Matches<ProbeDesigner.GappedProbeHit>(
                h => h.Start == 5),
                "The on-target site must NOT appear among off-target hits");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_IndelOffTarget_HasExactHandDerivedIdentity()
    {
        // MG3: Exact identity/coverage on the hand-derived indel alignment.
        // probe "ACGTAC-GTACGT" vs ref "ACGTACTGTACGT": 12/12 identical aligned columns,
        // one insertion gap → identity 1.0, coverage 1.0, HasGaps true. Site starts at 27.
        var gapped = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, IndelOffTargetReference);

        var off = gapped.OffTargetHits[0];

        Assert.Multiple(() =>
        {
            Assert.That(off.Start, Is.EqualTo(27),
                "Indel off-target begins at reference position 27");
            Assert.That(off.Identity, Is.EqualTo(1.0).Within(1e-10),
                "Indel off-target: 12 identical aligned columns / 12 probe length = 1.0");
            Assert.That(off.Coverage, Is.EqualTo(1.0).Within(1e-10),
                "Indel off-target: 12 ungapped columns / 12 probe length = 1.0");
            Assert.That(off.HasGaps, Is.True,
                "The off-target alignment contains the insertion gap");
            Assert.That(off.AlignedProbe, Is.EqualTo("ACGTAC-GTACGT"),
                "Probe side of the local alignment carries the gap at the insertion point");
            Assert.That(off.AlignedReference, Is.EqualTo("ACGTACTGTACGT"),
                "Reference side of the local alignment is the indel region");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_IndelPlusMismatch_IdentityIsHandDerivedFraction()
    {
        // MG4: An off-target with one insertion AND a trailing mismatch.
        // Region "ACGTACTGTACTT": SW (zero-floor) trims the mismatched "TT" tail →
        // probe "ACGTAC-GTAC" vs ref "ACGTACTGTAC" = 10 identical columns / 12 = 0.8333.
        var references = new[] { "NNNNNNNNNN" + "ACGTACTGTACTT" + "NNNNN" };

        var gapped = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, references, minIdentity: 0.75);

        Assert.Multiple(() =>
        {
            Assert.That(gapped.OnTargetHits, Is.Empty,
                "No perfect exact match in this reference → no on-target");
            Assert.That(gapped.OffTargetCount, Is.EqualTo(1),
                "The indel+mismatch site exceeds the 0.75 threshold → one off-target");
            Assert.That(gapped.OffTargetHits[0].Identity, Is.EqualTo(10.0 / 12.0).Within(1e-10),
                "10 identical aligned columns / 12 probe length = 0.8333... (SW trims the mismatched tail)");
            Assert.That(gapped.OffTargetHits[0].HasGaps, Is.True,
                "The hit required an insertion");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_IdentityThreshold_GatesHits()
    {
        // SG1: The minIdentity threshold (Kane et al. 2000, default 0.75) gates hits.
        // The 0.8333-identity indel+mismatch site is admitted at 0.75 but rejected at 0.90.
        var references = new[] { "NNNNNNNNNN" + "ACGTACTGTACTT" + "NNNNN" };

        var lenient = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, references, minIdentity: 0.75);
        var strict = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, references, minIdentity: 0.90);

        Assert.Multiple(() =>
        {
            Assert.That(lenient.OffTargetCount, Is.EqualTo(1),
                "0.8333 identity >= 0.75 → admitted");
            Assert.That(strict.OffTargetCount, Is.EqualTo(0),
                "0.8333 identity < 0.90 → rejected");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_SpecificProbe_IsSpecificTrue()
    {
        // SG2: A probe with exactly one on-target and no off-targets reports IsSpecific.
        var references = new[] { "GGGGG" + GappedProbe + "GGGGG" };

        var gapped = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, references);

        Assert.Multiple(() =>
        {
            Assert.That(gapped.OnTargetHits, Has.Count.EqualTo(1),
                "Single exact on-target");
            Assert.That(gapped.OffTargetCount, Is.EqualTo(0),
                "No off-targets in flanking poly-G");
            Assert.That(gapped.IsSpecific, Is.True,
                "Exactly one on-target and zero off-targets → specific");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_OnTargetPlusHighIdentityAndIndelOffTargets_SeparatedAndLowIdentityExcluded()
    {
        // SG3: full Kane (2000) separation cross-check on a hand-constructed target set:
        //   ref0 = exact on-target           → 1 on-target (identity 1.0, no gap)
        //   ref1 = 17/20 = 0.85 ungapped     → off-target ABOVE the 0.75 Kane threshold (mismatched, no gap)
        //   ref2 = ~scrambled (<0.75)        → must NOT be called (best local block far below 0.75)
        //   ref3 = one inserted base         → off-target reachable ONLY via a gap (identity 1.0, HasGaps)
        // Hand-derived with an independent Smith-Waterman (BlastDna +2/-3, gap -2):
        //   ref1 alignment ACGTGGCATTACGGCATTCA / ACATGGCATAACGGCAATCA → 17 identical / 20 = 0.85
        //   ref3 alignment ACGTGGCATT-ACGGCATTCA / ACGTGGCATTAACGGCATTCA → 20 identical / 20 = 1.0, one gap
        //   ref2 best local block only 4 identical columns / 20 = 0.20 (< 0.75) → rejected
        const string probe = "ACGTGGCATTACGGCATTCA"; // 20 nt
        var references = new[]
        {
            "TTTTT" + probe + "TTTTT",                          // exact on-target
            "GGGGG" + "ACATGGCATAACGGCAATCA" + "GGGGG",         // 0.85 mismatched off-target
            "CCCCC" + "TTAATTAATTAATTAATTAA" + "CCCCC",         // low-identity, must be excluded
            "AAAAA" + "ACGTGGCATTAACGGCATTCA" + "AAAAA",        // indel-only off-target (1 insertion)
        };

        var gapped = ProbeDesigner.ScanOffTargetsGapped(probe, references, minIdentity: 0.75);

        Assert.Multiple(() =>
        {
            Assert.That(gapped.OnTargetHits, Has.Count.EqualTo(1),
                "Exactly one intended on-target (the perfect exact match in ref0)");
            Assert.That(gapped.OnTargetHits[0].ReferenceIndex, Is.EqualTo(0),
                "The on-target lives in reference 0");

            Assert.That(gapped.OffTargetCount, Is.EqualTo(2),
                "Two genuine off-targets: the 0.85 mismatched site and the indel site");

            var hi = gapped.OffTargetHits.Single(h => h.ReferenceIndex == 1);
            Assert.That(hi.Identity, Is.EqualTo(17.0 / 20.0).Within(1e-10),
                "Mismatched off-target identity = 17/20 = 0.85 (Kane 2000: >0.75 cross-hybridizes)");
            Assert.That(hi.HasGaps, Is.False, "The 0.85 off-target is ungapped (pure substitutions)");

            var indel = gapped.OffTargetHits.Single(h => h.ReferenceIndex == 3);
            Assert.That(indel.Identity, Is.EqualTo(1.0).Within(1e-10),
                "Indel off-target: 20 identical aligned columns / 20 = 1.0");
            Assert.That(indel.HasGaps, Is.True,
                "Indel off-target reachable only via a gap (insertion) — missed by the ungapped Hamming scan");

            Assert.That(gapped.OffTargetHits.Any(h => h.ReferenceIndex == 2), Is.False,
                "Low-identity reference (best block 0.20 < 0.75) must NOT be called an off-target");
            Assert.That(gapped.IsSpecific, Is.False,
                "Two off-targets present → the probe is not specific");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_TwoPerfectCopies_CountsIndependentOfReferenceOrder()
    {
        // A5-3: the scan has no target coordinate, so exactly one perfect, ungapped, full-coverage copy (the first in
        // reference order) is the intended site; every further perfect copy is listed as on-target-class AND counted as
        // an off-target. Identical copies are indistinguishable, so the counts and IsSpecific must not depend on which
        // copy comes first.
        const string probe = "ACGTGGCATTACGGCATTCA"; // 20 nt
        string copyA = "TTTTT" + probe + "TTTTT";
        string copyB = "GGGGG" + probe + "GGGGG";

        var forward = ProbeDesigner.ScanOffTargetsGapped(probe, new[] { copyA, copyB });
        var reversed = ProbeDesigner.ScanOffTargetsGapped(probe, new[] { copyB, copyA });

        Assert.Multiple(() =>
        {
            foreach (var result in new[] { forward, reversed })
            {
                Assert.That(result.OnTargetHits, Has.Count.EqualTo(2), "Both perfect copies are on-target-class");
                Assert.That(result.OffTargetCount, Is.EqualTo(1), "Only the second perfect copy is an off-target");
                Assert.That(result.OnTargetHits[0].ReferenceIndex, Is.EqualTo(0), "First copy = intended site");
                Assert.That(result.OffTargetHits[0].ReferenceIndex, Is.EqualTo(1), "Extra copy = off-target");
                Assert.That(result.OffTargetHits[0].Start, Is.EqualTo(5), "Extra copy site after the 5-nt flank");
                Assert.That(result.IsSpecific, Is.False, "A second binding site makes the probe non-specific");
            }
        });
    }

    [Test]
    public void ScanOffTargetsGapped_NullProbe_ThrowsArgumentNullException()
    {
        // Guard: null probe must throw.
        Assert.Throws<ArgumentNullException>(() =>
            ProbeDesigner.ScanOffTargetsGapped(null!, IndelOffTargetReference),
            "Null probe should throw ArgumentNullException");
    }

    [Test]
    public void ScanOffTargetsGapped_NullReferences_ThrowsArgumentNullException()
    {
        // Guard: null references must throw.
        Assert.Throws<ArgumentNullException>(() =>
            ProbeDesigner.ScanOffTargetsGapped(GappedProbe, null!),
            "Null references should throw ArgumentNullException");
    }

    [Test]
    public void ScanOffTargetsGapped_EmptyProbe_ReturnsNoHits()
    {
        // Guard: empty probe yields no on/off-target hits.
        var gapped = ProbeDesigner.ScanOffTargetsGapped("", IndelOffTargetReference);

        Assert.Multiple(() =>
        {
            Assert.That(gapped.OnTargetHits, Is.Empty, "Empty probe → no on-target hits");
            Assert.That(gapped.OffTargetHits, Is.Empty, "Empty probe → no off-target hits");
        });
    }



    // KA1–KA7 — Karlin–Altschul statistics (Karlin & Altschul 1990, PNAS 87:2264;
    // Altschul et al. 1990, J Mol Biol 215:403).
    // Verbatim formulas (retrieved 2026-06-24):
    //   E = K·m·n·e^{−λS};  S' = (λS − ln K)/ln 2;  E = m·n·2^{−S'};
    //   λ = unique positive root of Σ p_i p_j e^{λ s_ij} = 1.
    //   Sources: NCBI "The Statistics of Sequence Similarity Scores" (Altschul),
    //   https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html ; Durand,
    //   "BLAST (Karlin–Altschul) Statistics", CMU 03-711 (cites both 1990 papers).

    // +1/−3 nucleotide scoring scheme (the scheme NCBI blastn reports λ ≈ 1.37, K ≈ 0.711 for).
    private static readonly Seqeron.Genomics.Infrastructure.ScoringMatrix MatchMismatch1_3 =
        new(Match: 1, Mismatch: -3, GapOpen: -5, GapExtend: -2);

    // KA1 — λ for the +1/−3, uniform-0.25 scheme equals the published NCBI blastn value ≈ 1.374.
    // This pins the numeric root-solver to a sourced value: a wrong solver fails this assertion.
    [Test]
    public void ComputeLambdaNucleotide_Plus1Minus3_UniformFrequencies_MatchesPublishedValue()
    {
        // 0.25·e^{λ·1} + 0.75·e^{λ·(−3)} = 1 → λ ≈ 1.3740631 (NCBI blastn +1/−3).
        double lambda = ProbeDesigner.ComputeLambdaNucleotide(match: 1, mismatch: -3);

        Assert.That(lambda, Is.EqualTo(1.3740631224599755).Within(1e-6),
            "λ for +1/−3 with uniform 0.25 base frequencies must equal the published NCBI blastn value ≈ 1.374");
    }

    // KA2 — the solved λ actually satisfies the defining equation Σ p_i p_j e^{λ s_ij} = 1.
    [Test]
    public void ComputeLambdaNucleotide_SolvedRoot_SatisfiesDefiningEquation()
    {
        double lambda = ProbeDesigner.ComputeLambdaNucleotide(match: 1, mismatch: -3);

        // p(match)=0.25, p(mismatch)=0.75 for four equiprobable bases.
        double f = 0.25 * Math.Exp(lambda * 1) + 0.75 * Math.Exp(lambda * -3);

        Assert.That(f, Is.EqualTo(1.0).Within(1e-9),
            "The returned λ must be a root of Σ p_i p_j e^{λ s_ij} = 1 (Karlin & Altschul 1990)");
    }

    // KA3 — bit score and E-value for a hand-derived (S, m, n) with the +1/−3 scheme.
    // S=30, m=20, n=1000, K=0.711, λ=1.3740631224599755:
    //   S' = (λ·30 − ln 0.711)/ln 2 = 59.9627001142850
    //   E  = 0.711·20·1000·e^{−λ·30} = 1.78015836860839e-14
    [Test]
    public void ComputeKarlinAltschul_HandDerivedExample_MatchesBitScoreAndEValue()
    {
        var stats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 30, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Lambda, Is.EqualTo(1.3740631224599755).Within(1e-6),
                "λ must be the +1/−3 published value");
            Assert.That(stats.BitScore, Is.EqualTo(59.962700114285006).Within(1e-6),
                "S' = (λS − ln K)/ln 2 (Altschul et al. 1990)");
            Assert.That(stats.EValue, Is.EqualTo(1.7801583686083893e-14).Within(1e-24),
                "E = K·m·n·e^{−λS} (Karlin & Altschul 1990)");
        });
    }

    // KA4 — the two equivalent E-value forms agree: E = K·m·n·e^{−λS} = m·n·2^{−S'}.
    [Test]
    public void ComputeKarlinAltschul_EValue_EqualsSearchSpaceTimesTwoToMinusBitScore()
    {
        var stats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 18, queryLength: 25, databaseLength: 5000, scoring: MatchMismatch1_3, k: 0.711);

        double fromBits = (double)stats.QueryLength * stats.DatabaseLength * Math.Pow(2.0, -stats.BitScore);

        Assert.That(stats.EValue, Is.EqualTo(fromBits).Within(stats.EValue * 1e-9),
            "E = m·n·2^{−S'} must equal E = K·m·n·e^{−λS} (Altschul et al. 1990)");
    }

    // KA5 — E-value strictly decreases as the raw score increases (a better hit is less likely by chance).
    [Test]
    public void ComputeKarlinAltschul_EValue_DecreasesAsScoreIncreases()
    {
        var low = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 20, queryLength: 30, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);
        var high = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 21, queryLength: 30, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);

        Assert.That(high.EValue, Is.LessThan(low.EValue),
            "Higher score → lower E (E = K·m·n·e^{−λS} is monotonically decreasing in S)");
    }

    // KA6 — E-value scales linearly with the search space m·n (double n → double E).
    [Test]
    public void ComputeKarlinAltschul_EValue_ScalesLinearlyWithSearchSpace()
    {
        var baseStats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);
        var doubledN = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 20, databaseLength: 2000, scoring: MatchMismatch1_3, k: 0.711);

        Assert.That(doubledN.EValue, Is.EqualTo(2.0 * baseStats.EValue).Within(baseStats.EValue * 1e-9),
            "Doubling n doubles E (E is linear in the search space m·n; Karlin & Altschul 1990)");
    }

    // KA7 (guards) — the Karlin–Altschul preconditions and argument validation.
    [Test]
    public void ComputeLambdaNucleotide_NonPositiveMatch_Throws()
    {
        // No positive score → λ undefined (Altschul et al. 1990).
        Assert.Throws<ArgumentOutOfRangeException>(
            () => ProbeDesigner.ComputeLambdaNucleotide(match: 0, mismatch: -3),
            "A scheme with no positive score has no λ");
    }

    [Test]
    public void ComputeLambdaNucleotide_NonNegativeExpectedScore_Throws()
    {
        // match=3, mismatch=−1 → expected = 0.25·3 + 0.75·(−1) = 0 → not negative → λ undefined.
        Assert.Throws<ArgumentOutOfRangeException>(
            () => ProbeDesigner.ComputeLambdaNucleotide(match: 3, mismatch: -1),
            "Expected per-pair score must be negative for λ to be defined");
    }

    [Test]
    public void ComputeKarlinAltschul_NonPositiveLength_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ProbeDesigner.ComputeKarlinAltschul(10, 0, 1000, MatchMismatch1_3),
                "Query length m must be positive");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ProbeDesigner.ComputeKarlinAltschul(10, 20, 0, MatchMismatch1_3),
                "Database length n must be positive");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ProbeDesigner.ComputeKarlinAltschul(10, 20, 1000, MatchMismatch1_3, k: 0),
                "K must be positive");
        });
    }

    // KA8 — score S = 0 boundary: E collapses to the raw search space scaled by K (E = K·m·n·e^0 = K·m·n),
    // and the bit score is the pure normalization offset S' = −ln K / ln 2 (Karlin & Altschul 1990).
    // m=20, n=1000, K=0.711 → E = 0.711·20·1000 = 14220; S' = −ln(0.711)/ln 2 = 0.4920785350426718.
    [Test]
    public void ComputeKarlinAltschul_ScoreZero_EValueEqualsKTimesSearchSpace()
    {
        var stats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 0, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);

        Assert.Multiple(() =>
        {
            Assert.That(stats.EValue, Is.EqualTo(14220.0).Within(1e-6),
                "At S=0, E = K·m·n·e^0 = K·m·n (Karlin & Altschul 1990)");
            Assert.That(stats.BitScore, Is.EqualTo(0.4920785350426718).Within(1e-9),
                "At S=0, S' = (0 − ln K)/ln 2 = −ln K/ln 2 (Altschul et al. 1990)");
        });
    }

    // KA9 — the K parameter: E is linear in K (E = K·m·n·e^{−λS}), and the bit score shifts by −log2(K),
    // so doubling K doubles E and lowers the bit score by exactly 1 bit (Karlin & Altschul 1990; Altschul et al. 1990).
    [Test]
    public void ComputeKarlinAltschul_DoublingK_DoublesEValueAndLowersBitByOneBit()
    {
        var baseStats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);
        var doubledK = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 1.422);

        Assert.Multiple(() =>
        {
            Assert.That(doubledK.EValue, Is.EqualTo(2.0 * baseStats.EValue).Within(baseStats.EValue * 1e-9),
                "E is linear in K: doubling K doubles E (E = K·m·n·e^{−λS})");
            Assert.That(doubledK.BitScore, Is.EqualTo(baseStats.BitScore - 1.0).Within(1e-9),
                "S' = (λS − ln K)/ln 2: doubling K lowers the bit score by exactly log2(2)=1 bit");
        });
    }

    // KA10 — E increases with the query length m as well (E linear in the search space m·n; KA6 covered n).
    // Doubling m must double E, the symmetric counterpart of doubling n (Karlin & Altschul 1990).
    [Test]
    public void ComputeKarlinAltschul_DoublingQueryLength_DoublesEValue()
    {
        var baseStats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);
        var doubledM = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 40, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);

        Assert.That(doubledM.EValue, Is.EqualTo(2.0 * baseStats.EValue).Within(baseStats.EValue * 1e-9),
            "Doubling m doubles E (E is linear in the search space m·n; Karlin & Altschul 1990)");
    }

    // KA11 — λ root-finder convergence is precise: f(λ) = Σ p_i p_j e^{λ s_ij} − 1 must vanish to
    // near machine precision (the bisection runs 200 iterations on [0,100], far past double resolution).
    // Independently re-solved oracle: λ(1,−3, uniform 0.25) = 1.3740631224599755.
    [Test]
    public void ComputeLambdaNucleotide_RootFinder_ConvergesToMachinePrecision()
    {
        double lambda = ProbeDesigner.ComputeLambdaNucleotide(match: 1, mismatch: -3);
        double residual = 0.25 * Math.Exp(lambda * 1) + 0.75 * Math.Exp(lambda * -3) - 1.0;

        Assert.Multiple(() =>
        {
            Assert.That(lambda, Is.EqualTo(1.3740631224599755).Within(1e-15),
                "Bisection must converge to the exact double-precision root λ = 1.3740631224599755");
            Assert.That(residual, Is.EqualTo(0.0).Within(1e-12),
                "Residual of the defining equation must vanish to near machine precision");
        });
    }

    // KA12 — the DEFAULT scheme path (no scoring arg → BlastDna +2/−3) computes λ from that matrix
    // under uniform 0.25 frequencies: λ(2,−3) = 0.6337314430979077, the ungapped Lambda NCBI blastn 2.12.0+
    // prints for -reward 2 -penalty 3 -ungapped ("0.634"). (blast_stat.c's blastn_values_2_3 {0,0,0.55,0.21}
    // row is the non-affine megablast entry, not the ungapped value.) K is now computed for the same scheme:
    // BlastKarlinLHtoK → 0.4081456625463167 (blastn prints 0.408), no longer the +1/−3 constant 0.711.
    [Test]
    public void ComputeKarlinAltschul_DefaultScheme_UsesBlastDna2_3UngappedLambdaAndK()
    {
        var stats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 30, queryLength: 20, databaseLength: 1000);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Lambda, Is.EqualTo(0.6337314430979077).Within(1e-9),
                "Default scoring is BlastDna (+2/−3); λ is its uniform-0.25 root 0.6337314430979077");
            Assert.That(stats.K, Is.EqualTo(0.4081456625463167).Within(1e-9),
                "K is computed for the +2/−3 scheme (NCBI BlastKarlinLHtoK; blastn -ungapped prints 0.408)");
            Assert.That(stats.EValue, Is.EqualTo(0.4081456625463167 * 20 * 1000 * Math.Exp(-0.6337314430979077 * 30))
                .Within(1e-12), "E = K·m·n·e^{−λS} with the scheme's own K");
        });
    }



    // Oracles: NCBI C++ Toolkit BLAST+ core blast_stat.c / blast_setup.c / blast_hits.c (raw.githubusercontent.com,
    // ncbi/ncbi-cxx-toolkit-public master, retrieved 2026-10-01), ported line-by-line to Python
    // (scratch ka_ref.py: Blast_KarlinLambdaNR, BlastKarlinLtoH, BlastKarlinLHtoK, BLAST_ComputeLengthAdjustment),
    // and NCBI blastn 2.12.0+ (Debian ncbi-blast+) runs whose Lambda/K/H footers, effective search spaces and
    // E-values the Python port reproduces exactly.

    // KA13 — ungapped λ, K, H for match/mismatch schemes under uniform composition (blastn -ungapped footers:
    // 1/−3 → 1.37 0.711 1.31; 2/−3 → 0.634 0.408 0.912; 1/−2 → 1.33 0.621 1.12). Values from the Python port.
    [TestCase(1, -3, 1.3740631224599753, 0.7106027952162398, 1.3072466039090012)]
    [TestCase(2, -3, 0.6337314430979075, 0.4081456625463167, 0.9124383922742278)]
    [TestCase(1, -2, 1.3327057628202603, 0.6209911172603866, 1.1240918464926624)]
    [TestCase(1, -1, 1.0986122886681096, 0.3333333333333333, 0.5493061443340547)]   // closed form low=−1, high=1
    [TestCase(1, -5, 1.3855589708750102, 0.7465088564177568, 1.3794476589446871)]   // closed form high=1
    [TestCase(3, -4, 0.40960018945921306, 0.3400126968378049, 0.8109980590982684)]  // series
    [TestCase(4, -5, 0.30105238556335095, 0.30626406181359495, 0.7531655560686317)] // series
    [TestCase(5, -4, 0.191529283390431, 0.17578794214241147, 0.35672385105693066)]  // series
    [TestCase(3, -2, 0.27117894897553707, 0.13043656940239723, 0.22232354260180037)] // series (low −2, high 3)
    [TestCase(2, -7, 0.6901463137128534, 0.5479023050006323, 1.3431256028187626)]
    public void ComputeUngappedKarlinParameters_MatchesNcbiBlastStat(
        int match, int mismatch, double lambda, double k, double h)
    {
        var p = ProbeDesigner.ComputeUngappedKarlinParameters(match, mismatch);

        Assert.Multiple(() =>
        {
            Assert.That(p.Lambda, Is.EqualTo(lambda).Within(1e-9));
            Assert.That(p.K, Is.EqualTo(k).Within(1e-9));
            Assert.That(p.H, Is.EqualTo(h).Within(1e-9));
            Assert.That(p.Alpha, Is.EqualTo(lambda / h).Within(1e-9), "ungapped α = λ/H (Blast_GetNuclAlphaBeta)");
            Assert.That(p.Gapped, Is.False);
            Assert.That(p.RoundDown, Is.False);
        });
    }

    // KA14 — K is invariant under scaling all scores by an integer, λ scales by its inverse (E is unchanged).
    [Test]
    public void ComputeUngappedKarlinParameters_ScaledScheme_SameKHalfLambda()
    {
        var p23 = ProbeDesigner.ComputeUngappedKarlinParameters(2, -3);
        var p46 = ProbeDesigner.ComputeUngappedKarlinParameters(4, -6);

        Assert.Multiple(() =>
        {
            Assert.That(p46.K, Is.EqualTo(0.4081456625463167).Within(1e-9), "K is a lattice-invariant (Karlin & Altschul 1990)");
            Assert.That(p46.Lambda, Is.EqualTo(p23.Lambda / 2).Within(1e-12));
            Assert.That(p46.H, Is.EqualTo(p23.H).Within(1e-9));
        });
    }

    // KA15 — explicit composition: A,C,G,T = 0.3,0.2,0.2,0.3 → p(match) = Σp² = 0.26 (Python port:
    // λ 1.333431429944318, K 0.6965155054507803, H 1.261161661613095); a uniform vector equals the scalar overload.
    [Test]
    public void ComputeUngappedKarlinParameters_Composition_MatchesNcbiPort()
    {
        var p = ProbeDesigner.ComputeUngappedKarlinParameters(1, -3, new[] { 0.3, 0.2, 0.2, 0.3 });
        var uniform = ProbeDesigner.ComputeUngappedKarlinParameters(1, -3, new[] { 1.0, 1.0, 1.0, 1.0 });

        Assert.Multiple(() =>
        {
            Assert.That(p.Lambda, Is.EqualTo(1.333431429944318).Within(1e-9));
            Assert.That(p.K, Is.EqualTo(0.6965155054507803).Within(1e-9));
            Assert.That(p.H, Is.EqualTo(1.261161661613095).Within(1e-9));
            Assert.That(uniform, Is.EqualTo(ProbeDesigner.ComputeUngappedKarlinParameters(1, -3)));
        });
    }

    [Test]
    public void ComputeUngappedKarlinParameters_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.ComputeUngappedKarlinParameters(1, -3, 0.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.ComputeUngappedKarlinParameters(1, -3, 0.0));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.ComputeUngappedKarlinParameters(1, -3, double.NaN));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.ComputeLambdaNucleotide(1, -3, 0.6));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.ComputeUngappedKarlinParameters(3, -1));
            Assert.Throws<ArgumentException>(() => ProbeDesigner.ComputeUngappedKarlinParameters(1, -3, new[] { 0.5, 0.5, 0.0 }));
            Assert.Throws<ArgumentException>(() => ProbeDesigner.ComputeUngappedKarlinParameters(1, -3, new[] { 0.5, 0.5, -0.1, 0.1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.ComputeUngappedKarlinParameters(1, -3, new[] { 1.0, 0, 0, 0 }));
        });
    }

    // KA16 — gapped parameters from blast_stat.c blastn_values_* (blastn "Gapped" footers: 2/−3 5/2 → 0.625 0.410 0.780;
    // 1/−3 2/2 → 1.37 0.700 1.20; 1/−2 2/2 → 1.33 0.620 1.10; 2/−3 4/4 → 0.630 0.420 0.840; 1/−3 5/5 → ungapped).
    [TestCase(2, -3, 5, 2, 0.625, 0.41, 0.78, 0.8, -2.0, true)]
    [TestCase(1, -3, 2, 2, 1.37, 0.70, 1.2, 1.1, 0.0, false)]
    [TestCase(1, -2, 2, 2, 1.33, 0.62, 1.1, 1.2, 0.0, false)]
    [TestCase(2, -3, 4, 4, 0.63, 0.42, 0.84, 0.75, -2.0, true)]
    [TestCase(2, -3, 0, 0, 0.55, 0.21, 0.46, 1.2, -5.0, true)]   // non-affine (megablast) row
    [TestCase(4, -6, 10, 4, 0.3125, 0.41, 0.78, 0.4, -2.0, true)] // gcd 2: gaps ×2, λ and α ÷2
    [TestCase(3, -4, 6, 3, 0.389, 0.25, 0.56, 0.7, -5.0, true)]
    public void GetBlastnGappedKarlinParameters_MatchesBlastStatTables(
        int reward, int penalty, int open, int extend, double lambda, double k, double h, double alpha, double beta, bool roundDown)
    {
        var p = ProbeDesigner.GetBlastnGappedKarlinParameters(reward, penalty, open, extend);

        Assert.That(p, Is.EqualTo(new ProbeDesigner.KarlinAltschulParameters(lambda, k, h, alpha, beta, roundDown, true)));
    }

    [Test]
    public void GetBlastnGappedKarlinParameters_InfiniteGapDomain_UsesUngappedValues()
    {
        var p = ProbeDesigner.GetBlastnGappedKarlinParameters(1, -3, 5, 5);
        var u = ProbeDesigner.ComputeUngappedKarlinParameters(1, -3);

        Assert.That(p, Is.EqualTo(u with { Gapped = true }),
            "gap costs ≥ gap_open_max/gap_extend_max (2/2) copy the ungapped Karlin block (blastn prints 1.37 0.711 1.31)");
    }

    [Test]
    public void GetBlastnGappedKarlinParameters_Unsupported_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => ProbeDesigner.GetBlastnGappedKarlinParameters(2, -3, 1, 1),
                "1/1 is not a 2/−3 table entry and is below the infinite domain (6/4)");
            Assert.Throws<ArgumentException>(() => ProbeDesigner.GetBlastnGappedKarlinParameters(7, -11, 5, 2),
                "7/−11 has no blastn_values table");
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.GetBlastnGappedKarlinParameters(2, -3, -1, 2));
        });
    }

    // KA17 — length adjustment / effective search space / E-value = blastn 2.12.0+ (query 40 nt, subject 3079 nt):
    // 2/−3 gap 5/2 → "Effective search space used: 88972", score 80 → 73.4 bits, E 7e-18; score 15 → E 5.8 (as 14).
    [Test]
    public void ComputeBlastnStatistics_Bl2seq_MatchesBlastn()
    {
        var hit = ProbeDesigner.ComputeBlastnStatistics(80, 40, 3079);
        var odd = ProbeDesigner.ComputeBlastnStatistics(15, 40, 3079);

        Assert.Multiple(() =>
        {
            Assert.That(hit.LengthAdjustment, Is.EqualTo(11));
            Assert.That(hit.EffectiveSearchSpace, Is.EqualTo(88972));
            Assert.That(hit.EValue, Is.EqualTo(7.035793990394873e-18).Within(1e-27));
            Assert.That(hit.BitScore, Is.EqualTo(73.42105622960482).Within(1e-9));
            Assert.That(odd.EValueScore, Is.EqualTo(14), "2/−3 table: odd scores are rounded down to even");
            Assert.That(odd.EValue, Is.EqualTo(5.780434617461434).Within(1e-9));
        });
    }

    // KA18 — database of N = 5 sequences, 6040 letters: blastn "Effective search space used: 167440", E 1.32e-17.
    [Test]
    public void ComputeBlastnStatistics_MultiSequenceDatabase_MatchesBlastn()
    {
        var s = ProbeDesigner.ComputeBlastnStatistics(80, 40, 6040, databaseSequenceCount: 5);

        Assert.Multiple(() =>
        {
            Assert.That(s.LengthAdjustment, Is.EqualTo(12));
            Assert.That(s.EffectiveSearchSpace, Is.EqualTo(167440));
            Assert.That(s.EValue, Is.EqualTo(1.3240944856266213e-17).Within(1e-26));
        });
    }

    // KA19 — other schemes: 1/−3 gap 2/2 (space 98272, score 31 → E 2e-14) and ungapped 2/−3 (space 95170,
    // score 80 → 74.4 bits, E 4e-18; ungapped β = −2 from s_GetUngappedBeta).
    [Test]
    public void ComputeBlastnStatistics_OtherSchemes_MatchBlastn()
    {
        var m13 = new Seqeron.Genomics.Infrastructure.ScoringMatrix(Match: 1, Mismatch: -3, GapOpen: -2, GapExtend: -2);
        var g = ProbeDesigner.ComputeBlastnStatistics(31, 40, 3079, scoring: m13);
        var u = ProbeDesigner.ComputeBlastnStatistics(80, 40, 3079, gapped: false);

        Assert.Multiple(() =>
        {
            Assert.That(g.EffectiveSearchSpace, Is.EqualTo(98272));
            Assert.That(g.EValue, Is.EqualTo(2.4719585736391905e-14).Within(1e-23));
            Assert.That(u.EffectiveSearchSpace, Is.EqualTo(95170));
            Assert.That(u.EValueScore, Is.EqualTo(80));
            Assert.That(u.EValue, Is.EqualTo(3.725887650102598e-18).Within(1e-24));
            Assert.That(u.BitScore, Is.EqualTo(74.4353407865069).Within(1e-6));
        });
    }

    // KA20 — end-to-end through the canonical affine aligner: probe vs a 279-nt subject carrying the probe with one
    // mismatch and a 1-nt deletion. blastn -task blastn (2/−3, 5/2) bl2seq: score 66, 60.8 bits, E 4.33e-15,
    // effective search space 8672; Biopython PairwiseAligner local (2, −3, open −7, extend −2) score 66.0.
    private const string EValueProbe = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACT";
    private const string EValueSubject =
        "CTTGTCTCCAAGTACCCATTTAGTAGACAAATCGTTCCATCACCAATTCGCTGGTTGTTGAACTATACGACCGGGGCACACTGCACTCAGTTCCCATTTAGAGGATCC" +
        "TAGCCTAGCTACGCTAAAGACAATGACATAACATACAGTCAGCACGAAACTGCGTTTGCGCATCAGGCTGTCCCATACATCAAGCGGTTCCCCTCAAATTATCCGGAC" +
        "TCGGTAAGGGCAGCGAGTAAATATTTTACAATACGTTTCTTGTCAATCTGCTGCTTTGTACGC";

    [Test]
    public void ComputeLengthAdjustment_TinySearchSpace_IsZero()
    {
        // c = n·m − max(m, n)/K < 0 → BLAST returns 0.
        Assert.Multiple(() =>
        {
            Assert.That(ProbeDesigner.ComputeLengthAdjustment(0.41, 1.28, -2, 1, 1), Is.EqualTo(0));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.ComputeLengthAdjustment(0, 1, 0, 10, 10));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.ComputeBlastnStatistics(10, 20, 100, databaseSequenceCount: 0));
        });
    }

    // KA21 — ComputeKarlinAltschul without k computes K for the scheme (+1/−3 → 0.7106027952162398 ≈ blastn 0.711).
    [Test]
    public void ComputeKarlinAltschul_NoK_ComputesSchemeK()
    {
        var stats = ProbeDesigner.ComputeKarlinAltschul(10, 100, 1000, MatchMismatch1_3);

        Assert.That(stats.K, Is.EqualTo(0.7106027952162398).Within(1e-9));
    }

    // KA22 — gcd > 1 schemes, K as BLAST+ computes it (KarlinKMethod.NcbiBlast): BlastKarlinLHtoK reduces low/high/λ by
    // δ = gcd but its series reads probArrayStartLow[j] from the UNREDUCED array. blastn 2.12.0+ -task blastn -ungapped
    // footers (λ K H): 4/−6 0.317 1.17 0.912; 6/−9 0.211 1.17 0.912; 8/−12 0.158 1.17 0.912; 4/−10 0.340 1.06 1.24;
    // 6/−4 0.136 1.63 0.222; 6/−10 0.216 1.03 0.997; 8/−10 0.151 1.07 0.753; 10/−8 0.0958 1.31 0.357; closed forms
    // 4/−2 0.132 0.0532 0.0722, 2/−2 0.549 0.333 0.549, 2/−4 0.666 0.621 1.12. Exact values: line-by-line Python port of
    // blast_stat.c (Blast_KarlinBlkUngappedCalc / BlastKarlinLtoH / BlastKarlinLHtoK, verbatim indexing).
    [TestCase(4, -6, 0.31686572154895387, 1.1666856431064105, 0.9124383922742287)]
    [TestCase(6, -9, 0.21124381436596928, 1.1666856431064097, 0.9124383922742297)]
    [TestCase(8, -12, 0.15843286077447694, 1.1666856431064099, 0.9124383922742294)]
    [TestCase(4, -10, 0.34025265263757676, 1.0551429674627688, 1.2420803627364985)]
    [TestCase(6, -4, 0.1355894744877687, 1.633485999782881, 0.22232354260180115)]
    [TestCase(6, -10, 0.21596619831296326, 1.030863998562112, 0.9968202121665337)]
    [TestCase(8, -10, 0.15052619278167545, 1.0682924577828306, 0.7531655560686317)]
    [TestCase(10, -8, 0.09576464169521554, 1.3071030460951232, 0.356723851056931)]
    [TestCase(4, -2, 0.1322485471578545, 0.05322292075469216, 0.07218608985058102)] // closed form (reduced low −1)
    [TestCase(2, -2, 0.549306144334055, 0.3333333333333333, 0.5493061443340556)]    // closed form (reduced ±1)
    [TestCase(2, -4, 0.6663528814101303, 0.6209911172603868, 1.1240918464926628)]   // closed form (reduced high +1)
    [TestCase(2, -3, 0.6337314430979077, 0.4081456625463164, 0.9124383922742287)]   // δ = 1: same as ReducedLattice
    public void ComputeUngappedKarlinParameters_NcbiBlastK_MatchesBlastn(int match, int mismatch, double lambda, double k, double h)
    {
        var ncbi = ProbeDesigner.ComputeUngappedKarlinParameters(
            match, mismatch, 0.25, ProbeDesigner.KarlinKMethod.NcbiBlast);

        Assert.Multiple(() =>
        {
            Assert.That(ncbi.Lambda, Is.EqualTo(lambda).Within(1e-9));
            Assert.That(ncbi.K, Is.EqualTo(k).Within(1e-9));
            Assert.That(ncbi.H, Is.EqualTo(h).Within(1e-9));
        });
    }

    // KA23 — the default stays the scale-invariant reduced lattice (Karlin & Altschul 1990): K(δ·a, −δ·b) = K(a, −b);
    // for closed-form and δ = 1 schemes both methods agree.
    [TestCase(4, -6, 2, -3)]
    [TestCase(4, -10, 2, -5)]
    [TestCase(6, -4, 3, -2)]
    [TestCase(10, -8, 5, -4)]
    [TestCase(4, -2, 2, -1)]
    public void ComputeUngappedKarlinParameters_ReducedLattice_IsScaleInvariant(int match, int mismatch, int rMatch, int rMismatch)
    {
        var scaled = ProbeDesigner.ComputeUngappedKarlinParameters(match, mismatch);
        var reduced = ProbeDesigner.ComputeUngappedKarlinParameters(rMatch, rMismatch);
        var explicitReduced = ProbeDesigner.ComputeUngappedKarlinParameters(
            match, mismatch, 0.25, ProbeDesigner.KarlinKMethod.ReducedLattice);

        Assert.Multiple(() =>
        {
            Assert.That(scaled.K, Is.EqualTo(reduced.K).Within(1e-12));
            Assert.That(explicitReduced, Is.EqualTo(scaled));
        });
    }

    [Test]
    public void ComputeUngappedKarlinParameters_KMethods_AgreeForCoprimeAndClosedForm_DifferOtherwise()
    {
        Assert.Multiple(() =>
        {
            foreach (var (m, mm) in new[] { (1, -3), (2, -3), (3, -4), (5, -4), (2, -2), (2, -4), (4, -2) })
                Assert.That(ProbeDesigner.ComputeUngappedKarlinParameters(m, mm, 0.25, ProbeDesigner.KarlinKMethod.NcbiBlast).K,
                    Is.EqualTo(ProbeDesigner.ComputeUngappedKarlinParameters(m, mm).K).Within(1e-12), $"{m}/{mm}");
            Assert.That(ProbeDesigner.ComputeUngappedKarlinParameters(4, -6).K, Is.EqualTo(0.4081456625463167).Within(1e-9));
            Assert.That(ProbeDesigner.ComputeUngappedKarlinParameters(1, -3, new[] { 0.3, 0.2, 0.2, 0.3 },
                ProbeDesigner.KarlinKMethod.NcbiBlast).K, Is.EqualTo(0.6965155054507803).Within(1e-9), "δ = 1 composition (KA15)");
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                ProbeDesigner.ComputeUngappedKarlinParameters(4, -6, 0.25, (ProbeDesigner.KarlinKMethod)7));
        });
    }

    // KA24 — ComputeBlastnStatistics (BLAST semantics) uses BLAST+'s K by default. blastn 2.12.0+ -task blastn, 200-nt
    // query vs 3100-nt subject: 4/−6 -ungapped → "Effective search space used: 573996", score 400 → 182 bits, E 6.03e-50;
    // 6/−4 -ungapped → 425600, score 602 → 117 bits, E 2.47e-30; 4/−6 gap 12/8 (≥ 2/−3 infinite domain 6/4 × 2) →
    // "Gapped 0.317 1.17 0.912", 573996, 6.03e-50. Exact E/bits from the Python port (E = space·K·e^{−λS}).
    [Test]
    public void ComputeBlastnStatistics_GcdScheme_UsesBlastK_MatchesBlastn()
    {
        var m46 = new Seqeron.Genomics.Infrastructure.ScoringMatrix(Match: 4, Mismatch: -6, GapOpen: -12, GapExtend: -8);
        var m64 = new Seqeron.Genomics.Infrastructure.ScoringMatrix(Match: 6, Mismatch: -4, GapOpen: -10, GapExtend: -6);
        var u46 = ProbeDesigner.ComputeBlastnStatistics(400, 200, 3100, scoring: m46, gapped: false);
        var g46 = ProbeDesigner.ComputeBlastnStatistics(400, 200, 3100, scoring: m46);
        var u64 = ProbeDesigner.ComputeBlastnStatistics(602, 200, 3100, scoring: m64, gapped: false);
        var reduced = ProbeDesigner.ComputeBlastnStatistics(400, 200, 3100, scoring: m46, gapped: false,
            kMethod: ProbeDesigner.KarlinKMethod.ReducedLattice);

        Assert.Multiple(() =>
        {
            Assert.That(u46.Parameters.K, Is.EqualTo(1.1666856431064105).Within(1e-9));
            Assert.That(u46.LengthAdjustment, Is.EqualTo(14));
            Assert.That(u46.EffectiveSearchSpace, Is.EqualTo(573996));
            Assert.That(u46.EValue, Is.EqualTo(6.034606696310044e-50).Within(1e-58));
            Assert.That(u46.BitScore, Is.EqualTo(182.63382615522124).Within(1e-6));
            Assert.That(g46.Parameters, Is.EqualTo(u46.Parameters with { Gapped = true, RoundDown = true }),
                "infinite domain copies BLAST's ungapped block (2/−3 table: round_down)");
            Assert.That(g46.EffectiveSearchSpace, Is.EqualTo(573996));
            Assert.That(g46.EValue, Is.EqualTo(6.034606696310044e-50).Within(1e-58));
            Assert.That(u64.Parameters.K, Is.EqualTo(1.633485999782881).Within(1e-9));
            Assert.That(u64.LengthAdjustment, Is.EqualTo(60));
            Assert.That(u64.EffectiveSearchSpace, Is.EqualTo(425600));
            Assert.That(u64.EValue, Is.EqualTo(2.471093453808567e-30).Within(1e-38));
            Assert.That(u64.BitScore, Is.EqualTo(117.05183189919178).Within(1e-6));
            Assert.That(reduced.Parameters.K, Is.EqualTo(0.4081456625463167).Within(1e-9), "opt-out: scale-invariant K");
            Assert.That(ProbeDesigner.GetBlastnGappedKarlinParameters(4, -6, 12, 8, ProbeDesigner.KarlinKMethod.ReducedLattice).K,
                Is.EqualTo(0.4081456625463167).Within(1e-9));
            Assert.That(ProbeDesigner.ComputeKarlinAltschul(400, 200, 3100, m46, kMethod: ProbeDesigner.KarlinKMethod.NcbiBlast).K,
                Is.EqualTo(1.1666856431064105).Within(1e-9));
            Assert.That(ProbeDesigner.ComputeKarlinAltschul(400, 200, 3100, m46).K,
                Is.EqualTo(0.4081456625463167).Within(1e-9), "ComputeKarlinAltschul default unchanged (reduced lattice)");
        });
    }

    // KA25 — ungapped schemes with no BLAST+ blastn_values_* table: Blast_GetNuclAlphaBeta returns the
    // s_GetNuclValuesArray error without setting α/β, BLAST_CalcEffLengths keeps α = β = 0 → ℓ = 0, search space m·n.
    // NCBI blastn 2.12.0+ -task blastn -ungapped -word_size 7 -dust no, 73-nt query vs 446-nt subject (query with 3
    // substitutions embedded, seed 71): "Effective search space used: 32558" (= 73·446) for every unsupported scheme;
    // top-HSP raw score, E-value and bit score as printed (-outfmt "6 score evalue bitscore"). Auditor repro (same
    // lengths): 2/−1 raw 98 → E 9.58e-09.
    [TestCase(3, -5, 195, "3.40e-33", 122.0)]
    [TestCase(1, -6, 52, "1.21e-27", 104.0)]
    [TestCase(3, -7, 189, "1.34e-33", 124.0)]
    [TestCase(5, -7, 329, "2.68e-32", 119.0)]
    [TestCase(2, -1, 137, "3.17e-13", 56.5)]
    [TestCase(6, -10, 390, "8.84e-33", 121.0)]
    [TestCase(2, -1, 98, "9.58e-09", double.NaN)]
    public void ComputeBlastnStatistics_UngappedSchemeWithoutBlastTable_NoLengthAdjustment_MatchesBlastn(
        int reward, int penalty, int rawScore, string blastnEValue, double blastnBits)
    {
        var m = new Seqeron.Genomics.Infrastructure.ScoringMatrix(Match: reward, Mismatch: penalty, GapOpen: -5, GapExtend: -2);
        var st = ProbeDesigner.ComputeBlastnStatistics(rawScore, 73, 446, scoring: m, gapped: false);
        var ka = ProbeDesigner.ComputeUngappedKarlinParameters(reward, penalty, kMethod: ProbeDesigner.KarlinKMethod.NcbiBlast);

        Assert.Multiple(() =>
        {
            Assert.That(st.LengthAdjustment, Is.EqualTo(0));
            Assert.That(st.EffectiveSearchSpace, Is.EqualTo(32558));
            Assert.That(st.Parameters.Alpha, Is.EqualTo(0.0));
            Assert.That(st.Parameters.Beta, Is.EqualTo(0.0));
            Assert.That(st.Parameters, Is.EqualTo(ka with { Alpha = 0, Beta = 0 }), "λ, K, H unchanged");
            Assert.That(st.EValue, Is.EqualTo(32558 * ka.K * Math.Exp(-ka.Lambda * rawScore)).Within(1e-12).Percent);
            Assert.That(st.EValue.ToString("0.00e+00", System.Globalization.CultureInfo.InvariantCulture),
                Is.EqualTo(blastnEValue));
            // Tabular output: bit scores > 99.9 are printed truncated to an integer, smaller ones with one decimal.
            if (!double.IsNaN(blastnBits))
                Assert.That(blastnBits >= 100 ? Math.Floor(st.BitScore) : Math.Round(st.BitScore, 1), Is.EqualTo(blastnBits));
        });
    }

    // KA25 control — tabulated schemes keep BLAST's α = λ/H, β = s_GetUngappedBeta and ℓ > 0 (same blastn runs:
    // 2/−3 → 28470, 1/−3 → 28974, 1/−2 → 28470, 4/−6 (gcd 2 → 2/−3 table, β 0) → 26970).
    [TestCase(2, -3, 131, 28470, "1.02e-32")]
    [TestCase(1, -3, 61, 28974, "8.17e-33")]
    [TestCase(1, -2, 64, 28470, "1.60e-33")]
    [TestCase(4, -6, 262, 26970, "2.77e-32")]
    public void ComputeBlastnStatistics_UngappedTabulatedScheme_KeepsLengthAdjustment_MatchesBlastn(
        int reward, int penalty, int rawScore, double blastnSpace, string blastnEValue)
    {
        var m = new Seqeron.Genomics.Infrastructure.ScoringMatrix(Match: reward, Mismatch: penalty, GapOpen: -5, GapExtend: -2);
        var st = ProbeDesigner.ComputeBlastnStatistics(rawScore, 73, 446, scoring: m, gapped: false);

        Assert.Multiple(() =>
        {
            Assert.That(st.LengthAdjustment, Is.GreaterThan(0));
            Assert.That(st.EffectiveSearchSpace, Is.EqualTo(blastnSpace));
            Assert.That(st.Parameters.Alpha, Is.EqualTo(st.Parameters.Lambda / st.Parameters.H));
            Assert.That(st.EValue.ToString("0.00e+00", System.Globalization.CultureInfo.InvariantCulture),
                Is.EqualTo(blastnEValue));
        });
    }



    // Kane et al. (2000) NAR 28:4552: non-targets > 75 % similar over the probe, or sharing a stretch of
    // > 15 contiguous identical bases, may cross-hybridize. Reference values: Biopython 1.88
    // PairwiseAligner(mode='local', match 2, mismatch −3, open_gap_score −7, extend_gap_score −2) = BLAST+ blastn
    // scoring (existence 5, extension 2) — score and the identical-column count of the (unique) optimal alignment;
    // longest common substring by dynamic programming.
    private const string KaneProbe = "TATGCCTCCGGTACATCAACTACAGTTAGCCTTAAGAGAAAAATCCCAAA";
    // Probe with a substitution at every 5th position (40/50 identical), random flanks.
    private const string KaneNonTargetA = "CCGCACCATGAGACTGTTTCTATGGCTCCTGTACCTCAAGTACATTTAGGCTTACGAGACAAATGCCAACCACATCGGCTTCGCACGTCT";
    // Probe[10..28) (18 nt) embedded in random sequence.
    private const string KaneNonTargetB = "GGTCCCACTGATAACGTGTTACCGGCTCTAGTACATCAACTACAGTTACCTAATGCAAAAAACTGTTAACACTTTAAA";
    // Unrelated random sequence.
    private const string KaneNonTargetC = "AATGTGATAGGATGTTAAAAGCGCCGAGACGGCGGTCTGCGATGTACCCCGCAACTGGTTCTTCCCCAGCCGCGGGGGTA";
    // Reverse complement of the probe embedded in random sequence.
    private const string KaneNonTargetD = "CCCCCGGCATTGTTCTTTGGGATTTTTCTCTTAAGGCTAACTGTAGTTGATGTACCGGAGGCATACGGCGCGGAATGACG";
    // Probe[0..15) (exactly 15 nt) embedded in random sequence — at, not above, the contiguous threshold.
    private const string KaneNonTargetE = "AAATTCCCGAAGAAGCGACTTGGTAGGGGATATGCCTCCGGTACACTAGTTCGTTATCCGTCTCGTATTCTGCTC";



    [Test]
    public void CheckSpecificity_BothStrands_CountsReverseComplementSites()
    {
        // Probe once on the indexed strand and once as its reverse complement (a site on the other strand).
        string genome = "TTTT" + UniqueProbe + "TTTT" + DnaSequence.GetReverseComplementString("ACCGTTAGGCATCGATGCAA") + "TTTT";
        var tree = global::SuffixTree.SuffixTree.Build(genome);
        const string probe = "ACCGTTAGGCATCGATGCAA";

        Assert.Multiple(() =>
        {
            Assert.That(ProbeDesigner.CheckSpecificity(probe, tree), Is.EqualTo(0.0), "indexed strand only: no site");
            Assert.That(ProbeDesigner.CheckSpecificity(probe, tree, bothStrands: true), Is.EqualTo(1.0));
            // Reverse-palindromic probe (GAATTC-like): its reverse complement is itself → counted once.
            var pal = global::SuffixTree.SuffixTree.Build("TTTGAATTCTTT");
            Assert.That(ProbeDesigner.CheckSpecificity("GAATTC", pal, bothStrands: true), Is.EqualTo(1.0));
            // Probe and its reverse complement each present once → two binding sites.
            var two = global::SuffixTree.SuffixTree.Build("AAAA" + probe + "AAAA" + DnaSequence.GetReverseComplementString(probe) + "AAAA");
            Assert.That(ProbeDesigner.CheckSpecificity(probe, two, bothStrands: true), Is.EqualTo(0.5));
            Assert.That(ProbeDesigner.CheckSpecificity(probe, two), Is.EqualTo(1.0));
        });
    }



    // Fixtures (random.seed(20261008)); expected values from an independent Python brute-force Hamming scan
    // (identity (L - d) / L, longest run of identical positions; Kane et al. 2000: > 0.75 or > 15 nt).
    private const string KaneSiteProbe12 = "GATCCGACGCTA";
    private const string KaneSiteProbe40 = "TATGCCGTACAGTTTTAAGATAGAGCGAAAGCGCAGACAA";

    // Exact site at 10; probe with positions 0..11 substituted at 60 (identity 0.70, run 28); probe with positions
    // 1, 4, …, 34 substituted at 110 (identity 0.70, longest run 5).
    private const string KaneSiteReference40 =
        "AAAAAAAAAATATGCCGTACAGTTTTAAGATAGAGCGAAAGCGCAGACAAAAAAAAAAAAACATGGTACGCTTTTTAAGATAGAGCGAAAGCGCAGACAA"
        + "AAAAAAAAAATCTGGCGAACCGTATTCAGCTATAGGGACAGGGCCGACAAAAAAAAAAAA";



    // Expected values: independent Python brute force over both strands (oracle.py — every reference window compared
    // position by position with the probe and with its Biopython reverse complement; the union of the hit positions per
    // reference is the site set, a position cross-hybridizing when either orientation meets a Kane criterion:
    // identity (L - d) / L > 0.75 or a run of identical positions > 15). Fixture random.seed(20261008):
    // reference = A + P + B + rc(P with positions 3, 11 substituted) + C + rc(P with 1, 5, 9, 13, 17 substituted) + D.
    private const string BothStrandsProbe = "GATCCGACGCTATATGCCGT";
    private const string BothStrandsReference =
        "ACAGTTTTAAGATAGGATCCGACGCTATATGCCGTAGCGAAAGCGCAGACACGGCATAGAGCGTCGCATCAATAAATAATCCGTAACCGCAGATACCGTAGGAGCGGGAGACCTGGCACA";

    // Random cases (Python random.Random(7): probe 8–30 nt, 1–3 references of 0–120 random nt with 0–3 planted copies of
    // the probe or its reverse complement carrying 0–4 substitutions; radius 0–4) where the two modes differ. Columns:
    // probe, references, radius, (hits, Kane sites) given strand, (hits, Kane sites) both strands — from oracle.py.
    private static IEnumerable<TestCaseData> BothStrandsOracleCases()
    {
        yield return new TestCaseData("CCCCCCAATGCCCCGC", new[] { "TAGGGCGGTGCGCGCCCAATGCGCGGGGCATTGGGGGGGCCGGATTTGGTGGGTAA" }, 4, 0, 0, 1, 1);
        yield return new TestCaseData("ACTAAAGCAAGCTCCCTTGGACTA", new[] { "TTCCGTTCCCTAGCAGTCGGCGCTAGTCCAAGTGCGCTTGCTTTAGTTAACGAGAAG" }, 4, 0, 0, 1, 1);
        yield return new TestCaseData("AGCGGAGACGGTAG", new[] { "GAACGGCTATAATCTACCGTCTCCGCTAAGCCGTCGGTAAGCTTAAACTTCTTCAGGCG" }, 1, 0, 0, 1, 1);
        yield return new TestCaseData("AGGTTCTAAAGGCTATGC", new[] { "GTGAGTAACATTGCATAGCCTTTAGAACCTCGCGCCACAGGTTCAAAAGGCAATGCATGAGCACG" }, 1, 0, 0, 1, 1);
        yield return new TestCaseData("CGATGCAAATTCCTCTGTTTCTAGT", new[] { "ACTAGAAACAGAGGAATTTGCATCGACATGACTAGAAACCGAGGAATTTGCATCGTCTACTTTAACTATTCGT" }, 1, 0, 0, 2, 2);
        yield return new TestCaseData("CTGCCCACCAGTCGCGAGGCAA", new[] { "TCCACTAACAGTACAGGCACGATCTCTATTCATTCACCAACAGCAGTCCCGAAGCCTTGCCTCGCGACTGGTGGGCAG" }, 0, 0, 0, 1, 1);
        yield return new TestCaseData("ACGGCGCTTTTATTTCGGGGTC", new[] { "AATACCCCGAAATCAAAGCGCCGTGGTCGTCCAAGGAGTGCAGCTATATTCATTTGCTTCAAAAAGTAGTCATTCCGGT" }, 3, 0, 0, 1, 1);
        yield return new TestCaseData("TATGAGAAAAGTTG", new[] { "CTTATTAATCCAACTTTTCTCATATCATGTAGCCGGCCCGCAGAAGCAGCCGGTTTTTGTTAGATATTAGAAAAGTAGCG" }, 2, 1, 1, 2, 2);
        yield return new TestCaseData("CTTGGCAT", new[] { "TTCGAGGTTTATTTCTTGGCATGTGAGCAGCATCGATAAGTATGCCAAGTATGGGCCTTGGCATAAAATTACGGGGGTAACGCCACCAGTTC" }, 4, 11, 2, 20, 3);
        yield return new TestCaseData("TGTTGAGC", new[] { "CGCTATGTCTAAACGCCGCGCTTAAGGCACAAGAGTTTCTTTTGAGGAGAAGTTCTATGAGTTTGTCGAGCACGGCACTCGCAAGAGAGACTCG" }, 2, 2, 1, 3, 1);
        yield return new TestCaseData("TCGCGATATAA", new[] { "TTATATTTATATCGGGACGGGAGCTGGCGATAAACTTTC", "TAATTAGTTATATCGGGAATACTCGCGATATAAGCGCGCGCTCCATTTACAGCCCAACGCTACAAGG" }, 4, 2, 1, 6, 3);
    }



    // 74-nt probe with a 10-bp perfect stem GCGCATCCAG·CTGGATGCGC closing the 4-nt loop TTTT;
    // the A/C-only flanks cannot pair. Auditor's repro: the former private scan (loop of exactly 3 nt) missed it.
    private const string StemLoop74 =
        "ACCAACCACACAACCCAAACCCAACGCGCATCCAGTTTTCTGGATGCGCAAAAACCCAACCCACCCCACCACCC";

    // 73-nt: a 5-bp "stem" ACCAC·GTAGT (4 of 5 complementary, central C·A mismatch) closing a 3-nt loop in an A/C-only
    // background, no exact 4-bp stem anywhere: flagged only by the former unsourced "≥ 80 % matched" tolerance.
    private const string MismatchedStem73 =
        "ACCAACCACACAACCCAAACCCAACACCAAACCACAAAGTAGTACCAACCACACAACCCAAACCCAACACCAA";

    // F55's random 61-mer (PrimerDesigner_NtthalMaxAlign_Tests.Random61): thal.c -DTHAL_MAX_ALIGN=10000 at
    // 50 mM / 0 / 0 / 50 nM: hairpin 33.980529935122263 °C, ANY −0.968 °C, END1 −79.6 °C (all ≤ 47 °C: no sourced
    // structure). It contains exact 4-bp stems closing ≥ 3-nt loops (independent Python scan: True), as 98 % of random
    // 61-mers do (F61 coverage).
    private const string Random61 = "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGT";


    // Audit round 7 (F68 follow-up): mismatch = int.MinValue overflowed the private gcd / series span. BLAST+
    // (blast_stat.h) limits one-letter scores to BLAST_SCORE_MIN = INT2_MIN … BLAST_SCORE_MAX = INT2_MAX.
    [TestCase(1, int.MinValue, "mismatch")]
    [TestCase(1, short.MinValue - 1, "mismatch")]
    [TestCase(1, short.MinValue, "mismatch")]   // blastn 2.12: 1/−32768 "Could not calculate ungapped Karlin-Altschul parameters"
    [TestCase(2, short.MinValue, "mismatch")]   // blastn 2.12: 2/−32768 likewise
    [TestCase(int.MaxValue, -3, "match")]
    [TestCase(short.MaxValue + 1, -3, "match")]
    [TestCase(short.MaxValue, -32767, "match")] // blastn 2.12: 32767/−32767 likewise
    public void ComputeUngappedKarlinParameters_ScoreOutsideBlastScoreRange_Throws(int match, int mismatch, string param)
    {
        foreach (var method in new[] { ProbeDesigner.KarlinKMethod.ReducedLattice, ProbeDesigner.KarlinKMethod.NcbiBlast })
            Assert.That(() => ProbeDesigner.ComputeUngappedKarlinParameters(match, mismatch, kMethod: method),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>().With.Property("ParamName").EqualTo(param));
    }

    [Test]
    public void ComputeUngappedKarlinParameters_BlastScoreRangeLimits_Accepted()
    {
        // blastn 2.12.0+ -ungapped -reward 1 -penalty -32767: lambda 1.39, K 0.750, H 1.39.
        var p = ProbeDesigner.ComputeUngappedKarlinParameters(1, short.MinValue + 1);
        Assert.That(p.Lambda, Is.EqualTo(1.39).Within(0.005));
        Assert.That(p.K, Is.EqualTo(0.750).Within(0.0005));
        Assert.That(p.H, Is.EqualTo(1.39).Within(0.005));
    }

    // Audit round 9 (A9-1): the BLAST score-range guard (and the other scheme preconditions) report the PUBLIC
    // method's own parameter name — match/mismatch, reward/penalty, or scoring for a ScoringMatrix argument.
    [TestCase("ComputeLambdaNucleotide", 1, int.MinValue, "mismatch")]
    [TestCase("ComputeLambdaNucleotide", short.MaxValue, -3, "match")]
    [TestCase("ComputeUngappedKarlinParameters", 1, short.MinValue, "mismatch")]
    [TestCase("ComputeUngappedKarlinParameters(frequencies)", 1, int.MinValue, "mismatch")]
    [TestCase("ComputeUngappedKarlinParameters(frequencies)", int.MaxValue, -3, "match")]
    [TestCase("GetBlastnGappedKarlinParameters", 2, int.MinValue, "penalty")]
    [TestCase("GetBlastnGappedKarlinParameters", short.MaxValue, -3, "reward")]
    [TestCase("GetBlastnGappedKarlinParameters", 0, -3, "reward")]
    [TestCase("GetBlastnGappedKarlinParameters", 5, -1, "reward")]   // unsupported scheme (ArgumentException)
    [TestCase("ComputeKarlinAltschul", 1, int.MinValue, "scoring")]
    [TestCase("ComputeKarlinAltschul", short.MaxValue, -3, "scoring")]
    [TestCase("ComputeKarlinAltschul", 4, -1, "scoring")]            // expected score ≥ 0
    [TestCase("ComputeBlastnStatistics", 2, int.MinValue, "scoring")]
    [TestCase("ComputeBlastnStatistics", short.MaxValue, -3, "scoring")]
    [TestCase("ComputeBlastnStatistics", 5, -1, "scoring")]          // unsupported gapped scheme (ArgumentException)
    [TestCase("ComputeBlastnStatistics(ungapped)", 2, int.MinValue, "scoring")]
    [TestCase("ComputeBlastnStatistics(ungapped)", 4, -1, "scoring")] // expected score ≥ 0
    [TestCase("ComputeBlastnStatistics(sequences)", 2, int.MinValue, "scoring")]
    public void KarlinAltschulEntryPoints_InvalidScheme_ThrowWithOwnParamName(string method, int match, int mismatch, string param)
    {
        var m = new Seqeron.Genomics.Infrastructure.ScoringMatrix(Match: match, Mismatch: mismatch, GapOpen: -5, GapExtend: -2);
        TestDelegate call = method switch
        {
            "ComputeLambdaNucleotide" => () => ProbeDesigner.ComputeLambdaNucleotide(match, mismatch),
            "ComputeUngappedKarlinParameters" => () => ProbeDesigner.ComputeUngappedKarlinParameters(match, mismatch),
            "ComputeUngappedKarlinParameters(frequencies)" =>
                () => ProbeDesigner.ComputeUngappedKarlinParameters(match, mismatch, new[] { 0.25, 0.25, 0.25, 0.25 }),
            "GetBlastnGappedKarlinParameters" => () => ProbeDesigner.GetBlastnGappedKarlinParameters(match, mismatch, 5, 2),
            "ComputeKarlinAltschul" => () => ProbeDesigner.ComputeKarlinAltschul(20, 20, 1000, m),
            "ComputeBlastnStatistics" => () => ProbeDesigner.ComputeBlastnStatistics(20, 20, 1000, 1, m),
            "ComputeBlastnStatistics(ungapped)" => () => ProbeDesigner.ComputeBlastnStatistics(20, 20, 1000, 1, m, gapped: false),
            "ComputeBlastnStatistics(sequences)" => () => ProbeDesigner.ComputeBlastnStatistics("ACGTACGT", "ACGTACGT", m),
            _ => throw new ArgumentException(method),
        };
        Assert.That(call, NUnit.Framework.Throws.InstanceOf<ArgumentException>().With.Property("ParamName").EqualTo(param));
    }

    // Audit round 8 (A8-2): Beer–Lambert c = A/(ε·l) in µM; ε, l ≤ 0 or non-finite arguments are undefined.
    [TestCase(0.5, 200000.0, 1.0, 2.5)]
    [TestCase(1.0, 100000.0, 0.5, 20.0)]
    [TestCase(-0.01, 100000.0, 1.0, -0.1)]
    public void CalculateConcentration_BeerLambert_Micromolar(double a, double eps, double l, double expected)
    {
        Assert.That(ProbeDesigner.CalculateConcentration(a, eps, l), Is.EqualTo(expected).Within(1e-12));
    }

    [TestCase(0.5, 0.0, 1.0, "extinctionCoefficient")]
    [TestCase(0.5, -10000.0, 1.0, "extinctionCoefficient")]
    [TestCase(0.5, double.NaN, 1.0, "extinctionCoefficient")]
    [TestCase(0.5, 10000.0, 0.0, "pathLength")]
    [TestCase(0.5, 10000.0, -1.0, "pathLength")]
    [TestCase(0.5, 10000.0, double.PositiveInfinity, "pathLength")]
    [TestCase(double.NaN, 10000.0, 1.0, "absorbance260")]
    public void CalculateConcentration_UndefinedArguments_Throw(double a, double eps, double l, string param)
    {
        Assert.That(() => ProbeDesigner.CalculateConcentration(a, eps, l),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>().With.Property("ParamName").EqualTo(param));
    }
}
