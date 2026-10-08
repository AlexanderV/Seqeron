// CODON-RARE-001 — Rare-codon cluster / run detection (%MinMax + Sherlocc RCC)
// Evidence: docs/Evidence/CODON-RARE-001-Evidence.md
// TestSpec: tests/TestSpecs/CODON-RARE-001.md
// Source: Clarke TF, Clark PL (2008). PLoS ONE 3(10):e3412 (%MinMax).
// Source: Chartier M, Gaudreault F, Najmanovich R (2012). Bioinformatics 28(11):1438-1445,
//         doi:10.1093/bioinformatics/bts149 (Sherlocc rare-codon cluster rule: 7-codon window,
//         >=4 rare/pause positions).

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the opt-in rare-codon cluster / run detection added to
/// <see cref="CodonOptimizer"/>: <c>CalculateMinMaxProfile</c> (Clarke &amp; Clark 2008 %MinMax)
/// and <c>FindRareCodonClusters</c> (Chartier et al. 2012 Sherlocc rule). The per-codon
/// <c>FindRareCodons</c> behaviour is unchanged and covered separately.
/// </summary>
[TestFixture]
public class CodonOptimizer_RareCodonClusters_Tests
{
    // E. coli K12 Arginine family (Kazusa species=316407, as embedded in CodonOptimizer.EColiK12):
    //   CGU=0.38, CGC=0.40, CGA=0.06, CGG=0.10, AGA=0.04, AGG=0.02.
    //   Sum = 1.00 -> Xavg = 1.00/6 = 0.16666666..., Xmax = 0.40 (CGC), Xmin = 0.02 (AGG).

    #region CalculateMinMaxProfile — %MinMax (Clarke & Clark 2008)

    // MM1 — Pure rarest-run window yields a strongly negative %Min.
    // Window of 3 AGA codons (Xij=0.04 < Xavg=0.16667): %Min =
    //   (Xavg-Xij)/(Xavg-Xmin) * 100 = (0.16667-0.04)/(0.16667-0.02)*100 = 86.363636..%,
    // returned negative (rare side).
    [Test]
    public void CalculateMinMaxProfile_AllRareArginineWindow_ReturnsExactNegativeMin()
    {
        // 3 codons, window 3 -> exactly one window.
        string sequence = "AGAAGAAGA";

        var profile = CodonOptimizer.CalculateMinMaxProfile(sequence, CodonOptimizer.EColiK12, windowSize: 3);

        Assert.Multiple(() =>
        {
            Assert.That(profile, Has.Count.EqualTo(1), "One window for 3 codons with windowSize 3");
            Assert.That(profile[0].WindowStartCodon, Is.EqualTo(0));
            Assert.That(profile[0].PercentMinMax, Is.EqualTo(-86.36363636363637).Within(1e-10),
                "%Min for an all-AGA window = -(0.16667-0.04)/(0.16667-0.02)*100");
        });
    }

    // MM2 — Pure most-common-codon run yields +100 (%Max). Window of 3 CGC (Xij=0.40=Xmax).
    // %Max = (Xij-Xavg)/(Xmax-Xavg)*100 = (0.40-0.16667)/(0.40-0.16667)*100 = 100.
    [Test]
    public void CalculateMinMaxProfile_AllMostCommonCodonWindow_ReturnsExactPlus100()
    {
        string sequence = "CGCCGCCGC"; // 3x CGC (Arg, most common)

        var profile = CodonOptimizer.CalculateMinMaxProfile(sequence, CodonOptimizer.EColiK12, windowSize: 3);

        Assert.Multiple(() =>
        {
            Assert.That(profile, Has.Count.EqualTo(1));
            Assert.That(profile[0].PercentMinMax, Is.EqualTo(100.0).Within(1e-10),
                "A window encoded only with the most common synonymous codon is +100 %Max");
        });
    }

    // MM3 — Mixed window: CUG (Leu, 0.50) + AGA (Arg, 0.04), window 2.
    // Leu family sum=1.00, avg=0.16667, max=0.50. Arg avg=0.16667, max=0.40.
    // sumXij=0.54 > sumXavg=0.33333 -> %Max = (0.54-0.33333)/((0.50-0.16667)+(0.40-0.16667))*100
    //   = 0.206667/0.566667*100 = 36.470588235..%.
    [Test]
    public void CalculateMinMaxProfile_MixedCommonRareWindow_ReturnsExactMax()
    {
        string sequence = "CUGAGA"; // Leu(CUG) + Arg(AGA)

        var profile = CodonOptimizer.CalculateMinMaxProfile(sequence, CodonOptimizer.EColiK12, windowSize: 2);

        Assert.Multiple(() =>
        {
            Assert.That(profile, Has.Count.EqualTo(1));
            Assert.That(profile[0].PercentMinMax, Is.EqualTo(36.470588235294116).Within(1e-10),
                "Mixed window where summed usage exceeds the average resolves to %Max");
        });
    }

    // MM4 — Window slides one codon at a time: codonCount - windowSize + 1 windows, in order.
    [Test]
    public void CalculateMinMaxProfile_SlidingWindow_ProducesOneWindowPerStart()
    {
        // 5 codons, window 3 -> 3 windows at starts 0,1,2.
        string sequence = "CGCCGCAGAAGAAGA";

        var profile = CodonOptimizer.CalculateMinMaxProfile(sequence, CodonOptimizer.EColiK12, windowSize: 3);

        Assert.Multiple(() =>
        {
            Assert.That(profile, Has.Count.EqualTo(3), "5 codons, window 3 -> 3 sliding windows");
            Assert.That(profile.Select(w => w.WindowStartCodon), Is.EqualTo(new[] { 0, 1, 2 }));
            // Window 2 is all-AGA -> the validated -86.3636.. %Min value.
            Assert.That(profile[2].PercentMinMax, Is.EqualTo(-86.36363636363637).Within(1e-10),
                "Last window (3x AGA) reproduces the all-rare %Min value");
        });
    }

    // MM5 — Sequence shorter than the window produces no windows.
    [Test]
    public void CalculateMinMaxProfile_SequenceShorterThanWindow_ReturnsEmpty()
    {
        string sequence = "AGAAGA"; // 2 codons

        var profile = CodonOptimizer.CalculateMinMaxProfile(sequence, CodonOptimizer.EColiK12, windowSize: 18);

        Assert.That(profile, Is.Empty, "Fewer complete codons than the window -> no %MinMax values");
    }

    // MM6 — Empty / null input returns empty.
    [Test]
    public void CalculateMinMaxProfile_EmptyOrNull_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CodonOptimizer.CalculateMinMaxProfile("", CodonOptimizer.EColiK12), Is.Empty);
            Assert.That(CodonOptimizer.CalculateMinMaxProfile(null!, CodonOptimizer.EColiK12), Is.Empty);
        });
    }

    // MM7 — Single-codon amino acids (Met/Trp) have no synonymous spread and contribute 0 to
    // both numerator and denominator, so a pure Met/Trp window is 0% (not a divide-by-zero / NaN).
    [Test]
    public void CalculateMinMaxProfile_SingleCodonAminoAcids_ReturnsZeroNotNaN()
    {
        string sequence = "AUGAUGUGG"; // Met, Met, Trp

        var profile = CodonOptimizer.CalculateMinMaxProfile(sequence, CodonOptimizer.EColiK12, windowSize: 3);

        Assert.Multiple(() =>
        {
            Assert.That(profile, Has.Count.EqualTo(1));
            Assert.That(double.IsNaN(profile[0].PercentMinMax), Is.False, "Must not be NaN");
            Assert.That(profile[0].PercentMinMax, Is.EqualTo(0.0).Within(1e-10),
                "No synonymous spread -> %MinMax is 0");
        });
    }

    // MM8 — %MinMax is bounded in [-100, 100] for any valid window (INV).
    [Test]
    public void CalculateMinMaxProfile_Invariant_Bounded()
    {
        string sequence = "AUGCGCCGUAGAAGGCGAUUAGCCGGCAGUGAAACCGGUGCUGAUCGC";

        var profile = CodonOptimizer.CalculateMinMaxProfile(sequence, CodonOptimizer.EColiK12, windowSize: 5);

        Assert.That(profile.All(w => w.PercentMinMax >= -100.0 - 1e-9 && w.PercentMinMax <= 100.0 + 1e-9),
            Is.True, "Every %MinMax value must lie in [-100, 100]");
    }

    // MM9 — DNA input (T) is normalised to RNA (U) before scoring.
    [Test]
    public void CalculateMinMaxProfile_DnaInput_NormalisedToRna()
    {
        // DNA spelling of 3x AGA: AGAAGAAGA has no T; use CGA (DNA) vs CGA(RNA identical).
        // Use a T-bearing codon: CTG (DNA) == CUG (RNA) Leu common.
        string dna = "CTGCTGCTG"; // 3x Leu CUG

        var profile = CodonOptimizer.CalculateMinMaxProfile(dna, CodonOptimizer.EColiK12, windowSize: 3);
        var rna = CodonOptimizer.CalculateMinMaxProfile("CUGCUGCUG", CodonOptimizer.EColiK12, windowSize: 3);

        Assert.Multiple(() =>
        {
            Assert.That(profile, Has.Count.EqualTo(1));
            Assert.That(profile[0].PercentMinMax, Is.EqualTo(rna[0].PercentMinMax).Within(1e-12),
                "DNA and RNA spellings of the same window score identically");
        });
    }

    // MM10 — windowSize < 1 throws.
    [Test]
    public void CalculateMinMaxProfile_WindowSizeZero_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => CodonOptimizer.CalculateMinMaxProfile("AGAAGA", CodonOptimizer.EColiK12, windowSize: 0));
    }

    // MM11 — windowSize = 1 is the SMALLEST legal window and must be accepted (pins "windowSize < 1"
    // against the "<= 1" off-by-one; a 1-codon window yields one profile point per codon).
    [Test]
    public void CalculateMinMaxProfile_WindowSizeOne_IsAcceptedAndPerCodon()
    {
        IReadOnlyList<CodonOptimizer.MinMaxWindow> profile = null!;
        Assert.DoesNotThrow(
            () => profile = CodonOptimizer.CalculateMinMaxProfile("AGAAGA", CodonOptimizer.EColiK12, windowSize: 1));
        Assert.That(profile, Has.Count.EqualTo(2), "2 codons, window 1 -> one window per codon");
    }

    #endregion

    #region FindRareCodonClusters — Sherlocc RCC rule (Chartier et al. 2012)

    // C1 — Canonical Sherlocc default: 7-codon window needs >=4 rare codons.
    // 7 consecutive AGA (all rare, freq 0.04 < 0.15) -> exactly one cluster, codons 0..6, 7 rare.
    [Test]
    public void FindRareCodonClusters_SevenRareCodons_DetectedAsOneCluster()
    {
        string sequence = string.Concat(Enumerable.Repeat("AGA", 7)); // 7 rare Arg codons

        var clusters = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(1));
            Assert.That(clusters[0].StartCodon, Is.EqualTo(0));
            Assert.That(clusters[0].EndCodon, Is.EqualTo(6), "7-codon window covers codons 0..6");
            Assert.That(clusters[0].RareCount, Is.EqualTo(7));
        });
    }

    // C2 — Below threshold: only 3 rare codons in a 7-codon window -> no cluster (4 required).
    [Test]
    public void FindRareCodonClusters_ThreeRareInSeven_NoCluster()
    {
        // 3 rare (AGA) + 4 common (CGC, freq 0.40) = 3 < 4.
        string sequence = "AGAAGAAGACGCCGCCGCCGC"; // 7 codons

        var clusters = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();

        Assert.That(clusters, Is.Empty, "3 rare codons < 4 required by the Sherlocc rule");
    }

    // C3 — Exactly 4 rare out of 7 in a window IS a cluster (boundary, >=4).
    [Test]
    public void FindRareCodonClusters_ExactlyFourRareInSeven_IsCluster()
    {
        // 4 rare (AGA) + 3 common (CGC) within a single 7-codon window.
        string sequence = "AGAAGAAGAAGACGCCGCCGC"; // 7 codons, 4 rare

        var clusters = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(1));
            Assert.That(clusters[0].RareCount, Is.EqualTo(4), "Exactly the minimum of 4 rare positions");
            Assert.That(clusters[0].StartCodon, Is.EqualTo(0));
            Assert.That(clusters[0].EndCodon, Is.EqualTo(6));
        });
    }

    // C4 — Isolated rare codons (per-codon FindRareCodons would flag them) do NOT form a cluster.
    // This is the exact capability gap the unit closes: clusters vs single rare codons.
    [Test]
    public void FindRareCodonClusters_IsolatedRareCodons_NoCluster()
    {
        // Rare codons separated by common ones; no 7-window ever holds 4 rare.
        // AGA C C C C C C  AGA C C C C C C  AGA  (rare at 0, 7, 14)
        string sequence = "AGACGCCGCCGCCGCCGCCGCAGACGCCGCCGCCGCCGCCGCAGA";

        var perCodon = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12).ToList();
        var clusters = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(perCodon, Has.Count.EqualTo(3), "Per-codon detection still flags the 3 isolated rare codons");
            Assert.That(clusters, Is.Empty, "But none of them form a Sherlocc rare-codon cluster");
        });
    }

    // C5 — Overlapping qualifying windows merge into a single maximal cluster.
    // 10 consecutive AGA: windows at starts 0..3 each have >=4 rare; merge to codons 0..9.
    [Test]
    public void FindRareCodonClusters_LongRareRun_MergedIntoOneCluster()
    {
        string sequence = string.Concat(Enumerable.Repeat("AGA", 10)); // 10 rare codons

        var clusters = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(1), "One merged cluster, not four overlapping windows");
            Assert.That(clusters[0].StartCodon, Is.EqualTo(0));
            Assert.That(clusters[0].EndCodon, Is.EqualTo(9), "Cluster spans all 10 rare codons");
            Assert.That(clusters[0].RareCount, Is.EqualTo(10));
        });
    }

    // C6 — Two separated rare runs produce two distinct clusters.
    // Layout: run1 = codons 0..6 (AGA rare), gap = codons 7..14 (CGC common), run2 = codons 15..21.
    // A 7-codon window needs >=4 rare. For run1 the last qualifying window starts at codon 3
    // (covers 3..9: codons 3..6 rare = 4), so cluster 1 spans codons 0..9. For run2 the first
    // qualifying window starts at codon 12 (covers 12..18: codons 15..18 rare = 4), and the last
    // valid window start is 15 (max start = 22-7), so cluster 2 spans codons 12..21. The two
    // clusters are distinct because no window between them holds 4 rare codons.
    [Test]
    public void FindRareCodonClusters_TwoSeparatedRuns_TwoClusters()
    {
        string run = string.Concat(Enumerable.Repeat("AGA", 7));
        string gap = string.Concat(Enumerable.Repeat("CGC", 8)); // 8 common codons
        string sequence = run + gap + run; // 22 codons total

        var clusters = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(2), "Two separated rare runs -> two clusters");
            Assert.That(clusters[0].StartCodon, Is.EqualTo(0));
            Assert.That(clusters[0].EndCodon, Is.EqualTo(9), "Cluster 1 = codons 0..9 (last qualifying window starts at 3)");
            Assert.That(clusters[0].RareCount, Is.EqualTo(7), "Only the 7 run-1 codons are rare in 0..9");
            Assert.That(clusters[1].StartCodon, Is.EqualTo(12), "Cluster 2 = first qualifying window starts at codon 12");
            Assert.That(clusters[1].EndCodon, Is.EqualTo(21));
            Assert.That(clusters[1].RareCount, Is.EqualTo(7), "Only the 7 run-2 codons are rare in 12..21");
        });
    }

    // C7 — Default rare threshold is the same 0.15 used by FindRareCodons.
    // AUA (Ile, 0.07) is rare at 0.15; build a 7-codon all-AUA run.
    [Test]
    public void FindRareCodonClusters_DefaultThreshold_MatchesPerCodonCriterion()
    {
        string sequence = string.Concat(Enumerable.Repeat("AUA", 7)); // Ile AUA freq 0.07 < 0.15

        var clusters = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(clusters, Has.Count.EqualTo(1));
            Assert.That(clusters[0].RareCount, Is.EqualTo(7), "AUA (0.07) is rare under the default 0.15 cutoff");
        });
    }

    // C8 — Organism-specific: AGA is rare in E. coli but common in yeast (0.48), so the same
    // 7-AGA run is a cluster in E. coli and not in yeast.
    [Test]
    public void FindRareCodonClusters_OrganismSpecific_DiffersByTable()
    {
        string sequence = string.Concat(Enumerable.Repeat("AGA", 7));

        var ecoli = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();
        var yeast = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.Yeast).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(ecoli, Has.Count.EqualTo(1), "AGA rare in E. coli (0.04) -> cluster");
            Assert.That(yeast, Is.Empty, "AGA common in yeast (0.48) -> no cluster");
        });
    }

    // C9 — Sequence shorter than the cluster window -> no clusters.
    [Test]
    public void FindRareCodonClusters_SequenceShorterThanWindow_ReturnsEmpty()
    {
        string sequence = string.Concat(Enumerable.Repeat("AGA", 6)); // 6 < default window 7

        var clusters = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();

        Assert.That(clusters, Is.Empty);
    }

    // C10 — Empty / null input -> no clusters.
    [Test]
    public void FindRareCodonClusters_EmptyOrNull_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CodonOptimizer.FindRareCodonClusters("", CodonOptimizer.EColiK12), Is.Empty);
            Assert.That(CodonOptimizer.FindRareCodonClusters(null!, CodonOptimizer.EColiK12), Is.Empty);
        });
    }

    // C11 — DNA input is normalised (T->U) before classification.
    [Test]
    public void FindRareCodonClusters_DnaInput_Normalised()
    {
        // ATA (DNA) == AUA (RNA) Ile rare. 7 codons.
        string dna = string.Concat(Enumerable.Repeat("ATA", 7));

        var clusters = CodonOptimizer.FindRareCodonClusters(dna, CodonOptimizer.EColiK12).ToList();

        Assert.That(clusters, Has.Count.EqualTo(1), "ATA -> AUA after normalisation, rare cluster detected");
    }

    // C12 — Tunable parameters: a custom 3-codon window with min 3 detects a short rare run that
    // the Sherlocc default (7/4) would miss.
    [Test]
    public void FindRareCodonClusters_CustomWindowAndThreshold_DetectsShortRun()
    {
        string sequence = string.Concat(Enumerable.Repeat("AGA", 3)); // 3 rare codons

        var defaultRule = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();
        var custom = CodonOptimizer.FindRareCodonClusters(
            sequence, CodonOptimizer.EColiK12, rareThreshold: 0.15, windowSize: 3, minRareCodons: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(defaultRule, Is.Empty, "3 codons < default window 7");
            Assert.That(custom, Has.Count.EqualTo(1), "Custom 3/3 rule detects the short run");
            Assert.That(custom[0].RareCount, Is.EqualTo(3));
        });
    }

    // C13 — Invalid parameters throw.
    [Test]
    public void FindRareCodonClusters_InvalidParameters_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CodonOptimizer.FindRareCodonClusters("AGAAGA", CodonOptimizer.EColiK12, windowSize: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CodonOptimizer.FindRareCodonClusters("AGAAGA", CodonOptimizer.EColiK12, minRareCodons: 0));
        });
    }

    // C13b — windowSize = 1 and minRareCodons = 1 are the SMALLEST legal values and must be accepted,
    // pinning both "windowSize < 1" and "minRareCodons < 1" guards against their "<= 1" off-by-one mutants.
    [Test]
    public void FindRareCodonClusters_WindowAndMinRareOfOne_AreAccepted()
    {
        Assert.DoesNotThrow(() =>
            CodonOptimizer.FindRareCodonClusters(
                "AGAAGA", CodonOptimizer.EColiK12, windowSize: 1, minRareCodons: 1));
    }

    // C14 — Determinism: same input -> same clusters.
    [Test]
    public void FindRareCodonClusters_Deterministic()
    {
        string sequence = string.Concat(Enumerable.Repeat("AGA", 7)) + "CGCCGC" + string.Concat(Enumerable.Repeat("AUA", 7));

        var r1 = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();
        var r2 = CodonOptimizer.FindRareCodonClusters(sequence, CodonOptimizer.EColiK12).ToList();

        Assert.That(r1.SequenceEqual(r2), Is.True, "Cluster detection is deterministic");
    }

    #endregion

    #region Review 2026-09 (B02) — reference-implementation locks

    // S. cerevisiae codon usage per thousand codons (Kazusa species=4932), verbatim from the
    // Clark-lab CHARMING repository file ScerCUB.txt (github.com/wrightgs/CHARMING); its
    // per-amino-acid ratios round to the python-codon-tables s_cerevisiae_4932 fractions
    // (stop codons within 0.01, as their per-thousand values are rounded to 0.1).
    private static readonly CodonOptimizer.CodonUsageTable YeastPerThousand = new(
        "S. cerevisiae (per thousand)",
        new Dictionary<string, double>
        {
            { "AAA", 41.9 }, { "AAC", 24.8 }, { "AAG", 30.8 }, { "AAU", 35.7 }, { "ACA", 17.8 }, { "ACC", 12.7 }, { "ACG", 8.0 }, { "ACU", 20.3 },
            { "AGA", 21.3 }, { "AGC", 9.8 }, { "AGG", 9.2 }, { "AGU", 14.2 }, { "AUA", 17.8 }, { "AUC", 17.2 }, { "AUG", 20.9 }, { "AUU", 30.1 },
            { "CAA", 27.3 }, { "CAC", 7.8 }, { "CAG", 12.1 }, { "CAU", 13.6 }, { "CCA", 18.3 }, { "CCC", 6.8 }, { "CCG", 5.3 }, { "CCU", 13.5 },
            { "CGA", 3.0 }, { "CGC", 2.6 }, { "CGG", 1.7 }, { "CGU", 6.4 }, { "CUA", 13.4 }, { "CUC", 5.4 }, { "CUG", 10.5 }, { "CUU", 12.3 },
            { "GAA", 45.6 }, { "GAC", 20.2 }, { "GAG", 19.2 }, { "GAU", 37.6 }, { "GCA", 16.2 }, { "GCC", 12.6 }, { "GCG", 6.2 }, { "GCU", 21.2 },
            { "GGA", 10.9 }, { "GGC", 9.8 }, { "GGG", 6.0 }, { "GGU", 23.9 }, { "GUA", 11.8 }, { "GUC", 11.8 }, { "GUG", 10.8 }, { "GUU", 22.1 },
            { "UAA", 1.1 }, { "UAC", 14.8 }, { "UAG", 0.5 }, { "UAU", 18.8 }, { "UCA", 18.7 }, { "UCC", 14.2 }, { "UCG", 8.6 }, { "UCU", 23.5 },
            { "UGA", 0.7 }, { "UGC", 4.8 }, { "UGG", 10.4 }, { "UGU", 8.1 }, { "UUA", 26.2 }, { "UUC", 18.4 }, { "UUG", 27.2 }, { "UUU", 26.1 },
        },
        new Dictionary<string, string>());

    // %MinMax with the reference's per-thousand input: values from the Clark-lab reference
    // code calculateMinMax (CHARMING.py, Wright et al. 2022), which prints them rounded to
    // 2 dp (0, 33.99 for CUG·AGA; 12.69, 17.03, 15.99 for the 20-codon gene); the unrounded
    // values below come from the same formula in the Python reference.
    [Test]
    public void CalculateMinMaxProfile_PerThousandTable_MatchesClarkLabReference()
    {
        var two = CodonOptimizer.CalculateMinMaxProfile("CUGAGA", YeastPerThousand, windowSize: 2);
        Assert.That(two.Select(w => w.PercentMinMax), Is.EqualTo(new[] { 33.99209486166008 }).Within(1e-9));

        const string gene = "AUGCUGAGAAAAGAAGCUCUGAGGUUUCCCGGAAUAUCUCGAUGGAAGACCCUAGUUCAG";
        var profile = CodonOptimizer.CalculateMinMaxProfile(gene, YeastPerThousand, windowSize: 18);
        Assert.That(profile.Select(w => w.WindowStartCodon), Is.EqualTo(new[] { 0, 1, 2 }));
        Assert.That(profile.Select(w => w.PercentMinMax),
            Is.EqualTo(new[] { 12.68791340950089, 17.028571428571425, 15.991062879029688 }).Within(1e-9));
    }

    // The same CUG·AGA window with the built-in per-amino-acid relative fractions gives a
    // different value (58.62 vs 33.99): relative fractions weight every residue equally, so
    // they reproduce the reference only for single-amino-acid windows (documented contract).
    [Test]
    public void CalculateMinMaxProfile_RelativeFractionPreset_WeightsResiduesEqually()
    {
        var rel = CodonOptimizer.CalculateMinMaxProfile("CUGAGA", CodonOptimizer.Yeast, windowSize: 2);
        Assert.That(rel[0].PercentMinMax, Is.EqualTo(58.620689655172406).Within(1e-9));
    }

    // The stop family {UAA, UAG, UGA} is a synonymous family in the reference (aaDict['*']):
    // UAG·UAG in E. coli (UAG = family minimum 0.07) is −100 %Min.
    [Test]
    public void CalculateMinMaxProfile_StopFamily_IsSynonymousFamily()
    {
        var p = CodonOptimizer.CalculateMinMaxProfile("UAGUAG", CodonOptimizer.EColiK12, windowSize: 2);
        Assert.That(p[0].PercentMinMax, Is.EqualTo(-100.0).Within(1e-9));
    }

    // Ambiguous triplets are never pause positions (no usage frequency). Previously NNN
    // counted as frequency 0 < threshold, so a run of seven NNN was reported as a cluster.
    [Test]
    public void FindRareCodonClusters_AmbiguousTriplets_AreNotPausePositions()
    {
        Assert.That(CodonOptimizer.FindRareCodonClusters(string.Concat(Enumerable.Repeat("NNN", 7)), CodonOptimizer.EColiK12),
            Is.Empty);

        var mixed = CodonOptimizer.FindRareCodonClusters(
            string.Concat(Enumerable.Repeat("AGA", 4)) + "NNNNNNNNN", CodonOptimizer.EColiK12);
        Assert.That(mixed, Is.EqualTo(new[] { new CodonOptimizer.RareCodonCluster(0, 6, 4) }));
    }

    #endregion
}
