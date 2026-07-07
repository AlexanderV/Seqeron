// ONCO-ARTIFACT-001 — Sequencing Artifact Detection (OxoG / FFPE deamination / strand bias)
// Evidence: docs/Evidence/ONCO-ARTIFACT-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-ARTIFACT-001.md
// Source: Chen L. et al. (2017). Science 355(6326):752-756. https://www.science.org/doi/10.1126/science.aai8690
//         Nature Methods (2017) 14:330 (GIV thresholds). https://www.nature.com/articles/nmeth.4254
//         Do H., Dobrovic A. (2015). FFPE deamination C:G>T:A. https://www.sciencedirect.com/science/article/pii/S152515781630188X
//         Broad GATK FisherStrand / StrandBiasTest. https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/StrandBiasTest.java

using Obs = Seqeron.Genomics.Oncology.OncologyAnalyzer.ArtifactObservation;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_FilterArtifacts_Tests
{
    // Balanced strand evidence shared by classification tests (strand bias is not the subject there).
    private static Obs Variant(char refBase, char altBase, int r1 = 0, int r2 = 0) =>
        new(refBase, altBase, RefForward: 10, RefReverse: 10, AltForward: 10, AltReverse: 10,
            AltReadsR1: r1, AltReadsR2: r2);

    #region ClassifyArtifact (substitution classes)

    // M1 — FFPE cytosine deamination: C>T (Do & Dobrovic 2015, C:G>T:A).
    [Test]
    public void ClassifyArtifact_CtoT_ClassifiesAsFfpeDeamination()
    {
        var call = OncologyAnalyzer.ClassifyArtifact(Variant('C', 'T'));

        Assert.Multiple(() =>
        {
            Assert.That(call.Type, Is.EqualTo(OncologyAnalyzer.ArtifactType.FfpeDeamination),
                "C>T is the canonical FFPE cytosine-deamination substitution (uracil pairs with adenine).");
            Assert.That(call.IsArtifact, Is.True, "FFPE deamination is flagged by substitution class alone.");
        });
    }

    // M2 — FFPE deamination on the antisense strand: G>A (Do & Dobrovic 2015).
    [Test]
    public void ClassifyArtifact_GtoA_ClassifiesAsFfpeDeamination()
    {
        var call = OncologyAnalyzer.ClassifyArtifact(Variant('G', 'A'));

        Assert.That(call.Type, Is.EqualTo(OncologyAnalyzer.ArtifactType.FfpeDeamination),
            "G>A is the antisense-strand read-out of C>T deamination (C:G>T:A).");
    }

    // M3 — OxoG oxidation: G>T (Chen 2017, G>T excess in read 1).
    [Test]
    public void ClassifyArtifact_GtoT_ClassifiesAsOxoG()
    {
        // GIV = 200/100 = 2.0 > 1.5 (damaged), so the OxoG class is confirmed and flagged.
        var call = OncologyAnalyzer.ClassifyArtifact(Variant('G', 'T', r1: 200, r2: 100));

        Assert.Multiple(() =>
        {
            Assert.That(call.Type, Is.EqualTo(OncologyAnalyzer.ArtifactType.OxoG),
                "G>T is the canonical 8-oxoguanine oxidation substitution.");
            Assert.That(call.IsArtifact, Is.True, "OxoG with GIV 2.0 (> 1.5) is a damaged-library artifact.");
        });
    }

    // M4 — OxoG oxidation reverse complement: C>A (Chen 2017, C>A excess in read 2).
    [Test]
    public void ClassifyArtifact_CtoA_ClassifiesAsOxoG()
    {
        var call = OncologyAnalyzer.ClassifyArtifact(Variant('C', 'A', r1: 200, r2: 100));

        Assert.That(call.Type, Is.EqualTo(OncologyAnalyzer.ArtifactType.OxoG),
            "C>A is the reverse complement of G>T, the read-2 read-out of OxoG.");
    }

    // M5 — A>G is neither deamination nor oxidation (the artifact classes are substitution-specific).
    [Test]
    public void ClassifyArtifact_AtoG_ClassifiesAsNoneAndNotArtifact()
    {
        var call = OncologyAnalyzer.ClassifyArtifact(Variant('A', 'G', r1: 200, r2: 100));

        Assert.Multiple(() =>
        {
            Assert.That(call.Type, Is.EqualTo(OncologyAnalyzer.ArtifactType.None),
                "A>G is outside both artifact classes {C>T,G>A} and {G>T,C>A}.");
            Assert.That(call.IsArtifact, Is.False, "A non-artifact substitution is never flagged.");
        });
    }

    // M5 (extension) — OxoG class but undamaged GIV (≤ 1.5) is NOT flagged: class alone is insufficient.
    [Test]
    public void ClassifyArtifact_GtoT_LowGiv_NotFlaggedAsArtifact()
    {
        // GIV = 100/100 = 1.0 (undamaged): OxoG class present but no read-orientation imbalance.
        var call = OncologyAnalyzer.ClassifyArtifact(Variant('G', 'T', r1: 100, r2: 100));

        Assert.Multiple(() =>
        {
            Assert.That(call.Type, Is.EqualTo(OncologyAnalyzer.ArtifactType.OxoG), "Still the OxoG substitution class.");
            Assert.That(call.GivScore, Is.EqualTo(1.0).Within(1e-10), "GIV = 100/100 = 1.0 (undamaged).");
            Assert.That(call.IsArtifact, Is.False, "OxoG below the 1.5 GIV threshold is not a damaged-library artifact.");
        });
    }

    // INV-04 — disjoint classes (lowercase input also accepted; classification is case-insensitive).
    [Test]
    public void ClassifyArtifact_LowercaseBases_ClassifiesSameAsUppercase()
    {
        var call = OncologyAnalyzer.ClassifyArtifact(Variant('c', 't'));

        Assert.That(call.Type, Is.EqualTo(OncologyAnalyzer.ArtifactType.FfpeDeamination),
            "Base classification must be case-insensitive (c>t == C>T).");
    }

    #endregion

    #region CalculateGivScore

    // M6 — GIV = R1/R2 = 200/100 = 2.0, above the 1.5 damaged threshold (Chen 2017 / Nature Methods).
    [Test]
    public void CalculateGivScore_R1TwiceR2_ReturnsTwoAndExceedsDamagedThreshold()
    {
        double giv = OncologyAnalyzer.CalculateGivScore(200, 100);

        Assert.Multiple(() =>
        {
            Assert.That(giv, Is.EqualTo(2.0).Within(1e-10), "GIV = read1/read2 = 200/100 = 2.0.");
            Assert.That(giv, Is.GreaterThan(OncologyAnalyzer.DamagedGivThreshold),
                "GIV 2.0 exceeds the documented damaged threshold of 1.5.");
        });
    }

    // M7 — GIV = R1/R2 = 100/100 = 1.0 (undamaged; balanced read 1 and read 2).
    [Test]
    public void CalculateGivScore_BalancedReads_ReturnsOne()
    {
        double giv = OncologyAnalyzer.CalculateGivScore(100, 100);

        Assert.That(giv, Is.EqualTo(OncologyAnalyzer.UndamagedGivScore).Within(1e-10),
            "Balanced read1 = read2 gives GIV = 1.0, the undamaged baseline (Chen 2017).");
    }

    // C2 — GIV is directional: R1 < R2 gives GIV < 1 (no OxoG imbalance on this strand orientation).
    [Test]
    public void CalculateGivScore_R1HalfR2_ReturnsHalf()
    {
        double giv = OncologyAnalyzer.CalculateGivScore(100, 200);

        Assert.That(giv, Is.EqualTo(0.5).Within(1e-10),
            "GIV = 100/200 = 0.5 (< 1): the imbalance favours read 2, not the OxoG (read-1) direction.");
    }

    // S1 — read-2 count zero with read-1 support: maximal one-sided imbalance, no exception.
    [Test]
    public void CalculateGivScore_ZeroR2_ReturnsPositiveInfinity()
    {
        double giv = OncologyAnalyzer.CalculateGivScore(50, 0);

        Assert.That(double.IsPositiveInfinity(giv), Is.True,
            "All alt reads in read 1 and none in read 2 is a maximal imbalance (GIV = +inf).");
    }

    // S2 — both counts zero: no imbalance evidence, GIV defined as 1.0 (not NaN/exception).
    [Test]
    public void CalculateGivScore_BothZero_ReturnsOne()
    {
        double giv = OncologyAnalyzer.CalculateGivScore(0, 0);

        Assert.That(giv, Is.EqualTo(1.0).Within(1e-10),
            "No alt reads in either mate means no imbalance evidence: GIV = 1.0 (undamaged).");
    }

    [Test]
    public void CalculateGivScore_NegativeR1_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateGivScore(-1, 10),
            "Negative read counts are invalid.");
    }

    // Negative read-2 count is rejected too (mirror branch of the read-1 validation).
    [Test]
    public void CalculateGivScore_NegativeR2_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateGivScore(10, -1),
            "Negative read-2 count is invalid.");
    }

    #endregion

    #region CalculateStrandBias (FisherStrand FS)

    // M8 — Balanced table [10,10,10,10]: two-sided Fisher p = 1 => FS = -10*log10(1) = 0 (GATK).
    [Test]
    public void CalculateStrandBias_BalancedTable_ReturnsZero()
    {
        double fs = OncologyAnalyzer.CalculateStrandBias(10, 10, 10, 10);

        Assert.That(fs, Is.EqualTo(0.0).Within(1e-9),
            "A perfectly balanced strand table has Fisher p = 1, so FS = -10*log10(1) = 0.");
    }

    // M9 — Fully segregated table [20,0,0,20]: exact two-sided Fisher p = 1.4508889103849754e-11,
    // so FS = -10*log10(p) = 108.38365838736458. The expected p/FS were derived independently from the
    // hypergeometric two-sided Fisher formula (not from this implementation's output).
    [Test]
    public void CalculateStrandBias_FullySegregatedTable_ReturnsExactPhredScore()
    {
        double fs = OncologyAnalyzer.CalculateStrandBias(20, 0, 0, 20);

        Assert.That(fs, Is.EqualTo(108.38365838736458).Within(1e-6),
            "Independently computed: two-sided Fisher p = 1.4508889e-11 -> FS = -10*log10(p) = 108.3836584.");
    }

    // M9 (extension) — Partially biased table [15,5,5,15]: exact p = 0.0038475273083775634,
    // FS = 24.148182890180962 (independently derived from the hypergeometric formula).
    [Test]
    public void CalculateStrandBias_PartiallyBiasedTable_ReturnsExactPhredScore()
    {
        double fs = OncologyAnalyzer.CalculateStrandBias(15, 5, 5, 15);

        Assert.That(fs, Is.EqualTo(24.148182890180962).Within(1e-6),
            "Independently computed: two-sided Fisher p = 0.0038475273 -> FS = 24.1481829.");
    }

    // Edge case — an all-zero strand table (no reads at all) provides no evidence of strand bias:
    // p = 1 => FS = 0. Documented in the spec (§3.3: "An all-zero strand table yields p = 1, FS = 0").
    [Test]
    public void CalculateStrandBias_EmptyTable_ReturnsZero()
    {
        double fs = OncologyAnalyzer.CalculateStrandBias(0, 0, 0, 0);

        Assert.That(fs, Is.EqualTo(0.0).Within(1e-9),
            "A table with no reads has no evidence of strand bias: two-sided Fisher p = 1, FS = 0.");
    }

    // Edge case — a zero-margin table (no reverse reads at all): the two-sided Fisher p = 1 (no bias),
    // so FS = 0. Independently confirmed with scipy.stats.fisher_exact([[10,0],[10,0]]) = 1.0.
    [Test]
    public void CalculateStrandBias_ZeroMarginTable_ReturnsZero()
    {
        double fs = OncologyAnalyzer.CalculateStrandBias(10, 0, 10, 0);

        Assert.That(fs, Is.EqualTo(0.0).Within(1e-9),
            "All reads on the forward strand for both alleles gives no strand bias: p = 1, FS = 0.");
    }

    // C1 / INV-05 — monotonicity: increasing strand segregation does not decrease FS.
    [Test]
    public void CalculateStrandBias_IncreasingSegregation_FsIsNonDecreasing()
    {
        double balanced = OncologyAnalyzer.CalculateStrandBias(10, 10, 10, 10);
        double partial = OncologyAnalyzer.CalculateStrandBias(15, 5, 5, 15);
        double full = OncologyAnalyzer.CalculateStrandBias(20, 0, 0, 20);

        Assert.Multiple(() =>
        {
            Assert.That(partial, Is.GreaterThanOrEqualTo(balanced),
                "More strand segregation cannot lower the strand-bias score (Fisher p decreases).");
            Assert.That(full, Is.GreaterThan(partial),
                "A fully segregated table is more extreme than a partially biased one.");
        });
    }

    // S5 — negative strand count is invalid.
    [Test]
    public void CalculateStrandBias_NegativeCount_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateStrandBias(-1, 10, 10, 10),
            "Negative strand read counts are invalid.");
    }

    #endregion

    #region DetectOxoGArtifacts

    // M10 — only the OxoG-class variant with GIV > 1.5 is returned.
    [Test]
    public void DetectOxoGArtifacts_MixedVariants_ReturnsOnlyDamagedOxoG()
    {
        var variants = new[]
        {
            Variant('G', 'T', r1: 200, r2: 100), // OxoG, GIV 2.0 -> detected
            Variant('C', 'T', r1: 200, r2: 100), // FFPE deamination -> not OxoG
            Variant('G', 'T', r1: 100, r2: 100), // OxoG class but GIV 1.0 -> not damaged
            Variant('A', 'G', r1: 200, r2: 100), // non-artifact
        };

        var oxoG = OncologyAnalyzer.DetectOxoGArtifacts(variants);

        Assert.Multiple(() =>
        {
            Assert.That(oxoG, Has.Count.EqualTo(1), "Only one variant is a damaged OxoG artifact.");
            Assert.That(oxoG[0].Type, Is.EqualTo(OncologyAnalyzer.ArtifactType.OxoG), "The returned call is OxoG.");
            Assert.That(oxoG[0].GivScore, Is.EqualTo(2.0).Within(1e-10), "Its GIV is 2.0 (> 1.5).");
        });
    }

    // Empty variant set yields an empty OxoG result (nothing to detect).
    [Test]
    public void DetectOxoGArtifacts_Empty_ReturnsEmpty()
    {
        var oxoG = OncologyAnalyzer.DetectOxoGArtifacts(System.Array.Empty<Obs>());

        Assert.That(oxoG, Is.Empty, "An empty variant set produces no OxoG calls.");
    }

    [Test]
    public void DetectOxoGArtifacts_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.DetectOxoGArtifacts(null!));
    }

    #endregion

    #region FilterArtifacts

    // M11 / INV-01 — artifacts removed, real variant kept; result is a subset in input order.
    [Test]
    public void FilterArtifacts_RemovesArtifactsKeepsRealVariant()
    {
        var real = Variant('A', 'G', r1: 50, r2: 50);       // non-artifact, must survive
        var variants = new[]
        {
            Variant('C', 'T', r1: 30, r2: 30),               // FFPE deamination -> removed
            real,                                            // kept
            Variant('G', 'T', r1: 200, r2: 100),             // OxoG GIV 2.0 -> removed
        };

        var kept = OncologyAnalyzer.FilterArtifacts(variants);

        Assert.Multiple(() =>
        {
            Assert.That(kept, Has.Count.EqualTo(1), "Both artifacts (FFPE C>T and damaged OxoG G>T) are removed.");
            Assert.That(kept[0], Is.EqualTo(real), "The surviving variant is the non-artifact A>G.");
            Assert.That(kept.All(v => variants.Contains(v)), Is.True, "INV-01: result is a subset of the input.");
        });
    }

    // FilterArtifacts keeps an OxoG-class variant whose GIV is undamaged (≤ 1.5).
    [Test]
    public void FilterArtifacts_KeepsOxoGClassVariantWithUndamagedGiv()
    {
        var variants = new[] { Variant('G', 'T', r1: 100, r2: 100) }; // OxoG class, GIV 1.0 (undamaged)

        var kept = OncologyAnalyzer.FilterArtifacts(variants);

        Assert.That(kept, Has.Count.EqualTo(1),
            "An OxoG-class substitution without read-orientation imbalance is not filtered as an artifact.");
    }

    // S4 — empty input yields empty output.
    [Test]
    public void FilterArtifacts_Empty_ReturnsEmpty()
    {
        var kept = OncologyAnalyzer.FilterArtifacts(System.Array.Empty<Obs>());

        Assert.That(kept, Is.Empty, "Filtering an empty variant set returns an empty set.");
    }

    // S3 — null input throws.
    [Test]
    public void FilterArtifacts_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.FilterArtifacts(null!));
    }

    #endregion
}
