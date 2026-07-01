// ONCO-ANNOT-001 — Cancer-Specific Variant Annotation (AMP/ASCO/CAP 2017 tiers)
// Evidence: docs/Evidence/ONCO-ANNOT-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-ANNOT-001.md
// Source: Li MM et al. (2017). Standards and Guidelines for the Interpretation and Reporting of
//         Sequence Variants in Cancer. J Mol Diagn 19(1):4-23. https://doi.org/10.1016/j.jmoldx.2016.10.002
//         Tate JG et al. (2019). COSMIC. Nucleic Acids Res 47(D1):D941-D947. https://doi.org/10.1093/nar/gky1015

using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Oncology;
using Level = Seqeron.Genomics.Oncology.OncologyAnalyzer.ClinicalEvidenceLevel;
using Tier = Seqeron.Genomics.Oncology.OncologyAnalyzer.VariantTier;
using Input = Seqeron.Genomics.Oncology.OncologyAnalyzer.CancerVariantAnnotationInput;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_AnnotateCancerVariants_Tests
{
    // Helper: build a variant evidence record with sensible defaults; tests override what matters.
    private static Input Variant(
        Level level = Level.None,
        double maf = 0.0,
        bool cancerAssociation = false,
        string gene = "GENE",
        string protein = "p.X1Y") =>
        new(gene, protein, level, maf, cancerAssociation);

    #region ClassifyVariantTier — evidence-level mapping (Tier I / II)

    // M1 — Level A ⇒ Tier I. Li (2017) Figure 2: "tier I ... (level A and B evidence)".
    [Test]
    public void ClassifyVariantTier_LevelA_ReturnsTierI()
    {
        // Even with MAF 0 and a cancer association, Level A alone fixes Tier I.
        var tier = OncologyAnalyzer.ClassifyVariantTier(Variant(Level.A, maf: 0.0, cancerAssociation: true));

        Assert.That(tier, Is.EqualTo(Tier.TierI_StrongClinicalSignificance),
            "Level A evidence maps to Tier I (strong clinical significance) per Li 2017 Figure 2.");
    }

    // M2 — Level B ⇒ Tier I. Li (2017) Figure 2.
    [Test]
    public void ClassifyVariantTier_LevelB_ReturnsTierI()
    {
        var tier = OncologyAnalyzer.ClassifyVariantTier(Variant(Level.B, maf: 0.0, cancerAssociation: true));

        Assert.That(tier, Is.EqualTo(Tier.TierI_StrongClinicalSignificance),
            "Level B evidence maps to Tier I per Li 2017 Figure 2 (level A and B = Tier I).");
    }

    // M3 — Level C ⇒ Tier II. Li (2017) Figure 2: "tier II ... (level C or D evidence)".
    [Test]
    public void ClassifyVariantTier_LevelC_ReturnsTierII()
    {
        var tier = OncologyAnalyzer.ClassifyVariantTier(Variant(Level.C, maf: 0.0, cancerAssociation: true));

        Assert.That(tier, Is.EqualTo(Tier.TierII_PotentialClinicalSignificance),
            "Level C evidence maps to Tier II (potential clinical significance) per Li 2017 Figure 2.");
    }

    // M4 — Level D ⇒ Tier II. Li (2017) Figure 2.
    [Test]
    public void ClassifyVariantTier_LevelD_ReturnsTierII()
    {
        var tier = OncologyAnalyzer.ClassifyVariantTier(Variant(Level.D, maf: 0.0, cancerAssociation: true));

        Assert.That(tier, Is.EqualTo(Tier.TierII_PotentialClinicalSignificance),
            "Level D evidence maps to Tier II per Li 2017 Figure 2 (level C or D = Tier II).");
    }

    // M8 — Clinical evidence takes priority over a high population frequency.
    // Li (2017) categorizes by evidence level (Figure 2); a Level A biomarker is Tier I even if common.
    [Test]
    public void ClassifyVariantTier_LevelAWithHighMaf_StillReturnsTierI()
    {
        // MAF 0.30 (> 1% benign cutoff) AND no cancer association would be Tier IV WITHOUT evidence;
        // Level A must override that to Tier I.
        var tier = OncologyAnalyzer.ClassifyVariantTier(Variant(Level.A, maf: 0.30, cancerAssociation: false));

        Assert.That(tier, Is.EqualTo(Tier.TierI_StrongClinicalSignificance),
            "A Level A biomarker is Tier I regardless of population frequency (evidence-level-first, Li 2017).");
    }

    #endregion

    #region ClassifyVariantTier — benign / unknown (Tier IV / III) with no evidence level

    // M5 — Common polymorphism (MAF 0.25, no evidence) ⇒ Tier IV. Li (2017) Table 7 (MAF >= 1%).
    [Test]
    public void ClassifyVariantTier_CommonPolymorphism_ReturnsTierIV()
    {
        // Cancer association true so the ONLY reason for Tier IV is the high MAF (>= 1%).
        var tier = OncologyAnalyzer.ClassifyVariantTier(Variant(Level.None, maf: 0.25, cancerAssociation: true));

        Assert.That(tier, Is.EqualTo(Tier.TierIV_BenignOrLikelyBenign),
            "MAF 0.25 (>= 1% cutoff) with no clinical evidence is a benign polymorphism (Tier IV), Li 2017 Table 7.");
    }

    // M6 — Rare, no evidence, no cancer association ⇒ Tier IV.
    // Li (2017) Figure 2 Tier IV box: "No existing published evidence of cancer association".
    [Test]
    public void ClassifyVariantTier_RareNoCancerAssociation_ReturnsTierIV()
    {
        var tier = OncologyAnalyzer.ClassifyVariantTier(Variant(Level.None, maf: 0.0001, cancerAssociation: false));

        Assert.That(tier, Is.EqualTo(Tier.TierIV_BenignOrLikelyBenign),
            "Rare variant with no cancer association is Tier IV (no published cancer association), Li 2017 Figure 2.");
    }

    // M7 — Rare, no evidence, WITH cancer association ⇒ Tier III (unknown significance).
    // Li (2017) Table 6: absent/extremely-low MAF, no clinical evidence, but cannot be called benign.
    [Test]
    public void ClassifyVariantTier_RareWithCancerAssociation_ReturnsTierIII()
    {
        var tier = OncologyAnalyzer.ClassifyVariantTier(Variant(Level.None, maf: 0.0001, cancerAssociation: true));

        Assert.That(tier, Is.EqualTo(Tier.TierIII_UnknownClinicalSignificance),
            "Rare variant with a cancer association but no clinical evidence is Tier III (VUS), Li 2017 Table 6.");
    }

    // M9 — MAF exactly at the 1% cutoff ⇒ Tier IV (cutoff is ">= 1%", inclusive).
    [Test]
    public void ClassifyVariantTier_MafExactlyAtOnePercent_ReturnsTierIV()
    {
        // 0.01 exactly; cancer association true so high MAF is the sole driver. A strict ">" implementation
        // would wrongly return Tier III here, so this value distinguishes correct (>=) from off-by-one code.
        var tier = OncologyAnalyzer.ClassifyVariantTier(
            Variant(Level.None, maf: OncologyAnalyzer.BenignPopulationMafThreshold, cancerAssociation: true));

        Assert.That(tier, Is.EqualTo(Tier.TierIV_BenignOrLikelyBenign),
            "MAF exactly 0.01 is benign (Tier IV): the 1% primary cutoff is inclusive (Li 2017, '>= 1%').");
    }

    // M10 — MAF just below the cutoff (0.0099) WITH cancer association ⇒ Tier III.
    [Test]
    public void ClassifyVariantTier_MafJustBelowOnePercent_ReturnsTierIII()
    {
        var tier = OncologyAnalyzer.ClassifyVariantTier(Variant(Level.None, maf: 0.0099, cancerAssociation: true));

        Assert.That(tier, Is.EqualTo(Tier.TierIII_UnknownClinicalSignificance),
            "MAF 0.0099 (< 1%) with a cancer association stays Tier III (not a common polymorphism), Li 2017.");
    }

    #endregion

    #region AnnotateCancerVariants — batch

    // M11 — Mixed batch: tiers and order preserved, one annotation per input. INV-4.
    [Test]
    public void AnnotateCancerVariants_MixedBatch_PreservesOrderAndTiers()
    {
        var batch = new[]
        {
            Variant(Level.A, maf: 0.0, cancerAssociation: true,  gene: "BRAF", protein: "p.V600E"), // Tier I
            Variant(Level.C, maf: 0.0, cancerAssociation: true,  gene: "JAK2", protein: "p.V617F"), // Tier II
            Variant(Level.None, maf: 0.25, cancerAssociation: false, gene: "SNP1", protein: "p.A1A"), // Tier IV
            Variant(Level.None, maf: 0.0001, cancerAssociation: true, gene: "VUS1", protein: "p.R2H"), // Tier III
        };

        var result = OncologyAnalyzer.AnnotateCancerVariants(batch);

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(4), "One annotation is produced per input variant (INV-4).");
            Assert.That(result[0].Tier, Is.EqualTo(Tier.TierI_StrongClinicalSignificance), "Index 0: Level A ⇒ Tier I.");
            Assert.That(result[1].Tier, Is.EqualTo(Tier.TierII_PotentialClinicalSignificance), "Index 1: Level C ⇒ Tier II.");
            Assert.That(result[2].Tier, Is.EqualTo(Tier.TierIV_BenignOrLikelyBenign), "Index 2: common SNP ⇒ Tier IV.");
            Assert.That(result[3].Tier, Is.EqualTo(Tier.TierIII_UnknownClinicalSignificance), "Index 3: rare VUS ⇒ Tier III.");
            Assert.That(result[0].Variant.Gene, Is.EqualTo("BRAF"), "Input order/content is preserved on the annotation.");
        });
    }

    // S6 — Empty batch ⇒ empty list.
    [Test]
    public void AnnotateCancerVariants_EmptyBatch_ReturnsEmpty()
    {
        var result = OncologyAnalyzer.AnnotateCancerVariants(Array.Empty<Input>());

        Assert.That(result, Is.Empty, "An empty input batch yields an empty annotation list.");
    }

    // S4 — Null variants ⇒ ArgumentNullException.
    [Test]
    public void AnnotateCancerVariants_NullVariants_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.AnnotateCancerVariants(null!),
            "Null variant collection is rejected per the API contract.");
    }

    #endregion

    #region GetCOSMICAnnotation — caller-supplied catalog lookup

    // M12 — Variant present in the supplied COSMIC catalog ⇒ returns the catalog value.
    [Test]
    public void GetCOSMICAnnotation_VariantInCatalog_ReturnsCosmicId()
    {
        var catalog = new Dictionary<(string Gene, string ProteinChange), string>
        {
            [("BRAF", "p.V600E")] = "COSV56056643",
        };
        var variant = Variant(Level.A, gene: "BRAF", protein: "p.V600E");

        string? id = OncologyAnalyzer.GetCOSMICAnnotation(variant, catalog);

        Assert.That(id, Is.EqualTo("COSV56056643"),
            "An exact (gene, protein change) match returns the caller-supplied COSMIC identifier (Tate 2019).");
    }

    // M13 — Variant absent from the catalog ⇒ null (do not fabricate).
    [Test]
    public void GetCOSMICAnnotation_VariantNotInCatalog_ReturnsNull()
    {
        var catalog = new Dictionary<(string Gene, string ProteinChange), string>
        {
            [("BRAF", "p.V600E")] = "COSV56056643",
        };
        var variant = Variant(gene: "TP53", protein: "p.R175H");

        string? id = OncologyAnalyzer.GetCOSMICAnnotation(variant, catalog);

        Assert.That(id, Is.Null, "A catalog miss returns null; COSMIC content is external and not fabricated.");
    }

    // S5 — Null catalog ⇒ ArgumentNullException.
    [Test]
    public void GetCOSMICAnnotation_NullCatalog_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.GetCOSMICAnnotation(Variant(), null!),
            "Null COSMIC catalog is rejected per the API contract.");
    }

    #endregion

    #region Input validation

    // S1 — Negative MAF ⇒ ArgumentOutOfRangeException.
    [Test]
    public void ClassifyVariantTier_NegativeMaf_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.ClassifyVariantTier(Variant(Level.None, maf: -0.1)),
            "A negative population MAF is invalid (must be in [0, 1]).");
    }

    // S2 — MAF > 1 ⇒ ArgumentOutOfRangeException.
    [Test]
    public void ClassifyVariantTier_MafAboveOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.ClassifyVariantTier(Variant(Level.None, maf: 1.5)),
            "A population MAF above 1 is invalid (must be in [0, 1]).");
    }

    // S3 — NaN MAF ⇒ ArgumentOutOfRangeException.
    [Test]
    public void ClassifyVariantTier_NaNMaf_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.ClassifyVariantTier(Variant(Level.None, maf: double.NaN)),
            "A NaN population MAF is invalid (must be a real number in [0, 1]).");
    }

    #endregion

    #region Property — totality (INV-1)

    // C1 — Every (evidence level × MAF band × cancer association) combination yields a defined tier.
    [Test]
    public void ClassifyVariantTier_AllEvidenceAndMafCombinations_ProduceDefinedTier()
    {
        Level[] levels = { Level.None, Level.A, Level.B, Level.C, Level.D };
        double[] mafs = { 0.0, 0.0099, 0.01, 0.25, 1.0 };
        bool[] associations = { false, true };

        Assert.Multiple(() =>
        {
            foreach (var level in levels)
            foreach (var maf in mafs)
            foreach (var assoc in associations)
            {
                var tier = OncologyAnalyzer.ClassifyVariantTier(Variant(level, maf, assoc));
                Assert.That(Enum.IsDefined(typeof(Tier), tier), Is.True,
                    $"Every input maps to one of the four defined tiers (level={level}, maf={maf}, assoc={assoc}) — INV-1.");
            }
        });
    }

    #endregion
}
