// ONCO-CHIP-001 — Clonal Hematopoiesis (CHIP) Filtering for cfDNA Liquid Biopsy
// Evidence: docs/Evidence/ONCO-CHIP-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-CHIP-001.md
// Source: Steensma DP et al. (2015). Blood 126(1):9-16 — CHIP = driver-gene somatic mutation at VAF >= 2% (0.02).
//         Genovese G et al. (2014). NEJM 371(26):2477-2487 — DNMT3A/TET2/ASXL1/PPM1D/JAK2/SF3B1 CH driver genes.
//         Razavi P et al. (2019). Nat Med 25:1928-1937 — matched-WBC subtraction identifies CH (non-tumour) cfDNA variants.
//         Wan JCM et al. (2020). Sci Transl Med 12(548):eaaz8084 — per-locus alt-read evidence.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_FilterCHIP_Tests
{
    // Helper: a cfDNA variant with gene + plasma VAF; locus defaults are display-only unless WBC matching is tested.
    private static OncologyAnalyzer.ChipVariant Var(
        string gene, double vaf, string chrom = "1", int pos = 100,
        string refA = "A", string altA = "T", int altReads = 0) =>
        new(chrom, pos, refA, altA, gene, vaf, altReads);

    #region IdentifyCHIPVariants

    // M1 — Steensma 2015: a driver-gene somatic mutation at VAF >= 0.02 meets the CHIP definition.
    [Test]
    public void IdentifyCHIPVariants_DriverGeneAboveThreshold_FlagsChip()
    {
        // Arrange: DNMT3A is a canonical CH driver gene (Genovese 2014), VAF 0.05 >= 0.02.
        var variants = new[] { Var("DNMT3A", 0.05) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> chip = OncologyAnalyzer.IdentifyCHIPVariants(variants);

        // Assert
        Assert.That(chip, Has.Count.EqualTo(1),
            "A CHIP driver-gene variant at VAF >= 2% is CHIP per Steensma 2015.");
    }

    // M2 — Steensma 2015: VAF below 2% does NOT meet the CHIP definition.
    [Test]
    public void IdentifyCHIPVariants_DriverGeneBelowThreshold_NotChip()
    {
        // Arrange: DNMT3A but VAF 0.01 < 0.02.
        var variants = new[] { Var("DNMT3A", 0.01) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> chip = OncologyAnalyzer.IdentifyCHIPVariants(variants);

        // Assert
        Assert.That(chip, Is.Empty,
            "VAF 0.01 < 0.02 fails the Steensma 2015 CHIP threshold, so the variant is not flagged.");
    }

    // M3 — Steensma 2015: gene not in the CHIP driver panel is not CHIP regardless of VAF.
    [Test]
    public void IdentifyCHIPVariants_NonChipGeneHighVaf_NotChip()
    {
        // Arrange: EGFR is a solid-tumour driver, not a CHIP gene; high VAF must not flag.
        var variants = new[] { Var("EGFR", 0.30) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> chip = OncologyAnalyzer.IdentifyCHIPVariants(variants);

        // Assert
        Assert.That(chip, Is.Empty,
            "EGFR is not a CHIP driver gene; high VAF alone does not make a variant CHIP (Steensma 2015).");
    }

    // M4 — Steensma 2015: threshold is inclusive ("must be >=2%"), so VAF exactly 0.02 is CHIP.
    [Test]
    public void IdentifyCHIPVariants_VafExactlyAtThreshold_FlagsChip()
    {
        // Arrange: TET2 at VAF exactly 0.02 (boundary). A '>' implementation would wrongly drop this.
        var variants = new[] { Var("TET2", 0.02) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> chip = OncologyAnalyzer.IdentifyCHIPVariants(variants);

        // Assert
        Assert.That(chip, Has.Count.EqualTo(1),
            "VAF == 0.02 is CHIP because the Steensma 2015 threshold is inclusive (>= 2%).");
    }

    // M5 — Genovese 2014 / Steensma 2015: every canonical default driver gene is recognized.
    [Test]
    public void IdentifyCHIPVariants_AllCanonicalGenes_AllFlagged()
    {
        // Arrange: one variant per canonical default gene, all at VAF 0.10.
        var genes = new[] { "DNMT3A", "TET2", "ASXL1", "TP53", "JAK2", "SF3B1", "SRSF2", "PPM1D" };
        var variants = genes.Select(g => Var(g, 0.10)).ToArray();

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> chip = OncologyAnalyzer.IdentifyCHIPVariants(variants);

        // Assert: every canonical gene must be returned (not merely the right count) — a bug that
        // dropped one gene but duplicated another would still pass a count-only check.
        Assert.That(chip.Select(v => v.Gene), Is.EquivalentTo(genes),
            "All canonical CH driver genes (Steensma 2015 Fig 2A / Genovese 2014) at VAF >= 0.02 are CHIP.");
    }

    // M6 — Razavi 2019 / Framework: a caller-supplied panel governs membership.
    [Test]
    public void IdentifyCHIPVariants_CallerSuppliedPanel_FlagsCustomGene()
    {
        // Arrange: custom panel contains ABCX (a gene absent from the default set).
        var variants = new[] { Var("ABCX", 0.05) };
        var panel = new[] { "ABCX" };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> chip = OncologyAnalyzer.IdentifyCHIPVariants(variants, panel);

        // Assert
        Assert.That(chip, Has.Count.EqualTo(1),
            "The CHIP panel is caller-supplied; ABCX in the custom panel at VAF >= 0.02 is flagged.");
    }

    // M7 — Steensma 2015 set: a gene outside the default panel is not CHIP.
    [Test]
    public void IdentifyCHIPVariants_DefaultPanel_ExcludesNonMember()
    {
        // Arrange: ABCX is not in the canonical default panel.
        var variants = new[] { Var("ABCX", 0.05) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> chip = OncologyAnalyzer.IdentifyCHIPVariants(variants);

        // Assert
        Assert.That(chip, Is.Empty,
            "ABCX is absent from the canonical default panel, so it is not flagged CHIP.");
    }

    // S1 — HGNC symbols are upper-case; gene comparison must be case-insensitive.
    [Test]
    public void IdentifyCHIPVariants_LowercaseGene_FlagsChip()
    {
        // Arrange: lower-case "dnmt3a".
        var variants = new[] { Var("dnmt3a", 0.05) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> chip = OncologyAnalyzer.IdentifyCHIPVariants(variants);

        // Assert
        Assert.That(chip, Has.Count.EqualTo(1),
            "Gene membership is case-insensitive; 'dnmt3a' matches the canonical DNMT3A.");
    }

    // S3 — threshold is configurable: a custom minVaf above the variant VAF excludes it.
    [Test]
    public void IdentifyCHIPVariants_CustomMinVaf_ExcludesBelowCustom()
    {
        // Arrange: DNMT3A VAF 0.05 with a stricter caller minVaf 0.10.
        var variants = new[] { Var("DNMT3A", 0.05) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> chip = OncologyAnalyzer.IdentifyCHIPVariants(variants, minVaf: 0.10);

        // Assert
        Assert.That(chip, Is.Empty,
            "VAF 0.05 < custom minVaf 0.10, so the variant is below the configured CHIP threshold.");
    }

    // V1 — repository convention: null input throws.
    [Test]
    public void IdentifyCHIPVariants_NullVariants_Throws()
    {
        Assert.That(() => OncologyAnalyzer.IdentifyCHIPVariants(null!),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(),
            "Null variant collection is rejected.");
    }

    // V4 — domain: minVaf must lie in (0, 1].
    [Test]
    public void IdentifyCHIPVariants_MinVafOutOfRange_Throws()
    {
        var variants = new[] { Var("DNMT3A", 0.05) };
        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.IdentifyCHIPVariants(variants, minVaf: 0.0),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(), "minVaf 0 is out of (0, 1].");
            Assert.That(() => OncologyAnalyzer.IdentifyCHIPVariants(variants, minVaf: 1.5),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(), "minVaf 1.5 is out of (0, 1].");
        });
    }

    #endregion

    #region FilterCHIP

    // M8 — Razavi 2019: a cfDNA variant present in the matched WBC is WBC/CH-derived and is removed.
    [Test]
    public void FilterCHIP_VariantInMatchedWbc_Removed()
    {
        // Arrange: EGFR (a non-CHIP gene) at the same locus in both cfDNA and matched WBC (WBC has alt reads).
        var cfDna = new[] { Var("EGFR", 0.30, chrom: "7", pos: 55259515, refA: "T", altA: "G") };
        var wbc = new[] { Var("EGFR", 0.30, chrom: "7", pos: 55259515, refA: "T", altA: "G", altReads: 8) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> kept = OncologyAnalyzer.FilterCHIP(cfDna, wbc);

        // Assert
        Assert.That(kept, Is.Empty,
            "A cfDNA variant matched in WBC is CH-derived and removed regardless of gene (Razavi 2019).");
    }

    // M9 — Razavi 2019: a cfDNA variant absent from matched WBC (and not a CHIP gene) is retained as tumour candidate.
    [Test]
    public void FilterCHIP_VariantAbsentFromWbc_Retained()
    {
        // Arrange: EGFR variant; matched WBC has a different locus only.
        var cfDna = new[] { Var("EGFR", 0.30, chrom: "7", pos: 55259515, refA: "T", altA: "G") };
        var wbc = new[] { Var("DNMT3A", 0.04, chrom: "2", pos: 25234374, refA: "C", altA: "T", altReads: 5) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> kept = OncologyAnalyzer.FilterCHIP(cfDna, wbc);

        // Assert
        Assert.That(kept, Has.Count.EqualTo(1),
            "An EGFR variant absent from matched WBC is retained as a candidate tumour variant (Razavi 2019).");
    }

    // M10 — Steensma 2015 fallback: a CHIP-gene variant meeting gene+VAF is removed even with no WBC evidence.
    [Test]
    public void FilterCHIP_ChipGeneVariantNotInWbc_RemovedByHeuristic()
    {
        // Arrange: DNMT3A VAF 0.05, no matched WBC evidence at all.
        var cfDna = new[] { Var("DNMT3A", 0.05, chrom: "2", pos: 25234374, refA: "C", altA: "T") };
        var wbc = Array.Empty<OncologyAnalyzer.ChipVariant>();

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> kept = OncologyAnalyzer.FilterCHIP(cfDna, wbc);

        // Assert
        Assert.That(kept, Is.Empty,
            "Without WBC evidence, a DNMT3A VAF>=0.02 variant is removed by the gene+VAF CHIP heuristic (Steensma 2015).");
    }

    // M11 — Steensma 2015: a sub-threshold CHIP-gene variant absent from WBC is retained (fails both rules).
    [Test]
    public void FilterCHIP_SubThresholdChipGeneAbsentFromWbc_Retained()
    {
        // Arrange: DNMT3A but VAF 0.01 (< 0.02) and not in WBC.
        var cfDna = new[] { Var("DNMT3A", 0.01, chrom: "2", pos: 25234374, refA: "C", altA: "T") };
        var wbc = Array.Empty<OncologyAnalyzer.ChipVariant>();

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> kept = OncologyAnalyzer.FilterCHIP(cfDna, wbc);

        // Assert
        Assert.That(kept, Has.Count.EqualTo(1),
            "VAF 0.01 fails the CHIP threshold and there is no WBC match, so the variant is retained.");
    }

    // M12 — combined: matched-WBC subtraction + gene heuristic leave only the tumour candidate.
    [Test]
    public void FilterCHIP_MixedPanel_OnlyTumorVariantRetained()
    {
        // Arrange:
        //   EGFR  (tumour, absent from WBC)        -> retained
        //   DNMT3A VAF 0.06 (CHIP gene+VAF)        -> removed by heuristic
        //   KRAS  matched in WBC (alt reads)       -> removed by matched-WBC rule
        var cfDna = new[]
        {
            Var("EGFR", 0.30, chrom: "7", pos: 55259515, refA: "T", altA: "G"),
            Var("DNMT3A", 0.06, chrom: "2", pos: 25234374, refA: "C", altA: "T"),
            Var("KRAS", 0.20, chrom: "12", pos: 25245350, refA: "C", altA: "A"),
        };
        var wbc = new[]
        {
            Var("KRAS", 0.18, chrom: "12", pos: 25245350, refA: "C", altA: "A", altReads: 6),
        };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> kept = OncologyAnalyzer.FilterCHIP(cfDna, wbc);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(kept, Has.Count.EqualTo(1), "Only the EGFR tumour candidate survives both CHIP rules.");
            Assert.That(kept[0].Gene, Is.EqualTo("EGFR"),
                "EGFR is absent from WBC and not a CHIP gene, so it is the only retained variant.");
        });
    }

    // S2 — Wan 2020: a WBC locus with 0 alt reads is not 'present', so the cfDNA variant is not subtracted.
    [Test]
    public void FilterCHIP_WbcLocusZeroAltReads_TreatedAbsent()
    {
        // Arrange: same locus in WBC but with 0 alt reads (no mutant evidence).
        var cfDna = new[] { Var("EGFR", 0.30, chrom: "7", pos: 55259515, refA: "T", altA: "G") };
        var wbc = new[] { Var("EGFR", 0.0, chrom: "7", pos: 55259515, refA: "T", altA: "G", altReads: 0) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> kept = OncologyAnalyzer.FilterCHIP(cfDna, wbc);

        // Assert
        Assert.That(kept, Has.Count.EqualTo(1),
            "A WBC locus with 0 alt reads carries no mutant evidence, so the cfDNA variant is not subtracted (Wan 2020).");
    }

    // C1 — INV-03: retained variants preserve input order.
    [Test]
    public void FilterCHIP_PreservesInputOrder()
    {
        // Arrange: three tumour-candidate variants (none CHIP genes, none in WBC), in a known order.
        var cfDna = new[]
        {
            Var("EGFR", 0.30, chrom: "7", pos: 55259515, refA: "T", altA: "G"),
            Var("KRAS", 0.25, chrom: "12", pos: 25245350, refA: "C", altA: "A"),
            Var("BRAF", 0.20, chrom: "7", pos: 140753336, refA: "A", altA: "T"),
        };
        var wbc = Array.Empty<OncologyAnalyzer.ChipVariant>();

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> kept = OncologyAnalyzer.FilterCHIP(cfDna, wbc);

        // Assert
        Assert.That(kept.Select(v => v.Gene), Is.EqualTo(new[] { "EGFR", "KRAS", "BRAF" }),
            "FilterCHIP preserves the input order of retained variants (INV-03).");
    }

    // V2 — repository convention: null cfDNA input throws.
    [Test]
    public void FilterCHIP_NullVariants_Throws()
    {
        var wbc = Array.Empty<OncologyAnalyzer.ChipVariant>();
        Assert.That(() => OncologyAnalyzer.FilterCHIP(null!, wbc),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(), "Null cfDNA collection is rejected.");
    }

    // V3 — repository convention: null matched-WBC collection throws.
    [Test]
    public void FilterCHIP_NullWbc_Throws()
    {
        var cfDna = new[] { Var("EGFR", 0.30) };
        Assert.That(() => OncologyAnalyzer.FilterCHIP(cfDna, null!),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(), "Null matched-WBC collection is rejected.");
    }

    // V4b — domain: FilterCHIP minVaf must lie in (0, 1] (same contract as IdentifyCHIPVariants).
    [Test]
    public void FilterCHIP_MinVafOutOfRange_Throws()
    {
        var cfDna = new[] { Var("DNMT3A", 0.05) };
        var wbc = Array.Empty<OncologyAnalyzer.ChipVariant>();
        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.FilterCHIP(cfDna, wbc, minVaf: 0.0),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(), "minVaf 0 is out of (0, 1].");
            Assert.That(() => OncologyAnalyzer.FilterCHIP(cfDna, wbc, minVaf: 1.5),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(), "minVaf 1.5 is out of (0, 1].");
        });
    }

    // V7 — domain (contract §3.3): minWbcAltReads must be >= 1; a value below 1 is rejected.
    [Test]
    public void FilterCHIP_MinWbcAltReadsBelowOne_Throws()
    {
        var cfDna = new[] { Var("EGFR", 0.30) };
        var wbc = Array.Empty<OncologyAnalyzer.ChipVariant>();
        Assert.That(() => OncologyAnalyzer.FilterCHIP(cfDna, wbc, minWbcAltReads: 0),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(),
            "minWbcAltReads must be at least 1 (Wan 2020 per-locus alt-read evidence); 0 is rejected.");
    }

    // M8b — Razavi 2019 boundary: the WBC alt-read cutoff is inclusive (>=). A WBC locus with exactly
    // minWbcAltReads alt reads counts as present, so the matched cfDNA variant is subtracted.
    [Test]
    public void FilterCHIP_WbcAltReadsExactlyAtCutoff_Removed()
    {
        // Arrange: matched locus, WBC has exactly 1 alt read (default cutoff = 1).
        var cfDna = new[] { Var("EGFR", 0.30, chrom: "7", pos: 55259515, refA: "T", altA: "G") };
        var wbc = new[] { Var("EGFR", 0.30, chrom: "7", pos: 55259515, refA: "T", altA: "G", altReads: 1) };

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> kept = OncologyAnalyzer.FilterCHIP(cfDna, wbc);

        // Assert
        Assert.That(kept, Is.Empty,
            "A WBC locus with alt reads == cutoff (1) is 'present' (inclusive >=), so the cfDNA variant is removed.");
    }

    // V5 — empty cfDNA input yields an empty result.
    [Test]
    public void FilterCHIP_EmptyVariants_ReturnsEmpty()
    {
        var kept = OncologyAnalyzer.FilterCHIP(
            Array.Empty<OncologyAnalyzer.ChipVariant>(), Array.Empty<OncologyAnalyzer.ChipVariant>());
        Assert.That(kept, Is.Empty, "Filtering an empty cfDNA set yields an empty result.");
    }

    // V6 — Razavi 2019: with no WBC evidence, only the gene+VAF heuristic applies.
    [Test]
    public void FilterCHIP_EmptyWbc_AppliesGeneHeuristicOnly()
    {
        // Arrange: one CHIP-gene variant (removed by heuristic) + one tumour candidate (retained).
        var cfDna = new[]
        {
            Var("TET2", 0.10, chrom: "4", pos: 105234567, refA: "G", altA: "A"),
            Var("EGFR", 0.30, chrom: "7", pos: 55259515, refA: "T", altA: "G"),
        };
        var wbc = Array.Empty<OncologyAnalyzer.ChipVariant>();

        // Act
        IReadOnlyList<OncologyAnalyzer.ChipVariant> kept = OncologyAnalyzer.FilterCHIP(cfDna, wbc);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(kept, Has.Count.EqualTo(1), "Only the non-CHIP tumour candidate survives the gene heuristic.");
            Assert.That(kept[0].Gene, Is.EqualTo("EGFR"), "TET2 is removed by the CHIP gene+VAF heuristic.");
        });
    }

    #endregion

    #region IsCanonicalChipGene

    // Supporting: case-insensitive membership and null/empty handling (used by both canonical methods).
    [Test]
    public void IsCanonicalChipGene_MembershipAndNullHandling()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsCanonicalChipGene("ASXL1"), Is.True,
                "ASXL1 is a canonical CH driver gene (Genovese 2014).");
            Assert.That(OncologyAnalyzer.IsCanonicalChipGene("asxl1"), Is.True,
                "Membership is case-insensitive.");
            Assert.That(OncologyAnalyzer.IsCanonicalChipGene("EGFR"), Is.False,
                "EGFR is not a CHIP driver gene.");
            Assert.That(OncologyAnalyzer.IsCanonicalChipGene(null), Is.False,
                "A null gene is not a CHIP gene.");
            Assert.That(OncologyAnalyzer.IsCanonicalChipGene(""), Is.False,
                "An empty gene symbol is not a CHIP gene.");
        });
    }

    #endregion
}
