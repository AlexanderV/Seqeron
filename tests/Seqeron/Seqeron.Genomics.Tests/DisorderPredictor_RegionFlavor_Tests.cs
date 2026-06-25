// DISORDER-REGION-001 — Disordered Region Detection (MobiDB-lite flavor labelling)
// Evidence: docs/Evidence/DISORDER-REGION-001-Evidence.md
// TestSpec: tests/TestSpecs/DISORDER-REGION-001.md
// Source: Necci M, Piovesan D, Clementel D, Dosztányi Z, Tosatto SCE (2020).
//         "MobiDB-lite 3.0: fast consensus annotation of intrinsic disorder flavors in proteins".
//         Bioinformatics 36(22-23):5533-5534. DOI 10.1093/bioinformatics/btaa1045. PMID 33325498.
//         Charge classes: Das RK, Pappu RV (2013). PNAS 110(33):13392-13397. PMID 23901099.
// Reference implementation (constants taken verbatim): BioComputingUP/MobiDB-lite (branch v3),
//         mdblib/states.py (get_disorder_class, is_enriched threshold=0.32) and
//         mdblib/consensus.py (get_region_features priority order).

using System;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class DisorderPredictor_RegionFlavor_Tests
{
    #region ClassifyRegionFlavorMobiDbLite — charge classes (Das & Pappu 2013 diagram of states)

    // F1 — PA: f+ = (R+K)/L = 0.5, f- = (D+E)/L = 0.5, FCR = 1.0 > 0.35, NCPR = 0 <= 0.35 → Polyampholyte.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_BalancedStrongCharge_ReturnsPolyampholyte()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("RKDERKDE");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.Polyampholyte),
            "FCR=1.0>0.35 with NCPR=0<=0.35 is the polyampholyte branch of get_disorder_class.");
    }

    // F2 — PPE: all R/K → f+ = 1.0, f- = 0, FCR = 1.0 > 0.35, NCPR = 1.0 > 0.35, f+ > 0.35 → Positive PE.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_AllPositiveCharge_ReturnsPositivePolyelectrolyte()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("RKRKRKRKRR");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.PositivePolyelectrolyte),
            "FCR>0.35, NCPR>0.35, f+ >0.35 selects the positive polyelectrolyte branch.");
    }

    // F3 — NPE: all D/E → f- = 1.0, f+ = 0, FCR = 1.0 > 0.35, NCPR = 1.0 > 0.35, f- > 0.35 → Negative PE.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_AllNegativeCharge_ReturnsNegativePolyelectrolyte()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("DEDEDEDEDD");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.NegativePolyelectrolyte),
            "FCR>0.35, NCPR>0.35, f- >0.35 selects the negative polyelectrolyte branch.");
    }

    // F4 — charge gate precedes composition: f+ = 0.4 (RK x2), FCR = 0.4 > 0.35, NCPR = 0.4 > 0.35,
    //      f+ > 0.35 → PPE even though P fraction (0.6) would otherwise be proline-rich.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_ChargedAndProlineRich_ChargeWinsAsPpe()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("RKRKPPPPPP");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.PositivePolyelectrolyte),
            "get_region_features tests the charge class first; PPE must override the enriched-P composition.");
    }

    // F5 — FCR at the strong-charge boundary: f+ = 0.35 exactly (FCR = 0.35, NOT > 0.35) → weakly charged,
    //      then no composition enriched → WeaklyCharged. (Strict '>' in get_disorder_class.)
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_FcrExactlyThreshold_IsNotStrongCharge()
    {
        // 7 of 20 R/K → f+ = 0.35; remaining 13 are A (no composition class).
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("RKRKRKRAAAAAAAAAAAAA");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.WeaklyCharged),
            "FCR=0.35 is NOT > 0.35, so the strong-charge gate fails and no composition class is enriched.");
    }

    // F5b — histidine counts as a positive charge (MobiDB-lite v3 states.py translation table
    //       intab='RKDEACFGHILMNPQSTVWY' / outab='PPNN____P___________' maps R,K,H → "P").
    //       8 H of 10 → f+ = 0.8 > 0.35, FCR = 0.8 > 0.35, NCPR = 0.8 > 0.35 → PPE.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_HistidineRich_ReturnsPositivePolyelectrolyte()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("HHHHHHHHAA");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.PositivePolyelectrolyte),
            "MobiDB-lite v3 maps H to the positive charge token; f+=0.8>0.35 → positive polyelectrolyte.");
    }

    // F5c — histidine contributes to the polyampholyte balance: 4 H (positive) + 4 D (negative)
    //       of 10 → f+ = 0.4, f- = 0.4, FCR = 0.8 > 0.35, NCPR = 0 <= 0.35 → Polyampholyte.
    //       (If H were ignored this would be NPE; the v3 source counts H as positive.)
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_HistidineBalancesNegative_ReturnsPolyampholyte()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("HHHHDDDDAA");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.Polyampholyte),
            "H (positive) + D (negative) balance to NCPR=0 with FCR=0.8 → polyampholyte per v3 states.py.");
    }

    #endregion

    #region ClassifyRegionFlavorMobiDbLite — composition classes (is_enriched threshold = 0.32)

    // F6 — Cysteine-rich: weakly charged, C fraction = 0.4 >= 0.32.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_CysteineEnriched_ReturnsCysteineRich()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("CCCCAAAAAA");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.CysteineRich),
            "C fraction 0.4 >= 0.32 with no strong charge yields the cysteine-rich class.");
    }

    // F7 — Proline-rich: weakly charged, P fraction = 0.4 >= 0.32.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_ProlineEnriched_ReturnsProlineRich()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("PPPPAAAAAA");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.ProlineRich),
            "P fraction 0.4 >= 0.32 yields the proline-rich class.");
    }

    // F8 — Glycine-rich: weakly charged, G fraction = 0.4 >= 0.32.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_GlycineEnriched_ReturnsGlycineRich()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("GGGGAAAAAA");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.GlycineRich),
            "G fraction 0.4 >= 0.32 yields the glycine-rich class.");
    }

    // F9 — Polar: weakly charged, {S,T,N,Q} fraction = 0.8 >= 0.32.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_PolarEnriched_ReturnsPolar()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("SSTTNNQQAA");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.Polar),
            "{S,T,N,Q} fraction 0.8 >= 0.32 yields the polar class (is_enriched(['S','T','N','Q'])).");
    }

    // F10 — composition priority C → P → G → polar: both C and P at 0.4, C must win.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_CysteineAndProlineEnriched_CysteineWinsByPriority()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("CCCCPPPPAA");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.CysteineRich),
            "consensus.py tests C before P; cysteine-rich must take priority over proline-rich.");
    }

    // F11 — enrichment threshold is inclusive (>= 0.32): 8 of 25 C = 0.32 exactly → CysteineRich.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_FractionExactlyThreshold_IsEnriched()
    {
        // 8 C + 17 A = length 25 → C fraction = 8/25 = 0.32 exactly.
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("CCCCCCCCAAAAAAAAAAAAAAAAA");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.CysteineRich),
            "is_enriched uses 's >= threshold'; a fraction of exactly 0.32 counts as enriched.");
    }

    // F12 — just below threshold: 7 of 25 C = 0.28 < 0.32 → not enriched → WeaklyCharged.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_FractionJustBelowThreshold_NotEnriched()
    {
        // 7 C + 18 A = length 25 → C fraction = 7/25 = 0.28 < 0.32.
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("CCCCCCCAAAAAAAAAAAAAAAAAA");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.WeaklyCharged),
            "C fraction 0.28 < 0.32 is not enriched; with no charge class the region is weakly charged.");
    }

    #endregion

    #region ClassifyRegionFlavorMobiDbLite — fallback and edge cases

    // F13 — no charge and no composition class enriched → WeaklyCharged (no MobiDB-lite subregion).
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_NoEnrichment_ReturnsWeaklyCharged()
    {
        var flavor = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("ALIVMFWALIVMFW");

        Assert.That(flavor, Is.EqualTo(DisorderPredictor.DisorderFlavor.WeaklyCharged),
            "A hydrophobic stretch has FCR=0 and no enriched composition class → weakly charged.");
    }

    // F14 — case-insensitive: lowercase input classifies identically to uppercase.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_LowercaseInput_SameAsUppercase()
    {
        var upper = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("RKDERKDE");
        var lower = DisorderPredictor.ClassifyRegionFlavorMobiDbLite("rkderkde");

        Assert.That(lower, Is.EqualTo(upper),
            "Input is upper-cased before classification, so case must not change the flavor.");
    }

    // F15 — null / empty input is rejected.
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_NullOrEmpty_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(
                () => DisorderPredictor.ClassifyRegionFlavorMobiDbLite(null!),
                "Null sequence must throw ArgumentException.");
            Assert.Throws<ArgumentException>(
                () => DisorderPredictor.ClassifyRegionFlavorMobiDbLite(""),
                "Empty sequence must throw ArgumentException.");
        });
    }

    // F16 — boundaries unchanged: the opt-in flavor labelling does not alter the validated
    //       TOP-IDP region boundaries reported by PredictDisorder (a poly-P region stays [0, L-1]).
    [Test]
    public void ClassifyRegionFlavorMobiDbLite_DoesNotAffectRegionBoundaries()
    {
        // 30 prolines → single disordered region spanning the whole sequence (validated TOP-IDP boundary).
        string polyP = new string('P', 30);
        var result = DisorderPredictor.PredictDisorder(polyP);

        Assert.Multiple(() =>
        {
            Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
                "TOP-IDP grouping yields exactly one region for a 30-residue poly-proline.");
            var region = result.DisorderedRegions[0];
            Assert.That(region.Start, Is.EqualTo(0), "Validated boundary Start is unaffected by flavor labelling.");
            Assert.That(region.End, Is.EqualTo(29), "Validated boundary End is unaffected by flavor labelling.");

            // The opt-in flavor of that exact region sequence is proline-rich (P fraction 1.0 >= 0.32).
            string regionSeq = polyP.Substring(region.Start, region.End - region.Start + 1);
            Assert.That(DisorderPredictor.ClassifyRegionFlavorMobiDbLite(regionSeq),
                Is.EqualTo(DisorderPredictor.DisorderFlavor.ProlineRich),
                "The MobiDB-lite flavor is computed from the region sequence; boundaries come from TOP-IDP.");
        });
    }

    #endregion
}
