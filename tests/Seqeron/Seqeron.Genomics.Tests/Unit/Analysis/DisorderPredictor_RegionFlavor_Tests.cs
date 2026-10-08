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

namespace Seqeron.Genomics.Tests.Unit.Analysis;

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

    #region PredictFlavorSubregionsMobiDbLite — full windowed MobiDB-lite v3 feature step

    // Expected values below were produced by running the VERBATIM MobiDB-lite v3 code
    // (raw.githubusercontent.com/BioComputingUP/MobiDB-lite/v3/mdblib/{states,consensus,prediction}.py,
    // MobidbLiteConsensus.get_region_features(merge=True, only_in_idr=True)) with the given IDRs
    // and SEG mask; coordinates converted from its 1-based output to 0-based inclusive.
    // Additionally, 400 random fixtures (seed 20260928) matched the reference 400/400 during review.

    private static readonly bool[] NoLc58 = new bool[58];

    private const string AlphaSynuclein =
        "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA";

    // FS1 — two feature blocks separated by neutral residues: PPE and proline-rich sub-regions.
    [Test]
    public void PredictFlavorSubregions_TwoFeatureBlocks_MatchesMobiDbLiteReference()
    {
        string seq = new string('A', 10) + string.Concat(Enumerable.Repeat("RK", 7)) + new string('A', 11)
                     + new string('P', 13) + new string('A', 10);
        var got = DisorderPredictor.PredictFlavorSubregionsMobiDbLite(seq, new[] { (0, 57) }, NoLc58);
        Assert.That(got.Select(r => (r.Start, r.End, r.Flavor)), Is.EqualTo(new[]
        {
            (9, 24, DisorderPredictor.DisorderFlavor.PositivePolyelectrolyte),
            (34, 48, DisorderPredictor.DisorderFlavor.ProlineRich)
        }));
    }

    // FS2 — α-synuclein (P37840) with its TOP-IDP regions (10–43, 47–66, 94–139): the acidic C-terminal
    // tail splits into PA (94–104) and NPE (111–139); the region-level label of 94–139 is PA
    // (f+ = 3/46, f− = 15/46, FCR = 0.391 > 0.35, NCPR = 0.261 ≤ 0.35) — the windowed algorithm
    // resolves the NPE tail that the whole-region label cannot.
    [Test]
    public void PredictFlavorSubregions_AlphaSynuclein_MatchesMobiDbLiteReference()
    {
        var regions = DisorderPredictor.PredictDisorderRegions(AlphaSynuclein).DisorderedRegions
            .Select(r => (r.Start, r.End)).ToList();
        var got = DisorderPredictor.PredictFlavorSubregionsMobiDbLite(
            AlphaSynuclein, regions, new bool[AlphaSynuclein.Length]);
        Assert.Multiple(() =>
        {
            Assert.That(regions, Is.EqualTo(new[] { (10, 43), (47, 66), (94, 139) }));
            Assert.That(got.Select(r => (r.Start, r.End, r.Flavor)), Is.EqualTo(new[]
            {
                (94, 104, DisorderPredictor.DisorderFlavor.Polyampholyte),
                (111, 139, DisorderPredictor.DisorderFlavor.NegativePolyelectrolyte)
            }));
            Assert.That(DisorderPredictor.ClassifyRegionFlavorMobiDbLite(AlphaSynuclein[94..140]),
                Is.EqualTo(DisorderPredictor.DisorderFlavor.Polyampholyte));
        });
    }

    // FS3 — the SEG low-complexity class (code 7) outranks polar (code 8); without a mask the
    // same (SQ)10 block is polar, and morphology/windowing widens it to 9–30.
    [Test]
    public void PredictFlavorSubregions_LowComplexityBeatsPolar_MatchesMobiDbLiteReference()
    {
        string seq = new string('W', 10) + string.Concat(Enumerable.Repeat("SQ", 10)) + new string('W', 10);
        var lc = Enumerable.Range(0, 40).Select(i => i >= 10 && i < 30).ToArray();
        var withLc = DisorderPredictor.PredictFlavorSubregionsMobiDbLite(seq, new[] { (0, 39) }, lc);
        var noLc = DisorderPredictor.PredictFlavorSubregionsMobiDbLite(seq, new[] { (0, 39) }, new bool[40]);
        Assert.Multiple(() =>
        {
            Assert.That(withLc.Select(r => (r.Start, r.End, r.Flavor)),
                Is.EqualTo(new[] { (10, 29, DisorderPredictor.DisorderFlavor.LowComplexity) }));
            Assert.That(noLc.Select(r => (r.Start, r.End, r.Flavor)),
                Is.EqualTo(new[] { (9, 30, DisorderPredictor.DisorderFlavor.Polar) }));
        });
    }

    // FS4 — math_morphology(rmax=5): a 3-residue neutral gap inside NPE is closed (one run 0–26);
    // an 8-residue gap is not (two runs).
    [Test]
    public void PredictFlavorSubregions_Morphology_ClosesShortGapsOnly()
    {
        string shortGap = new string('E', 12) + "AAA" + new string('E', 12);
        string longGap = new string('E', 12) + new string('A', 8) + new string('E', 12);
        var a = DisorderPredictor.PredictFlavorSubregionsMobiDbLite(shortGap, null, new bool[shortGap.Length]);
        var b = DisorderPredictor.PredictFlavorSubregionsMobiDbLite(longGap, null, new bool[longGap.Length]);
        Assert.Multiple(() =>
        {
            Assert.That(a.Select(r => (r.Start, r.End, r.Flavor)),
                Is.EqualTo(new[] { (0, 26, DisorderPredictor.DisorderFlavor.NegativePolyelectrolyte) }));
            Assert.That(b.Select(r => (r.Start, r.End, r.Flavor)), Is.EqualTo(new[]
            {
                (0, 12, DisorderPredictor.DisorderFlavor.NegativePolyelectrolyte),
                (19, 31, DisorderPredictor.DisorderFlavor.NegativePolyelectrolyte)
            }));
        });
    }

    // FS5 — feature_len_thr = 10 is applied to the run clipped to the IDR: 9 residues → none, 10 → one.
    [Test]
    public void PredictFlavorSubregions_MinimumLengthTen_AppliedWithinRegion()
    {
        string polyE = new string('E', 30);
        var nine = DisorderPredictor.PredictFlavorSubregionsMobiDbLite(polyE, new[] { (5, 13) }, new bool[30]);
        var ten = DisorderPredictor.PredictFlavorSubregionsMobiDbLite(polyE, new[] { (5, 14) }, new bool[30]);
        Assert.Multiple(() =>
        {
            Assert.That(nine, Is.Empty);
            Assert.That(ten.Select(r => (r.Start, r.End, r.Flavor)),
                Is.EqualTo(new[] { (5, 14, DisorderPredictor.DisorderFlavor.NegativePolyelectrolyte) }));
        });
    }

    // FS6 — per-residue windows separate a glycine block from an acidic block (G12 D12), whereas the
    // whole-region label is NPE (f− = 0.5).
    [Test]
    public void PredictFlavorSubregions_GlycineThenAcidic_SplitsIntoTwoFlavors()
    {
        string seq = new string('G', 12) + new string('D', 12);
        var got = DisorderPredictor.PredictFlavorSubregionsMobiDbLite(seq, null, new bool[24]);
        Assert.Multiple(() =>
        {
            Assert.That(got.Select(r => (r.Start, r.End, r.Flavor)), Is.EqualTo(new[]
            {
                (0, 10, DisorderPredictor.DisorderFlavor.GlycineRich),
                (11, 23, DisorderPredictor.DisorderFlavor.NegativePolyelectrolyte)
            }));
            Assert.That(DisorderPredictor.ClassifyRegionFlavorMobiDbLite(seq),
                Is.EqualTo(DisorderPredictor.DisorderFlavor.NegativePolyelectrolyte));
        });
    }

    // FS7 — short / edge inputs: L = 2 (window shrinks to L, no run ≥ 10); lowercase input.
    [Test]
    public void PredictFlavorSubregions_ShortAndLowercaseInputs()
    {
        Assert.Multiple(() =>
        {
            Assert.That(DisorderPredictor.PredictFlavorSubregionsMobiDbLite("RK", null, new bool[2]), Is.Empty);
            Assert.That(DisorderPredictor.PredictFlavorSubregionsMobiDbLite(new string('e', 12), null, new bool[12])
                    .Select(r => (r.Start, r.End, r.Flavor)),
                Is.EqualTo(new[] { (0, 11, DisorderPredictor.DisorderFlavor.NegativePolyelectrolyte) }));
        });
    }

    // FS8 — input validation.
    [Test]
    public void PredictFlavorSubregions_InvalidInputs_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => DisorderPredictor.PredictFlavorSubregionsMobiDbLite(""));
            Assert.Throws<ArgumentException>(() => DisorderPredictor.PredictFlavorSubregionsMobiDbLite(null!));
            Assert.Throws<ArgumentException>(() =>
                DisorderPredictor.PredictFlavorSubregionsMobiDbLite("EEEE", null, new bool[3]));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                DisorderPredictor.PredictFlavorSubregionsMobiDbLite("EEEE", new[] { (2, 4) }, new bool[4]));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                DisorderPredictor.PredictFlavorSubregionsMobiDbLite("EEEE", new[] { (3, 2) }, new bool[4]));
        });
    }

    #endregion
}
