// ONCO-DRIVER-001 — Driver Mutation Detection (20/20 rule)
// Evidence: docs/Evidence/ONCO-DRIVER-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-DRIVER-001.md
// Source: Vogelstein B et al. (2013). Science 339(6127):1546–1558. https://doi.org/10.1126/science.1235122
//         Tokheim C, Karchin R (2020). Bioinformatics 36(6):1712–1719. https://doi.org/10.1093/bioinformatics/btz759
//         Schroeder MP et al. (2014). Bioinformatics 30(17):i549–i555. https://doi.org/10.1093/bioinformatics/btu467
//         Miller ML et al. (2017). Oncotarget 8(20):33321–33333. https://doi.org/10.18632/oncotarget.15514

using GM = Seqeron.Genomics.Oncology.OncologyAnalyzer.GeneMutation;
using Cons = Seqeron.Genomics.Oncology.OncologyAnalyzer.MutationConsequence;
using Role = Seqeron.Genomics.Oncology.OncologyAnalyzer.DriverGeneRole;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_IdentifyDriverMutations_Tests
{
    // IDH1 archetype: 10 missense, all at codon 132 (Vogelstein 2013 / Miller 2017).
    private static GM[] Idh1Recurrent() =>
        Enumerable.Range(0, 10).Select(_ => new GM("IDH1", 132, Cons.Missense)).ToArray();

    // Tumor-suppressor archetype: 8 truncating at distinct positions + 2 missense at distinct positions.
    private static GM[] DispersedTsg()
    {
        var muts = new List<GM>();
        for (int i = 0; i < 5; i++) muts.Add(new GM("TP53", 100 + i, Cons.Nonsense));
        for (int i = 0; i < 2; i++) muts.Add(new GM("TP53", 200 + i, Cons.Frameshift));
        muts.Add(new GM("TP53", 300, Cons.SpliceSite));
        muts.Add(new GM("TP53", 400, Cons.Missense));
        muts.Add(new GM("TP53", 401, Cons.Missense));
        return muts.ToArray(); // 8 truncating, 2 distinct missense, N = 10
    }

    #region ClassifyGene Tests

    // M1 — IDH1: all missense at one codon → Oncogene; f_OG = 1.00, f_TSG = 0.00.
    // Evidence: Vogelstein 2013 (oncogene >20% recurrent missense); Miller 2017 (IDH1 codon 132).
    [Test]
    public void ClassifyGene_AllRecurrentMissense_ClassifiesOncogene()
    {
        var c = OncologyAnalyzer.ClassifyGene(Idh1Recurrent());

        Assert.Multiple(() =>
        {
            Assert.That(c.Role, Is.EqualTo(Role.Oncogene), "f_OG = 1.00 > 0.20 → oncogene");
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(1.00).Within(1e-10), "10/10 missense at codon 132 are recurrent");
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.00).Within(1e-10), "no truncating mutations");
            Assert.That(c.MutationCount, Is.EqualTo(10), "denominator N = 10");
        });
    }

    // M2 — Dispersed truncating → TumorSuppressor; f_TSG = 0.80.
    // Evidence: Vogelstein 2013 (TSG >20% truncating); Schroeder 2014 (nonsense/frameshift/splice).
    [Test]
    public void ClassifyGene_DispersedTruncating_ClassifiesTumorSuppressor()
    {
        var c = OncologyAnalyzer.ClassifyGene(DispersedTsg());

        Assert.Multiple(() =>
        {
            Assert.That(c.Role, Is.EqualTo(Role.TumorSuppressor), "f_TSG = 0.80 > 0.20 → tumor suppressor");
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.80).Within(1e-10), "8/10 truncating");
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(0.00).Within(1e-10), "2 missense at distinct positions are not recurrent");
        });
    }

    // M3 — Truncating fraction exactly 0.20 is NOT a TSG (strict > threshold).
    // Evidence: Vogelstein 2013 "more than 20%"; Tokheim 2020 ">20%".
    [Test]
    public void ClassifyGene_TruncatingFractionAtThreshold_IsNotTumorSuppressor()
    {
        var muts = new List<GM> { new("BRCA1", 10, Cons.Nonsense), new("BRCA1", 11, Cons.Frameshift) };
        for (int i = 0; i < 8; i++) muts.Add(new GM("BRCA1", 100 + i, Cons.Missense)); // distinct → not recurrent
        // f_TSG = 2/10 = 0.20 (not > 0.20); f_OG = 0.00

        var c = OncologyAnalyzer.ClassifyGene(muts);

        Assert.Multiple(() =>
        {
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.20).Within(1e-10), "2/10 truncating");
            Assert.That(c.Role, Is.Not.EqualTo(Role.TumorSuppressor), "exactly 0.20 does not exceed the strict threshold");
            Assert.That(c.Role, Is.EqualTo(Role.Ambiguous), "neither criterion exceeds 0.20");
        });
    }

    // M4 — Truncating fraction 0.30 (just above) → TumorSuppressor.
    // Evidence: Vogelstein 2013.
    [Test]
    public void ClassifyGene_TruncatingFractionAboveThreshold_ClassifiesTumorSuppressor()
    {
        var muts = new List<GM>();
        for (int i = 0; i < 3; i++) muts.Add(new GM("APC", 10 + i, Cons.Nonsense));
        for (int i = 0; i < 7; i++) muts.Add(new GM("APC", 100 + i, Cons.Missense)); // distinct
        // f_TSG = 3/10 = 0.30

        var c = OncologyAnalyzer.ClassifyGene(muts);

        Assert.Multiple(() =>
        {
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.30).Within(1e-10), "3/10 truncating");
            Assert.That(c.Role, Is.EqualTo(Role.TumorSuppressor), "0.30 > 0.20 → tumor suppressor");
        });
    }

    // S2 — Singleton missense positions are not recurrent.
    // Evidence: Miller 2017 — recurrence requires ≥ 2 at the same position.
    [Test]
    public void ClassifyGene_AllMissenseDistinctPositions_NotRecurrentNotOncogene()
    {
        var muts = new[]
        {
            new GM("EGFR", 700, Cons.Missense),
            new GM("EGFR", 800, Cons.Missense),
            new GM("EGFR", 900, Cons.Missense),
        };

        var c = OncologyAnalyzer.ClassifyGene(muts);

        Assert.Multiple(() =>
        {
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(0.00).Within(1e-10), "no position has ≥ 2 missense");
            Assert.That(c.Role, Is.EqualTo(Role.Ambiguous), "no recurrent missense and no truncating → ambiguous");
        });
    }

    // S3 — Recurrent-missense fraction exactly 0.20 is NOT an oncogene (strict >).
    // Evidence: Miller 2017 (≥2 = recurrent) + Vogelstein ">20%".
    [Test]
    public void ClassifyGene_RecurrentMissenseFractionAtThreshold_IsNotOncogene()
    {
        var muts = new List<GM>
        {
            new("KRAS", 12, Cons.Missense), // codon 12 seen twice → recurrent
            new("KRAS", 12, Cons.Missense),
        };
        for (int i = 0; i < 8; i++) muts.Add(new GM("KRAS", 200 + i, Cons.Missense)); // distinct, not recurrent
        // recurrent missense = 2/10 = 0.20 (not > 0.20)

        var c = OncologyAnalyzer.ClassifyGene(muts);

        Assert.Multiple(() =>
        {
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(0.20).Within(1e-10), "2/10 missense at codon 12");
            Assert.That(c.Role, Is.Not.EqualTo(Role.Oncogene), "exactly 0.20 does not exceed the strict threshold");
            Assert.That(c.Role, Is.EqualTo(Role.Ambiguous), "neither criterion exceeds 0.20 → ambiguous");
        });
    }

    // S1 — Low-recurrence gene → Ambiguous (neither criterion > 0.20).
    // Evidence: OncodriveROLE 2014 — lowly recurrent genes left unclassified.
    [Test]
    public void ClassifyGene_LowRecurrence_ClassifiesAmbiguous()
    {
        var muts = new List<GM> { new("XYZ1", 10, Cons.Nonsense) }; // 1 truncating → 0.20
        for (int i = 0; i < 4; i++) muts.Add(new GM("XYZ1", 50 + i, Cons.Missense)); // distinct
        // f_TSG = 1/5 = 0.20 (not >); f_OG = 0.00

        var c = OncologyAnalyzer.ClassifyGene(muts);

        Assert.That(c.Role, Is.EqualTo(Role.Ambiguous), "neither fraction exceeds 0.20 → ambiguous");
    }

    // Dual-pass tie-break — both criteria exceed 0.20; dominant fraction wins (Assumption #1; INV-02/03).
    // Construction: pos X carries 3 missense (recurrent) and pos Y 3 missense (recurrent) → f_OG = 6/10 = 0.60;
    // 3 nonsense → f_TSG = 3/10 = 0.30; 1 Other pads the denominator. Both > 0.20, f_OG dominant → Oncogene.
    // Evidence: Vogelstein 2013 well-documented genes "far surpass" one criterion; dominant-signal resolution.
    [Test]
    public void ClassifyGene_DualPassOncogeneDominant_ClassifiesOncogene()
    {
        var muts = new List<GM>
        {
            new("MIX", 10, Cons.Missense), new("MIX", 10, Cons.Missense), new("MIX", 10, Cons.Missense),
            new("MIX", 20, Cons.Missense), new("MIX", 20, Cons.Missense), new("MIX", 20, Cons.Missense),
            new("MIX", 30, Cons.Nonsense), new("MIX", 31, Cons.Nonsense), new("MIX", 32, Cons.Nonsense),
            new("MIX", 40, Cons.Other),
        };

        var c = OncologyAnalyzer.ClassifyGene(muts);

        Assert.Multiple(() =>
        {
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(0.60).Within(1e-10), "6/10 recurrent missense");
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.30).Within(1e-10), "3/10 truncating");
            Assert.That(c.Role, Is.EqualTo(Role.Oncogene), "both > 0.20; f_OG (0.60) dominates f_TSG (0.30)");
        });
    }

    // Dual-pass exact tie — both fractions equal and > 0.20 → Ambiguous (Assumption #1: exact tie is ambiguous).
    // Construction: pos X carries 3 missense (recurrent) → f_OG = 3/10 = 0.30; 3 nonsense → f_TSG = 3/10 = 0.30;
    // 4 Other pad the denominator. f_OG == f_TSG, both > 0.20 → genuinely ambiguous.
    [Test]
    public void ClassifyGene_DualPassExactTie_ClassifiesAmbiguous()
    {
        var muts = new List<GM>
        {
            new("TIE", 10, Cons.Missense), new("TIE", 10, Cons.Missense), new("TIE", 10, Cons.Missense),
            new("TIE", 30, Cons.Nonsense), new("TIE", 31, Cons.Nonsense), new("TIE", 32, Cons.Nonsense),
            new("TIE", 40, Cons.Other), new("TIE", 41, Cons.Other), new("TIE", 42, Cons.Other), new("TIE", 43, Cons.Other),
        };

        var c = OncologyAnalyzer.ClassifyGene(muts);

        Assert.Multiple(() =>
        {
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(0.30).Within(1e-10), "3/10 recurrent missense");
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.30).Within(1e-10), "3/10 truncating");
            Assert.That(c.Role, Is.EqualTo(Role.Ambiguous), "exact tie of two passing criteria → ambiguous");
        });
    }

    // Other-consequence counts toward the denominator only (neither truncating nor missense).
    // Evidence: Vogelstein 2013 — fractions are over ALL recorded mutations in the gene; synonymous/other
    // changes dilute both fractions. Here 2 nonsense among N=10 → f_TSG = 0.20 (not > 0.20) → not a TSG.
    [Test]
    public void ClassifyGene_OtherConsequence_CountsInDenominatorOnly()
    {
        var muts = new List<GM>
        {
            new("DEN", 10, Cons.Nonsense), new("DEN", 11, Cons.Nonsense),
            new("DEN", 20, Cons.Missense), new("DEN", 21, Cons.Missense), // distinct → not recurrent
        };
        for (int i = 0; i < 6; i++) muts.Add(new GM("DEN", 50 + i, Cons.Other));
        // N = 10; truncating = 2 → f_TSG = 0.20 (NOT > 0.20); recurrent missense = 0

        var c = OncologyAnalyzer.ClassifyGene(muts);

        Assert.Multiple(() =>
        {
            Assert.That(c.MutationCount, Is.EqualTo(10), "Other mutations are part of the denominator");
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.20).Within(1e-10), "2 truncating / 10 total");
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(0.00).Within(1e-10), "no recurrent missense");
            Assert.That(c.Role, Is.EqualTo(Role.Ambiguous), "0.20 does not exceed the strict threshold");
        });
    }

    // C1 — Empty gene mutations → Ambiguous, fractions 0.
    [Test]
    public void ClassifyGene_EmptyMutations_ReturnsAmbiguousZeroFractions()
    {
        var c = OncologyAnalyzer.ClassifyGene(System.Array.Empty<GM>());

        Assert.Multiple(() =>
        {
            Assert.That(c.Role, Is.EqualTo(Role.Ambiguous), "no evidence → ambiguous");
            Assert.That(c.MutationCount, Is.EqualTo(0), "denominator 0");
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.0).Within(1e-10), "no truncating fraction");
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(0.0).Within(1e-10), "no recurrent missense fraction");
        });
    }

    // C2 — Null input rejected.
    [Test]
    public void ClassifyGene_NullInput_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => OncologyAnalyzer.ClassifyGene(null!), "Null mutation collection is rejected");
    }

    #endregion

    #region ScoreDriverPotential Tests

    // M8 — Score equals the recurrent-missense fraction for the IDH1 set.
    // Evidence: Vogelstein 2013 20/20 driver-signal = max(f_TSG, f_OG).
    [Test]
    public void ScoreDriverPotential_RecurrentMissenseGene_ReturnsRecurrentFraction()
    {
        double score = OncologyAnalyzer.ScoreDriverPotential(Idh1Recurrent());

        Assert.That(score, Is.EqualTo(1.00).Within(1e-10), "max(0.00, 1.00) = 1.00");
    }

    // M9 — Score equals the truncating fraction for the dispersed-TSG set.
    [Test]
    public void ScoreDriverPotential_TruncatingGene_ReturnsTruncatingFraction()
    {
        double score = OncologyAnalyzer.ScoreDriverPotential(DispersedTsg());

        Assert.That(score, Is.EqualTo(0.80).Within(1e-10), "max(0.80, 0.00) = 0.80");
    }

    [Test]
    public void ScoreDriverPotential_NullInput_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => OncologyAnalyzer.ScoreDriverPotential(null!), "Null mutation collection is rejected");
    }

    #endregion

    #region MatchCancerHotspots Tests

    // M5 — A mutation at a known (gene, position) hotspot matches.
    // Evidence: Miller 2017 recurrent-position / caller-supplied catalog.
    [Test]
    public void MatchCancerHotspots_KnownHotspot_ReturnsTrue()
    {
        var hotspots = new HashSet<(string, int)> { ("BRAF", 600), ("KRAS", 12) };
        bool hit = OncologyAnalyzer.MatchCancerHotspots(new GM("BRAF", 600, Cons.Missense), hotspots);

        Assert.That(hit, Is.True, "(BRAF, 600) is in the supplied hotspot set");
    }

    // M6 — A mutation not in the set does not match.
    [Test]
    public void MatchCancerHotspots_NotInSet_ReturnsFalse()
    {
        var hotspots = new HashSet<(string, int)> { ("BRAF", 600) };

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.MatchCancerHotspots(new GM("BRAF", 601, Cons.Missense), hotspots),
                Is.False, "wrong position is not a hotspot");
            Assert.That(OncologyAnalyzer.MatchCancerHotspots(new GM("EGFR", 600, Cons.Missense), hotspots),
                Is.False, "wrong gene at same position is not a hotspot");
        });
    }

    [Test]
    public void MatchCancerHotspots_NullSet_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => OncologyAnalyzer.MatchCancerHotspots(new GM("BRAF", 600, Cons.Missense), null!),
            "Null hotspot set is rejected");
    }

    #endregion

    #region IdentifyDriverMutations Tests

    // M7 — Returns only mutations in driver genes; result ⊆ input (INV-01).
    // Evidence: Vogelstein 2013 per-gene 20/20 classification.
    [Test]
    public void IdentifyDriverMutations_MixedPanel_ReturnsOnlyDriverGeneMutations()
    {
        var input = new List<GM>();
        input.AddRange(Idh1Recurrent());       // IDH1 → Oncogene (driver)
        input.AddRange(DispersedTsg());        // TP53 → TumorSuppressor (driver)
        // PASS gene: 5 missense at distinct positions → Ambiguous (non-driver)
        for (int i = 0; i < 5; i++) input.Add(new GM("PASS", 10 + i, Cons.Missense));

        var drivers = OncologyAnalyzer.IdentifyDriverMutations(input);

        Assert.Multiple(() =>
        {
            Assert.That(drivers, Has.Count.EqualTo(20), "10 IDH1 + 10 TP53 driver mutations; PASS excluded");
            Assert.That(drivers.All(m => m.Gene is "IDH1" or "TP53"), Is.True, "only driver-gene mutations retained");
            Assert.That(drivers.Any(m => m.Gene == "PASS"), Is.False, "ambiguous gene mutations excluded");
            Assert.That(drivers.All(input.Contains), Is.True, "INV-01: result ⊆ input");
        });
    }

    // M5/M6 integration — a hotspot rescues a mutation in an otherwise non-driver gene.
    // Evidence: Miller 2017 — known hotspot is a driver signal independent of the 20/20 gene call.
    [Test]
    public void IdentifyDriverMutations_HotspotInNonDriverGene_IsRetained()
    {
        var input = new List<GM>();
        for (int i = 0; i < 5; i++) input.Add(new GM("PASS", 10 + i, Cons.Missense)); // Ambiguous gene
        input.Add(new GM("PASS", 600, Cons.Missense)); // but this position is a known hotspot
        var hotspots = new HashSet<(string, int)> { ("PASS", 600) };

        var drivers = OncologyAnalyzer.IdentifyDriverMutations(input, hotspots);

        Assert.Multiple(() =>
        {
            Assert.That(drivers, Has.Count.EqualTo(1), "only the hotspot mutation is a driver");
            Assert.That(drivers[0].ProteinPosition, Is.EqualTo(600), "the retained mutation is the hotspot");
        });
    }

    // C3 — Property: every returned mutation is present in the input (INV-01), and order is preserved.
    [Test]
    public void IdentifyDriverMutations_Property_ResultIsSubsetInInputOrder()
    {
        var input = new List<GM>();
        input.AddRange(DispersedTsg()); // driver
        for (int i = 0; i < 4; i++) input.Add(new GM("Z", 1 + i, Cons.Missense)); // ambiguous

        var drivers = OncologyAnalyzer.IdentifyDriverMutations(input);

        // Verify subset + order: driver list is the input filtered, so it must appear as an in-order subsequence.
        int idx = -1;
        bool inOrder = true;
        foreach (var d in drivers)
        {
            idx = input.IndexOf(d, idx + 1);
            if (idx < 0) { inOrder = false; break; }
        }
        Assert.Multiple(() =>
        {
            Assert.That(drivers.All(input.Contains), Is.True, "INV-01: every driver is an input mutation");
            Assert.That(inOrder, Is.True, "drivers preserve input order");
        });
    }

    [Test]
    public void IdentifyDriverMutations_EmptyInput_ReturnsEmpty()
    {
        var drivers = OncologyAnalyzer.IdentifyDriverMutations(System.Array.Empty<GM>());

        Assert.That(drivers, Is.Empty, "no mutations in, no drivers out");
    }

    [Test]
    public void IdentifyDriverMutations_NullInput_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => OncologyAnalyzer.IdentifyDriverMutations(null!), "Null mutation collection is rejected");
    }

    #endregion
}
