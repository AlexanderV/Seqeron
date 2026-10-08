using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// MIRNA-* mutation killers: exact-value tests that pin the documented formulas of the
/// public helper methods whose canonical tests only used range assertions
/// (<c>GreaterThan</c>/<c>GreaterThanOrEqualTo</c>), leaving boundary, arithmetic and
/// logical mutants alive under Stryker. Each assertion below reproduces the published
/// rule exactly so an injected operator change diverges from the asserted value.
///
/// Evidence:
///  - Target-site context (3'UTR AU enrichment, positional bias): Grimson et al. (2007)
///    Mol Cell 27:91–105; Agarwal et al. (2015) eLife 4:e05005 (TargetScan context++).
///  - Local accessibility: McCaskill region-unpaired probability (RNAplfold-style local
///    window; Bernhart et al. 2006; TargetScan SA), via RnaSecondaryStructure.
///  - Seed families (shared 7-nt seed, single-mismatch neighbours): Bartel (2009) Cell
///    136:215–233; Lewis et al. (2005) Cell 120:15–20.
/// </summary>
[TestFixture]
public class MiRnaAnalyzerMutationTests
{
    private const double Tol = 1e-9;

    #region AnalyzeTargetContext — Grimson (2007) site context

    // Model (MiRnaAnalyzer.AnalyzeTargetContext XML doc; Grimson et al. 2007 Mol Cell 27:91;
    // TargetScan targetscan_70_context_scores.pl getLocalAU_contribution / $MIN_DIST_TO_CDS = 15):
    //   AuContent  = Σ_{A/U flank nt} 1/d ÷ Σ_{flank nt} 1/d, d = distance from the site edge
    //                (up to W = 30 nt each side, site excluded)
    //   NearStart  = 1-based start < 15 (ribosome-occluded, TargetScan "too_close")
    //   NearEnd    = end > len*0.85 (descriptive)
    //   ContextScore = NearStart ? 0 : 0.5*AuContent + 0.5*(1 − min(d5,d3)/((d5+d3)/2))

    private static readonly string PolyA100 = new string('A', 100);

    [Test]
    public void AnalyzeTargetContext_MiddleSite_ScoresLowerThanEndSite()
    {
        // Grimson 2007: sites in the middle of the UTR are the LEAST effective; sites near the
        // ends are more effective. (40,50): d5=40, d3=49 ⇒ EndProximity = 1 − 40/44.5.
        var mid = AnalyzeTargetContext(PolyA100, 40, 50);
        // (90,95): d5=90, d3=4 ⇒ EndProximity = 1 − 4/47.
        var end = AnalyzeTargetContext(PolyA100, 90, 95);

        Assert.Multiple(() =>
        {
            Assert.That(mid.AuContent, Is.EqualTo(1.0).Within(Tol));
            Assert.That(mid.NearStart, Is.False);
            Assert.That(mid.NearEnd, Is.False);
            Assert.That(mid.ContextScore, Is.EqualTo(0.5 + 0.5 * (1 - 40 / 44.5)).Within(Tol)); // 0.5505617977528090
            Assert.That(end.NearEnd, Is.True);
            Assert.That(end.ContextScore, Is.EqualTo(0.5 + 0.5 * (1 - 4.0 / 47.0)).Within(Tol)); // 0.9574468085106382
            Assert.That(end.ContextScore, Is.GreaterThan(mid.ContextScore));
        });
    }

    [Test]
    public void AnalyzeTargetContext_WithinFirst15nt_NearStart_ContextScoreZero()
    {
        // Grimson 2007: sites within ~15 nt of the stop codon are cleared by the ribosome.
        // TargetScan: `if ($utrStart < $MIN_DIST_TO_CDS)` (1-based, MIN = 15) ⇒ no context score.
        var ctx = AnalyzeTargetContext(PolyA100, 5, 10);
        Assert.That(ctx.NearStart, Is.True);
        Assert.That(ctx.ContextScore, Is.EqualTo(0.0).Within(Tol));
    }

    [Test]
    public void AnalyzeTargetContext_Start15ntBoundary()
    {
        // 0-based 13 ⇒ 1-based 14 < 15 ⇒ too close; 0-based 14 ⇒ 1-based 15 ⇒ scored.
        Assert.That(AnalyzeTargetContext(PolyA100, 13, 20).NearStart, Is.True);
        var ok = AnalyzeTargetContext(PolyA100, 14, 20);
        Assert.That(ok.NearStart, Is.False);
        Assert.That(ok.ContextScore, Is.GreaterThan(0.0));
    }

    [Test]
    public void AnalyzeTargetContext_EndExactlyAtThreshold_IsNotNearEnd()
    {
        // end = 85 == len*0.85. Strict '>' ⇒ NOT nearEnd.
        var ctx = AnalyzeTargetContext(PolyA100, 50, 85);
        Assert.That(ctx.NearEnd, Is.False);
    }

    [Test]
    public void AnalyzeTargetContext_NoAu_OnlyPositionTerm()
    {
        string polyGc = string.Concat(Enumerable.Repeat("GC", 50)); // 100 nt, no A/U
        var ctx = AnalyzeTargetContext(polyGc, 40, 50);

        Assert.That(ctx.AuContent, Is.EqualTo(0.0).Within(Tol));
        Assert.That(ctx.ContextScore, Is.EqualTo(0.5 * (1 - 40 / 44.5)).Within(Tol)); // 0.050561797752809
    }

    [Test]
    public void AnalyzeTargetContext_AuContent_IsInverseDistanceWeighted_SiteExcluded()
    {
        // TargetScan getLocalAU_contribution: weight 1/(i+1) for the i-th flank nt outward.
        // Site [10..15] = CCCCCC; upstream (outward) G, A×9; downstream G, U×9.
        // Each flank: Σ weights = H10, A/U weight = H10 − 1 ⇒ fraction = (H10 − 1)/H10
        // = 0.6585828478525945 (an unweighted count would give 18/20 = 0.9; counting the
        // site itself, as the old code did, gives even less).
        string mrna = "AAAAAAAAAG" + "CCCCCC" + "GUUUUUUUUU";
        double h10 = Enumerable.Range(1, 10).Sum(k => 1.0 / k);
        var ctx = AnalyzeTargetContext(mrna, 10, 15);
        Assert.That(ctx.AuContent, Is.EqualTo((h10 - 1) / h10).Within(Tol));

        // contextWindow = 1 ⇒ only the adjacent G on each side ⇒ 0.
        Assert.That(AnalyzeTargetContext(mrna, 10, 15, contextWindow: 1).AuContent, Is.EqualTo(0.0).Within(Tol));
    }

    [Test]
    public void AnalyzeTargetContext_InvalidCoordinates_ReturnZeros()
    {
        Assert.That(AnalyzeTargetContext(PolyA100, -1, 5), Is.EqualTo((0.0, false, false, 0.0)));
        Assert.That(AnalyzeTargetContext(PolyA100, 10, 100), Is.EqualTo((0.0, false, false, 0.0)));
        Assert.That(AnalyzeTargetContext(PolyA100, 20, 10), Is.EqualTo((0.0, false, false, 0.0)));
    }

    #endregion

    #region CalculateSiteAccessibility — McCaskill region-unpaired probability (delegates to RnaSecondaryStructure)

    // Accessibility = P(site entirely unpaired) = Z_open/Z (Turner 2004 McCaskill), via the canonical
    // RnaSecondaryStructure.CalculateRegionUnpairedProbability over an RNAplfold-style W = 80 local
    // context. For sequences ≤ 80 nt the context is the whole sequence.
    // ViennaRNA 2.x cross-check (RNA.fold_compound pf with hc_add_up over the site; Z_c/Z):
    //   seq                                              site    Vienna d0   Vienna d2   Seqeron
    //   GGGGGCUACCUCAGGGGG                               5..12   1.872e-3    1.28e-4     3.225e-4
    //   GGGAAACCCAAAGGGUUUCCCAAGCUACCUCAAA               23..30  0.993777    0.976954    0.972098
    //   AUGCUACCUCAAAAAAAAAAAAAAAAA                      3..10   0.993008    0.970119    0.963046
    // The superseded pair-density heuristic returned 0 for all three.

    [TestCase("GGGAAACCCAAAGGGUUUCCCAAGCUACCUCAAA", 23, 30, 0.976954)]
    [TestCase("AUGCUACCUCAAAAAAAAAAAAAAAAA", 3, 10, 0.970119)]
    public void CalculateSiteAccessibility_UnstructuredSite_MatchesViennaRna(string seq, int s, int e, double vienna)
    {
        double acc = CalculateSiteAccessibility(seq, s, e);
        Assert.That(acc, Is.EqualTo(vienna).Within(0.01));
    }

    [Test]
    public void CalculateSiteAccessibility_SiteLockedInHelix_NearZero_MatchesViennaRna()
    {
        // GGGGG·CUACCUCA·GGGGG: the site's C's pair with the G flanks. Vienna: 1.9e-3 (d0) / 1.3e-4 (d2).
        double acc = CalculateSiteAccessibility("GGGGGCUACCUCAGGGGG", 5, 12);
        Assert.That(acc, Is.LessThan(0.002));
    }

    [Test]
    public void CalculateSiteAccessibility_DelegatesToCanonicalRegionUnpairedProbability()
    {
        const string seq = "GGGAAACCCAAAGGGUUUCCCAAGCUACCUCAAA"; // 34 nt ⇒ context = whole sequence
        double canonical = Seqeron.Genomics.Analysis.RnaSecondaryStructure
            .CalculateRegionUnpairedProbability(seq, windowEnd: 30, windowLength: 8);
        Assert.That(CalculateSiteAccessibility(seq, 23, 30), Is.EqualTo(canonical));
        // DNA / lower-case input is the same molecule.
        Assert.That(CalculateSiteAccessibility(seq.ToLowerInvariant().Replace('u', 't'), 23, 30), Is.EqualTo(canonical));
    }

    [Test]
    public void CalculateSiteAccessibility_LongSequence_UsesLocal80ntContext()
    {
        // 200-nt sequence: the site [100..107] is folded within the 80-nt window centred on it
        // (RNAplfold -W 80), i.e. context [64..143], local site end = 107 − 64 = 43.
        var rnd = new Random(7);
        string seq = new string(Enumerable.Range(0, 200).Select(_ => "ACGU"[rnd.Next(4)]).ToArray());
        double expected = Seqeron.Genomics.Analysis.RnaSecondaryStructure
            .CalculateRegionUnpairedProbability(seq.Substring(64, 80), windowEnd: 43, windowLength: 8);
        Assert.That(CalculateSiteAccessibility(seq, 100, 107), Is.EqualTo(expected));
    }

    [Test]
    public void CalculateSiteAccessibility_PolyA_NoPairs_FullyAccessible()
    {
        Assert.That(CalculateSiteAccessibility(new string('A', 120), 60, 65), Is.EqualTo(1.0).Within(Tol));
    }

    [Test]
    public void CalculateSiteAccessibility_InvalidCoordinates_ReturnZero()
    {
        const string seq = "GAAAAUAAAC";
        Assert.Multiple(() =>
        {
            Assert.That(CalculateSiteAccessibility(seq, 2, seq.Length), Is.EqualTo(0.0));
            Assert.That(CalculateSiteAccessibility(seq, -1, 5), Is.EqualTo(0.0));
            Assert.That(CalculateSiteAccessibility(seq, 6, 2), Is.EqualTo(0.0));
            Assert.That(CalculateSiteAccessibility("", 0, 0), Is.EqualTo(0.0));
        });
    }

    #endregion

    #region FindSimilarMiRnas — seed-family Hamming neighbours

    // Seeds are miRNA positions 2-8 (7 nt). Family membership = seed Hamming distance ≤ maxMismatches.
    //   query seed = GGGGGGG
    //   m1   seed  = GGGGGCC  (2 mismatches)
    //   m2   seed  = GGGGCCC  (3 mismatches)
    private static readonly MiRna SimQuery = CreateMiRna("q", "AGGGGGGG");
    private static readonly MiRna SimM1 = CreateMiRna("m1", "AGGGGGCC"); // seed GGGGGCC, 2 mm
    private static readonly MiRna SimM2 = CreateMiRna("m2", "AGGGGCCC"); // seed GGGGCCC, 3 mm

    [Test]
    public void FindSimilarMiRnas_SeedMismatchAtThreshold_IsIncluded()
    {
        // m1 has exactly 2 mismatches; with maxMismatches=2 the inclusive '<=' keeps it,
        // a strict '<' mutant would drop it; an 'i>min' loop-bound mutant (mismatches≡0)
        // or a '>' mutant would additionally let m2 (3 mm) in.
        var hits = FindSimilarMiRnas(SimQuery, new[] { SimM1, SimM2 }, maxMismatches: 2)
            .Select(m => m.Name).ToList();

        Assert.That(hits, Does.Contain("m1"));
        Assert.That(hits, Does.Not.Contain("m2"));
    }

    [Test]
    public void FindSimilarMiRnas_BelowThreshold_TwoMismatchExcluded()
    {
        // With maxMismatches=1, m1's 2 mismatches exceed the cutoff ⇒ excluded.
        var hits = FindSimilarMiRnas(SimQuery, new[] { SimM1, SimM2 }, maxMismatches: 1)
            .Select(m => m.Name).ToList();

        Assert.That(hits, Is.Empty);
    }

    [Test]
    public void FindSimilarMiRnas_SeedMismatchCountIsExact()
    {
        // Sanity check on the reference seeds so the threshold tests rest on known distances.
        Assert.That(SimQuery.SeedSequence, Is.EqualTo("GGGGGGG"));
        Assert.That(SimM1.SeedSequence, Is.EqualTo("GGGGGCC"));
        Assert.That(SimM2.SeedSequence, Is.EqualTo("GGGGCCC"));
    }

    #endregion
}
