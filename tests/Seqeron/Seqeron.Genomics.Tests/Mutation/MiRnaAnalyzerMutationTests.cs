using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Annotation;
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
///  - Local accessibility (less local structure ⇒ more accessible site): Kertesz et al.
///    (2007) Nat Genet 39:1278–1284 (PITA); Watson–Crick pairing per Crick (1966).
///  - Seed families (shared 7-nt seed, single-mismatch neighbours): Bartel (2009) Cell
///    136:215–233; Lewis et al. (2005) Cell 120:15–20.
/// </summary>
[TestFixture]
public class MiRnaAnalyzerMutationTests
{
    private const double Tol = 1e-9;

    #region AnalyzeTargetContext — positional thresholds, AU bonus, context score

    // Model (algorithm doc §"Target Context Analysis"):
    //   window  = mrna[max(0,start-W) .. min(len,end+W))     (W = contextWindow = 30)
    //   auContent = #(A|U in window) / window.Length
    //   nearStart = start < len*0.15 ;  nearEnd = end > len*0.85
    //   contextScore = auContent*0.5  (+0.3 only when NOT nearStart AND NOT nearEnd)
    //                  clamped to ≤ 1.0
    // A 100-nt poly(A) carrier makes auContent ≡ 1.0 for every window (auContent term is
    // isolated), and len*0.15 = 15.0, len*0.85 = 85.0 are exact integers so the thresholds
    // are unambiguous.

    private static readonly string PolyA100 = new string('A', 100);

    [Test]
    public void AnalyzeTargetContext_MiddleSite_GetsAuBonus_ExactScore()
    {
        // start=40,end=50: nearStart=40<15=false, nearEnd=50>85=false ⇒ +0.3 bonus.
        var ctx = AnalyzeTargetContext(PolyA100, 40, 50);

        Assert.That(ctx.AuContent, Is.EqualTo(1.0).Within(Tol));
        Assert.That(ctx.NearStart, Is.False);   // kills start>… and start/0.15 mutants
        Assert.That(ctx.NearEnd, Is.False);
        Assert.That(ctx.ContextScore, Is.EqualTo(0.8).Within(Tol)); // 1.0*0.5 + 0.3
    }

    [Test]
    public void AnalyzeTargetContext_StartProximal_NoBonus_ExactScore()
    {
        // start=5<15 ⇒ nearStart=true; end=10>85=false. !nearEnd && !nearStart = false ⇒ no bonus.
        var ctx = AnalyzeTargetContext(PolyA100, 5, 10);

        Assert.That(ctx.NearStart, Is.True);    // kills start>… mutant (5>15=false)
        Assert.That(ctx.NearEnd, Is.False);
        Assert.That(ctx.ContextScore, Is.EqualTo(0.5).Within(Tol)); // 1.0*0.5, no bonus
        // With the && replaced by ||, !nearEnd||!nearStart = true ⇒ bonus would push this to 0.8.
    }

    [Test]
    public void AnalyzeTargetContext_EndProximal_NoBonus_ExactScore()
    {
        // start=90 ⇒ nearStart=false; end=95>85 ⇒ nearEnd=true ⇒ no bonus.
        var ctx = AnalyzeTargetContext(PolyA100, 90, 95);

        Assert.That(ctx.NearStart, Is.False);
        Assert.That(ctx.NearEnd, Is.True);      // kills end/0.85 mutant (95>117.6=false)
        Assert.That(ctx.ContextScore, Is.EqualTo(0.5).Within(Tol));
    }

    [Test]
    public void AnalyzeTargetContext_StartExactlyAtThreshold_IsNotNearStart()
    {
        // start = 15 == len*0.15. Strict '<' ⇒ NOT nearStart; a '<=' mutant would flip it.
        var ctx = AnalyzeTargetContext(PolyA100, 15, 20);

        Assert.That(ctx.NearStart, Is.False);   // kills start<=len*0.15 mutant
        Assert.That(ctx.NearEnd, Is.False);
        Assert.That(ctx.ContextScore, Is.EqualTo(0.8).Within(Tol)); // bonus applies
    }

    [Test]
    public void AnalyzeTargetContext_EndExactlyAtThreshold_IsNotNearEnd()
    {
        // end = 85 == len*0.85. Strict '>' ⇒ NOT nearEnd; a '>=' mutant would flip it.
        var ctx = AnalyzeTargetContext(PolyA100, 50, 85);

        Assert.That(ctx.NearStart, Is.False);
        Assert.That(ctx.NearEnd, Is.False);     // kills end>=len*0.85 mutant
        Assert.That(ctx.ContextScore, Is.EqualTo(0.8).Within(Tol));
    }

    [Test]
    public void AnalyzeTargetContext_AuContentDrivesHalfWeight()
    {
        // GC-only window ⇒ auContent = 0 ⇒ contextScore = 0*0.5 (+0.3 bonus, middle site) = 0.3.
        // This pins the auContent*0.5 term (separates *0.5 from /0.5: 0/0.5 == 0 either way,
        // but the additive bonus is isolated to exactly 0.3 here).
        string polyGc = string.Concat(Enumerable.Repeat("GC", 50)); // 100 nt, no A/U
        var ctx = AnalyzeTargetContext(polyGc, 40, 50);

        Assert.That(ctx.AuContent, Is.EqualTo(0.0).Within(Tol));
        Assert.That(ctx.ContextScore, Is.EqualTo(0.3).Within(Tol)); // 0*0.5 + 0.3
    }

    #endregion

    #region CalculateSiteAccessibility — windowed structure density

    // Model (algorithm doc §"Site accessibility"):
    //   guard: empty || siteStart<0 || siteEnd>=len ⇒ 0
    //   window  = mrna[max(0,siteStart-50) .. min(len,siteEnd+50))
    //   structureScore = #{(i,j): j>=i+4, CanPair(w_i,w_j) AND NOT G:U-wobble}
    //   maxPairs = (W*(W-4))/2  ;  density = structureScore / max(1,maxPairs)
    //   accessibility = max(0, 1 - density*10)
    //
    // Reference sequence GAAAAUAAAC (len 10). The window covers the whole sequence.
    // Watson–Crick (non-wobble) pairs with j>=i+4:
    //   (0=G,9=C) G-C ; (1=A,5=U) A-U  ⇒ structureScore = 2  (G0:U5 is wobble, excluded).
    //   maxPairs = (10*6)/2 = 30 ;  accessibility = 1 - (2/30)*10 = 1 - 20/30.
    private const string AccSeq = "GAAAAUAAAC";
    private static readonly double AccExpected = 1.0 - 20.0 / 30.0; // ≈ 0.333333…

    [Test]
    public void CalculateSiteAccessibility_KnownWindow_ExactValue()
    {
        double acc = CalculateSiteAccessibility(AccSeq, 2, 7);

        Assert.That(acc, Is.EqualTo(AccExpected).Within(Tol));
        // Pins: pair-count loop bounds, the CanPair AND !wobble predicate, the structureScore++
        // block, maxPairs = (W*(W-4))/2, density division, and the 1 - density*10 form.
    }

    [Test]
    public void CalculateSiteAccessibility_SiteStartZero_StillComputes()
    {
        // siteStart == 0 is valid (guard is strict siteStart < 0). A '<=' mutant returns 0.
        double acc = CalculateSiteAccessibility(AccSeq, 0, 7);
        Assert.That(acc, Is.EqualTo(AccExpected).Within(Tol));
    }

    [Test]
    public void CalculateSiteAccessibility_SiteEndEqualsLength_ReturnsZeroByGuard()
    {
        // siteEnd == len triggers the siteEnd >= len guard ⇒ 0. A '>' mutant would compute instead.
        double acc = CalculateSiteAccessibility(AccSeq, 2, AccSeq.Length);
        Assert.That(acc, Is.EqualTo(0.0).Within(Tol));
    }

    [Test]
    public void CalculateSiteAccessibility_NegativeStart_WithRoomAfter_ReturnsZero()
    {
        // siteStart < 0 with siteEnd < len: only the (siteStart<0) disjunct fires.
        // The '||'→'&&' mutant on the first guard term would let this fall through and compute.
        double acc = CalculateSiteAccessibility(AccSeq, -1, 5);
        Assert.That(acc, Is.EqualTo(0.0).Within(Tol));
    }

    [Test]
    public void CalculateSiteAccessibility_OffsetWindow_NoStructure_FullyAccessible()
    {
        // siteStart=60 ⇒ window start = 60-50 = 10 (non-zero), exercising end-start length.
        // poly(A) ⇒ no base pairs ⇒ density 0 ⇒ accessibility = 1.0. An 'end+start' mutant
        // requests an out-of-range substring length and throws.
        string polyA = new string('A', 120);
        double acc = CalculateSiteAccessibility(polyA, 60, 65);
        Assert.That(acc, Is.EqualTo(1.0).Within(Tol));
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
