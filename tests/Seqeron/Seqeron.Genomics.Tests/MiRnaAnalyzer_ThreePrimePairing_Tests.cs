// MIRNA-CONTEXT-001 (08_DIFFERENTIAL strategy REF) — spec-based killers for the 3'-supplementary
// pairing CORE SCORER, asserted against the PUBLISHED Grimson et al. (2007) rule (NOT the perl prep):
//   over an alignment, score only runs of >= 2 consecutive Watson-Crick pairs; each pair adds 1.0 when
//   it falls in the supplementary window (offset-adjusted position 4..7) and 0.5 otherwise; take the
//   best such run. Base codes A=1,U=2,C=3,G=4,N=5; a pair is A:U (1*2=2) or G:C (3*4=12).
// The expected values are computed BY HAND from that rule on tiny arrays the test fully controls, so a
// mutated index / threshold / product / run-commit in PairingRunScore diverges -> killed. This is the
// model-vs-implementation differential, independent of the (idiosyncratic) subseq-extraction prep.

using NUnit.Framework;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

[TestFixture]
[Category("MIRNA-CONTEXT-001")]
public class MiRnaAnalyzer_ThreePrimePairing_Tests
{
    private const double Tol = 1e-9;
    // Codes: G=4, C=3, A=1, U=2, N=5, mismatch=0.

    [Test]
    public void PairingRunScore_ContiguousRun_WindowWeighting()
    {
        // UTR=GCGCA, MIRNA=CGCGU, offset 0, overhang 0, gapOnTop. All 5 pair (G:C,C:G,G:C,C:G,A:U).
        // posCheck = i; i=0..3 → 0.5 each (outside [4,7]); i=4 → 1.0 (in [4,7]). Total = 4*0.5 + 1.0 = 3.0.
        var utr = new[] { 4, 3, 4, 3, 1 };
        var mir = new[] { 3, 4, 3, 4, 2 };
        Assert.That(MiRnaAnalyzer.PairingRunScore(utr, mir, 0, 0, true), Is.EqualTo(3.0).Within(Tol));
    }

    [Test]
    public void PairingRunScore_RunMustBeAtLeastTwo_AndBreakCommits()
    {
        // UTR=G C X G C, MIRNA=C G C G G. i0:G:C pair, i1:C:G pair (run=2, 0.5+0.5=1.0), i2: 0 vs 3 no
        // pair → commit 1.0; i3: G:G (16) no pair; i4: C:G pair but lone (run=1, not committed).
        var utr = new[] { 4, 3, 0, 4, 3 };
        var mir = new[] { 3, 4, 3, 4, 4 };
        Assert.That(MiRnaAnalyzer.PairingRunScore(utr, mir, 0, 0, true), Is.EqualTo(1.0).Within(Tol),
            "the length-2 run scores 1.0; the lone trailing pair is not counted");
    }

    [Test]
    public void PairingRunScore_OverhangShiftsTheWindow()
    {
        // Same as the first case but overhang=2 → posCheck = i-2 = -2..2, none in [4,7] → all 0.5.
        var utr = new[] { 4, 3, 4, 3, 1 };
        var mir = new[] { 3, 4, 3, 4, 2 };
        Assert.That(MiRnaAnalyzer.PairingRunScore(utr, mir, 0, 2, true), Is.EqualTo(2.5).Within(Tol),
            "5 pairs * 0.5 = 2.5 when the supplementary window is shifted out");
    }

    [Test]
    public void PairingRunScore_GapOnBottom_IndexesUtrByOffset()
    {
        // gapOnTop=false: u = UTR[i+offset], m = MIRNA[i], posCheck = i - overhang. offset 1:
        // UTR=N G C G C, MIRNA=C G C G (len 4). i0:UTR[1]=G:C pair, i1:UTR[2]=C:G, i2:UTR[3]=G:C,
        // i3:UTR[4]=C:G. limit=min(UTR.Len-1-1, MIRNA.Len-1)=min(3,3)=3. posCheck=i; i=0..3 → 0.5 each
        // = 2.0.
        var utr = new[] { 5, 4, 3, 4, 3 };
        var mir = new[] { 3, 4, 3, 4 };
        Assert.That(MiRnaAnalyzer.PairingRunScore(utr, mir, 1, 0, false), Is.EqualTo(2.0).Within(Tol));
    }

    [Test]
    public void PairingRunScore_WindowUpperBoundIsInclusiveAtSeven()
    {
        // An 8-pair run (G:C x8), offset 0, overhang 0 → posCheck = 0..7. Positions 4,5,6,7 are inside
        // the inclusive [4,7] window (1.0 each); 0,1,2,3 are outside (0.5 each). Total = 4*0.5 + 4*1.0 = 6.0.
        // A mutant making the upper bound exclusive (`< 7`) would score position 7 at 0.5 → 5.5.
        var utr = new[] { 4, 3, 4, 3, 4, 3, 4, 3 };
        var mir = new[] { 3, 4, 3, 4, 3, 4, 3, 4 };
        Assert.That(MiRnaAnalyzer.PairingRunScore(utr, mir, 0, 0, true), Is.EqualTo(6.0).Within(Tol));
    }

    [Test]
    public void PairingRunScore_GapOnTop_PositionUsesOffset()
    {
        // gapOnTop=true, offset=2, overhang=0: u=UTR[i], m=MIRNA[i+2], posCheck = i + offset - overhang
        // = i + 2. UTR=GCGC; MIRNA="AA"+CGCG so MIRNA[i+2] pairs with UTR[i] (G:C). limit = min(6-1-2,3)=3.
        // posCheck: i=0→2 (0.5), i=1→3 (0.5), i=2→4 (1.0), i=3→5 (1.0) → 3.0. A mutant dropping the
        // +offset term (or flipping the overhang sign) would mis-window every position.
        var utr = new[] { 4, 3, 4, 3 };
        var mir = new[] { 1, 1, 3, 4, 3, 4 };
        Assert.That(MiRnaAnalyzer.PairingRunScore(utr, mir, 2, 0, true), Is.EqualTo(3.0).Within(Tol));
    }

    [Test]
    public void PairingRunScore_NoRun_IsZero()
    {
        // Alternating single pairs separated by mismatches → no run of >= 2 → 0.
        var utr = new[] { 4, 0, 4, 0, 4 };
        var mir = new[] { 3, 0, 3, 0, 3 };
        Assert.That(MiRnaAnalyzer.PairingRunScore(utr, mir, 0, 0, true), Is.EqualTo(0.0).Within(Tol));
    }
}
