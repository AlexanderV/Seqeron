// SEQ-GC-PROFILE-001 — GC Content Profile (sliding-window GC content)
// Evidence: docs/Evidence/SEQ-GC-PROFILE-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-GC-PROFILE-001.md
// Source: Wikipedia, GC-content (citing primary literature),
//         https://en.wikipedia.org/wiki/GC-content (accessed 2026-06-14).
//         Biopython Bio.SeqUtils.gc_fraction, Cock P.J.A. et al. (2009)
//         Bioinformatics 25(11):1422-1423, doi:10.1093/bioinformatics/btp163.

using System.Linq;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceStatistics_CalculateGcContentProfile_Tests
{
    // Expected GC% values are derived by hand from GC% = (G+C)/(A+T+G+C)×100,
    // computed independently of the implementation, per the Evidence datasets.
    private const double Tolerance = 1e-10;

    #region CalculateGcContentProfile (sliding window)

    // M1 — all-GC window: 10/10×100 = 100.0 (not the fraction 1.0).
    // Evidence: Wikipedia (G+C)/(A+T+G+C)×100.
    [Test]
    public void CalculateGcContentProfile_AllGcWindow_Returns100Percent()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGGGGGGGGG", 10).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(1),
                "n=10,w=10,step=1 → exactly one window (INV-03)");
            Assert.That(profile[0], Is.EqualTo(100.0).Within(Tolerance),
                "GGGGGGGGGG → 10/10×100 = 100.0 GC% (INV-01)");
        });
    }

    // M2 — half-GC windows: every ATGC window is 2/4×100 = 50.0.
    // Evidence: Wikipedia formula; Biopython gc_fraction("ACTG")=0.50 ×100.
    [Test]
    public void CalculateGcContentProfile_AtgcRepeats_Returns50PercentPerWindow()
    {
        // "ATGCATGCATGC", window 4, step 4 → three disjoint ATGC windows.
        var profile = SequenceStatistics.CalculateGcContentProfile("ATGCATGCATGC", 4, 4).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(3),
                "n=12,w=4,step=4 → offsets 0,4,8 = 3 windows (INV-03)");
            Assert.That(profile.All(p => System.Math.Abs(p - 50.0) <= Tolerance), Is.True,
                "each ATGC window → 2 GC / 4 bases × 100 = 50.0 GC%");
        });
    }

    // M3 — exact mixed profile: GGGAAATGCC, w=4, step=3 → GGGA, AAAT, TGCC.
    // Evidence: Wikipedia formula per window. Values 75/0/75 are not what a
    // fraction (0.75/0/0.75) or an off-by-one window would produce.
    [Test]
    public void CalculateGcContentProfile_MixedSequence_ReturnsExactProfile()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 3).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(3),
                "n=10,w=4,step=3 → offsets 0,3,6 = 3 windows (INV-03)");
            Assert.That(profile[0], Is.EqualTo(75.0).Within(Tolerance),
                "GGGA → 3 GC / 4 × 100 = 75.0 GC%");
            Assert.That(profile[1], Is.EqualTo(0.0).Within(Tolerance),
                "AAAT → 0 GC / 4 × 100 = 0.0 GC%");
            Assert.That(profile[2], Is.EqualTo(75.0).Within(Tolerance),
                "TGCC → 3 GC / 4 × 100 = 75.0 GC%");
        });
    }

    // M4 — N excluded from the denominator: GGAN → 2 GC / 3 standard bases × 100.
    // Evidence: Biopython gc_fraction("ACTGN")=0.50 under default "remove" (N removed
    // from length). Counting N in the denominator would give 2/4×100 = 50.0 (wrong).
    [Test]
    public void CalculateGcContentProfile_AmbiguousN_ExcludedFromDenominator()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGAN", 4).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(1), "single window (INV-03)");
            Assert.That(profile[0], Is.EqualTo(200.0 / 3.0).Within(Tolerance),
                "GGAN → 2 GC / 3 standard bases (N excluded) × 100 = 66.66… GC% (INV-02)");
        });
    }

    // M5 — RNA U is a non-GC base equivalent to T: GGAU → 2 GC / 4 × 100 = 50.0.
    // Evidence: Biopython gc_fraction("GGAUCUUCGGAUCU")=0.50 (U non-GC).
    [Test]
    public void CalculateGcContentProfile_RnaUracil_TreatedAsNonGc()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGAU", 4).ToArray();

        Assert.That(profile[0], Is.EqualTo(50.0).Within(Tolerance),
            "GGAU → G,G are GC, A and U are non-GC → 2/4×100 = 50.0 GC% (INV-04)");
    }

    // M6 — window count obeys INV-03 for several step sizes.
    // Evidence: sliding-window enumeration count = ⌊(n-w)/step⌋+1.
    [Test]
    public void CalculateGcContentProfile_WindowCount_MatchesInvariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 1).Count(),
                Is.EqualTo(7), "n=10,w=4,step=1 → 7 windows (INV-03)");
            Assert.That(SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 3).Count(),
                Is.EqualTo(3), "n=10,w=4,step=3 → offsets 0,3,6 = 3 windows (INV-03)");
            Assert.That(SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 2).Count(),
                Is.EqualTo(4), "n=10,w=4,step=2 → offsets 0,2,4,6 = 4 windows (INV-03)");
        });
    }

    #endregion

    #region Edge cases and invariants

    // S1 — windowSize greater than sequence length yields an empty profile.
    // Evidence: INV-03 (no full window when W > n).
    [Test]
    public void CalculateGcContentProfile_WindowLargerThanSequence_ReturnsEmpty()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("ATGC", 100).ToArray();

        Assert.That(profile, Is.Empty,
            "W=100 > n=4 → no full window → empty profile (INV-03)");
    }

    // S2 — windowSize equal to length yields exactly one whole-sequence window.
    // Evidence: INV-03; GGCC → 4/4×100 = 100.0.
    [Test]
    public void CalculateGcContentProfile_WindowEqualsLength_ReturnsSingleValue()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGCC", 4).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(1), "W == n → one window (INV-03)");
            Assert.That(profile[0], Is.EqualTo(100.0).Within(Tolerance),
                "GGCC → 4 GC / 4 × 100 = 100.0 GC%");
        });
    }

    // S3 — null and empty input yield empty profiles.
    // Evidence: guarded input (§3.3).
    [Test]
    public void CalculateGcContentProfile_NullOrEmpty_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateGcContentProfile(null!, 4), Is.Empty,
                "null sequence → empty profile");
            Assert.That(SequenceStatistics.CalculateGcContentProfile("", 4), Is.Empty,
                "empty sequence → empty profile");
        });
    }

    // S4 — window with no standard base (all-N) → 0 (zero-division convention).
    // Evidence: Assumption A1 (repository convention, matches SEQ-GC-ANALYSIS-001).
    [Test]
    public void CalculateGcContentProfile_AllAmbiguousWindow_ReturnsZero()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("NNNN", 4).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(profile.Length, Is.EqualTo(1), "single window (INV-03)");
            Assert.That(profile[0], Is.EqualTo(0.0).Within(Tolerance),
                "no standard base → 0 (zero-division convention, INV-05)");
        });
    }

    // C1 — case-insensitivity: lowercase input produces the same profile as uppercase.
    // Evidence: implementation case-folds before counting (§3.3).
    [Test]
    public void CalculateGcContentProfile_LowercaseInput_MatchesUppercase()
    {
        var upper = SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", 4, 1).ToArray();
        var lower = SequenceStatistics.CalculateGcContentProfile("gggaaatgcc", 4, 1).ToArray();

        Assert.That(lower, Is.EqualTo(upper).Within(Tolerance),
            "case-folded counting → lowercase profile equals uppercase profile");
    }

    // C2 — INV-01: every window value is bounded in [0, 100].
    // Evidence: numerator ≤ denominator, then ×100 (Wikipedia formula).
    [Test]
    public void CalculateGcContentProfile_AllWindows_BoundedZeroToHundred()
    {
        var profile = SequenceStatistics
            .CalculateGcContentProfile("ATGCGGGGAAAACCCCATGCATGC", 6, 1).ToArray();

        Assert.That(profile, Is.Not.Empty, "profile must contain windows");
        Assert.That(profile.All(p => p >= -Tolerance && p <= 100.0 + Tolerance), Is.True,
            "every GC% value in [0,100] (INV-01)");
    }

    #endregion
}
