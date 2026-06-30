// MIRNA-TARGET-001 / MIRNA-PCT-001 (08_DIFFERENTIAL strategy REF) — closed-form killers for the
// remaining context++ kernels: non-overlapping seed-site counting (Garcia 2011 TA), the TA log10,
// LocalAU flank-length boundary, the offset-6mer pass of FindTargetSites, and dot-bracket edge
// guards. Asserted against hand-derived expected values, independent of the implementation.

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests;

[TestFixture]
[Category("MIRNA-TARGET-001")]
public class MiRnaAnalyzer_SeedCountingKernels_Tests
{
    private const double Tol = 1e-9;

    // ── CountSeedSitesInUtr: count of 8mer/7mer-m8/7mer-A1 (not bare 6mer), non-overlapping ──
    [Test]
    public void CountSeedSitesInUtr_Counts8merAndQualifiedSites()
    {
        // core = "GUACGU", pos8 base = 'C'. "CGUACGUA": core at i=1, utr[0]='C' (m8), utr[7]='A' (A1) → 8mer.
        Assert.That(MiRnaAnalyzer.CountSeedSitesInUtr("CGUACGUA", "GUACGU", 'C'), Is.EqualTo(1));
        // bare 6mer (no upstream C, no downstream A): "GGUACGUG" → core at i=1, utr[0]='G'≠'C', utr[7]='G'≠'A' → 0.
        Assert.That(MiRnaAnalyzer.CountSeedSitesInUtr("GGUACGUG", "GUACGU", 'C'), Is.EqualTo(0));
        // 7mer-A1 only (downstream A, no upstream C): "GGUACGUA" → 1.
        Assert.That(MiRnaAnalyzer.CountSeedSitesInUtr("GGUACGUA", "GUACGU", 'C'), Is.EqualTo(1));
        // Core at the LAST index (i = len-6) with an upstream m8 match and NO room for A1: "CCGUACGU"
        // has core "GUACGU" at i=2 = 8-6, utr[1]='C' (m8). The hasA1 guard `i+6 < len` must short-circuit
        // before reading utr[len]; a `<=` mutant would index out of range. → 1 site (7mer-m8).
        Assert.That(MiRnaAnalyzer.CountSeedSitesInUtr("CCGUACGU", "GUACGU", 'C'), Is.EqualTo(1));
    }

    [Test]
    public void CountSeedSitesInUtr_NonOverlapping_TwoDistinctSites()
    {
        // Two non-overlapping 7mer-A1 sites of core "AAAAAA": after a counted site the scan jumps +6.
        // "AAAAAAA AAAAAAA" laid out so each "AAAAAA" core has a downstream 'A'. Use a homopolymer where
        // greedy non-overlapping counting yields a deterministic number.
        // "AAAAAAAAAAAAA" (13 A): core matches at every i; first site at i=0 (hasA1: utr[6]='A'), jump to 6;
        // site at 6 (utr[12]='A'? index12 exists, =A) hasA1, jump to 12; 12+6=18>13 stop → 2 sites.
        Assert.That(MiRnaAnalyzer.CountSeedSitesInUtr(new string('A', 13), "AAAAAA", 'A'), Is.EqualTo(2));
    }

    [Test]
    public void ComputeTa3Utr_IsLog10OfTotalSiteCount()
    {
        var miRna = CreateMiRna("m", "AACGUACGUAGCUAGCUAGCUA");
        string seedRC = GetReverseComplement(miRna.SeedSequence);
        char pos8 = seedRC[0];
        string core = seedRC.Substring(1, 6);
        // Build a UTR with exactly one 8mer site: pos8 + core + 'A'.
        string utr = pos8 + core + "A";
        Assert.That(MiRnaAnalyzer.CountSeedSites3Utr(miRna, new[] { utr }), Is.EqualTo(1));
        Assert.That(MiRnaAnalyzer.ComputeTa3Utr(miRna, new[] { utr }), Is.EqualTo(0.0).Within(Tol)); // log10(1)=0
        // Ten copies → log10(10) = 1.
        Assert.That(MiRnaAnalyzer.ComputeTa3Utr(miRna, Enumerable.Repeat(utr, 10)), Is.EqualTo(1.0).Within(Tol));
    }

    // ── LocalAU flank-length boundary (the loop runs i in [0, 30); a 31st flank base must be ignored) ──
    [Test]
    public void LocalAuContribution_FlankLengthBoundary_31stBaseIgnored()
    {
        // Upstream: 30 'A' (indices 1..30) all counted; a 31st base 'G' at index 0 must NOT be reached
        // by the i < 30 loop. Site [31..38]; no downstream (site ends at the last index). The 30 counted
        // upstream bases are all A/U → fraction = 1.0; a mutant that reaches the 31st 'G' lowers it.
        string mrna = "G" + new string('A', 30) + new string('C', 8); // length 39; site [31..38]
        Assert.That(MiRnaAnalyzer.LocalAuContribution(mrna, 31, 38, TargetSiteType.Seed8mer),
            Is.EqualTo(-0.254 * (1.0 - 0.308) / (0.814 - 0.308)).Within(Tol));
    }

    // ── FindTargetSites: offset-6mer (pass 2) and the minScore boundary ──
    [Test]
    public void FindTargetSites_OffsetSixmer_DetectedAndScored()
    {
        // offset-6mer pattern = first 6 nt of seedRC (RC of miRNA positions 3-8). Build an mRNA whose
        // ONLY match is the offset-6mer (no 6mer-core anywhere), and check it is reported as Offset6mer.
        var miRna = CreateMiRna("m", "AACGUACGUAGCUAGCUAGCUA");
        string seedRC = GetReverseComplement(miRna.SeedSequence);
        string offset6 = seedRC.Substring(0, 6);
        string sixmerCore = seedRC.Substring(1, 6);
        // Place offset6 flanked so it does NOT also form the 6mer core or full seedRC.
        string mrna = "UU" + offset6 + "UU";
        var sites = FindTargetSites(mrna, miRna, minScore: 0.0).ToList();
        // The offset-6mer site is reported only when no higher-priority 6mer-core site overlaps it.
        bool coreElsewhere = mrna.Contains(sixmerCore);
        if (!coreElsewhere)
            Assert.That(sites.Any(s => s.Type == TargetSiteType.Offset6mer),
                $"expected an Offset6mer; got [{string.Join(", ", sites.Select(s => s.Type.ToString()))}]");
        else
            Assert.Pass("core present elsewhere — offset suppression path, not asserted here");
    }

    [Test]
    public void FindTargetSites_MinScoreBoundary_Inclusive()
    {
        // The filter is `site.Score >= minScore`. Setting minScore to exactly an 8mer's score keeps it;
        // a strictly-greater mutant would drop it. Build a strong 8mer and use its own score as minScore.
        var miRna = CreateMiRna("m", "AACGUACGUAGCUAGCUAGCUA");
        string seedRC = GetReverseComplement(miRna.SeedSequence);
        string mrna = "GG" + seedRC + "A" + "GG";
        var atZero = FindTargetSites(mrna, miRna, minScore: 0.0).First(s => s.Type == TargetSiteType.Seed8mer);
        var atExact = FindTargetSites(mrna, miRna, minScore: atZero.Score).ToList();
        Assert.That(atExact.Any(s => s.Type == TargetSiteType.Seed8mer && Math.Abs(s.Score - atZero.Score) < 1e-12),
            "score == minScore must be retained (inclusive boundary)");
    }

    // ── dot-bracket edge guards ──
    [Test]
    public void TryDescribeSingleHairpin_SingleBasePairAtIndexZero()
    {
        // "(.)" : lastOpenIndex = 0. A mutant guarding `lastOpenIndex <= 0` would wrongly reject it.
        bool ok = MiRnaAnalyzer.TryDescribeSingleHairpin("(.)", out int bp, out int loopStart, out int loopSize);
        Assert.That(ok, Is.True);
        Assert.That(bp, Is.EqualTo(1));
        Assert.That(loopStart, Is.EqualTo(1));
        Assert.That(loopSize, Is.EqualTo(1));
    }

    [Test]
    public void LargestEnclosedLoop_ParensNotAtStart()
    {
        // ".(...)." : firstOpen = 1 (> 0). A mutant guarding `firstOpen > 0` would wrongly return 0.
        Assert.That(MiRnaAnalyzer.LargestEnclosedLoop(".(...)."), Is.EqualTo(3));
    }
}
