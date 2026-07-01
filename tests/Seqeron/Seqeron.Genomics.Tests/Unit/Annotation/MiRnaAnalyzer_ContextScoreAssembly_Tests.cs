// MIRNA-CONTEXT-001 / MIRNA-CLEAVAGE-001 — killers for the context++ score ASSEMBLY (the partial-sum
// of all per-feature contributions), the SA window-fit boundary, the CNNC max-offset boundary, and the
// short-input guards. Runs under the Permissive LimitationPolicy (test bootstrap) so the full
// context++ score returns its partial value instead of throwing on omitted features.

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

[TestFixture]
[Category("MIRNA-CONTEXT-001")]
public class MiRnaAnalyzer_ContextScoreAssembly_Tests
{
    private const double Tol = 1e-9;
    private static readonly MiRna Mir = CreateMiRna("m", "AACGUACGUAGCUAGCUAGCUA");

    [Test]
    public void ContextScorePartial_EqualsSumOfAllContributions()
    {
        string seedRC = GetReverseComplement(Mir.SeedSequence);
        // Long flanks so SA and local-AU windows fit; a single 8mer in the middle.
        string mrna = new string('G', 20) + seedRC + "A" + new string('U', 20);
        var site = FindTargetSites(mrna, Mir, minScore: 0.0).First(s => s.Type == TargetSiteType.Seed8mer);

        var sc = ScoreTargetSiteContextPlusPlus(mrna, Mir, site);

        double sum = sc.Intercept + sc.LocalAuContribution + sc.SRna1Contribution + sc.SRna8Contribution
            + sc.Site8Contribution + sc.SaContribution + sc.ThreePrimePairingContribution + sc.MinDistContribution
            + sc.Len3UtrContribution + sc.Off6mContribution + sc.SpsContribution + sc.TaContribution
            + sc.LenOrfContribution + sc.Orf8mContribution + sc.PctContribution;
        Assert.That(sc.ContextScorePartial, Is.EqualTo(sum).Within(Tol),
            "ContextScorePartial must be the exact sum of every per-feature contribution");
    }

    // SA is included iff the 14-nt window fits: windowStart0 = siteStart-6 >= 0 and windowEnd0 =
    // siteStart+7 < length. siteStart=6 → windowStart0=0 (included); siteStart=5 → windowStart0=-1
    // (omitted, returns 0). Pins the `windowStart0 < 0` boundary (a `<= 0` mutant rejects siteStart=6).
    [Test]
    public void SaContribution_WindowStartBoundary()
    {
        string mrna = new string('A', 40);
        MiRnaAnalyzer.SaContribution(mrna, 6, TargetSiteType.Seed8mer, out bool includedAt6);
        MiRnaAnalyzer.SaContribution(mrna, 5, TargetSiteType.Seed8mer, out bool includedAt5);
        Assert.That(includedAt6, Is.True, "siteStart=6 → windowStart0=0 → window fits");
        Assert.That(includedAt5, Is.False, "siteStart=5 → windowStart0=-1 → omitted");
    }

    [Test]
    public void SaContribution_WindowDoesNotFit_ReturnsZero()
    {
        double sa = MiRnaAnalyzer.SaContribution(new string('A', 10), 1, TargetSiteType.Seed8mer, out bool inc);
        Assert.That(inc, Is.False);
        Assert.That(sa, Is.EqualTo(0.0));
    }

    // CNNC detection window = offsets {16,17,18} downstream of the Drosha cut. Placing a CNNC ONLY at
    // the MAX offset (18) pins the inclusive `offset <= CnncMaxNtDownstreamOfDroshaCut` bound: a
    // `offset < max` mutant never tests offset 18 and reports no motif.
    [Test]
    public void Cnnc_AtMaxOffset18_IsDetected()
    {
        // basalJunction=5 → droshaCut5Prime=16; offset 18 → start=34; CNNC = C at 34 and 37.
        var s = new char[60];
        Array.Fill(s, 'A');
        s[34] = 'C'; s[37] = 'C'; // CNNC only at offset 18
        var c = MiRnaAnalyzer.PredictDroshaDicerCleavage(new string(s), 5);
        Assert.That(c!.Value.HasCnncMotif, Is.True, "CNNC at the max offset (18) must be detected");
    }

    // ── short-input guards ──
    [Test]
    public void SRna1Contribution_EmptyMiRna_IsZero()
        => Assert.That(MiRnaAnalyzer.SRna1Contribution("", TargetSiteType.Seed8mer), Is.EqualTo(0.0));

    [Test]
    public void Off6mContribution_MiRnaShorterThan8_IsZero()
        => Assert.That(MiRnaAnalyzer.Off6mContribution("AAAAAAAAAA", "ACGUACG", TargetSiteType.Seed8mer),
            Is.EqualTo(0.0), "miRNA length 7 (< 8) → no offset-6mer pattern → 0");

    [Test]
    public void CountSeedSites3Utr_SeedShorterThan7_IsZero()
    {
        // A miRNA whose seed is shorter than 7 → no site set → 0.
        var shortMir = CreateMiRna("s", "ACGUA"); // seed (pos 2-8) shorter than 7
        Assert.That(MiRnaAnalyzer.CountSeedSites3Utr(shortMir, new[] { new string('A', 50) }), Is.EqualTo(0));
    }
}
