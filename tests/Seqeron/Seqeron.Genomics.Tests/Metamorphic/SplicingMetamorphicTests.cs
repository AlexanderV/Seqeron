using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Splicing area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SPLICE-DONOR-001 — 5' donor splice-site prediction (Splicing).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 77.
///
/// API under test (SpliceSitePredictor.FindDonorSites):
///   A canonical donor requires the invariant GU dinucleotide; its score is the average
///   match of the surrounding bases (offsets −3..+5 around the G of GU) to the consensus
///   PWM "MAG|GURAGU". Bases outside that 9-position window do not enter the score.
///
/// Relations (derived from the PWM window, NOT from output):
///   • COMP (non-GT ⇒ no donor): without a GU dinucleotide there is no canonical donor site.
///   • MON  (closer to consensus ⇒ higher score): restoring consensus positions (keeping the
///          invariant GU fixed) monotonically raises the donor score, peaking at 1.0.
///   • INV  (downstream change ⇒ same score): edits beyond offset +5 lie outside the scoring
///          window and leave a donor's score unchanged.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class SplicingMetamorphicTests
{
    /// <summary>Score of the donor site anchored at index <paramref name="pos"/> (the G of GU).</summary>
    private static double DonorScoreAt(string sequence, int pos) =>
        SpliceSitePredictor.FindDonorSites(sequence, minScore: 0.0).Single(s => s.Position == pos).Score;

    #region SPLICE-DONOR-001 COMP — no GU dinucleotide ⇒ no canonical donor

    [Test]
    [Description("COMP: a canonical donor requires the invariant GU, so a sequence with no GU dinucleotide yields no donor sites.")]
    public void FindDonorSites_NoGu_ReturnsNothing()
    {
        // G's are always followed by C or A here — never U — so no canonical GU donor exists.
        foreach (var seq in new[] { "CAGCAGCAGCAGCAG", "AAGAAGAAGAAGAAG", "CCCCCCCCCCCC" })
            SpliceSitePredictor.FindDonorSites(seq, minScore: 0.0)
                .Should().BeEmpty(because: $"'{seq}' contains no GU dinucleotide to anchor a canonical donor");
    }

    #endregion

    #region SPLICE-DONOR-001 MON — closer to the consensus raises the score

    [Test]
    [Description("MON: with the invariant GU fixed at one locus, restoring consensus positions monotonically raises the donor score up to the perfect-consensus 1.0.")]
    public void FindDonorSites_CloserToConsensus_HigherScore()
    {
        // Donor anchored at index 3 (G of GU). Consensus −3..+5 = MAG|GURAGU → "CAGGUAAGU" scores 1.0.
        // Each variant restores one more consensus position while keeping GU (indices 3,4) intact.
        double perfect = DonorScoreAt("CAGGUAAGU", 3); // all 9 positions match consensus
        double oneOff  = DonorScoreAt("CAGGUACGU", 3); // +3 A→C (one mismatch)
        double twoOff  = DonorScoreAt("CAGGUACGC", 3); // +3 A→C, +5 U→C (two mismatches)
        double threeOff = DonorScoreAt("GAGGUACGC", 3); // −3 also broken (three mismatches)

        perfect.Should().Be(1.0, because: "every PWM position matches the consensus MAG|GURAGU");
        perfect.Should().BeGreaterThan(oneOff, because: "breaking a consensus position lowers the averaged PWM score");
        oneOff.Should().BeGreaterThan(twoOff, because: "a second consensus mismatch lowers the score further");
        twoOff.Should().BeGreaterThan(threeOff, because: "a third consensus mismatch lowers the score further");
    }

    #endregion

    #region SPLICE-DONOR-001 INV — edits beyond the scoring window don't change the score

    [Test]
    [Description("INV: the donor score uses only offsets −3..+5, so editing bases beyond offset +5 (downstream) leaves the donor's score unchanged.")]
    public void FindDonorSites_DownstreamEdits_DoNotChangeScore()
    {
        // Donor at index 3; its window is indices 0..8. Indices ≥9 are downstream of the window.
        const string core = "CAGGUAAGU";
        double baseScore = DonorScoreAt(core + "AAAAAA", 3);

        foreach (var tail in new[] { "CCCCCC", "GUGUGU", "UUUUUU", "ACGUAC" })
            DonorScoreAt(core + tail, 3).Should().Be(baseScore,
                because: "bases beyond offset +5 are outside the PWM window and cannot change the donor score");
    }

    #endregion
}
