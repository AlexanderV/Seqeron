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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SPLICE-ACCEPTOR-001 — 3' acceptor splice-site prediction (Splicing).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 78.
    //
    // API under test (SpliceSitePredictor.FindAcceptorSites):
    //   A canonical acceptor requires the invariant AG dinucleotide; its score combines a
    //   polypyrimidine-tract count over positions i−15..i−4 with a consensus PWM over the AG
    //   context (i−13..i+2). Bases outside [i−15, i+2] do not enter the score.
    //
    // Relations (derived from the scoring window, NOT from output):
    //   • COMP (non-AG ⇒ no acceptor): without an AG dinucleotide there is no canonical acceptor.
    //   • MON  (closer to consensus ⇒ higher score): a strong polypyrimidine tract and the
    //          consensus C at −3 raise the acceptor score.
    //   • INV  (far-upstream change ⇒ same score): edits before position i−15 lie outside the
    //          scoring window and leave the acceptor's score unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    /// <summary>Score of the acceptor site reported at index <paramref name="reportedPos"/> (the G of AG).</summary>
    private static double AcceptorScoreAt(string sequence, int reportedPos) =>
        SpliceSitePredictor.FindAcceptorSites(sequence, minScore: 0.0).Single(s => s.Position == reportedPos).Score;

    #region SPLICE-ACCEPTOR-001 COMP — no AG dinucleotide ⇒ no canonical acceptor

    [Test]
    [Description("COMP: a canonical acceptor requires the invariant AG, so a sequence with no AG dinucleotide yields no acceptor sites.")]
    public void FindAcceptorSites_NoAg_ReturnsNothing()
    {
        foreach (var seq in new[] { "CCCCCCCCCCCCCCCCCCCC", "CACACACACACACACACACA", "UUUUUUUUUUUUUUUUUUUU" })
            SpliceSitePredictor.FindAcceptorSites(seq, minScore: 0.0)
                .Should().BeEmpty(because: $"'{seq}' contains no AG dinucleotide to anchor a canonical acceptor");
    }

    #endregion

    #region SPLICE-ACCEPTOR-001 MON — closer to the consensus raises the score

    [Test]
    [Description("MON: a strong polypyrimidine tract and the consensus C at position −3 raise the acceptor score; weakening them lowers it.")]
    public void FindAcceptorSites_CloserToConsensus_HigherScore()
    {
        // AG anchored at indices 16/17 (reported Position 17). 15-nt tract precedes it.
        const string strongPpt = "UCUCUCUCUCUCUCU"; // pyrimidine-rich
        const string weakPpt = "AGAGAGAGAGAGAGA";   // purine-rich (no C/U)

        double best   = AcceptorScoreAt(strongPpt + "CAGG" + "AAAA", 17); // strong PPT + consensus C at −3
        double midC   = AcceptorScoreAt(strongPpt + "AAGG" + "AAAA", 17); // strong PPT, −3 broken (A)
        double weakest = AcceptorScoreAt(weakPpt + "AAGG" + "AAAA", 17);   // weak PPT, −3 broken

        best.Should().BeGreaterThan(midC, because: "the consensus C at position −3 scores higher than a non-consensus A");
        midC.Should().BeGreaterThan(weakest, because: "a pyrimidine-rich tract scores higher than a purine-rich one");
    }

    #endregion

    #region SPLICE-ACCEPTOR-001 INV — far-upstream edits don't change the score

    [Test]
    [Description("INV: the acceptor score uses only positions i−15..i+2, so editing bases before i−15 (far upstream) leaves the acceptor's score unchanged.")]
    public void FindAcceptorSites_FarUpstreamEdits_DoNotChangeScore()
    {
        // 10-nt pad, then a 15-nt tract and CAGG. AG at indices 26/27 (Position 27);
        // its scoring window starts at index 11, so the pad (indices 0..9) is far upstream.
        const string rest = "UCUCUCUCUCUCUCU" + "CAGG" + "AAAA";
        double baseScore = AcceptorScoreAt("AAAAAAAAAA" + rest, 27);

        foreach (var pad in new[] { "CCCCCCCCCC", "GGGGGGGGGG", "UUUUUUUUUU", "GAGAGAGAGA" })
            AcceptorScoreAt(pad + rest, 27).Should().Be(baseScore,
                because: "bases before position i−15 are outside the scoring window and cannot change the acceptor score");
    }

    #endregion
}
