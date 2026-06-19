using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the ProteinMotif area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-FIND-001 — pattern-based protein motif finding (ProteinMotif).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 82.
///
/// API under test (ProteinMotifFinder.FindMotifByPattern):
///   Reports every (overlapping) occurrence of a regex motif in a protein, as (Start, End,
///   matched substring). A match is a purely local property of the residues it spans.
///
/// Relations (derived from local regex matching, NOT from output):
///   • SHIFT (prepend flank): prepending a non-matching flank of length k shifts every match
///          start by exactly k and preserves the matched substrings.
///   • MON  (broader pattern ⇒ ≥ matches): a pattern that generalises another matches a
///          superset of positions.
///   • INV  (flank change ⇒ same motif): editing residues outside a motif occurrence (without
///          creating/removing matches) leaves that occurrence detected at the same position.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class ProteinMotifMetamorphicTests
{
    private static (int Start, int End, string Seq)[] Matches(string protein, string pattern) =>
        ProteinMotifFinder.FindMotifByPattern(protein, pattern).Select(m => (m.Start, m.End, m.Sequence)).ToArray();

    #region PROTMOTIF-FIND-001 SHIFT — prepending a flank shifts match starts

    [Test]
    [Description("SHIFT: prepending a non-matching flank of length k shifts every match start by exactly k and preserves the matched substrings.")]
    public void FindMotifByPattern_PrependFlank_ShiftsStarts()
    {
        const string protein = "AAACGGCDDD"; // contains one C-x(2)-C motif "CGGC" at index 3
        const string pattern = "C.{2}C";

        var baseMatches = Matches(protein, pattern);
        baseMatches.Should().NotBeEmpty(because: "the protein contains a C-x(2)-C motif");

        foreach (int k in new[] { 1, 3, 5 })
        {
            string flank = new string('W', k); // tryptophans — no cysteine, cannot match or extend the motif
            var shifted = Matches(flank + protein, pattern);

            shifted.Select(m => (m.Start, m.Seq)).Should().Equal(
                baseMatches.Select(m => (m.Start + k, m.Seq)),
                because: $"a non-matching {k}-residue flank shifts every match by {k} while preserving the matched residues");
        }
    }

    #endregion

    #region PROTMOTIF-FIND-001 MON — a broader pattern matches a superset

    [Test]
    [Description("MON: a generalised pattern matches a superset of the positions of the more specific pattern.")]
    public void FindMotifByPattern_BroaderPattern_MatchesSuperset()
    {
        const string protein = "CAACGGCAAS"; // C-x(2)-C twice, plus a C-x(2)-S

        var narrow = Matches(protein, "C.{2}C").Select(m => m.Start).ToHashSet();
        var broad = Matches(protein, "C.{2}[CS]").Select(m => m.Start).ToHashSet();

        broad.IsSupersetOf(narrow).Should().BeTrue(
            because: "C-x(2)-[CS] generalises C-x(2)-C, so it matches at least everywhere the narrower pattern does");
        broad.Count.Should().BeGreaterThan(narrow.Count,
            because: "the broader pattern additionally matches the C-x(2)-S occurrence");
    }

    #endregion

    #region PROTMOTIF-FIND-001 INV — flank edits don't change the detected motif

    [Test]
    [Description("INV: editing residues in the flanks (no cysteine, so no new/removed motif) leaves the motif occurrence detected at the same position with the same residues.")]
    public void FindMotifByPattern_FlankChange_SameMotif()
    {
        const string pattern = "C.{2}C";
        // Motif "CAAC" sits at a fixed index 5 behind a 5-residue prefix; flanks carry no cysteine.
        var baseline = Matches("WWWWW" + "CAAC" + "WWWWW", pattern);
        baseline.Should().ContainSingle().Which.Should().Be((5, 8, "CAAC"));

        foreach (var (prefix, suffix) in new[]
                 {
                     ("AAAAA", "GGGGG"),
                     ("DDDDD", "HHHHH"),
                     ("GGGGG", "AAAAA"),
                 })
        {
            Matches(prefix + "CAAC" + suffix, pattern).Should().Equal(baseline,
                because: "cysteine-free flank edits keep the lengths fixed and create no new motif, so the single occurrence is unchanged");
        }
    }

    #endregion
}
