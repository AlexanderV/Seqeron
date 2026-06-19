using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the MolTools area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs
/// of multiple runs under an input transformation, with no hardcoded oracle. The
/// relations are derived from the PAM-site-detection *definition*, not from the
/// current implementation's output: a PAM site exists wherever the strand-local
/// PAM motif occurs and the guide-length target fits within sequence bounds.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CRISPR-PAM-001 — CRISPR PAM-site finding (MolTools).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 18.
/// Relations:
///   • MON   — appending bases cannot remove existing PAM sites ⇒ count is
///             non-decreasing (the original region's sites persist).
///   • INV   — appending a region that creates NO new PAM motif (and no
///             junction-spanning motif) leaves the PAM-site count unchanged and
///             preserves the original sites exactly.
///   • SHIFT — prepending a non-PAM-creating flank F shifts every PAM-site
///             forward-strand Position by |F| and preserves the site set.
///
/// Source (semantics, pinned from spec + source, NOT from observed output):
///   docs/algorithms/Molecular_Tools/PAM_Site_Detection.md §2.2, §4 (both strands
///   scanned; forward NGG, reverse strand = revComp NGG = forward-strand CCN);
///   src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs
///   (FindPamSites / FindPamSitesCore). PamSite.Position is ALWAYS a forward-strand
///   coordinate (record doc-comment, CrisprDesigner.cs); 0-based.
///
/// Why the relations hold (SpCas9, PAM = NGG, 20-nt guide, PAM 3' of target):
///   A forward-strand PAM is an "NGG" (a GG dinucleotide) at index i, reported only
///   when its 20-nt target fits UPSTREAM (i ≥ 20). A reverse-strand PAM is an "NGG"
///   on the reverse complement, i.e. a "CCN" (a CC dinucleotide) on the forward
///   strand; its 20-nt target lies on the revComp, which maps to forward positions
///   DOWNSTREAM of the CC (forward [p+3, p+22] for a CC at forward index p), so a
///   reverse site needs room toward the 3' end.
///   • A "non-PAM-creating" flank/append must introduce neither a GG nor a CC
///     dinucleotide — internally OR across any junction. We use the alternating A/T
///     region "ATAT…" (no G, no C at all), which starts and ends in A/T, so every
///     junction with it is "<base>A" / "T<base>" → never GG or CC.
///   • MON (append, any region): appending only ADDS bases downstream. No upstream
///     motif is removed and no forward target (upstream) is shortened, so every
///     original site persists ⇒ count is non-decreasing. (This holds for ANY append,
///     even a PAM-rich one — hence MON is tested on arbitrary suffixes.)
///   • INV / SHIFT require EXACT preservation, so we additionally guard the relevant
///     end: appending downstream could only *add* a reverse (CC) site that was
///     previously starved of 3' room, and prepending upstream could only *add* a
///     forward (GG) site previously starved of 5' room. Padding every test body with
///     ≥22 A/T bases on BOTH ends means every in-body motif already has its full
///     20-nt target room, so neither prepend nor append can create a new in-bounds
///     site ⇒ the site set is preserved exactly.
///   • SHIFT positions: Position is reported on the FORWARD strand for both strands,
///     so a length-|F| prepend advances every site's Position by exactly |F|.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class MolToolsMetamorphicTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed so random inputs are reproducible.</summary>
    private static readonly Random Rng = new(20260619);

    /// <summary>
    /// A flank/append region guaranteed to create NO PAM motif on either strand:
    /// it contains only A and T, so it has neither a GG (forward NGG) nor a CC
    /// (reverse-strand CCN) dinucleotide, and—starting and ending in A/T—it forms
    /// none across any junction either.
    /// </summary>
    private static string NonPamRegion(int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (i % 2 == 0) ? 'A' : 'T';
        return new string(chars);
    }

    /// <summary>Generates a random DNA string of the given length over {A,C,G,T}.</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// SpCas9 PAM sites as a comparable set of (Position, PamSequence, IsForwardStrand)
    /// tuples — the strand-aware identity of a site, independent of run order.
    /// </summary>
    private static HashSet<(int Position, string Pam, bool Fwd)> SiteSet(string sequence)
        => CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9)
            .Select(s => (s.Position, s.PamSequence, s.IsForwardStrand))
            .ToHashSet();

    private static int SiteCount(string sequence)
        => CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).Count();

    /// <summary>
    /// Sequences with KNOWN PAM placement plus a few fixed-seed random ones.
    /// Bodies embed explicit "…AGG…" / "…CCA…" motifs at predictable offsets while
    /// keeping the 20-nt guide windows free of GG/CC so the placement is unambiguous.
    /// </summary>
    private static IEnumerable<string> SampleSequences()
    {
        // 20-nt GG/CC-free guide, then a forward NGG ("AGG"), then more GG/CC-free body.
        yield return "ATATATATATATATATATAT" + "AGG" + "ATATATATAT";
        // Two overlapping forward PAMs ("AGGTGG") after a 20-nt guide.
        yield return "ATATATATATATATATATAT" + "AGGTGG" + "ATATAT";
        // A forward PAM and a reverse-strand motif ("CCA") downstream with its own room.
        yield return "ATATATATATATATATATAT" + "AGG" + "ATATATATATATATATATAT" + "CCA" + "ATAT";
        // PAM-dense block.
        yield return "ATATATATATATATATATAT" + "AGGCGGTGGGGG" + "ATATATATAT";
        // No PAM at all (no GG, no CC anywhere).
        yield return "ATATATATATATATATATATATATATATATAT";
        // Fixed-seed random sequences (relations must hold for arbitrary input too).
        yield return RandomDna(60);
        yield return RandomDna(120);
        yield return RandomDna(200);
    }

    /// <summary>22-nt A/T pad (&gt; 20-nt guide), with no GG/CC, used to saturate body ends.</summary>
    private const string Pad = "ATATATATATATATATATATAT";

    /// <summary>
    /// Bodies padded with a ≥22-nt A/T region on BOTH ends so every embedded motif
    /// already has its full 20-nt target room within the body. This makes the EXACT
    /// preservation relations (INV append, SHIFT prepend) hold: a downstream append
    /// cannot newly enable a 3'-starved reverse site, and an upstream prepend cannot
    /// newly enable a 5'-starved forward site, because none exist.
    /// </summary>
    private static IEnumerable<string> PaddedSequences()
    {
        yield return Pad + "AGG" + Pad;                       // one forward PAM
        yield return Pad + "AGGTGG" + Pad;                    // two overlapping forward PAMs
        yield return Pad + "CCA" + Pad;                       // one reverse-strand motif (CC)
        yield return Pad + "AGG" + Pad + "CCA" + Pad;         // forward + reverse, both with room
        yield return Pad + "AGGCGGTGG" + Pad + "CCACCT" + Pad; // PAM-dense, both strands
        yield return Pad + Pad;                               // no PAM at all
        yield return Pad + RandomDna(40) + Pad;               // random core, saturated ends
        yield return Pad + RandomDna(80) + Pad;
    }

    #endregion

    #region MON — appending bases cannot remove PAM sites (count non-decreasing)

    [Test]
    [Description("MON: appending ANY region (even one rich in new PAMs) never lowers the PAM-site count.")]
    public void FindPamSites_AppendArbitrarySuffix_CountIsNonDecreasing()
    {
        foreach (var body in SampleSequences())
        {
            int baseCount = SiteCount(body);

            // Append several arbitrary suffixes, including PAM-rich ones.
            foreach (var suffix in new[] { "AGGTGGCCACCT", "GGGGGGGGGG", RandomDna(40), NonPamRegion(30) })
            {
                int extendedCount = SiteCount(body + suffix);
                extendedCount.Should().BeGreaterThanOrEqualTo(baseCount,
                    because: $"appending '{suffix[..Math.Min(6, suffix.Length)]}…' to a length-{body.Length} body " +
                             "can only add room/motifs downstream; sites in the original region persist, so count cannot drop");
            }
        }
    }

    [Test]
    [Description("MON: every PAM site of the original sequence is still present after appending a suffix (set is a subset).")]
    public void FindPamSites_AppendSuffix_OriginalSitesPersist()
    {
        foreach (var body in SampleSequences())
        {
            var baseSites = SiteSet(body);
            var extendedSites = SiteSet(body + NonPamRegion(24));

            extendedSites.IsSupersetOf(baseSites).Should().BeTrue(
                because: "appending downstream removes no upstream motif and shortens no target, " +
                         "so the original site set survives as a subset of the extended set");
        }
    }

    #endregion

    #region INV — appending a non-PAM region preserves the count and the sites exactly

    [Test]
    [Description("INV: appending an A/T-only region (no GG, no CC, no junction motif) leaves the PAM-site count unchanged.")]
    public void FindPamSites_AppendNonPamRegion_CountUnchanged()
    {
        foreach (var body in PaddedSequences())
        {
            int baseCount = SiteCount(body);

            foreach (var len in new[] { 1, 4, 13, 50 })
            {
                int extendedCount = SiteCount(body + NonPamRegion(len));
                extendedCount.Should().Be(baseCount,
                    because: $"an A/T-only append of length {len} introduces neither an NGG (GG) nor a reverse CCN (CC) " +
                             "motif, internally or across the junction, so the PAM-site count is exactly preserved");
            }
        }
    }

    [Test]
    [Description("INV: appending an A/T-only region preserves the EXACT set of original PAM sites (identity, not just count).")]
    public void FindPamSites_AppendNonPamRegion_SiteSetIdentical()
    {
        foreach (var body in PaddedSequences())
        {
            var baseSites = SiteSet(body);
            var extendedSites = SiteSet(body + NonPamRegion(40));

            extendedSites.SetEquals(baseSites).Should().BeTrue(
                because: "with no new motif created the extended sequence yields exactly the original sites — " +
                         "same positions, same PAMs, same strands");
        }
    }

    #endregion

    #region SHIFT — prepending a non-PAM flank shifts every Position by |F|, preserving the set

    [Test]
    [Description("SHIFT: prepending an A/T-only flank F shifts every PAM-site Position by exactly |F| and preserves the count.")]
    public void FindPamSites_PrependNonPamFlank_ShiftsPositionsByFlankLength()
    {
        foreach (var body in PaddedSequences())
        {
            var baseSites = CrisprDesigner.FindPamSites(body, CrisprSystemType.SpCas9).ToList();

            foreach (var flankLen in new[] { 1, 5, 20, 37 })
            {
                string flank = NonPamRegion(flankLen);
                var shifted = CrisprDesigner.FindPamSites(flank + body, CrisprSystemType.SpCas9).ToList();

                shifted.Count.Should().Be(baseSites.Count,
                    because: $"an A/T-only prefix of length {flankLen} creates no new motif and removes none, " +
                             "so the site count is preserved");

                // The shifted set must equal the original set with every Position advanced by |F|.
                var expected = baseSites
                    .Select(s => (s.Position + flankLen, s.PamSequence, s.IsForwardStrand))
                    .ToHashSet();
                var actual = shifted
                    .Select(s => (s.Position, s.PamSequence, s.IsForwardStrand))
                    .ToHashSet();

                actual.SetEquals(expected).Should().BeTrue(
                    because: $"Position is a forward-strand coordinate for BOTH strands, so a length-{flankLen} prepend " +
                             "advances every site's Position by exactly that amount while keeping PAM and strand identical");
            }
        }
    }

    [Test]
    [Description("SHIFT: a known forward AGG at a fixed offset moves by exactly |F| when an A/T-only flank is prepended.")]
    public void FindPamSites_PrependFlank_KnownForwardPam_ShiftsExactly()
    {
        // 20-nt GG/CC-free guide, then forward "AGG" PAM at index 20.
        const string body = "ATATATATATATATATATAT" + "AGG" + "ATATAT";

        var basePam = CrisprDesigner.FindPamSites(body, CrisprSystemType.SpCas9)
            .Single(s => s.IsForwardStrand && s.PamSequence == "AGG");
        basePam.Position.Should().Be(20,
            because: "the AGG PAM follows a 20-nt target, so its 0-based forward position is 20");

        const int flankLen = 12;
        string shifted = NonPamRegion(flankLen) + body;
        var shiftedPam = CrisprDesigner.FindPamSites(shifted, CrisprSystemType.SpCas9)
            .Single(s => s.IsForwardStrand && s.PamSequence == "AGG");

        shiftedPam.Position.Should().Be(20 + flankLen,
            because: $"prepending {flankLen} A/T bases shifts the forward AGG PAM to position {20 + flankLen}");
    }

    #endregion
}
