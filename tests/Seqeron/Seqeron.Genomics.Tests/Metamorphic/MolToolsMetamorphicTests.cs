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
/// Units: CRISPR-PAM-001 (PAM-site finding), CRISPR-GUIDE-001 (guide-RNA design),
///        CRISPR-OFF-001 (off-target scoring).
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

    #region CRISPR-GUIDE-001 — guide RNA design

    // ─────────────────────────────────────────────────────────────────────────
    // Unit: CRISPR-GUIDE-001 — CRISPR guide-RNA design (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 19.
    //
    // API under test (CrisprDesigner.cs):
    //   DesignGuideRnas(DnaSequence, regionStart, regionEnd, systemType, GuideRnaParameters?)
    //     1. FindPamSitesCore over the WHOLE sequence,
    //     2. keep sites whose cut-site lies in [regionStart, regionEnd] (IsInRegion),
    //     3. EvaluateGuideRna(pamSite.TargetSequence) → a heuristic Score in [0,100],
    //     4. yield candidates with Score >= effectiveParams.MinScore.
    //   A designed guide's identity is (Sequence, Position, IsForwardStrand) where
    //   Sequence == pamSite.TargetSequence (the 20-nt protospacer adjacent to the PAM),
    //   Position == pamSite.TargetStart, IsForwardStrand == pamSite.IsForwardStrand.
    //
    // Relations (derived from the design definition, NOT from observed output):
    //   • SUB — every designed guide corresponds to a real PAM site: its protospacer
    //           is the TargetSequence of an actual site reported by FindPamSites, and
    //           the guide set is ⊆ the candidate PAM-adjacent windows in the region.
    //           Filtering (MinScore + cut-site-in-region) only REMOVES, never invents;
    //           hence guide count ≤ in-region PAM-site count ≤ total PAM-site count.
    //   • MON — MinScore is the real filtering threshold. A guide's Score is fixed by
    //           its protospacer sequence alone (EvaluateGuideRna is composition-only),
    //           so raising MinScore can only DROP guides whose score < threshold:
    //           the higher-threshold guide set is a SUBSET of the lower one and the
    //           count is non-increasing along a threshold chain.
    //   • INV — Score, Sequence, Position and IsForwardStrand of a guide depend ONLY
    //           on the local PAM-adjacent window (EvaluateGuideRna never scans the
    //           genome; the protospacer/PAM are read from a fixed local slice). A
    //           sequence edit FAR downstream — introducing no new PAM motif and whose
    //           cut-site is outside the design region — therefore leaves the in-region
    //           guide set IDENTICAL (same sequences, positions, strands, scores).
    //           Scope encoded: the LOCAL guide IDENTITY (Sequence/Position/Strand) and
    //           its Score, for a design region pinned away from the edited tail.
    //
    // Sequence construction reuses the CRISPR-PAM-001 helpers above:
    //   - Pad / NonPamRegion(): A/T-only, GG-free and CC-free, so they create no
    //     SpCas9 PAM (forward NGG = GG; reverse NGG = forward CC) on either strand and
    //     none across any junction (they start and end in A/T).
    //   - "AGG"/"CGG"/"TGG" embed forward NGG PAMs at known offsets.
    //   GG/CC-free 20-nt guide windows keep each protospacer unambiguous.
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// A 20-nt protospacer with no GG/CC dinucleotide and a comfortable GC content,
    /// so when followed by an NGG PAM it forms a single clean, high-scoring candidate.
    /// "ACGTACGTACGTACGTACGT" = 50% GC, no poly-T, seed GC in range ⇒ score 100.
    /// </summary>
    private const string CleanGuide20 = "ACGTACGTACGTACGTACGT";

    /// <summary>
    /// Sequences carrying one or more forward-strand SpCas9 PAMs, each preceded by a
    /// 20-nt protospacer, plus fixed-seed random sequences. Used as design templates.
    /// The cut site for an SpCas9 forward PAM at index p is p−3, i.e. inside the
    /// protospacer, so passing the whole sequence as the region captures every PAM.
    /// </summary>
    private static IEnumerable<string> GuideTemplates()
    {
        // Single clean PAM: 20-nt guide + AGG.
        yield return CleanGuide20 + "AGG";
        // Three tandem clean PAMs (distinct N in NGG).
        yield return CleanGuide20 + "AGG" + CleanGuide20 + "CGG" + CleanGuide20 + "TGG";
        // Mixed-quality protospacers so MinScore actually discriminates:
        //   high-GC guide (penalised) then a clean guide.
        yield return "GCGCGCGCGCGCGCGCGCGC" + "AGG" + CleanGuide20 + "CGG";
        // Fixed-seed random templates with appended PAMs (relations hold for any input).
        yield return RandomDna(40) + "AGG" + RandomDna(20) + "TGG";
        yield return RandomDna(60) + "CGG" + RandomDna(30) + "AGG" + RandomDna(20) + "TGG";
        yield return RandomDna(100) + "AGG";
    }

    private static List<GuideRnaCandidate> Design(string sequence, GuideRnaParameters? parameters = null)
        => CrisprDesigner.DesignGuideRnas(
                new DnaSequence(sequence), 0, sequence.Length - 1, CrisprSystemType.SpCas9, parameters)
            .ToList();

    /// <summary>Run-order-independent identity of a designed guide.</summary>
    private static (int Position, string Sequence, bool Fwd) GuideId(GuideRnaCandidate g)
        => (g.Position, g.Sequence, g.IsForwardStrand);

    private static HashSet<(int, string, bool)> GuideIdSet(IEnumerable<GuideRnaCandidate> guides)
        => guides.Select(GuideId).ToHashSet();

    #region SUB — every guide comes from a real PAM site; guides ⊆ candidate windows

    [Test]
    [Description("SUB: each designed guide's protospacer is the TargetSequence of an actual PAM site found by FindPamSites.")]
    public void DesignGuideRnas_EveryGuide_MapsToARealPamSite()
    {
        foreach (var template in GuideTemplates())
        {
            var pamTargets = CrisprDesigner.FindPamSites(template, CrisprSystemType.SpCas9)
                .Select(p => (p.TargetSequence, p.TargetStart, p.IsForwardStrand))
                .ToHashSet();

            foreach (var guide in Design(template))
            {
                pamTargets.Should().Contain((guide.Sequence, guide.Position, guide.IsForwardStrand),
                    because: "a designed guide is never invented — its protospacer/position/strand are copied " +
                             "verbatim from a PAM site that FindPamSites reports for the same sequence");
            }
        }
    }

    [Test]
    [Description("SUB: guide count never exceeds the number of PAM sites — filtering only removes candidates.")]
    public void DesignGuideRnas_GuideCount_DoesNotExceedPamSiteCount()
    {
        foreach (var template in GuideTemplates())
        {
            int pamCount = CrisprDesigner.FindPamSites(template, CrisprSystemType.SpCas9).Count();
            int guideCount = Design(template).Count;

            guideCount.Should().BeLessThanOrEqualTo(pamCount,
                because: "DesignGuideRnas draws guides from PAM sites and then DROPS some (cut-site-in-region " +
                         "and MinScore filters); it can never produce more guides than there are PAM sites");
        }
    }

    [Test]
    [Description("SUB: dropping the MinScore floor to 0 makes guides ⊆ in-region PAM-adjacent windows (one guide per kept site).")]
    public void DesignGuideRnas_NoScoreFloor_GuidesSubsetOfInRegionPamWindows()
    {
        var permissive = GuideRnaParameters.Default with { MinScore = 0 };

        foreach (var template in GuideTemplates())
        {
            var allTargets = CrisprDesigner.FindPamSites(template, CrisprSystemType.SpCas9)
                .Select(p => (p.TargetStart, p.TargetSequence, p.IsForwardStrand))
                .ToHashSet();

            var guideIds = GuideIdSet(Design(template, permissive));

            allTargets.IsSupersetOf(guideIds).Should().BeTrue(
                because: "with no score floor every kept guide still maps onto a distinct PAM-adjacent candidate " +
                         "window — the guide set is a subset of the candidate windows, never a superset");
        }
    }

    #endregion

    #region MON — stricter scoring (higher MinScore) → subset of guides, non-increasing count

    [Test]
    [Description("MON: along an increasing MinScore chain the guide set shrinks monotonically (each stricter set ⊆ the looser).")]
    public void DesignGuideRnas_RaisingMinScore_YieldsSubset_CountNonIncreasing()
    {
        // A protospacer's Score is fixed by its sequence, so MinScore is a pure cutoff.
        double[] thresholdChain = { 0, 25, 50, 75, 90, 100, 101 };

        foreach (var template in GuideTemplates())
        {
            HashSet<(int, string, bool)>? previousSet = null;
            int previousCount = int.MaxValue;

            foreach (var minScore in thresholdChain)
            {
                var current = Design(template, GuideRnaParameters.Default with { MinScore = minScore });
                var currentSet = GuideIdSet(current);

                if (previousSet is not null)
                {
                    previousSet.IsSupersetOf(currentSet).Should().BeTrue(
                        because: $"raising MinScore to {minScore} can only DROP guides whose fixed score falls below it, " +
                                 "so the stricter guide set is a subset of the looser one");
                    current.Count.Should().BeLessThanOrEqualTo(previousCount,
                        because: $"a stricter MinScore ({minScore}) removes-or-keeps each candidate, never adds one — count is non-increasing");
                }

                previousSet = currentSet;
                previousCount = current.Count;
            }
        }
    }

    [Test]
    [Description("MON: a MinScore strictly above every candidate's score yields no guides; a floor of 0 yields the most.")]
    public void DesignGuideRnas_ThresholdAboveAllScores_Empties_WhileZeroFloorIsMaximal()
    {
        foreach (var template in GuideTemplates())
        {
            int floorZero = Design(template, GuideRnaParameters.Default with { MinScore = 0 }).Count;
            int above100 = Design(template, GuideRnaParameters.Default with { MinScore = 101 }).Count;

            above100.Should().Be(0,
                because: "the heuristic score is clamped to [0,100], so a MinScore of 101 excludes every candidate");
            floorZero.Should().BeGreaterThanOrEqualTo(above100,
                because: "the most permissive floor keeps a superset of any stricter floor's guides");
        }
    }

    #endregion

    #region INV — a distant downstream change preserves each guide's local identity & score

    [Test]
    [Description("INV: editing a far-downstream non-PAM tail, with the design region pinned upstream, leaves the in-region guide set identical.")]
    public void DesignGuideRnas_DistantDownstreamChange_PreservesInRegionGuides()
    {
        // Build: [protospacer + PAM core] anchored at the 5' end, then a long A/T spacer,
        // then a downstream tail we will MUTATE. The design region covers only the core,
        // so the tail is outside it. The A/T spacer keeps any tail PAM ≥ 20 nt away from
        // the core, and the tail itself never reaches back into the core protospacer.
        const string core = CleanGuide20 + "AGG"; // one clean forward PAM, cut-site inside the protospacer
        string spacer = NonPamRegion(40);          // ≥ guide length of A/T padding, no PAM
        int regionEnd = core.Length - 1;           // design region = the core only

        DnaSequence Build(string tail) => new(core + spacer + tail);

        var baseGuides = CrisprDesigner
            .DesignGuideRnas(Build(NonPamRegion(30)), 0, regionEnd, CrisprSystemType.SpCas9)
            .ToList();
        var baseSet = GuideIdSet(baseGuides);
        baseSet.Should().NotBeEmpty(
            because: "the clean 20-nt protospacer + AGG yields a high-scoring guide inside the design region");

        // Several distinct downstream tails, including PAM-rich and fixed-seed random ones.
        foreach (var tail in new[]
                 {
                     NonPamRegion(60),                       // different length / content, still no PAM
                     "AGGTGGCCACCT" + NonPamRegion(20),      // PAM-rich tail (outside the region)
                     "GGGGGGGGGGGG" + NonPamRegion(20),
                     RandomDna(50),
                     RandomDna(120),
                 })
        {
            var changedGuides = CrisprDesigner
                .DesignGuideRnas(Build(tail), 0, regionEnd, CrisprSystemType.SpCas9)
                .ToList();

            GuideIdSet(changedGuides).SetEquals(baseSet).Should().BeTrue(
                because: "the in-region guides' identities depend only on the upstream local windows, " +
                         "which the distant-tail edit does not touch — same sequences, positions, strands");

            // And the local INV is exact down to the heuristic Score (composition-only metric).
            foreach (var g in changedGuides)
            {
                var twin = baseGuides.Single(b => GuideId(b) == GuideId(g));
                g.Score.Should().Be(twin.Score,
                    because: "EvaluateGuideRna scores the protospacer alone and never inspects the genome, " +
                             "so a distant change cannot move a guide's score");
            }
        }
    }

    [Test]
    [Description("INV: a single guide's identity & score are byte-for-byte stable as the downstream context is rewritten.")]
    public void DesignGuideRnas_SingleGuide_IdentityAndScoreStableUnderDownstreamRewrite()
    {
        const string core = CleanGuide20 + "AGG";
        string spacer = NonPamRegion(40);
        int regionEnd = core.Length - 1;

        var reference = CrisprDesigner
            .DesignGuideRnas(new DnaSequence(core + spacer + NonPamRegion(10)), 0, regionEnd, CrisprSystemType.SpCas9)
            .Single();

        reference.Sequence.Should().Be(CleanGuide20,
            because: "the protospacer is the 20 nt immediately 5' of the AGG PAM");
        reference.IsForwardStrand.Should().BeTrue();

        foreach (var tail in new[] { NonPamRegion(80), "TGGAGGCGG" + RandomDna(40), RandomDna(70) })
        {
            var again = CrisprDesigner
                .DesignGuideRnas(new DnaSequence(core + spacer + tail), 0, regionEnd, CrisprSystemType.SpCas9)
                .Single();

            again.Sequence.Should().Be(reference.Sequence,
                because: "the protospacer is a fixed local slice, independent of the downstream tail");
            again.Position.Should().Be(reference.Position,
                because: "the guide's position is anchored to its local PAM, untouched by the distant edit");
            again.IsForwardStrand.Should().Be(reference.IsForwardStrand);
            again.Score.Should().Be(reference.Score,
                because: "the composition-only score reads only the protospacer, so it is invariant to distant context");
        }
    }

    #endregion

    #endregion

    #region CRISPR-OFF-001 — off-target scoring

    // ─────────────────────────────────────────────────────────────────────────
    // Unit: CRISPR-OFF-001 — CRISPR off-target scoring (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 20.
    //
    // API under test (CrisprDesigner.cs):
    //   FindOffTargets(guide, genome, maxMismatches, systemType)
    //     → OffTargetSite.OffTargetScore: a POSITION-WEIGHTED PENALTY (higher = worse;
    //       grows as mismatches accumulate). Per the spec, seed mismatch = 5, non-seed = 2.
    //   CalculateSpecificityScore(guide, genome, systemType)
    //     → an off-target SCORE in [0,100] (higher = better / more specific):
    //       100 − Σ off-target penalties, clamped to [0,100]; 100 when no off-target.
    //
    // The MR phrase "off-target score" denotes the higher-is-better SCORE for which a
    // documented MAXIMUM exists and which DECREASES with off-target burden — i.e. the
    // specificity score (max = 100; the spec's "No off-targets found ⇒ 100"). The
    // complementary penalty (OffTargetScore) is the dual: it is NON-DECREASING as
    // mismatches accumulate. Both directions are asserted.
    //
    // Relation DIRECTIONS are derived from the cited source/spec, NOT from code output:
    //   Source: docs/algorithms/MolTools/Off_Target_Analysis.md
    //     §2.1 — "PAM-proximal mismatches are especially important for specificity"
    //            (Hsu 2013, Fu 2013): a seed mismatch costs MORE specificity than a
    //            PAM-distal one.
    //     §2.2 / §4.2 — seed region = the 12 PAM-proximal positions (SpCas9: last 12 =
    //            positions 8..19); penalty(seed) > penalty(non-seed); the score is
    //            max(0, 100 − Σ penalties).
    //     §6.1 — "No off-targets found ⇒ specificity 100" (the MAX endpoint).
    //
    //   • COMP (0 mismatches → max score): a candidate identical to the guide is the
    //     on-target (0 mismatches), which is NOT an off-target (INV-01), so it adds no
    //     penalty; with no off-target the specificity score is the documented MAXIMUM 100.
    //   • MON (more mismatches → lower score): along a chain that adds mismatches at
    //     fixed positions, the off-target PENALTY is non-decreasing and the specificity
    //     SCORE is non-increasing; strictly so when the added position carries a
    //     non-zero penalty (every position does: seed=5, non-seed=2 > 0).
    //   • MON (seed mismatch penalized more): a single mismatch placed in the
    //     PAM-PROXIMAL seed region yields a LOWER specificity score (HIGHER penalty)
    //     than the same single mismatch placed at a PAM-distal position — derived from
    //     the source's "PAM-proximal mismatches matter more" structure, not from W values.
    //
    // Construction (SpCas9: NGG PAM 3' of target, 20-nt guide, seed = positions 8..19):
    //   A genome of the form  TARGET(20) + "AGG"  has exactly one forward PAM at index
    //   20, whose protospacer is TARGET. We start from a clean 20-nt guide with no GG/CC
    //   dinucleotide (so the only PAM is the appended AGG and the target window is
    //   unambiguous), then introduce mismatches at chosen positions by an A↔T / C↔G flip
    //   (transversion within the same {weak|strong} pair) so a flip never creates a new
    //   GG/CC PAM motif inside the window. Several fixed-seed random guides are included.
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>SpCas9 seed region = the 12 PAM-proximal positions of a 20-nt guide: indices 8..19.</summary>
    private const int SeedStartSpCas9 = 8;

    /// <summary>
    /// A 20-nt guide with no GG and no CC dinucleotide, so appending "AGG" yields a
    /// single unambiguous forward PAM and a single protospacer window.
    /// </summary>
    private const string OffGuide20 = "ACGTACGTACGTACGTACGT";

    /// <summary>
    /// Flip a base to a guaranteed mismatch WITHOUT creating a new G/C run: A↔T and
    /// C↔G keep the base within its weak/strong pair, so neither a GG nor a CC
    /// dinucleotide can be introduced by the flip.
    /// </summary>
    private static char FlipBase(char b) => b switch
    {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => throw new ArgumentException($"Unexpected base '{b}'.", nameof(b)),
    };

    /// <summary>
    /// Returns a copy of <paramref name="target"/> with a mismatch (vs itself) introduced
    /// at each position in <paramref name="positions"/> via <see cref="FlipBase"/>.
    /// </summary>
    private static string MutateAt(string target, params int[] positions)
    {
        var chars = target.ToCharArray();
        foreach (int p in positions)
            chars[p] = FlipBase(chars[p]);
        return new string(chars);
    }

    /// <summary>Builds the SpCas9 genome TARGET(20)+"AGG" whose sole forward protospacer is <paramref name="target"/>.</summary>
    private static DnaSequence OffGenome(string target) => new(target + "AGG");

    /// <summary>Off-target specificity score (higher = better; 100 = no off-target) for guide vs genome TARGET+AGG.</summary>
    private static double Specificity(string guide, string target)
        => CrisprDesigner.CalculateSpecificityScore(guide, OffGenome(target), CrisprSystemType.SpCas9);

    /// <summary>The single off-target penalty (OffTargetScore) for guide against TARGET+AGG, or 0 if no off-target (exact match).</summary>
    private static double OffTargetPenalty(string guide, string target)
        => CrisprDesigner.FindOffTargets(guide, OffGenome(target), 5, CrisprSystemType.SpCas9)
            .Select(ot => ot.OffTargetScore)
            .DefaultIfEmpty(0)
            .Sum();

    /// <summary>
    /// 20-nt guides with no GG/CC dinucleotide (clean single-PAM windows) plus fixed-seed
    /// random guides — relations must hold for arbitrary input. Random guides are sanitized
    /// to remove GG/CC so the appended AGG stays the only PAM.
    /// </summary>
    private static IEnumerable<string> OffGuides()
    {
        yield return OffGuide20;
        yield return "ATATATATATATATATATAT";          // weak-only
        yield return "ACATACATACATACATACAT";          // mixed, GG/CC-free
        yield return "TGTATGTATGTATGTATGTA";          // GG/CC-free
        for (int i = 0; i < 3; i++)
            yield return SanitizeNoGgCc(RandomDna(20)); // fixed-seed random, sanitized
    }

    /// <summary>Removes GG and CC dinucleotides by flipping the second base of any such pair (keeps it a valid 20-nt guide).</summary>
    private static string SanitizeNoGgCc(string s)
    {
        var c = s.ToCharArray();
        for (int i = 1; i < c.Length; i++)
            if ((c[i] == 'G' && c[i - 1] == 'G') || (c[i] == 'C' && c[i - 1] == 'C'))
                c[i] = c[i] == 'G' ? 'A' : 'T';
        return new string(c);
    }

    #region COMP — 0 mismatches → maximum off-target score (specificity = 100)

    [Test]
    [Description("COMP: a candidate identical to the guide is the on-target (0 mismatches), not an off-target, so the off-target score is the documented MAXIMUM (100).")]
    public void CalculateSpecificityScore_GuideIdenticalToTarget_IsMaximum100()
    {
        foreach (var guide in OffGuides())
        {
            // The genome's protospacer equals the guide → 0 mismatches → excluded as off-target.
            double score = Specificity(guide, target: guide);

            score.Should().Be(100.0,
                because: "a 0-mismatch candidate is the on-target (INV-01: not an off-target), so no penalty " +
                         "accrues and the specificity score is the documented maximum endpoint of 100");

            // Dual: with no off-target the accumulated penalty is exactly 0 (the penalty's minimum).
            OffTargetPenalty(guide, target: guide).Should().Be(0.0,
                because: "an exact (0-mismatch) match is filtered out by the mismatches>0 rule, so it contributes no off-target penalty");
        }
    }

    #endregion

    #region MON — more mismatches → lower off-target score (penalty non-decreasing, specificity non-increasing)

    [Test]
    [Description("MON: along a chain that adds mismatches at fixed positions, the off-target score is non-increasing — and strictly decreasing, since every added position carries a positive penalty.")]
    public void CalculateSpecificityScore_AddingMismatches_ScoreStrictlyDecreasing()
    {
        // A fixed chain of distinct positions; each step ADDS one more mismatch to the
        // off-target window. The set is the same for every guide, so the relation is
        // a property of the SCORING MODEL, not of any one input.
        //
        // The chain is kept to 4 positions: CalculateSpecificityScore evaluates off-targets
        // with a FIXED mismatch cap of 4 (spec §3.3 / §5.2). A candidate with >4 mismatches
        // leaves the off-target search window entirely (INV-02) — it is no longer an
        // off-target, so its penalty vanishes and the score returns to the no-off-target
        // maximum (100). The strict monotone-decrease relation is therefore stated where it
        // is theoretically guaranteed: while the candidate remains an in-window off-target
        // (1..4 mismatches). (The dual penalty, probed via FindOffTargets with cap 5, is
        // separately covered up to its own bound.)
        int[] chainPositions = { 0, 7, 12, 19 };

        foreach (var guide in OffGuides())
        {
            double previousScore = 100.0;   // 0 mismatches ⇒ on-target ⇒ specificity 100
            double previousPenalty = 0.0;
            var applied = new List<int>();

            foreach (int pos in chainPositions)
            {
                applied.Add(pos);
                string target = MutateAt(guide, applied.ToArray());

                double score = Specificity(guide, target);
                double penalty = OffTargetPenalty(guide, target);

                score.Should().BeLessThan(previousScore,
                    because: $"adding a mismatch at position {pos} adds a positive penalty (seed=5 or non-seed=2), " +
                             "so the specificity score must strictly drop relative to the previous step in the chain");
                penalty.Should().BeGreaterThan(previousPenalty,
                    because: $"the dual off-target penalty is non-decreasing in mismatches and strictly increases at position {pos} (penalty > 0)");

                previousScore = score;
                previousPenalty = penalty;
            }
        }
    }

    [Test]
    [Description("MON: a candidate with a SUPERSET of another's mismatch positions never scores higher (penalty never lower) — non-strict monotonicity over arbitrary position sets.")]
    public void CalculateSpecificityScore_SupersetOfMismatches_NeverScoresHigher()
    {
        // Pairs (subset ⊂ superset) of mismatch-position sets, fixed across guides.
        var cases = new[]
        {
            (Subset: new[] { 5 },           Superset: new[] { 5, 12 }),
            (Subset: new[] { 0, 1 },        Superset: new[] { 0, 1, 18 }),
            (Subset: new[] { 9, 14 },       Superset: new[] { 2, 9, 14, 17 }),
            (Subset: new[] { 19 },          Superset: new[] { 4, 8, 19 }),
        };

        foreach (var guide in OffGuides())
        {
            foreach (var (subset, superset) in cases)
            {
                double subScore = Specificity(guide, MutateAt(guide, subset));
                double supScore = Specificity(guide, MutateAt(guide, superset));

                supScore.Should().BeLessThanOrEqualTo(subScore,
                    because: "a superset of mismatch positions accumulates at least the subset's penalty plus more, " +
                             "so the specificity score is non-increasing as the mismatch set grows");

                double subPenalty = OffTargetPenalty(guide, MutateAt(guide, subset));
                double supPenalty = OffTargetPenalty(guide, MutateAt(guide, superset));
                supPenalty.Should().BeGreaterThanOrEqualTo(subPenalty,
                    because: "adding positions to a mismatch set can only add penalty, never remove it");
            }
        }
    }

    #endregion

    #region MON — a seed-region mismatch is penalized more than a PAM-distal one

    [Test]
    [Description("MON: a single mismatch in the PAM-proximal SEED region yields a LOWER off-target score (higher penalty) than the same single mismatch at a PAM-distal position.")]
    public void CalculateSpecificityScore_SeedMismatch_LowerScoreThanDistalMismatch()
    {
        // Distal positions: 0..7 (outside seed 8..19). Seed positions: 8..19.
        // Derived from the source: PAM-proximal (seed) mismatches matter MORE for
        // specificity ⇒ score(seed-mismatch) < score(distal-mismatch). We pair one
        // distal and one seed position so each candidate carries EXACTLY ONE mismatch.
        var pairs = new[]
        {
            (Distal: 0,  Seed: 8),
            (Distal: 3,  Seed: 12),
            (Distal: 7,  Seed: 19),
            (Distal: 2,  Seed: 15),
        };

        foreach (var guide in OffGuides())
        {
            foreach (var (distal, seed) in pairs)
            {
                distal.Should().BeLessThan(SeedStartSpCas9, because: "the distal position must lie OUTSIDE the seed (indices 8..19)");
                seed.Should().BeGreaterThanOrEqualTo(SeedStartSpCas9, because: "the seed position must lie INSIDE the seed (indices 8..19)");

                double distalScore = Specificity(guide, MutateAt(guide, distal));
                double seedScore = Specificity(guide, MutateAt(guide, seed));

                seedScore.Should().BeLessThan(distalScore,
                    because: $"a single seed mismatch (pos {seed}) is weighted more heavily than a single distal mismatch " +
                             $"(pos {distal}) per the source (PAM-proximal mismatches matter more), so it lowers the score more");

                // Dual on the penalty: the seed mismatch's penalty strictly exceeds the distal one's.
                double distalPenalty = OffTargetPenalty(guide, MutateAt(guide, distal));
                double seedPenalty = OffTargetPenalty(guide, MutateAt(guide, seed));
                seedPenalty.Should().BeGreaterThan(distalPenalty,
                    because: "the seed-region mismatch penalty exceeds the PAM-distal mismatch penalty (seed weighted heavier)");
            }
        }
    }

    [Test]
    [Description("MON: ordering is consistent across the whole guide — moving a single mismatch from any distal position into the seed never raises the off-target score.")]
    public void CalculateSpecificityScore_MovingMismatchIntoSeed_NeverRaisesScore()
    {
        foreach (var guide in OffGuides())
        {
            // Every distal position (0..7) paired against every seed position (8..19):
            // a single-mismatch in the seed is never better-scored than one PAM-distal.
            for (int distal = 0; distal < SeedStartSpCas9; distal++)
            {
                double distalScore = Specificity(guide, MutateAt(guide, distal));
                for (int seed = SeedStartSpCas9; seed < 20; seed++)
                {
                    double seedScore = Specificity(guide, MutateAt(guide, seed));
                    seedScore.Should().BeLessThanOrEqualTo(distalScore,
                        because: $"relocating a single mismatch from distal {distal} into seed {seed} weights it more heavily, " +
                                 "so the specificity score cannot increase");
                }
            }
        }
    }

    #endregion

    #endregion
}
