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
///        CRISPR-OFF-001 (off-target scoring), PRIMER-TM-001 (primer melting temperature),
///        PRIMER-DESIGN-001 (PCR primer design),
///        PRIMER-STRUCT-001 (primer secondary structure: self-dimer / hairpin).
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

    #region PRIMER-TM-001 — melting temperature

    // ─────────────────────────────────────────────────────────────────────────
    // Unit: PRIMER-TM-001 — primer melting temperature (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 21.
    //
    // API under test (PrimerDesigner.cs / ThermoConstants.cs):
    //   PrimerDesigner.CalculateMeltingTemperature(string primer) → double (°C)
    //     1. Uppercases the input; counts AT = #{A,T} and GC = #{G,C}; ignores any
    //        non-ACGT character. validLength = AT + GC.
    //     2. validLength == 0  ⇒ 0.
    //     3. validLength <  14  ⇒ Wallace rule:  Tm = 2·AT + 4·GC
    //                             (ThermoConstants.CalculateWallaceTm).
    //     4. validLength >= 14  ⇒ Marmur-Doty:   Tm = max(0, 64.9 + 41·(GC − 16.4)/N)
    //                             (ThermoConstants.CalculateMarmurDotyTm, clamped ≥ 0).
    //   The branch threshold is ThermoConstants.WallaceMaxLength = 14 (counted bases).
    //
    // Source (formulas pinned from the spec, NOT from observed output):
    //   docs/algorithms/Molecular_Tools/Melting_Temperature.md §2.2, §4.2, §2.4.
    //     Wallace: A/T contributes +2 °C, G/C contributes +4 °C (INV-01, exact constants).
    //     Marmur-Doty: depends only on GC count and counted length N (INV-02).
    //     Longer-primer branch is clamped to ≥ 0 (INV-04).
    //   ThermoConstants.WallaceAtContribution = 2, WallaceGcContribution = 4,
    //   MarmurDotyBase = 64.9, MarmurDotyGcCoefficient = 41, MarmurDotyGcOffset = 16.4.
    //
    // ── Reconciling the checklist MR with the implemented formula ──
    //   The checklist row reads: "MON: add GC → Tm increases; MON: add AT → Tm
    //   decreases; INV: same sequence → same Tm". The literal "add AT → Tm decreases"
    //   is FALSE under BOTH implemented branches if "add" means APPEND a base:
    //     • Wallace: appending ANY base raises Tm (+2 for A/T, +4 for G/C) — never
    //       decreases.
    //     • Marmur-Doty: appending A/T raises N at fixed GC, so it can change Tm either
    //       way depending on the current GC/N ratio — it is NOT monotone-down in general.
    //   The defensible, formula-guaranteed readings encoded here are:
    //     • INV — deterministic & case-insensitive: identical (up to case / ignored
    //       non-ACGT padding) input ⇒ identical Tm. (Holds for both branches.)
    //     • MON (GC contributes more than AT, APPEND) — appending a strong (G/C) base
    //       raises Tm strictly MORE than appending a weak (A/T) base:
    //       Tm(s+'G') > Tm(s+'A') and Tm(s+'C') > Tm(s+'T'). This is the honest reading
    //       of "add GC ↑ vs add AT". (Both branches; see per-branch proof below.)
    //     • MON (raise GC CONTENT at fixed length) — substituting an A/T position with a
    //       G/C base (length unchanged) STRICTLY increases Tm; the converse G/C→A/T
    //       STRICTLY decreases it. THIS is the precise, true sense of "add AT → Tm
    //       decreases": at fixed length, swapping a strong base out for a weak one (i.e.
    //       making the composition more AT-rich) lowers Tm. (Both branches.)
    //     • Wallace EXACT increments — the doc pins the constants 2 (A/T) and 4 (G/C),
    //       so in the short-oligo regime we additionally assert the exact deltas:
    //       append G/C ⇒ +4, append A/T ⇒ +2, substitute A/T→G/C ⇒ +2.
    //
    //   Per-branch proof of the comparative-append relation Tm(s+strong) > Tm(s+weak),
    //   with both extended strings kept inside the SAME branch:
    //     • Wallace (|s|+1 < 14): Δ = +4 (strong) vs +2 (weak) ⇒ strict, gap exactly 2.
    //     • Marmur-Doty (|s|+1 ≥ 14, same N = |s|+1 for both, both kept above the 0-clamp):
    //         Tm(s+strong) − Tm(s+weak) = 41·(GC+1 − GC)/N = 41/N > 0 ⇒ strict.
    //   The substitution relation at fixed length N:
    //     • Wallace: one base flips A/T(+2)→G/C(+4) ⇒ Δ = +2 (strict).
    //     • Marmur-Doty: N fixed, GC ↑ by 1 ⇒ Δ = 41/N > 0 (strict; the unclamped value
    //       only rises, so max(0,·) preserves the strict order whenever the larger one is
    //       positive — guaranteed when we build GC-rich enough inputs).
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>Counted-base branch threshold: validLength &lt; 14 ⇒ Wallace, else Marmur-Doty.</summary>
    private const int WallaceMaxLen = 14; // == ThermoConstants.WallaceMaxLength

    private static double Tm(string primer) => PrimerDesigner.CalculateMeltingTemperature(primer);

    /// <summary>A random DNA string over {A,C,G,T} of the given length (fixed-seed <see cref="Rng"/>).</summary>
    private static string RandomDnaTm(int length) => RandomDna(length);

    /// <summary>
    /// SHORT bodies (extended length stays &lt; 14 counted bases) — the Wallace regime,
    /// where exact +2 / +4 increments hold. Each is GC/AT-balanced enough that an append
    /// keeps the result inside the short-oligo branch.
    /// </summary>
    private static IEnumerable<string> ShortBodies()
    {
        yield return "A";
        yield return "ACGT";
        yield return "AATT";
        yield return "GCGC";
        yield return "ACGTACG";       // 7-mer
        yield return "ACGTACGTAC";    // 10-mer
        yield return "AAATTTGGGCC";   // 11-mer, mixed
    }

    /// <summary>
    /// LONG bodies whose appended/substituted variants stay in the Marmur-Doty regime
    /// (≥ 14 counted bases) AND keep Tm above the 0-clamp, so the strict comparative
    /// relations are theoretically guaranteed. They are GC-bearing (GC ≥ ~50%) and
    /// length ≥ 14, giving 64.9 + 41·(GC−16.4)/N well above 0.
    /// </summary>
    private static IEnumerable<string> LongBodies()
    {
        yield return "ACGTACGTACGTAC";                    // 14-mer, 50% GC
        yield return "ACGTACGTACGTACGT";                  // 16-mer
        yield return "GCGCGCGCGCGCGCGCGC";                // 18-mer, GC-rich
        yield return "ACGTACGTACGTACGTACGT";              // 20-mer, 50% GC
        yield return "GGCCGGCCGGCCGGCCGGCCGGCC";          // 24-mer, GC-rich
        yield return SanitizeGcRich(RandomDnaTm(30));     // fixed-seed random, GC-bearing
        yield return SanitizeGcRich(RandomDnaTm(50));
    }

    /// <summary>
    /// Forces a sequence to carry a healthy GC fraction (every A/T at an even index is
    /// turned into G/C) so its Marmur-Doty Tm sits safely above the 0-clamp — keeping the
    /// strict comparative relations in their guaranteed region.
    /// </summary>
    private static string SanitizeGcRich(string s)
    {
        var c = s.ToCharArray();
        for (int i = 0; i < c.Length; i += 2)
            if (c[i] is 'A' or 'T')
                c[i] = (i % 4 == 0) ? 'G' : 'C';
        return new string(c);
    }

    #region INV — same sequence → same Tm (deterministic, case-insensitive)

    [Test]
    [Description("INV: CalculateMeltingTemperature is deterministic — identical input yields identical Tm across repeated calls.")]
    public void CalculateMeltingTemperature_SameSequence_IsDeterministic()
    {
        foreach (var seq in ShortBodies().Concat(LongBodies()))
        {
            double first = Tm(seq);
            for (int i = 0; i < 5; i++)
                Tm(seq).Should().Be(first,
                    because: "Tm is a pure function of base counts and length, so repeated calls on the same input must agree exactly");
        }
    }

    [Test]
    [Description("INV: Tm is case-insensitive — the implementation uppercases input, so lower/mixed case gives the same Tm.")]
    public void CalculateMeltingTemperature_CaseFolding_PreservesTm()
    {
        foreach (var seq in ShortBodies().Concat(LongBodies()))
        {
            double upper = Tm(seq.ToUpperInvariant());

            Tm(seq.ToLowerInvariant()).Should().Be(upper,
                because: "the implementation calls ToUpperInvariant before counting bases, so case carries no information");

            // Mixed case: alternate the case of every other character.
            var mixed = new string(seq.Select((ch, i) => i % 2 == 0 ? char.ToUpperInvariant(ch) : char.ToLowerInvariant(ch)).ToArray());
            Tm(mixed).Should().Be(upper,
                because: "case folding is total, so any case pattern of the same letters yields the same counted A/T and G/C");
        }
    }

    [Test]
    [Description("INV: non-ACGT characters are ignored, so interleaving them (same A/C/G/T in order) leaves Tm unchanged.")]
    public void CalculateMeltingTemperature_NonAcgtPadding_PreservesTm()
    {
        // The spec (§3.3, §5.2) states only A/C/G/T contribute to the counted length;
        // every other character is ignored. So padding with N/-/spaces is a no-op for Tm.
        foreach (var seq in ShortBodies().Concat(LongBodies()))
        {
            double clean = Tm(seq);
            string padded = string.Concat(seq.Select(ch => ch + "N-")); // sprinkle ignored chars
            Tm(padded).Should().Be(clean,
                because: "only A/C/G/T are counted, so non-DNA characters cannot change the counted A/T, G/C, or length");
        }
    }

    #endregion

    #region MON — appending a G/C base raises Tm strictly more than appending an A/T base

    [Test]
    [Description("MON (Wallace): in the short-oligo regime, appending a G/C base adds exactly +4 °C and an A/T base exactly +2 °C, so Tm(s+strong) > Tm(s+weak) by exactly 2.")]
    public void CalculateMeltingTemperature_AppendStrongVsWeak_Wallace_ExactIncrements()
    {
        foreach (var body in ShortBodies())
        {
            // Keep the extended string strictly inside the Wallace branch (< 14 counted).
            if (body.Length + 1 >= WallaceMaxLen)
                continue;

            double baseTm = Tm(body);

            foreach (char gc in new[] { 'G', 'C' })
            {
                Tm(body + gc).Should().Be(baseTm + ThermoConstants.WallaceGcContribution,
                    because: $"the Wallace rule adds exactly +{ThermoConstants.WallaceGcContribution} °C for a G/C base appended to '{body}'");
            }

            foreach (char at in new[] { 'A', 'T' })
            {
                Tm(body + at).Should().Be(baseTm + ThermoConstants.WallaceAtContribution,
                    because: $"the Wallace rule adds exactly +{ThermoConstants.WallaceAtContribution} °C for an A/T base appended to '{body}'");
            }

            // The comparative MR: strong append strictly beats weak append (gap = 4 − 2 = 2).
            Tm(body + 'G').Should().BeGreaterThan(Tm(body + 'A'),
                because: "appending a strong (G/C) base raises Tm more than appending a weak (A/T) base — the honest reading of 'add GC ↑ vs add AT'");
            Tm(body + 'C').Should().BeGreaterThan(Tm(body + 'T'),
                because: "G/C contributes +4 vs A/T's +2 in the Wallace rule, so the strong append wins");
            (Tm(body + 'G') - Tm(body + 'A')).Should().Be(
                ThermoConstants.WallaceGcContribution - ThermoConstants.WallaceAtContribution,
                because: "the exact Wallace gap between a strong and a weak append is 4 − 2 = 2 °C");
        }
    }

    [Test]
    [Description("MON (Marmur-Doty): in the long-oligo regime, appending a G/C base raises Tm strictly more than appending an A/T base (Δ = 41/N > 0).")]
    public void CalculateMeltingTemperature_AppendStrongVsWeak_MarmurDoty_StrictOrder()
    {
        foreach (var body in LongBodies())
        {
            // Both extended strings share N = body.Length + 1 ≥ 14 (same Marmur-Doty branch).
            (body.Length + 1).Should().BeGreaterThanOrEqualTo(WallaceMaxLen,
                because: "the long-body fixture is constructed so the appended variant stays in the Marmur-Doty branch");

            double withG = Tm(body + 'G');
            double withC = Tm(body + 'C');
            double withA = Tm(body + 'A');
            double withT = Tm(body + 'T');

            // GC-rich bodies keep all four extensions above the 0-clamp, so the strict
            // unclamped order survives the Math.Max(0, ·).
            withG.Should().BeGreaterThan(withA,
                because: "at the same length N, the G append has one more G/C than the A append, adding 41/N > 0 °C");
            withC.Should().BeGreaterThan(withT,
                because: "at the same length N, the C append has one more G/C than the T append, adding 41/N > 0 °C");

            // A/T appends are interchangeable, and so are G/C appends (Tm depends on GC count, not identity).
            withG.Should().Be(withC, because: "Marmur-Doty Tm depends only on GC COUNT and N, so G and C appends are equivalent");
            withA.Should().Be(withT, because: "Marmur-Doty Tm depends only on GC COUNT and N, so A and T appends are equivalent");
        }
    }

    #endregion

    #region MON — raising GC CONTENT at fixed length raises Tm (substitution); lowering it (→AT) decreases Tm

    [Test]
    [Description("MON (Wallace): substituting an A/T position with a G/C base at fixed length raises Tm by exactly +2 °C; the reverse G/C→A/T lowers it by 2 — the precise sense of 'add AT → Tm decreases'.")]
    public void CalculateMeltingTemperature_SubstituteAtToGc_Wallace_RaisesByTwo()
    {
        // Short, in-Wallace bodies that contain at least one A/T to upgrade.
        foreach (var body in new[] { "AATT", "ACGTACG", "ACGTACGTAC", "AAATTTGGGCC" })
        {
            body.Length.Should().BeLessThan(WallaceMaxLen, because: "this sub-test asserts exact Wallace increments");

            int idx = body.IndexOfAny(new[] { 'A', 'T' });
            idx.Should().BeGreaterThanOrEqualTo(0, because: "the body must contain a weak base to upgrade");

            double baseTm = Tm(body);
            string upgraded = body[..idx] + 'G' + body[(idx + 1)..]; // A/T → G, length unchanged

            Tm(upgraded).Should().Be(baseTm + (ThermoConstants.WallaceGcContribution - ThermoConstants.WallaceAtContribution),
                because: "swapping a weak (A/T, +2) base for a strong (G/C, +4) one at fixed length raises Wallace Tm by exactly +2 °C");

            // The DUAL, encoding 'add AT → Tm decreases' literally: G/C → A/T strictly lowers Tm.
            int gcIdx = body.IndexOfAny(new[] { 'G', 'C' });
            if (gcIdx >= 0)
            {
                string atified = body[..gcIdx] + 'A' + body[(gcIdx + 1)..]; // G/C → A, length unchanged
                Tm(atified).Should().Be(baseTm - (ThermoConstants.WallaceGcContribution - ThermoConstants.WallaceAtContribution),
                    because: "making the composition more AT-rich at fixed length (G/C → A/T) lowers Wallace Tm by exactly 2 °C");
                Tm(atified).Should().BeLessThan(baseTm,
                    because: "this is the literal, formula-true reading of 'add AT → Tm decreases': swapping a strong base for a weak one at fixed length");
            }
        }
    }

    [Test]
    [Description("MON (Marmur-Doty): at fixed length, raising GC content (A/T → G/C) strictly raises Tm and lowering it (G/C → A/T) strictly lowers Tm.")]
    public void CalculateMeltingTemperature_RaiseGcContentAtFixedLength_MarmurDoty_StrictlyMonotone()
    {
        foreach (var body in LongBodies())
        {
            body.Length.Should().BeGreaterThanOrEqualTo(WallaceMaxLen, because: "fixed-length substitution must stay in the Marmur-Doty branch");

            double baseTm = Tm(body);

            // Upgrade one A/T position to G/C (length fixed) ⇒ Tm strictly UP.
            int weakIdx = body.IndexOfAny(new[] { 'A', 'T' });
            if (weakIdx >= 0)
            {
                string upgraded = body[..weakIdx] + 'G' + body[(weakIdx + 1)..];
                Tm(upgraded).Should().BeGreaterThan(baseTm,
                    because: "at fixed N, one more G/C raises GC by 1 and adds 41/N > 0 °C — higher GC content ⇒ higher Tm");
            }

            // Downgrade one G/C position to A/T (length fixed) ⇒ Tm strictly DOWN.
            int strongIdx = body.IndexOfAny(new[] { 'G', 'C' });
            if (strongIdx >= 0)
            {
                string downgraded = body[..strongIdx] + 'A' + body[(strongIdx + 1)..];
                Tm(downgraded).Should().BeLessThan(baseTm,
                    because: "at fixed N, one fewer G/C lowers GC by 1 and subtracts 41/N > 0 °C — the precise sense of 'add AT → Tm decreases'");
            }
        }
    }

    [Test]
    [Description("MON: along a chain that flips A/T positions to G/C one at a time at FIXED length, Tm is strictly increasing (monotone in GC content).")]
    public void CalculateMeltingTemperature_GcContentChain_FixedLength_StrictlyIncreasing()
    {
        // Start from an all-AT 20-mer (Marmur-Doty branch, GC = 0) and flip positions to
        // G one by one. Each flip keeps N = 20 and raises GC by 1, so Tm strictly rises.
        var chars = "ATATATATATATATATATAT".ToCharArray(); // 20 nt, all weak
        chars.Length.Should().BeGreaterThanOrEqualTo(WallaceMaxLen);

        double previous = Tm(new string(chars));
        for (int i = 0; i < chars.Length; i++)
        {
            if (chars[i] is not ('A' or 'T'))
                continue;
            chars[i] = 'G'; // upgrade weak → strong, length unchanged
            double current = Tm(new string(chars));
            current.Should().BeGreaterThan(previous,
                because: $"flipping position {i} from A/T to G raises GC content at fixed length, so Tm strictly increases (+41/N)");
            previous = current;
        }
    }

    #endregion

    #endregion

    #region PRIMER-DESIGN-001 — primer design

    // ─────────────────────────────────────────────────────────────────────────
    // Unit: PRIMER-DESIGN-001 — PCR primer design (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 22.
    //
    // API under test (PrimerDesigner.cs):
    //   GeneratePrimerCandidates(DnaSequence template, int regionStart, int regionEnd,
    //                            bool forward, PrimerParameters? parameters)
    //     → enumerates EVERY (start, len) window in [regionStart, regionEnd] with
    //       MinLength ≤ len ≤ MaxLength, evaluates each via EvaluatePrimer, and yields
    //       the resulting PrimerCandidate (valid AND invalid alike — it does NOT pre-filter).
    //   EvaluatePrimer(string seq, int position, bool isForward, PrimerParameters?)
    //     → sets IsValid = true iff EVERY screen passes; the relevant filters here are:
    //          MinLength ≤ Length ≤ MaxLength
    //          MinGcContent ≤ GcContent ≤ MaxGcContent
    //          MinTm ≤ MeltingTemperature ≤ MaxTm
    //          homopolymer ≤ MaxHomopolymer, dinuc-repeat ≤ MaxDinucleotideRepeats,
    //          !HasHairpin, (Check3PrimeStability ⇒ ΔG ≥ −9), (Avoid3PrimeGC ⇒ GC clamp).
    //
    // The MR phrase "primers" / "results" denotes the set of VALID candidates a parameter
    // set admits over a fixed region — i.e. GeneratePrimerCandidates(...).Where(IsValid).
    // A candidate's run-order-independent identity is (Position, Sequence, IsForward); its
    // metrics (GcContent, Tm, Length, IsValid) are a PURE function of that local window,
    // independent of every other candidate and of the rest of the template.
    //
    // Relation DIRECTIONS are derived from the FILTER STRUCTURE (definition), NOT from
    // observed output. Every filter is an independent conjunct, so widening exactly ONE
    // window while holding ALL other parameters fixed can only flip a candidate
    // invalid→valid (never valid→invalid):
    //
    //   • MON (wider Tm range → ≥ primers): EvaluatePrimer accepts a candidate's Tm iff
    //     MinTm ≤ Tm ≤ MaxTm. Widening the window to [MinTm', MaxTm'] ⊇ [MinTm, MaxTm]
    //     (MinTm' ≤ MinTm, MaxTm' ≥ MaxTm) keeps every previously-passing Tm passing and
    //     may admit more. With all other params fixed, the valid set is a SUPERSET and the
    //     count is non-decreasing: Valid(narrow Tm) ⊆ Valid(wide Tm).
    //
    //   • SUB (stricter GC% → ⊆ results): acceptance requires MinGc ≤ GC ≤ MaxGc.
    //     TIGHTENING the window to [MinGc', MaxGc'] ⊆ [MinGc, MaxGc] can only DROP
    //     candidates whose fixed GC% falls outside the narrower band; none are added.
    //     Along a tightening chain the valid set shrinks monotonically:
    //     Valid(strict GC) ⊆ Valid(loose GC); count non-increasing.
    //
    //   • MON (longer template → ≥ candidates): GeneratePrimerCandidates enumerates every
    //     window fitting inside [regionStart, regionEnd]. Extend the template T → T + X and
    //     the design region [0, |T|−1] → [0, |T+X|−1]. Every window that fit in the original
    //     region still fits in the extended one and—because the first |T| bases are
    //     byte-for-byte unchanged—evaluates to the IDENTICAL candidate (same Position,
    //     Sequence, IsForward, IsValid). So the original valid set is preserved EXACTLY as a
    //     subset and the count is non-decreasing: Valid(T) ⊆ Valid(T+X).
    //
    // For each relation only ONE parameter (or the region) is varied; ALL others are held
    // fixed — the essence of the metamorphic relation. Templates include hand-built
    // GC-balanced sequences that actually yield primers plus fixed-seed random ones.
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Design templates rich enough to yield several valid primers under DefaultParameters:
    /// 40–60% GC, length ≥ 60, no long homopolymer/dinucleotide runs, plus fixed-seed random.
    /// </summary>
    private static IEnumerable<string> DesignTemplates()
    {
        // ~50% GC, varied composition, 60 nt.
        yield return "ACGTGACTGACTGGATCAGTCAGTACGATCGATGCATGCATCGTAGCATGCATGCATGCA";
        // Another balanced 72-nt template.
        yield return "TGCATGCAGTCAGTACGTACGATCGATCGTAGCTAGCATGCATGCATCGATCGATCAGTCAGTACGTACGTA";
        // GC-leaning but still mostly in-range, 64 nt.
        yield return "GCGATCGATGCATCGATCGTAGCTAGCATCGATCGATGCATGCATCGATCGATGCATCGATCGT";
        // Fixed-seed random templates (relations hold for arbitrary input too).
        yield return RandomDna(80);
        yield return RandomDna(120);
        yield return RandomDna(160);
    }

    /// <summary>Run-order-independent identity of a primer candidate.</summary>
    private static (int Position, string Sequence, bool Fwd) PrimerId(PrimerCandidate c)
        => (c.Position, c.Sequence, c.IsForward);

    /// <summary>
    /// The set of VALID primer candidates a parameter set admits over [regionStart, regionEnd],
    /// keyed by identity (Position, Sequence, IsForward) so the comparison ignores run order.
    /// </summary>
    private static HashSet<(int, string, bool)> ValidPrimerIds(
        string template, int regionStart, int regionEnd, bool forward, PrimerParameters param)
        => PrimerDesigner
            .GeneratePrimerCandidates(new DnaSequence(template), regionStart, regionEnd, forward, param)
            .Where(c => c.IsValid)
            .Select(PrimerId)
            .ToHashSet();

    /// <summary>Count of VALID primer candidates over the whole template, both orientations.</summary>
    private static int ValidPrimerCount(string template, PrimerParameters param)
        => ValidPrimerIds(template, 0, template.Length, forward: true, param).Count
         + ValidPrimerIds(template, 0, template.Length, forward: false, param).Count;

    #region MON — widening the Tm range yields a superset of primers (count non-decreasing)

    [Test]
    [Description("MON: along a chain that widens [MinTm, MaxTm] (all other params fixed) the valid-primer set grows monotonically — each narrower set is a subset of the wider one.")]
    public void GeneratePrimerCandidates_WideningTmRange_YieldsSuperset_CountNonDecreasing()
    {
        // Increasingly wide Tm windows centred on the default optimum (60 °C). Each ⊇ the prior.
        (double Min, double Max)[] tmChain =
        {
            (59.5, 60.5),
            (58, 62),
            (55, 65),
            (50, 70),
            (0, 200),   // accept any Tm — the Tm filter is effectively disabled
        };

        foreach (var template in DesignTemplates())
        {
            foreach (var forward in new[] { true, false })
            {
                HashSet<(int, string, bool)>? previousSet = null;
                int previousCount = -1;

                foreach (var (min, max) in tmChain)
                {
                    var param = PrimerDesigner.DefaultParameters with { MinTm = min, MaxTm = max };
                    var currentSet = ValidPrimerIds(template, 0, template.Length, forward, param);

                    if (previousSet is not null)
                    {
                        currentSet.IsSupersetOf(previousSet).Should().BeTrue(
                            because: $"widening the Tm window to [{min},{max}] keeps every Tm that already passed " +
                                     "passing (all other filters fixed), so the valid set is a superset of the narrower one");
                        currentSet.Count.Should().BeGreaterThanOrEqualTo(previousCount,
                            because: $"a wider Tm window [{min},{max}] admits-or-keeps each candidate, never drops one — count is non-decreasing");
                    }

                    previousSet = currentSet;
                    previousCount = currentSet.Count;
                }
            }
        }
    }

    [Test]
    [Description("MON: disabling the Tm filter (accept any Tm) yields at least as many valid primers as the default narrow window, holding all else fixed.")]
    public void GeneratePrimerCandidates_DisablingTmFilter_NeverFewerThanNarrowDefault()
    {
        foreach (var template in DesignTemplates())
        {
            var narrow = PrimerDesigner.DefaultParameters; // 57–63 °C
            var wide = PrimerDesigner.DefaultParameters with { MinTm = 0, MaxTm = 200 };

            int narrowCount = ValidPrimerCount(template, narrow);
            int wideCount = ValidPrimerCount(template, wide);

            wideCount.Should().BeGreaterThanOrEqualTo(narrowCount,
                because: "the only changed parameter is the Tm window, widened to accept everything; " +
                         "every primer valid under the narrow window stays valid, so the count cannot drop");
        }
    }

    #endregion

    #region SUB — tightening the GC% bounds yields a subset of primers (count non-increasing)

    [Test]
    [Description("SUB: along a chain that TIGHTENS [MinGcContent, MaxGcContent] (all other params fixed) the valid-primer set shrinks monotonically — each stricter set is a subset of the looser one.")]
    public void GeneratePrimerCandidates_TighteningGcBounds_YieldsSubset_CountNonIncreasing()
    {
        // Increasingly STRICT GC windows, each ⊆ the previous one (chain from loose to strict).
        (double Min, double Max)[] gcChain =
        {
            (0, 100),     // accept any GC% — filter disabled
            (30, 70),
            (40, 60),     // the default band
            (45, 55),
            (48, 52),
        };

        foreach (var template in DesignTemplates())
        {
            foreach (var forward in new[] { true, false })
            {
                HashSet<(int, string, bool)>? previousSet = null;
                int previousCount = int.MaxValue;

                foreach (var (min, max) in gcChain)
                {
                    var param = PrimerDesigner.DefaultParameters with { MinGcContent = min, MaxGcContent = max };
                    var currentSet = ValidPrimerIds(template, 0, template.Length, forward, param);

                    if (previousSet is not null)
                    {
                        previousSet.IsSupersetOf(currentSet).Should().BeTrue(
                            because: $"tightening the GC window to [{min},{max}] can only DROP candidates whose fixed GC% " +
                                     "falls outside the narrower band, so the stricter set is a subset of the looser one");
                        currentSet.Count.Should().BeLessThanOrEqualTo(previousCount,
                            because: $"a stricter GC window [{min},{max}] removes-or-keeps each candidate, never adds one — count is non-increasing");
                    }

                    previousSet = currentSet;
                    previousCount = currentSet.Count;
                }
            }
        }
    }

    [Test]
    [Description("SUB: every primer valid under a stricter GC band is also valid under a looser band that contains it (membership is preserved when loosening).")]
    public void GeneratePrimerCandidates_StrictGcBand_SubsetOfLooseBand()
    {
        foreach (var template in DesignTemplates())
        {
            foreach (var forward in new[] { true, false })
            {
                var strict = ValidPrimerIds(template, 0, template.Length, forward,
                    PrimerDesigner.DefaultParameters with { MinGcContent = 45, MaxGcContent = 55 });
                var loose = ValidPrimerIds(template, 0, template.Length, forward,
                    PrimerDesigner.DefaultParameters with { MinGcContent = 40, MaxGcContent = 60 });

                loose.IsSupersetOf(strict).Should().BeTrue(
                    because: "[45,55] ⊆ [40,60]; a GC% inside the tighter band is inside the wider band too, " +
                             "and no other filter changed, so the strict valid set is a subset of the loose one");
            }
        }
    }

    #endregion

    #region MON — extending the template (and design region) never removes candidates anchored in the original region

    [Test]
    [Description("MON: extending template T → T+X and the design region to cover it preserves every valid primer anchored in the original region EXACTLY (subset), so the count is non-decreasing.")]
    public void GeneratePrimerCandidates_LongerTemplate_OriginalCandidatesPreserved_CountNonDecreasing()
    {
        var param = PrimerDesigner.DefaultParameters;

        foreach (var body in DesignTemplates())
        {
            foreach (var forward in new[] { true, false })
            {
                var baseIds = ValidPrimerIds(body, 0, body.Length, forward, param);

                // Append several distinct extensions, including GC-rich and fixed-seed random ones.
                foreach (var ext in new[] { "ACGTACGTACGT", "GGGGCCCCAAAA", RandomDna(30), RandomDna(60) })
                {
                    string extended = body + ext;
                    var extendedIds = ValidPrimerIds(extended, 0, extended.Length, forward, param);

                    extendedIds.IsSupersetOf(baseIds).Should().BeTrue(
                        because: $"the first {body.Length} bases are unchanged, so every window anchored in the original " +
                                 "region yields the IDENTICAL candidate; the extended region only ADDS windows, never removes one");
                    extendedIds.Count.Should().BeGreaterThanOrEqualTo(baseIds.Count,
                        because: "extending the template/region can only add candidate windows downstream — the original valid set survives, so the count is non-decreasing");
                }
            }
        }
    }

    [Test]
    [Description("MON: appending bases after the design region's end leaves the in-region valid set IDENTICAL — the appended tail lies outside [regionStart, regionEnd].")]
    public void GeneratePrimerCandidates_AppendOutsideFixedRegion_InRegionSetUnchanged()
    {
        var param = PrimerDesigner.DefaultParameters;

        foreach (var body in DesignTemplates())
        {
            foreach (var forward in new[] { true, false })
            {
                int regionEnd = body.Length; // pin the design region to the original body
                var baseIds = ValidPrimerIds(body, 0, regionEnd, forward, param);

                foreach (var ext in new[] { "ACGTACGTACGT", RandomDna(40) })
                {
                    // Hold the region FIXED at [0, regionEnd]; the appended tail is outside it.
                    var sameRegionIds = ValidPrimerIds(body + ext, 0, regionEnd, forward, param);

                    sameRegionIds.SetEquals(baseIds).Should().BeTrue(
                        because: "GeneratePrimerCandidates only enumerates windows inside [regionStart, regionEnd]; " +
                                 "bases appended beyond regionEnd are never read, so the in-region valid set is exactly preserved");
                }
            }
        }
    }

    #endregion

    #endregion

    #region PRIMER-STRUCT-001 — primer secondary structure (self-dimer / hairpin)

    // ─────────────────────────────────────────────────────────────────────────
    // Unit: PRIMER-STRUCT-001 — primer secondary structure (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 23.
    //
    // API under test (PrimerDesigner.cs):
    //   HasPrimerDimer(primer1, primer2, minComplementarity = 4) → bool
    //     1. seq1 = primer1.upper();  seq2 = revComp(primer2.upper()).
    //     2. checkLength = min(8, len1, len2).
    //     3. end1 = last `checkLength` bases of seq1; end2 = first `checkLength` bases of seq2.
    //     4. complementary = #{ i : IsComplementary(end1[i], end2[i]) }.
    //     5. return complementary >= minComplementarity.
    //   HasHairpinPotential(seq, minStemLength = 4, minLoopLength = 3) → bool
    //     true ⇔ ∃ a length-minStemLength window at i complementary (reversed) to a
    //     window at j with j ≥ i + minStemLength + minLoopLength (i.e. a stem-loop
    //     with stem ≥ minStem and loop ≥ minLoop). false below 2·minStem + minLoop.
    //
    // Source (semantics pinned from spec, NOT from observed output):
    //   docs/algorithms/Molecular_Tools/Primer_Structure_Analysis.md §2.1–§2.2, §4.2;
    //   docs/algorithms/MolTools/Primer_Design.md.
    //   "Primer-dimers arise when primers have complementary 3' ends"; the dimer signal
    //   is the COUNT of complementary base pairs in the compared 3' window (more
    //   self-complementarity ⇒ stronger dimer). "Hairpins form when a primer contains
    //   self-complementary regions separated by a loop" — the boolean depends ONLY on
    //   the existence of such a stem-loop.
    //
    // ── Self-dimer score (the MON subject) ──
    //   For a SELF-dimer we evaluate HasPrimerDimer(primer, primer, ·). The hidden
    //   score is the complementary-pair count of step 4. Because HasPrimerDimer exposes
    //   only the boolean `count >= minComplementarity`, we RECOVER the exact score by a
    //   threshold sweep: SelfDimerScore = the largest threshold t in [0, checkLen] for
    //   which HasPrimerDimer(primer, primer, t) is still true; that t equals `count`.
    //   This is the genuine score the implementation computes, read through its only
    //   public surface — not an output-fitted proxy.
    //
    //   Geometry of the SELF window (derived, not observed): for self-dimer
    //   seq2 = revComp(primer), so end2[i] = complement(w[checkLen-1-i]) where w is the
    //   terminal checkLen-mer and end1[i] = w[i]. IsComplementary(w[i], complement(x))
    //   holds iff w[i] == x, so a comparison position i pairs iff w[i] == w[checkLen-1-i].
    //   So the self-dimer score COUNTS the positions where the terminal checkLen-mer
    //   equals its own REVERSE (a reverse-EQUAL palindrome — A·T/C·G pairing across the
    //   antiparallel self-fold maps "complementary" back to "equal"):
    //     • a terminal window that is a perfect reverse-palindrome (w == reverse(w))
    //       scores checkLen (the MAXIMUM, fully self-complementary 3' end);
    //     • a window whose every mirror pair (i, checkLen-1-i) differs scores 0.
    //   (The matched count is always even: a matched mirror pair contributes 2 positions.)
    //   MON: a primer with MORE self-complementary terminal pairs has a score ≥ one with
    //   fewer; strictly higher when the added pairing is real.
    //
    // ── Hairpin INV/MON dependency (the exact thing tested) ──
    //   The hairpin boolean depends ONLY on whether some self-complementary stem pair
    //   (stem ≥ minStem, loop ≥ minLoop) exists inside the sequence. Appending bases at
    //   the 3' end preserves every original substring/position, so it can NEVER remove a
    //   stem ⇒ detection is MONOTONE: once true, any 3' append keeps it true (false→true
    //   flips are real — an append can fabricate a new stem, including a JUNCTION window
    //   that mixes the body's tail bases with appended bases). We therefore assert exact
    //   "same result" only where it is provably sound: appending a HOMOPOLYMER X^k whose
    //   complement is ABSENT from the body. Then no stem arm overlapping the appended X's —
    //   all-X, all-body, or a mixed junction window — can find a reverse-complementary
    //   partner (a partner needs complement(X) at the X positions, and complement(X) occurs
    //   nowhere in body+append because the append is pure X and X ≠ complement(X)). With no
    //   stem added and none removable, HasHairpinPotential is preserved exactly, both ways.
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>The terminal comparison window cap used by HasPrimerDimer (min(8, len1, len2)).</summary>
    private const int DimerCheckCap = 8;

    private const int HairpinMinStem = 4;
    private const int HairpinMinLoop = 3;

    /// <summary>
    /// Recovers the implementation's hidden self-dimer score — the complementary-pair
    /// count in the 3' window — via the only public surface (HasPrimerDimer's boolean
    /// `count &gt;= threshold`): the largest threshold still returning true equals the count.
    /// </summary>
    private static int SelfDimerScore(string primer)
    {
        int checkLen = Math.Min(DimerCheckCap, primer.Length);
        int score = 0;
        for (int t = 1; t <= checkLen; t++)
            if (PrimerDesigner.HasPrimerDimer(primer, primer, minComplementarity: t))
                score = t;
        return score;
    }

    /// <summary>DNA complement of a single base (uppercase ACGT).</summary>
    private static char Complement(char b) => b switch
    {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => throw new ArgumentException($"Unexpected base '{b}'.", nameof(b)),
    };

    /// <summary>
    /// Picks an append base X whose homopolymer run is guaranteed to add NO hairpin stem
    /// to <paramref name="body"/>: the body must contain NO base equal to complement(X) at
    /// all. Then no stem arm — whether all-X, all-body, or a junction window mixing body
    /// tail bases with appended X's — can find a reverse-complementary partner that pairs
    /// with the X positions (a partner needs complement(X) where X sits, and complement(X)
    /// occurs nowhere in body+append, since the append is pure X and X ≠ complement(X)).
    /// Returns null if every base's complement appears in the body (then the body is not
    /// usable for the exact-invariance relation and is skipped).
    /// </summary>
    private static char? SafeHairpinAppendBase(string body)
    {
        string upper = body.ToUpperInvariant();
        foreach (char x in "ACGT")
            if (!upper.Contains(Complement(x)))
                return x;
        return null;
    }

    #region MON — more self-complementary 3' end → higher self-dimer score

    [Test]
    [Description("MON: making one more terminal mirror-pair self-complementary STRICTLY raises the self-dimer score; the chain is monotonically non-decreasing up to the window cap.")]
    public void HasPrimerDimer_ExtendSelfComplementaryCore_ScoreStrictlyIncreasing()
    {
        // A fixed 5' leader keeps the discriminating signal in the 3' terminal 8-base
        // window. Each tail below is the previous one with one MORE mirror pair
        // (i, 7-i) made equal — turning two more comparison positions into self-paired
        // ones. So the recovered self-dimer score (reverse-palindrome position count of
        // the last 8 bases) rises by exactly 2 at each step: 0, 2, 4, 6, 8.
        const string leader = "AAAA"; // primer = leader + 8-base tail ⇒ window == tail
        string[] tailsByRisingSelfComp =
        {
            "ACGTACGT", // 0 matched mirror pairs: (A,T)(C,G)(G,C)(T,A) all differ
            "ACGTTCGT", // + pair (3,4)=T,T            ⇒ 2
            "ACGTTGGT", // + pair (2,5)=G,G            ⇒ 4
            "ACGTTGCT", // + pair (1,6)=C,C            ⇒ 6
            "ACGTTGCA", // + pair (0,7)=A,A (reverse-palindrome) ⇒ 8 (window cap)
        };

        int previous = -1;
        foreach (var tail in tailsByRisingSelfComp)
        {
            int score = SelfDimerScore(leader + tail);

            score.Should().BeGreaterThan(previous,
                because: $"making one more terminal mirror pair self-complementary (tail '{tail}') turns two more " +
                         "comparison positions into pairing ones, so the self-dimer score must strictly rise");
            previous = score;
        }

        // The fully self-complementary (reverse-palindromic) terminal window saturates the cap.
        SelfDimerScore(leader + "ACGTTGCA").Should().Be(DimerCheckCap,
            because: "the 8-base 3' window 'ACGTTGCA' equals its own reverse, so every comparison position pairs — the maximum");
    }

    [Test]
    [Description("MON: a primer with strictly more self-complementary terminal pairs never scores below one with fewer — checked across several primers incl. fixed-seed random.")]
    public void HasPrimerDimer_MoreSelfComplementaryTerminus_ScoresAtLeastAsHigh()
    {
        // For each body we compare its self-dimer score to the same body with one MORE
        // mirror pair made self-complementary at the 3' end. We make the terminal window
        // self-palindromic by overwriting its outermost mismatching mirror pair, which
        // can only turn a non-paired terminal position into a paired one.
        foreach (var body in StructureSamples())
        {
            int baseScore = SelfDimerScore(body);

            // Force a perfectly self-complementary 8-base 3' end: append a window that
            // equals its own reverse. This is the self-complementarity MAXIMUM.
            string maximallySelfComp = body + "ACGTTGCA"; // 8-base reverse-palindrome tail
            int maxScore = SelfDimerScore(maximallySelfComp);

            maxScore.Should().BeGreaterThanOrEqualTo(baseScore,
                because: "replacing the 3' terminus with a fully self-complementary window can only add complementary " +
                         "pairs to the compared window, so its self-dimer score is ≥ the original");
            maxScore.Should().Be(DimerCheckCap,
                because: "a reverse-palindromic 3' window (window == reverse(window)) pairs at every comparison position — the documented maximum");
        }
    }

    [Test]
    [Description("MON anchors: a fully self-complementary (palindromic) 3' end scores the maximum; a 3' end with no self-complementary mirror pair scores the minimum (0).")]
    public void HasPrimerDimer_KnownExtremes_MaxForPalindromeZeroForNonComplementary()
    {
        // Maximum: the 8-base 3' window equals its own reverse (every mirror pair matches).
        SelfDimerScore("AAAAACGTTGCA").Should().Be(DimerCheckCap,
            because: "the 8-base 3' window 'ACGTTGCA' equals reverse('ACGTTGCA'), so all 8 comparison positions pair — the maximum");
        // An all-same 3' end is the trivial reverse-palindrome and also hits the maximum.
        SelfDimerScore("GCGCAAAAAAAA").Should().Be(DimerCheckCap,
            because: "an all-A 8-base 3' window is its own reverse, so every comparison position pairs — the maximum");

        // Minimum: build an 8-base 3' window whose every mirror pair (i, 7-i) DIFFERS,
        // so no comparison position is self-paired → score 0.
        //   window 'ACGTACGT': pairs (A,T)(C,G)(G,C)(T,A) — every mirror pair differs.
        SelfDimerScore("AAAAACGTACGT").Should().Be(0,
            because: "in the 3' window 'ACGTACGT' every mirror pair (i, 7-i) differs, so none is self-complementary — the minimum 0");
    }

    #endregion

    #region INV — a 3' append preserves a detected hairpin (monotone); a complement-absent homopolymer preserves the result exactly

    // NOTE on theory: an exact "append → same hairpin" relation is FALSE for an arbitrary
    // append, because a 3' extension can fabricate a NEW stem — including a junction window
    // that mixes the body's tail bases with appended bases and happens to be reverse-
    // complementary to an existing body window. So we assert only what the algorithm
    // actually guarantees:
    //   (1) MONOTONE preservation — a 3' append NEVER destroys a detected hairpin
    //       (HasHairpinPotentialSimple scans all i<j windows; extending the string only
    //       ADDS (i,j) pairs and never shortens an existing complementary pair). Rigorous
    //       for ANY appended region.
    //   (2) EXACT invariance under a complement-absent homopolymer — if the body contains
    //       no base equal to complement(X), then appending X^k adds no stem at all (no arm
    //       overlapping the X's can find a partner), and removes none, so the boolean is
    //       byte-for-byte preserved in BOTH directions.

    [Test]
    [Description("MONOTONE: a 3' append NEVER destroys a detected hairpin — once HasHairpinPotential is true for the body it stays true after ANY appended region.")]
    public void HasHairpinPotential_ThreePrimeAppend_NeverDestroysDetectedHairpin()
    {
        bool sawHairpin = false;

        // Includes "GCGCTTTTTGCGC" (a guaranteed hairpin) so the relation is exercised,
        // not vacuous, regardless of the random fixtures.
        foreach (var body in StructureSamples())
        {
            bool baseResult = PrimerDesigner.HasHairpinPotential(body, HairpinMinStem, HairpinMinLoop);
            if (!baseResult)
                continue;
            sawHairpin = true;

            foreach (var ext in new[] { "A", "C", "G", "T", "AAAA", "GCGC", NonPamRegion(6), RandomDna(8) })
            {
                PrimerDesigner.HasHairpinPotential(body + ext, HairpinMinStem, HairpinMinLoop).Should().BeTrue(
                    because: $"extending the 3' end with '{ext[..Math.Min(4, ext.Length)]}…' only ADDS candidate stem-pair windows " +
                             "and never shortens an existing one, so a hairpin detected in the body is still detected");
            }
        }

        sawHairpin.Should().BeTrue(because: "at least one fixture (the GCGC…GCGC hairpin) forms a hairpin, so the monotone relation is actually tested");
    }

    [Test]
    [Description("EXACT INV: appending a homopolymer X whose complement is ABSENT from the body leaves HasHairpinPotential byte-for-byte unchanged — no arm overlapping the X's can ever pair, so no stem is added or removed.")]
    public void HasHairpinPotential_AppendComplementAbsentHomopolymer_ResultUnchanged()
    {
        // Deterministic bodies over three bases (so one base's complement is absent),
        // covering BOTH a hairpin-positive and a hairpin-negative case:
        //   "GGCCAAAGGCC" — stem GGCC … (AAA loop) … GGCC, no T  ⇒ hairpin = true.
        //   "ACGCAGCAGCAG" — no inverted-repeat stem, no T       ⇒ hairpin = false.
        var bodies = new[] { "GGCCAAAGGCC", "ACGCAGCAGCAG" };
        bool sawTrue = false, sawFalse = false;

        foreach (var body in bodies)
        {
            char? x = SafeHairpinAppendBase(body);
            x.Should().NotBeNull(because: $"'{body}' omits a base, so the omitted base is a complement-absent append target");

            bool baseResult = PrimerDesigner.HasHairpinPotential(body, HairpinMinStem, HairpinMinLoop);
            sawTrue |= baseResult;
            sawFalse |= !baseResult;

            foreach (int runLen in new[] { 1, 3, 5, 8 })
            {
                string ext = new string(x!.Value, runLen);
                PrimerDesigner.HasHairpinPotential(body + ext, HairpinMinStem, HairpinMinLoop).Should().Be(baseResult,
                    because: $"complement('{x}') occurs nowhere in '{body}', so a '{x}'×{runLen} append (pure {x}) gives every " +
                             "X-overlapping window no possible reverse-complementary partner — no stem is added, none removed, the result is preserved exactly");
            }
        }

        sawTrue.Should().BeTrue(because: "the GGCC…GGCC body forms a hairpin, exercising exact true→true invariance");
        sawFalse.Should().BeTrue(because: "the non-stem body forms none, exercising exact false→false invariance");
    }

    #endregion

    /// <summary>
    /// Primers for the secondary-structure relations: deliberate palindromic / non-
    /// palindromic 3' ends, structured bodies, and fixed-seed random primers (relations
    /// must hold for arbitrary input too).
    /// </summary>
    private static IEnumerable<string> StructureSamples()
    {
        yield return "ACGTACGTACGTACGTACGT";      // 20-mer, mixed
        yield return "GCGCTTTTTGCGC";             // a hairpin-forming primer
        yield return "AAAAAAAAAAAA";              // flat, no structure
        yield return "ATATATATATAT";              // weak alternating
        yield return "GACGTCAAAAAA";              // self-palindromic head, flat tail
        yield return "TTTTTTTTGACATGTC";          // palindromic 3' end
        yield return RandomDna(20);               // fixed-seed random
        yield return RandomDna(24);
        yield return RandomDna(30);
    }

    #endregion

    #region PROBE-DESIGN-001 — hybridization probe design

    // ─────────────────────────────────────────────────────────────────────────
    // Unit: PROBE-DESIGN-001 — hybridization probe design (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 24.
    //
    // API under test (ProbeDesigner.cs):
    //   DesignProbes(target, ProbeParameters?, maxProbes)
    //     1. reject if target shorter than MinLength; uppercase target.
    //     2. enumerate EVERY (start, len) window with MinLength ≤ len ≤ MaxLength.
    //     3. EvaluateProbeWithGc: raw score starts at 1.0 and SOFT penalties subtract
    //        (doc §2.2): GC out of [MinGc,MaxGc] −0.3, Tm out of [MinTm,MaxTm] −0.3,
    //        homopolymer>Max −0.2, selfComp>Max −0.2, structure −0.15, repeats −0.1,
    //        5'/3' G·C −0.02 each. A window is a CANDIDATE iff raw score > 0
    //        (EvaluateProbeWithGc returns null at score ≤ 0) — except GC beyond the
    //        early-reject band [MinGc−0.1, MaxGc+0.1] is dropped BEFORE scoring.
    //     4. return the top-`maxProbes` candidates by descending score.
    //   DesignProbes(target, ISuffixTree genomeIndex, params, maxProbes, requireUnique)
    //     forms the same raw-score shortlist, then (requireUnique) SKIPS any probe with
    //     CheckSpecificity < 1.0, i.e. any probe occurring more than once in the genome
    //     index (INV-03: spec = 0 / 1 / 1·hits⁻¹).
    //
    // A probe's run-order-independent identity is (Start, Sequence); its full record
    // (Tm, GcContent, Score, Warnings) is a PURE function of that one window's substring
    // — independent of every OTHER candidate and of the rest of the target.
    //
    // Relation DIRECTIONS are derived from this filter/scoring structure (definition),
    // NOT from observed output. To isolate each metamorphic relation we lift the
    // top-K cap by designing with a very large `maxProbes` (AllProbes) so the returned
    // set is exactly {candidates : raw score > 0} (no truncation churn), then vary ONE
    // thing and hold the rest fixed.
    //
    //   • MON (wider Tm → ≥ probes): Tm is one independent SOFT penalty. Widening
    //     [MinTm,MaxTm] removes the −0.3 Tm penalty from every window whose Tm enters
    //     the wider band and changes NOTHING else, so each window's raw score is
    //     NON-DECREASING ⇒ {raw score > 0} grows ⇒ probe count is non-decreasing and the
    //     valid set is a SUPERSET. The exact, non-vacuous mechanism is asserted directly:
    //     for a probe shared by the narrow- and wide-Tm designs, every field but Score is
    //     identical and Score rises by EXACTLY the documented 0.3 iff the probe's Tm lay
    //     outside the narrow window (the salt-adjusted formula puts a balanced 50–60mer
    //     near 69 °C, below the 75–85 °C Microarray window, so the toggle is real).
    //
    //   • SUB (stricter uniqueness → ⊆ results): `requireUnique=true` applies an EXTRA
    //     conjunct (skip specificity < 1) over the SAME ordered candidate stream that
    //     `requireUnique=false` yields. With the cap lifted it therefore yields a SUBSET:
    //     ids(requireUnique) ⊆ ids(¬requireUnique), count non-increasing. A genome that
    //     fully duplicates the target makes EVERY probe non-unique, so the strict filter
    //     empties the result (the SUB endpoint), while the lenient design stays non-empty.
    //
    //   • INV (unrelated region append → same probes): a probe's whole record depends
    //     only on its own window. Appending a region downstream leaves every window that
    //     lies entirely in the original target BYTE-IDENTICAL, so each original probe is
    //     preserved EXACTLY (same Start/Sequence/Tm/GcContent/Score/Warnings) and the
    //     extended valid set is a SUPERSET. For EXACT set equality we append after a
    //     poly-G tail of length ≥ MaxLength using a poly-G extension: then every window
    //     that reaches the appended bases starts inside the tail (|tail| ≥ MaxLength), so
    //     it is all-G ⇒ GC = 1.0 > MaxGc+0.1 ⇒ early-rejected ⇒ the append creates NO new
    //     candidate and design(T) == design(T+X) exactly.
    //
    // Source (semantics pinned from spec, NOT from observed output):
    //   docs/algorithms/MolTools/Hybridization_Probe_Design.md §2.2 (penalty table),
    //   §2.4 (INV-01 raw-score-positive shortlist, INV-03 specificity 0/1/1·hits⁻¹),
    //   §3.1/§3.3 (requireUnique filters specificity < 1.0; top-maxProbes after ranking),
    //   §4.2 (Microarray default: length 50–60, Tm 75–85 °C, GC 0.40–0.60).
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>Microarray default probe parameters (length 50–60, Tm 75–85, GC 0.40–0.60).</summary>
    private static ProbeDesigner.ProbeParameters Microarray => ProbeDesigner.Defaults.Microarray;

    /// <summary>Documented Tm-window penalty (doc §2.2): a probe whose Tm is out of [MinTm,MaxTm] loses 0.3.</summary>
    private const double ProbeTmPenalty = 0.3;

    /// <summary>
    /// Large probe cap that lifts the top-K truncation for these fixtures, so a design
    /// returns exactly {candidate windows : raw score &gt; 0} — making the set relations
    /// (superset / subset / equality) reflect the SCORING/FILTER definition, not the cut.
    /// </summary>
    private const int AllProbes = 100_000;

    /// <summary>Run-order-independent identity of a designed probe.</summary>
    private static (int Start, string Sequence) ProbeId(ProbeDesigner.Probe p) => (p.Start, p.Sequence);

    /// <summary>All probes a parameter set admits over the target (cap lifted ⇒ the full score-positive set).</summary>
    private static List<ProbeDesigner.Probe> DesignAll(string target, ProbeDesigner.ProbeParameters param)
        => ProbeDesigner.DesignProbes(target, param, AllProbes).ToList();

    private static Dictionary<(int, string), ProbeDesigner.Probe> DesignById(string target, ProbeDesigner.ProbeParameters param)
        => DesignAll(target, param).ToDictionary(ProbeId);

    /// <summary>
    /// Targets ≥ MinLength with ~50% GC and varied composition, so several Microarray
    /// probes (50–60 nt) fit, plus fixed-seed random targets — relations must hold for
    /// arbitrary input too.
    /// </summary>
    private static IEnumerable<string> ProbeTargets()
    {
        yield return "ACGTGACTGACTGGATCAGTCAGTACGATCGATGCATGCATCGTAGCATGCATGCATGCAACGTGACTGACTGGATCAGT";
        yield return "TGCATGCAGTCAGTACGTACGATCGATCGTAGCTAGCATGCATGCATCGATCGATCAGTCAGTACGTACGTAGCATCGAT";
        yield return "GCGATCGATGCATCGATCGTAGCTAGCATCGATCGATGCATGCATCGATCGATGCATCGATCGTACGATCGTAGCTAGCA";
        yield return RandomDna(90);
        yield return RandomDna(120);
    }

    #region MON — widening the Tm window never removes probes; the Tm penalty toggles by exactly 0.3

    [Test]
    [Description("MON: along a chain that widens [MinTm,MaxTm] (all else fixed) the valid-probe set grows monotonically — each narrower set is a subset of the wider one, count non-decreasing.")]
    public void DesignProbes_WideningTmWindow_YieldsSuperset_CountNonDecreasing()
    {
        // Increasingly wide Tm windows, each ⊇ the prior (Microarray default 75–85 outward).
        (double Min, double Max)[] tmChain =
        {
            (75, 85),   // Microarray default
            (70, 90),
            (50, 100),
            (0, 200),   // Tm filter effectively disabled
        };

        foreach (var target in ProbeTargets())
        {
            HashSet<(int, string)>? previousSet = null;
            int previousCount = -1;

            foreach (var (min, max) in tmChain)
            {
                var ids = DesignAll(target, Microarray with { MinTm = min, MaxTm = max })
                    .Select(ProbeId).ToHashSet();

                if (previousSet is not null)
                {
                    ids.IsSupersetOf(previousSet).Should().BeTrue(
                        because: $"widening the Tm window to [{min},{max}] only removes the −{ProbeTmPenalty} Tm penalty " +
                                 "(all other penalties fixed), so every window that already scored > 0 still does — the valid set is a superset");
                    ids.Count.Should().BeGreaterThanOrEqualTo(previousCount,
                        because: $"a wider Tm window [{min},{max}] can only raise raw scores, never lower one past 0 — probe count is non-decreasing");
                }

                previousSet = ids;
                previousCount = ids.Count;
            }
        }
    }

    [Test]
    [Description("MON mechanism: a probe shared by the narrow- and wide-Tm designs keeps every field but Score; its Score rises by EXACTLY the documented 0.3 iff its Tm was outside the narrow window.")]
    public void DesignProbes_WideningTmWindow_RaisesSharedProbeScoreByExactlyTheTmPenalty()
    {
        bool sawTmToggle = false;

        foreach (var target in ProbeTargets())
        {
            // narrow = Microarray default (75–85 °C); wide = Tm filter disabled.
            var narrow = DesignById(target, Microarray with { MinTm = 75, MaxTm = 85 });
            var wide = DesignById(target, Microarray with { MinTm = 0, MaxTm = 200 });

            foreach (var (id, probe) in narrow)
            {
                // A probe valid under the narrow Tm window stays valid when the window widens
                // (its raw score only rises), so the wide design must contain its twin.
                wide.Should().ContainKey(id,
                    because: "removing the Tm penalty cannot drop a probe that already scored > 0 under the narrow window");
                var twin = wide[id];

                twin.Sequence.Should().Be(probe.Sequence);
                twin.Tm.Should().Be(probe.Tm, because: "Tm is a pure function of the probe window, independent of the Tm parameter range");
                twin.GcContent.Should().Be(probe.GcContent, because: "GC content depends only on the window, not on the Tm range");

                // The ONLY score difference is the Tm penalty: present in narrow iff the probe's
                // Tm is outside [75,85]; never present in the wide (0–200) window.
                bool tmOutsideNarrow = probe.Tm < 75 || probe.Tm > 85;
                double expectedDelta = tmOutsideNarrow ? ProbeTmPenalty : 0.0;
                if (tmOutsideNarrow) sawTmToggle = true;

                twin.Score.Should().BeApproximately(probe.Score + expectedDelta, 1e-9,
                    because: tmOutsideNarrow
                        ? "the probe's Tm is outside the narrow window, so widening it removes exactly the documented 0.3 Tm penalty"
                        : "the probe's Tm is already inside the narrow window, so widening the window leaves its score unchanged");
            }
        }

        sawTmToggle.Should().BeTrue(
            because: "the salt-adjusted formula puts balanced 50–60mers near 69 °C (below 75–85 °C), so at least one probe's Tm toggles — the relation is exercised, not vacuous");
    }

    #endregion

    #region SUB — requiring genome uniqueness yields a subset of the lenient design

    [Test]
    [Description("SUB: requireUnique=true applies an extra 'specificity = 1' filter over the same candidate stream, so its probe set is a subset of the requireUnique=false design (count non-increasing).")]
    public void DesignProbes_RequireUnique_YieldsSubsetOfLenientDesign()
    {
        foreach (var target in ProbeTargets())
        {
            string upper = target.ToUpperInvariant();
            // Genome where the FIRST part of the target is duplicated: probes drawn from that
            // region occur twice (non-unique), the rest occur once (unique). 'N' separators
            // carry no probe base, so they break any cross-junction match.
            int dupLen = Math.Min(upper.Length, Microarray.MinLength + 15);
            var genome = global::SuffixTree.SuffixTree.Build(upper + "NNNNNNNN" + upper[..dupLen]);

            var unique = ProbeDesigner.DesignProbes(target, genome, Microarray, AllProbes, requireUnique: true)
                .Select(ProbeId).ToHashSet();
            var lenient = ProbeDesigner.DesignProbes(target, genome, Microarray, AllProbes, requireUnique: false)
                .Select(ProbeId).ToHashSet();

            lenient.IsSupersetOf(unique).Should().BeTrue(
                because: "requiring genome uniqueness only SKIPS candidates with specificity < 1 from the same ordered stream, " +
                         "so the unique-only probe set is a subset of the lenient one");
            unique.Count.Should().BeLessThanOrEqualTo(lenient.Count,
                because: "a stricter uniqueness requirement removes-or-keeps each candidate, never adds one — count is non-increasing");
        }
    }

    [Test]
    [Description("SUB endpoint: against a genome that fully duplicates the target every probe is non-unique, so requireUnique empties the result while the lenient design stays non-empty.")]
    public void DesignProbes_FullyDuplicatedGenome_RequireUnique_EmptiesWhileLenientKeepsProbes()
    {
        foreach (var target in ProbeTargets())
        {
            string upper = target.ToUpperInvariant();
            // Two full copies of the target ⇒ every probe (a substring of the target) occurs ≥ 2×.
            var genome = global::SuffixTree.SuffixTree.Build(upper + "NNNNNNNN" + upper);

            int uniqueCount = ProbeDesigner.DesignProbes(target, genome, Microarray, AllProbes, requireUnique: true).Count();
            int lenientCount = ProbeDesigner.DesignProbes(target, genome, Microarray, AllProbes, requireUnique: false).Count();

            lenientCount.Should().BeGreaterThan(0,
                because: "the target yields probes, and the lenient design keeps non-unique candidates");
            uniqueCount.Should().Be(0,
                because: "every probe occurs in both target copies (specificity = 0.5 < 1), so the strict uniqueness filter removes them all");
        }
    }

    #endregion

    #region INV — appending an unrelated region preserves the originally designed probes

    [Test]
    [Description("INV: appending any unrelated downstream region preserves every original probe EXACTLY (same Start/Sequence/Tm/GC/Score/Warnings); the extended valid set is a superset.")]
    public void DesignProbes_AppendUnrelatedRegion_PreservesOriginalProbesExactly()
    {
        foreach (var target in ProbeTargets())
        {
            var baseProbes = DesignAll(target, Microarray);
            baseProbes.Should().NotBeEmpty(because: "each target is constructed to yield Microarray probes");

            foreach (var ext in new[] { RandomDna(30), RandomDna(60), "GGGGCCCCAAAATTTT", NonPamRegion(40) })
            {
                var extById = DesignById(target + ext, Microarray);

                foreach (var bp in baseProbes)
                {
                    extById.Should().ContainKey(ProbeId(bp),
                        because: "windows lying entirely in the original target are byte-identical after a downstream append, " +
                                 "so every original probe still appears");
                    var twin = extById[ProbeId(bp)];

                    twin.Sequence.Should().Be(bp.Sequence);
                    twin.Tm.Should().Be(bp.Tm, because: "Tm depends only on the probe window, which the downstream append does not touch");
                    twin.GcContent.Should().Be(bp.GcContent);
                    twin.Score.Should().Be(bp.Score, because: "the heuristic score reads only the probe's own window, so a distant append cannot move it");
                    twin.Warnings.Should().Equal(bp.Warnings, because: "warnings are recorded from the probe window alone");
                }
            }
        }
    }

    [Test]
    [Description("INV exact: with a poly-G tail ≥ MaxLength, a poly-G append creates no new candidate (every append-touching window is all-G ⇒ GC-rejected), so the designed probe set is exactly unchanged.")]
    public void DesignProbes_PolyGTailThenPolyGAppend_DesignSetExactlyUnchanged()
    {
        // good region (~50% GC, yields probes) + poly-G tail of length ≥ MaxLength.
        const string goodRegion = "ACGTGACTGACTGGATCAGTCAGTACGATCGATGCATGCATCGTAGCATGCATGCATGCAACGTGACTGACTGGATCAGT";
        string tail = new string('G', Microarray.MaxLength + 4); // ≥ MaxLength ⇒ append windows start inside the tail
        string baseSeq = goodRegion + tail;

        var baseById = DesignById(baseSeq, Microarray);
        baseById.Should().NotBeEmpty(because: "the balanced good region yields Microarray probes; the poly-G tail itself is GC-rejected");

        foreach (int appendLen in new[] { 1, 10, 30 })
        {
            string extended = baseSeq + new string('G', appendLen);
            var extById = DesignById(extended, Microarray);

            extById.Keys.Should().BeEquivalentTo(baseById.Keys,
                because: $"a poly-G×{appendLen} append after a poly-G tail of length ≥ MaxLength adds only all-G windows " +
                         "(GC = 100% > MaxGc+0.1), which are early-rejected — so no probe is added or removed");

            foreach (var (id, bp) in baseById)
            {
                var twin = extById[id];
                twin.Sequence.Should().Be(bp.Sequence);
                twin.Tm.Should().Be(bp.Tm);
                twin.GcContent.Should().Be(bp.GcContent);
                twin.Score.Should().Be(bp.Score, because: "the surviving probes are byte-identical windows, so every metric is preserved exactly");
            }
        }
    }

    #endregion

    #endregion

    #region PROBE-VALID-001 — probe validation

    // ─────────────────────────────────────────────────────────────────────────
    // Unit: PROBE-VALID-001 — hybridization probe validation (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 25.
    //
    // API under test (ProbeDesigner.ValidateProbe):
    //   ValidateProbe(probe, references, maxMismatches = 3, selfComplementarityThreshold = 0.3)
    //     → ProbeValidation(IsValid, SpecificityScore, OffTargetHits, SelfComplementarity,
    //                       HasSecondaryStructure, Issues).
    //   offTargetHits = Σ_refs |{ i : Hamming(ref[i..i+|probe|], probe) ≤ maxMismatches }|
    //                   (ungapped fixed-length substitution scan, both sides uppercased).
    //   IsValid = issues.Count == 0  ∨  (offTargetHits ≤ 1 ∧ selfComp ≤ 0.4), where issues
    //   holds an off-target item iff offTargetHits > 1, a self-comp item iff
    //   selfComp > selfComplementarityThreshold, and a structure item iff hasStructure.
    //
    //   Algebra (definition, NOT observed output): substituting the issue conditions,
    //     IsValid = (offTargetHits ≤ 1) ∧ R,   R = ((selfComp ≤ thr ∧ ¬hasStructure) ∨ selfComp ≤ 0.4)
    //   and R is INDEPENDENT of maxMismatches (selfComp / hasStructure read the probe alone).
    //
    // Relation DIRECTIONS (derived from that algebra and the spec, NOT from output):
    //   Source: docs/algorithms/MolTools/Probe_Validation.md §2.2 (specificity map),
    //           §2.4 (INV-01..04), §4.1 step 5 (the IsValid rule), §5.2.
    //
    //   • MON (lower specificity threshold → more pass): `maxMismatches` is the off-target
    //     detection stringency — the mismatch tolerance below which a window counts as a
    //     cross-hybridization site. FindApproximateMatches is monotone in it:
    //     matches(k) ⊆ matches(k+1) (a window within k mismatches is within k+1), so
    //     offTargetHits(k) is NON-DECREASING in k. Since IsValid = (offTargetHits ≤ 1) ∧ R
    //     with R fixed, IsValid(k) is NON-INCREASING: lowering the threshold can only ADD
    //     passing probes (the passing set is downward-closed in k). The dual SpecificityScore
    //     (1 / hits for hits > 1) is likewise non-increasing as hits grow.
    //
    //   • INV (same input → same result): ValidateProbe is a pure function whose result
    //     depends only on the MULTISET of references (offTargetHits is a SUM over them) and
    //     on case-folded sequences (probe and each reference are uppercased). Hence the full
    //     record is invariant under (a) repeated calls, (b) permuting the reference order,
    //     and (c) changing the case of the probe and/or references.
    //
    // Construction: a "clean" probe has selfComp ≤ threshold AND no secondary structure, so
    // R = true and IsValid reduces to (offTargetHits ≤ 1) — isolating the maxMismatches MR.
    // References are exact-length copies of the probe mutated at a fixed number of positions
    // (via the same A↔T / C↔G flip used above), so each reference contributes one window at
    // a KNOWN Hamming distance and the hit count crosses 1 as the tolerance rises. The clean
    // probes' R-precondition is asserted from the returned fields, so the setup is self-checking.
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>Default self-complementarity issue threshold used by ValidateProbe.</summary>
    private const double SelfCompThreshold = 0.3;

    private static ProbeDesigner.ProbeValidation Validate(
        string probe, IEnumerable<string> references, int maxMismatches, double selfThr = SelfCompThreshold)
        => ProbeDesigner.ValidateProbe(probe, references, maxMismatches, selfThr);

    /// <summary>
    /// 20-nt probes that are NOT self-complementary and form no hairpin, so the validity
    /// rule reduces to (offTargetHits ≤ 1) — letting the maxMismatches relation be tested
    /// in isolation. (Each probe's R-precondition is verified in-test from the result.)
    /// </summary>
    private static IEnumerable<string> CleanProbes()
    {
        yield return "AAAAAAAAAAGGGGGGGGGG";
        yield return "AAAAAAAAAACCCCCCCCCC";
        yield return "TTTTTTTTTTCCCCCCCCCC";
    }

    /// <summary>
    /// References for a clean probe: an EXACT copy (0 mismatches) plus copies at 2, 4 and 6
    /// mismatches. Each is the probe's length, so it contributes exactly one window whose
    /// Hamming distance to the probe is the stated value; the hit count therefore steps up as
    /// the mismatch tolerance rises, crossing the validity boundary (offTargetHits ≤ 1).
    /// </summary>
    private static List<string> MismatchLadder(string probe) => new()
    {
        probe,                                  // 0 mismatches (exact on-target)
        MutateAt(probe, 0, 1),                  // 2 mismatches
        MutateAt(probe, 0, 1, 2, 3),            // 4 mismatches
        MutateAt(probe, 0, 1, 2, 3, 4, 5),      // 6 mismatches
    };

    private static void AssertSameValidation(ProbeDesigner.ProbeValidation a, ProbeDesigner.ProbeValidation b, string why)
    {
        a.IsValid.Should().Be(b.IsValid, because: why);
        a.SpecificityScore.Should().Be(b.SpecificityScore, because: why);
        a.OffTargetHits.Should().Be(b.OffTargetHits, because: why);
        a.SelfComplementarity.Should().Be(b.SelfComplementarity, because: why);
        a.HasSecondaryStructure.Should().Be(b.HasSecondaryStructure, because: why);
        a.Issues.Should().Equal(b.Issues, because: why);
    }

    #region MON — lowering the mismatch (specificity) threshold never removes a passing probe

    [Test]
    [Description("MON: along an increasing maxMismatches chain the off-target hit count is non-decreasing and IsValid is non-increasing — lowering the threshold can only ADD passing probes, never remove one.")]
    public void ValidateProbe_RaisingMismatchTolerance_HitsNonDecreasing_ValidityNonIncreasing()
    {
        int[] toleranceChain = { 0, 1, 2, 3, 4, 5, 6 };

        foreach (var probe in CleanProbes())
        {
            var references = MismatchLadder(probe);

            // Precondition: the probe is "clean" (R = true), so IsValid ⇔ offTargetHits ≤ 1.
            var clean = Validate(probe, references, maxMismatches: 0);
            clean.SelfComplementarity.Should().BeLessThanOrEqualTo(SelfCompThreshold,
                because: "the construction requires a non-self-complementary probe so validity tracks only the off-target count");
            clean.HasSecondaryStructure.Should().BeFalse(
                because: "the construction requires a probe with no hairpin so validity tracks only the off-target count");

            int previousHits = -1;
            bool? previousValid = null;
            bool sawValid = false, sawInvalid = false;

            foreach (int k in toleranceChain)
            {
                var result = Validate(probe, references, maxMismatches: k);

                result.OffTargetHits.Should().BeGreaterThanOrEqualTo(previousHits,
                    because: $"raising the tolerance to {k} can only ADD approximate matches (matches(k) ⊆ matches(k+1)), so the hit count is non-decreasing");

                // IsValid ⇔ offTargetHits ≤ 1 for a clean probe — the algebra, checked exactly.
                result.IsValid.Should().Be(result.OffTargetHits <= 1,
                    because: "for a clean probe (R = true) the validity rule reduces to offTargetHits ≤ 1");

                if (previousValid is not null)
                    (result.IsValid && !previousValid.Value).Should().BeFalse(
                        because: $"a higher mismatch tolerance ({k}) never turns an invalid probe valid — validity is non-increasing in the threshold");

                sawValid |= result.IsValid;
                sawInvalid |= !result.IsValid;
                previousHits = result.OffTargetHits;
                previousValid = result.IsValid;
            }

            sawValid.Should().BeTrue(because: "at low tolerance only the exact on-target is found (1 hit) ⇒ the probe passes");
            sawInvalid.Should().BeTrue(because: "at high tolerance the near-duplicates are also found (>1 hit) ⇒ the probe fails — the relation is exercised, not vacuous");
        }
    }

    [Test]
    [Description("MON (set): the set of probes that PASS validation is downward-closed in maxMismatches — lowering the threshold yields a superset of passing probes (count non-decreasing).")]
    public void ValidateProbe_LoweringThreshold_PassingSetIsSuperset()
    {
        // Each probe carries its own mismatch ladder; pass/fail is decided per probe.
        var cases = CleanProbes().Select(p => (Probe: p, Refs: MismatchLadder(p))).ToList();
        int[] descendingThresholds = { 6, 4, 3, 2, 1, 0 };

        HashSet<string>? higherSet = null;
        foreach (int k in descendingThresholds)
        {
            var passing = cases.Where(c => Validate(c.Probe, c.Refs, k).IsValid).Select(c => c.Probe).ToHashSet();

            if (higherSet is not null)
                passing.IsSupersetOf(higherSet).Should().BeTrue(
                    because: $"lowering the threshold to {k} keeps every probe that already passed and may add more — the passing set is a superset");

            higherSet = passing;
        }
    }

    [Test]
    [Description("MON dual: the specificity score is non-increasing as the off-target hit count grows (1 hit → 1.0, N hits → 1/N).")]
    public void ValidateProbe_MoreOffTargetHits_SpecificityScoreNonIncreasing()
    {
        foreach (var probe in CleanProbes())
        {
            var references = MismatchLadder(probe);
            double previousScore = double.PositiveInfinity;

            foreach (int k in new[] { 0, 1, 2, 3, 4, 5, 6 })
            {
                var result = Validate(probe, references, k);

                result.SpecificityScore.Should().BeLessThanOrEqualTo(previousScore,
                    because: $"as the tolerance rises to {k} the hit count grows and specificity (1/hits for hits>1) is non-increasing");
                if (result.OffTargetHits >= 1)
                    result.SpecificityScore.Should().BeApproximately(1.0 / result.OffTargetHits, 1e-12,
                        because: "the spec maps hits>0 to 1/hits (and 1 hit to 1.0), independent of which references produced them");

                previousScore = result.SpecificityScore;
            }
        }
    }

    #endregion

    #region INV — same input → same result (deterministic, reference-order- and case-independent)

    [Test]
    [Description("INV: ValidateProbe is deterministic — repeated calls on identical input return an identical validation record.")]
    public void ValidateProbe_SameInput_IsDeterministic()
    {
        foreach (var probe in CleanProbes())
        {
            var references = MismatchLadder(probe);
            foreach (int k in new[] { 0, 2, 4 })
            {
                var first = Validate(probe, references, k);
                for (int i = 0; i < 4; i++)
                    AssertSameValidation(Validate(probe, references, k), first,
                        "ValidateProbe is a pure function of (probe, references, params), so repeated calls must agree exactly");
            }
        }
    }

    [Test]
    [Description("INV: permuting the reference order leaves the validation record unchanged — offTargetHits is a sum over references, which is order-independent.")]
    public void ValidateProbe_ReferenceOrderPermutation_ResultUnchanged()
    {
        foreach (var probe in CleanProbes())
        {
            var references = MismatchLadder(probe);
            var reversed = Enumerable.Reverse(references).ToList();
            var rotated = references.Skip(1).Concat(references.Take(1)).ToList();

            foreach (int k in new[] { 0, 2, 4, 6 })
            {
                var baseline = Validate(probe, references, k);
                AssertSameValidation(Validate(probe, reversed, k), baseline,
                    "the total hit count sums over references regardless of order, and every other field reads the probe alone");
                AssertSameValidation(Validate(probe, rotated, k), baseline,
                    "rotating the reference list cannot change the summed hit count or the probe-only metrics");
            }
        }
    }

    [Test]
    [Description("INV: ValidateProbe is case-insensitive — lower/mixed casing the probe and references gives the same record, since both sides are uppercased before analysis.")]
    public void ValidateProbe_CaseFolding_ResultUnchanged()
    {
        foreach (var probe in CleanProbes())
        {
            var references = MismatchLadder(probe);
            foreach (int k in new[] { 0, 2, 4 })
            {
                var baseline = Validate(probe, references, k);

                AssertSameValidation(
                    Validate(probe.ToLowerInvariant(), references.Select(r => r.ToLowerInvariant()).ToList(), k),
                    baseline,
                    "the probe and each reference are uppercased before matching, so casing carries no information");

                var mixedRefs = references
                    .Select(r => new string(r.Select((ch, i) => i % 2 == 0 ? char.ToLowerInvariant(ch) : ch).ToArray()))
                    .ToList();
                AssertSameValidation(Validate(probe, mixedRefs, k), baseline,
                    "case folding is total, so any case pattern of the same letters yields the same validation");
            }
        }
    }

    #endregion

    #endregion

    #region RESTR-FIND-001 — restriction enzyme site finding

    // ─────────────────────────────────────────────────────────────────────────
    // Unit: RESTR-FIND-001 — restriction enzyme site finding (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 26.
    //
    // API under test (RestrictionAnalyzer.FindSitesCore via the FindSites overloads):
    //   FindSites(seq, enzymeName)        → sites of ONE enzyme.
    //   FindSites(seq, params names[])    → Σ over names of FindSites(seq, name).
    //   FindAllSites(seq)                 → Σ over the whole enzyme database.
    //   A site is found at every forward-strand match of the recognition pattern (IUPAC-
    //   aware), AND at every match on the reverse complement, whose coordinate is mapped
    //   BACK to the forward strand: for a reverse hit at revComp index i,
    //   forwardPos = |seq| − i − L. So RestrictionSite.Position and .CutPosition are ALWAYS
    //   forward-strand coordinates (0-based), for both strands. A palindromic recognition
    //   site therefore yields TWO entries (forward + reverse) at the same Position.
    //
    // A site's run-order-independent identity is
    //   (Position, EnzymeName, IsForwardStrand, CutPosition, RecognizedSequence).
    //
    // Relation DIRECTIONS (derived from the search definition, NOT from observed output):
    //   Source: docs/algorithms/MolTools/Restriction_Site_Detection.md; RestrictionAnalyzer.cs
    //   (FindSitesCore: forward scan + reverse-complement scan with forwardPos remap).
    //
    //   • MON (more enzymes → ≥ total sites): FindSites over a name set is the UNION of the
    //     per-enzyme site sets (the params overload concatenates them; FindAllSites unions
    //     the whole database). Enlarging the enzyme set can only ADD sites, so the site set
    //     is a SUPERSET and the count is non-decreasing; every subset's sites ⊆ FindAllSites.
    //
    //   • SHIFT (prepend flank shifts positions): Position/CutPosition are forward-strand
    //     coordinates for BOTH strands. Prepending a flank F that creates NO new site (and
    //     destroys none) moves every forward match from i to i+|F|; a reverse match keeps
    //     its place in revComp(F+seq) = revComp(seq)+revComp(F) (the revComp(seq) prefix is
    //     unchanged), so forwardPos = |F+seq|−i−L = old forwardPos + |F|. Thus EVERY site's
    //     Position and CutPosition advance by exactly |F|, with enzyme/strand/recognized
    //     sequence preserved.
    //
    //   • INV (non-site append → same sites): appending X leaves every forward match inside
    //     the original body byte-identical (same Position), and a reverse match maps back to
    //     forwardPos = |seq+X| − (i+|X|) − L = old forwardPos (forward coordinates are
    //     anchored at the 5' end, which the append does not move). So an append that creates
    //     no new site leaves the site set EXACTLY unchanged.
    //
    // Neutral flank/append construction (so SHIFT/INV preserve sites exactly): for a
    // PALINDROMIC pattern P a homopolymer of base b ∉ {P[0], P[last]} introduces no match —
    // internally (a homopolymer can't contain a ≥2-letter pattern) and across either
    // junction on either strand. Forward junction safety needs P[last] ≠ b (a suffix of P
    // overlapping the X-run would have to be all-b); reverse junction safety needs
    // P[0] ≠ comp(b); for a palindrome comp(P[0]) = P[last], so both reduce to
    // b ∉ {P[0], P[last]}. Each test also GUARDS this empirically (the flank alone yields no
    // site) and the set-equality assertions would catch any stray junction site.
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>Run-order-independent identity of a restriction site (all forward-strand coordinates).</summary>
    private static (int Position, string Enzyme, bool Fwd, int Cut, string Rec) RestrSiteId(RestrictionSite s)
        => (s.Position, s.Enzyme.Name, s.IsForwardStrand, s.CutPosition, s.RecognizedSequence);

    private static HashSet<(int, string, bool, int, string)> RestrSiteSet(string seq, params string[] enzymes)
        => RestrictionAnalyzer.FindSites(new DnaSequence(seq), enzymes).Select(RestrSiteId).ToHashSet();

    /// <summary>Palindromic, pure-ACGT Type-II enzymes — each recognition site yields a forward + reverse entry.</summary>
    private static readonly string[] PalindromicEnzymes = { "EcoRI", "BamHI", "HindIII", "XhoI", "SalI", "NcoI" };

    /// <summary>
    /// A homopolymer base provably neutral for a palindromic pattern: b ∉ {P[0], P[last]}
    /// adds no site internally or across either junction on either strand.
    /// </summary>
    private static char SafeFlankBase(string palindromicPattern)
    {
        foreach (char b in "ACGT")
            if (b != palindromicPattern[0] && b != palindromicPattern[^1])
                return b;
        throw new InvalidOperationException($"No neutral flank base for pattern '{palindromicPattern}'.");
    }

    /// <summary>A body embedding two copies of <paramref name="pattern"/> separated by a neutral spacer.</summary>
    private static string RestrBody(string pattern)
    {
        const string spacer = "ACTGAC";
        return spacer + pattern + spacer + pattern + spacer;
    }

    /// <summary>A body embedding one site each for EcoRI, BamHI, HindIII and XhoI.</summary>
    private const string RestrMultiBody =
        "ACTGAC" + "GAATTC" + "ACTGAC" + "GGATCC" + "ACTGAC" + "AAGCTT" + "ACTGAC" + "CTCGAG" + "ACTGAC";

    #region MON — enlarging the enzyme set never removes sites (count non-decreasing, superset)

    [Test]
    [Description("MON: along a growing chain of enzyme sets the site set grows monotonically — each smaller set's sites are a subset of the larger set's, total count non-decreasing.")]
    public void FindSites_AddingEnzymes_YieldsSuperset_CountNonDecreasing()
    {
        string[][] enzymeChain =
        {
            new[] { "EcoRI" },
            new[] { "EcoRI", "BamHI" },
            new[] { "EcoRI", "BamHI", "HindIII" },
            new[] { "EcoRI", "BamHI", "HindIII", "XhoI" },
        };

        HashSet<(int, string, bool, int, string)>? previous = null;
        foreach (var enzymes in enzymeChain)
        {
            var sites = RestrSiteSet(RestrMultiBody, enzymes);

            if (previous is not null)
            {
                sites.IsSupersetOf(previous).Should().BeTrue(
                    because: $"FindSites over {enzymes.Length} enzymes is the UNION of the per-enzyme site sets, " +
                             "so adding an enzyme keeps every existing site and can only add more");
                sites.Count.Should().BeGreaterThanOrEqualTo(previous.Count,
                    because: "enlarging the enzyme set unions in more sites — the total count is non-decreasing");
            }

            previous = sites;
        }
    }

    [Test]
    [Description("MON: every subset of enzymes finds a subset of the sites that FindAllSites (the whole database) reports.")]
    public void FindAllSites_IsSupersetOfAnyEnzymeSubset()
    {
        var allForward = RestrictionAnalyzer.FindAllSites(new DnaSequence(RestrMultiBody))
            .Select(RestrSiteId).ToHashSet();

        foreach (var enzymes in new[]
                 {
                     new[] { "EcoRI" },
                     new[] { "BamHI", "HindIII" },
                     new[] { "EcoRI", "BamHI", "HindIII", "XhoI" },
                 })
        {
            var subset = RestrSiteSet(RestrMultiBody, enzymes);
            subset.Should().NotBeEmpty(because: "the multi-site body contains a site for each of these enzymes");
            allForward.IsSupersetOf(subset).Should().BeTrue(
                because: "FindAllSites scans the full enzyme database, so it reports every site any subset of those enzymes finds");
        }
    }

    #endregion

    #region SHIFT — prepending a neutral flank advances every site by exactly the flank length

    [Test]
    [Description("SHIFT: prepending a neutral homopolymer flank F advances every site's Position and CutPosition by exactly |F| (both strands), preserving enzyme, strand and recognized sequence.")]
    public void FindSites_PrependNeutralFlank_ShiftsAllSitesByFlankLength()
    {
        foreach (var name in PalindromicEnzymes)
        {
            string pattern = RestrictionAnalyzer.GetEnzyme(name)!.RecognitionSequence;
            char b = SafeFlankBase(pattern);
            string body = RestrBody(pattern);

            var baseSites = RestrictionAnalyzer.FindSites(new DnaSequence(body), name).ToList();
            baseSites.Should().NotBeEmpty(because: $"the body embeds two {name} recognition sites");

            foreach (int flankLen in new[] { 1, 5, 17 })
            {
                string flank = new string(b, flankLen);
                RestrictionAnalyzer.FindSites(flank, name).Should().BeEmpty(
                    because: $"a '{b}'-homopolymer flank contains no {name} site (b ∉ {{{pattern[0]},{pattern[^1]}}})");

                var shifted = RestrictionAnalyzer.FindSites(new DnaSequence(flank + body), name).Select(RestrSiteId).ToHashSet();
                var expected = baseSites
                    .Select(s => (s.Position + flankLen, s.Enzyme.Name, s.IsForwardStrand, s.CutPosition + flankLen, s.RecognizedSequence))
                    .ToHashSet();

                shifted.SetEquals(expected).Should().BeTrue(
                    because: $"Position and CutPosition are forward-strand coordinates for both strands, so a length-{flankLen} neutral prefix " +
                             "advances every site by exactly that amount while preserving enzyme, strand and recognized sequence");
            }
        }
    }

    #endregion

    #region INV — appending a neutral (non-site) region leaves the site set exactly unchanged

    [Test]
    [Description("INV: appending a neutral homopolymer region (no new site on either strand) leaves the restriction-site set EXACTLY unchanged — forward coordinates are anchored at the 5' end.")]
    public void FindSites_AppendNeutralRegion_SiteSetUnchanged()
    {
        foreach (var name in PalindromicEnzymes)
        {
            string pattern = RestrictionAnalyzer.GetEnzyme(name)!.RecognitionSequence;
            char b = SafeFlankBase(pattern);
            string body = RestrBody(pattern);

            var baseSites = RestrictionAnalyzer.FindSites(new DnaSequence(body), name).Select(RestrSiteId).ToHashSet();
            baseSites.Should().NotBeEmpty(because: $"the body embeds two {name} recognition sites");

            foreach (int appendLen in new[] { 1, 5, 17 })
            {
                string ext = new string(b, appendLen);
                RestrictionAnalyzer.FindSites(ext, name).Should().BeEmpty(
                    because: $"a '{b}'-homopolymer append contains no {name} site");

                var appended = RestrictionAnalyzer.FindSites(new DnaSequence(body + ext), name).Select(RestrSiteId).ToHashSet();

                appended.SetEquals(baseSites).Should().BeTrue(
                    because: "appending a non-site region downstream creates no new match and shifts no forward-anchored coordinate, " +
                             "so the site set is preserved exactly (same positions, cut sites, strands, recognized sequences)");
            }
        }
    }

    #endregion

    #endregion

    #region RESTR-DIGEST-001 — restriction digest simulation

    // ─────────────────────────────────────────────────────────────────────────
    // Unit: RESTR-DIGEST-001 — restriction digest simulation (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 27.
    //
    // API under test (RestrictionAnalyzer.Digest, LINEAR topology = default overload):
    //   Digest(seq, params enzymeNames)
    //     1. collect DISTINCT forward-strand cut positions over all enzymes (a SortedSet,
    //        forward strand only to avoid double-counting palindromes).
    //     2. zero cuts ⇒ a single fragment = the whole sequence.
    //     3. otherwise fragments span consecutive boundaries of {0} ∪ cuts ∪ {|seq|};
    //        zero-length pieces (a cut at 0 or |seq|) are dropped.
    //   k distinct interior cut positions ⇒ k+1 fragments; the fragments TILE the sequence.
    //
    // Relation DIRECTIONS (derived from the partition definition, NOT from observed output):
    //   Source: docs/algorithms/MolTools/Restriction_Digest_Simulation.md; RestrictionAnalyzer.Digest.
    //
    //   • COMP (0 sites → 1 fragment = full seq): with no cut position the source returns the
    //     entire molecule as one fragment (Sequence = seq, Length = |seq|, Start = 0, no
    //     flanking enzymes). The documented endpoint.
    //   • MON (more enzymes → ≥ fragments): cut positions are a SET unioned over enzymes, so
    //     adding an enzyme can only ADD cut positions; each new interior cut splits one
    //     fragment into two and coincident/boundary cuts leave the count unchanged — the
    //     fragment count is NON-DECREASING in the enzyme set.
    //   • INV (fragment sum = seq length): the cut boundaries partition [0,|seq|], so the
    //     fragment lengths telescope to exactly |seq| and the fragments concatenate back to
    //     the original sequence — for ANY enzyme set, with or without cuts.
    // ─────────────────────────────────────────────────────────────────────────

    private static List<DigestFragment> Digest(string body, params string[] enzymes)
        => RestrictionAnalyzer.Digest(new DnaSequence(body), enzymes).ToList();

    /// <summary>Bodies for digest relations: multi-site, single-enzyme, no-site, and fixed-seed random.</summary>
    private static IEnumerable<string> DigestBodies()
    {
        yield return RestrMultiBody;                 // sites for EcoRI, BamHI, HindIII, XhoI
        yield return RestrBody("GAATTC");            // two EcoRI sites
        yield return "AAAAAAAAAAAAAAAAAAAA";          // no site for the enzymes used below
        yield return RandomDna(80);
        yield return RandomDna(120);
    }

    private static readonly string[][] DigestEnzymeSets =
    {
        new[] { "EcoRI" },
        new[] { "EcoRI", "BamHI" },
        new[] { "EcoRI", "BamHI", "HindIII", "XhoI" },
    };

    #region COMP — a sequence with no recognition site digests to a single full-length fragment

    [Test]
    [Description("COMP: digesting a sequence that contains no recognition site yields exactly ONE fragment equal to the whole input sequence.")]
    public void Digest_NoRecognitionSite_ReturnsSingleFullLengthFragment()
    {
        const string noSiteBody = "AAAAAAAAAAAAAAAAAAAA"; // contains no GAATTC / GGATCC / AAGCTT

        foreach (var enzymes in new[] { new[] { "EcoRI" }, new[] { "EcoRI", "BamHI", "HindIII" } })
        {
            var fragments = Digest(noSiteBody, enzymes);

            fragments.Should().ContainSingle(because: "with no cut site the molecule is returned intact as one fragment");
            var only = fragments[0];
            only.Sequence.Should().Be(noSiteBody, because: "the single fragment is the entire input sequence");
            only.Length.Should().Be(noSiteBody.Length);
            only.StartPosition.Should().Be(0);
            only.LeftEnzyme.Should().BeNull(because: "an uncut fragment has no flanking enzyme on the left");
            only.RightEnzyme.Should().BeNull(because: "an uncut fragment has no flanking enzyme on the right");
        }
    }

    #endregion

    #region MON — adding enzymes never reduces the fragment count

    [Test]
    [Description("MON: along a growing enzyme chain the fragment count is non-decreasing — each added cutter only splits fragments, never merges them.")]
    public void Digest_AddingEnzymes_FragmentCountNonDecreasing()
    {
        string[][] chain =
        {
            new[] { "EcoRI" },
            new[] { "EcoRI", "BamHI" },
            new[] { "EcoRI", "BamHI", "HindIII" },
            new[] { "EcoRI", "BamHI", "HindIII", "XhoI" },
        };

        int previousCount = -1;
        int firstCount = -1;
        foreach (var enzymes in chain)
        {
            int count = Digest(RestrMultiBody, enzymes).Count;
            if (firstCount < 0) firstCount = count;

            count.Should().BeGreaterThanOrEqualTo(previousCount,
                because: $"adding cutters to {string.Join("+", enzymes)} unions in more cut positions, each splitting a fragment — the count cannot drop");
            previousCount = count;
        }

        previousCount.Should().BeGreaterThan(firstCount,
            because: "the body carries a distinct site per enzyme, so the four-enzyme digest yields strictly more fragments than the single-enzyme one — the relation is exercised");
    }

    #endregion

    #region INV — fragments partition the sequence: lengths sum to |seq| and concatenate back

    [Test]
    [Description("INV: for any enzyme set the fragment lengths sum to exactly the sequence length, and the fragments concatenate back to the original sequence — a digest is a partition.")]
    public void Digest_Fragments_TileTheSequence()
    {
        foreach (var body in DigestBodies())
        {
            foreach (var enzymes in DigestEnzymeSets)
            {
                var fragments = Digest(body, enzymes);

                fragments.Should().NotBeEmpty(because: "a linear digest always yields at least the whole molecule");
                fragments.Should().OnlyContain(f => f.Length > 0, because: "zero-length boundary pieces are dropped");

                fragments.Sum(f => f.Length).Should().Be(body.Length,
                    because: "the cut boundaries partition [0,|seq|], so the fragment lengths telescope to the full length");

                string reconstructed = string.Concat(fragments.OrderBy(f => f.StartPosition).Select(f => f.Sequence));
                reconstructed.Should().Be(body,
                    because: "the fragments are the consecutive pieces between cuts, so concatenating them in position order rebuilds the input exactly");
            }
        }
    }

    [Test]
    [Description("INV: the fragment count equals the number of distinct forward-strand cut positions plus one (interior cuts), confirming the k cuts → k+1 fragments partition rule.")]
    public void Digest_FragmentCount_EqualsInteriorCutsPlusOne()
    {
        foreach (var body in DigestBodies())
        {
            foreach (var enzymes in DigestEnzymeSets)
            {
                var forwardCuts = enzymes
                    .SelectMany(e => RestrictionAnalyzer.FindSites(new DnaSequence(body), e))
                    .Where(s => s.IsForwardStrand)
                    .Select(s => s.CutPosition)
                    .Where(c => c > 0 && c < body.Length)   // interior cuts produce real fragment boundaries
                    .Distinct()
                    .Count();

                Digest(body, enzymes).Count.Should().Be(forwardCuts + 1,
                    because: "k distinct interior forward-strand cut positions split a linear molecule into k+1 fragments");
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RESTR-FILTER-001 — restriction-enzyme catalogue filtering (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 224.
    //
    // API under test (RestrictionAnalyzer.GetEnzymesByCutLength / GetBluntCutters / GetStickyCutters):
    //   Select enzymes from the fixed catalogue by recognition-site length (inclusive [min,max]) or by
    //   end type (blunt vs sticky).
    //
    // Relations (derived from the predicate-filter semantics, NOT from output):
    //   • SUB (filtered ⊆ all): every filter returns a subset of the full catalogue; blunt and sticky
    //         partition it (disjoint, union = all).
    //   • MON (stricter criteria ⇒ subset): narrowing the inclusive length window yields a subset of
    //         the wider window's enzymes.
    // ───────────────────────────────────────────────────────────────────────────

    #region RESTR-FILTER-001 — Helpers

    private static HashSet<string> Names(IEnumerable<RestrictionEnzyme> enzymes) =>
        enzymes.Select(e => e.Name).ToHashSet();

    private static readonly HashSet<string> AllEnzymeNames = RestrictionAnalyzer.Enzymes.Values.Select(e => e.Name).ToHashSet();

    #endregion

    #region RESTR-FILTER-001 SUB — every filtered set is a subset of the catalogue

    [Test]
    [Description("SUB: each catalogue filter returns a subset of all enzymes; blunt and sticky partition the catalogue.")]
    public void RestrictionFilter_FilteredSets_AreSubsetsOfCatalogue()
    {
        Names(RestrictionAnalyzer.GetEnzymesByCutLength(4, 8)).Should().BeSubsetOf(AllEnzymeNames,
            because: "a length filter selects from the catalogue, adding nothing");

        var blunt = Names(RestrictionAnalyzer.GetBluntCutters());
        var sticky = Names(RestrictionAnalyzer.GetStickyCutters());

        blunt.Should().BeSubsetOf(AllEnzymeNames);
        sticky.Should().BeSubsetOf(AllEnzymeNames);
        blunt.Should().NotIntersectWith(sticky, because: "an enzyme is either blunt or sticky, never both");
        blunt.Union(sticky).Should().BeEquivalentTo(AllEnzymeNames, because: "blunt and sticky exhaust the catalogue");
    }

    #endregion

    #region RESTR-FILTER-001 MON — narrowing the length window yields a subset

    [Test]
    [Description("MON: narrowing the inclusive recognition-length window selects a subset of the wider window's enzymes.")]
    public void RestrictionFilter_NarrowerLengthWindow_IsSubset()
    {
        var wide = Names(RestrictionAnalyzer.GetEnzymesByCutLength(4, 8));
        var mid = Names(RestrictionAnalyzer.GetEnzymesByCutLength(5, 7));
        var narrow = Names(RestrictionAnalyzer.GetEnzymesByCutLength(6, 6));

        narrow.Should().BeSubsetOf(mid, because: "[6,6] ⊆ [5,7] as length predicates");
        mid.Should().BeSubsetOf(wide, because: "[5,7] ⊆ [4,8] as length predicates");

        // The single-length overload agrees with the degenerate window [n,n].
        Names(RestrictionAnalyzer.GetEnzymesByCutLength(6))
            .Should().BeEquivalentTo(narrow, because: "GetEnzymesByCutLength(6) is the window [6,6]");

        // Non-vacuous: the canonical 6-cutter EcoRI (GAATTC) lands in every window that includes 6.
        narrow.Should().Contain("EcoRI", because: "EcoRI has a 6-bp recognition site");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PRIMER-NNTM-001 — SantaLucia nearest-neighbour, salt-corrected Tm (MolTools)
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (SantaLucia 1998 PNAS 95:1460; SantaLucia & Hicks 2004; Owczarzy 2004 Biochemistry
    //   43:3537; docs/algorithms/MolTools/NearestNeighbor_Salt_Corrected_Tm.md):
    //   CalculateMeltingTemperatureNN sums per-dinucleotide ΔH°/ΔS° nearest-neighbour parameters
    //   plus terminal-A·T and (self-complementary) symmetry corrections, computes the bimolecular
    //   Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) and applies the Owczarzy (2004) monovalent-salt
    //   correction. Two metamorphic relations (checklist row 240):
    //
    //   • MON (raising monovalent salt raises Tm): higher [Na⁺] stabilises the duplex; the
    //     Owczarzy 2004 correction is monotone increasing in Tm over the physiological salt range,
    //     so Tm is strictly increasing in [Na⁺] for a fixed primer.
    //   • INV (reverse-complement has equal duplex Tm): the SantaLucia parameters are
    //     strand-symmetric — the multiset of nearest-neighbour stacks of revcomp(S) equals that of
    //     S, and revcomp preserves each terminus's A·T status — so ΔH°, ΔS° and hence Tm are
    //     identical for a sequence and its reverse complement (the same duplex read from the other
    //     strand).
    //
    // API under test: PrimerDesigner.CalculateMeltingTemperatureNN.

    #region PRIMER-NNTM-001 — nearest-neighbour salt-corrected Tm

    // Deterministic, non-self-complementary primers spanning a range of GC content and length.
    // (Literal fixtures — NOT drawn from the shared Rng, which would perturb every other fixture.)
    private static readonly string[] NnTmPrimers =
    {
        "CAGGTGGCACCTTAACG",
        "GCTAGCATCGGATCCAA",
        "TTGACCAGTCCATGGCA",
        "ATGCGGTCAATTGCAACGT",
        "GGGCGCGGCACCGTCCA",
    };

    [Test]
    [Description("MON: higher monovalent [Na⁺] stabilises the duplex; the Owczarzy 2004 salt correction is monotone, so the nearest-neighbour Tm strictly increases with [Na⁺] for a fixed primer.")]
    public void NnTm_RaisingMonovalentSalt_RaisesTm()
    {
        double[] sodium = { 0.02, 0.05, 0.10, 0.25, 0.50 }; // 20 mM → 500 mM, ascending

        foreach (string primer in NnTmPrimers)
        {
            double previous = double.NegativeInfinity;
            foreach (double na in sodium)
            {
                double tm = PrimerDesigner.CalculateMeltingTemperatureNN(primer, sodiumMolar: na);
                tm.Should().NotBe(double.NaN, because: $"'{primer}' is a valid ACGT primer");
                tm.Should().BeGreaterThan(previous,
                    because: $"raising [Na⁺] to {na} M stabilises the duplex, so the NN Tm of '{primer}' strictly increases");
                previous = tm;
            }
        }
    }

    [Test]
    [Description("INV: the SantaLucia NN parameters are strand-symmetric and revcomp preserves each terminus's A·T status, so a primer and its reverse complement have identical ΔH°/ΔS° and hence identical NN Tm.")]
    public void NnTm_ReverseComplement_HasEqualDuplexTm()
    {
        double[] sodium = { 0.05, 0.50 };

        foreach (string primer in NnTmPrimers)
        {
            string rc = DnaSequence.GetReverseComplementString(primer);

            // Non-vacuity: the primer is not its own reverse complement, so this is a real transform.
            rc.Should().NotBe(primer, because: $"'{primer}' must be non-self-complementary for the INV to be non-trivial");

            foreach (double na in sodium)
            {
                double tm = PrimerDesigner.CalculateMeltingTemperatureNN(primer, sodiumMolar: na);
                double tmRc = PrimerDesigner.CalculateMeltingTemperatureNN(rc, sodiumMolar: na);

                tm.Should().BeGreaterThan(0.0, because: $"'{primer}' has a well-defined positive Tm");
                tmRc.Should().BeApproximately(tm, 1e-9,
                    because: $"the reverse complement is the same duplex read from the other strand, so its NN Tm equals that of '{primer}' at [Na⁺]={na} M");
            }
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PRIMER-HAIRPIN-001 — most-stable DNA hairpin ΔG (MolTools)
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (SantaLucia 1998; SantaLucia & Hicks 2004 Table 1 + Table 4;
    //   docs/algorithms/MolTools/DNA_Hairpin_Folding_Tm.md):
    //   FindMostStableHairpin finds the minimum-ΔG°37 intramolecular hairpin (one Watson-Crick
    //   stem closing one loop): ΔG°37 = Σ stem NN stacks (each ≤ 0) + loop initiation (≥ 0). It
    //   returns null when no stem of ≥ 2 bp can close a loop of ≥ 3 nt. Two metamorphic relations
    //   (checklist row 241):
    //
    //   • MON (lengthening a complementary stem lowers ΔG): each extra base pair adds one more
    //     stabilising (negative) nearest-neighbour stack while the loop term is unchanged, so the
    //     most-stable hairpin's ΔG°37 is strictly DECREASING (more negative) as the stem grows.
    //   • INV (no stem possible → no hairpin): a sequence that admits no ≥2-bp Watson-Crick stem
    //     closing a ≥3-nt loop (a homopolymer, or a lone WC pair) has no hairpin — the method
    //     returns null — and that holds however long the non-pairing sequence is.
    //
    // API under test: PrimerDesigner.FindMostStableHairpin (HairpinResult?).

    #region PRIMER-HAIRPIN-001 — most-stable hairpin ΔG

    [Test]
    [Description("MON: each added stem base pair contributes one more stabilising NN stack while the loop term is fixed, so lengthening a complementary stem strictly lowers (makes more negative) the hairpin ΔG°37.")]
    public void Hairpin_LengtheningStem_LowersDeltaG()
    {
        // A fixed GC-rich arm; prefixes give perfect stems of increasing length around a fixed loop.
        const string fullArm = "GCAGTCAGGTC";
        const string loop = "TTTT";

        double previous = double.PositiveInfinity;
        for (int stem = 3; stem <= 8; stem++)
        {
            string leftArm = fullArm.Substring(0, stem);
            string sequence = leftArm + loop + DnaSequence.GetReverseComplementString(leftArm);

            var hairpin = PrimerDesigner.FindMostStableHairpin(sequence);
            hairpin.Should().NotBeNull(because: $"a perfect {stem}-bp stem closing a 4-nt loop must form a hairpin ('{sequence}')");
            hairpin!.Value.StemLength.Should().BeGreaterThanOrEqualTo(stem,
                because: $"the most stable hairpin uses the full perfect {stem}-bp stem");

            hairpin.Value.DeltaG37.Should().BeLessThan(previous - 1e-9,
                because: $"extending the stem to {stem} bp adds another stabilising NN stack, so ΔG°37 strictly decreases");
            previous = hairpin.Value.DeltaG37;
        }
    }

    [Test]
    [Description("INV: a sequence that admits no ≥2-bp Watson-Crick stem closing a ≥3-nt loop has no hairpin — FindMostStableHairpin returns null — regardless of how long the non-pairing sequence is.")]
    public void Hairpin_NoStemPossible_YieldsNoHairpin()
    {
        // Homopolymers can form no Watson-Crick pair at all → no hairpin, at any length.
        foreach (char b in "ACGT")
        {
            foreach (int len in new[] { 6, 12, 30 })
            {
                string homo = new string(b, len);
                PrimerDesigner.FindMostStableHairpin(homo).Should().BeNull(
                    because: $"a poly-{b} homopolymer has no complementary bases, so no stem and no hairpin ('{homo}')");
            }
        }

        // A single Watson-Crick pair (A…T) cannot make a ≥2-bp stem → still no hairpin.
        PrimerDesigner.FindMostStableHairpin("AAAATAAAA").Should().BeNull(
            because: "one A·T pair closing an all-A loop is only a 1-bp stem (< 2 bp), so no hairpin forms");

        // Non-vacuity contrast: a genuine stem-loop DOES yield a hairpin.
        PrimerDesigner.FindMostStableHairpin("GCGCTTTTGCGC").Should().NotBeNull(
            because: "a 4-bp GC stem closing a 4-nt loop is a real hairpin — the relation is non-vacuous");
    }

    #endregion

    #endregion
}
