using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;
using Site = Seqeron.Genomics.Oncology.OncologyAnalyzer.BreakpointSite;
using Frame = Seqeron.Genomics.Oncology.OncologyAnalyzer.BreakpointFrameStatus;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology fusion-BREAKPOINT-ANALYSIS area — ONCO-FUSION-003.
/// The units under test are <see cref="OncologyAnalyzer.AnalyzeBreakpoint"/> (per-partner
/// breakpoint SITE classification + junction reading-frame consequence) and its protein
/// companion <see cref="OncologyAnalyzer.PredictFusionProtein"/> (chimeric-CDS assembly by
/// breakpoint OFFSET), both implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// This file is scoped to FUSION-003 BREAKPOINT analysis ONLY. Fusion DETECTION
/// (<c>DetectFusions</c>, ONCO-FUSION-001) and known-fusion DB lookup
/// (<c>MatchKnownFusions</c>, ONCO-FUSION-002) are separate checklist rows, covered by
/// OncologyFusionDetectionFuzzTests / OncologyFusionDatabaseFuzzTests, and are NOT
/// re-exercised here.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime exception
/// (IndexOutOfRange / NullReference / Overflow). Every input must resolve to EITHER
/// a well-defined, theory-correct value OR a *documented, intentional* outcome
/// (here, <see cref="ArgumentOutOfRangeException"/> for a breakpoint OFFSET outside its
/// CDS or a coding-junction phase outside {0,1,2}, and <see cref="ArgumentNullException"/>
/// for a null CDS string).
/// For a breakpoint analyzer the headline hazards are:
///   • an IndexOutOfRangeException leaking from the chimeric-CDS slice when a breakpoint
///     OFFSET is past the gene/CDS end or negative ("out-of-bounds") — the documented
///     contract is a guarded ArgumentOutOfRangeException, never an unhandled slice crash;
///   • an OFF-BY-ONE at the gene/exon BOUNDARY: an offset EXACTLY at 0 (gene start) or
///     EXACTLY at CDS.Length (gene end) is the INCLUSIVE valid boundary (empty prefix /
///     empty suffix), while Length+1 is the first invalid position — both must be honoured
///     to the base;
///   • a reading-frame call FABRICATED for an "intronic" (or any non-CDS) breakpoint — a
///     frame is defined ONLY for a coding-to-coding junction; an intronic/UTR/intergenic
///     break must classify as NotPredicted with BreakpointInCoding=false (Arriba
///     reading_frame = '.'), never InFrame/OutOfFrame and never crash.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-FUSION-003 — Fusion breakpoint analysis (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 102.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення (0, -1, MaxInt, empty).
///     Targets (checklist row 102): "breakpoint at gene boundary, intronic, out-of-bounds".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Fusion_Breakpoint_Analysis.md (docs/algorithms/Oncology/Fusion_Breakpoint_Analysis.md):
///   • Site categories = {5'UTR, 3'UTR, UTR, CDS, exon, intron, intergenic}    (§2.2)
///   • A reading-frame call (InFrame/OutOfFrame) is made ONLY when BOTH breakpoints are
///     CDS; else NotPredicted (Arriba reading_frame = '.')                     (INV-01, §6.1)
///   • InFrame ⟺ (fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0      (INV-02)
///   • BreakpointInCoding == (Site5==CDS && Site3==CDS); partners carried through (§3.2)
///   • AnalyzeBreakpoint validates the frame quantities ONLY when both sites are CDS,
///     delegating to IsInFrame (ArgumentOutOfRangeException: phase ∉ {0,1,2})  (§3.3, §6.1)
///   • ChimericCds = 5' CDS prefix [0:junction5] ++ 3' CDS suffix [junction3:] (INV-04, §4.1)
///   • Breakpoint OFFSETS: junction5 = FivePrimeCodingBases, junction3 = ThreePrimeStartPhase;
///     each must be within [0, CDS.Length] — out of range ⇒ ArgumentOutOfRangeException (§3.3)
///   • Out-of-frame chimeric CDS trimmed to whole codons before translation    (INV-05, §4.1)
///   • Peptide truncated at the first stop codon; HasPrematureStop flags it     (INV-03, §6.1)
///   • Null CDS ⇒ ArgumentNullException                                        (§3.3)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyFusionBreakpointFuzzTests
{
    private const int CodonLength = 3;

    // The non-CDS site categories — every one of these on EITHER partner must force
    // a NotPredicted frame call (INV-01). "Intron" is the headline BE target.
    private static readonly Site[] NonCdsSites =
    {
        Site.FivePrimeUtr,
        Site.ThreePrimeUtr,
        Site.Utr,
        Site.Exon,
        Site.Intron,
        Site.Intergenic,
    };

    // ── FusionBreakpoint builder (mirrors the source record signature) ────────
    private static FusionBreakpoint Bp(
        Site site5, Site site3, int fivePrimeCodingBases, int threePrimeStartPhase,
        string g5 = "EML4", string g3 = "ALK") =>
        new(g5, g3, site5, site3, fivePrimeCodingBases, threePrimeStartPhase);

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented structural contract on EVERY BreakpointAnalysis: the partner
    // symbols and site categories are carried through verbatim, BreakpointInCoding is
    // EXACTLY (both sites CDS), a frame is called IFF both sites are CDS (INV-01), and the
    // frame label is never a value outside the documented enum. This is what stops a fuzz
    // test from rubber-stamping a fabricated frame call (e.g. InFrame on an intron) green.
    private static void AssertWellFormedAnalysis(BreakpointAnalysis a, FusionBreakpoint src)
    {
        a.Gene5Prime.Should().Be(src.Gene5Prime, "partners are carried through unchanged (§3.2)");
        a.Gene3Prime.Should().Be(src.Gene3Prime, "partners are carried through unchanged (§3.2)");
        a.Site5Prime.Should().Be(src.Site5Prime, "the 5' site category is carried through (§3.2)");
        a.Site3Prime.Should().Be(src.Site3Prime, "the 3' site category is carried through (§3.2)");

        bool bothCoding = src.Site5Prime == Site.Cds && src.Site3Prime == Site.Cds;
        a.BreakpointInCoding.Should().Be(bothCoding,
            "BreakpointInCoding ⟺ both breakpoints are CDS (§3.2)");

        if (bothCoding)
        {
            a.FrameStatus.Should().BeOneOf(new[] { Frame.InFrame, Frame.OutOfFrame },
                "a coding-to-coding junction is called in- or out-of-frame (INV-01/02)");
        }
        else
        {
            a.FrameStatus.Should().Be(Frame.NotPredicted,
                "a non-coding junction has no frame call (INV-01, Arriba reading_frame = '.')");
        }
    }

    #region ONCO-FUSION-003 — positive sanity (CDS boundary IS classified in-frame; an intron IS NotPredicted)

    // The worked example from §7.1: a CDS::CDS junction with (FivePrimeCodingBases=6,
    // ThreePrimeStartPhase=0) — (6−0) mod 3 == 0 — is in-frame and coding; the chimeric
    // CDS "ATGAAA"++"GATGGT" translates to "MKDG". A fuzz suite that never asserts a TRUE
    // positive proves nothing, so this pins the documented behaviour end-to-end.
    [Test]
    public void AnalyzeBreakpoint_DocumentedWorkedExample_CdsToCds_IsInFrameAndCoding()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 6, threePrimeStartPhase: 0);

        var analysis = AnalyzeBreakpoint(bp);

        analysis.BreakpointInCoding.Should().BeTrue("CDS::CDS is a coding-to-coding junction");
        analysis.FrameStatus.Should().Be(Frame.InFrame, "(6−0) mod 3 == 0 (INV-02)");
        analysis.Gene5Prime.Should().Be("EML4");
        analysis.Gene3Prime.Should().Be("ALK");
        AssertWellFormedAnalysis(analysis, bp);

        var protein = PredictFusionProtein(bp, ("ATGAAA", "GATGGT"));
        protein.ChimericCds.Should().Be("ATGAAAGATGGT", "5' prefix [0:6] ++ 3' suffix [0:] (INV-04)");
        protein.Peptide.Should().Be("MKDG", "ATG(M) AAA(K) GAT(D) GGT(G)");
        protein.Effect.Should().Be(Frame.InFrame);
        protein.HasPrematureStop.Should().BeFalse();
    }

    // An intronic breakpoint (the BE "intronic" target) on the 3' partner must be
    // NotPredicted — a frame is undefined off the coding sequence (INV-01). This is the
    // direct counterpart of the in-frame positive: the SAME frame quantities that would
    // be in-frame in CDS produce NO frame call once a partner is intronic.
    [Test]
    public void AnalyzeBreakpoint_IntronicThreePrimeBreakpoint_IsNotPredicted_NotCoding()
    {
        var bp = Bp(Site.Cds, Site.Intron, fivePrimeCodingBases: 6, threePrimeStartPhase: 0);

        var analysis = AnalyzeBreakpoint(bp);

        analysis.BreakpointInCoding.Should().BeFalse("an intronic breakpoint is not coding (§2.2)");
        analysis.FrameStatus.Should().Be(Frame.NotPredicted,
            "frame is undefined off CDS (INV-01, Arriba reading_frame = '.')");
        AssertWellFormedAnalysis(analysis, bp);
    }

    #endregion

    #region ONCO-FUSION-003 / BE — breakpoint at gene boundary (CDS::CDS, no off-by-one in the codon-phase rule)

    // BE "breakpoint at gene boundary": the codon-phase boundary itself. With phase 0 a
    // 5' coding-base count that is EXACTLY a multiple of 3 sits on a codon boundary and is
    // in-frame; one base either side is out-of-frame. This pins the boundary to the base —
    // no off-by-one in (b − p) mod 3.
    [TestCase(0, 0, Frame.InFrame,    TestName = "b=0  phase=0 → on codon boundary → InFrame")]
    [TestCase(3, 0, Frame.InFrame,    TestName = "b=3  phase=0 → on codon boundary → InFrame")]
    [TestCase(6, 0, Frame.InFrame,    TestName = "b=6  phase=0 → on codon boundary → InFrame")]
    [TestCase(2, 0, Frame.OutOfFrame, TestName = "b=2  phase=0 → 1 base short → OutOfFrame")]
    [TestCase(4, 0, Frame.OutOfFrame, TestName = "b=4  phase=0 → 1 base over → OutOfFrame")]
    public void AnalyzeBreakpoint_CodonBoundaryPhase0_InFrameOnlyOnMultipleOfThree(
        int fivePrimeCodingBases, int phase, Frame expected)
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases, phase);

        var analysis = AnalyzeBreakpoint(bp);

        analysis.FrameStatus.Should().Be(expected,
            "InFrame ⟺ (b − phase) mod 3 == 0, to the base (INV-02)");
        AssertWellFormedAnalysis(analysis, bp);
    }

    // The phase boundary: phase ∈ {0,1,2} is the inclusive valid set; the in-frame base
    // count shifts with the phase. b = phase (any of 0,1,2) is always in-frame
    // ((phase−phase) mod 3 == 0); b = phase+1 is always out-of-frame.
    [TestCase(0, 0, Frame.InFrame,    TestName = "phase=0 b=0 → InFrame")]
    [TestCase(0, 1, Frame.OutOfFrame, TestName = "phase=0 b=1 → OutOfFrame")]
    [TestCase(1, 1, Frame.InFrame,    TestName = "phase=1 b=1 → InFrame")]
    [TestCase(1, 2, Frame.OutOfFrame, TestName = "phase=1 b=2 → OutOfFrame")]
    [TestCase(2, 2, Frame.InFrame,    TestName = "phase=2 b=2 → InFrame")]
    [TestCase(2, 5, Frame.InFrame,    TestName = "phase=2 b=5 → InFrame (5−2=3)")]
    public void AnalyzeBreakpoint_PhaseBoundary_FrameTracksCodonPhase(
        int phase, int fivePrimeCodingBases, Frame expected)
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases, phase);

        AnalyzeBreakpoint(bp).FrameStatus.Should().Be(expected,
            "the in-frame base count is shifted by the 3' coding-start phase (INV-02)");
    }

    // An extreme but in-bounds 5' coding-base count (int.MaxValue is a multiple of 3:
    // 2147483647 mod 3 == 1, so use a constructed large multiple) must not overflow the
    // mod-3 arithmetic. b is a non-negative dividend, so (b − phase) % 3 is well-defined.
    [Test]
    public void AnalyzeBreakpoint_HugeInBoundsCodingBases_NoOverflow_FrameByModThree()
    {
        const int hugeMultipleOfThree = 2147483646; // int.MaxValue − 1, divisible by 3
        var bp = Bp(Site.Cds, Site.Cds, hugeMultipleOfThree, 0);

        AnalyzeBreakpoint(bp).FrameStatus.Should().Be(Frame.InFrame,
            "a huge multiple-of-three coding-base count is in-frame; no overflow in mod 3");
    }

    #endregion

    #region ONCO-FUSION-003 / BE — intronic / non-CDS breakpoint (no fabricated frame call)

    // BE "intronic": an intron (or any non-CDS site) on EITHER partner forces
    // NotPredicted, whatever the frame quantities would otherwise compute. Spans the full
    // non-CDS vocabulary on the 5' side, with 3' fixed CDS.
    [Test]
    public void AnalyzeBreakpoint_NonCdsFivePrimeSite_IsAlwaysNotPredicted([ValueSource(nameof(NonCdsSites))] Site site5)
    {
        // Frame quantities that WOULD be in-frame in CDS::CDS — proving the site, not the
        // arithmetic, drives the NotPredicted call.
        var bp = Bp(site5, Site.Cds, fivePrimeCodingBases: 9, threePrimeStartPhase: 0);

        var analysis = AnalyzeBreakpoint(bp);

        analysis.FrameStatus.Should().Be(Frame.NotPredicted,
            "a non-CDS 5' breakpoint ({0}) has no frame call (INV-01)", site5);
        analysis.BreakpointInCoding.Should().BeFalse();
        AssertWellFormedAnalysis(analysis, bp);
    }

    [Test]
    public void AnalyzeBreakpoint_NonCdsThreePrimeSite_IsAlwaysNotPredicted([ValueSource(nameof(NonCdsSites))] Site site3)
    {
        var bp = Bp(Site.Cds, site3, fivePrimeCodingBases: 9, threePrimeStartPhase: 0);

        var analysis = AnalyzeBreakpoint(bp);

        analysis.FrameStatus.Should().Be(Frame.NotPredicted,
            "a non-CDS 3' breakpoint ({0}) has no frame call (INV-01)", site3);
        analysis.BreakpointInCoding.Should().BeFalse();
        AssertWellFormedAnalysis(analysis, bp);
    }

    // A non-CDS partner short-circuits the frame quantities: AnalyzeBreakpoint validates
    // them ONLY for a coding-to-coding junction (§3.3). So even an INVALID phase (3) or a
    // NEGATIVE coding-base count must NOT throw when a site is intronic — it is simply
    // NotPredicted, never an ArgumentOutOfRangeException leak.
    [TestCase(-1, 0, TestName = "intron + negative coding bases → NotPredicted, no throw")]
    [TestCase(9, 3, TestName = "intron + invalid phase 3 → NotPredicted, no throw")]
    [TestCase(9, -1, TestName = "intron + negative phase → NotPredicted, no throw")]
    [TestCase(int.MinValue, int.MaxValue, TestName = "intron + extreme quantities → NotPredicted, no throw")]
    public void AnalyzeBreakpoint_IntronicSite_DoesNotValidateFrameQuantities(int b, int phase)
    {
        var bp = Bp(Site.Cds, Site.Intron, b, phase);

        Action act = () =>
        {
            var analysis = AnalyzeBreakpoint(bp);
            analysis.FrameStatus.Should().Be(Frame.NotPredicted,
                "frame quantities are not validated for a non-coding junction (§3.3)");
        };

        act.Should().NotThrow("a non-CDS junction short-circuits the frame check (§3.3)");
    }

    // Fuzz: random site pairs over the FULL vocabulary with random frame quantities. A
    // frame call must appear IFF both sites are CDS; otherwise NotPredicted. The analyzer
    // must never throw on a non-CDS junction (it does not reach the frame validation), and
    // the result must always be well-formed.
    [Test]
    [CancelAfter(10000)]
    public void AnalyzeBreakpoint_FuzzSitePairs_FrameCalledIffBothCds([Random(1, 100000, 50)] int seed)
    {
        var rng = new Random(seed);
        Site[] allSites = (Site[])Enum.GetValues(typeof(Site));
        Site site5 = allSites[rng.Next(allSites.Length)];
        Site site3 = allSites[rng.Next(allSites.Length)];
        // Keep frame quantities VALID so a CDS::CDS draw exercises the real frame call.
        int b = rng.Next(0, 100);
        int phase = rng.Next(0, CodonLength);

        var bp = Bp(site5, site3, b, phase);
        var analysis = AnalyzeBreakpoint(bp);

        bool bothCds = site5 == Site.Cds && site3 == Site.Cds;
        if (bothCds)
        {
            Frame expected = (b - phase) % CodonLength == 0 ? Frame.InFrame : Frame.OutOfFrame;
            analysis.FrameStatus.Should().Be(expected,
                "CDS::CDS is called by the codon-phase rule (INV-02, seed {0})", seed);
        }
        else
        {
            analysis.FrameStatus.Should().Be(Frame.NotPredicted,
                "a non-CDS junction is NotPredicted (INV-01, seed {0})", seed);
        }

        AssertWellFormedAnalysis(analysis, bp);
    }

    #endregion

    #region ONCO-FUSION-003 / BE — invalid coding-junction phase (documented ArgumentOutOfRangeException)

    // For a CDS::CDS junction the frame quantities ARE validated (delegated to IsInFrame).
    // A phase outside {0,1,2} is the documented ArgumentOutOfRangeException (§3.3, §6.1) —
    // the boundary just past the valid set (3) and the negative side (−1).
    [TestCase(9, 3, TestName = "phase=3 (just past {0,1,2}) → throws")]
    [TestCase(9, -1, TestName = "phase=−1 → throws")]
    [TestCase(9, int.MaxValue, TestName = "phase=MaxInt → throws")]
    public void AnalyzeBreakpoint_CdsJunctionWithInvalidPhase_ThrowsArgumentOutOfRange(int b, int phase)
    {
        var bp = Bp(Site.Cds, Site.Cds, b, phase);

        Action act = () => AnalyzeBreakpoint(bp);

        act.Should().Throw<ArgumentOutOfRangeException>(
            "a coding junction validates the phase ∈ {0,1,2} via IsInFrame (§3.3)");
    }

    // A negative coding-base count on a CDS::CDS junction is also the documented throw.
    [TestCase(-1, 0)]
    [TestCase(int.MinValue, 0)]
    public void AnalyzeBreakpoint_CdsJunctionWithNegativeCodingBases_ThrowsArgumentOutOfRange(int b, int phase)
    {
        var bp = Bp(Site.Cds, Site.Cds, b, phase);

        Action act = () => AnalyzeBreakpoint(bp);

        act.Should().Throw<ArgumentOutOfRangeException>(
            "a coding junction rejects a negative coding-base count via IsInFrame (§3.3)");
    }

    #endregion

    #region ONCO-FUSION-003 / BE — breakpoint OFFSET at the CDS boundary (PredictFusionProtein, inclusive [0, Length])

    // BE "breakpoint at gene boundary" for the protein offset: the 5' prefix length
    // junction5 may be EXACTLY 0 (empty prefix, gene start) or EXACTLY the CDS length
    // (whole 5' CDS, gene end) — both INCLUSIVE valid boundaries. No off-by-one.
    [Test]
    public void PredictFusionProtein_FivePrimePrefixAtZeroBoundary_EmptyPrefix_NoCrash()
    {
        // junction5 = 0 (gene start); phase 0 ⇒ in-frame; chimeric = "" ++ full 3' CDS.
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 0, threePrimeStartPhase: 0);

        var p = PredictFusionProtein(bp, ("ATGAAA", "GATGGT"));

        p.ChimericCds.Should().Be("GATGGT", "an empty 5' prefix [0:0] leaves only the 3' suffix (INV-04)");
        p.Effect.Should().Be(Frame.InFrame, "(0−0) mod 3 == 0");
    }

    [Test]
    public void PredictFusionProtein_FivePrimePrefixAtLengthBoundary_WholeCds_NoOffByOne()
    {
        // junction5 == CDS.Length (6) is the INCLUSIVE upper boundary: the WHOLE 5' CDS is
        // taken as the prefix. Length+1 would be the first out-of-bounds offset (tested
        // separately). phase 0 keeps it in-frame (6 mod 3 == 0).
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 6, threePrimeStartPhase: 0);

        var p = PredictFusionProtein(bp, ("ATGAAA", "GGGTTT"));

        p.ChimericCds.Should().Be("ATGAAAGGGTTT",
            "junction5 == CDS.Length takes the whole 5' CDS (inclusive boundary, no off-by-one)");
    }

    // The 3' suffix start junction3 may be EXACTLY the 3' CDS length — an empty suffix
    // (the whole 3' CDS is upstream of the breakpoint). The peptide then comes from the
    // 5' prefix alone (§6.1 "Empty 3' CDS suffix"). Here phase must be in {0,1,2}, so the
    // 3' CDS is built to length 2 so that junction3 == 2 is both a valid phase AND the
    // exact length boundary.
    [Test]
    public void PredictFusionProtein_ThreePrimeSuffixAtLengthBoundary_EmptySuffix_PeptideFromFivePrime()
    {
        // 3' CDS length 2, junction3 = 2 == length ⇒ empty suffix; 5' prefix length 3.
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 3, threePrimeStartPhase: 2);

        var p = PredictFusionProtein(bp, ("ATGAAA", "GG"));

        p.ChimericCds.Should().Be("ATG", "5' prefix [0:3] ++ 3' suffix [2:] (empty) (§6.1)");
        p.Peptide.Should().Be("M", "ATG(M); the peptide comes from the 5' prefix alone");
    }

    #endregion

    #region ONCO-FUSION-003 / BE — breakpoint OFFSET out-of-bounds (documented guard, NO IndexOutOfRange leak)

    // BE "out-of-bounds": a 5' prefix length PAST the CDS end (Length+1, the first invalid
    // offset) must be the documented ArgumentOutOfRangeException — NEVER an unhandled
    // IndexOutOfRange/ArgumentException from the Substring slice. This is the headline
    // boundary-exploitation hazard for an offset-driven slice.
    [Test]
    public void PredictFusionProtein_FivePrimeOffsetPastCdsEnd_ThrowsArgumentOutOfRange_NoIndexLeak()
    {
        // CDS length 6, junction5 = 7 = Length+1 (first invalid). phase 1 in {0,1,2}.
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 7, threePrimeStartPhase: 1);

        Action act = () => PredictFusionProtein(bp, ("ATGAAA", "GATGGT"));

        act.Should().Throw<ArgumentOutOfRangeException>(
            "a 5' offset past CDS.Length is guarded, not an IndexOutOfRange slice crash (§3.3)");
    }

    [Test]
    public void PredictFusionProtein_FivePrimeOffsetNegative_ThrowsArgumentOutOfRange()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: -1, threePrimeStartPhase: 0);

        Action act = () => PredictFusionProtein(bp, ("ATGAAA", "GATGGT"));

        act.Should().Throw<ArgumentOutOfRangeException>(
            "a negative 5' offset is guarded (§3.3)");
    }

    // The 3' suffix start past the 3' CDS end is likewise guarded. (junction3 is also the
    // phase, but the OFFSET range check fires when it exceeds the 3' CDS length.)
    [Test]
    public void PredictFusionProtein_ThreePrimeOffsetPastCdsEnd_ThrowsArgumentOutOfRange_NoIndexLeak()
    {
        // 3' CDS length 1, junction3 = 2 = Length+1 (first invalid for THIS short CDS).
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 3, threePrimeStartPhase: 2);

        Action act = () => PredictFusionProtein(bp, ("ATGAAA", "G"));

        act.Should().Throw<ArgumentOutOfRangeException>(
            "a 3' offset past CDS.Length is guarded, not an IndexOutOfRange slice crash (§3.3)");
    }

    // A null CDS string is the documented ArgumentNullException (§3.3), on either side.
    [Test]
    public void PredictFusionProtein_NullFivePrimeCds_ThrowsArgumentNull()
    {
        var bp = Bp(Site.Cds, Site.Cds, 3, 0);

        Action act = () => PredictFusionProtein(bp, (null!, "GATGGT"));

        act.Should().Throw<ArgumentNullException>("a null 5' CDS is the documented throw (§3.3)");
    }

    [Test]
    public void PredictFusionProtein_NullThreePrimeCds_ThrowsArgumentNull()
    {
        var bp = Bp(Site.Cds, Site.Cds, 3, 0);

        Action act = () => PredictFusionProtein(bp, ("ATGAAA", null!));

        act.Should().Throw<ArgumentNullException>("a null 3' CDS is the documented throw (§3.3)");
    }

    // Fuzz: random offsets against fixed-length CDS strings. Whatever the offsets, the
    // result is EITHER a well-formed prediction whose ChimericCds is exactly
    // prefix[0:j5]++suffix[j3:] (when BOTH offsets are in [0,Length]) OR a guarded
    // ArgumentOutOfRangeException — never an IndexOutOfRange/ArgumentException slice leak.
    [Test]
    [CancelAfter(15000)]
    public void PredictFusionProtein_FuzzOffsets_InRangeBuildsExactChimera_OutOfRangeGuarded(
        [Random(1, 100000, 60)] int seed)
    {
        var rng = new Random(seed);
        const string cds5 = "ATGAAAGGGTTT"; // length 12
        const string cds3 = "GATGGTCCC";    // length 9

        // Draw offsets that straddle the boundaries on both sides (negative, in-range, past).
        int j5 = rng.Next(-3, cds5.Length + 4);
        int j3 = rng.Next(-3, cds3.Length + 4);
        var bp = Bp(Site.Cds, Site.Cds, j5, j3);

        bool j5InRange = j5 >= 0 && j5 <= cds5.Length;
        bool j3InRange = j3 >= 0 && j3 <= cds3.Length;
        // For a CDS::CDS junction the source ALSO validates the phase ∈ {0,1,2} via the
        // frame effect path AFTER the offset checks, but the offset checks fire first; a
        // valid prediction additionally needs j3 to be a valid phase only for the Effect
        // label — the offset range alone governs whether a slice throws. j3 is bounded by
        // cds3.Length (9) so j3 ∈ {0,1,2} is the only in-frame-eligible subset; we assert
        // structural slicing on the offsets and accept either frame label.

        Action act = () =>
        {
            if (j5InRange && j3InRange)
            {
                var p = PredictFusionProtein(bp, (cds5, cds3));
                string expected = cds5.Substring(0, j5) + cds3.Substring(j3);
                p.ChimericCds.Should().Be(expected,
                    "ChimericCds = 5' prefix [0:j5] ++ 3' suffix [j3:] (INV-04, seed {0})", seed);
                p.Effect.Should().BeOneOf(new[] { Frame.InFrame, Frame.OutOfFrame },
                    "the frame effect is always in- or out-of-frame (seed {0})", seed);
                int translatable = p.ChimericCds.Length - (p.ChimericCds.Length % CodonLength);
                p.Peptide.Length.Should().BeLessThanOrEqualTo(translatable / CodonLength,
                    "the peptide cannot be longer than the translatable codon count (seed {0})", seed);
            }
            else
            {
                Action call = () => PredictFusionProtein(bp, (cds5, cds3));
                call.Should().Throw<ArgumentOutOfRangeException>(
                    "an out-of-range offset is guarded, never an IndexOutOfRange leak (seed {0})", seed);
            }
        };

        act.Should().NotThrow("the only thrown exception is the documented ArgumentOutOfRangeException (seed {0})", seed);
    }

    #endregion

    #region ONCO-FUSION-003 / BE — out-of-frame trimming and premature-stop boundary

    // INV-05: an out-of-frame chimeric CDS is trimmed to whole codons before translation.
    // With junction5=5 (phase 0 ⇒ out-of-frame, since 5 mod 3 == 2) the chimera is read in
    // its shifted frame; a trailing partial codon is dropped, never an IndexOutOfRange.
    [Test]
    public void PredictFusionProtein_OutOfFrameJunction_TrimsToWholeCodons_NoCrash()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 5, threePrimeStartPhase: 0);

        var p = PredictFusionProtein(bp, ("ATGAA", "GATGGT")); // 5' CDS length 5

        p.Effect.Should().Be(Frame.OutOfFrame, "(5−0) mod 3 != 0 (INV-02)");
        // chimeric "ATGAA" ++ "GATGGT" = "ATGAAGATGGT" (length 11) → trimmed to 9 → 3 codons.
        p.ChimericCds.Should().Be("ATGAAGATGGT");
        p.Peptide.Length.Should().BeLessThanOrEqualTo(p.ChimericCds.Length / CodonLength,
            "translation reads whole codons only (INV-05)");
    }

    // INV-03: a stop codon in the chimeric ORF truncates the peptide at the first stop and
    // flags HasPrematureStop — the §7.1 walk-through (chimera reaches TAA → "MKD").
    [Test]
    public void PredictFusionProtein_PrematureStop_TruncatesAtFirstStop_FlagsIt()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 6, threePrimeStartPhase: 0);

        var p = PredictFusionProtein(bp, ("ATGAAA", "GATTAAGGT"));

        p.Peptide.Should().Be("MKD", "ATG(M) AAA(K) GAT(D) TAA(stop) (INV-03, §7.1)");
        p.HasPrematureStop.Should().BeTrue("a stop codon was reached before the ORF end (INV-03)");
    }

    #endregion
}
