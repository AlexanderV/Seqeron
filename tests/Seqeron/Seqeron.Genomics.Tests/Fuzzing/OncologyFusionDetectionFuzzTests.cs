using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology gene-fusion-DETECTION area — ONCO-FUSION-001.
/// The unit under test is the deterministic, threshold-driven STAR-Fusion-rule
/// fusion caller <see cref="OncologyAnalyzer.DetectFusions"/> (and its building
/// block <see cref="OncologyAnalyzer.ComputeTotalSupport"/>), implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// This file is scoped to FUSION-001 DETECTION ONLY. Partner-DB lookup
/// (<c>MatchKnownFusions</c>, ONCO-FUSION-002) and breakpoint/protein analysis
/// (<c>AnalyzeBreakpoint</c>/<c>PredictFusionProtein</c>, ONCO-FUSION-003) are
/// separate checklist rows and are NOT exercised here.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime exception
/// (DivideByZero / NullReference / Overflow). Every input must resolve to EITHER
/// a well-defined, theory-correct value OR a *documented, intentional* outcome
/// (here, <see cref="ArgumentNullException"/> for a null candidate enumerable,
/// and <see cref="ArgumentException"/> for a negative supporting-read count).
/// For threshold-driven fusion detection the headline hazards are:
///   • a NullReferenceException while reading a malformed candidate record (a
///     null gene symbol) — the documented contract treats a null/empty 5' or 3'
///     symbol as data, never a crash;
///   • a self-fusion or identical-gene candidate (gene5p == gene3p) being
///     FALSELY reported — the documented rule (INV-01) requires two DISTINCT
///     partner genes, so such candidates are silently SKIPPED, never called and
///     never crashed on;
///   • an empty / no-chimeric-read input producing anything other than an empty
///     list (no DivideByZeroException in any support fraction, no NullReference);
///   • a false positive: a candidate BELOW the STAR-Fusion support threshold must
///     yield NO call (INV-03), and a candidate AT/above it must yield exactly one
///     call with the documented partners and TotalSupport (INV-02/04).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-FUSION-001 — Gene fusion detection (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 100.
/// Fuzz strategy exercised for THIS unit:
///   • MC = Malformed Content — невалідний контент. Targets (checklist row 100):
///     "no chimeric reads, self-fusion, identical genes, empty reads".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Fusion_Gene_Detection.md (docs/algorithms/Oncology/Fusion_Gene_Detection.md):
///   • junction_reads = split_reads1 + split_reads2                          (§2.2)
///   • TotalSupport   = split_reads1 + split_reads2 + discordant_mates       (§2.2, INV-02)
///   • Detection rule (STAR-Fusion defaults, §2.2 / §4.2, INV-03):
///       - if junction_reads ≥ 1: report iff junction_reads ≥ MIN_JUNCTION_READS(=1)
///           AND total_support ≥ MIN_SUM_FRAGS(=2);
///       - if junction_reads == 0: report iff discordant_mates ≥ MIN_SPANNING_FRAGS_ONLY(=5).
///   • gene5p == gene3p (case-insensitive) ⇒ candidate SKIPPED, NOT a fusion    (§3.3, INV-01, §6.1)
///   • Empty input ⇒ empty list                                               (§3.3, §6.1)
///   • Negative supporting-read count ⇒ ArgumentException                      (§3.3)
///   • Null candidates ⇒ ArgumentNullException                                (§3.3)
///   • Results ordered by DESCENDING TotalSupport, gene pair as tie-break      (§2.4 INV-04, §4.1)
///   • Reading frame: InFrame ⇔ (fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0;
///       phase fields unset (-1) ⇒ Unknown (never guessed)                     (§2.2 INV-05, §3.3, §6.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyFusionDetectionFuzzTests
{
    private const int MinJunctionReads = DefaultMinJunctionReads;     // 1
    private const int MinSumFrags = DefaultMinSumFrags;               // 2
    private const int MinSpanningFragsOnly = DefaultMinSpanningFragsOnly; // 5

    // ── Candidate builders (mirror the source FusionCandidate signature) ─────
    private static FusionCandidate Cand(
        string g5, string g3, int split5, int split3, int discordant,
        int fivePrimeCodingBases = -1, int threePrimeStartPhase = -1) =>
        new(g5, g3, split5, split3, discordant, fivePrimeCodingBases, threePrimeStartPhase);

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented numeric contract on EVERY emitted call: TotalSupport is
    // the exact Arriba sum (INV-02), JunctionReads/DiscordantMates are
    // non-negative, the partners are DISTINCT (INV-01), and the STAR-Fusion
    // pass-rule (INV-03) genuinely held. This is what stops a fuzz test from
    // rubber-stamping a mislabelled / impossible call green.
    private static void AssertWellFormedCall(FusionCall c)
    {
        c.JunctionReads.Should().BeGreaterThanOrEqualTo(0, "junction reads is a count");
        c.DiscordantMates.Should().BeGreaterThanOrEqualTo(0, "discordant mates is a count");
        c.TotalSupport.Should().Be(c.JunctionReads + c.DiscordantMates,
            "TotalSupport = split1+split2+discordant (INV-02, §2.2)");

        // INV-01: a reported fusion has two DISTINCT partners (case-insensitive).
        string.Equals(c.Gene5Prime, c.Gene3Prime, StringComparison.OrdinalIgnoreCase)
            .Should().BeFalse("a reported fusion needs two distinct partner genes (INV-01)");

        // INV-03: the STAR-Fusion threshold rule genuinely passed.
        bool passes = c.JunctionReads >= MinJunctionReads
            ? c.TotalSupport >= MinSumFrags
            : c.DiscordantMates >= MinSpanningFragsOnly;
        passes.Should().BeTrue("a reported fusion must satisfy the STAR-Fusion rule (INV-03)");
    }

    // Confirms the result list is ordered by descending TotalSupport (INV-04).
    private static void AssertDescendingSupport(IReadOnlyList<FusionCall> calls)
    {
        for (int i = 1; i < calls.Count; i++)
        {
            calls[i - 1].TotalSupport.Should().BeGreaterThanOrEqualTo(calls[i].TotalSupport,
                "results are ordered by descending TotalSupport (INV-04, §2.4)");
        }
    }

    #region ONCO-FUSION-001 — positive sanity (a real fusion IS called; a sub-threshold one is NOT)

    // The worked example from §7.1: EML4::ALK (junction=5, total=9, in-frame) and
    // CD74::ROS1 (spanning-only, total=5) are called; NCOA4::RET (spanning-only,
    // 4 < 5) is rejected. A fuzz suite that never asserts a TRUE positive proves
    // nothing, so this pins the documented behaviour end-to-end.
    [Test]
    public void DetectFusions_DocumentedWorkedExample_CallsExactlyTheTwoSupportedFusions()
    {
        var candidates = new[]
        {
            Cand("EML4", "ALK", 3, 2, 4, 300, 0),  // junction=5, total=9, in-frame
            Cand("CD74", "ROS1", 0, 0, 5),          // spanning-only, total=5 → detected
            Cand("NCOA4", "RET", 0, 0, 4),          // spanning-only, 4 < 5 → rejected
        };

        var calls = DetectFusions(candidates);

        calls.Should().HaveCount(2, "two candidates pass the STAR-Fusion rule, one is sub-threshold");
        calls.Select(c => (c.Gene5Prime, c.Gene3Prime))
            .Should().NotContain(("NCOA4", "RET"), "NCOA4::RET has 4 < MIN_SPANNING_FRAGS_ONLY (5)");

        var eml4Alk = calls.Single(c => c.Gene5Prime == "EML4");
        eml4Alk.Gene3Prime.Should().Be("ALK");
        eml4Alk.JunctionReads.Should().Be(5, "split1+split2 = 3+2");
        eml4Alk.DiscordantMates.Should().Be(4);
        eml4Alk.TotalSupport.Should().Be(9, "3+2+4 (INV-02)");
        eml4Alk.ReadingFrame.Should().Be(FusionReadingFrame.InFrame, "(300−0) mod 3 == 0 (INV-05)");

        var cd74Ros1 = calls.Single(c => c.Gene5Prime == "CD74");
        cd74Ros1.Gene3Prime.Should().Be("ROS1");
        cd74Ros1.TotalSupport.Should().Be(5);

        AssertDescendingSupport(calls);
        foreach (var c in calls)
        {
            AssertWellFormedCall(c);
        }
    }

    // The exact STAR-Fusion threshold boundaries, on both branches (§6.1 edge table).
    [TestCase(1, 0, 0, false, TestName = "Junction=1 total=1 < MIN_SUM_FRAGS(2) → rejected")]
    [TestCase(1, 1, 0, true,  TestName = "Junction=1 total=2 == MIN_SUM_FRAGS(2) → detected")]
    [TestCase(0, 0, 4, false, TestName = "Spanning-only 4 < MIN_SPANNING_FRAGS_ONLY(5) → rejected")]
    [TestCase(0, 0, 5, true,  TestName = "Spanning-only 5 == MIN_SPANNING_FRAGS_ONLY(5) → detected")]
    public void DetectFusions_AtAndBelowSupportThreshold_CallsIffRuleSatisfied(
        int split5, int split3, int discordant, bool expectCalled)
    {
        var calls = DetectFusions(new[] { Cand("A", "B", split5, split3, discordant) });

        calls.Should().HaveCount(expectCalled ? 1 : 0,
            "the STAR-Fusion threshold rule decides detection (INV-03)");
        foreach (var c in calls)
        {
            AssertWellFormedCall(c);
        }
    }

    #endregion

    #region ONCO-FUSION-001 / MC — empty reads (no candidates)

    // Target "empty reads": no candidates ⇒ empty list, never a DivideByZero in a
    // support fraction nor a NullReference (§3.3, §6.1).
    [Test]
    public void DetectFusions_EmptyCandidateList_ReturnsEmpty_NoCrash()
    {
        var calls = DetectFusions(Array.Empty<FusionCandidate>());

        calls.Should().NotBeNull();
        calls.Should().BeEmpty("empty input ⇒ empty list (§6.1)");
    }

    [Test]
    public void DetectFusions_NullCandidates_ThrowsArgumentNullException()
    {
        Action act = () => DetectFusions(null!);

        act.Should().Throw<ArgumentNullException>("null candidates is the documented throw (§3.3)");
    }

    // A lazily-evaluated empty enumerable (a different code path from an array)
    // must still terminate with an empty list and no exception.
    [Test]
    [CancelAfter(5000)]
    public void DetectFusions_LazyEmptyEnumerable_ReturnsEmpty_NoHang()
    {
        IEnumerable<FusionCandidate> Empty()
        {
            yield break;
        }

        var calls = DetectFusions(Empty());

        calls.Should().BeEmpty();
    }

    #endregion

    #region ONCO-FUSION-001 / MC — no chimeric reads (candidates present, zero support)

    // Target "no chimeric reads": candidates are present but carry NO supporting
    // reads (all classes zero). No candidate passes the rule ⇒ empty result, no
    // crash. This is the all-zero / no-evidence input.
    [Test]
    public void DetectFusions_CandidatesWithZeroSupport_ReturnsEmpty_NoCrash()
    {
        var candidates = new[]
        {
            Cand("EML4", "ALK", 0, 0, 0),
            Cand("BCR", "ABL1", 0, 0, 0),
            Cand("TMPRSS2", "ERG", 0, 0, 0),
        };

        var calls = DetectFusions(candidates);

        calls.Should().BeEmpty("zero support fails the STAR-Fusion rule for every candidate (INV-03)");
    }

    // Fuzz: many random ALL-ZERO-support distinct-gene candidates. Whatever the
    // partner symbols, with no reads none can ever pass; the result must stay
    // empty and the call must never throw.
    [Test]
    [CancelAfter(10000)]
    public void DetectFusions_FuzzManyZeroSupportCandidates_AlwaysEmpty([Random(1, 100000, 25)] int seed)
    {
        var rng = new Random(seed);
        int n = rng.Next(0, 40);
        var candidates = new List<FusionCandidate>(n);
        for (int i = 0; i < n; i++)
        {
            candidates.Add(Cand($"G{rng.Next(0, 1000)}", $"H{rng.Next(0, 1000)}", 0, 0, 0));
        }

        var calls = DetectFusions(candidates);

        calls.Should().BeEmpty("no supporting reads ⇒ no fusion, for any partner symbols (seed {0})", seed);
    }

    #endregion

    #region ONCO-FUSION-001 / MC — self-fusion / identical genes (gene5p == gene3p)

    // Target "self-fusion" / "identical genes": a candidate whose two partners are
    // the SAME gene must be SKIPPED (INV-01) — a fusion needs two distinct genes.
    // It must never be reported and never crash, even with overwhelming support.
    [Test]
    public void DetectFusions_SelfFusionWithHeavySupport_IsNeverCalled()
    {
        var candidates = new[]
        {
            Cand("ALK", "ALK", 1000, 1000, 1000, 30, 0), // identical genes, massive support
        };

        var calls = DetectFusions(candidates);

        calls.Should().BeEmpty("gene5p == gene3p is not a fusion (INV-01, §6.1)");
    }

    // Identical genes are matched case-INSENSITIVELY (source uses OrdinalIgnoreCase),
    // so "alk"/"ALK" is still a self-fusion and must be skipped.
    [TestCase("ALK", "ALK")]
    [TestCase("alk", "ALK")]
    [TestCase("Alk", "aLK")]
    [TestCase("BCR", "bcr")]
    public void DetectFusions_IdenticalGenesCaseInsensitive_AreSkipped(string g5, string g3)
    {
        var calls = DetectFusions(new[] { Cand(g5, g3, 10, 10, 10) });

        calls.Should().BeEmpty("same gene (case-insensitive) is not a fusion (INV-01)");
    }

    // A self-fusion mixed into a batch must be the ONLY thing dropped — the real,
    // distinct-gene fusion alongside it is still called correctly.
    [Test]
    public void DetectFusions_SelfFusionAmongRealFusions_DropsOnlyTheSelfFusion()
    {
        var candidates = new[]
        {
            Cand("ALK", "ALK", 50, 50, 50),     // self-fusion → skipped
            Cand("EML4", "ALK", 3, 2, 4),       // real fusion → called
        };

        var calls = DetectFusions(candidates);

        calls.Should().ContainSingle("only the distinct-gene candidate is a fusion");
        calls[0].Gene5Prime.Should().Be("EML4");
        calls[0].Gene3Prime.Should().Be("ALK");
        AssertWellFormedCall(calls[0]);
    }

    #endregion

    #region ONCO-FUSION-001 / MC — malformed read records (null / empty / whitespace symbols)

    // A null partner symbol is malformed CONTENT, not a crash trigger. With
    // distinct (null vs non-null) symbols and sufficient support, the candidate is
    // a valid fusion and must be reported; the call must never NullReference.
    [Test]
    public void DetectFusions_NullPartnerSymbolWithSupport_DoesNotCrash_AndIsTreatedAsDistinct()
    {
        Action act = () =>
        {
            var calls = DetectFusions(new[] { Cand(null!, "ALK", 3, 2, 4) });
            // null vs "ALK" are distinct ⇒ a fusion is callable; assert no corruption.
            foreach (var c in calls)
            {
                AssertWellFormedCall(c);
            }
        };

        act.Should().NotThrow("a null partner symbol is malformed data, not a crash (MC)");
    }

    // Two null symbols compare equal (string.Equals(null, null) == true), so this
    // is a self-fusion and must be SKIPPED, not crash.
    [Test]
    public void DetectFusions_BothPartnersNull_AreEqual_SoSkipped_NoCrash()
    {
        Action act = () =>
        {
            var calls = DetectFusions(new[] { Cand(null!, null!, 10, 10, 10) });
            calls.Should().BeEmpty("null == null ⇒ self-fusion ⇒ skipped (INV-01)");
        };

        act.Should().NotThrow("two null symbols are malformed data, not a crash (MC)");
    }

    // Empty / whitespace symbols are distinct from a real gene and must not crash.
    [TestCase("", "ALK")]
    [TestCase("ALK", "")]
    [TestCase("   ", "ALK")]
    [TestCase("\t", "\n")]
    public void DetectFusions_EmptyOrWhitespaceSymbols_DoNotCrash(string g5, string g3)
    {
        Action act = () =>
        {
            var calls = DetectFusions(new[] { Cand(g5, g3, 3, 2, 4) });
            foreach (var c in calls)
            {
                AssertWellFormedCall(c);
            }
        };

        act.Should().NotThrow("empty/whitespace partner symbols are malformed data, not a crash (MC)");
    }

    #endregion

    #region ONCO-FUSION-001 / MC — negative read counts (documented ArgumentException)

    // Negative supporting-read counts are the documented ArgumentException (§3.3).
    [TestCase(-1, 0, 0)]
    [TestCase(0, -1, 0)]
    [TestCase(0, 0, -1)]
    [TestCase(int.MinValue, 0, 0)]
    public void DetectFusions_NegativeReadCount_ThrowsArgumentException(int split5, int split3, int disc)
    {
        Action act = () => DetectFusions(new[] { Cand("A", "B", split5, split3, disc) });

        act.Should().Throw<ArgumentException>("negative supporting-read counts are rejected (§3.3)");
    }

    #endregion

    #region ONCO-FUSION-001 / MC — mixed batch invariants under fuzz

    // Broad fuzz over mixed batches: distinct-gene candidates with random support,
    // self-fusions, and zero-support entries interleaved. Whatever the mix, the
    // result must (a) never throw, (b) contain ONLY distinct-gene calls satisfying
    // the STAR-Fusion rule, (c) be ordered by descending TotalSupport, and
    // (d) never contain any self-fusion partner pair.
    [Test]
    [CancelAfter(15000)]
    public void DetectFusions_FuzzMixedBatch_HoldsAllInvariants([Random(1, 100000, 40)] int seed)
    {
        var rng = new Random(seed);
        int n = rng.Next(0, 30);
        var candidates = new List<FusionCandidate>(n);

        for (int i = 0; i < n; i++)
        {
            int kind = rng.Next(0, 3);
            int s5 = rng.Next(0, 20);
            int s3 = rng.Next(0, 20);
            int dm = rng.Next(0, 20);
            string gene = $"G{rng.Next(0, 8)}";
            candidates.Add(kind switch
            {
                0 => Cand(gene, gene, s5, s3, dm),                       // self-fusion (must drop)
                1 => Cand(gene, $"H{rng.Next(0, 8)}", 0, 0, 0),          // zero support (must drop)
                _ => Cand(gene, $"P{rng.Next(0, 8)}", s5, s3, dm),       // candidate fusion
            });
        }

        IReadOnlyList<FusionCall> calls = DetectFusions(candidates);

        calls.Should().NotBeNull();
        AssertDescendingSupport(calls);
        foreach (var c in calls)
        {
            AssertWellFormedCall(c); // includes distinct-partner + threshold checks
        }
    }

    // ComputeTotalSupport is the Arriba sum and must equal split1+split2+discordant
    // for arbitrary non-negative counts (INV-02), with no overflow surprises.
    [Test]
    public void ComputeTotalSupport_MatchesArribaSum([Random(0, 100000, 30)] int seed)
    {
        var rng = new Random(seed);
        int s5 = rng.Next(0, 10000);
        int s3 = rng.Next(0, 10000);
        int dm = rng.Next(0, 10000);

        ComputeTotalSupport(Cand("A", "B", s5, s3, dm))
            .Should().Be(s5 + s3 + dm, "TotalSupport = split1+split2+discordant (INV-02)");
    }

    #endregion
}
