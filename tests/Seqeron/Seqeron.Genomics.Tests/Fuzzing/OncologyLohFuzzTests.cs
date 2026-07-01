using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology loss-of-heterozygosity area — ONCO-LOH-001.
/// The units under test are the deterministic HRD-LOH entry points
/// <see cref="OncologyAnalyzer.DetectLOH(IEnumerable{OncologyAnalyzer.AlleleSpecificSegment})"/>
/// (qualifying regions + HRD-LOH score),
/// <see cref="OncologyAnalyzer.CalculateHrdLohScore(IEnumerable{OncologyAnalyzer.AlleleSpecificSegment})"/>
/// (the score directly), and
/// <see cref="OncologyAnalyzer.CalculateLOHFraction(IEnumerable{OncologyAnalyzer.AlleleSpecificSegment}, string)"/>
/// (per-chromosome length-weighted LOH burden ∈ [0,1]),
/// implemented in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / extreme inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense
/// output, and no *unhandled* runtime fault. Every input must resolve to EITHER
/// a well-defined, theory-correct value OR a *documented, intentional* outcome
/// (here, an ArgumentNullException for null inputs or an ArgumentException for an
/// invalid segment). The reported HRD-LOH score must always be a non-negative
/// count (INV-01) and the LOH fraction must always be a FINITE value in [0,1]
/// (INV-02) — never NaN / ±Inf from a divide-by-zero over an empty or single-SNP
/// region. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// Note on "BAF" framing. The checklist row phrases the targets in B-allele-
/// frequency terms ("BAF=0.5 everywhere, BAF=0/1, single SNP"); the documented
/// unit operates one step downstream on already-segmented allele-specific copy
/// number (Loss_Of_Heterozygosity.md §5.2: "upstream B-allele-frequency modelling
/// is out of scope"). The BAF targets map exactly onto the segment contract:
///   • BAF ≈ 0.5 everywhere  ⇔ perfectly heterozygous, minor CN ≠ 0
///                              ⇒ NOT LOH, score 0, fraction 0   (§6.1 "Heterozygous retained");
///   • BAF = 0 or 1 everywhere ⇔ complete allelic imbalance, minor CN == 0
///                              ⇒ LOH segments (saturating LOH fraction → 1.0)  (§2.2, INV-02);
///   • single SNP             ⇔ a region described by a single segment
///                              ⇒ no divide-by-zero / variance-of-one crash, FINITE result (§6.1 "Empty input").
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-LOH-001 — HRD loss of heterozygosity detection (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 95.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 95): "BAF=0.5 everywhere, BAF=0/1, single SNP".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (Loss_Of_Heterozygosity.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • A segment is LOH iff minor CN == 0 AND major CN != 0
///       (homozygous deletion minor==0,major==0 is NOT LOH)              (§2.2, INV-03, §6.1)
///   • HRD-LOH score = number of LOH regions with length STRICTLY > 15 Mb,
///       after excluding whole-chromosome-LOH chromosomes               (§2.2, §4.2, INV-04, INV-05)
///   • Length exactly 15,000,000 bp ⇒ NOT counted (strict >)            (§6.1)
///   • A chromosome whose every segment is LOH ⇒ whole-chromosome LOH,
///       its regions excluded                                           (§2.2, §6.1, INV-05)
///   • HRD-LOH score ≥ 0                                                 (INV-01)
///   • Score is independent of input segment order                      (INV-06)
///   • CalculateLOHFraction ∈ [0,1]; absent chromosome ⇒ 0.0            (§3.2, §3.3, INV-02)
///   • Null segments / null chromosome ⇒ ArgumentNullException          (§3.3)
///   • End ≤ Start, or a negative copy number ⇒ ArgumentException        (§3.3, §6.1)
///   • Empty input ⇒ score 0; fraction 0                                (§6.1)
///   • HrdLohMinRegionLengthBp == 15,000,000                            (§4.2)
///
/// No source bug was found; no test was weakened.
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyLohFuzzTests
{
    private const long SizeLimit = HrdLohMinRegionLengthBp; // 15,000,000 bp (Abkevich 2012; scarHRD)

    // ── Well-formed-result assertion helpers ─────────────────────────────────
    // Pin the documented numeric contract on EVERY returned value so a fuzz test
    // cannot rubber-stamp a malformed result (negative count / NaN fraction) green.
    private static void AssertWellFormedResult(LohResult result)
    {
        result.Score.Should().BeGreaterThanOrEqualTo(0, "HRD-LOH score is a region count (INV-01)");
        result.Regions.Should().NotBeNull("a result always carries its region list");
        result.Regions.Count.Should().Be(result.Score, "Score == Regions.Count by construction (§3.2)");
        foreach (LohRegion region in result.Regions)
        {
            // INV-03/INV-04: every counted region is a long, valid LOH stretch.
            region.Length.Should().Be(region.End - region.Start, "region length = End − Start");
            region.Length.Should().BeGreaterThan(SizeLimit, "a counted region is strictly > 15 Mb (INV-04)");
        }
    }

    private static void AssertWellFormedFraction(double fraction)
    {
        double.IsNaN(fraction).Should().BeFalse("LOH fraction must never be NaN (INV-02)");
        double.IsInfinity(fraction).Should().BeFalse("LOH fraction must never be ±Inf (INV-02)");
        fraction.Should().BeInRange(0.0, 1.0, "LOH fraction ∈ [0,1] (INV-02, §3.2)");
    }

    // Convenience constructors that name the BAF interpretation at each call site.
    private static AlleleSpecificSegment Het(string chr, long start, long end)
        => new(chr, start, end, MajorCopyNumber: 1, MinorCopyNumber: 1);   // BAF ≈ 0.5

    private static AlleleSpecificSegment Loh(string chr, long start, long end)
        => new(chr, start, end, MajorCopyNumber: 1, MinorCopyNumber: 0);   // BAF = 0 or 1

    #region ONCO-LOH-001 — positive sanity (hand-computed score & fraction)

    [Test]
    public void DetectLOH_DocWorkedExample_ScoreIsOne()
    {
        // Docs §7.1: a 20 Mb LOH on chr1 (counted), a het run on chr1, and a
        // 10 Mb LOH on chr2 (≤ 15 Mb, not counted) ⇒ HRD-LOH score 1.
        var segments = new[]
        {
            new AlleleSpecificSegment("1", 0, 20_000_000, 1, 0),            // 20 Mb LOH → counted
            new AlleleSpecificSegment("1", 20_000_000, 60_000_000, 1, 1),  // het, not LOH
            new AlleleSpecificSegment("2", 0, 10_000_000, 2, 0),           // LOH but ≤ 15 Mb
        };

        LohResult result = DetectLOH(segments);

        AssertWellFormedResult(result);
        result.Score.Should().Be(1, "only the 20 Mb chr1 LOH region qualifies (§7.1)");
        result.Regions.Should().ContainSingle()
            .Which.Should().BeEquivalentTo(new { Chromosome = "1", Start = 0L, End = 20_000_000L, Length = 20_000_000L });
        CalculateHrdLohScore(segments).Should().Be(1, "wrapper agrees with DetectLOH (§5.1)");
    }

    [Test]
    public void CalculateLOHFraction_DocWorkedExample_IsOneThird()
    {
        // Docs §7.1: chr1 LOH 20 Mb of 60 Mb covered ⇒ 20M / 60M = 0.3333…
        var segments = new[]
        {
            new AlleleSpecificSegment("1", 0, 20_000_000, 1, 0),
            new AlleleSpecificSegment("1", 20_000_000, 60_000_000, 1, 1),
            new AlleleSpecificSegment("2", 0, 10_000_000, 2, 0),
        };

        double frac = CalculateLOHFraction(segments, "1");

        AssertWellFormedFraction(frac);
        frac.Should().BeApproximately(20_000_000.0 / 60_000_000.0, 1e-12, "20M LOH of 60M covered (§7.1)");
    }

    [Test]
    public void DetectLOH_HetRegionVsLohRegion_OnlyLohIsCalled()
    {
        // POSITIVE sanity contrast (checklist intent): a heterozygous (BAF≈0.5)
        // region is NOT LOH; a region skewed to one allele (BAF=0/1) IS LOH —
        // provided the chromosome is not WHOLLY LOH (a het anchor keeps it).
        var hetOnly = new[] { Het("1", 0, 30_000_000) };
        var lohWithHetAnchor = new[]
        {
            Loh("1", 0, 30_000_000),                 // 30 Mb LOH → counted
            Het("1", 30_000_000, 31_000_000),        // het anchor so chr1 is not whole-chromosome LOH
        };

        DetectLOH(hetOnly).Score.Should().Be(0, "BAF≈0.5 everywhere ⇒ no LOH (§6.1)");
        DetectLOH(lohWithHetAnchor).Score.Should().Be(1, "a > 15 Mb LOH region IS called (§2.2)");
    }

    #endregion

    #region ONCO-LOH-001 — BE: BAF = 0.5 everywhere (perfectly heterozygous ⇒ no LOH)

    [Test]
    public void DetectLOH_AllHeterozygous_NoLohCalled()
    {
        // BAF ≈ 0.5 everywhere ⇔ minor CN ≠ 0 everywhere ⇒ zero LOH, even for
        // arbitrarily long segments (the size filter is irrelevant when nothing is LOH).
        var segments = new[]
        {
            Het("1", 0, 100_000_000),
            Het("2", 0, 200_000_000),
            new AlleleSpecificSegment("3", 0, 50_000_000, 2, 1), // balanced gain, minor ≠ 0 → not LOH
        };

        LohResult result = DetectLOH(segments);

        AssertWellFormedResult(result);
        result.Score.Should().Be(0, "no segment has minor CN 0 ⇒ no LOH (§6.1 het retained)");
        result.Regions.Should().BeEmpty();
    }

    [Test]
    public void CalculateLOHFraction_AllHeterozygous_IsZero()
    {
        var segments = new[] { Het("X", 0, 75_000_000), Het("X", 75_000_000, 90_000_000) };

        double frac = CalculateLOHFraction(segments, "X");

        AssertWellFormedFraction(frac);
        frac.Should().Be(0.0, "no LOH length over a fully heterozygous chromosome (INV-02)");
    }

    [Test]
    public void DetectLOH_FuzzedHeterozygousGenomes_NeverCallLoh()
    {
        // BE: scatter many het segments (minor CN ≥ 1) across chromosomes of random
        // lengths; the score is invariant at 0 no matter how large or how many.
        var rng = new Random(95_0001);
        for (int trial = 0; trial < 400; trial++)
        {
            int n = rng.Next(1, 12);
            var segments = new List<AlleleSpecificSegment>(n);
            long cursor = 0;
            for (int i = 0; i < n; i++)
            {
                long len = rng.Next(1, 250_000_000);
                int major = rng.Next(1, 6);
                int minor = rng.Next(1, major + 1);     // minor ≥ 1 ⇒ heterozygous (never LOH)
                segments.Add(new AlleleSpecificSegment(((i % 3) + 1).ToString(), cursor, cursor + len, major, minor));
                cursor += len + 1;
            }

            LohResult result = DetectLOH(segments);
            AssertWellFormedResult(result);
            result.Score.Should().Be(0, "minor CN ≥ 1 everywhere ⇒ no LOH (INV-03)");
        }
    }

    #endregion

    #region ONCO-LOH-001 — BE: BAF = 0/1 everywhere (complete allelic imbalance)

    [Test]
    public void DetectLOH_WholeChromosomeLoh_IsExcluded()
    {
        // BAF=0/1 EVERYWHERE on a chromosome ⇒ every segment LOH ⇒ whole-chromosome
        // LOH ⇒ EXCLUDED from the score (Abkevich: "< whole chromosome"; INV-05),
        // even though each segment is individually a huge LOH stretch.
        var segments = new[]
        {
            Loh("1", 0, 100_000_000),
            Loh("1", 100_000_000, 250_000_000),
        };

        LohResult result = DetectLOH(segments);

        AssertWellFormedResult(result);
        result.Score.Should().Be(0, "all segments LOH ⇒ whole-chromosome LOH excluded (§6.1, INV-05)");
    }

    [Test]
    public void CalculateLOHFraction_AllLohSaturatesAtOne_NoOverflow()
    {
        // The fraction has NO whole-chromosome exclusion (§5.2): an all-LOH
        // chromosome saturates the burden at exactly 1.0 — finite, never > 1.
        var segments = new[]
        {
            Loh("7", 0, 1_000_000),
            Loh("7", 1_000_000, 240_000_000),
        };

        double frac = CalculateLOHFraction(segments, "7");

        AssertWellFormedFraction(frac);
        frac.Should().Be(1.0, "all covered length is LOH ⇒ fraction saturates at 1.0 (INV-02)");
    }

    [Test]
    public void DetectLOH_LohWithHetAnchor_CountedWhenStrictlyOver15Mb()
    {
        // A LOH region wins ONLY when the chromosome is not wholly LOH and the
        // region is strictly > 15 Mb. Three sizes around the boundary on the same
        // not-wholly-LOH chromosome.
        long anchorStart = 300_000_000;
        var segments = new[]
        {
            Loh("9", 0, SizeLimit),                                   // exactly 15 Mb → NOT counted
            Loh("9", 20_000_000, 20_000_000 + SizeLimit + 1),        // 15 Mb + 1 bp → counted
            Loh("9", 40_000_000, 40_000_000 + 50_000_000),           // 50 Mb → counted
            Het("9", anchorStart, anchorStart + 1_000_000),          // het anchor (not whole-chromosome LOH)
        };

        LohResult result = DetectLOH(segments);

        AssertWellFormedResult(result);
        result.Score.Should().Be(2, "15 Mb fails strict >, 15 Mb+1 and 50 Mb pass (§6.1, INV-04)");
        result.Regions.Should().OnlyContain(r => r.Length > SizeLimit);
    }

    [Test]
    public void DetectLOH_HomozygousDeletionIsNotLoh()
    {
        // BAF degenerate sibling: minor==0 AND major==0 is a homozygous deletion,
        // explicitly NOT LOH (§2.2, §6.1) — the major != 0 clause excludes it.
        var segments = new[]
        {
            new AlleleSpecificSegment("4", 0, 40_000_000, 0, 0),           // homozygous deletion → not LOH
            new AlleleSpecificSegment("4", 40_000_000, 41_000_000, 1, 1),  // het anchor
        };

        LohResult result = DetectLOH(segments);

        AssertWellFormedResult(result);
        result.Score.Should().Be(0, "minor==0 & major==0 is a homozygous deletion, not LOH (§6.1)");
    }

    [Test]
    public void DetectLOH_FuzzedAllelicImbalance_AlwaysWellFormedNoOverflow()
    {
        // BE: build chromosomes of many LOH segments at huge bp coordinates (near
        // genome scale) plus a het anchor; whatever the verdict, the result is
        // well-formed (non-negative count, every region > 15 Mb), no Int64 overflow
        // and no NaN/Inf anywhere.
        var rng = new Random(95_0002);
        for (int trial = 0; trial < 400; trial++)
        {
            var segments = new List<AlleleSpecificSegment>();
            long cursor = rng.Next(0, 1_000_000);
            int n = rng.Next(1, 8);
            for (int i = 0; i < n; i++)
            {
                long len = rng.Next(1, 60_000_000);
                segments.Add(Loh("12", cursor, cursor + len));
                cursor += len + rng.Next(0, 3);   // sometimes adjacent (≤1bp) → merged
            }
            // Het anchor so the chromosome is not whole-chromosome LOH.
            segments.Add(Het("12", cursor + 10, cursor + 1_000_010));

            LohResult result = DetectLOH(segments);
            AssertWellFormedResult(result);

            double frac = CalculateLOHFraction(segments, "12");
            AssertWellFormedFraction(frac);
        }
    }

    #endregion

    #region ONCO-LOH-001 — BE: single SNP / single-segment region (no divide-by-zero)

    [Test]
    public void DetectLOH_SingleSegment_NoCrash_WellFormed()
    {
        // "Single SNP" ⇔ a region described by a single segment. A solitary LOH
        // segment IS whole-chromosome LOH ⇒ excluded ⇒ score 0; no divide-by-zero
        // or variance-of-one fault on a 1-element per-chromosome stat.
        var loneLoh = new[] { Loh("1", 0, 50_000_000) };
        var loneHet = new[] { Het("1", 0, 50_000_000) };

        DetectLOH(loneLoh).Score.Should().Be(0, "single all-LOH chromosome is whole-chromosome LOH (INV-05)");
        DetectLOH(loneHet).Score.Should().Be(0, "single het segment is not LOH (§6.1)");
    }

    [Test]
    public void CalculateLOHFraction_SingleSegment_NoDivideByZero()
    {
        // A single-segment chromosome must yield a clean 0/1 fraction, never a
        // 0/0 → NaN from a stat computed over one element.
        CalculateLOHFraction(new[] { Loh("1", 0, 1_000_000) }, "1").Should().Be(1.0);
        CalculateLOHFraction(new[] { Het("1", 0, 1_000_000) }, "1").Should().Be(0.0);

        // A single-bp segment (the smallest possible) is still finite and in range.
        double tiny = CalculateLOHFraction(new[] { Loh("1", 0, 1) }, "1");
        AssertWellFormedFraction(tiny);
        tiny.Should().Be(1.0);
    }

    [Test]
    public void DetectLOH_EmptyInput_ScoreZero()
    {
        // §6.1: empty input ⇒ score 0 (empty domain), never a crash.
        LohResult result = DetectLOH(Array.Empty<AlleleSpecificSegment>());

        AssertWellFormedResult(result);
        result.Score.Should().Be(0);
        result.Regions.Should().BeEmpty();
    }

    [Test]
    public void CalculateLOHFraction_EmptyOrAbsentChromosome_IsZero()
    {
        // §3.3 / §6.1: an absent chromosome (no covered length) ⇒ 0.0, the guarded
        // totalLength==0 branch — the single place a divide-by-zero could occur.
        CalculateLOHFraction(Array.Empty<AlleleSpecificSegment>(), "1").Should().Be(0.0);
        CalculateLOHFraction(new[] { Loh("1", 0, 1_000_000) }, "ABSENT").Should().Be(0.0);
    }

    [Test]
    public void CalculateLOHFraction_FuzzedSingleSegment_AlwaysCleanZeroOrOne()
    {
        // BE: random single segments — the per-chromosome fraction is exactly 0
        // (het) or 1 (LOH), always finite, never NaN/Inf from a one-element ratio.
        var rng = new Random(95_0003);
        for (int i = 0; i < 600; i++)
        {
            long start = rng.Next(0, 100_000_000);
            long len = rng.Next(1, 100_000_000);
            bool isLoh = rng.Next(2) == 0;
            var seg = isLoh ? Loh("1", start, start + len) : Het("1", start, start + len);

            double frac = CalculateLOHFraction(new[] { seg }, "1");
            AssertWellFormedFraction(frac);
            frac.Should().Be(isLoh ? 1.0 : 0.0, "a single-segment fraction is a clean 0 or 1");
        }
    }

    #endregion

    #region ONCO-LOH-001 — INV-06 (order independence) & INV-04 boundary

    [Test]
    public void DetectLOH_ScoreIsOrderIndependent_AcrossShuffledFuzzedInputs()
    {
        // INV-06: per-chromosome aggregation is set-based ⇒ the score is invariant
        // under input reordering. Build a mixed genome and shuffle it repeatedly.
        var rng = new Random(95_0004);
        for (int trial = 0; trial < 200; trial++)
        {
            var segments = new List<AlleleSpecificSegment>();
            for (int c = 1; c <= 3; c++)
            {
                long cursor = 0;
                int n = rng.Next(2, 6);
                for (int i = 0; i < n; i++)
                {
                    long len = rng.Next(1, 60_000_000);
                    bool isLoh = rng.Next(2) == 0;
                    var seg = isLoh ? Loh(c.ToString(), cursor, cursor + len)
                                    : Het(c.ToString(), cursor, cursor + len);
                    segments.Add(seg);
                    cursor += len + rng.Next(1, 5);
                }
            }

            int baseline = DetectLOH(segments).Score;
            var shuffled = segments.OrderBy(_ => rng.Next()).ToList();
            DetectLOH(shuffled).Score.Should().Be(baseline, "score is order-independent (INV-06)");
        }
    }

    [Test]
    public void HrdLohMinRegionLengthBp_Is15Mb()
        => SizeLimit.Should().Be(15_000_000L, "Abkevich 2012 / scarHRD sizelimitLOH = 15e6 (§4.2)");

    #endregion

    #region ONCO-LOH-001 — BE: documented exceptions (null / invalid segment)

    [Test]
    public void DetectLOH_NullSegments_ThrowsArgumentNull()
    {
        Action act = () => DetectLOH(null!);
        act.Should().Throw<ArgumentNullException>().Which.ParamName.Should().Be("segments");
    }

    [Test]
    public void CalculateHrdLohScore_NullSegments_ThrowsArgumentNull()
    {
        Action act = () => CalculateHrdLohScore(null!);
        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void CalculateLOHFraction_NullSegmentsOrChromosome_ThrowsArgumentNull()
    {
        Action nullSeg = () => CalculateLOHFraction(null!, "1");
        nullSeg.Should().Throw<ArgumentNullException>();

        Action nullChr = () => CalculateLOHFraction(new[] { Het("1", 0, 1_000_000) }, null!);
        nullChr.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void DetectLOH_NonPositiveLength_ThrowsArgumentException()
    {
        // §3.3 / §6.1: End ≤ Start is an invalid segment.
        Action equal = () => DetectLOH(new[] { new AlleleSpecificSegment("1", 100, 100, 1, 0) });
        equal.Should().Throw<ArgumentException>();

        Action inverted = () => DetectLOH(new[] { new AlleleSpecificSegment("1", 200, 100, 1, 0) });
        inverted.Should().Throw<ArgumentException>();
    }

    [Test]
    public void DetectLOH_NegativeCopyNumber_ThrowsArgumentException()
    {
        // §3.3 / §6.1: a negative copy number is rejected — copy numbers ≥ 0.
        Action negMinor = () => DetectLOH(new[] { new AlleleSpecificSegment("1", 0, 1_000_000, 1, -1) });
        negMinor.Should().Throw<ArgumentException>();

        Action negMajor = () => DetectLOH(new[] { new AlleleSpecificSegment("1", 0, 1_000_000, -1, 0) });
        negMajor.Should().Throw<ArgumentException>();
    }

    [Test]
    public void DetectLOH_FuzzedInvalidSegments_AlwaysThrowArgumentException()
    {
        // BE: inject one out-of-contract segment (bad length or negative CN) into
        // an otherwise valid batch and confirm the documented throw every time.
        var rng = new Random(95_0005);
        for (int i = 0; i < 400; i++)
        {
            var segments = new List<AlleleSpecificSegment>
            {
                Het("1", 0, 1_000_000),
            };

            if (rng.Next(2) == 0)
            {
                long start = rng.Next(0, 1000);
                segments.Add(new AlleleSpecificSegment("2", start, start - rng.Next(0, 1000), 1, 0)); // End ≤ Start
            }
            else
            {
                int badCn = -rng.Next(1, int.MaxValue);
                segments.Add(rng.Next(2) == 0
                    ? new AlleleSpecificSegment("2", 0, 1_000_000, badCn, 0)
                    : new AlleleSpecificSegment("2", 0, 1_000_000, 1, badCn));
            }

            Action act = () => DetectLOH(segments);
            act.Should().Throw<ArgumentException>("an invalid segment is documented out of contract (§3.3)");
        }
    }

    #endregion
}
