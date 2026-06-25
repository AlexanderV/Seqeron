using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Transcriptome alternative-splicing unit — Percent Spliced In
/// (<see cref="TranscriptomeAnalyzer.CalculatePSI"/>) and event classification
/// (<see cref="TranscriptomeAnalyzer.DetectAlternativeSplicing"/>).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// NaN/Infinity leaking where a real number is contracted, and no *unhandled*
/// runtime exception (IndexOutOfRange, NullReference, DivideByZero, Overflow, …).
/// Every input must yield EITHER a well-defined, theory-correct result, OR a
/// *documented, intentional* convention (NaN for the 0/0 undefined PSI case, an
/// ArgumentOutOfRangeException on a negative read count, an empty event sequence
/// for &lt; 2 isoforms / null input). A raw runtime exception, a hang, a PSI outside
/// [0,1], or a fabricated event on a single isoform is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: TRANS-SPLICE-001 — alternative splicing / PSI quantification (Transcriptome)
/// Checklist: docs/checklists/03_FUZZING.md, row 200.
/// Fuzz strategy for THIS unit: BE = Boundary Exploitation (0, -1, MaxInt, empty)
///   — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Fuzz targets (checklist row 200): "single isoform, no junction reads, empty".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (docs/algorithms/Transcriptome/Alternative_Splicing.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • Ψ = I/(I+S), unnormalized read-count ratio. — §2.2, §4.1.
///   • rMATS length-normalized ψ̂ = (I/l_I)/(I/l_I + S/l_S) when BOTH effective
///     lengths are strictly positive; otherwise the unnormalized ratio. — §2.2, §3.3, §4.1.
///   • DetectAlternativeSplicing groups by gene, skips genes with &lt; 2 isoforms,
///     emits no event for identical isoform pairs, treats null as empty. — §3.3, §4.1.
///
/// Documented invariants this fixture pins (§2.4):
///   • INV-01: 0 ≤ Ψ ≤ 1 when I,S ≥ 0 and I+S &gt; 0.
///   • INV-02: I+S = 0 ⇒ Ψ = NaN (0/0 undefined — no divide-by-zero, no Infinity).
///   • INV-03: S=0,I&gt;0 ⇒ Ψ = 1 ; I=0,S&gt;0 ⇒ Ψ = 0.
///   • INV-04: both lengths &gt; 0 ⇒ Ψ = (I/l_I)/(I/l_I + S/l_S).
///   • INV-05: a detected event references two isoforms of one gene differing in structure.
///
/// Boundary handling fixed by the doc (§3.3, §6.1) and pinned here so the contract
/// can never silently drift:
///   • SINGLE ISOFORM (BE): a gene with &lt; 2 isoforms ⇒ NO event (an AS event needs
///     two isoforms). — §6.1 ("&lt; 2 isoforms for a gene"), INV-05.
///   • NO JUNCTION READS (BE, "zero-coverage"): I = S = 0 ⇒ Ψ = NaN by the documented
///     0/0 undefined convention; no DivideByZero, no Infinity. The rMATS form with
///     positive lengths but zero reads is ALSO 0/0 ⇒ NaN. — §6.1 ("I+S = 0"), INV-02.
///   • EMPTY / NULL isoform enumerable (BE) ⇒ empty event result. — §3.3, §6.1.
///   • NEGATIVE reads (BE, "-1") ⇒ ArgumentOutOfRangeException (counts are
///     non-negative; thrown, not swallowed). — §3.3, §6.1 ("negative reads").
///   • MaxInt-scale reads ⇒ a finite Ψ ∈ [0,1], no Overflow. — BE.
///
/// Positive sanity (worked example, §7.1, derived INDEPENDENTLY from the formulae,
/// NOT echoed off the implementation):
///   • Ψ(80, 20) = 80/(80+20) = 0.80.
///   • rMATS Ψ(80, 20, l_I=200, l_S=100) = (80/200)/((80/200)+(20/100))
///       = 0.4/(0.4+0.2) = 0.4/0.6 = 2/3 = 0.6666…  (corrects the longer-isoform bias).
///   • A skipped-exon isoform pair (one isoform missing one internal exon) is
///     classified as exactly one SkippedExon event.
/// These pin a non-degenerate contract: a stub returning a constant Ψ (e.g. always
/// 0.5) or never emitting an event would fail here even though it "never crashes".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Determinism
/// ───────────────────────────────────────────────────────────────────────────
/// All inputs are hand-built or generated from a LOCALLY fixed-seed
/// `new Random(seed)` (never a shared static Rng), so every run is reproducible.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class TranscriptomeSplicingFuzzTests
{
    private const double Tol = 1e-12;

    // Worked-example constants (§7.1), derived from the formula, NOT echoed off code.
    private const double PsiWorked = 0.80;                 // 80/(80+20)
    private const double PsiWorkedNormalized = 2.0 / 3.0;  // (80/200)/((80/200)+(20/100))

    private static TranscriptomeAnalyzer.TranscriptIsoform Iso(
        string transcriptId, string geneId, params (int Start, int End)[] exons)
        => new(transcriptId, geneId,
               Length: exons.Sum(e => e.End - e.Start + 1),
               ExonCount: exons.Length, Expression: 1.0, IsProteinCoding: true,
               Exons: exons.ToList());

    #region TRANS-SPLICE-001 — positive sanity (worked example must be exact)

    // ════════════════════════════════════════════════════════════════════════
    //  Guards against a degenerate analyzer (constant Ψ / never-emit) that would
    //  pass every boundary test below. Expected values derived from the §7.1
    //  formula, cross-checked by hand, never read off the implementation.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void CalculatePSI_WorkedExample_EqualsZeroPointEight()
    {
        double psi = TranscriptomeAnalyzer.CalculatePSI(80, 20);

        psi.Should().BeApproximately(PsiWorked, Tol,
            "Ψ = I/(I+S) = 80/(80+20) = 0.80 — §7.1 unnormalized read-count definition");
    }

    [Test]
    public void CalculatePSI_WorkedExample_LengthNormalized_EqualsTwoThirds()
    {
        // rMATS: I=80, S=20, l_I=200, l_S=100 ⇒ (80/200)/((80/200)+(20/100)) = 0.4/0.6 = 2/3.
        double psi = TranscriptomeAnalyzer.CalculatePSI(80, 20, 200, 100);

        psi.Should().BeApproximately(PsiWorkedNormalized, Tol,
            "ψ̂ = (I/l_I)/(I/l_I + S/l_S) = 0.4/0.6 = 2/3 — §7.1 rMATS length-normalized form (§2.2 INV-04)");
    }

    // INV-03: full inclusion / full exclusion are exactly 1 and 0.
    [Test]
    public void CalculatePSI_FullInclusionOrExclusion_IsExactlyOneOrZero()
    {
        TranscriptomeAnalyzer.CalculatePSI(50, 0).Should().Be(1.0,
            "S = 0, I > 0 ⇒ Ψ = I/I = 1 (INV-03)");
        TranscriptomeAnalyzer.CalculatePSI(0, 50).Should().Be(0.0,
            "I = 0, S > 0 ⇒ Ψ = 0/S = 0 (INV-03)");
    }

    // INV-05 / §4.1: a clean skipped-exon pair yields exactly one SkippedExon event.
    [Test]
    public void DetectAlternativeSplicing_SkippedExonPair_ClassifiesAsSingleSkippedExon()
    {
        // Inclusion isoform has the middle exon [200,300]; skipping isoform lacks it.
        var inclusion = Iso("t-inc", "geneX", (100, 150), (200, 300), (400, 450));
        var skipping  = Iso("t-skip", "geneX", (100, 150), (400, 450));

        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(new[] { inclusion, skipping }).ToList();

        events.Should().ContainSingle("one structural difference between the two isoforms (§4.1)");
        events[0].EventType.Should().Be("SkippedExon", "an exon present in one isoform and absent in the other is a cassette/skipped exon (§4.2)");
        events[0].GeneId.Should().Be("geneX", "the event references the shared gene (INV-05)");
    }

    #endregion

    #region TRANS-SPLICE-001 — BE boundary: SINGLE ISOFORM (< 2 isoforms ⇒ no event)

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: SINGLE ISOFORM (BE). An AS event is defined per isoform PAIR
    // of a gene, so a gene with fewer than two isoforms must produce NO event —
    // never a fabricated self-event. — §6.1 ("< 2 isoforms for a gene"), INV-05.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void DetectAlternativeSplicing_SingleIsoform_ProducesNoEvent()
    {
        var only = Iso("solo", "geneSolo", (1, 100), (200, 300), (400, 500));

        TranscriptomeAnalyzer.DetectAlternativeSplicing(new[] { only })
            .Should().BeEmpty("a single isoform has no pair to compare ⇒ no AS event (§6.1, INV-05)");
    }

    // Two genes, each with exactly one isoform: still no event for either.
    [Test]
    public void DetectAlternativeSplicing_OneIsoformPerGene_ProducesNoEvent()
    {
        var a = Iso("a1", "geneA", (1, 100), (200, 300));
        var b = Iso("b1", "geneB", (1, 100), (200, 300));

        TranscriptomeAnalyzer.DetectAlternativeSplicing(new[] { a, b })
            .Should().BeEmpty("events are within a gene; one isoform per gene ⇒ no pair (§4.1)");
    }

    // Two structurally IDENTICAL isoforms of one gene ⇒ still no event.
    [Test]
    public void DetectAlternativeSplicing_IdenticalIsoformPair_ProducesNoEvent()
    {
        var a = Iso("x1", "geneI", (1, 100), (200, 300), (400, 500));
        var b = Iso("x2", "geneI", (1, 100), (200, 300), (400, 500));

        TranscriptomeAnalyzer.DetectAlternativeSplicing(new[] { a, b })
            .Should().BeEmpty("identical exon structure ⇒ no structural difference ⇒ no event (§6.1)");
    }

    #endregion

    #region TRANS-SPLICE-001 — BE boundary: NO JUNCTION READS (zero-coverage I=S=0)

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: NO JUNCTION READS (BE, "zero-coverage"). I + S = 0 is the
    // documented 0/0 undefined case ⇒ Ψ = NaN. Pinned EXPLICITLY because the
    // alternative (a swallowed 0 or a thrown DivideByZero) would corrupt the
    // contract. — §6.1 ("I+S = 0"), INV-02. No Infinity, no exception.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculatePSI_NoJunctionReads_IsNaNNotDivideByZero()
    {
        double psi = TranscriptomeAnalyzer.CalculatePSI(0, 0);

        double.IsNaN(psi).Should().BeTrue("I+S = 0 ⇒ 0/0 undefined ⇒ Ψ = NaN (INV-02, §6.1)");
        double.IsInfinity(psi).Should().BeFalse("0/0 is NaN, never ±Infinity");
    }

    // rMATS form with positive lengths but zero reads is ALSO 0/0 ⇒ NaN, not 0/L.
    [Test]
    public void CalculatePSI_NoJunctionReads_LengthNormalized_IsNaN()
    {
        double psi = TranscriptomeAnalyzer.CalculatePSI(0, 0, 200, 100);

        double.IsNaN(psi).Should().BeTrue(
            "rMATS rates 0/200 and 0/100 both 0 ⇒ denominator 0 ⇒ ψ̂ = NaN (INV-02)");
        double.IsInfinity(psi).Should().BeFalse();
    }

    #endregion

    #region TRANS-SPLICE-001 — BE boundary: EMPTY / NULL isoform enumerable

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY / NULL isoform enumerable (BE) ⇒ empty event result,
    // no NullReferenceException. — §3.3 ("treats null input as an empty
    // sequence"), §6.1 ("null isoforms ⇒ empty result").
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void DetectAlternativeSplicing_EmptyOrNullIsoforms_YieldEmptyResult()
    {
        TranscriptomeAnalyzer.DetectAlternativeSplicing(
            Array.Empty<TranscriptomeAnalyzer.TranscriptIsoform>())
            .Should().BeEmpty("no isoforms ⇒ no events (§6.1)");

        TranscriptomeAnalyzer.DetectAlternativeSplicing(null!)
            .Should().BeEmpty("null isoform enumerable ⇒ empty result, no NRE (§3.3, §6.1)");
    }

    // An isoform whose Exons list itself is null must not crash classification.
    [Test]
    public void DetectAlternativeSplicing_IsoformWithNullExons_DoesNotCrash()
    {
        var withExons = Iso("e1", "geneN", (1, 100), (200, 300));
        var nullExons = new TranscriptomeAnalyzer.TranscriptIsoform(
            "e2", "geneN", Length: 0, ExonCount: 0, Expression: 1.0,
            IsProteinCoding: true, Exons: null!);

        Action act = () => TranscriptomeAnalyzer
            .DetectAlternativeSplicing(new[] { withExons, nullExons }).ToList();

        act.Should().NotThrow("null Exons is treated as the empty exon set, not a crash (§3.3)");
    }

    #endregion

    #region TRANS-SPLICE-001 — BE boundary: NEGATIVE reads (-1) ⇒ thrown

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: NEGATIVE reads (BE, "-1"). Read counts are non-negative;
    // a negative count is rejected with ArgumentOutOfRangeException — a thrown,
    // documented precondition violation, NOT a silent NaN or a garbage ratio.
    // — §3.3 ("throws ArgumentOutOfRangeException if either read count is
    // negative"), §6.1 ("negative reads").
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculatePSI_NegativeInclusion_Throws()
    {
        Action act = () => TranscriptomeAnalyzer.CalculatePSI(-1, 20);
        act.Should().Throw<ArgumentOutOfRangeException>("negative read counts violate the non-negativity precondition (§3.3)");
    }

    [Test]
    public void CalculatePSI_NegativeExclusion_Throws()
    {
        Action act = () => TranscriptomeAnalyzer.CalculatePSI(20, -1);
        act.Should().Throw<ArgumentOutOfRangeException>("negative read counts violate the non-negativity precondition (§3.3)");
    }

    #endregion

    #region TRANS-SPLICE-001 — randomized boundary sweep

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: random boundary batch (BE) under a time budget.
    // A deterministic, locally-seeded generator draws read counts and effective
    // lengths spanning the boundaries (0, MaxInt-scale, mixed zero) and feeds
    // both the unnormalized and rMATS forms. Every call must satisfy:
    //   • non-negative I,S with I+S > 0 ⇒ Ψ finite ∈ [0,1] (INV-01);
    //   • I+S = 0 (no junction reads) ⇒ Ψ = NaN, never ±Infinity (INV-02);
    //   • negative count ⇒ ArgumentOutOfRangeException (thrown, §3.3);
    //   • MaxInt-scale counts ⇒ finite Ψ ∈ [0,1], no Overflow.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void CalculatePSI_RandomBoundaryBatch_NeverCrashesAndStaysWellFormed()
    {
        var rng = new Random(20260620); // locally fixed seed — deterministic

        for (int trial = 0; trial < 5000; trial++)
        {
            // Read counts: include exact 0, ordinary, and MaxInt-scale magnitudes.
            double i = PickCount(rng);
            double s = PickCount(rng);
            // Effective lengths: 0 (⇒ unnormalized path) or positive (⇒ rMATS path).
            double li = rng.Next(2) == 0 ? 0 : rng.NextDouble() * 1000 + 1;
            double ls = rng.Next(2) == 0 ? 0 : rng.NextDouble() * 1000 + 1;

            double psi = TranscriptomeAnalyzer.CalculatePSI(i, s, li, ls);

            bool noReads = (li > 0 && ls > 0)
                ? (i / li) + (s / ls) == 0   // rMATS denominator zero
                : i + s == 0;                // unnormalized denominator zero

            if (noReads)
            {
                double.IsNaN(psi).Should().BeTrue("no junction reads ⇒ 0/0 ⇒ Ψ = NaN (INV-02)");
            }
            else
            {
                double.IsNaN(psi).Should().BeFalse("with supporting reads Ψ is defined (INV-01)");
                double.IsInfinity(psi).Should().BeFalse("Ψ is a bounded ratio, never Infinity (INV-01)");
                psi.Should().BeInRange(0.0, 1.0, "Ψ is a part/whole fraction ⇒ 0 ≤ Ψ ≤ 1 (INV-01)");
            }
        }
    }

    // Negative counts in the sweep must always throw, regardless of the other args.
    [Test]
    [CancelAfter(30000)]
    public void CalculatePSI_RandomNegativeCounts_AlwaysThrow()
    {
        var rng = new Random(987654321);

        for (int trial = 0; trial < 2000; trial++)
        {
            double mag = rng.NextDouble() * int.MaxValue + double.Epsilon;
            bool negI = rng.Next(2) == 0;
            double i = negI ? -mag : rng.NextDouble() * 1000;
            double s = negI ? rng.NextDouble() * 1000 : -mag;

            Action act = () => TranscriptomeAnalyzer.CalculatePSI(i, s);
            act.Should().Throw<ArgumentOutOfRangeException>(
                "any negative read count is rejected (§3.3, BE -1 boundary)");
        }
    }

    // Pick a count that is exactly 0 sometimes, otherwise ordinary or MaxInt-scale.
    private static double PickCount(Random rng) => rng.Next(4) switch
    {
        0 => 0.0,                                  // zero-coverage boundary
        1 => rng.NextDouble() * 100.0,             // ordinary
        2 => (double)int.MaxValue,                 // MaxInt boundary (BE)
        _ => rng.NextDouble() * int.MaxValue,      // large magnitude
    };

    #endregion
}
