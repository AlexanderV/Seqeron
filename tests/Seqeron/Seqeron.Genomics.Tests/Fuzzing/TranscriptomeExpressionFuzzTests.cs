using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Transcriptome expression-quantification unit — TPM
/// (<see cref="TranscriptomeAnalyzer.CalculateTPM"/>), FPKM/RPKM
/// (<see cref="TranscriptomeAnalyzer.CalculateFPKM"/>) and quantile normalization
/// (<see cref="TranscriptomeAnalyzer.QuantileNormalize"/>).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// NaN/Infinity leaking into a TPM/FPKM/normalized value, and no *unhandled*
/// runtime exception (IndexOutOfRange, NullReference, DivideByZero, Overflow, …).
/// Every input must yield EITHER a well-defined, theory-correct result, OR a
/// *documented, intentional* convention (empty sequence, TPM = 0, FPKM = 0). A raw
/// runtime exception, a hang, a TPM that does not sum to 10^6 on a non-degenerate
/// sample, or a NaN/Infinity is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: TRANS-EXPR-001 — expression quantification (Transcriptome)
/// Checklist: docs/checklists/03_FUZZING.md, row 199.
/// Fuzz strategy for THIS unit: BE = Boundary Exploitation (0, -1, MaxInt, empty)
///   — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Fuzz targets (checklist row 199): "zero reads, single transcript, all-multimapped".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (docs/algorithms/Transcriptome/Expression_Quantification.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • TPM_i = (X_i / l_i) / Σ_j(X_j / l_j) · 10^6. — §2.2 (verbatim from [2][3]).
///   • FPKM_i = X_i · 10^9 / (l_i · N); 0 when l ≤ 0 or N ≤ 0. — §2.2, §3.3, §6.1.
///   • Quantile normalization: replace each value by the across-sample mean of the
///     sorted value at its rank, preserving original positions; tied values inside a
///     sample share the average of the rank means they span (Bolstad tie rule).
///     — §2.2, §4.1.
///
/// Documented invariants this fixture pins (§2.4):
///   • INV-01: Σ_i TPM_i = 10^6 for any sample with Σ(X/l) > 0.
///   • INV-02: TPM_i ≥ 0; equal X/l ⇒ equal TPM.
///   • INV-03: FPKM_i ≥ 0; FPKM = 0 when l ≤ 0 or N ≤ 0.
///   • INV-04: quantile normalization preserves within-column rank order.
///   • INV-05: untied columns map to a permutation of the same rank-mean multiset.
///
/// Boundary handling fixed by the doc (§3.3, §6.1, §5.4) and pinned here so the
/// contract can never silently drift (the row's checklist targets map onto the
/// raw-count model of this unit — it consumes already-assigned counts, §3.1):
///   • ZERO READS (BE): every transcript count 0 ⇒ Σ(X/l) = 0 (0/0 undefined) ⇒
///     TPM = 0 and FPKM = 0 for every gene; no DivideByZero, no NaN. — §3.3, §6.1
///     ("All-zero counts"), §5.4 Assumption 2.
///   • SINGLE TRANSCRIPT (BE): one gene with a positive count ⇒ its rate is the
///     whole denominator ⇒ TPM = 10^6 (the full per-million budget). — INV-01.
///   • ALL-MULTIMAPPED (BE): reads that map equally to every transcript appear as a
///     uniform count across genes; equal X with equal l ⇒ equal X/l ⇒ every TPM is
///     10^6 / n (the budget split evenly), summing to 10^6. — INV-01/INV-02, §6.1.
///   • EMPTY geneCounts / EMPTY samples ⇒ empty sequence, no crash. — §3.3, §6.1.
///   • FPKM with length ≤ 0 or N ≤ 0 (BE: 0, −1) ⇒ 0 (formula undefined). — INV-03.
///
/// Positive sanity (§7.1 worked example, derived INDEPENDENTLY from the formulae,
/// NOT echoed off the implementation):
///   • TPM: counts/lengths (10/2000, 20/4000, 30/1000) ⇒ rates (0.005, 0.005, 0.030),
///       Σ = 0.04 ⇒ TPM = rate/0.04·10^6 = (125000, 125000, 750000), Σ = 10^6.
///   • FPKM(1000, 2000, 1_000_000) = 1000·10^9 / (2000·10^6) = 10^12 / (2·10^9) = 500.
/// A real, hand-checkable quantification result must therefore appear, so a passing
/// "no crash" result cannot come from a degenerate analyzer that returns 0 for all.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Determinism
/// ───────────────────────────────────────────────────────────────────────────
/// All inputs are hand-built or generated from a LOCALLY fixed-seed
/// `new Random(seed)` (never a shared static Rng), so every run is reproducible.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class TranscriptomeExpressionFuzzTests
{
    private const double TpmTotal = 1_000_000.0;

    private static (string GeneId, double RawCount, int Length) Gene(string id, double count, int length)
        => (id, count, length);

    #region TRANS-EXPR-001 — positive sanity (§7.1 worked example)

    // ════════════════════════════════════════════════════════════════════════
    //  The §7.1 worked example must be reproduced EXACTLY from the formula.
    //  Guards against a degenerate analyzer (constant 0 / empty) that would pass
    //  every boundary test below.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void CalculateTPM_WorkedExample_MatchesHandComputedValuesAndSumsToMillion()
    {
        // rates: 10/2000 = 0.005, 20/4000 = 0.005, 30/1000 = 0.030; Σ = 0.04.
        // TPM = rate/0.04·10^6 ⇒ A = 125000, B = 125000, C = 750000.
        var tpm = TranscriptomeAnalyzer.CalculateTPM(new[]
        {
            Gene("A", 10.0, 2000),
            Gene("B", 20.0, 4000),
            Gene("C", 30.0, 1000),
        }).ToList();

        tpm.Single(g => g.GeneId == "A").TPM.Should().BeApproximately(125_000.0, 1e-6,
            "0.005/0.04·10^6 = 125000 — §7.1, derived from TPM_i = (X/l)/Σ(X/l)·10^6");
        tpm.Single(g => g.GeneId == "B").TPM.Should().BeApproximately(125_000.0, 1e-6, "0.005/0.04·10^6");
        tpm.Single(g => g.GeneId == "C").TPM.Should().BeApproximately(750_000.0, 1e-6, "0.030/0.04·10^6");

        tpm.Sum(g => g.TPM).Should().BeApproximately(TpmTotal, 1e-6, "Σ TPM = 10^6 (INV-01)");
    }

    [Test]
    public void CalculateFPKM_WorkedExample_EqualsFiveHundred()
    {
        // 1000·10^9 / (2000·10^6) = 10^12 / (2·10^9) = 500.
        TranscriptomeAnalyzer.CalculateFPKM(1000, 2000, 1_000_000)
            .Should().BeApproximately(500.0, 1e-9,
                "X·10^9/(l·N) = 1000·10^9/(2000·10^6) = 500 — §7.1, derived from the FPKM formula");
    }

    // INV-05: an untied two-column quantile normalization reproduces the canonical
    // Bolstad worked result; rank means are a permutation of the same multiset.
    [Test]
    public void QuantileNormalize_TwoUntiedColumns_AssignsRankMeansByRank()
    {
        // Sample 1 = (5, 2, 3, 4); Sample 2 = (4, 1, 4, 2) (cols).
        // Sorted col1 = (2,3,4,5); sorted col2 = (1,2,4,4).
        // rank means = ((2+1)/2,(3+2)/2,(4+4)/2,(5+4)/2) = (1.5, 2.5, 4.0, 4.5).
        // col1: 5→rank3→4.5, 2→rank0→1.5, 3→rank1→2.5, 4→rank2→4.0.
        var result = TranscriptomeAnalyzer.QuantileNormalize(new[]
        {
            new double[] { 5, 2, 3, 4 },
            new double[] { 4, 1, 4, 2 },
        }).ToList();

        result[0].Should().Equal(new double[] { 4.5, 1.5, 2.5, 4.0 },
            "untied column maps each value to the rank mean at its rank (INV-04/INV-05)");
    }

    #endregion

    #region TRANS-EXPR-001 — BE boundary: ZERO READS (all-zero counts)

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: ZERO READS (BE). Every transcript count 0 ⇒ Σ(X/l) = 0 (0/0
    // undefined) ⇒ TPM = 0 AND FPKM = 0 for every gene; no DivideByZero, no NaN.
    // — §3.3, §6.1 ("All-zero counts"), §5.4 Assumption 2.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateTPM_ZeroReads_AllTpmAndFpkmZeroNoNaN()
    {
        var tpm = TranscriptomeAnalyzer.CalculateTPM(new[]
        {
            Gene("A", 0.0, 2000),
            Gene("B", 0.0, 1500),
            Gene("C", 0.0, 800),
        }).ToList();

        tpm.Should().HaveCount(3, "one record per gene even with zero counts");
        foreach (var g in tpm)
        {
            g.TPM.Should().Be(0.0, "Σ(X/l) = 0 ⇒ 0/0 undefined ⇒ convention TPM = 0 (§6.1)");
            g.FPKM.Should().Be(0.0, "total mapped reads N = 0 ⇒ FPKM = 0 (INV-03)");
            double.IsNaN(g.TPM).Should().BeFalse("no 0/0 NaN leaks out");
            double.IsInfinity(g.TPM).Should().BeFalse("no divide-by-zero Infinity");
        }
    }

    [Test]
    public void CalculateFPKM_ZeroReads_IsZero()
    {
        TranscriptomeAnalyzer.CalculateFPKM(0.0, 2000, 1_000_000)
            .Should().Be(0.0, "X = 0 ⇒ FPKM = 0 (numerator 0, l and N positive)");
    }

    #endregion

    #region TRANS-EXPR-001 — BE boundary: SINGLE TRANSCRIPT

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: SINGLE TRANSCRIPT (BE). One gene with a positive count owns the
    // whole denominator ⇒ TPM = 10^6 (the full per-million budget). — INV-01.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateTPM_SingleTranscript_GetsFullMillion()
    {
        var tpm = TranscriptomeAnalyzer.CalculateTPM(new[] { Gene("only", 42.0, 1337) }).ToList();

        tpm.Should().ContainSingle();
        tpm[0].TPM.Should().BeApproximately(TpmTotal, 1e-6,
            "single non-zero gene is the entire Σ(X/l) ⇒ TPM = 10^6 (INV-01)");
        // FPKM for the single gene with N = its own count: X·10^9/(l·X) = 10^9/l.
        tpm[0].FPKM.Should().BeApproximately(FpkmScaling / 1337.0, 1e-3,
            "N = own count ⇒ FPKM = 10^9/l (INV-03 form)");
    }

    [Test]
    public void CalculateTPM_SingleZeroTranscript_IsZeroNotMillion()
    {
        // A lone gene with count 0 is the all-zero degenerate case ⇒ TPM 0, not 10^6.
        var tpm = TranscriptomeAnalyzer.CalculateTPM(new[] { Gene("only", 0.0, 500) }).ToList();

        tpm.Should().ContainSingle();
        tpm[0].TPM.Should().Be(0.0, "Σ(X/l) = 0 even for one gene ⇒ TPM = 0 (§6.1)");
        tpm[0].FPKM.Should().Be(0.0, "N = 0 ⇒ FPKM = 0 (INV-03)");
    }

    #endregion

    #region TRANS-EXPR-001 — BE boundary: ALL-MULTIMAPPED (uniform counts)

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: ALL-MULTIMAPPED (BE). Reads mapping equally to every transcript
    // surface as a uniform count across genes. Equal X with equal l ⇒ equal X/l ⇒
    // every TPM = 10^6/n (the budget split evenly), still summing to 10^6.
    // — INV-01/INV-02, §6.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateTPM_AllMultimappedUniformEqualLength_SplitsMillionEvenly()
    {
        // 5 transcripts, identical count and identical length ⇒ identical X/l ⇒ equal TPM.
        var tpm = TranscriptomeAnalyzer.CalculateTPM(Enumerable.Range(0, 5)
            .Select(i => Gene($"t{i}", 17.0, 900))).ToList();

        tpm.Should().HaveCount(5);
        double expected = TpmTotal / 5.0; // 200000 each.
        foreach (var g in tpm)
            g.TPM.Should().BeApproximately(expected, 1e-6,
                "equal X/l ⇒ equal TPM = 10^6/n (INV-02), each 200000");
        tpm.Sum(g => g.TPM).Should().BeApproximately(TpmTotal, 1e-6, "Σ TPM = 10^6 (INV-01)");
    }

    [Test]
    public void CalculateTPM_AllMultimappedUniformCountVaryingLength_WeightsByInverseLength()
    {
        // Same count everywhere but different lengths: X/l is larger for shorter l, so
        // shorter transcripts receive proportionally more TPM (length normalization).
        var tpm = TranscriptomeAnalyzer.CalculateTPM(new[]
        {
            Gene("short", 12.0, 1000), // rate 0.012
            Gene("long",  12.0, 3000), // rate 0.004
        }).ToList();

        // Σrate = 0.016 ⇒ short = 0.012/0.016·10^6 = 750000, long = 250000.
        tpm.Single(g => g.GeneId == "short").TPM.Should().BeApproximately(750_000.0, 1e-6,
            "shorter transcript gets more TPM per the X/l weighting (INV-02)");
        tpm.Single(g => g.GeneId == "long").TPM.Should().BeApproximately(250_000.0, 1e-6,
            "longer transcript gets less TPM (length normalization)");
        tpm.Sum(g => g.TPM).Should().BeApproximately(TpmTotal, 1e-6, "Σ TPM = 10^6 (INV-01)");
    }

    #endregion

    #region TRANS-EXPR-001 — BE boundary: EMPTY / non-positive inputs

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY geneCounts / EMPTY samples (BE) ⇒ empty sequence, no
    // crash. — §3.3, §6.1 ("Empty geneCounts / samples").
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateTPM_EmptyGeneCounts_YieldsEmpty()
    {
        TranscriptomeAnalyzer.CalculateTPM(
            Array.Empty<(string, double, int)>())
            .Should().BeEmpty("no genes ⇒ empty sequence (§6.1)");
    }

    [Test]
    public void QuantileNormalize_EmptySamplesOrEmptyColumns_YieldsEmpty()
    {
        TranscriptomeAnalyzer.QuantileNormalize(Array.Empty<IEnumerable<double>>())
            .Should().BeEmpty("no samples ⇒ empty sequence (§6.1)");

        TranscriptomeAnalyzer.QuantileNormalize(new[] { Enumerable.Empty<double>() })
            .Should().BeEmpty("zero genes per sample ⇒ empty sequence (§3.3)");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: FPKM with non-positive length / N (BE: 0, −1) ⇒ 0 (the formula
    // is undefined for l ≤ 0 or N ≤ 0). — INV-03, §3.3, §6.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateFPKM_NonPositiveLengthOrDepth_ReturnsZero()
    {
        TranscriptomeAnalyzer.CalculateFPKM(1000, 0, 1_000_000)
            .Should().Be(0.0, "length = 0 ⇒ undefined ⇒ 0 (INV-03)");
        TranscriptomeAnalyzer.CalculateFPKM(1000, -1, 1_000_000)
            .Should().Be(0.0, "length = −1 ⇒ undefined ⇒ 0 (BE −1, INV-03)");
        TranscriptomeAnalyzer.CalculateFPKM(1000, 2000, 0)
            .Should().Be(0.0, "N = 0 ⇒ undefined ⇒ 0 (INV-03)");
        TranscriptomeAnalyzer.CalculateFPKM(1000, 2000, -1)
            .Should().Be(0.0, "N = −1 ⇒ undefined ⇒ 0 (BE −1, INV-03)");
    }

    // QN with a single sample: rank means equal that sample's own sorted values, so
    // the output must equal the input exactly (mean of one value is the value).
    [Test]
    public void QuantileNormalize_SingleSample_ReturnsInputUnchanged()
    {
        var input = new double[] { 9, 3, 7, 1 };
        var result = TranscriptomeAnalyzer.QuantileNormalize(new[] { input }).ToList();

        result.Should().ContainSingle();
        result[0].Should().Equal(input, "one sample ⇒ rank means = own values ⇒ identity (§6.1)");
    }

    // QN with identical columns must return the columns unchanged (§6.1).
    [Test]
    public void QuantileNormalize_IdenticalColumns_ReturnsThemUnchanged()
    {
        var col = new double[] { 4, 1, 9, 2 };
        var result = TranscriptomeAnalyzer.QuantileNormalize(new[]
        {
            col.ToArray(), col.ToArray(), col.ToArray(),
        }).ToList();

        foreach (var r in result)
            r.Should().Equal(col, "mean of equal values is the value ⇒ output = input (§6.1)");
    }

    #endregion

    #region TRANS-EXPR-001 — randomized boundary sweep

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: random boundary batch (BE) under a time budget. A deterministic,
    // locally-seeded generator builds gene panels spanning the boundaries (zero
    // reads, single transcript, all-multimapped uniform, ±MaxInt-scale magnitudes,
    // empty panels). CalculateTPM must process every case without crashing or
    // hanging, and EVERY result must be well-formed:
    //   • TPM and FPKM finite (no NaN/Infinity) and ≥ 0 (INV-02/INV-03);
    //   • a non-degenerate sample (Σ(X/l) > 0) ⇒ Σ TPM = 10^6 (INV-01);
    //   • an all-zero sample ⇒ every TPM = 0 (§6.1).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void CalculateTPM_RandomBoundaryBatch_NeverCrashesAndStaysWellFormed()
    {
        var rng = new Random(20260621); // locally fixed seed — deterministic

        for (int trial = 0; trial < 400; trial++)
        {
            int geneCount = rng.Next(0, 9);
            var genes = new List<(string, double, int)>(geneCount);
            bool forcedAllZero = geneCount > 0 && rng.Next(4) == 0;

            for (int g = 0; g < geneCount; g++)
            {
                int length = rng.Next(6) switch
                {
                    0 => 1,                       // shortest meaningful length
                    1 => int.MaxValue,            // BE: MaxInt length
                    _ => rng.Next(1, 50_000),     // ordinary length
                };

                double count = forcedAllZero
                    ? 0.0
                    : rng.Next(5) switch
                    {
                        0 => 0.0,                                 // zero reads for this gene
                        1 => (double)int.MaxValue,                // BE: MaxInt count
                        2 => 17.0,                                // uniform (all-multimapped flavour)
                        _ => rng.NextDouble() * 1e6,              // ordinary count
                    };

                genes.Add(($"g{g}", count, length));
            }

            var results = TranscriptomeAnalyzer.CalculateTPM(genes).ToList();

            results.Should().HaveCount(geneCount, "exactly one record per gene");

            foreach (var r in results)
            {
                double.IsNaN(r.TPM).Should().BeFalse("TPM never NaN on boundary input");
                double.IsInfinity(r.TPM).Should().BeFalse("TPM never Infinity (no divide-by-zero)");
                r.TPM.Should().BeGreaterThanOrEqualTo(0.0, "TPM ≥ 0 (INV-02)");
                double.IsNaN(r.FPKM).Should().BeFalse("FPKM never NaN");
                double.IsInfinity(r.FPKM).Should().BeFalse("FPKM never Infinity");
                r.FPKM.Should().BeGreaterThanOrEqualTo(0.0, "FPKM ≥ 0 (INV-03)");
            }

            if (geneCount == 0)
                continue;

            double sumRates = genes.Sum(g => g.Item2 / Math.Max(g.Item3, 1));
            if (sumRates == 0)
            {
                results.Should().OnlyContain(r => r.TPM == 0.0,
                    "Σ(X/l) = 0 ⇒ every TPM = 0 (§6.1 all-zero convention)");
            }
            else
            {
                // Σ TPM = 10^6 (INV-01); relative tolerance for MaxInt-scale magnitudes.
                results.Sum(r => r.TPM).Should().BeApproximately(TpmTotal, TpmTotal * 1e-6,
                    "Σ TPM = 10^6 for any sample with Σ(X/l) > 0 (INV-01)");
            }
        }
    }

    #endregion

    // FPKM 10^9 scaling factor, derived from the doc (§2.2: 10^3 per-kb · 10^6
    // per-million), NOT echoed off the implementation constant.
    private const double FpkmScaling = 1_000_000_000.0;
}
