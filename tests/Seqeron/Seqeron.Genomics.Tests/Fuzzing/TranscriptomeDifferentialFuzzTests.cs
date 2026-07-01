using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Transcriptome differential-expression unit — log2 fold change
/// (<see cref="TranscriptomeAnalyzer.CalculateFoldChange"/>) plus the full per-gene
/// pipeline (<see cref="TranscriptomeAnalyzer.FindDifferentiallyExpressed"/>): Welch
/// unequal-variance t-test + Benjamini-Hochberg FDR + two-criterion significance gate.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// NaN/Infinity leaking into a fold change or p-value, and no *unhandled* runtime
/// exception (IndexOutOfRange, NullReference, DivideByZero, Overflow, …). Every
/// input must yield EITHER a well-defined, theory-correct result, OR a *documented,
/// intentional* convention (empty result, log2FC = 0, p = 1). A raw runtime
/// exception, a hang, a p-value outside [0, 1], a NaN/Infinity, or a *fabricated*
/// significant call on degenerate input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: TRANS-DIFF-001 — differential expression analysis (Transcriptome)
/// Checklist: docs/checklists/03_FUZZING.md, row 198.
/// Fuzz strategy for THIS unit: BE = Boundary Exploitation (0, -1, MaxInt, empty)
///   — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Fuzz targets (checklist row 198): "A=B, zero counts, single replicate".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (docs/algorithms/Transcriptome/Differential_Expression.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • log2FC = log2((mean2 + c) / (mean1 + c)) with pseudocount c = 1; positive =
///     up in condition 2 (treatment / numerator). — §2.2, §4.2.
///   • Raw p-value = exact two-sided Welch t tail
///       t = (m2 − m1)/√(s1²/N1 + s2²/N2),  ν = Welch-Satterthwaite df,
///       p = I_{ν/(ν+t²)}(ν/2, ½)  (regularized incomplete beta). — §2.2.
///   • Adjusted p = Benjamini-Hochberg step-up (R p.adjust BH). — §2.2, §4.2.
///   • Significant ⇔ |log2FC| ≥ threshold AND adjusted p &lt; α. — §2.2, §3.2.
///
/// Documented invariants this fixture pins (§2.4):
///   • INV-01: CalculateFoldChange(a,b) = −CalculateFoldChange(b,a).
///   • INV-02: equal group means ⇒ log2FC = 0 AND p = 1 ⇒ not significant.
///   • INV-03: adjusted p ≥ raw p and ≤ 1.
///   • INV-05: significant ⇔ |log2FC| ≥ threshold AND adjusted p &lt; α.
///
/// Boundary handling fixed by the doc (§3.3, §6.1, §5.4) and pinned here so the
/// contract can never silently drift:
///   • A = B (BE): identical condition vectors ⇒ ratio = 1 ⇒ log2FC = 0, t = 0 ⇒
///     p = 1, never significant. — INV-02, §6.1 ("Identical groups").
///   • ZERO COUNTS (BE): an all-zero gene ⇒ mean1 = mean2 = 0 ⇒ finite log2FC via
///     the pseudocount = log2(1/1) = 0, p = 1, never significant; no DivideByZero,
///     no NaN. — §3.3, §6.1 ("Zero mean expression"), §5.4 Deviation 1.
///   • SINGLE REPLICATE (BE): a group with &lt; 2 replicates has no unbiased variance,
///     so its raw p-value is set to 1.0 (not testable) ⇒ never significant, even
///     when |log2FC| is large. — §3.3, §6.1 ("Group with &lt;2 replicates"),
///     §5.4 Deviation 2, source MinReplicatesForTTest = 2.
///   • EMPTY / NULL gene enumerable ⇒ empty result. — §3.3, §6.1.
///   • NULL / empty expression list in CalculateFoldChange ⇒ fold change 0. — §3.1, §3.3.
///
/// Positive sanity (worked example, derived INDEPENDENTLY from the formulae and
/// cross-checked against SciPy ttest_ind(equal_var=False) / false_discovery_control,
/// NOT echoed off the implementation):
///   • §7.1: control = {1,2,3}, treatment = {7,8,9}
///       means 2 and 8; s1² = s2² = 1; se = √(1/3+1/3) = 0.8164966;
///       t = 6/0.8164966 = 7.3484692; Welch ν = 4;
///       p = I_{4/(4+t²)}(2, ½) = 0.0018262607  (SciPy-confirmed);
///       log2FC = log2((8+1)/(2+1)) = log2(3) = 1.5849625.
///   • §7.1 BH: raw p-values (0.001, 0.4, 0.5, 0.9) ⇒ adjusted (0.004, 0.66667,
///       0.66667, 0.9)  (SciPy false_discovery_control-confirmed).
/// A genuine large-effect, low-p gene must therefore be flagged Upregulated &amp;
/// significant, so a passing "no crash" result cannot be a degenerate analyzer that
/// returns log2FC = 0 / p = 1 / not-significant for everything.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Determinism
/// ───────────────────────────────────────────────────────────────────────────
/// All inputs are hand-built or generated from a LOCALLY fixed-seed
/// `new Random(seed)` (never a shared static Rng), so every run is reproducible.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class TranscriptomeDifferentialFuzzTests
{
    // Independently-derived worked-example constants (§7.1), NOT echoed off the code.
    //   log2FC = log2((8+1)/(2+1)) = log2(3).
    private const double WorkedLog2Fc = 1.5849625007211562; // Math.Log2(3)
    //   p = I_{4/(4+t²)}(2, ½) = 0.0018262607 (SciPy ttest_ind equal_var=False).
    private const double WorkedRawP = 0.0018262606682599833;

    private static (string, IReadOnlyList<double>, IReadOnlyList<double>)
        Gene(string id, IReadOnlyList<double> c1, IReadOnlyList<double> c2) => (id, c1, c2);

    #region TRANS-DIFF-001 — differential expression (log2FC + Welch t-test + BH)

    // ════════════════════════════════════════════════════════════════════════
    //  Positive sanity — the §7.1 worked example must be reproduced EXACTLY.
    //  Guards against a degenerate analyzer (constant log2FC = 0 / p = 1 /
    //  not-significant) that would pass every boundary test below.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void CalculateFoldChange_WorkedExample_EqualsLog2Of3()
    {
        // control {1,2,3} mean 2, treatment {7,8,9} mean 8 ⇒ log2((8+1)/(2+1)) = log2 3.
        double fc = TranscriptomeAnalyzer.CalculateFoldChange(
            new double[] { 1, 2, 3 }, new double[] { 7, 8, 9 });

        fc.Should().BeApproximately(WorkedLog2Fc, 1e-12,
            "log2((mean2+1)/(mean1+1)) = log2(9/3) = log2 3 — §7.1, derived from the formula");
    }

    [Test]
    public void FindDifferentiallyExpressed_WorkedExample_ReproducesLog2FcRawPAndSignificance()
    {
        var genes = new[] { Gene("G1", new double[] { 1, 2, 3 }, new double[] { 7, 8, 9 }) };

        var r = TranscriptomeAnalyzer.FindDifferentiallyExpressed(genes, alpha: 0.05, log2FoldChangeThreshold: 1.0)
            .Should().ContainSingle().Subject;

        r.Log2FoldChange.Should().BeApproximately(WorkedLog2Fc, 1e-12, "§7.1 log2 3");
        r.PValue.Should().BeApproximately(WorkedRawP, 1e-10,
            "exact two-sided Welch t tail I_{4/(4+t²)}(2,½) = 0.0018262607 — §7.1 (SciPy-confirmed)");
        // Single gene ⇒ BH m/rank = 1 ⇒ adjusted p == raw p.
        r.AdjustedPValue.Should().BeApproximately(WorkedRawP, 1e-10,
            "one gene ⇒ BH factor m/rank = 1 ⇒ adjusted = raw (INV-03 lower bound met)");
        r.IsSignificant.Should().BeTrue(
            "|log2 3| = 1.585 ≥ 1 AND adjusted p 0.0018 < 0.05 — two-criterion gate (INV-05)");
        r.Regulation.Should().Be("Upregulated", "log2FC > 0 ⇒ up in condition 2 (treatment)");
    }

    // INV-01: antisymmetry of the fold change under swapping the two conditions.
    [Test]
    public void CalculateFoldChange_IsAntisymmetricUnderConditionSwap()
    {
        var a = new double[] { 1, 2, 3 };
        var b = new double[] { 7, 8, 9 };

        double ab = TranscriptomeAnalyzer.CalculateFoldChange(a, b);
        double ba = TranscriptomeAnalyzer.CalculateFoldChange(b, a);

        ab.Should().BeApproximately(-ba, 1e-12, "log2((m2+c)/(m1+c)) = −log2((m1+c)/(m2+c)) — INV-01");
    }

    // §7.1 BH walk-through: raw (0.001,0.4,0.5,0.9) ⇒ adjusted (0.004,0.66667,0.66667,0.9).
    // Built so each gene produces (approximately) one of those raw p-values is hard
    // analytically; instead pin BH through the public pipeline using identical groups
    // for the high-p genes and a strong gene for the low-p, then assert INV-03/INV-04.
    [Test]
    public void FindDifferentiallyExpressed_AdjustedPValuesObeyBenjaminiHochbergInvariants()
    {
        // One strongly DE gene + three null genes (A=B ⇒ raw p = 1).
        var genes = new[]
        {
            Gene("strong", new double[] { 1, 2, 3 }, new double[] { 7, 8, 9 }), // raw p ≈ 0.00183
            Gene("null1",  new double[] { 5, 6, 7 }, new double[] { 5, 6, 7 }),  // raw p = 1
            Gene("null2",  new double[] { 4, 4, 4 }, new double[] { 4, 4, 4 }),  // raw p = 1
            Gene("null3",  new double[] { 9, 9, 9 }, new double[] { 9, 9, 9 }),  // raw p = 1
        };

        var results = TranscriptomeAnalyzer.FindDifferentiallyExpressed(genes, alpha: 0.05).ToList();

        foreach (var r in results)
        {
            r.AdjustedPValue.Should().BeGreaterThanOrEqualTo(r.PValue - 1e-12,
                "BH multiplies by m/rank ≥ 1 ⇒ adjusted ≥ raw — INV-03");
            r.AdjustedPValue.Should().BeLessThanOrEqualTo(1.0 + 1e-12, "BH is clamped to 1 — INV-03");
        }

        // Strongest gene: raw 0.00183, m = 4, smallest rank 1 ⇒ adjusted = 4·0.00183 = 0.0073.
        var strong = results.Single(r => r.GeneId == "strong");
        strong.AdjustedPValue.Should().BeApproximately(WorkedRawP * 4.0, 1e-9,
            "BH adjusted = m·p(1)/1 = 4·0.00183 — Benjamini-Hochberg step-up (§2.2)");
        strong.IsSignificant.Should().BeTrue("adjusted 0.0073 < 0.05 AND |log2 3| ≥ 1");
    }

    #endregion

    #region TRANS-DIFF-001 — BE boundary: A = B (identical conditions)

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: A = B (BE). Identical condition vectors ⇒ ratio = 1 ⇒
    // log2FC = 0, t = 0 ⇒ p = 1 ⇒ never significant; regulation "Unchanged".
    // — INV-02, §6.1 ("Identical groups").
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateFoldChange_IdenticalConditions_IsExactlyZero()
    {
        var v = new double[] { 3, 5, 8, 13 };
        TranscriptomeAnalyzer.CalculateFoldChange(v, v)
            .Should().Be(0.0, "equal means ⇒ ratio = 1 ⇒ log2 1 = 0 (INV-02)");
    }

    [Test]
    public void FindDifferentiallyExpressed_IdenticalConditions_NotSignificantUnchanged()
    {
        var genes = new[] { Gene("AeqB", new double[] { 2, 4, 6 }, new double[] { 2, 4, 6 }) };

        var r = TranscriptomeAnalyzer.FindDifferentiallyExpressed(genes).Should().ContainSingle().Subject;

        r.Log2FoldChange.Should().Be(0.0, "A = B ⇒ log2FC = 0 (INV-02)");
        r.PValue.Should().Be(1.0, "A = B ⇒ t = 0 ⇒ p = I_x(ν/2,½) at x = 1 = 1 (INV-02)");
        r.AdjustedPValue.Should().Be(1.0, "raw 1 ⇒ adjusted 1 (clamped, INV-03)");
        r.IsSignificant.Should().BeFalse("no differential expression when A = B (INV-02/INV-05)");
        r.Regulation.Should().Be("Unchanged", "log2FC = 0 ⇒ neither up nor down");
    }

    #endregion

    #region TRANS-DIFF-001 — BE boundary: ZERO COUNTS (all-zero genes)

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: ZERO COUNTS (BE). An all-zero gene ⇒ mean1 = mean2 = 0 ⇒
    // finite log2FC via the pseudocount = log2((0+1)/(0+1)) = 0; t = 0 ⇒ p = 1;
    // never significant. No DivideByZero, no NaN. — §3.3, §6.1, §5.4 Deviation 1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateFoldChange_AllZeroCounts_IsFiniteZeroViaPseudocount()
    {
        var zeros = new double[] { 0, 0, 0 };
        double fc = TranscriptomeAnalyzer.CalculateFoldChange(zeros, zeros);

        fc.Should().Be(0.0, "log2((0+1)/(0+1)) = log2 1 = 0 — pseudocount keeps it finite (§5.4)");
        double.IsNaN(fc).Should().BeFalse("pseudocount prevents 0/0 ⇒ no NaN");
        double.IsInfinity(fc).Should().BeFalse("pseudocount prevents log2(0) ⇒ no Infinity");
    }

    // One side zero, other non-zero: still finite via pseudocount, sign correct.
    [Test]
    public void CalculateFoldChange_ZeroControlNonZeroTreatment_FiniteAndPositive()
    {
        // mean1 = 0, mean2 = 7 ⇒ log2((7+1)/(0+1)) = log2 8 = 3.
        double fc = TranscriptomeAnalyzer.CalculateFoldChange(
            new double[] { 0, 0, 0 }, new double[] { 7, 7, 7 });

        fc.Should().BeApproximately(3.0, 1e-12, "log2((7+1)/(0+1)) = log2 8 = 3 — regularized ratio (§6.1)");
        double.IsNaN(fc).Should().BeFalse();
        double.IsInfinity(fc).Should().BeFalse();
    }

    [Test]
    public void FindDifferentiallyExpressed_AllZeroGene_NotSignificantNoNaN()
    {
        var genes = new[] { Gene("zero", new double[] { 0, 0, 0 }, new double[] { 0, 0, 0 }) };

        var r = TranscriptomeAnalyzer.FindDifferentiallyExpressed(genes).Should().ContainSingle().Subject;

        r.Log2FoldChange.Should().Be(0.0, "all-zero ⇒ log2FC = 0 (pseudocount)");
        double.IsNaN(r.PValue).Should().BeFalse("zero-variance identical means ⇒ p = 1, not NaN");
        r.PValue.Should().Be(1.0, "se = 0 with equal means ⇒ p = 1 (§5.4 Deviation 3)");
        r.IsSignificant.Should().BeFalse("an all-zero gene is never differentially expressed");
    }

    #endregion

    #region TRANS-DIFF-001 — BE boundary: SINGLE REPLICATE per condition

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: SINGLE REPLICATE (BE). A group with < 2 replicates has no
    // unbiased variance, so its raw p-value is set to 1.0 (not testable) ⇒ never
    // flagged significant, EVEN when |log2FC| is large. — §3.3, §6.1, §5.4
    // Deviation 2, source MinReplicatesForTTest = 2.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindDifferentiallyExpressed_SingleReplicatePerCondition_PIsOneNeverSignificant()
    {
        // Huge effect (mean1 = 5, mean2 = 100 ⇒ log2((100+1)/(5+1)) ≈ 4.07) but N = 1
        // per group ⇒ not testable ⇒ p = 1 ⇒ not significant.
        var genes = new[] { Gene("single", new double[] { 5 }, new double[] { 100 }) };

        var r = TranscriptomeAnalyzer.FindDifferentiallyExpressed(genes).Should().ContainSingle().Subject;

        r.Log2FoldChange.Should().BeApproximately(Math.Log2(101.0 / 6.0), 1e-12,
            "fold change is still computed from the single replicate means (§3.2)");
        r.Log2FoldChange.Should().BeGreaterThan(1.0, "the effect is large, threshold = 1");
        r.PValue.Should().Be(1.0, "N < 2 ⇒ variance undefined ⇒ p = 1 (Welch precondition, §5.4 Deviation 2)");
        r.AdjustedPValue.Should().Be(1.0, "raw 1 ⇒ adjusted 1");
        r.IsSignificant.Should().BeFalse(
            "single-replicate genes are NEVER significant despite a large |log2FC| (INV-05 fails on adj p)");
    }

    // Single replicate on only ONE side is also not testable ⇒ p = 1.
    [Test]
    public void FindDifferentiallyExpressed_OneSideSingleReplicate_PIsOne()
    {
        var genes = new[] { Gene("mixed", new double[] { 1 }, new double[] { 7, 8, 9 }) };

        var r = TranscriptomeAnalyzer.FindDifferentiallyExpressed(genes).Should().ContainSingle().Subject;

        r.PValue.Should().Be(1.0, "either group with N < 2 ⇒ not testable ⇒ p = 1 (§6.1)");
        r.IsSignificant.Should().BeFalse();
    }

    #endregion

    #region TRANS-DIFF-001 — BE boundary: EMPTY / NULL inputs

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY / NULL gene enumerable (BE) ⇒ empty result, no crash.
    // — §3.3, §6.1 ("Empty gene enumerable / null").
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindDifferentiallyExpressed_EmptyOrNullGenes_YieldEmptyResult()
    {
        TranscriptomeAnalyzer.FindDifferentiallyExpressed(
            Array.Empty<(string, IReadOnlyList<double>, IReadOnlyList<double>)>())
            .Should().BeEmpty("no genes to test ⇒ empty result (§6.1)");

        TranscriptomeAnalyzer.FindDifferentiallyExpressed(null!)
            .Should().BeEmpty("null gene enumerable ⇒ empty result (§3.3)");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Null / empty expression list in CalculateFoldChange ⇒ fold change 0,
    // no NullReferenceException, no crash. — §3.1, §3.3.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateFoldChange_NullOrEmptyExpression_ReturnsZero()
    {
        TranscriptomeAnalyzer.CalculateFoldChange(null!, new double[] { 1, 2 })
            .Should().Be(0.0, "null condition ⇒ 0 (§3.1)");
        TranscriptomeAnalyzer.CalculateFoldChange(new double[] { 1, 2 }, null!)
            .Should().Be(0.0, "null condition ⇒ 0 (§3.1)");
        TranscriptomeAnalyzer.CalculateFoldChange(Array.Empty<double>(), new double[] { 1, 2 })
            .Should().Be(0.0, "empty condition ⇒ 0 (§3.3)");
        TranscriptomeAnalyzer.CalculateFoldChange(new double[] { 1, 2 }, Array.Empty<double>())
            .Should().Be(0.0, "empty condition ⇒ 0 (§3.3)");
    }

    #endregion

    #region TRANS-DIFF-001 — randomized boundary sweep

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: random boundary batch (BE) under a time budget.
    // A deterministic, locally-seeded generator builds gene panels spanning the
    // boundaries (A = B, all-zero, single replicate, ±1/MaxInt-scale magnitudes,
    // empty groups). FindDifferentiallyExpressed must process every case without
    // crashing or hanging, and EVERY result must be well-formed:
    //   • Log2FoldChange finite (no NaN/Infinity);
    //   • raw & adjusted p ∈ [0, 1], finite;
    //   • adjusted p ≥ raw p − ε (INV-03);
    //   • A = B / all-zero / single-replicate ⇒ NOT significant;
    //   • Regulation matches the sign of Log2FoldChange (INV-05 labelling).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void FindDifferentiallyExpressed_RandomBoundaryBatch_NeverCrashesAndStaysWellFormed()
    {
        var rng = new Random(20260620); // locally fixed seed — deterministic

        for (int trial = 0; trial < 300; trial++)
        {
            int geneCount = rng.Next(0, 8);
            var genes = new List<(string, IReadOnlyList<double>, IReadOnlyList<double>)>(geneCount);
            // Track which genes are forced degenerate (must end up not-significant).
            var degenerate = new HashSet<string>();

            for (int g = 0; g < geneCount; g++)
            {
                string id = $"g{g}";
                IReadOnlyList<double> c1, c2;

                switch (rng.Next(6))
                {
                    case 0: // A = B (identical) — must be not significant
                        c1 = RandomVector(rng, rng.Next(1, 5));
                        c2 = c1.ToArray();
                        degenerate.Add(id);
                        break;
                    case 1: // all-zero — must be not significant
                        c1 = new double[rng.Next(1, 5)];
                        c2 = new double[rng.Next(1, 5)];
                        degenerate.Add(id);
                        break;
                    case 2: // single replicate per side — not testable ⇒ not significant
                        c1 = new[] { rng.NextDouble() * 1e6 };
                        c2 = new[] { rng.NextDouble() * 1e6 };
                        degenerate.Add(id);
                        break;
                    case 3: // extreme magnitudes near int.MaxValue scale (BE)
                        c1 = RandomVector(rng, 3, scale: int.MaxValue);
                        c2 = RandomVector(rng, 3, scale: int.MaxValue);
                        break;
                    case 4: // one side empty ⇒ fold change 0, p = 1 (N < 2)
                        c1 = Array.Empty<double>();
                        c2 = RandomVector(rng, rng.Next(2, 5));
                        degenerate.Add(id);
                        break;
                    default: // ordinary replicated gene
                        c1 = RandomVector(rng, rng.Next(2, 6));
                        c2 = RandomVector(rng, rng.Next(2, 6));
                        break;
                }

                genes.Add((id, c1, c2));
            }

            var results = TranscriptomeAnalyzer
                .FindDifferentiallyExpressed(genes, alpha: 0.05, log2FoldChangeThreshold: 1.0)
                .ToList();

            results.Should().HaveCount(geneCount, "exactly one result per gene");

            foreach (var r in results)
            {
                double.IsNaN(r.Log2FoldChange).Should().BeFalse("log2FC never NaN on boundary input");
                double.IsInfinity(r.Log2FoldChange).Should().BeFalse("log2FC never Infinity (pseudocount)");

                r.PValue.Should().BeInRange(0.0, 1.0, "raw p is a probability (INV)");
                r.AdjustedPValue.Should().BeInRange(0.0, 1.0 + 1e-12, "adjusted p is clamped to 1 (INV-03)");
                double.IsNaN(r.PValue).Should().BeFalse("no NaN p-value");
                double.IsNaN(r.AdjustedPValue).Should().BeFalse("no NaN adjusted p-value");
                r.AdjustedPValue.Should().BeGreaterThanOrEqualTo(r.PValue - 1e-9,
                    "BH multiplies by m/rank ≥ 1 ⇒ adjusted ≥ raw (INV-03)");

                string expectedReg = r.Log2FoldChange > 0 ? "Upregulated"
                    : r.Log2FoldChange < 0 ? "Downregulated" : "Unchanged";
                r.Regulation.Should().Be(expectedReg, "regulation labels the sign of log2FC");

                if (degenerate.Contains(r.GeneId))
                    r.IsSignificant.Should().BeFalse(
                        "A = B / all-zero / single-replicate / empty genes are never significant (INV-02/§6.1)");
            }
        }
    }

    private static double[] RandomVector(Random rng, int n, double scale = 100.0) =>
        Enumerable.Range(0, n).Select(_ => rng.NextDouble() * scale).ToArray();

    #endregion
}
