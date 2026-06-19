using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Transcriptome area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: TRANS-DIFF-001 — differential expression / fold change (Transcriptome).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 198.
///
/// API under test (TranscriptomeAnalyzer.CalculateFoldChange / FindDifferentiallyExpressed):
///   log2 fold change between conditions; Welch t-test with Benjamini-Hochberg FDR.
///
/// Relations (derived from the log-ratio and rank-based FDR, NOT from output):
///   • SYM  (FC(A,B) = −FC(B,A)): the log2 ratio swaps sign when the conditions are swapped.
///   • INV  (gene order independent): BH adjustment is rank-based, so per-gene results do not depend
///          on the order of the input genes.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class TranscriptomeMetamorphicTests
{
    #region TRANS-DIFF-001 SYM — fold change is antisymmetric

    [Test]
    [Description("SYM: log2 fold change is a log ratio, so swapping the two conditions negates it: FC(A,B) = −FC(B,A).")]
    public void FoldChange_SwapConditions_NegatesSign()
    {
        var pairs = new[]
        {
            (new double[] { 10, 12, 11 }, new double[] { 40, 44, 42 }),
            (new double[] { 5, 6, 7 }, new double[] { 5, 6, 7 }),
            (new double[] { 100, 90 }, new double[] { 10, 12 }),
        };

        foreach (var (a, b) in pairs)
            TranscriptomeAnalyzer.CalculateFoldChange(b, a)
                .Should().BeApproximately(-TranscriptomeAnalyzer.CalculateFoldChange(a, b), 1e-12,
                    because: "log2(meanA/meanB) = −log2(meanB/meanA)");
    }

    #endregion

    #region TRANS-DIFF-001 INV — differential expression is independent of gene order

    [Test]
    [Description("INV: Benjamini-Hochberg FDR is rank-based, so reordering the input genes leaves every gene's fold change, adjusted p-value and significance unchanged.")]
    public void DifferentialExpression_GeneOrder_Invariant()
    {
        var genes = new (string, IReadOnlyList<double>, IReadOnlyList<double>)[]
        {
            ("up",   new double[] { 10, 11, 9 },   new double[] { 80, 82, 78 }),
            ("down", new double[] { 90, 88, 92 },  new double[] { 9, 11, 10 }),
            ("flat", new double[] { 50, 51, 49 },  new double[] { 50, 50, 51 }),
            ("mild", new double[] { 20, 22, 19 },  new double[] { 26, 25, 27 }),
        };

        Dictionary<string, (double Fc, double Adj, bool Sig)> Run(
            IEnumerable<(string, IReadOnlyList<double>, IReadOnlyList<double>)> g) =>
            TranscriptomeAnalyzer.FindDifferentiallyExpressed(g)
                .ToDictionary(d => d.GeneId, d => (d.Log2FoldChange, d.AdjustedPValue, d.IsSignificant));

        var forward = Run(genes);
        var reversed = Run(genes.Reverse());

        foreach (var id in forward.Keys)
        {
            reversed[id].Fc.Should().BeApproximately(forward[id].Fc, 1e-12, because: $"{id}'s fold change is per-gene");
            reversed[id].Adj.Should().BeApproximately(forward[id].Adj, 1e-12, because: $"{id}'s BH-adjusted p-value is rank-based, not order-based");
            reversed[id].Sig.Should().Be(forward[id].Sig, because: $"{id}'s significance call is order-independent");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: TRANS-EXPR-001 — TPM normalization (Transcriptome).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 199.
    //
    // API under test (TranscriptomeAnalyzer.CalculateTPM):
    //   TPM_i = (X_i/l_i) / Σ_j(X_j/l_j) × 10^6.
    //
    // Relations (derived from the within-sample normalization, NOT from output):
    //   • HOMO (scaling depth preserves TPM): multiplying every raw count by a constant cancels in
    //          the ratio, so TPM is unchanged.
    //   • INV  (read/gene order independent): each gene's TPM depends on its rate and the order-
    //          independent sum of rates, so reordering the genes preserves every TPM.
    // ───────────────────────────────────────────────────────────────────────────

    private static readonly (string, double, int)[] TpmGenes =
    {
        ("g1", 100, 1000),
        ("g2", 50, 500),
        ("g3", 200, 2000),
    };

    private static Dictionary<string, double> Tpm(IEnumerable<(string, double, int)> genes) =>
        TranscriptomeAnalyzer.CalculateTPM(genes).ToDictionary(g => g.GeneId, g => g.TPM);

    #region TRANS-EXPR-001 HOMO — scaling sequencing depth preserves TPM

    [Test]
    [Description("HOMO: multiplying every raw count by a constant scales the numerator and the sum-of-rates equally, so the TPM of every gene is unchanged.")]
    public void Tpm_ScaleDepth_Preserved()
    {
        var baseline = Tpm(TpmGenes);

        foreach (double k in new[] { 2.0, 10.0 })
        {
            var scaled = Tpm(TpmGenes.Select(g => (g.Item1, g.Item2 * k, g.Item3)));
            foreach (var id in baseline.Keys)
                scaled[id].Should().BeApproximately(baseline[id], 1e-9,
                    because: $"scaling all counts by {k} cancels in the TPM ratio");
        }
    }

    #endregion

    #region TRANS-EXPR-001 INV — TPM is independent of gene order

    [Test]
    [Description("INV: each gene's TPM depends on its rate and the order-independent sum of rates, so reversing the gene order yields the same TPMs.")]
    public void Tpm_GeneOrder_Invariant()
    {
        var baseline = Tpm(TpmGenes);
        var reversed = Tpm(TpmGenes.Reverse());

        foreach (var id in baseline.Keys)
            reversed[id].Should().BeApproximately(baseline[id], 1e-12, because: "TPM does not depend on the order of the genes");
    }

    #endregion
}
