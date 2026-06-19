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
}
