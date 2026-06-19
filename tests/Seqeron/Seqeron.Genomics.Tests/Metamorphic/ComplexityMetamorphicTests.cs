using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Complexity area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-COMPLEX-COMPRESS-001 — compression-based complexity (Complexity).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 228.
///
/// API under test (SequenceComplexity.EstimateCompressionRatio):
///   Normalized Lempel–Ziv complexity c / (n / log_b n); lower ⇒ more repetitive/less complex.
///
/// Relations (derived from the LZ parse, NOT from output):
///   • INV   (case change preserves ratio): counting is case-folded, so upper/lower/mixed case give
///           the same complexity.
///   • ORDER (concatenating repeats lowers the ratio): a tandem repetition of a sequence is highly
///           compressible — the LZ parse reuses earlier components — so its normalized complexity is
///           strictly below the single copy's and does not increase as more copies are appended.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class ComplexityMetamorphicTests
{
    // A moderately complex DNA sequence (no obvious short period), so the single copy is far from the
    // repetitive floor and repetition can visibly lower the complexity.
    private const string ComplexSeq = "ACGTGACTTGACATGCGTAACGTTAGC";

    #region SEQ-COMPLEX-COMPRESS-001 INV — case folding preserves the ratio

    [Test]
    [Description("INV: complexity counting is case-folded, so upper, lower and mixed case give identical normalized LZ complexity.")]
    public void Compression_CaseChange_PreservesRatio()
    {
        double upper = SequenceComplexity.EstimateCompressionRatio(ComplexSeq);

        SequenceComplexity.EstimateCompressionRatio(ComplexSeq.ToLowerInvariant())
            .Should().Be(upper, because: "the sequence is upper-cased before the LZ parse");

        string mixed = new string(ComplexSeq.Select((c, i) => i % 2 == 0 ? char.ToLowerInvariant(c) : c).ToArray());
        SequenceComplexity.EstimateCompressionRatio(mixed)
            .Should().Be(upper, because: "case folding makes the LZ parse independent of letter case");
    }

    #endregion

    #region SEQ-COMPLEX-COMPRESS-001 ORDER — tandem repetition lowers the ratio

    [Test]
    [Description("ORDER: as a sequence is repeated more times it becomes increasingly compressible, so its normalized LZ complexity falls — strictly decreasing over well-separated repetition counts (the asymptotic trend; adjacent counts can wobble through the n/log_b n normalisation).")]
    public void Compression_ConcatenatingRepeats_LowersRatio()
    {
        double previous = double.MaxValue;
        foreach (int copies in new[] { 1, 8, 64 })
        {
            string repeated = string.Concat(Enumerable.Repeat(ComplexSeq, copies));
            double ratio = SequenceComplexity.EstimateCompressionRatio(repeated);

            ratio.Should().BeLessThan(previous,
                because: $"a {copies}× tandem repeat reuses earlier LZ components, lowering normalized complexity relative to fewer copies");
            previous = ratio;
        }
    }

    #endregion
}
