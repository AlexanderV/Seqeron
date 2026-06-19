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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-COMPLEX-DUST-001 — DUST low-complexity score (Complexity).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 229.
    //
    // API under test (SequenceComplexity.CalculateDustScore):
    //   Σ_t c_t·(c_t−1)/2 over overlapping triplets, divided by the word count; HIGHER ⇒ LOWER
    //   complexity (Morgulis et al. 2006).
    //
    // Relations (derived from the triplet-count sum, NOT from output):
    //   • INV  (complement preserves DUST): complement maps each base A↔T/C↔G, a bijection on the
    //          triplet alphabet, so the multiset of triplet counts — and hence the score — is unchanged.
    //   • MONO (adding a homopolymer run raises score): appending a homopolymer accumulates one triplet
    //          quadratically while the length grows linearly, so the DUST score increases.
    // ───────────────────────────────────────────────────────────────────────────

    private static string Complement(string seq) =>
        new string(seq.Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c }).ToArray());

    #region SEQ-COMPLEX-DUST-001 INV — complement preserves the DUST score

    [Test]
    [Description("INV: complement is a bijection on the triplet alphabet, so it preserves the multiset of triplet counts and hence the DUST score.")]
    public void Dust_Complement_PreservesScore()
    {
        double original = SequenceComplexity.CalculateDustScore(ComplexSeq);
        original.Should().BeGreaterThan(0, because: "the test sequence has at least one repeated triplet — a non-vacuous fixture");

        SequenceComplexity.CalculateDustScore(Complement(ComplexSeq))
            .Should().BeApproximately(original, 1e-12,
                because: "A↔T / C↔G relabels every triplet bijectively, leaving the count distribution unchanged");
    }

    #endregion

    #region SEQ-COMPLEX-DUST-001 MONO — appending a homopolymer run raises the score

    [Test]
    [Description("MONO: appending a longer homopolymer run accumulates one triplet quadratically against a linear length, so the DUST score increases.")]
    public void Dust_AddingHomopolymerRun_RaisesScore()
    {
        double previous = double.MinValue;
        foreach (int runLength in new[] { 0, 10, 20, 40 })
        {
            string seq = ComplexSeq + new string('A', runLength);
            double score = SequenceComplexity.CalculateDustScore(seq);

            score.Should().BeGreaterThan(previous,
                because: $"a longer poly-A run adds more identical 'AAA' triplets, lowering complexity (raising DUST) — run length {runLength}");
            previous = score;
        }
    }

    #endregion
}
