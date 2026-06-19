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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-COMPLEX-KMER-001 — k-mer entropy (Complexity).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 230.
    //
    // API under test (SequenceComplexity.CalculateKmerEntropy):
    //   Shannon entropy (bits) of the overlapping k-mer frequency distribution.
    //
    // Relations (derived from the entropy of the k-mer distribution, NOT from output):
    //   • INV  (reverse preserves k-mer entropy): reversing the sequence maps each k-mer to its
    //          reverse — a bijection — and preserves the position-for-position counts, so the count
    //          distribution (and its entropy) is unchanged.
    //   • MONO (more distinct k-mers ⇒ higher entropy): with near-uniform usage, increasing the number
    //          of distinct k-mers raises the entropy of the distribution.
    // ───────────────────────────────────────────────────────────────────────────

    private static int DistinctKmers(string seq, int k) =>
        Enumerable.Range(0, seq.Length - k + 1).Select(i => seq.Substring(i, k)).Distinct().Count();

    #region SEQ-COMPLEX-KMER-001 INV — reversing the sequence preserves k-mer entropy

    [Test]
    [Description("INV: reversal maps each k-mer to its reverse (a bijection) and preserves the per-position counts, so the k-mer entropy is unchanged.")]
    public void KmerEntropy_Reverse_PreservesEntropy()
    {
        foreach (int k in new[] { 2, 3 })
            foreach (var seq in new[] { ComplexSeq, "ACGTACGTTTGCA", "AACCGGTTACGT" })
            {
                string reversed = new string(seq.Reverse().ToArray());
                SequenceComplexity.CalculateKmerEntropy(reversed, k)
                    .Should().BeApproximately(SequenceComplexity.CalculateKmerEntropy(seq, k), 1e-12,
                        because: $"reversal is a bijection on {k}-mers and keeps the count distribution, so the entropy is invariant");
            }
    }

    #endregion

    #region SEQ-COMPLEX-KMER-001 MONO — more distinct k-mers gives higher entropy

    [Test]
    [Description("MONO: with near-uniform usage, increasing the number of distinct k-mers raises the k-mer entropy.")]
    public void KmerEntropy_MoreDistinctKmers_HigherEntropy()
    {
        const int k = 2;
        // 1, 2, 3, 4 distinct near-uniform dinucleotides, by construction.
        var sequences = new[]
        {
            "AAAAAAAAAAAA", // {AA}
            "ATATATATATAT", // {AT, TA}
            "ACGACGACGACG", // {AC, CG, GA}
            "ACGTACGTACGT", // {AC, CG, GT, TA}
        };

        int previousDistinct = 0;
        double previousEntropy = double.MinValue;
        foreach (var seq in sequences)
        {
            int distinct = DistinctKmers(seq, k);
            double entropy = SequenceComplexity.CalculateKmerEntropy(seq, k);

            distinct.Should().BeGreaterThan(previousDistinct, because: "the fixtures are ordered by increasing distinct dinucleotide count");
            entropy.Should().BeGreaterThan(previousEntropy, because: $"more distinct near-uniform {k}-mers spread the distribution, raising its entropy");
            previousDistinct = distinct;
            previousEntropy = entropy;
        }
    }

    #endregion
}
