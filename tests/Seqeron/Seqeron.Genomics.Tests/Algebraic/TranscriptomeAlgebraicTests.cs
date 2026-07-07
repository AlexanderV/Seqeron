using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Transcriptome area (TPM normalization).
///
/// Algebraic testing pins the sequencing-depth scale-invariance of TPM (a
/// within-sample relative measure) and its empty-input identity.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, row 199.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Transcriptome")]
public class TranscriptomeAlgebraicTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: TRANS-EXPR-001 — TPM expression (Transcriptome), row 199.
    //
    // Model: TPM_i = (X_i/l_i) / Σ_j(X_j/l_j) × 10^6; because the per-gene rate is
    //        normalized by the sum of rates, multiplying every raw count by a common
    //        factor (deeper sequencing) leaves every TPM unchanged.
    //   — docs/algorithms/Transcriptome; TranscriptomeAnalyzer.CalculateTPM.
    //
    // Laws (row 199): HOMO — scaling all counts by k → identical TPM.
    //                 ID — empty input → no expression rows.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Tpm_Identity_EmptyYieldsNoRows()
    {
        TranscriptomeAnalyzer.CalculateTPM(
            System.Array.Empty<(string, double, int)>()).Should().BeEmpty();
    }

    [FsCheck.NUnit.Property]
    public Property Tpm_Homogeneous_DepthScaleInvariant()
    {
        var gen = (from n in Gen.Choose(1, 6)
                   from counts in Gen.Choose(1, 1000).ArrayOf(n)
                   from lengths in Gen.Choose(100, 5000).ArrayOf(n)
                   from k in Gen.Choose(2, 50)
                   select (counts, lengths, k)).ToArbitrary();
        return Prop.ForAll(gen, t =>
        {
            var baseGenes = t.counts.Zip(t.lengths, (c, l) => (id: $"g{c}_{l}", (double)c, l)).ToList();
            // Distinct gene ids to keep them as separate rows.
            var genes = baseGenes.Select((g, i) => ($"g{i}", g.Item2, g.l)).ToList();
            var scaled = genes.Select(g => (g.Item1, g.Item2 * t.k, g.l)).ToList();

            var baseTpm = TranscriptomeAnalyzer.CalculateTPM(genes).OrderBy(e => e.GeneId)
                .Select(e => e.TPM).ToList();
            var scaledTpm = TranscriptomeAnalyzer.CalculateTPM(scaled).OrderBy(e => e.GeneId)
                .Select(e => e.TPM).ToList();

            bool ok = baseTpm.Count == scaledTpm.Count
                && baseTpm.Zip(scaledTpm, (a, b) => System.Math.Abs(a - b) < 1e-6).All(x => x);
            return ok.Label($"TPM changed under depth scaling by {t.k}");
        });
    }
}
