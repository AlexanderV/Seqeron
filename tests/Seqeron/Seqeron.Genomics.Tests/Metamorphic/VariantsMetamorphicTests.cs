using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Variants area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: VARIANT-ANNOT-001 — variant consequence annotation (Variants).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 186.
///
/// API under test (VariantAnnotator.AnnotateVariant / AnnotateVariants):
///   Predicts a variant's consequence on transcripts from the relative position of the variant
///   within the transcript's exon/CDS geometry.
///
/// Relations (derived from the position-relative consequence model, NOT from output):
///   • SHIFT (coordinate shift shifts annotations): shifting the variant AND the transcript by the
///          same offset preserves the relative geometry, so the consequence is unchanged.
///   • INV  (variant order independent): each variant is annotated against the transcripts
///          independently, so the per-variant consequences do not depend on input order.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class VariantsMetamorphicTests
{
    private static VariantAnnotator.Transcript Transcript(int offset = 0) =>
        new("T1", "G1", "GENE", "chr1",
            100 + offset, 200 + offset, '+',
            new[] { (100 + offset, 200 + offset) },
            new[] { (100 + offset, 200 + offset) },
            100 + offset, 200 + offset);

    private static VariantAnnotator.Variant Snv(int position) =>
        new("chr1", position, "A", "G", VariantAnnotator.ClassifyVariant("A", "G"));

    #region VARIANT-ANNOT-001 SHIFT — a uniform coordinate shift preserves the consequence

    [Test]
    [Description("SHIFT: shifting both the variant and the transcript by the same offset preserves their relative geometry, so the predicted consequence is unchanged.")]
    public void Annotation_UniformCoordinateShift_PreservesConsequence()
    {
        var baseConsequence = VariantAnnotator.AnnotateVariant(Snv(150), new[] { Transcript() }).Single().Consequence;

        foreach (int offset in new[] { 1000, 50000 })
        {
            var shifted = VariantAnnotator.AnnotateVariant(Snv(150 + offset), new[] { Transcript(offset) }).Single().Consequence;
            shifted.Should().Be(baseConsequence,
                because: $"shifting variant and transcript together by {offset} keeps the variant at the same position within the transcript");
        }
    }

    #endregion

    #region VARIANT-ANNOT-001 INV — annotation is independent of variant order

    [Test]
    [Description("INV: each variant is annotated against the transcripts independently, so reordering the input variants yields the same per-variant consequences.")]
    public void Annotation_VariantOrder_Invariant()
    {
        var transcripts = new[] { Transcript() };
        var variants = new[] { Snv(120), Snv(180), Snv(9000) };

        Dictionary<int, List<VariantAnnotator.ConsequenceType>> Consequences(IEnumerable<VariantAnnotator.Variant> vs) =>
            VariantAnnotator.AnnotateVariants(vs, transcripts)
                .ToDictionary(g => g.Key.Position, g => g.Select(a => a.Consequence).OrderBy(c => c).ToList());

        var forward = Consequences(variants);
        var reversed = Consequences(variants.Reverse());

        foreach (var pos in forward.Keys)
            reversed[pos].Should().Equal(forward[pos],
                because: $"the consequences of the variant at {pos} do not depend on the order of the variant list");
    }

    #endregion
}
