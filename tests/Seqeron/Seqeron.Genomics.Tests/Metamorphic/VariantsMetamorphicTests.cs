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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: VARIANT-CALL-001 — variant detection from reference↔query (Variants).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 187.
    //
    // API under test (VariantCaller.CallVariantsFromAlignment):
    //   Reports one variant per differing alignment column (SNP / indel).
    //
    // Relations (derived from "a variant is a difference from reference", NOT from output):
    //   • MON  (more differences ⇒ superset of calls): each additional mismatch is called, so adding
    //          substitutions yields a superset of variant calls (the analog of "deeper evidence →
    //          superset of confident calls" for this difference-based caller).
    //   • INV  (identity ⇒ no variants; deterministic): identical sequences produce no variants, and
    //          the call set is a pure function of the inputs.
    // ───────────────────────────────────────────────────────────────────────────

    private const string CallReference = "ACGTACGTAC";

    private static string SubstituteFirst(string seq, int count)
    {
        char[] arr = seq.ToCharArray();
        for (int i = 0; i < count; i++) arr[i] = arr[i] == 'A' ? 'C' : 'A';
        return new string(arr);
    }

    #region VARIANT-CALL-001 MON — more differences yield a superset of calls

    [Test]
    [Description("MON: each additional substitution is called as a variant, so introducing more substitutions yields a superset of variant positions.")]
    public void VariantCalls_MoreDifferences_Superset()
    {
        int previous = -1;
        var prevPositions = new HashSet<int>();
        foreach (int subs in new[] { 1, 2, 4 })
        {
            var positions = VariantCaller.CallVariantsFromAlignment(CallReference, SubstituteFirst(CallReference, subs))
                .Select(v => v.Position).ToHashSet();
            if (previous >= 0)
                prevPositions.IsSubsetOf(positions).Should().BeTrue(because: $"the {previous}-substitution calls are a subset of the {subs}-substitution calls");
            positions.Count.Should().BeGreaterThan(previous, because: $"{subs} nested substitutions are all called");
            previous = subs;
            prevPositions = positions;
        }
    }

    #endregion

    #region VARIANT-CALL-001 INV — identity gives no variants and is deterministic

    [Test]
    [Description("INV: identical sequences differ nowhere, so no variant is called; and the call is a pure function returning the same variants on repeat.")]
    public void VariantCalls_Identity_NoVariants_Deterministic()
    {
        VariantCaller.CallVariantsFromAlignment(CallReference, CallReference).Should().BeEmpty(
            because: "a variant is a difference from the reference, and there is none");

        string query = SubstituteFirst(CallReference, 3);
        VariantCaller.CallVariantsFromAlignment(CallReference, query).Select(v => v.Position).ToList()
            .Should().Equal(VariantCaller.CallVariantsFromAlignment(CallReference, query).Select(v => v.Position).ToList(),
                because: "variant calling has no hidden state");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: VARIANT-INDEL-001 — indel detection (Variants).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 188.
    //
    // API under test (VariantCaller.FindIndels):
    //   Aligns reference and query and reports the insertion/deletion variants.
    //
    // Relations (derived from coordinate-relative calling, NOT from output):
    //   • SHIFT (prepend flank shifts indel positions): prepending an identical flank to both
    //          sequences shifts the indel positions by the flank length.
    //   • INV  (deterministic): the indel call set is a pure function of the inputs.
    // ───────────────────────────────────────────────────────────────────────────

    private const string IndelReference = "ACGTCATGCATAG";
    private const string IndelQuery = "ACGTCAGTGCATAG"; // inserted 'G' between the A and T (unambiguous)

    #region VARIANT-INDEL-001 SHIFT — a prepended flank shifts indel positions

    [Test]
    [Description("SHIFT: prepending an identical flank to both reference and query shifts every indel position by the flank length.")]
    public void Indels_PrependFlank_ShiftsPositions()
    {
        var original = VariantCaller.FindIndels(new DnaSequence(IndelReference), new DnaSequence(IndelQuery))
            .Select(v => v.Position).ToList();
        original.Should().NotBeEmpty(because: "the query has a one-base insertion relative to the reference");

        foreach (var flank in new[] { "TTTT", "GGCCAA" })
        {
            var shifted = VariantCaller.FindIndels(new DnaSequence(flank + IndelReference), new DnaSequence(flank + IndelQuery))
                .Select(v => v.Position).ToList();
            shifted.Should().Equal(original.Select(p => p + flank.Length),
                because: $"prepending an identical {flank.Length}-base flank shifts the indel by {flank.Length}");
        }
    }

    #endregion

    #region VARIANT-INDEL-001 INV — indel detection is deterministic

    [Test]
    [Description("INV: indel detection is a pure function, so repeated calls return the identical indel positions.")]
    public void Indels_Deterministic()
    {
        VariantCaller.FindIndels(new DnaSequence(IndelReference), new DnaSequence(IndelQuery)).Select(v => v.Position).ToList()
            .Should().Equal(VariantCaller.FindIndels(new DnaSequence(IndelReference), new DnaSequence(IndelQuery)).Select(v => v.Position).ToList(),
                because: "indel calling has no hidden state");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: VARIANT-SNP-001 — SNP detection (Variants).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 189.
    //
    // API under test (VariantCaller.FindSnpsDirect):
    //   Reports a SNP at every position where the aligned reference and query bases differ.
    //
    // Relations (derived from positional comparison, NOT from output):
    //   • SHIFT (prepend flank shifts SNP positions): prepending an identical flank to both shifts
    //          every SNP position by the flank length.
    //   • INV  (deterministic): the SNP set is a pure function of the inputs.
    // ───────────────────────────────────────────────────────────────────────────

    private const string SnpReference = "ACGTACGTAC";
    private const string SnpQuery = "ACTTACGAAC"; // substitutions at positions 2 and 7

    #region VARIANT-SNP-001 SHIFT — a prepended flank shifts SNP positions

    [Test]
    [Description("SHIFT: prepending an identical flank to both sequences shifts every SNP position by the flank length.")]
    public void Snps_PrependFlank_ShiftsPositions()
    {
        var original = VariantCaller.FindSnpsDirect(SnpReference, SnpQuery).Select(v => v.Position).ToList();
        original.Should().Equal(2, 7);

        foreach (var flank in new[] { "GG", "TTTTT" })
        {
            var shifted = VariantCaller.FindSnpsDirect(flank + SnpReference, flank + SnpQuery).Select(v => v.Position).ToList();
            shifted.Should().Equal(original.Select(p => p + flank.Length),
                because: $"the identical {flank.Length}-base flank shifts every SNP by {flank.Length}");
        }
    }

    #endregion

    #region VARIANT-SNP-001 INV — SNP detection is deterministic

    [Test]
    [Description("INV: SNP detection is a pure function, so repeated calls return the identical SNP positions.")]
    public void Snps_Deterministic()
    {
        VariantCaller.FindSnpsDirect(SnpReference, SnpQuery).Select(v => v.Position).ToList()
            .Should().Equal(VariantCaller.FindSnpsDirect(SnpReference, SnpQuery).Select(v => v.Position).ToList(),
                because: "SNP calling has no hidden state");
    }

    #endregion
}
