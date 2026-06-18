using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for variant calling and annotation (VariantCaller).
///
/// Test Units: VARIANT-ANNOT-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Variants")]
public class VariantProperties
{
    /// <summary>Two independent DNA sequences of length 12..18.</summary>
    private static Arbitrary<(string reference, string query)> SeqPairArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            return (RandDna(rng, 12 + rng.Next(7)), RandDna(rng, 12 + rng.Next(7)));
        }).ToArbitrary();

    private static string RandDna(Random rng, int len)
    {
        const string bases = "ACGT";
        var c = new char[len];
        for (int i = 0; i < len; i++) c[i] = bases[rng.Next(4)];
        return new string(c);
    }

    #region VARIANT-ANNOT-001: R: effect ∈ enum; P: annotation preserves the called variant; D: deterministic

    // AnnotateVariants wraps each called variant with a functional effect (VariantEffect) and a
    // mutation type. The underlying variant (position, alleles, type) is preserved.

    /// <summary>
    /// INV-1 (P): annotation preserves the called variants exactly — same count and same underlying
    /// Variant in order.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annotate_PreservesCalledVariants()
    {
        return Prop.ForAll(SeqPairArbitrary(), input =>
        {
            var (r, q) = input;
            var called = VariantCaller.CallVariants(new DnaSequence(r), new DnaSequence(q)).ToList();
            var annotated = VariantCaller.AnnotateVariants(new DnaSequence(r), new DnaSequence(q)).ToList();
            bool ok = annotated.Count == called.Count
                      && annotated.Zip(called).All(p => p.First.Variant == p.Second);
            return ok.Label("annotation did not preserve the called variants");
        });
    }

    /// <summary>
    /// INV-2 (R): every annotation has a defined VariantEffect and a mutation type consistent with the
    /// variant type (SNP → transition/transversion; otherwise Other).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annotate_EffectsAndMutationTypesAreValid()
    {
        return Prop.ForAll(SeqPairArbitrary(), input =>
        {
            var (r, q) = input;
            var annotated = VariantCaller.AnnotateVariants(new DnaSequence(r), new DnaSequence(q)).ToList();
            bool ok = annotated.All(a =>
                Enum.IsDefined(a.Effect) &&
                (a.Variant.Type == VariantType.SNP
                    ? a.MutationType is MutationType.Transition or MutationType.Transversion
                    : a.MutationType == MutationType.Other));
            return ok.Label("an annotation had an undefined effect or inconsistent mutation type");
        });
    }

    /// <summary>
    /// INV-3 (P, coding effects): in coding mode, indels are Frameshift and SNPs map to a coding effect.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annotate_CodingMode_EffectsAreConsistent()
    {
        return Prop.ForAll(SeqPairArbitrary(), input =>
        {
            var (r, q) = input;
            var annotated = VariantCaller.AnnotateVariants(new DnaSequence(r), new DnaSequence(q), isCodingSequence: true).ToList();
            bool ok = annotated.All(a => a.Variant.Type switch
            {
                VariantType.Insertion or VariantType.Deletion => a.Effect == VariantEffect.Frameshift,
                VariantType.SNP => a.Effect is VariantEffect.Synonymous or VariantEffect.Missense
                    or VariantEffect.Nonsense or VariantEffect.StopLoss or VariantEffect.Unknown,
                _ => true,
            });
            return ok.Label("a coding-mode effect was inconsistent with the variant type");
        });
    }

    /// <summary>
    /// INV-4 (D): Annotation is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annotate_IsDeterministic()
    {
        return Prop.ForAll(SeqPairArbitrary(), input =>
        {
            var (r, q) = input;
            var a = VariantCaller.AnnotateVariants(new DnaSequence(r), new DnaSequence(q)).ToList();
            var b = VariantCaller.AnnotateVariants(new DnaSequence(r), new DnaSequence(q)).ToList();
            return (a.Count == b.Count && a.Zip(b).All(p => p.First == p.Second))
                .Label("AnnotateVariants must be deterministic");
        });
    }

    #endregion
}
