using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for variant calling and annotation (VariantCaller).
///
/// Test Units: VARIANT-ANNOT-001, VARIANT-CALL-001, VARIANT-INDEL-001, VARIANT-SNP-001
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

    #region VARIANT-CALL-001: P: variants only where sequences differ; R: positions valid; D: deterministic

    // CallVariants aligns reference and query then emits a variant at each differing alignment column.
    // (There is no read pileup/depth model here, so the checklist's depth→confidence is N/A.)

    /// <summary>Two gap-free DNA sequences of equal length L (12..18).</summary>
    private static Arbitrary<(string a, string b)> EqualLengthPairArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int len = 12 + rng.Next(7);
            return (RandDna(rng, len), RandDna(rng, len));
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (P): a sequence compared to itself produces no variants.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Call_Identity_HasNoVariants()
    {
        return Prop.ForAll(EqualLengthPairArbitrary(), input =>
        {
            var (a, _) = input;
            return GenomicCallEmpty(a).Label("identical sequences produced variants");
        });
    }

    private static bool GenomicCallEmpty(string a) =>
        !VariantCaller.CallVariants(new DnaSequence(a), new DnaSequence(a)).Any();

    /// <summary>
    /// INV-2 (P + R): over a gap-free alignment, variants are exactly the differing columns — each is a
    /// SNP at the differing position with the correct reference/alternate bases.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Call_FromAlignment_OnlyAtDifferences()
    {
        return Prop.ForAll(EqualLengthPairArbitrary(), input =>
        {
            var (a, b) = input;
            var variants = VariantCaller.CallVariantsFromAlignment(a, b).ToList();
            var expectedPositions = Enumerable.Range(0, a.Length).Where(i => a[i] != b[i]).ToList();

            bool positionsMatch = variants.Select(v => v.Position).SequenceEqual(expectedPositions);
            bool wellFormed = variants.All(v =>
                v.Type == VariantType.SNP &&
                v.Position >= 0 && v.Position < a.Length &&
                v.ReferenceAllele == a[v.Position].ToString() &&
                v.AlternateAllele == b[v.Position].ToString());
            return (positionsMatch && wellFormed).Label("variants were not exactly the differing columns");
        });
    }

    /// <summary>
    /// INV-3 (R): every called variant has a position within the reference span.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Call_Positions_AreValid()
    {
        return Prop.ForAll(SeqPairArbitrary(), input =>
        {
            var (r, q) = input;
            var variants = VariantCaller.CallVariants(new DnaSequence(r), new DnaSequence(q)).ToList();
            return variants.All(v => v.Position >= 0 && v.Position <= r.Length)
                .Label("a variant position was out of range");
        });
    }

    /// <summary>
    /// INV-4 (D): Variant calling is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Call_IsDeterministic()
    {
        return Prop.ForAll(SeqPairArbitrary(), input =>
        {
            var (r, q) = input;
            var a = VariantCaller.CallVariants(new DnaSequence(r), new DnaSequence(q)).ToList();
            var b = VariantCaller.CallVariants(new DnaSequence(r), new DnaSequence(q)).ToList();
            return (a.Count == b.Count && a.Zip(b).All(p => p.First == p.Second))
                .Label("CallVariants must be deterministic");
        });
    }

    #endregion

    #region VARIANT-INDEL-001: P: indel has a gap allele on exactly one side; R: positions valid; D: deterministic

    // FindIndels returns the insertion/deletion variants. In this representation an indel carries the
    // gap allele "-" on exactly one side (Ref="-" for insertions, Alt="-" for deletions).

    private static bool IsWellFormedIndel(Variant v) => v.Type switch
    {
        VariantType.Insertion => v.ReferenceAllele == "-" && v.AlternateAllele.Length == 1 && v.AlternateAllele != "-",
        VariantType.Deletion => v.AlternateAllele == "-" && v.ReferenceAllele.Length == 1 && v.ReferenceAllele != "-",
        _ => false,
    };

    /// <summary>
    /// INV-1 (P + R): every reported indel is a well-formed gap-allele insertion or deletion at a
    /// valid position.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Indels_AreWellFormed()
    {
        return Prop.ForAll(SeqPairArbitrary(), input =>
        {
            var (r, q) = input;
            var indels = VariantCaller.FindIndels(new DnaSequence(r), new DnaSequence(q)).ToList();
            return indels.All(v => IsWellFormedIndel(v) && v.Position >= 0 && v.Position <= r.Length)
                .Label("a reported indel was malformed");
        });
    }

    /// <summary>
    /// INV-2 (P, positive controls): a single-base deletion and insertion are detected with the
    /// correct gap-allele orientation.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Indels_SingleBaseDelAndIns_AreDetected()
    {
        const string reference = "ACGTACGTAC";
        string deletion = reference.Remove(4, 1);            // query missing one base → deletion
        string insertion = reference.Insert(4, "G");          // query has an extra base → insertion

        var dels = VariantCaller.FindIndels(new DnaSequence(reference), new DnaSequence(deletion)).ToList();
        var ins = VariantCaller.FindIndels(new DnaSequence(reference), new DnaSequence(insertion)).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(dels.Any(v => v.Type == VariantType.Deletion), Is.True, "deletion expected");
            Assert.That(ins.Any(v => v.Type == VariantType.Insertion), Is.True, "insertion expected");
            Assert.That(dels.Concat(ins).All(IsWellFormedIndel), Is.True, "indels well-formed");
        });
    }

    /// <summary>
    /// INV-3 (D): Indel detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Indels_IsDeterministic()
    {
        return Prop.ForAll(SeqPairArbitrary(), input =>
        {
            var (r, q) = input;
            var a = VariantCaller.FindIndels(new DnaSequence(r), new DnaSequence(q)).ToList();
            var b = VariantCaller.FindIndels(new DnaSequence(r), new DnaSequence(q)).ToList();
            return (a.Count == b.Count && a.Zip(b).All(p => p.First == p.Second))
                .Label("FindIndels must be deterministic");
        });
    }

    #endregion

    #region VARIANT-SNP-001: P: ref/alt single bases with ref≠alt; R: positions valid; D: deterministic

    private static bool IsWellFormedSnp(Variant v) =>
        v.Type == VariantType.SNP &&
        v.ReferenceAllele.Length == 1 && v.AlternateAllele.Length == 1 &&
        v.ReferenceAllele != v.AlternateAllele &&
        "ACGT".Contains(v.ReferenceAllele[0]) && "ACGT".Contains(v.AlternateAllele[0]);

    /// <summary>
    /// INV-1 (P + completeness): FindSnpsDirect reports a single-base ref≠alt SNP at exactly each
    /// position where the equal-length sequences differ.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Snps_Direct_AreExactSingleBaseDifferences()
    {
        return Prop.ForAll(EqualLengthPairArbitrary(), input =>
        {
            var (a, b) = input;
            var snps = VariantCaller.FindSnpsDirect(a, b).ToList();
            var expected = Enumerable.Range(0, a.Length).Where(i => a[i] != b[i]).ToList();
            bool ok = snps.Select(v => v.Position).SequenceEqual(expected)
                      && snps.All(v => IsWellFormedSnp(v) && v.ReferenceAllele[0] == a[v.Position] && v.AlternateAllele[0] == b[v.Position]);
            return ok.Label("SNPs were not exactly the single-base differences");
        });
    }

    /// <summary>
    /// INV-2 (P + R): every alignment-based SNP is a single-base ref≠alt change at a valid position.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Snps_AreWellFormed()
    {
        return Prop.ForAll(SeqPairArbitrary(), input =>
        {
            var (r, q) = input;
            var snps = VariantCaller.FindSnps(new DnaSequence(r), new DnaSequence(q)).ToList();
            return snps.All(v => IsWellFormedSnp(v) && v.Position >= 0 && v.Position < r.Length)
                .Label("a SNP was malformed or out of range");
        });
    }

    /// <summary>
    /// INV-3 (D): SNP detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Snps_IsDeterministic()
    {
        return Prop.ForAll(EqualLengthPairArbitrary(), input =>
        {
            var (a, b) = input;
            return VariantCaller.FindSnpsDirect(a, b).SequenceEqual(VariantCaller.FindSnpsDirect(a, b))
                .Label("FindSnpsDirect must be deterministic");
        });
    }

    #endregion
}
