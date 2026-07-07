// VARIANT-CALL-001 — Variant Detection
// Evidence: docs/Evidence/VARIANT-CALL-001-Evidence.md
// TestSpec: tests/TestSpecs/VARIANT-CALL-001.md
// Source: VCFv4.3 spec (samtools/hts-specs); Danecek et al. (2011) Bioinformatics 27(15):2156-2158;
//         Collins & Jukes (1994) Genomics 20(3):386-396; Wikipedia Transition/Transversion (Futuyma 2013).

namespace Seqeron.Genomics.Tests.Unit.Annotation;

[TestFixture]
public class VariantCaller_CallVariants_Tests
{
    #region CallVariants

    // M1 — Source 1 (variant = difference from reference); INV-01.
    [Test]
    public void CallVariants_IdenticalSequences_ReturnsNoVariants()
    {
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("ATGCATGC");

        var variants = VariantCaller.CallVariants(reference, query).ToList();

        Assert.That(variants, Is.Empty,
            "Identical reference and query have no differences, so zero variants must be reported (Danecek 2011: variant = difference from reference).");
    }

    // M3 — Source 2 (simple SNP G->A is a single-base substitution); INV-02/04.
    [Test]
    public void CallVariants_SingleSubstitution_ReturnsOneSnpWithExactAlleles()
    {
        var reference = new DnaSequence("ATGC");
        var query = new DnaSequence("ATTC");

        var variants = VariantCaller.CallVariants(reference, query).ToList();

        Assert.That(variants, Has.Count.EqualTo(1), "A single substituted base is exactly one SNP.");
        Assert.Multiple(() =>
        {
            Assert.That(variants[0].Type, Is.EqualTo(VariantType.SNP), "A mismatch column is a SNP.");
            Assert.That(variants[0].Position, Is.EqualTo(2), "The substituted base is at 0-based reference position 2.");
            Assert.That(variants[0].ReferenceAllele, Is.EqualTo("G"), "Reference base at position 2 is G.");
            Assert.That(variants[0].AlternateAllele, Is.EqualTo("T"), "Query base at position 2 is T.");
        });
    }

    // S3 — input validation.
    [Test]
    public void CallVariants_NullReference_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantCaller.CallVariants(null!, new DnaSequence("ATGC")).ToList(),
            "A null reference is invalid input and must throw ArgumentNullException.");
    }

    // S4 — input validation.
    [Test]
    public void CallVariants_NullQuery_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantCaller.CallVariants(new DnaSequence("ATGC"), null!).ToList(),
            "A null query is invalid input and must throw ArgumentNullException.");
    }

    #endregion

    #region CallVariantsFromAlignment

    // M2 — Source 2 (simple SNP); INV-02/04.
    [Test]
    public void CallVariantsFromAlignment_SingleMismatchColumn_ReturnsExactSnp()
    {
        var variants = VariantCaller.CallVariantsFromAlignment("ATGC", "ATTC").ToList();

        Assert.That(variants, Has.Count.EqualTo(1), "One mismatched column yields exactly one SNP.");
        Assert.Multiple(() =>
        {
            Assert.That(variants[0].Type, Is.EqualTo(VariantType.SNP));
            Assert.That(variants[0].Position, Is.EqualTo(2), "0-based reference position of the mismatch.");
            Assert.That(variants[0].ReferenceAllele, Is.EqualTo("G"));
            Assert.That(variants[0].AlternateAllele, Is.EqualTo("T"));
        });
    }

    // M4 — Source 2 (microsatellite insertion of one base T); INV-04, ASM-01/02.
    [Test]
    public void CallVariantsFromAlignment_RefGapColumn_ReturnsInsertionWithGapRefAllele()
    {
        // Reference-side gap = a base present in the query but absent in the reference (insertion).
        var variants = VariantCaller.CallVariantsFromAlignment("AT-GC", "ATTGC").ToList();

        Assert.That(variants, Has.Count.EqualTo(1), "A single reference-gap column is exactly one insertion.");
        Assert.Multiple(() =>
        {
            Assert.That(variants[0].Type, Is.EqualTo(VariantType.Insertion), "A reference-side gap is an insertion in the query.");
            Assert.That(variants[0].Position, Is.EqualTo(2), "Insertion sits after the 2 consumed reference bases (0-based refPos = 2).");
            Assert.That(variants[0].ReferenceAllele, Is.EqualTo("-"), "Inserted base has no reference allele (gap sentinel).");
            Assert.That(variants[0].AlternateAllele, Is.EqualTo("T"), "The inserted base is T.");
            Assert.That(variants[0].QueryPosition, Is.EqualTo(2), "Inserted base is the 3rd query base (0-based query position 2).");
        });
    }

    // M5 — Source 2 (microsatellite deletion); INV-04, ASM-01/02.
    [Test]
    public void CallVariantsFromAlignment_QueryGapColumn_ReturnsDeletionWithGapAltAllele()
    {
        var variants = VariantCaller.CallVariantsFromAlignment("ATTGC", "AT-GC").ToList();

        Assert.That(variants, Has.Count.EqualTo(1), "A single query-gap column is exactly one deletion.");
        Assert.Multiple(() =>
        {
            Assert.That(variants[0].Type, Is.EqualTo(VariantType.Deletion), "A query-side gap is a deletion from the query.");
            Assert.That(variants[0].Position, Is.EqualTo(2), "Deleted reference base is at 0-based position 2.");
            Assert.That(variants[0].ReferenceAllele, Is.EqualTo("T"), "The deleted reference base is T.");
            Assert.That(variants[0].AlternateAllele, Is.EqualTo("-"), "Deleted base has no query allele (gap sentinel).");
        });
    }

    // M6 — Source 2 §1.1 (microsatellite "deletion of 2 bases (TC)"); INV-04.
    [Test]
    public void CallVariantsFromAlignment_TwoBaseDeletion_ReturnsTwoConsecutiveDeletions()
    {
        // GTC (ref) vs G (query) deletes T then C — VCFv4.3 §1.1 microsatellite example.
        var variants = VariantCaller.CallVariantsFromAlignment("GTCAA", "G--AA").ToList();

        Assert.That(variants, Has.Count.EqualTo(2), "Deleting TC produces two single-base deletion columns.");
        Assert.Multiple(() =>
        {
            Assert.That(variants[0].Type, Is.EqualTo(VariantType.Deletion));
            Assert.That(variants[0].Position, Is.EqualTo(1), "First deleted base T is at reference position 1.");
            Assert.That(variants[0].ReferenceAllele, Is.EqualTo("T"), "First deleted base is T.");
            Assert.That(variants[1].Type, Is.EqualTo(VariantType.Deletion));
            Assert.That(variants[1].Position, Is.EqualTo(2), "Second deleted base C is at reference position 2.");
            Assert.That(variants[1].ReferenceAllele, Is.EqualTo("C"), "Second deleted base is C.");
        });
    }

    // S5 — INV-04 (each mismatch column an independent SNP).
    [Test]
    public void CallVariantsFromAlignment_MultipleMismatches_ReturnsSnpsAtExactPositions()
    {
        var variants = VariantCaller.CallVariantsFromAlignment("AAAA", "TGTA").ToList();

        Assert.That(variants, Has.Count.EqualTo(3), "Three differing columns (positions 0,1,2) yield three SNPs.");
        Assert.Multiple(() =>
        {
            Assert.That(variants.Select(v => v.Position), Is.EqualTo(new[] { 0, 1, 2 }), "SNP reference positions.");
            Assert.That(variants.All(v => v.Type == VariantType.SNP), "All three are SNPs.");
            Assert.That(variants.Select(v => v.AlternateAllele), Is.EqualTo(new[] { "T", "G", "T" }), "Query bases at the mismatches.");
        });
    }

    // S1 — documented contract.
    [Test]
    public void CallVariantsFromAlignment_UnequalLengths_ThrowsArgumentException()
    {
        Assert.Throws<ArgumentException>(
            () => VariantCaller.CallVariantsFromAlignment("ATGC", "ATG").ToList(),
            "Aligned sequences must have equal length; unequal columns are an error.");
    }

    // S2 — documented contract.
    [Test]
    public void CallVariantsFromAlignment_EmptyInput_ReturnsEmpty()
    {
        var variants = VariantCaller.CallVariantsFromAlignment("", "").ToList();

        Assert.That(variants, Is.Empty, "Empty alignment has no columns and therefore no variants.");
    }

    // C1 — property test (O(n×m) algorithm): INV-02 and INV-03 over a constructed alignment.
    [Test]
    public void CallVariantsFromAlignment_AnyAlignment_PositionsInBoundsAndSnpsDiffer()
    {
        // Deterministic mixed alignment: SNP, insertion, deletion, matches.
        const string alignedRef = "AC-GTAC";
        const string alignedQuery = "AGTG-AC";
        int refLength = alignedRef.Count(c => c != '-');

        var variants = VariantCaller.CallVariantsFromAlignment(alignedRef, alignedQuery).ToList();

        Assert.That(variants, Is.Not.Empty, "The constructed alignment contains differences.");
        Assert.Multiple(() =>
        {
            foreach (var v in variants)
            {
                Assert.That(v.Position, Is.InRange(0, refLength),
                    $"Every variant position must lie within reference bounds [0,{refLength}] (INV-03).");
                if (v.Type == VariantType.SNP)
                    Assert.That(v.ReferenceAllele, Is.Not.EqualTo(v.AlternateAllele),
                        "A SNP must have distinct reference and alternate alleles (INV-02).");
            }
        });
    }

    #endregion

    #region ClassifyMutation

    // M7 — Source 5 (transition A<->G); INV-05.
    [Test]
    public void ClassifyMutation_AtoG_ReturnsTransition()
    {
        var result = VariantCaller.ClassifyMutation(new Variant(0, "A", "G", VariantType.SNP, 0));
        Assert.That(result, Is.EqualTo(MutationType.Transition), "A->G is purine<->purine, a transition.");
    }

    // M8 — Source 5 (transition C<->T); INV-05.
    [Test]
    public void ClassifyMutation_CtoT_ReturnsTransition()
    {
        var result = VariantCaller.ClassifyMutation(new Variant(0, "C", "T", VariantType.SNP, 0));
        Assert.That(result, Is.EqualTo(MutationType.Transition), "C->T is pyrimidine<->pyrimidine, a transition.");
    }

    // M9 — Source 5; INV-05.
    [Test]
    public void ClassifyMutation_GtoAandTtoC_ReturnTransition()
    {
        Assert.Multiple(() =>
        {
            Assert.That(VariantCaller.ClassifyMutation(new Variant(0, "G", "A", VariantType.SNP, 0)),
                Is.EqualTo(MutationType.Transition), "G->A is a transition (purine<->purine).");
            Assert.That(VariantCaller.ClassifyMutation(new Variant(0, "T", "C", VariantType.SNP, 0)),
                Is.EqualTo(MutationType.Transition), "T->C is a transition (pyrimidine<->pyrimidine).");
        });
    }

    // M10 — Source 6 (the four transversions A<->C, A<->T, G<->C, G<->T); INV-05.
    [Test]
    public void ClassifyMutation_PurinePyrimidineSwaps_ReturnTransversion()
    {
        Assert.Multiple(() =>
        {
            Assert.That(VariantCaller.ClassifyMutation(new Variant(0, "A", "C", VariantType.SNP, 0)),
                Is.EqualTo(MutationType.Transversion), "A->C is purine<->pyrimidine, a transversion.");
            Assert.That(VariantCaller.ClassifyMutation(new Variant(0, "A", "T", VariantType.SNP, 0)),
                Is.EqualTo(MutationType.Transversion), "A->T transversion.");
            Assert.That(VariantCaller.ClassifyMutation(new Variant(0, "G", "C", VariantType.SNP, 0)),
                Is.EqualTo(MutationType.Transversion), "G->C transversion.");
            Assert.That(VariantCaller.ClassifyMutation(new Variant(0, "G", "T", VariantType.SNP, 0)),
                Is.EqualTo(MutationType.Transversion), "G->T transversion.");
        });
    }

    // M11 — Source 2 (REF/ALT case-insensitive); INV-05.
    [Test]
    public void ClassifyMutation_LowercaseBases_ClassifiedCaseInsensitively()
    {
        var result = VariantCaller.ClassifyMutation(new Variant(0, "a", "g", VariantType.SNP, 0));
        Assert.That(result, Is.EqualTo(MutationType.Transition),
            "a->g must classify identically to A->G; VCF base comparison is case-insensitive.");
    }

    // M12 — INV-05 (classification defined only for SNPs).
    [Test]
    public void ClassifyMutation_NonSnp_ReturnsOther()
    {
        var result = VariantCaller.ClassifyMutation(new Variant(0, "A", "-", VariantType.Deletion, 0));
        Assert.That(result, Is.EqualTo(MutationType.Other), "Transition/transversion is defined only for SNPs; indels are Other.");
    }

    #endregion

    #region CalculateTiTvRatio

    // M13 — definition; INV-06.
    [Test]
    public void CalculateTiTvRatio_OneTransitionOneTransversion_ReturnsOne()
    {
        var variants = new[]
        {
            new Variant(0, "A", "G", VariantType.SNP, 0), // transition
            new Variant(1, "A", "C", VariantType.SNP, 1)  // transversion
        };

        double ratio = VariantCaller.CalculateTiTvRatio(variants);

        Assert.That(ratio, Is.EqualTo(1.0).Within(1e-10), "1 transition / 1 transversion = 1.0.");
    }

    // M14 — definition; INV-06.
    [Test]
    public void CalculateTiTvRatio_TwoTransitionsOneTransversion_ReturnsTwo()
    {
        var variants = new[]
        {
            new Variant(0, "A", "G", VariantType.SNP, 0), // transition
            new Variant(1, "C", "T", VariantType.SNP, 1), // transition
            new Variant(2, "A", "C", VariantType.SNP, 2)  // transversion
        };

        double ratio = VariantCaller.CalculateTiTvRatio(variants);

        Assert.That(ratio, Is.EqualTo(2.0).Within(1e-10), "2 transitions / 1 transversion = 2.0.");
    }

    // M15 — ASM-03 (undefined denominator maps to 0); INV-06.
    [Test]
    public void CalculateTiTvRatio_NoTransversions_ReturnsZero()
    {
        var variants = new[]
        {
            new Variant(0, "A", "G", VariantType.SNP, 0),
            new Variant(1, "C", "T", VariantType.SNP, 1)
        };

        double ratio = VariantCaller.CalculateTiTvRatio(variants);

        Assert.That(ratio, Is.EqualTo(0.0).Within(1e-10),
            "With zero transversions the ratio is undefined; the contract returns 0 rather than infinity.");
    }

    // S (validation) — null input.
    [Test]
    public void CalculateTiTvRatio_NullVariants_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantCaller.CalculateTiTvRatio(null!),
            "Null variant collection is invalid input.");
    }

    #endregion
}
