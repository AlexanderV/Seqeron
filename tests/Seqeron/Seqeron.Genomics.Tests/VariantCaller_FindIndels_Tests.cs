// VARIANT-INDEL-001 — Indel Detection
// Evidence: docs/Evidence/VARIANT-INDEL-001-Evidence.md
// TestSpec: tests/TestSpecs/VARIANT-INDEL-001.md
// Source: VCFv4.3 spec (samtools/hts-specs), REF field & §1.1 indel examples;
//         Tan A, Abecasis GR, Kang HM (2015). Bioinformatics 31(13):2202-2204, doi:10.1093/bioinformatics/btv112;
//         ericminikel/minimal_representation normalize.py (Tan et al. 2015 Algorithm 1).

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class VariantCaller_FindIndels_Tests
{
    // Gap sentinel for the absent allele in the in-memory model (ASM-01).
    private const string Gap = "-";

    #region FindInsertions

    // M1 — VCF single-base insertion (REF=C, ALT=CA: ALT longer than REF); INV-03. Unique alignment.
    [Test]
    public void FindInsertions_SingleInsertedBase_ReturnsExactInsertion()
    {
        // ref "ATGCAT" vs query "ATGTCAT": a T is inserted after reference index 2 (between G and C).
        // The optimal global alignment ref "ATG-CAT" / qry "ATGTCAT" (6 matches, 1 gap) is unique,
        // so the indel position is deterministic.
        var reference = new DnaSequence("ATGCAT");
        var query = new DnaSequence("ATGTCAT");

        var insertions = VariantCaller.FindInsertions(reference, query).ToList();

        Assert.That(insertions, Has.Count.EqualTo(1),
            "One base present in the query but absent from the reference is exactly one insertion (VCFv4.3 insertion class).");
        Assert.Multiple(() =>
        {
            Assert.That(insertions[0].Type, Is.EqualTo(VariantType.Insertion),
                "A reference-gap column is an insertion in the query (VCFv4.3: insertion REF=C, ALT=CA).");
            Assert.That(insertions[0].ReferenceAllele, Is.EqualTo(Gap),
                "The reference has no base at an insertion column, so the in-memory reference allele is the gap sentinel (ASM-01).");
            Assert.That(insertions[0].AlternateAllele, Is.EqualTo("T"),
                "The inserted base is T; ALT is longer than REF for an insertion (VCFv4.3).");
            Assert.That(insertions[0].Position, Is.EqualTo(3),
                "After consuming reference ATG (indices 0-2), the insertion is reported at 0-based reference position 3 (ASM-02 unique alignment).");
        });
    }

    // M3 — FindInsertions filters to insertions only (no SNP, no deletion); INV-02.
    [Test]
    public void FindInsertions_InsertionAndSubstitutionInput_ReturnsInsertionsOnly()
    {
        // ref "ATGCATGC" vs query "ATGTCATGG": a T inserted after index 2 AND a substitution C->G near the end.
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("ATGTCATGG");

        var insertions = VariantCaller.FindInsertions(reference, query).ToList();

        Assert.That(insertions, Is.Not.Empty, "The input contains at least one insertion.");
        Assert.That(insertions.All(v => v.Type == VariantType.Insertion),
            "FindInsertions returns insertions only; SNPs and deletions are filtered out (insertion is a distinct VCF class, INV-02).");
    }

    // M5 — multi-base insertion = consecutive per-base columns; INV-05. Unique alignment.
    [Test]
    public void FindInsertions_MultiBaseInsertion_ReturnsConsecutiveColumns()
    {
        // ref "ATGCAT" vs query "ATGTTCAT": two T bases inserted after reference index 2.
        // Optimal alignment ref "ATG--CAT" / qry "ATGTTCAT" (6 matches, 2 gaps) is unique.
        var reference = new DnaSequence("ATGCAT");
        var query = new DnaSequence("ATGTTCAT");

        var insertions = VariantCaller.FindInsertions(reference, query).ToList();

        Assert.That(insertions, Has.Count.EqualTo(2),
            "A contiguous block of 2 inserted bases yields exactly 2 insertion columns (each extra base is one event, INV-05).");
        Assert.Multiple(() =>
        {
            Assert.That(insertions.All(v => v.Type == VariantType.Insertion), "Both reported variants are insertions.");
            Assert.That(insertions.All(v => v.ReferenceAllele == Gap), "Both insertion columns have the gap sentinel as reference allele (ASM-01).");
            Assert.That(insertions.Select(v => v.AlternateAllele), Is.EqualTo(new[] { "T", "T" }), "Both inserted bases are T.");
            Assert.That(insertions.Select(v => v.Position), Is.EqualTo(new[] { 3, 3 }),
                "Both inserted bases sit at the same reference junction (0-based reference position 3) since the reference is not consumed by an insertion column.");
        });
    }

    // M7 — identical sequences produce no insertions; INV-01.
    [Test]
    public void FindInsertions_IdenticalSequences_ReturnsEmpty()
    {
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("ATGCATGC");

        var insertions = VariantCaller.FindInsertions(reference, query).ToList();

        Assert.That(insertions, Is.Empty,
            "Identical sequences have no length-changing difference, so there are no insertions (VCFv4.3 indel class, INV-01).");
    }

    // S5 — a length-preserving substitution is not an insertion.
    [Test]
    public void FindInsertions_SubstitutionOnlyInput_ReturnsEmpty()
    {
        // Equal-length input differing only at index 3 (C->A): a substitution, not an indel.
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("ATGAATGC");

        var insertions = VariantCaller.FindInsertions(reference, query).ToList();

        Assert.That(insertions, Is.Empty,
            "A substitution preserves length and is a SNP, not an insertion (VCFv4.3 distinct classes).");
    }

    // S1 — null reference propagates ArgumentNullException from CallVariants.
    [Test]
    public void FindInsertions_NullReference_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantCaller.FindInsertions(null!, new DnaSequence("ATGC")).ToList(),
            "A null reference is invalid input and must throw ArgumentNullException.");
    }

    // S2 — null query propagates ArgumentNullException from CallVariants.
    [Test]
    public void FindInsertions_NullQuery_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantCaller.FindInsertions(new DnaSequence("ATGC"), null!).ToList(),
            "A null query is invalid input and must throw ArgumentNullException.");
    }

    #endregion

    #region FindDeletions

    // M2 — VCF single-base deletion (REF=TC, ALT=T: REF longer than ALT); INV-04. Unique alignment.
    [Test]
    public void FindDeletions_SingleDeletedBase_ReturnsExactDeletion()
    {
        // ref "ATGTCAT" vs query "ATGCAT": the T at reference index 3 is deleted.
        // Optimal alignment ref "ATGTCAT" / qry "ATG-CAT" (6 matches, 1 gap) is unique.
        var reference = new DnaSequence("ATGTCAT");
        var query = new DnaSequence("ATGCAT");

        var deletions = VariantCaller.FindDeletions(reference, query).ToList();

        Assert.That(deletions, Has.Count.EqualTo(1),
            "One reference base absent from the query is exactly one deletion (VCFv4.3 deletion class).");
        Assert.Multiple(() =>
        {
            Assert.That(deletions[0].Type, Is.EqualTo(VariantType.Deletion),
                "A query-gap column is a deletion in the query (VCFv4.3: deletion REF=TC, ALT=T).");
            Assert.That(deletions[0].ReferenceAllele, Is.EqualTo("T"),
                "The deleted reference base is T; REF is longer than ALT for a deletion (VCFv4.3).");
            Assert.That(deletions[0].AlternateAllele, Is.EqualTo(Gap),
                "The query has no base at a deletion column, so the in-memory alternate allele is the gap sentinel (ASM-01).");
            Assert.That(deletions[0].Position, Is.EqualTo(3),
                "After consuming reference ATG (indices 0-2), the deletion is reported at 0-based reference position 3 (ASM-02 unique alignment).");
        });
    }

    // M4 — FindDeletions filters to deletions only (no SNP, no insertion); INV-02.
    [Test]
    public void FindDeletions_DeletionAndSubstitutionInput_ReturnsDeletionsOnly()
    {
        // ref "ATGTCATGC" vs query "ATGCATGG": the T at index 3 deleted AND a substitution C->G near the end.
        var reference = new DnaSequence("ATGTCATGC");
        var query = new DnaSequence("ATGCATGG");

        var deletions = VariantCaller.FindDeletions(reference, query).ToList();

        Assert.That(deletions, Is.Not.Empty, "The input contains at least one deletion.");
        Assert.That(deletions.All(v => v.Type == VariantType.Deletion),
            "FindDeletions returns deletions only; SNPs and insertions are filtered out (deletion is a distinct VCF class, INV-02).");
    }

    // M6 — multi-base deletion = consecutive per-base columns; INV-05. Unique alignment.
    [Test]
    public void FindDeletions_MultiBaseDeletion_ReturnsConsecutiveColumns()
    {
        // ref "ATGTTCAT" vs query "ATGCAT": two T bases deleted at reference indices 3 and 4.
        // Optimal alignment ref "ATGTTCAT" / qry "ATG--CAT" (6 matches, 2 gaps) is unique.
        var reference = new DnaSequence("ATGTTCAT");
        var query = new DnaSequence("ATGCAT");

        var deletions = VariantCaller.FindDeletions(reference, query).ToList();

        Assert.That(deletions, Has.Count.EqualTo(2),
            "A contiguous block of 2 deleted bases yields exactly 2 deletion columns (each absent base is one event, INV-05).");
        Assert.Multiple(() =>
        {
            Assert.That(deletions.All(v => v.Type == VariantType.Deletion), "Both reported variants are deletions.");
            Assert.That(deletions.All(v => v.AlternateAllele == Gap), "Both deletion columns have the gap sentinel as alternate allele (ASM-01).");
            Assert.That(deletions.Select(v => v.ReferenceAllele), Is.EqualTo(new[] { "T", "T" }), "Both deleted reference bases are T.");
            Assert.That(deletions.Select(v => v.Position), Is.EqualTo(new[] { 3, 4 }),
                "The two deleted reference bases occupy consecutive 0-based reference positions 3 and 4 (the reference is consumed by each deletion column).");
        });
    }

    // M8 — identical sequences produce no deletions; INV-01.
    [Test]
    public void FindDeletions_IdenticalSequences_ReturnsEmpty()
    {
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("ATGCATGC");

        var deletions = VariantCaller.FindDeletions(reference, query).ToList();

        Assert.That(deletions, Is.Empty,
            "Identical sequences have no length-changing difference, so there are no deletions (VCFv4.3 indel class, INV-01).");
    }

    // S6 — a length-preserving substitution is not a deletion.
    [Test]
    public void FindDeletions_SubstitutionOnlyInput_ReturnsEmpty()
    {
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("ATGAATGC");

        var deletions = VariantCaller.FindDeletions(reference, query).ToList();

        Assert.That(deletions, Is.Empty,
            "A substitution preserves length and is a SNP, not a deletion (VCFv4.3 distinct classes).");
    }

    // S3 — null reference propagates ArgumentNullException from CallVariants.
    [Test]
    public void FindDeletions_NullReference_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantCaller.FindDeletions(null!, new DnaSequence("ATGC")).ToList(),
            "A null reference is invalid input and must throw ArgumentNullException.");
    }

    // S4 — null query propagates ArgumentNullException from CallVariants.
    [Test]
    public void FindDeletions_NullQuery_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantCaller.FindDeletions(new DnaSequence("ATGC"), null!).ToList(),
            "A null query is invalid input and must throw ArgumentNullException.");
    }

    #endregion

    #region FindIndels (delegate)

    // S7 — delegate smoke test: union of insertions and deletions, no SNPs.
    [Test]
    public void FindIndels_InsertionAndDeletionInput_ReturnsBothTypes()
    {
        // ref "ATGTCATGCAT" vs query "ATGCATGTCAT": one deletion (T at index 3) and one insertion (T) downstream,
        // constructed so the alignment yields both an insertion and a deletion column.
        var reference = new DnaSequence("ATGTCATGCAT");
        var query = new DnaSequence("ATGCATGTCAT");

        var indels = VariantCaller.FindIndels(reference, query).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(indels.All(v => v.Type == VariantType.Insertion || v.Type == VariantType.Deletion),
                "FindIndels returns only insertions and deletions (no SNPs).");
            Assert.That(indels.Any(v => v.Type == VariantType.Insertion), "At least one insertion is reported.");
            Assert.That(indels.Any(v => v.Type == VariantType.Deletion), "At least one deletion is reported.");
        });
    }

    #endregion

    #region Property-based (O(n*m) invariant)

    // C1 — property: a contiguous block of k inserted bases yields exactly k insertion columns,
    //      all of type Insertion and all within reference bounds (INV-05, INV-06). Unique alignment.
    [Test]
    public void FindInsertions_ContiguousBlock_CountEqualsBlockLengthAndPositionsInBounds()
    {
        // ref "ATGCAT" vs query "ATGTTTCAT": a contiguous block of k=3 T bases inserted after index 2.
        // Optimal alignment ref "ATG---CAT" / qry "ATGTTTCAT" (6 matches, 3 gaps) is unique.
        const int k = 3;
        var reference = new DnaSequence("ATGCAT");
        var query = new DnaSequence("ATGTTTCAT");

        var insertions = VariantCaller.FindInsertions(reference, query).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(insertions, Has.Count.EqualTo(k),
                "A contiguous block of k inserted bases yields exactly k insertion columns (INV-05).");
            Assert.That(insertions.All(v => v.Type == VariantType.Insertion), "Every reported variant is an insertion.");
            Assert.That(insertions.All(v => v.Position >= 0 && v.Position <= reference.Length),
                "Every reported indel position lies within [0, reference.Length] (INV-06).");
        });
    }

    #endregion
}
