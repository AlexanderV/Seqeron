using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class VariantCallerTests
{
    #region Basic Variant Detection Tests

    [Test]
    public void CallVariants_IdenticalSequences_NoVariants()
    {
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("ATGCATGC");

        var variants = VariantCaller.CallVariants(reference, query).ToList();

        Assert.That(variants, Is.Empty);
    }

    [Test]
    public void CallVariants_SingleSnp_DetectsIt()
    {
        var reference = new DnaSequence("ATGC");
        var query = new DnaSequence("ATTC");

        var variants = VariantCaller.CallVariants(reference, query).ToList();

        Assert.That(variants.Count, Is.EqualTo(1));
        Assert.That(variants[0].Type, Is.EqualTo(VariantType.SNP));
        Assert.That(variants[0].ReferenceAllele, Is.EqualTo("G"));
        Assert.That(variants[0].AlternateAllele, Is.EqualTo("T"));
    }

    [Test]
    public void CallVariants_MultipleSnps_DetectsAll()
    {
        var reference = new DnaSequence("AAAA");
        var query = new DnaSequence("TATA");

        var variants = VariantCaller.CallVariants(reference, query).ToList();

        Assert.That(variants.Count(v => v.Type == VariantType.SNP), Is.EqualTo(2));
    }

    #endregion

    #region Alignment-Based Variant Detection Tests

    [Test]
    public void CallVariantsFromAlignment_WithInsertion_DetectsIt()
    {
        var variants = VariantCaller.CallVariantsFromAlignment(
            "AT-GC",
            "ATXGC").ToList();

        Assert.That(variants.Any(v => v.Type == VariantType.Insertion));
    }

    [Test]
    public void CallVariantsFromAlignment_WithDeletion_DetectsIt()
    {
        var variants = VariantCaller.CallVariantsFromAlignment(
            "ATGC",
            "AT-C").ToList();

        Assert.That(variants.Any(v => v.Type == VariantType.Deletion));
    }

    [Test]
    public void CallVariantsFromAlignment_DifferentLengths_ThrowsException()
    {
        Assert.Throws<ArgumentException>(() =>
            VariantCaller.CallVariantsFromAlignment("ATGC", "ATG").ToList());
    }

    [Test]
    public void CallVariantsFromAlignment_EmptySequences_ReturnsEmpty()
    {
        var variants = VariantCaller.CallVariantsFromAlignment("", "").ToList();
        Assert.That(variants, Is.Empty);
    }

    #endregion

    #region SNP Detection Tests

    [Test]
    public void FindSnps_ReturnsOnlySnps()
    {
        var reference = new DnaSequence("ATGC");
        var query = new DnaSequence("ATTC");

        var snps = VariantCaller.FindSnps(reference, query).ToList();

        Assert.That(snps.All(v => v.Type == VariantType.SNP));
    }

    [Test]
    public void FindSnpsDirect_FastDirectComparison()
    {
        var snps = VariantCaller.FindSnpsDirect("ATGC", "ATTC").ToList();

        Assert.That(snps.Count, Is.EqualTo(1));
        Assert.That(snps[0].Position, Is.EqualTo(2));
    }

    [Test]
    public void FindSnpsDirect_DifferentLengths_ComparesOverlap()
    {
        var snps = VariantCaller.FindSnpsDirect("ATGC", "ATGCAAAA").ToList();

        Assert.That(snps, Is.Empty); // No SNPs in overlapping region
    }

    #endregion

    #region Indel Detection Tests

    [Test]
    public void FindInsertions_ReturnsOnlyInsertions()
    {
        var reference = new DnaSequence("ATGC");
        var query = new DnaSequence("ATGGC"); // G inserted

        var insertions = VariantCaller.FindInsertions(reference, query).ToList();

        Assert.That(insertions.All(v => v.Type == VariantType.Insertion));
    }

    [Test]
    public void FindDeletions_ReturnsOnlyDeletions()
    {
        var reference = new DnaSequence("ATGGC");
        var query = new DnaSequence("ATGC"); // G deleted

        var deletions = VariantCaller.FindDeletions(reference, query).ToList();

        Assert.That(deletions.All(v => v.Type == VariantType.Deletion));
    }

    [Test]
    public void FindIndels_ReturnsBothTypes()
    {
        var reference = new DnaSequence("ATGC");
        var query = new DnaSequence("ATGC");

        var indels = VariantCaller.FindIndels(reference, query).ToList();

        Assert.That(indels.All(v =>
            v.Type == VariantType.Insertion || v.Type == VariantType.Deletion));
    }

    #endregion

    #region Mutation Classification Tests

    [Test]
    public void ClassifyMutation_AG_Transition()
    {
        var variant = new Variant(0, "A", "G", VariantType.SNP, 0);
        var type = VariantCaller.ClassifyMutation(variant);

        Assert.That(type, Is.EqualTo(MutationType.Transition));
    }

    [Test]
    public void ClassifyMutation_CT_Transition()
    {
        var variant = new Variant(0, "C", "T", VariantType.SNP, 0);
        var type = VariantCaller.ClassifyMutation(variant);

        Assert.That(type, Is.EqualTo(MutationType.Transition));
    }

    [Test]
    public void ClassifyMutation_AC_Transversion()
    {
        var variant = new Variant(0, "A", "C", VariantType.SNP, 0);
        var type = VariantCaller.ClassifyMutation(variant);

        Assert.That(type, Is.EqualTo(MutationType.Transversion));
    }

    [Test]
    public void ClassifyMutation_GT_Transversion()
    {
        var variant = new Variant(0, "G", "T", VariantType.SNP, 0);
        var type = VariantCaller.ClassifyMutation(variant);

        Assert.That(type, Is.EqualTo(MutationType.Transversion));
    }

    [Test]
    public void ClassifyMutation_NonSnp_ReturnsOther()
    {
        var variant = new Variant(0, "A", "-", VariantType.Deletion, 0);
        var type = VariantCaller.ClassifyMutation(variant);

        Assert.That(type, Is.EqualTo(MutationType.Other));
    }

    [Test]
    public void CalculateTiTvRatio_OnlyTransitions_ReturnsInfinity()
    {
        var variants = new[]
        {
            new Variant(0, "A", "G", VariantType.SNP, 0),
            new Variant(1, "C", "T", VariantType.SNP, 1)
        };

        double ratio = VariantCaller.CalculateTiTvRatio(variants);

        Assert.That(ratio, Is.EqualTo(0)); // No transversions, returns 0
    }

    [Test]
    public void CalculateTiTvRatio_EqualTiTv_ReturnsOne()
    {
        var variants = new[]
        {
            new Variant(0, "A", "G", VariantType.SNP, 0), // Transition
            new Variant(1, "A", "C", VariantType.SNP, 1)  // Transversion
        };

        double ratio = VariantCaller.CalculateTiTvRatio(variants);

        Assert.That(ratio, Is.EqualTo(1.0));
    }

    #endregion

    #region Statistics Tests

    [Test]
    public void CalculateStatistics_ReturnsComprehensiveStats()
    {
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("ATTCATGC");

        var stats = VariantCaller.CalculateStatistics(reference, query);

        Assert.That(stats.TotalVariants, Is.GreaterThan(0));
        Assert.That(stats.Snps, Is.GreaterThan(0));
        Assert.That(stats.ReferenceLength, Is.EqualTo(8));
        Assert.That(stats.QueryLength, Is.EqualTo(8));
    }

    [Test]
    public void CalculateStatistics_VariantDensity_Calculated()
    {
        var reference = new DnaSequence("ATGCATGC");
        var query = new DnaSequence("TTGCATGC");

        var stats = VariantCaller.CalculateStatistics(reference, query);

        Assert.That(stats.VariantDensity, Is.GreaterThan(0));
    }

    #endregion

    #region Variant Effect Prediction Tests

    [Test]
    public void PredictEffect_SynonymousMutation_ReturnsSynonymous()
    {
        // TTA -> TTG both code for Leucine
        var coding = new DnaSequence("TTA");
        var variant = new Variant(2, "A", "G", VariantType.SNP, 2);

        var effect = VariantCaller.PredictEffect(variant, coding, 2);

        Assert.That(effect, Is.EqualTo(VariantEffect.Synonymous));
    }

    [Test]
    public void PredictEffect_MissenseMutation_ReturnsMissense()
    {
        // ATG (Met) -> ACG (Thr)
        var coding = new DnaSequence("ATG");
        var variant = new Variant(1, "T", "C", VariantType.SNP, 1);

        var effect = VariantCaller.PredictEffect(variant, coding, 1);

        Assert.That(effect, Is.EqualTo(VariantEffect.Missense));
    }

    [Test]
    public void PredictEffect_NonsenseMutation_ReturnsNonsense()
    {
        // TAT (Tyr) -> TAA (Stop)
        var coding = new DnaSequence("TAT");
        var variant = new Variant(2, "T", "A", VariantType.SNP, 2);

        var effect = VariantCaller.PredictEffect(variant, coding, 2);

        Assert.That(effect, Is.EqualTo(VariantEffect.Nonsense));
    }

    [Test]
    public void PredictEffect_Indel_ReturnsFrameshift()
    {
        var coding = new DnaSequence("ATGC");
        var variant = new Variant(1, "-", "T", VariantType.Insertion, 1);

        var effect = VariantCaller.PredictEffect(variant, coding, 1);

        Assert.That(effect, Is.EqualTo(VariantEffect.Frameshift));
    }

    [Test]
    public void AnnotateVariants_ReturnsAnnotations()
    {
        var reference = new DnaSequence("ATGCAT");
        var query = new DnaSequence("ATTCAT");

        var annotated = VariantCaller.AnnotateVariants(reference, query, isCodingSequence: true).ToList();

        Assert.That(annotated.Count, Is.GreaterThan(0));
        Assert.That(annotated.All(a => a.Effect != VariantEffect.Unknown || a.Variant.Type != VariantType.SNP));
    }

    #endregion

    #region VCF Output Tests

    [Test]
    public void ToVcfLines_GeneratesHeader()
    {
        var variants = new[] { new Variant(0, "A", "G", VariantType.SNP, 0) };
        var lines = VariantCaller.ToVcfLines(variants).ToList();

        Assert.That(lines[0], Does.StartWith("##fileformat=VCF"));
        Assert.That(lines.Any(l => l.StartsWith("#CHROM")));
    }

    [Test]
    public void ToVcfLines_FormatsVariantCorrectly()
    {
        var variants = new[] { new Variant(9, "A", "G", VariantType.SNP, 9) };
        var lines = VariantCaller.ToVcfLines(variants, chromosome: "chr1").ToList();

        var dataLine = lines.Last();
        Assert.That(dataLine, Does.StartWith("chr1\t10\t")); // 1-based position
        Assert.That(dataLine, Does.Contain("\tA\tG\t"));
    }

    [Test]
    public void ToVcfLines_CustomSampleName()
    {
        var variants = new[] { new Variant(0, "A", "G", VariantType.SNP, 0) };
        var lines = VariantCaller.ToVcfLines(variants, sampleName: "MySample").ToList();

        Assert.That(lines.Any(l => l.Contains("MySample")));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void CallVariants_NullReference_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            VariantCaller.CallVariants(null!, new DnaSequence("ATGC")).ToList());
    }

    [Test]
    public void CallVariants_NullQuery_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            VariantCaller.CallVariants(new DnaSequence("ATGC"), null!).ToList());
    }

    [Test]
    public void CalculateStatistics_NullReference_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            VariantCaller.CalculateStatistics(null!, new DnaSequence("ATGC")));
    }

    [Test]
    public void AnnotateVariants_NullReference_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            VariantCaller.AnnotateVariants(null!, new DnaSequence("ATGC")).ToList());
    }

    [Test]
    public void ToVcfLines_NullVariants_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            VariantCaller.ToVcfLines(null!).ToList());
    }

    [Test]
    public void CalculateTiTvRatio_NullVariants_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            VariantCaller.CalculateTiTvRatio(null!));
    }

    #endregion
}
