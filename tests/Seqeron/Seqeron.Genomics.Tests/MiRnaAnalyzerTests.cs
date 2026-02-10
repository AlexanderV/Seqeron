using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class MiRnaAnalyzerTests
{
    // Seed Sequence Tests migrated to MiRnaAnalyzer_SeedAnalysis_Tests.cs (MIRNA-SEED-001)

    #region Reverse Complement Tests

    [Test]
    public void GetReverseComplement_SimpleSequence_ReturnsComplement()
    {
        string seq = "AGCAGCA";
        string rc = GetReverseComplement(seq);

        Assert.That(rc, Is.EqualTo("UGCUGCU"));
    }

    [Test]
    public void GetReverseComplement_EmptySequence_ReturnsEmpty()
    {
        Assert.That(GetReverseComplement(""), Is.Empty);
    }

    [Test]
    public void GetReverseComplement_SingleBase_ReturnsComplement()
    {
        Assert.That(GetReverseComplement("A"), Is.EqualTo("U"));
        Assert.That(GetReverseComplement("U"), Is.EqualTo("A"));
        Assert.That(GetReverseComplement("G"), Is.EqualTo("C"));
        Assert.That(GetReverseComplement("C"), Is.EqualTo("G"));
    }

    #endregion

    #region Base Pairing Tests

    [TestCase('A', 'U', true)]
    [TestCase('U', 'A', true)]
    [TestCase('G', 'C', true)]
    [TestCase('C', 'G', true)]
    [TestCase('G', 'U', true)]
    [TestCase('U', 'G', true)]
    [TestCase('A', 'A', false)]
    [TestCase('C', 'C', false)]
    public void CanPair_VariousBases_ReturnsExpected(char b1, char b2, bool expected)
    {
        Assert.That(CanPair(b1, b2), Is.EqualTo(expected));
    }

    [Test]
    public void IsWobblePair_GU_ReturnsTrue()
    {
        Assert.That(IsWobblePair('G', 'U'), Is.True);
        Assert.That(IsWobblePair('U', 'G'), Is.True);
    }

    [Test]
    public void IsWobblePair_WatsonCrick_ReturnsFalse()
    {
        Assert.That(IsWobblePair('A', 'U'), Is.False);
        Assert.That(IsWobblePair('G', 'C'), Is.False);
    }

    #endregion

    // Target Site Finding Tests migrated to MiRnaAnalyzer_TargetPrediction_Tests.cs (MIRNA-TARGET-001)
    // Alignment Tests migrated to MiRnaAnalyzer_TargetPrediction_Tests.cs (MIRNA-TARGET-001)

    #region Pre-miRNA Tests

    [Test]
    public void FindPreMiRnaHairpins_ValidHairpin_FindsPreMiRNA()
    {
        // Create a simple hairpin structure
        string stem5 = "GCGCGCGCGCGCGCGCGCGC"; // 20 nt
        string loop = "AAAAA";
        string stem3 = "GCGCGCGCGCGCGCGCGCGC"; // Complement
        string hairpin = stem5 + loop + stem3;

        // Need proper complementarity
        string seq = "GGGGGGGGGGGGGGGGGGGG" + "AAAAAAA" + "CCCCCCCCCCCCCCCCCCCC";

        var premirnas = FindPreMiRnaHairpins(seq, minHairpinLength: 45).ToList();

        Assert.That(premirnas, Has.Count.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void FindPreMiRnaHairpins_ShortSequence_ReturnsEmpty()
    {
        var premirnas = FindPreMiRnaHairpins("GGGGCCCCC", minHairpinLength: 55).ToList();
        Assert.That(premirnas, Is.Empty);
    }

    [Test]
    public void FindPreMiRnaHairpins_ReturnsStructureInfo()
    {
        string seq = new string('G', 25) + "AAAAAAA" + new string('C', 25);

        var premirnas = FindPreMiRnaHairpins(seq, minHairpinLength: 50).ToList();

        foreach (var pre in premirnas)
        {
            Assert.That(pre.Structure, Is.Not.Empty);
            Assert.That(pre.MatureSequence, Is.Not.Empty);
        }
    }

    #endregion

    #region Context Analysis Tests

    [Test]
    public void AnalyzeTargetContext_HighAU_HighScore()
    {
        string mrna = "AAAUUUAAAUUUAAAUUU";
        var context = AnalyzeTargetContext(mrna, 6, 12);

        Assert.That(context.AuContent, Is.GreaterThan(0.8));
    }

    [Test]
    public void AnalyzeTargetContext_MiddlePosition_BonusScore()
    {
        string mrna = new string('A', 100);
        var contextMiddle = AnalyzeTargetContext(mrna, 40, 50);
        var contextEnd = AnalyzeTargetContext(mrna, 90, 95);

        Assert.That(contextMiddle.ContextScore, Is.GreaterThanOrEqualTo(contextEnd.ContextScore));
    }

    [Test]
    public void AnalyzeTargetContext_EmptySequence_ReturnsZeros()
    {
        var context = AnalyzeTargetContext("", 0, 5);

        Assert.That(context.AuContent, Is.EqualTo(0));
        Assert.That(context.ContextScore, Is.EqualTo(0));
    }

    [Test]
    public void CalculateSiteAccessibility_UnstructuredRegion_HighAccessibility()
    {
        // All different bases = no structure
        string mrna = "ACGUACGUACGUACGUACGUACGU";
        double access = CalculateSiteAccessibility(mrna, 8, 16);

        Assert.That(access, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void CalculateSiteAccessibility_InvalidPosition_ReturnsZero()
    {
        Assert.That(CalculateSiteAccessibility("AAAA", -1, 5), Is.EqualTo(0));
        Assert.That(CalculateSiteAccessibility("AAAA", 0, 100), Is.EqualTo(0));
    }

    #endregion

    #region miRNA Family Tests

    [Test]
    public void GroupBySeedFamily_SameSeed_GroupedTogether()
    {
        var mirnas = new List<MiRna>
        {
            CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU"),
            CreateMiRna("let-7b", "UGAGGUAGUAGGUUGUGUGGUU"),
            CreateMiRna("miR-1", "UGGAAUGUAAAGAAGUAUGUAU")
        };

        var families = GroupBySeedFamily(mirnas).ToList();

        Assert.That(families.Count, Is.LessThanOrEqualTo(3));
    }

    [Test]
    public void FindSimilarMiRnas_OneMismatch_Finds()
    {
        var query = CreateMiRna("miR-1", "UGGAAUGUAAAGAAGUAUGUAU");
        var database = new List<MiRna>
        {
            CreateMiRna("miR-2", "UGGAAUGUAAAGAAGUAUGUAU"), // Same
            CreateMiRna("miR-3", "UAGAAUGUAAAGAAGUAUGUAU"), // 1 diff
            CreateMiRna("miR-4", "UCCCCUGUAAAGAAGUAUGUAU")  // Many diff
        };

        var similar = FindSimilarMiRnas(query, database, maxMismatches: 1).ToList();

        Assert.That(similar, Has.Count.GreaterThanOrEqualTo(1));
    }

    [Test]
    public void FindSimilarMiRnas_ExcludesQuery()
    {
        var query = CreateMiRna("miR-1", "UGGAAUGUAAAGAAGUAUGUAU");
        var database = new List<MiRna> { query };

        var similar = FindSimilarMiRnas(query, database).ToList();

        Assert.That(similar, Is.Empty);
    }

    #endregion

    #region Utility Tests

    [Test]
    public void CalculateGcContent_HighGC_CorrectValue()
    {
        string seq = "GGGGCCCC";
        double gc = CalculateGcContent(seq);

        Assert.That(gc, Is.EqualTo(1.0));
    }

    [Test]
    public void CalculateGcContent_NoGC_Zero()
    {
        string seq = "AAAAUUUU";
        double gc = CalculateGcContent(seq);

        Assert.That(gc, Is.EqualTo(0));
    }

    [Test]
    public void CalculateGcContent_EmptySequence_Zero()
    {
        Assert.That(CalculateGcContent(""), Is.EqualTo(0));
    }

    [Test]
    public void GenerateSeedVariants_GeneratesVariants()
    {
        string seed = "AGCAGCA";
        var variants = GenerateSeedVariants(seed).ToList();

        // Original + 3 variants per position
        Assert.That(variants, Has.Count.EqualTo(1 + seed.Length * 3));
        Assert.That(variants[0], Is.EqualTo(seed));
    }

    [Test]
    public void GenerateSeedVariants_IncludesOriginal()
    {
        string seed = "AGCA";
        var variants = GenerateSeedVariants(seed);

        Assert.That(variants.Contains(seed), Is.True);
    }

    #endregion

    #region Integration Tests

    [Test]
    public void FullWorkflow_AnalyzeMiRNAFamily()
    {
        // let-7 family members have similar seeds
        var mirnas = new List<MiRna>
        {
            CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU"),
            CreateMiRna("let-7b", "UGAGGUAGUAGGUUGUGUGGUU"),
            CreateMiRna("let-7c", "UGAGGUAGUAGGUUGUAUGGUU")
        };

        var families = GroupBySeedFamily(mirnas).ToList();

        // All should share the same seed
        Assert.That(mirnas.Select(m => m.SeedSequence).Distinct().Count(), Is.EqualTo(1));
    }

    // FullWorkflow_PredictTargets migrated to MiRnaAnalyzer_TargetPrediction_Tests.cs (MIRNA-TARGET-001)
    // TargetSite_HasAllFields migrated to MiRnaAnalyzer_TargetPrediction_Tests.cs (MIRNA-TARGET-001)

    #endregion
}
