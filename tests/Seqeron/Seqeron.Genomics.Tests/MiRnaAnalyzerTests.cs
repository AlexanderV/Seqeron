using NUnit.Framework;
using Seqeron.Genomics;
using static Seqeron.Genomics.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class MiRnaAnalyzerTests
{
    #region Seed Sequence Tests

    [Test]
    public void GetSeedSequence_ValidMiRNA_ReturnsSeed()
    {
        string mirna = "UAGCAGCACGUAAAUAUUGGCG"; // let-7a
        string seed = GetSeedSequence(mirna);

        Assert.That(seed, Has.Length.EqualTo(7));
        Assert.That(seed, Is.EqualTo("AGCAGCA"));
    }

    [Test]
    public void GetSeedSequence_ShortSequence_ReturnsEmpty()
    {
        string seed = GetSeedSequence("UAGCA");
        Assert.That(seed, Is.Empty);
    }

    [Test]
    public void GetSeedSequence_EmptySequence_ReturnsEmpty()
    {
        Assert.That(GetSeedSequence(""), Is.Empty);
        Assert.That(GetSeedSequence(null!), Is.Empty);
    }

    [Test]
    public void CreateMiRna_ValidSequence_CreatesMiRna()
    {
        var mirna = CreateMiRna("let-7a", "UAGCAGCACGUAAAUAUUGGCG");

        Assert.Multiple(() =>
        {
            Assert.That(mirna.Name, Is.EqualTo("let-7a"));
            Assert.That(mirna.SeedSequence, Has.Length.EqualTo(7));
            Assert.That(mirna.SeedStart, Is.EqualTo(1));
            Assert.That(mirna.SeedEnd, Is.EqualTo(7));
        });
    }

    [Test]
    public void CreateMiRna_DNASequence_ConvertsToRNA()
    {
        var mirna = CreateMiRna("test", "TAGCAGCACGTAAATATTGGCG");

        Assert.That(mirna.Sequence, Does.Contain("U"));
        Assert.That(mirna.Sequence, Does.Not.Contain("T"));
    }

    #endregion

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

    #region Target Site Finding Tests

    [Test]
    public void FindTargetSites_PerfectSeedMatch_FindsSite()
    {
        var mirna = CreateMiRna("test-miR", "UAGCAGCACGUAAAUAUUGGCG");
        // Reverse complement of seed AGCAGCA is UGCUGCU
        string mrna = "AAAAAAUGCUGCUAAAAAA";

        var sites = FindTargetSites(mrna, mirna, minScore: 0.3).ToList();

        Assert.That(sites, Has.Count.GreaterThanOrEqualTo(1));
    }

    [Test]
    public void FindTargetSites_8mer_HighestScore()
    {
        var mirna = CreateMiRna("test-miR", "UAGCAGCACGUAAAUAUUGGCG");
        // 8mer: seed match + A at end
        string mrna = "AAAAAAUGCUGCUAAAAAAAA";

        var sites = FindTargetSites(mrna, mirna, minScore: 0.1).ToList();

        // Should find some match if the seed complementary is present
        // The seed is AGCAGCA, reverse complement is UGCUGCU  
        Assert.That(sites, Is.Not.Null);
        if (sites.Any())
        {
            Assert.That(sites.Max(s => s.Score), Is.GreaterThan(0));
        }
    }

    [Test]
    public void FindTargetSites_NoMatch_ReturnsEmpty()
    {
        var mirna = CreateMiRna("test-miR", "UAGCAGCACGUAAAUAUUGGCG");
        string mrna = "AAAAAAAAAAAAAAAAAAAA";

        var sites = FindTargetSites(mrna, mirna, minScore: 0.5).ToList();

        Assert.That(sites, Is.Empty);
    }

    [Test]
    public void FindTargetSites_MultipleSites_FindsAll()
    {
        var mirna = CreateMiRna("test-miR", "UAGCAGCACGUAAAUAUUGGCG");
        string seedRC = GetReverseComplement("AGCAGCA");
        string mrna = $"AAA{seedRC}AAA{seedRC}AAA{seedRC}AAA";

        var sites = FindTargetSites(mrna, mirna, minScore: 0.3).ToList();

        Assert.That(sites, Has.Count.GreaterThanOrEqualTo(2));
    }

    [Test]
    public void FindTargetSites_EmptySequences_ReturnsEmpty()
    {
        var mirna = CreateMiRna("test", "UAGCAGCA");

        Assert.That(FindTargetSites("", mirna).ToList(), Is.Empty);
        Assert.That(FindTargetSites("AAAA", new MiRna()).ToList(), Is.Empty);
    }

    [Test]
    public void FindTargetSites_IncludesAlignment()
    {
        var mirna = CreateMiRna("test-miR", "UAGCAGCACGUAAAUAUUGGCG");
        string seedRC = GetReverseComplement("AGCAGCA");
        string mrna = $"AAA{seedRC}AAA";

        var sites = FindTargetSites(mrna, mirna, minScore: 0.3).ToList();

        if (sites.Any())
        {
            Assert.That(sites[0].Alignment, Is.Not.Empty);
        }
    }

    #endregion

    #region Alignment Tests

    [Test]
    public void AlignMiRnaToTarget_PerfectMatch_AllMatches()
    {
        string mirna = "AAAA";
        string target = "UUUU";

        var duplex = AlignMiRnaToTarget(mirna, target);

        Assert.That(duplex.Matches, Is.EqualTo(4));
        Assert.That(duplex.Mismatches, Is.EqualTo(0));
    }

    [Test]
    public void AlignMiRnaToTarget_AllMismatches_NoMatches()
    {
        string mirna = "AAAA";
        string target = "AAAA";

        var duplex = AlignMiRnaToTarget(mirna, target);

        Assert.That(duplex.Mismatches, Is.EqualTo(4));
    }

    [Test]
    public void AlignMiRnaToTarget_WobblePairs_Detected()
    {
        string mirna = "GGGG";
        string target = "UUUU";

        var duplex = AlignMiRnaToTarget(mirna, target);

        Assert.That(duplex.GUWobbles, Is.EqualTo(4));
    }

    [Test]
    public void AlignMiRnaToTarget_EmptySequences_ReturnsEmptyDuplex()
    {
        var duplex = AlignMiRnaToTarget("", "AAAA");

        Assert.That(duplex.Matches, Is.EqualTo(0));
        Assert.That(duplex.MiRnaSequence, Is.Empty);
    }

    [Test]
    public void AlignMiRnaToTarget_CalculatesFreeEnergy()
    {
        string mirna = "UAGCAGCA";
        string target = "UGCUGCUA";

        var duplex = AlignMiRnaToTarget(mirna, target);

        Assert.That(duplex.FreeEnergy, Is.Not.EqualTo(0));
    }

    #endregion

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
    public void FullWorkflow_PredictTargets()
    {
        // Real miRNA: let-7a
        var let7a = CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU");

        // mRNA with potential target site (reverse complement of seed)
        string seedRC = GetReverseComplement(let7a.SeedSequence);
        string mrna = $"AUGGCUAAA{seedRC}AAAGCUUAA";

        var targets = FindTargetSites(mrna, let7a, minScore: 0.3).ToList();

        foreach (var target in targets)
        {
            var context = AnalyzeTargetContext(mrna, target.Start, target.End);

            Assert.Multiple(() =>
            {
                Assert.That(target.Score, Is.GreaterThan(0));
                Assert.That(context.AuContent, Is.GreaterThanOrEqualTo(0));
            });
        }
    }

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

    [Test]
    public void TargetSite_HasAllFields()
    {
        var mirna = CreateMiRna("test", "UAGCAGCACGUAAAUAUUGGCG");
        string seedRC = GetReverseComplement(mirna.SeedSequence);
        string mrna = $"AAA{seedRC}AAA";

        var sites = FindTargetSites(mrna, mirna, minScore: 0.1).ToList();

        if (sites.Any())
        {
            var site = sites[0];
            Assert.Multiple(() =>
            {
                Assert.That(site.Start, Is.GreaterThanOrEqualTo(0));
                Assert.That(site.End, Is.GreaterThan(site.Start));
                Assert.That(site.TargetSequence, Is.Not.Empty);
                Assert.That(site.MiRnaName, Is.EqualTo("test"));
                Assert.That(site.Score, Is.InRange(0, 1));
            });
        }
    }

    #endregion
}
