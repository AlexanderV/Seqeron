using static Seqeron.Genomics.Metagenomics.PanGenomeAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Metagenomics;

/// <summary>
/// PANGEN-* mutation killers (batch 2): empty-sequence identity handling, parsimony-informative
/// site counting, and the open/closed classification driven by the Heaps log-log regression slope.
/// </summary>
[TestFixture]
public class PanGenomeAnalyzer_MutationKillers2_Tests
{
    private const double Tol = 1e-9;

    private static IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> Genomes(
        params (string Genome, (string GeneId, string Sequence)[] Genes)[] entries)
    {
        var d = new Dictionary<string, IReadOnlyList<(string, string)>>();
        foreach (var (g, genes) in entries)
            d[g] = genes.Select(x => (x.GeneId, x.Sequence)).ToList();
        return d;
    }

    [Test]
    public void ClusterGenes_EmptyAndNonEmpty_IdentityIsZeroNotNaN()
    {
        // identityThreshold 0 forces "" to join "AAA"; the pair identity is exactly 0
        // (one-empty rule). An '||'→'&&' mutant in the empty guard yields 0/0 = NaN instead.
        var g = Genomes(
            ("g1", new[] { ("x", "AAA") }),
            ("g2", new[] { ("y", "") }));
        var cluster = ClusterGenes(g, identityThreshold: 0.0).Single();
        Assert.That(cluster.AverageIdentity, Is.EqualTo(0.0).Within(Tol));
    }

    [Test]
    public void CountParsimonyInformativeSites_ExactCount()
    {
        // Col0: A,A,T,T → two states each in ≥2 rows ⇒ informative; Col1: T,T,T,T monomorphic;
        // Col2: G,G,G,G monomorphic ⇒ exactly 1 informative site.
        Assert.That(CountParsimonyInformativeSites(new[] { "ATG", "ATG", "TTG", "TTG" }), Is.EqualTo(1));

        // A singleton variant (A,A,A,T) is NOT parsimony-informative (only one state has ≥2).
        Assert.That(CountParsimonyInformativeSites(new[] { "AA", "AA", "AA", "AT" }), Is.EqualTo(0));

        // Too few sequences ⇒ 0.
        Assert.That(CountParsimonyInformativeSites(new[] { "AT" }), Is.EqualTo(0));
    }

    [Test]
    public void ConstructPanGenome_SaturatingAccessory_IsClosed()
    {
        // 3 genomes, clusters {core, x, y, z}; genome i lacks exactly one accessory cluster.
        // Under every ordering the 2nd genome adds 1 new cluster and the 3rd adds 0 (curve 1,0).
        // micropan heaps() objective, Python reference (exhaustive 3! orderings): K = 3.3402,
        // alpha = 2.0 (upper bound) > 1 => Closed (Tettelin 2008). A mutant that ignores the
        // fit (or inverts the alpha < 1 test) flips this to Open.
        var g = Genomes(
            ("g1", new[] { ("g1a", "AAAAAAAAAA"), ("g1y", "GGGGGGGGGG"), ("g1z", "TTTTTTTTTT") }),
            ("g2", new[] { ("g2a", "AAAAAAAAAA"), ("g2x", "CCCCCCCCCC"), ("g2z", "TTTTTTTTTT") }),
            ("g3", new[] { ("g3a", "AAAAAAAAAA"), ("g3x", "CCCCCCCCCC"), ("g3y", "GGGGGGGGGG") }));
        var r = ConstructPanGenome(g, identityThreshold: 0.9, coreFraction: 0.99);
        Assert.That(r.Statistics.Type, Is.EqualTo(PanGenomeType.Closed));
    }

    [Test]
    public void ConstructPanGenome_SustainedNovelty_IsOpen()
    {
        // g2 and g3 each carry 4 strain-specific clusters: under random orderings the expected
        // new-gene curve is flat, n(2) = n(3) = 8/3 (Python, exhaustive 3! orderings) ⇒ micropan
        // heaps() global optimum alpha = 0 < 1.0 ⇒ Open. N = 3 is the minimum for the fit: a
        // mutant raising the minimum genome count falls back to Closed.
        var g = Genomes(
            ("g1", new[] { ("a", "AAAAAAAAAA") }),
            ("g2", new[] { ("a", "AAAAAAAAAA"), ("b", "CCCCCCCCCC"), ("c", "GGGGGGGGGG"), ("d", "TTTTTTTTTT"), ("e", "ACACACACAC") }),
            ("g3", new[] { ("a", "AAAAAAAAAA"), ("f", "GTGTGTGTGT"), ("h", "CGCGCGCGCG"), ("i", "ATATATATAT"), ("j", "TGTGTGTGTG") }));
        var r = ConstructPanGenome(g, identityThreshold: 0.9, coreFraction: 0.99);
        Assert.That(r.Statistics.Type, Is.EqualTo(PanGenomeType.Open));
    }
}
