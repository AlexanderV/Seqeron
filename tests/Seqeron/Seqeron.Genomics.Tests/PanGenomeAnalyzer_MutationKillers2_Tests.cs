using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using static Seqeron.Genomics.Metagenomics.PanGenomeAnalyzer;

namespace Seqeron.Genomics.Tests;

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
    public void ConstructPanGenome_DecayingNoveltyCurve_IsClosed()
    {
        // 3 genomes; novelty per added genome decays 3 → 1, giving a steep negative log-log slope
        // (alpha ≈ 2.7 > 1.0) ⇒ Closed. A broken regression-slope mutant would flip this to Open.
        var g = Genomes(
            ("g1", new[] { ("g1a", "AAAAAAAAAA"), ("g1b", "CCCCCCCCCC") }),
            ("g2", new[] { ("g2a", "AAAAAAAAAA"), ("g2c", "GGGGGGGGGG"), ("g2d", "TTTTTTTTTT"), ("g2e", "ACACACACAC") }),
            ("g3", new[] { ("g3a", "AAAAAAAAAA"), ("g3f", "GTGTGTGTGT") }));
        var r = ConstructPanGenome(g, identityThreshold: 0.9, coreFraction: 0.99);
        Assert.That(r.Statistics.Type, Is.EqualTo(PanGenomeType.Closed));
    }

    [Test]
    public void ConstructPanGenome_SustainedNovelty_IsOpen()
    {
        // Each added genome contributes the same large number of new clusters (4, 4) ⇒ flat
        // log-log curve, alpha ≈ 0 < 1.0 ⇒ Open. Requires BOTH curve points (k = 2 and k = 3):
        // a 'k > 2' or 'logK.Count <= 2' mutant drops to one point and falls back to Closed.
        var g = Genomes(
            ("g1", new[] { ("a", "AAAAAAAAAA") }),
            ("g2", new[] { ("a", "AAAAAAAAAA"), ("b", "CCCCCCCCCC"), ("c", "GGGGGGGGGG"), ("d", "TTTTTTTTTT"), ("e", "ACACACACAC") }),
            ("g3", new[] { ("a", "AAAAAAAAAA"), ("f", "GTGTGTGTGT"), ("h", "CGCGCGCGCG"), ("i", "ATATATATAT"), ("j", "TGTGTGTGTG") }));
        var r = ConstructPanGenome(g, identityThreshold: 0.9, coreFraction: 0.99);
        Assert.That(r.Statistics.Type, Is.EqualTo(PanGenomeType.Open));
    }
}
