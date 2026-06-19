using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Comparative-genomics area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-ANI-001 — Average Nucleotide Identity (ANIb, Goris et al. 2007).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 131.
///
/// API under test (ComparativeGenomics.CalculateANI):
///   The query genome is cut into consecutive fragments of fragmentLength; each fragment is
///   ungapped-best-placed in the reference; ANI is the mean identity of the qualifying fragments,
///   reported as a fraction in [0, 1] (1.0 ≡ 100 %).
///
/// Relations (derived from the ANIb definition, NOT from output):
///   • INV  (self-identity is maximal): every fragment of a genome is a perfect substring of
///          itself, so ANI(A,A)=1.0 regardless of fragmentation.
///   • MON  (more mutations ⇒ lower ANI): each substituted base lowers the matching count of its
///          fragment, so adding substitutions monotonically lowers ANI.
///   • SYM  (reciprocal symmetry under equal-length, full-length alignment): ANIb is asymmetric in
///          general because only the query is fragmented (Goris 2007 / pyani). For equal-length
///          genomes aligned as a single full-length fragment, the only placement is the diagonal
///          (offset 0), so both directions count the same matched positions and ANI(A,B)=ANI(B,A).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class ComparativeMetamorphicTests
{
    // A fixed pseudo-genome over {A,C,G,T}; its own .Length is used as the full-length fragment size.
    private const string Genome = "ACGTTGCAACGTGGATCCGTACGATCGATTACAGGCATTAGCATCGTA";

    // Substitutes the first <paramref name="count"/> bases with a guaranteed-different base,
    // producing a copy that differs from the original at exactly <paramref name="count"/> positions.
    private static string Substitute(string seq, int count)
    {
        char[] arr = seq.ToCharArray();
        for (int i = 0; i < count; i++)
            arr[i] = arr[i] == 'A' ? 'C' : 'A';
        return new string(arr);
    }

    #region COMPGEN-ANI-001 INV — identical genomes give the maximal ANI of 1.0

    [Test]
    [Description("INV: every fragment of a genome is a perfect substring of itself, so ANI(A,A)=1.0 (the [0,1] form of 100 %), here exercised over multiple fragments.")]
    public void Ani_IdenticalGenomes_IsOne()
    {
        // fragmentLength 12 over the length-48 genome → 4 fragments, each a perfect self-match.
        ComparativeGenomics.CalculateANI(Genome, Genome, fragmentLength: 12).Should().Be(1.0,
            because: "each consecutive fragment occurs verbatim in the genome, so every per-fragment identity is 1.0 and so is their mean");
    }

    #endregion

    #region COMPGEN-ANI-001 MON — more substitutions lower ANI

    [Test]
    [Description("MON: each substituted base reduces the matched-position count, so introducing progressively more substitutions monotonically lowers ANI.")]
    public void Ani_MoreSubstitutions_LowerAni()
    {
        double previous = double.MaxValue;
        foreach (int mutations in new[] { 0, 1, 3, 6 })
        {
            // Full-length single fragment: ANI = (L − mutations)/L exactly, so the trend is strict.
            double ani = ComparativeGenomics.CalculateANI(Genome, Substitute(Genome, mutations), fragmentLength: Genome.Length);
            ani.Should().BeLessThan(previous, because: $"{mutations} substitutions leave fewer matched bases than {previous} did");
            previous = ani;
        }
    }

    #endregion

    #region COMPGEN-ANI-001 SYM — reciprocal symmetry for equal-length full-length alignment

    [Test]
    [Description("SYM: for equal-length genomes aligned as a single full-length fragment the only placement is the diagonal, so the two reciprocal directions count the same matched positions and ANI(A,B)=ANI(B,A).")]
    public void Ani_EqualLengthFullFragment_Symmetric()
    {
        foreach (int mutations in new[] { 2, 4, 8 })
        {
            string a = Genome;
            string b = Substitute(Genome, mutations);

            ComparativeGenomics.CalculateANI(a, b, fragmentLength: Genome.Length)
                .Should().Be(ComparativeGenomics.CalculateANI(b, a, fragmentLength: Genome.Length),
                    because: "matched-position count is symmetric and, at equal length, only the offset-0 placement exists");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: COMPGEN-CLUSTER-001 — conserved gene clusters as common intervals (Comparative).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 132.
    //
    // API under test (ComparativeGenomics.FindConservedClusters):
    //   A conserved cluster is a COMMON INTERVAL of the ortholog-group permutations — a set of
    //   group labels that occupies a contiguous window in EVERY genome (Uno & Yagiura 2000;
    //   Heber & Stoye 2001; Bui-Xuan et al. 2013). The size cut-off is minClusterSize.
    //
    // Relations (derived from the common-interval definition, NOT from output):
    //   • MON  (lower size threshold ⇒ superset): the checklist's "identity threshold" maps to this
    //          model's only monotone parameter, minClusterSize — lowering it can only admit more
    //          (smaller) clusters while keeping every larger one, so the cluster set grows.
    //   • INV  (genome order independent): a common interval is an interval of every genome by
    //          definition, so permuting the genome list cannot change the set of common intervals.
    // ───────────────────────────────────────────────────────────────────────────

    // Three genomes (ortholog-group order) with known common intervals {g1,g2}, {g4,g5},
    // {g1,g2,g3} and the trivial whole set {g1..g5}:
    //   A: g1 g2 g3 g4 g5   B: g3 g2 g1 g5 g4   C: g2 g1 g3 g4 g5
    private static readonly string[][] GenomeGroupOrders =
    {
        new[] { "g1", "g2", "g3", "g4", "g5" },
        new[] { "g3", "g2", "g1", "g5", "g4" },
        new[] { "g2", "g1", "g3", "g4", "g5" },
    };

    // Builds the genome gene-lists and the gene-id → ortholog-group map for a given genome ordering.
    private static (List<IReadOnlyList<ComparativeGenomics.Gene>> Genomes, Dictionary<string, string> Groups)
        BuildScenario(IReadOnlyList<string[]> orders)
    {
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>>();
        var groups = new Dictionary<string, string>();
        for (int gi = 0; gi < orders.Count; gi++)
        {
            var genes = new List<ComparativeGenomics.Gene>();
            for (int p = 0; p < orders[gi].Length; p++)
            {
                string id = $"G{gi}_{p}";
                genes.Add(new ComparativeGenomics.Gene(id, $"genome{gi}", p, p, '+'));
                groups[id] = orders[gi][p];
            }
            genomes.Add(genes);
        }
        return (genomes, groups);
    }

    // Canonical, order-independent key set of the reported clusters.
    private static HashSet<string> ClusterKeys(IReadOnlyList<string[]> orders, int minClusterSize)
    {
        var (genomes, groups) = BuildScenario(orders);
        return ComparativeGenomics.FindConservedClusters(genomes, groups, minClusterSize)
            .Select(c => string.Join(",", c))
            .ToHashSet();
    }

    #region COMPGEN-CLUSTER-001 MON — lowering the size threshold yields a superset

    [Test]
    [Description("MON: lowering minClusterSize keeps every larger common interval and can only admit additional smaller ones, so the cluster set at a lower threshold is a superset of the set at a higher threshold.")]
    public void ConservedClusters_LowerSizeThreshold_Superset()
    {
        var size4 = ClusterKeys(GenomeGroupOrders, minClusterSize: 4);
        var size3 = ClusterKeys(GenomeGroupOrders, minClusterSize: 3);
        var size2 = ClusterKeys(GenomeGroupOrders, minClusterSize: 2);

        size4.IsSubsetOf(size3).Should().BeTrue(because: "every cluster surviving the stricter size-4 cut-off also survives size-3");
        size3.Count.Should().BeGreaterThan(size4.Count, because: "lowering the threshold from 4 to 3 admits the size-3 common interval {g1,g2,g3}");

        size3.IsSubsetOf(size2).Should().BeTrue(because: "every size-3 cluster also passes the size-2 cut-off");
        size2.Count.Should().BeGreaterThan(size3.Count, because: "lowering the threshold from 3 to 2 admits the size-2 common intervals {g1,g2} and {g4,g5}");
    }

    #endregion

    #region COMPGEN-CLUSTER-001 INV — the cluster set is independent of genome order

    [Test]
    [Description("INV: a common interval is an interval of every genome, so permuting the genome list leaves the set of conserved clusters unchanged.")]
    public void ConservedClusters_GenomeOrder_Invariant()
    {
        var original = GenomeGroupOrders;
        var permuted = new[] { GenomeGroupOrders[2], GenomeGroupOrders[0], GenomeGroupOrders[1] };

        foreach (int minSize in new[] { 2, 3 })
            ClusterKeys(permuted, minSize).Should().BeEquivalentTo(ClusterKeys(original, minSize),
                because: "the common-interval property is symmetric across the genome family, so input order is irrelevant");
    }

    #endregion
}
