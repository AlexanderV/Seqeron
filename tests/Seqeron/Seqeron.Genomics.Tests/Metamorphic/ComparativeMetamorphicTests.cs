using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Metamorphic;

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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: COMPGEN-COMPARE-001 — two-genome comparison (Comparative).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 133.
    //
    // API under test (ComparativeGenomics.CompareGenomes):
    //   Partitions genes into the core (shared, reciprocal-best-hit) set and each genome's
    //   genome-specific set (Tettelin pan-genome model); ConservedGenes counts the shared genes.
    //   Similarity is k-mer Jaccard, so identical gene sequences are reciprocal best hits.
    //
    // Relations (derived from the RBH/pan-genome definition, NOT from output):
    //   • MON  (more shared genes ⇒ higher similarity): adding an orthologous gene pair adds one
    //          reciprocal best hit, so ConservedGenes (the similarity measure) increases.
    //   • SYM  (order independent): RBH is an unordered matching, so swapping the two genomes keeps
    //          ConservedGenes and OverallSynteny identical and merely swaps the two genome-specific
    //          counts.
    // ───────────────────────────────────────────────────────────────────────────

    // Each ortholog group is encoded as a distinct repeated letter, so two genes share all their
    // 5-mers (identity 1.0 → reciprocal best hit) iff they carry the same group id, and share none
    // otherwise (no qualifying hit → genome-specific).
    private static ComparativeGenomics.Gene GeneOf(string genomeId, int position, int groupId) =>
        new($"{genomeId}_g{position}", genomeId, position, position, '+', new string((char)('A' + groupId), 10));

    // A genome whose genes carry the given ortholog-group ids, in order.
    private static List<ComparativeGenomics.Gene> GenomeOf(string genomeId, params int[] groupIds) =>
        groupIds.Select((grp, pos) => GeneOf(genomeId, pos, grp)).ToList();

    #region COMPGEN-COMPARE-001 MON — more shared genes raise the conserved-gene count

    [Test]
    [Description("MON: each added orthologous gene pair contributes one reciprocal best hit, so ConservedGenes grows as more shared genes are present.")]
    public void CompareGenomes_MoreSharedGenes_HigherConservedCount()
    {
        int previous = int.MinValue;
        foreach (int shared in new[] { 1, 2, 4, 7 })
        {
            int[] groups = Enumerable.Range(0, shared).ToArray();
            var result = ComparativeGenomics.CompareGenomes(GenomeOf("A", groups), GenomeOf("B", groups));

            result.ConservedGenes.Should().Be(shared, because: $"the two genomes share exactly {shared} ortholog groups");
            result.ConservedGenes.Should().BeGreaterThan(previous, because: "adding a shared gene increases the conserved-gene similarity");
            previous = result.ConservedGenes;
        }
    }

    #endregion

    #region COMPGEN-COMPARE-001 SYM — swapping the genomes only swaps the genome-specific counts

    [Test]
    [Description("SYM: RBH is an unordered matching, so CompareGenomes(A,B) and CompareGenomes(B,A) report the same ConservedGenes and OverallSynteny and merely exchange the two genome-specific counts.")]
    public void CompareGenomes_SwapGenomes_ConsistentResult()
    {
        // Three shared groups (0,1,2); A has two specific genes (20,21), B has one (24).
        var a = GenomeOf("A", 0, 1, 2, 20, 21);
        var b = GenomeOf("B", 0, 1, 2, 24);

        var ab = ComparativeGenomics.CompareGenomes(a, b);
        var ba = ComparativeGenomics.CompareGenomes(b, a);

        ab.ConservedGenes.Should().Be(ba.ConservedGenes, because: "the reciprocal best-hit matching is the same regardless of argument order");
        ab.Orthologs.Count.Should().Be(ba.Orthologs.Count, because: "ortholog pairs are unordered");
        ab.OverallSynteny.Should().BeApproximately(ba.OverallSynteny, 1e-12, because: "synteny is normalised by the smaller genome, which is order-independent");
        ab.GenomeSpecificGenes1.Should().Be(ba.GenomeSpecificGenes2, because: "genome A's specific genes are reported as genome-1-specific in (A,B) and genome-2-specific in (B,A)");
        ab.GenomeSpecificGenes2.Should().Be(ba.GenomeSpecificGenes1, because: "genome B's specific genes swap roles symmetrically");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: COMPGEN-DOTPLOT-001 — word-match dot plot (Comparative).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 134.
    //
    // API under test (ComparativeGenomics.GenerateDotPlot):
    //   Reports every (x, y) where the length-w word at sequence1[x] exactly equals the word at
    //   sequence2[y] (EMBOSS dottup; Gibbs & McIntyre 1970). x ranges over sequence1, y over sequence2.
    //
    // Relations (derived from exact word matching, NOT from output):
    //   • INV  (reverse-complement maps the diagonal): the dot plot detects only FORWARD exact
    //          matches, so reverse-complementing BOTH axes reflects the whole plot through its
    //          centre — dot (x,y) ↔ (L1−w−x, L2−w−y) — because revcomp(s1)[x..]=revcomp(s2)[y..]
    //          iff the mirror-position forward words are equal. The main diagonal maps onto the
    //          reflected diagonal. (A forward→anti-diagonal map under reverse-complementing a single
    //          axis does not hold for arbitrary sequences, since forward word matching cannot detect
    //          reverse-complement matches as a coordinate transform.)
    //   • SHIFT (prepend flank shifts dots): prepending an f-base flank to sequence1 shifts every
    //          existing match's x-coordinate by f, preserving all original dots at (x+f, y).
    // ───────────────────────────────────────────────────────────────────────────

    private static string RevComp(string dna) =>
        new(dna.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c }).ToArray());

    private static HashSet<(int x, int y)> DotSet(string s1, string s2, int wordSize) =>
        ComparativeGenomics.GenerateDotPlot(s1, s2, wordSize).ToHashSet();

    #region COMPGEN-DOTPLOT-001 INV — reverse-complementing both axes reflects the plot

    [Test]
    [Description("INV: forward word matching only detects forward matches, so reverse-complementing both sequences reflects the dot plot through its centre, mapping every dot (x,y) to (L1−w−x, L2−w−y).")]
    public void DotPlot_ReverseComplementBothAxes_ReflectsPlot()
    {
        const string s1 = "ACGTACGTTTGCA";
        const string s2 = "ACGTACGTAAGGC";
        const int w = 4;

        var original = DotSet(s1, s2, w);
        var reflected = DotSet(RevComp(s1), RevComp(s2), w);

        var expected = original.Select(d => (s1.Length - w - d.x, s2.Length - w - d.y)).ToHashSet();
        reflected.Should().BeEquivalentTo(expected,
            because: "revcomp(s1)[x..]=revcomp(s2)[y..] iff the mirror-position forward words match, so the plot is point-reflected through its centre");
    }

    #endregion

    #region COMPGEN-DOTPLOT-001 SHIFT — prepending a flank shifts the dots along x

    [Test]
    [Description("SHIFT: prepending an f-base flank to sequence1 shifts every existing match's x-coordinate by f, so all original dots reappear at (x+f, y).")]
    public void DotPlot_PrependFlankToSequence1_ShiftsDots()
    {
        const string s1 = "ACGTACGTTTGCA";
        const string s2 = "ACGTACGTAAGGC";
        const int w = 4;
        var original = DotSet(s1, s2, w);

        foreach (var flank in new[] { "TT", "GGGCCC" })
        {
            var shifted = DotSet(flank + s1, s2, w);
            var expected = original.Select(d => (d.x + flank.Length, d.y));
            shifted.Should().Contain(expected,
                because: $"words wholly inside the original sequence keep matching the same sequence2 words, only their x-coordinate moves by the {flank.Length}-base flank");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: COMPGEN-ORTHO-001 — ortholog detection (reciprocal best hits) (Comparative).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 135.
    //
    // API under test (ComparativeGenomics.FindOrthologs → FindReciprocalBestHits):
    //   Two genes are orthologs iff each is the other's best hit across the two genomes
    //   (Moreno-Hagelsieb & Latimer 2008). Similarity is k-mer Jaccard, so identical gene
    //   sequences are reciprocal best hits.
    //
    // Relations (derived from the reciprocal-best-hit definition, NOT from output):
    //   • SYM  (ortholog relation symmetric): RBH is reciprocal by construction, so a↔b is an
    //          ortholog pair in FindOrthologs(A,B) iff b↔a is one in FindOrthologs(B,A) — the
    //          unordered pair set is identical.
    //   • INV  (gene order independent): best-hit ties are broken deterministically, so permuting
    //          the gene order within each genome leaves the ortholog pair set unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    private static string UnorderedPair(string a, string b) =>
        string.CompareOrdinal(a, b) < 0 ? $"{a}|{b}" : $"{b}|{a}";

    private static HashSet<string> OrthologPairs(
        IReadOnlyList<ComparativeGenomics.Gene> g1, IReadOnlyList<ComparativeGenomics.Gene> g2) =>
        ComparativeGenomics.FindOrthologs(g1, g2).Select(o => UnorderedPair(o.Gene1Id, o.Gene2Id)).ToHashSet();

    #region COMPGEN-ORTHO-001 SYM — the ortholog relation is symmetric

    [Test]
    [Description("SYM: reciprocal best hits are reciprocal by definition, so swapping the two genomes yields the same set of (unordered) ortholog pairs.")]
    public void Orthologs_SwapGenomes_SamePairs()
    {
        var a = GenomeOf("A", 0, 1, 2, 20);
        var b = GenomeOf("B", 0, 1, 2, 24);

        OrthologPairs(b, a).Should().BeEquivalentTo(OrthologPairs(a, b),
            because: "a is b's best hit iff b is a's best hit, so the ortholog relation is symmetric");
    }

    #endregion

    #region COMPGEN-ORTHO-001 INV — ortholog set is independent of gene order

    [Test]
    [Description("INV: best-hit selection depends only on sequence similarity (ties broken deterministically), so reordering the genes within each genome leaves the ortholog pair set unchanged.")]
    public void Orthologs_PermuteGeneOrder_SamePairs()
    {
        var a = GenomeOf("A", 0, 1, 2, 20);
        var b = GenomeOf("B", 0, 1, 2, 24);

        var aShuffled = ((IEnumerable<ComparativeGenomics.Gene>)a).Reverse().ToList();
        var bShuffled = ((IEnumerable<ComparativeGenomics.Gene>)b).Reverse().ToList();

        OrthologPairs(aShuffled, bShuffled).Should().BeEquivalentTo(OrthologPairs(a, b),
            because: "orthology is a property of the gene sequences, not their chromosomal order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: COMPGEN-RBH-001 — reciprocal best hits with hit metrics (Comparative).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 136.
    //
    // API under test (ComparativeGenomics.FindReciprocalBestHits):
    //   The dedicated RBH entry point, returning each pair with its identity, coverage and
    //   alignment length. Distinct from COMPGEN-ORTHO-001 in also checking the reported metrics.
    //
    // Relations (derived from the symmetric similarity score, NOT from output):
    //   • SYM  (RBH symmetric): the k-mer similarity is symmetric, so swapping the genomes reverses
    //          each pair AND preserves its identity/coverage.
    //   • INV  (input order independent): the matching is unique with deterministic tie-breaking, so
    //          permuting the gene order in each genome preserves the pairs and their metrics.
    // ───────────────────────────────────────────────────────────────────────────

    // Maps each unordered gene pair to its rounded (identity, coverage) hit metrics.
    private static Dictionary<string, (double Id, double Cov)> RbhMetrics(
        IReadOnlyList<ComparativeGenomics.Gene> g1, IReadOnlyList<ComparativeGenomics.Gene> g2) =>
        ComparativeGenomics.FindReciprocalBestHits(g1, g2)
            .ToDictionary(o => UnorderedPair(o.Gene1Id, o.Gene2Id),
                         o => (System.Math.Round(o.Identity, 9), System.Math.Round(o.Coverage, 9)));

    #region COMPGEN-RBH-001 SYM — RBH is symmetric, metrics included

    [Test]
    [Description("SYM: the similarity score is symmetric, so FindReciprocalBestHits(A,B) and (B,A) report the same pairs (reversed) with identical identity and coverage.")]
    public void Rbh_SwapGenomes_SamePairsAndMetrics()
    {
        var a = GenomeOf("A", 0, 1, 2, 20);
        var b = GenomeOf("B", 0, 1, 2, 24);

        RbhMetrics(b, a).Should().BeEquivalentTo(RbhMetrics(a, b),
            because: "k-mer Jaccard identity and coverage do not depend on which sequence is the query");
    }

    #endregion

    #region COMPGEN-RBH-001 INV — RBH is independent of input order

    [Test]
    [Description("INV: with deterministic tie-breaking the reciprocal-best-hit matching depends only on the gene sequences, so permuting the input order preserves both the pairs and their metrics.")]
    public void Rbh_PermuteInputOrder_SamePairsAndMetrics()
    {
        var a = GenomeOf("A", 0, 1, 2, 20);
        var b = GenomeOf("B", 0, 1, 2, 24);

        var aShuffled = ((IEnumerable<ComparativeGenomics.Gene>)a).Reverse().ToList();
        var bShuffled = ((IEnumerable<ComparativeGenomics.Gene>)b).Reverse().ToList();

        RbhMetrics(aShuffled, bShuffled).Should().BeEquivalentTo(RbhMetrics(a, b),
            because: "the RBH matching is a function of the sequence set, not the input order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: COMPGEN-REARR-001 — rearrangement detection (Comparative).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 137.
    //
    // API under test (ComparativeGenomics.DetectRearrangements):
    //   Reads the genome-1 orthologs (in order) as a signed permutation of their genome-2 ranks and
    //   reports each breakpoint of the extended permutation (Bafna & Pevzner 1998).
    //
    // Relations (derived from the breakpoint model, NOT from output):
    //   • INV  (identity ⇒ no rearrangements): when the ortholog order and strands are identical the
    //          relative permutation is the identity, which has no breakpoints, so no events are reported.
    //   • SYM  ((A,B) consistent with (B,A)): the breakpoint distance of a permutation equals that of
    //          its inverse, so swapping the genomes (and the ortholog map) reports the same number of
    //          breakpoints.
    // ───────────────────────────────────────────────────────────────────────────

    // Builds a genome whose genes carry the given ortholog-group ids in list order (all '+' strand).
    private static List<ComparativeGenomics.Gene> OrderedGenome(string genomeId, params int[] groupOrder) =>
        groupOrder.Select((grp, pos) => new ComparativeGenomics.Gene($"{genomeId}{grp}", genomeId, pos, pos, '+'))
                  .ToList();

    // Ortholog map linking each group's gene in genome A to the same group's gene in genome B.
    private static Dictionary<string, string> GroupMap(string fromGenome, string toGenome, int groupCount) =>
        Enumerable.Range(0, groupCount).ToDictionary(g => $"{fromGenome}{g}", g => $"{toGenome}{g}");

    #region COMPGEN-REARR-001 INV — identical order yields no rearrangements

    [Test]
    [Description("INV: when genome 2 has the same ortholog order and strands as genome 1 the relative permutation is the identity, which has no breakpoints, so no rearrangements are reported.")]
    public void Rearrangements_IdenticalOrder_None()
    {
        var g1 = OrderedGenome("A", 0, 1, 2, 3, 4);
        var g2 = OrderedGenome("B", 0, 1, 2, 3, 4);

        ComparativeGenomics.DetectRearrangements(g1, g2, GroupMap("A", "B", 5))
            .Should().BeEmpty(because: "an identity permutation has no breakpoints");
    }

    #endregion

    #region COMPGEN-REARR-001 SYM — swapping the genomes reports the same breakpoint count

    [Test]
    [Description("SYM: the breakpoint distance of a permutation equals that of its inverse, so DetectRearrangements(A,B) and (B,A) report the same number of breakpoints.")]
    public void Rearrangements_SwapGenomes_SameBreakpointCount()
    {
        int[] permuted = { 2, 0, 4, 1, 3 };
        var g1 = OrderedGenome("A", 0, 1, 2, 3, 4);
        var g2 = OrderedGenome("B", permuted);

        int forward = ComparativeGenomics.DetectRearrangements(g1, g2, GroupMap("A", "B", 5)).Count();
        int reverse = ComparativeGenomics.DetectRearrangements(g2, g1, GroupMap("B", "A", 5)).Count();

        forward.Should().BeGreaterThan(0, because: "the chosen order is a non-trivial permutation, so breakpoints exist");
        reverse.Should().Be(forward, because: "a permutation and its inverse share the same breakpoint distance");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: COMPGEN-REVERSAL-001 — reversal distance (Comparative).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 138.
    //
    // API under test (ComparativeGenomics.CalculateReversalDistance):
    //   The unsigned breakpoint lower bound ⌈b/2⌉ on the reversal (sorting-by-reversals) distance
    //   between two gene orders (Bafna & Pevzner 1998).
    //
    // Relations (derived from the breakpoint bound, NOT from output):
    //   • INV  (identical permutation ⇒ 0): identical orders have no breakpoints, so the distance is 0.
    //   • SYM  (symmetric): the breakpoint count is invariant under inverting the relative
    //          permutation, so d(p1,p2)=d(p2,p1).
    //   • MON  (more reversals applied ⇒ ≥ distance): each reversal removes at most two breakpoints, so
    //          the returned lower bound never exceeds the number of reversals actually applied — the
    //          applied count is an upper bound (≥) on the distance.
    // ───────────────────────────────────────────────────────────────────────────

    // Reverses the inclusive sub-range [i, j] of a copy of the permutation (one reversal operation).
    private static int[] ApplyReversal(int[] perm, int i, int j)
    {
        var copy = (int[])perm.Clone();
        System.Array.Reverse(copy, i, j - i + 1);
        return copy;
    }

    #region COMPGEN-REVERSAL-001 INV — identical permutations are distance 0

    [Test]
    [Description("INV: an order compared with itself has no breakpoints, so its reversal distance is 0.")]
    public void ReversalDistance_IdenticalPermutation_Zero()
    {
        foreach (var p in new[] { new[] { 0, 1, 2, 3, 4 }, new[] { 3, 1, 4, 0, 2 } })
            ComparativeGenomics.CalculateReversalDistance(p, p).Should().Be(0,
                because: "identical orders share every adjacency, so there are no breakpoints");
    }

    #endregion

    #region COMPGEN-REVERSAL-001 SYM — the distance is symmetric

    [Test]
    [Description("SYM: the breakpoint count is invariant under inverting the relative permutation, so d(p1,p2)=d(p2,p1).")]
    public void ReversalDistance_Symmetric()
    {
        var p1 = new[] { 2, 0, 4, 1, 3 };
        var p2 = new[] { 0, 1, 2, 3, 4 };
        var p3 = new[] { 4, 3, 2, 1, 0 };

        ComparativeGenomics.CalculateReversalDistance(p2, p1)
            .Should().Be(ComparativeGenomics.CalculateReversalDistance(p1, p2), because: "reversal distance is symmetric");
        ComparativeGenomics.CalculateReversalDistance(p3, p1)
            .Should().Be(ComparativeGenomics.CalculateReversalDistance(p1, p3), because: "reversal distance is symmetric");
    }

    #endregion

    #region COMPGEN-REVERSAL-001 MON — the distance never exceeds the reversals applied

    [Test]
    [Description("MON: each reversal removes at most two breakpoints, so the breakpoint lower bound is at most the number of reversals applied to reach the permutation from the identity.")]
    public void ReversalDistance_AtMostReversalsApplied()
    {
        int[] identity = { 0, 1, 2, 3, 4, 5, 6, 7 };

        // A fixed sequence of reversal operations; after applying the first k of them, the distance
        // back to the identity must be ≤ k.
        var operations = new (int i, int j)[] { (1, 4), (0, 2), (3, 7), (2, 6), (1, 5) };

        int[] current = identity;
        for (int k = 0; k < operations.Length; k++)
        {
            current = ApplyReversal(current, operations[k].i, operations[k].j);
            int distance = ComparativeGenomics.CalculateReversalDistance(current, identity);
            distance.Should().BeLessThanOrEqualTo(k + 1,
                because: $"reaching this order took {k + 1} reversals, an upper bound on the true (and hence the lower-bound) distance");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: COMPGEN-SYNTENY-001 — syntenic block detection (Comparative).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 139.
    //
    // API under test (ComparativeGenomics.FindSyntenicBlocks):
    //   MCScanX-style collinear chaining of orthologous anchors (Wang et al. 2012). A chain is a
    //   block when its score ≥ 250 and it has ≥ minAnchors anchors (the block-size threshold).
    //
    // Relations (derived from the chaining/report rule, NOT from output):
    //   • MON  (lower minBlockSize ⇒ superset): lowering minAnchors keeps every larger block and can
    //          only admit additional smaller ones, so the block set grows.
    //   • INV  (reverse preserves block count): reverse-complementing both genomes reverses the gene
    //          order on both axes, which keeps each collinear chain collinear, so the block count is
    //          unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    // A two-block scenario: a forward block of 6 anchors (genome-2 positions 0..5) and a block of 5
    // anchors (genome-2 positions 40..44), separated by a genome-2 gap larger than maxGap so the
    // greedy chainer cannot merge them.
    private static (List<ComparativeGenomics.Gene> G1, List<ComparativeGenomics.Gene> G2, Dictionary<string, string> Map)
        BuildSyntenyScenario()
    {
        var g1 = new List<ComparativeGenomics.Gene>();
        for (int i = 0; i <= 10; i++)
            g1.Add(new ComparativeGenomics.Gene($"A{i}", "A", i, i, '+'));

        var g2 = new List<ComparativeGenomics.Gene>();
        int pos = 0;
        for (int i = 0; i <= 5; i++) g2.Add(new ComparativeGenomics.Gene($"B{i}", "B", pos, pos++, '+')); // B0..B5 at 0..5
        for (int f = 0; f < 34; f++) g2.Add(new ComparativeGenomics.Gene($"F{f}", "B", pos, pos++, '+')); // filler 6..39
        for (int i = 6; i <= 10; i++) g2.Add(new ComparativeGenomics.Gene($"B{i}", "B", pos, pos++, '+')); // B6..B10 at 40..44

        var map = Enumerable.Range(0, 11).ToDictionary(i => $"A{i}", i => $"B{i}");
        return (g1, g2, map);
    }

    private static HashSet<string> BlockKeys(
        IReadOnlyList<ComparativeGenomics.Gene> g1, IReadOnlyList<ComparativeGenomics.Gene> g2,
        Dictionary<string, string> map, int minAnchors) =>
        ComparativeGenomics.FindSyntenicBlocks(g1, g2, map, minAnchors)
            .Select(b => $"{b.Start1}-{b.End1}:{b.Start2}-{b.End2}:{b.GeneCount}")
            .ToHashSet();

    #region COMPGEN-SYNTENY-001 MON — lowering the block-size threshold yields a superset

    [Test]
    [Description("MON: lowering minAnchors keeps every larger syntenic block and admits smaller ones, so the block set at a lower threshold is a superset of the set at a higher threshold.")]
    public void SyntenicBlocks_LowerMinBlockSize_Superset()
    {
        var (g1, g2, map) = BuildSyntenyScenario();

        var size5 = BlockKeys(g1, g2, map, minAnchors: 5);
        var size6 = BlockKeys(g1, g2, map, minAnchors: 6);

        size6.IsSubsetOf(size5).Should().BeTrue(because: "every block passing the stricter size-6 cut-off also passes size-5");
        size5.Count.Should().BeGreaterThan(size6.Count, because: "the 5-anchor block is admitted at threshold 5 but not at threshold 6");
    }

    #endregion

    #region COMPGEN-SYNTENY-001 INV — reversing both genomes preserves the block count

    [Test]
    [Description("INV: reverse-complementing both genomes reverses the gene order on both axes, keeping each collinear chain collinear, so the number of syntenic blocks is unchanged.")]
    public void SyntenicBlocks_ReverseBothGenomes_SameBlockCount()
    {
        var (g1, g2, map) = BuildSyntenyScenario();

        int original = ComparativeGenomics.FindSyntenicBlocks(g1, g2, map, minAnchors: 5).Count();

        var g1Rev = ((IEnumerable<ComparativeGenomics.Gene>)g1).Reverse().ToList();
        var g2Rev = ((IEnumerable<ComparativeGenomics.Gene>)g2).Reverse().ToList();
        int reversed = ComparativeGenomics.FindSyntenicBlocks(g1Rev, g2Rev, map, minAnchors: 5).Count();

        original.Should().Be(2, because: "the scenario is constructed with exactly two collinear blocks");
        reversed.Should().Be(original, because: "reversing both axes keeps each chain collinear, preserving the block count");
    }

    #endregion
}
