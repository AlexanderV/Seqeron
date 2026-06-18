using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for comparative genomics algorithms.
/// Verifies invariants drawn from the literature each algorithm implements.
///
/// Test Units: COMPGEN-ANI-001, COMPGEN-CLUSTER-001, COMPGEN-COMPARE-001, COMPGEN-DOTPLOT-001, COMPGEN-ORTHO-001, COMPGEN-RBH-001, COMPGEN-REARR-001, COMPGEN-REVERSAL-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("ComparativeGenomics")]
public class ComparativeGenomicsProperties
{
    /// <summary>Generates a DNA string of exactly <paramref name="len"/> bases over {A,C,G,T}.</summary>
    private static Arbitrary<string> DnaArbitrary(int len) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= len)
            .Select(a => new string(a, 0, len))
            .ToArbitrary();

    #region COMPGEN-ANI-001: R: ANI ∈ [0,1]; S: ANI(A,B)=ANI(B,A); I: ANI(A,A)=1; D: deterministic

    // ANIb (Goris et al. 2007, Int J Syst Evol Microbiol 57:81-91) cuts the query genome into
    // consecutive fragments, places each ungapped on the reference, keeps matches passing the
    // identity / alignable-fraction cut-offs, and averages their identities. The implementation
    // returns ANI as a *fraction in [0,1]* (Goris' percentage / 100), so the perfect-identity
    // value is 1.0 rather than 100.

    // A fragment length that fits inside the short sequences used by the property generators, so
    // that whole-sequence fragmentation actually produces qualifying fragments.
    private const int FragLen = 20;

    /// <summary>
    /// INV-1 (R): ANI is always a fraction in [0,1] for arbitrary query/reference DNA.
    /// Evidence: each fragment's contribution is min(matches/fragLen, 1.0) ∈ [0,1]; ANI is the
    /// mean of such values, hence bounded by [0,1] (Goris et al. 2007).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Ani_InUnitInterval()
    {
        return Prop.ForAll(DnaArbitrary(3 * FragLen), DnaArbitrary(3 * FragLen), (query, reference) =>
        {
            double ani = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: FragLen);
            return (ani >= 0.0 && ani <= 1.0)
                .Label($"ANI={ani} outside [0,1]");
        });
    }

    /// <summary>
    /// INV-2 (I): Self-identity. ANI(A,A) == 1.0 whenever A admits at least one fragment.
    /// Evidence: every fragment of A is a substring of A, so its best ungapped placement is exact
    /// (identity 1.0) and qualifies; the mean of 1.0 values is 1.0. This is the biological anchor
    /// that a genome is 100% identical to itself.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Ani_SelfIdentity_IsOne()
    {
        return Prop.ForAll(DnaArbitrary(3 * FragLen), seq =>
        {
            double ani = ComparativeGenomics.CalculateANI(seq, seq, fragmentLength: FragLen);
            return (Math.Abs(ani - 1.0) < 1e-12)
                .Label($"ANI(A,A)={ani}, expected 1.0");
        });
    }

    /// <summary>
    /// INV-3 (S): Symmetry. When both genomes have length equal to the fragment length there is a
    /// single fragment placed at the single possible offset, so the match count A↔B equals B↔A and
    /// ANI(A,B) == ANI(B,A) exactly. (In general ANIb is only approximately symmetric because the
    /// query alone is fragmented; this configuration isolates the exact-symmetry regime.)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Ani_IsSymmetric_ForEqualLengthGenomes()
    {
        return Prop.ForAll(DnaArbitrary(FragLen), DnaArbitrary(FragLen), (a, b) =>
        {
            double ab = ComparativeGenomics.CalculateANI(a, b, fragmentLength: FragLen, minIdentity: 0.0);
            double ba = ComparativeGenomics.CalculateANI(b, a, fragmentLength: FragLen, minIdentity: 0.0);
            return (Math.Abs(ab - ba) < 1e-12)
                .Label($"ANI(A,B)={ab} ≠ ANI(B,A)={ba}");
        });
    }

    /// <summary>
    /// INV-4 (M): Divergence does not increase identity. Mutating more positions of a copy of A
    /// yields ANI no higher than mutating fewer positions: ANI(A, more-mutated) ≤ ANI(A, less-mutated).
    /// Evidence: with one full-length fragment the qualifying identity is matches/length; extra
    /// substitutions can only reduce the match count (Goris et al. 2007 identity definition).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Ani_MoreMutations_DoNotIncreaseIdentity()
    {
        const string baseSeq = "ACGTACGTACGTACGTACGT"; // length 20 == FragLen → single fragment
        string fewMutations = "TCGTACGTACGTACGTACGT"; // 1 substitution
        string moreMutations = "TCGAACGTTCGTACGTACGT"; // superset: positions 0,3,8 changed

        double aniFew = ComparativeGenomics.CalculateANI(baseSeq, fewMutations, fragmentLength: FragLen, minIdentity: 0.0);
        double aniMore = ComparativeGenomics.CalculateANI(baseSeq, moreMutations, fragmentLength: FragLen, minIdentity: 0.0);

        Assert.That(aniMore, Is.LessThanOrEqualTo(aniFew + 1e-12),
            $"More divergence raised ANI: ANI(more)={aniMore} > ANI(few)={aniFew}");
    }

    /// <summary>
    /// INV-5 (D): Determinism. The same query/reference pair yields the same ANI on repeated calls.
    /// Evidence: CalculateANI is a pure function of its inputs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Ani_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(3 * FragLen), DnaArbitrary(3 * FragLen), (query, reference) =>
        {
            double a1 = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: FragLen);
            double a2 = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: FragLen);
            return (a1 == a2).Label("CalculateANI must be deterministic");
        });
    }

    /// <summary>
    /// INV-6 (boundary): ANI is 0 when no full-length fragment fits — empty inputs or a reference
    /// shorter than the fragment leave no qualifying match, so the mean is defined as 0
    /// (Goris et al. 2007: only fragments with a placeable, qualifying match contribute).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Ani_ReturnsZero_WhenNoFragmentQualifies()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ComparativeGenomics.CalculateANI("", "ACGTACGTACGT", fragmentLength: FragLen),
                Is.EqualTo(0.0), "empty query");
            Assert.That(ComparativeGenomics.CalculateANI("ACGTACGTACGT", "", fragmentLength: FragLen),
                Is.EqualTo(0.0), "empty reference");
            // Reference shorter than the fragment: no full-length placement exists.
            Assert.That(ComparativeGenomics.CalculateANI(new string('A', FragLen), "ACGT", fragmentLength: FragLen),
                Is.EqualTo(0.0), "reference shorter than fragment");
        });
    }

    #endregion

    #region COMPGEN-CLUSTER-001: R: cluster size ≥ minClusterSize; P: cluster is a common interval of every genome; M: lower minClusterSize → ≥ clusters; D: deterministic

    // A conserved gene cluster is a COMMON INTERVAL of the ortholog-group permutations: a set of
    // group labels that forms a contiguous window in every genome (Uno & Yagiura 2000; Heber &
    // Stoye 2001; Bui-Xuan, Habib & Paul 2013, Def. 1). The checklist's "lower threshold → ≥
    // clusters" knob is minClusterSize here (this model has no identity threshold).

    private static readonly string[] Alphabet = { "1", "2", "3", "4", "5", "6" };

    /// <summary>A permutation of <see cref="Alphabet"/> produced deterministically from a seed.</summary>
    private static Arbitrary<string[]> PermutationArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(Permute).ToArbitrary();

    private static string[] Permute(int seed)
    {
        var arr = (string[])Alphabet.Clone();
        var rng = new Random(seed);
        for (int i = arr.Length - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (arr[i], arr[j]) = (arr[j], arr[i]);
        }
        return arr;
    }

    /// <summary>Builds gene-ordered genomes (one gene per group label) plus the gene→group map.</summary>
    private static (List<IReadOnlyList<ComparativeGenomics.Gene>> genomes, Dictionary<string, string> map)
        Build(params string[][] orders)
    {
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>>(orders.Length);
        var map = new Dictionary<string, string>(StringComparer.Ordinal);
        for (int gi = 0; gi < orders.Length; gi++)
        {
            var genes = new List<ComparativeGenomics.Gene>(orders[gi].Length);
            for (int i = 0; i < orders[gi].Length; i++)
            {
                string id = $"G{gi}_{i}";
                genes.Add(new ComparativeGenomics.Gene(id, $"G{gi}", i * 100, i * 100 + 50, '+'));
                map[id] = orders[gi][i];
            }
            genomes.Add(genes);
        }
        return (genomes, map);
    }

    /// <summary>
    /// Independent (definition-level) test that a label set occupies a contiguous window of a
    /// permutation: the indices of its members span a gap-free range of exactly its own size.
    /// </summary>
    private static bool IsCommonInterval(string[][] orders, IReadOnlyList<string> cluster)
    {
        var set = cluster.ToHashSet(StringComparer.Ordinal);
        foreach (var order in orders)
        {
            int min = int.MaxValue, max = int.MinValue, count = 0;
            for (int i = 0; i < order.Length; i++)
            {
                if (!set.Contains(order[i])) continue;
                count++;
                if (i < min) min = i;
                if (i > max) max = i;
            }
            if (count != set.Count || max - min != set.Count - 1)
                return false;
        }
        return true;
    }

    private static HashSet<string> KeySet(IEnumerable<IReadOnlyList<string>> clusters)
        => clusters.Select(c => string.Join(",", c.OrderBy(x => x, StringComparer.Ordinal)))
                   .ToHashSet(StringComparer.Ordinal);

    /// <summary>
    /// INV-1 (P): Every reported cluster is genuinely a common interval — a contiguous window in
    /// EVERY genome — verified by an independent definition-level check, not the implementation's
    /// own routine. The full alphabet is always such an interval, so the result is also non-empty
    /// (guards against vacuous success).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Clusters_AreCommonIntervalsOfEveryGenome()
    {
        return Prop.ForAll(PermutationArbitrary(), PermutationArbitrary(), PermutationArbitrary(),
            (p1, p2, p3) =>
        {
            var orders = new[] { p1, p2, p3 };
            var (genomes, map) = Build(orders);
            var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2).ToList();

            bool allValid = clusters.All(c => IsCommonInterval(orders, c));
            bool fullSetPresent = clusters.Any(c => c.Count == Alphabet.Length);
            return (allValid && fullSetPresent)
                .Label($"validAll={allValid}, fullSetPresent={fullSetPresent}, n={clusters.Count}");
        });
    }

    /// <summary>
    /// INV-2 (R): Every reported cluster has at least <c>minClusterSize</c> distinct groups
    /// (and never fewer than 2, the smallest meaningful interval).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Clusters_RespectMinClusterSize()
    {
        return Prop.ForAll(PermutationArbitrary(), PermutationArbitrary(), (p1, p2) =>
        {
            const int minSize = 3;
            var (genomes, map) = Build(p1, p2);
            var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: minSize).ToList();
            return clusters.All(c => c.Count >= minSize)
                .Label($"a cluster smaller than {minSize} was reported");
        });
    }

    /// <summary>
    /// INV-3 (M): Lowering minClusterSize cannot remove clusters — the size-3 result is a subset of
    /// the size-2 result. Evidence: shrinking the size threshold only admits more candidate windows;
    /// the common-interval test itself is independent of the threshold.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Clusters_LowerMinSize_IsSuperset()
    {
        return Prop.ForAll(PermutationArbitrary(), PermutationArbitrary(), PermutationArbitrary(),
            (p1, p2, p3) =>
        {
            var (genomes, map) = Build(p1, p2, p3);
            var loose = KeySet(ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2));
            var strict = KeySet(ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 3));
            return strict.IsSubsetOf(loose)
                .Label($"min=3 result ({strict.Count}) not ⊆ min=2 result ({loose.Count})");
        });
    }

    /// <summary>
    /// INV-4 (D): Cluster detection is deterministic in both content and order.
    /// Evidence: FindConservedClusters sorts its output by size then label, with no random state.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Clusters_AreDeterministic()
    {
        return Prop.ForAll(PermutationArbitrary(), PermutationArbitrary(), (p1, p2) =>
        {
            var (genomes, map) = Build(p1, p2);
            var first = ComparativeGenomics.FindConservedClusters(genomes, map)
                .Select(c => string.Join(",", c)).ToList();
            var second = ComparativeGenomics.FindConservedClusters(genomes, map)
                .Select(c => string.Join(",", c)).ToList();
            return first.SequenceEqual(second)
                .Label("FindConservedClusters must be deterministic in content and order");
        });
    }

    #endregion

    #region COMPGEN-COMPARE-001: S: shared/specific metrics symmetric; R: synteny ∈ [0,1]; P: core+specific = genome size; D: deterministic

    // CompareGenomes partitions genes into the conserved (RBH) core and each genome's dispensable
    // set (pan-genome model, Tettelin et al. 2005) and reports OverallSynteny as a fraction.
    // RBH is a symmetric (bidirectional best hit) matching (Moreno-Hagelsieb & Latimer 2008), so
    // the conserved count is order-independent and the genome-specific counts swap under swapping.

    private static Gen<string> DnaGen(int len) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= len).Select(a => new string(a, 0, len));

    /// <summary>1..5 genes, each a 20-base DNA sequence (long enough for the k=5 similarity).</summary>
    private static Arbitrary<string[]> GenomeSeqsArbitrary() =>
        DnaGen(20).ArrayOf().Where(a => a.Length is >= 1 and <= 5).ToArbitrary();

    private static IReadOnlyList<ComparativeGenomics.Gene> MakeGenome(string genomeId, string[] seqs) =>
        seqs.Select((s, i) => new ComparativeGenomics.Gene($"{genomeId}_{i}", genomeId, i * 100, i * 100 + s.Length, '+', s))
            .ToList();

    /// <summary>
    /// INV-1 (S): The conserved-gene count is symmetric and the genome-specific counts swap under
    /// argument swapping, because RBH is a bidirectional-best-hit matching independent of which
    /// genome is the query (Moreno-Hagelsieb &amp; Latimer 2008).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Compare_SharedAndSpecificCounts_AreSymmetric()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var g1 = MakeGenome("A", s1);
            var g2 = MakeGenome("B", s2);
            var ab = ComparativeGenomics.CompareGenomes(g1, g2);
            var ba = ComparativeGenomics.CompareGenomes(g2, g1);
            return (ab.ConservedGenes == ba.ConservedGenes
                    && ab.GenomeSpecificGenes1 == ba.GenomeSpecificGenes2
                    && ab.GenomeSpecificGenes2 == ba.GenomeSpecificGenes1)
                .Label($"asymmetry: conserved {ab.ConservedGenes}/{ba.ConservedGenes}, " +
                       $"specific ({ab.GenomeSpecificGenes1},{ab.GenomeSpecificGenes2}) vs ({ba.GenomeSpecificGenes1},{ba.GenomeSpecificGenes2})");
        });
    }

    /// <summary>
    /// INV-2 (R): OverallSynteny is a fraction in [0,1] (syntenic genes ÷ smaller genome, clamped).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Compare_OverallSynteny_InUnitInterval()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var result = ComparativeGenomics.CompareGenomes(MakeGenome("A", s1), MakeGenome("B", s2));
            return (result.OverallSynteny is >= 0.0 and <= 1.0)
                .Label($"OverallSynteny={result.OverallSynteny} outside [0,1]");
        });
    }

    /// <summary>
    /// INV-3 (P): Pan-genome partition is exact — for each genome, conserved (core) + genome-specific
    /// (dispensable) equals its gene count (Tettelin et al. 2005); RBH maps each gene at most once.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Compare_CorePlusSpecific_EqualsGenomeSize()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var g1 = MakeGenome("A", s1);
            var g2 = MakeGenome("B", s2);
            var r = ComparativeGenomics.CompareGenomes(g1, g2);
            return (r.ConservedGenes + r.GenomeSpecificGenes1 == g1.Count
                    && r.ConservedGenes + r.GenomeSpecificGenes2 == g2.Count)
                .Label($"partition broken: core={r.ConservedGenes}, spec1={r.GenomeSpecificGenes1}/{g1.Count}, spec2={r.GenomeSpecificGenes2}/{g2.Count}");
        });
    }

    /// <summary>
    /// INV-4 (D): Comparison is deterministic across all reported metrics.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Compare_IsDeterministic()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var g1 = MakeGenome("A", s1);
            var g2 = MakeGenome("B", s2);
            var a = ComparativeGenomics.CompareGenomes(g1, g2);
            var b = ComparativeGenomics.CompareGenomes(g1, g2);
            return (a.ConservedGenes == b.ConservedGenes
                    && a.GenomeSpecificGenes1 == b.GenomeSpecificGenes1
                    && a.GenomeSpecificGenes2 == b.GenomeSpecificGenes2
                    && a.OverallSynteny == b.OverallSynteny
                    && a.Orthologs.Count == b.Orthologs.Count
                    && a.SyntenicBlocks.Count == b.SyntenicBlocks.Count)
                .Label("CompareGenomes must be deterministic");
        });
    }

    /// <summary>
    /// INV-5 (non-vacuity, positive control): two genomes sharing two near-identical genes plus one
    /// unique gene each yield exactly 2 conserved genes and 1 genome-specific gene per genome.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Compare_SharedGenes_AreCountedAsConserved()
    {
        const string x = "AAAAACCCCCGGGGGTTTTT";
        const string y = "TTTTTGGGGGCCCCCAAAAA";
        const string z = "ACGTACGTACGTACGTACGT"; // genome-1 only
        const string w = "GATCGATCGATCGATCGATC"; // genome-2 only

        var g1 = MakeGenome("A", new[] { x, y, z });
        var g2 = MakeGenome("B", new[] { x, y, w });
        var r = ComparativeGenomics.CompareGenomes(g1, g2);

        Assert.Multiple(() =>
        {
            Assert.That(r.ConservedGenes, Is.EqualTo(2), "shared genes X,Y must be conserved (RBH)");
            Assert.That(r.GenomeSpecificGenes1, Is.EqualTo(1), "Z is genome-1 specific");
            Assert.That(r.GenomeSpecificGenes2, Is.EqualTo(1), "W is genome-2 specific");
        });
    }

    #endregion

    #region COMPGEN-DOTPLOT-001: R: dot positions valid; P: each dot is a true word match; P: full main diagonal for identical seqs; D: deterministic

    // GenerateDotPlot is an EMBOSS-dottup-style word-match dot matrix (Gibbs & McIntyre 1970): it
    // reports every (x,y) where the wordSize-mer starting at sequence1[x] exactly equals the one at
    // sequence2[y]. Regions of similarity form diagonal runs; identical sequences put the whole main
    // diagonal in the plot.

    private const int DotWordSize = 4;

    /// <summary>
    /// INV-1 (R): Every reported dot lies in range — x ∈ [0, len1−w] and y ∈ [0, len2−w] — so a
    /// full word fits at both coordinates.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DotPlot_Coordinates_AreInRange()
    {
        return Prop.ForAll(DnaArbitrary(20), DnaArbitrary(25), (s1, s2) =>
        {
            var dots = ComparativeGenomics.GenerateDotPlot(s1, s2, wordSize: DotWordSize).ToList();
            bool inRange = dots.All(d =>
                d.x >= 0 && d.x <= s1.Length - DotWordSize &&
                d.y >= 0 && d.y <= s2.Length - DotWordSize);
            return inRange.Label($"a dot fell outside the valid coordinate range (n={dots.Count})");
        });
    }

    /// <summary>
    /// INV-2 (P): Every reported dot is a genuine exact word match — sequence1[x..x+w] equals
    /// sequence2[y..y+w]. This is the defining property of a word-match dot plot.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DotPlot_EveryDot_IsAnExactWordMatch()
    {
        return Prop.ForAll(DnaArbitrary(20), DnaArbitrary(25), (s1, s2) =>
        {
            var dots = ComparativeGenomics.GenerateDotPlot(s1, s2, wordSize: DotWordSize).ToList();
            bool allMatch = dots.All(d =>
                string.Equals(s1.Substring(d.x, DotWordSize), s2.Substring(d.y, DotWordSize), StringComparison.Ordinal));
            return allMatch.Label($"a reported dot was not an exact {DotWordSize}-mer match");
        });
    }

    /// <summary>
    /// INV-3 (P): For identical sequences the entire main diagonal {(i,i)} is present, because every
    /// word of a sequence trivially matches itself at the same offset.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DotPlot_IdenticalSequences_ContainFullMainDiagonal()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var dots = ComparativeGenomics.GenerateDotPlot(seq, seq, wordSize: DotWordSize).ToHashSet();
            bool diagonalComplete = Enumerable.Range(0, seq.Length - DotWordSize + 1)
                .All(i => dots.Contains((i, i)));
            return diagonalComplete.Label("main diagonal incomplete for identical sequences");
        });
    }

    /// <summary>
    /// INV-4 (D): The dot set is deterministic for a given input.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DotPlot_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(20), DnaArbitrary(25), (s1, s2) =>
        {
            var a = ComparativeGenomics.GenerateDotPlot(s1, s2, wordSize: DotWordSize).ToList();
            var b = ComparativeGenomics.GenerateDotPlot(s1, s2, wordSize: DotWordSize).ToList();
            return a.SequenceEqual(b).Label("GenerateDotPlot must be deterministic");
        });
    }

    /// <summary>
    /// INV-5 (boundary): empty when a sequence is shorter than the word; invalid word/step sizes
    /// throw <see cref="ArgumentOutOfRangeException"/>.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DotPlot_Boundaries()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ComparativeGenomics.GenerateDotPlot("ACG", "ACGTACGT", wordSize: DotWordSize).Any(),
                Is.False, "sequence shorter than word must yield no dots");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ComparativeGenomics.GenerateDotPlot("ACGTACGT", "ACGTACGT", wordSize: 0).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ComparativeGenomics.GenerateDotPlot("ACGTACGT", "ACGTACGT", wordSize: DotWordSize, stepSize: 0).ToList());
        });
    }

    #endregion

    #region COMPGEN-ORTHO-001: P: ortholog pairs bidirectional (RBH); R: ids/metrics valid; P: matching (each gene ≤ 1 pair); D: deterministic

    // FindOrthologs is reciprocal-best-hit ortholog detection (Tatusov et al. 1997; Moreno-Hagelsieb
    // & Latimer 2008): a pair (g1,g2) is orthologous iff each is the other's best qualifying hit in
    // the opposite genome. The RBH condition is symmetric in its two arguments.

    private static HashSet<(string, string)> OrthoPairs(IEnumerable<ComparativeGenomics.OrthologPair> pairs) =>
        pairs.Select(p => (p.Gene1Id, p.Gene2Id)).ToHashSet();

    /// <summary>
    /// INV-1 (P): Orthology is bidirectional — FindOrthologs(A,B) and FindOrthologs(B,A) describe the
    /// same matching with the gene roles swapped. This is the reciprocal-best-hit symmetry.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Orthologs_AreBidirectional()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var g1 = MakeGenome("A", s1);
            var g2 = MakeGenome("B", s2);
            var ab = OrthoPairs(ComparativeGenomics.FindOrthologs(g1, g2));
            var baSwapped = ComparativeGenomics.FindOrthologs(g2, g1)
                .Select(p => (p.Gene2Id, p.Gene1Id)).ToHashSet();
            return ab.SetEquals(baSwapped)
                .Label($"RBH not bidirectional: AB={ab.Count}, BA(swapped)={baSwapped.Count}");
        });
    }

    /// <summary>
    /// INV-2 (R): Each pair references a real gene of each genome and carries metrics in range —
    /// identity and coverage in [threshold,1], alignment length positive.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Orthologs_HaveValidMembershipAndMetrics()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var g1 = MakeGenome("A", s1);
            var g2 = MakeGenome("B", s2);
            var ids1 = g1.Select(g => g.Id).ToHashSet();
            var ids2 = g2.Select(g => g.Id).ToHashSet();
            var pairs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();
            bool ok = pairs.All(p =>
                ids1.Contains(p.Gene1Id) && ids2.Contains(p.Gene2Id) &&
                p.Identity is >= 0.3 and <= 1.0 &&
                p.Coverage is >= 0.5 and <= 1.0 &&
                p.AlignmentLength > 0);
            return ok.Label("an ortholog pair had an unknown gene id or out-of-range metric");
        });
    }

    /// <summary>
    /// INV-3 (P): The result is a matching — every gene of either genome appears in at most one pair
    /// (RBH maps each gene at most once).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Orthologs_FormAMatching()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var pairs = ComparativeGenomics.FindOrthologs(MakeGenome("A", s1), MakeGenome("B", s2)).ToList();
            bool gene1Unique = pairs.Select(p => p.Gene1Id).Distinct().Count() == pairs.Count;
            bool gene2Unique = pairs.Select(p => p.Gene2Id).Distinct().Count() == pairs.Count;
            return (gene1Unique && gene2Unique)
                .Label("a gene appeared in more than one ortholog pair");
        });
    }

    /// <summary>
    /// INV-4 (D): Ortholog detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Orthologs_AreDeterministic()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var g1 = MakeGenome("A", s1);
            var g2 = MakeGenome("B", s2);
            var a = OrthoPairs(ComparativeGenomics.FindOrthologs(g1, g2));
            var b = OrthoPairs(ComparativeGenomics.FindOrthologs(g1, g2));
            return a.SetEquals(b).Label("FindOrthologs must be deterministic");
        });
    }

    /// <summary>
    /// INV-5 (non-vacuity, positive control): two shared near-identical genes are recovered as the
    /// exact ortholog pairs; unrelated genes are not paired.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Orthologs_RecoverSharedGenes()
    {
        const string x = "AAAAACCCCCGGGGGTTTTT";
        const string y = "TTTTTGGGGGCCCCCAAAAA";
        const string z = "ACGTACGTACGTACGTACGT"; // genome-1 only
        const string w = "GATCGATCGATCGATCGATC"; // genome-2 only

        var g1 = MakeGenome("A", new[] { x, y, z });
        var g2 = MakeGenome("B", new[] { x, y, w });
        var pairs = OrthoPairs(ComparativeGenomics.FindOrthologs(g1, g2));

        Assert.That(pairs, Is.EquivalentTo(new[] { ("A_0", "B_0"), ("A_1", "B_1") }),
            "X↔X and Y↔Y must be the only reciprocal best hits");
    }

    #endregion

    #region COMPGEN-RBH-001: S: RBH(A,B)=RBH(B,A); P: each gene ≤ 1 RBH pair; P: only mutual best hits; D: deterministic

    // FindReciprocalBestHits keeps a pair (g1,g2) iff g1's best qualifying hit in genome 2 is g2 AND
    // g2's best qualifying hit in genome 1 is g1 (Tatusov et al. 1997; Moreno-Hagelsieb & Latimer
    // 2008). A one-directional best hit is explicitly NOT an ortholog.

    /// <summary>
    /// INV-1 (S): RBH is symmetric — the pair set of RBH(A,B) equals that of RBH(B,A) with the gene
    /// roles swapped, since the mutual-best-hit condition is symmetric in its arguments.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Rbh_IsSymmetric()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var g1 = MakeGenome("A", s1);
            var g2 = MakeGenome("B", s2);
            var ab = OrthoPairs(ComparativeGenomics.FindReciprocalBestHits(g1, g2));
            var ba = ComparativeGenomics.FindReciprocalBestHits(g2, g1)
                .Select(p => (p.Gene2Id, p.Gene1Id)).ToHashSet();
            return ab.SetEquals(ba).Label($"RBH(A,B)≠RBH(B,A): {ab.Count} vs {ba.Count}");
        });
    }

    /// <summary>
    /// INV-2 (P): The RBH set is a matching — each gene of either genome occurs in at most one pair.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Rbh_FormsAMatching()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var pairs = ComparativeGenomics.FindReciprocalBestHits(MakeGenome("A", s1), MakeGenome("B", s2)).ToList();
            return (pairs.Select(p => p.Gene1Id).Distinct().Count() == pairs.Count
                    && pairs.Select(p => p.Gene2Id).Distinct().Count() == pairs.Count)
                .Label("a gene appeared in more than one RBH pair");
        });
    }

    /// <summary>
    /// INV-3 (D): RBH detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Rbh_IsDeterministic()
    {
        return Prop.ForAll(GenomeSeqsArbitrary(), GenomeSeqsArbitrary(), (s1, s2) =>
        {
            var g1 = MakeGenome("A", s1);
            var g2 = MakeGenome("B", s2);
            var a = OrthoPairs(ComparativeGenomics.FindReciprocalBestHits(g1, g2));
            var b = OrthoPairs(ComparativeGenomics.FindReciprocalBestHits(g1, g2));
            return a.SetEquals(b).Label("FindReciprocalBestHits must be deterministic");
        });
    }

    /// <summary>
    /// INV-4 (P, reciprocity): comparing a genome of distinct genes to itself pairs every gene with
    /// itself — the canonical mutual-best-hit fixpoint (each gene's best match is itself, identity 1).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Rbh_SelfComparison_PairsEachGeneWithItself()
    {
        var genome = MakeGenome("A", new[]
        {
            "AAAAACCCCCGGGGGTTTTT",
            "TTTTTGGGGGCCCCCAAAAA",
            "ACGTACGTACGTACGTACGT"
        });
        var pairs = OrthoPairs(ComparativeGenomics.FindReciprocalBestHits(genome, genome));
        Assert.That(pairs, Is.EquivalentTo(new[] { ("A_0", "A_0"), ("A_1", "A_1"), ("A_2", "A_2") }));
    }

    /// <summary>
    /// INV-5 (P, one-directional exclusion): when two genome-1 genes are equally the best hit of a
    /// single genome-2 gene, only the reciprocal one is kept; the other is a one-directional best hit
    /// and is excluded ("a one-directional best hit is NOT an ortholog").
    /// </summary>
    [Test]
    [Category("Property")]
    public void Rbh_OneDirectionalHit_IsExcluded()
    {
        const string seq = "AAAAACCCCCGGGGGTTTTT";
        var g1 = MakeGenome("A", new[] { seq, seq }); // A_0 and A_1 both identical to the genome-2 gene
        var g2 = MakeGenome("B", new[] { seq });      // B_0 best-hits A_0 (ordinal tie-break), not A_1

        var pairs = OrthoPairs(ComparativeGenomics.FindReciprocalBestHits(g1, g2));
        Assert.That(pairs, Is.EquivalentTo(new[] { ("A_0", "B_0") }),
            "only the reciprocal pair survives; A_1's hit is one-directional");
    }

    #endregion

    #region COMPGEN-REARR-001: P: collinear genomes have no rearrangement; R: events have valid coordinates/types; P: type encoding consistent; D: deterministic

    // DetectRearrangements reads the orthologous markers in genome-1 order, relabels them to a signed
    // genome-2 order permutation (sign = relative strand), and yields one RearrangementEvent per
    // breakpoint of the extended permutation (Bafna & Pevzner 1998; Tannier et al. 2009). Only the two
    // operation classes derivable from a single signed in-order permutation are produced: Inversion
    // (sign reversal across the boundary) and Transposition (orientation-preserving discontinuity).

    /// <summary>
    /// Generates (σ, flip): a random genome-2 ordering permutation of n markers and per-marker strand
    /// flips, deterministically from a seed.
    /// </summary>
    private static Arbitrary<(int[] perm, bool[] flip)> RearrArbitrary(int n = 6) =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            var perm = Enumerable.Range(0, n).ToArray();
            for (int i = n - 1; i > 0; i--)
            {
                int j = rng.Next(i + 1);
                (perm[i], perm[j]) = (perm[j], perm[i]);
            }
            var flip = new bool[n];
            for (int j = 0; j < n; j++) flip[j] = rng.Next(2) == 0;
            return (perm, flip);
        }).ToArbitrary();

    /// <summary>
    /// Builds two gene-ordered genomes plus the ortholog map. genome-1 gene i (always '+') is
    /// orthologous to genome-2 gene σ(i); genome-2 gene j has '-' strand iff flip[j].
    /// </summary>
    private static (List<ComparativeGenomics.Gene> g1, List<ComparativeGenomics.Gene> g2, Dictionary<string, string> map)
        BuildRearr(int[] perm, bool[] flip)
    {
        int n = perm.Length;
        var g1 = new List<ComparativeGenomics.Gene>(n);
        var g2 = new List<ComparativeGenomics.Gene>(n);
        var map = new Dictionary<string, string>(n);
        for (int i = 0; i < n; i++)
            g1.Add(new ComparativeGenomics.Gene($"G1_{i}", "genome1", i * 100, i * 100 + 50, '+'));
        for (int j = 0; j < n; j++)
            g2.Add(new ComparativeGenomics.Gene($"G2_{j}", "genome2", j * 100, j * 100 + 50, flip[j] ? '-' : '+'));
        for (int i = 0; i < n; i++)
            map[$"G1_{i}"] = $"G2_{perm[i]}";
        return (g1, g2, map);
    }

    /// <summary>
    /// INV-1 (P): Collinear genomes have no rearrangement. When the orthologs appear in the same
    /// order and orientation in both genomes, the relabelled permutation is the identity and there
    /// are zero breakpoints.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Rearrangements_CollinearGenomes_HaveNoBreakpoints()
    {
        int n = 6;
        var (g1, g2, map) = BuildRearr(Enumerable.Range(0, n).ToArray(), new bool[n]);
        Assert.That(ComparativeGenomics.DetectRearrangements(g1, g2, map), Is.Empty,
            "identity gene order must produce no rearrangement events");
    }

    /// <summary>
    /// INV-2 (R): Every event has a valid genome-1 anchor coordinate, non-negative length, and one of
    /// the two derivable rearrangement types (Inversion or Transposition).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Rearrangements_HaveValidCoordinatesAndTypes()
    {
        return Prop.ForAll(RearrArbitrary(), pf =>
        {
            var (g1, g2, map) = BuildRearr(pf.perm, pf.flip);
            var starts = g1.Select(g => g.Start).ToHashSet();
            var events = ComparativeGenomics.DetectRearrangements(g1, g2, map).ToList();
            bool ok = events.All(e =>
                starts.Contains(e.Position) &&
                e.Length >= 0 &&
                (e.Type == ComparativeGenomics.RearrangementType.Inversion ||
                 e.Type == ComparativeGenomics.RearrangementType.Transposition));
            return ok.Label("an event had an invalid coordinate, negative length, or unexpected type");
        });
    }

    /// <summary>
    /// INV-3 (P): The stored event type is consistent with its breakpoint encoding —
    /// re-deriving the type from the event via ClassifyRearrangement reproduces Type.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Rearrangements_TypeEncoding_IsConsistent()
    {
        return Prop.ForAll(RearrArbitrary(), pf =>
        {
            var (g1, g2, map) = BuildRearr(pf.perm, pf.flip);
            var events = ComparativeGenomics.DetectRearrangements(g1, g2, map).ToList();
            return events.All(e => ComparativeGenomics.ClassifyRearrangement(e) == e.Type)
                .Label("ClassifyRearrangement disagreed with the stored event Type");
        });
    }

    /// <summary>
    /// INV-4 (D): Rearrangement detection is deterministic in content and order.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Rearrangements_AreDeterministic()
    {
        return Prop.ForAll(RearrArbitrary(), pf =>
        {
            var (g1, g2, map) = BuildRearr(pf.perm, pf.flip);
            var a = ComparativeGenomics.DetectRearrangements(g1, g2, map)
                .Select(e => (e.Type, e.Position, e.Length, e.TargetPosition)).ToList();
            var b = ComparativeGenomics.DetectRearrangements(g1, g2, map)
                .Select(e => (e.Type, e.Position, e.Length, e.TargetPosition)).ToList();
            return a.SequenceEqual(b).Label("DetectRearrangements must be deterministic");
        });
    }

    /// <summary>
    /// INV-5 (P, positive control): a whole-genome reversal (orthologs in reverse order with flipped
    /// strands) yields only Inversion events — a reversal negates the signs of the block it reverses.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Rearrangements_WholeGenomeReversal_AreAllInversions()
    {
        int n = 6;
        var reversed = Enumerable.Range(0, n).Select(i => n - 1 - i).ToArray();
        var allFlipped = Enumerable.Repeat(true, n).ToArray();
        var (g1, g2, map) = BuildRearr(reversed, allFlipped);

        var events = ComparativeGenomics.DetectRearrangements(g1, g2, map).ToList();
        Assert.That(events, Is.Not.Empty, "a whole-genome reversal must register breakpoints");
        Assert.That(events.All(e => e.Type == ComparativeGenomics.RearrangementType.Inversion), Is.True,
            "every breakpoint of a strand-flipped reversal must classify as an inversion");
    }

    #endregion

    #region COMPGEN-REVERSAL-001: R: distance ≥ 0; S: d(A,B)=d(B,A); I: d(A,A)=0; D: deterministic

    // CalculateReversalDistance returns the unsigned breakpoint lower bound ⌈b/2⌉ on the reversal
    // distance (Bafna & Pevzner 1998): a single reversal removes at most two breakpoints, so
    // d(π) ≥ b(π)/2. Breakpoint distance is a metric-like quantity over permutations of a marker set.

    private static Arbitrary<int[]> IntPermArbitrary(int n = 7) =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            var perm = Enumerable.Range(0, n).ToArray();
            for (int i = n - 1; i > 0; i--)
            {
                int j = rng.Next(i + 1);
                (perm[i], perm[j]) = (perm[j], perm[i]);
            }
            return perm;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R): The reversal-distance lower bound is non-negative for any two orderings of a marker
    /// set (it is ⌈b/2⌉ with b ≥ 0).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReversalDistance_IsNonNegative()
    {
        return Prop.ForAll(IntPermArbitrary(), IntPermArbitrary(), (a, b) =>
        {
            int d = ComparativeGenomics.CalculateReversalDistance(a, b);
            return (d >= 0).Label($"reversal distance {d} must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-2 (I): Self-distance is zero — a gene order requires no reversals to reach itself
    /// (the relative permutation is the identity, with no breakpoints).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReversalDistance_SelfIsZero()
    {
        return Prop.ForAll(IntPermArbitrary(), a =>
        {
            int d = ComparativeGenomics.CalculateReversalDistance(a, a);
            return (d == 0).Label($"d(A,A)={d}, expected 0");
        });
    }

    /// <summary>
    /// INV-3 (S): The distance is symmetric — d(A,B) = d(B,A). The breakpoint count depends only on
    /// the (unordered) shared adjacencies and matching endpoints of the two orderings.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReversalDistance_IsSymmetric()
    {
        return Prop.ForAll(IntPermArbitrary(), IntPermArbitrary(), (a, b) =>
        {
            int ab = ComparativeGenomics.CalculateReversalDistance(a, b);
            int ba = ComparativeGenomics.CalculateReversalDistance(b, a);
            return (ab == ba).Label($"d(A,B)={ab} ≠ d(B,A)={ba}");
        });
    }

    /// <summary>
    /// INV-4 (D): Distance computation is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReversalDistance_IsDeterministic()
    {
        return Prop.ForAll(IntPermArbitrary(), IntPermArbitrary(), (a, b) =>
        {
            return (ComparativeGenomics.CalculateReversalDistance(a, b)
                    == ComparativeGenomics.CalculateReversalDistance(a, b))
                .Label("CalculateReversalDistance must be deterministic");
        });
    }

    /// <summary>
    /// INV-5 (positive control): a single whole-sequence reversal is one reversal away from the
    /// identity — d(identity, reverse) = 1 (two endpoint breakpoints, ⌈2/2⌉ = 1).
    /// </summary>
    [Test]
    [Category("Property")]
    public void ReversalDistance_WholeReversal_IsOne()
    {
        var identity = new[] { 0, 1, 2, 3, 4 };
        var reversed = new[] { 4, 3, 2, 1, 0 };
        Assert.That(ComparativeGenomics.CalculateReversalDistance(identity, reversed), Is.EqualTo(1));
    }

    /// <summary>
    /// INV-6 (boundary): orderings over marker sets of different sizes are not comparable.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ReversalDistance_LengthMismatch_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => ComparativeGenomics.CalculateReversalDistance(new[] { 0, 1, 2 }, new[] { 0, 1 }));
    }

    #endregion
}
