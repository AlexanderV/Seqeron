namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Comparative Genomics area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Comparative")]
public class ComparativeCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: COMPGEN-ANI-001 — Average Nucleotide Identity (ANIb) (Comparative)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 131.
    // Spec: tests/TestSpecs/COMPGEN-ANI-001.md (ComparativeGenomics.CalculateANI).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Goris et al. (2007) ANIb — the query is cut into consecutive non-overlapping fragments of
    // fragmentLength nt; each is ungapped-aligned to the reference for its best match; fragments with
    // identity > minIdentity (0.30) and alignable fraction ≥ minAlignableFraction (0.70) contribute; ANI is
    // the mean of qualifying fragment identities.
    //
    // Checklist axes kmerSize(3) × genomeLen(3) × divergence(3) map onto the real knobs:
    //   • kmerSize   → fragmentLength ∈ {5, 10, 20}.
    //   • genomeLen  → sequence length ∈ {40, 80, 120}.
    //   • divergence → the fraction of query bases mutated from the reference ∈ {0.0, 0.1, 0.3}.
    // Grid = 3³ = 27 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: ANI is a JOINT function of fragment size, genome length and divergence. Each
    // cell re-derives the ANIb procedure independently from the constructed sequences and checks production
    // matches, plus the structural bounds (ANI ∈ [0,1], identical → 1.0).
    // ═══════════════════════════════════════════════════════════════════════

    private const string Bases = "ACGT";

    /// <summary>
    /// For every (fragment length, genome length, divergence) the production ANI equals the ANIb procedure
    /// re-derived independently from the constructed query/reference, and lies in [0, 1].
    /// </summary>
    [Test, Combinatorial]
    public void CalculateANI_FragmentLengthGenomeDivergenceGrid_MatchesAnibProcedure(
        [Values(5, 10, 20)] int fragmentLength,
        [Values(40, 80, 120)] int genomeLen,
        [Values(0.0, 0.1, 0.3)] double divergence)
    {
        string reference = BuildVariedDna(genomeLen);
        string query = MutateAtDensity(reference, divergence);

        double expected = GroundTruthAni(query, reference, fragmentLength, 0.30, 0.70);

        double ani = ComparativeGenomics.CalculateANI(query, reference, fragmentLength, 0.30, 0.70);

        ani.Should().BeApproximately(expected, 1e-9, "ANI = mean qualifying-fragment identity (Goris ANIb)");
        ani.Should().BeInRange(0.0, 1.0, "[INV-1] ANI ∈ [0, 1]");
        if (divergence == 0.0)
            ani.Should().BeApproximately(1.0, 1e-9, "[INV-2] identical sequences → ANI = 1.0");
    }

    /// <summary>
    /// Interaction witness (divergence axis monotonicity): at fixed fragment and genome length, higher
    /// query–reference divergence gives a lower (or equal) ANI. Source: Goris et al. (2007).
    /// </summary>
    [Test]
    public void CalculateANI_DivergenceAxis_IsMonotoneDecreasing()
    {
        string reference = BuildVariedDna(80);

        double identical = ComparativeGenomics.CalculateANI(reference, reference, 10);
        double low = ComparativeGenomics.CalculateANI(MutateAtDensity(reference, 0.1), reference, 10);
        double high = ComparativeGenomics.CalculateANI(MutateAtDensity(reference, 0.3), reference, 10);

        identical.Should().Be(1.0);
        identical.Should().BeGreaterThanOrEqualTo(low);
        low.Should().BeGreaterThanOrEqualTo(high);
    }

    /// <summary>
    /// Interaction witness (fragmentation, INV-4): the query is cut into floor(N / fragmentLength)
    /// consecutive non-overlapping fragments; a trailing partial fragment is ignored. With an identical
    /// reference every fragment is perfect so ANI = 1.0 regardless of the trailing remainder. Source:
    /// Goris et al. (2007) "consecutive fragments".
    /// </summary>
    [Test]
    public void CalculateANI_TrailingPartialFragment_IsIgnored()
    {
        string reference = BuildVariedDna(45); // 45 = 4×10 + 5 (trailing partial)
        ComparativeGenomics.CalculateANI(reference, reference, 10).Should().Be(1.0, "4 full fragments, partial ignored");
    }

    /// <summary>Witness (INV-5): empty input yields ANI 0.</summary>
    [Test]
    public void CalculateANI_EmptyInput_IsZero()
    {
        ComparativeGenomics.CalculateANI("", BuildVariedDna(40), 10).Should().Be(0.0);
        ComparativeGenomics.CalculateANI(BuildVariedDna(40), "", 10).Should().Be(0.0);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: COMPGEN-CLUSTER-001 — Conserved gene clusters (common intervals) (Comparative)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 132.
    // Spec: tests/TestSpecs/COMPGEN-CLUSTER-001.md (ComparativeGenomics.FindConservedClusters).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Uno & Yagiura (2000); Heber & Stoye (2001); Bui-Xuan et al. (2013) — a conserved cluster is
    // a COMMON INTERVAL: a set of ortholog groups that occupies a contiguous window in EVERY genome.
    //
    // Checklist axes nGenomes(3) × identityThreshold(3) map onto the real knobs (this is a synteny model,
    // not an identity one): nGenomes → number of genomes K ∈ {2,3,4}; identityThreshold → minClusterSize ∈
    // {2,3,4} (the cluster-size stringency). Grid = 3 × 3 = 9 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: the returned clusters are a JOINT function of how many genomes constrain the
    // synteny and the minimum cluster size — adding genomes can only remove common intervals, and raising
    // minClusterSize filters by size. Each cell re-derives the common intervals by brute force from the
    // permutation definition and checks production returns exactly that size-filtered set.
    // ═══════════════════════════════════════════════════════════════════════

    // Four genomes over labels {A,B,C,D,E} in which the block A,B,C stays contiguous (a synteny block).
    private static readonly string[][] ClusterGenomeLabels =
    {
        new[] { "A", "B", "C", "D", "E" },
        new[] { "D", "A", "B", "C", "E" },
        new[] { "E", "D", "A", "B", "C" },
        new[] { "A", "B", "C", "E", "D" },
    };

    /// <summary>
    /// For every (genome count, minClusterSize) the reported conserved clusters are exactly the common
    /// intervals of the permutations (re-derived by brute force) that meet the size threshold.
    /// </summary>
    [Test, Combinatorial]
    public void FindConservedClusters_GenomeCountMinSizeGrid_MatchesCommonIntervals(
        [Values(2, 3, 4)] int nGenomes,
        [Values(2, 3, 4)] int minClusterSize)
    {
        var labelSets = ClusterGenomeLabels.Take(nGenomes).ToList();
        var (genomes, map) = BuildGenomes(labelSets);

        var expected = BruteForceCommonIntervals(labelSets, minClusterSize);

        var actual = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize)
            .Select(c => string.Join(",", c.OrderBy(x => x, StringComparer.Ordinal)))
            .ToHashSet(StringComparer.Ordinal);

        actual.Should().BeEquivalentTo(expected, "reported clusters = size-filtered common intervals");
    }

    /// <summary>
    /// Interaction witness (genome-count axis tightens synteny): a block (D,E) that is contiguous in the
    /// first genome but split apart once a genome that separates them is added is no longer a common
    /// interval — adding genomes can only remove clusters. Source: Bui-Xuan et al. (2013) Def. 1.
    /// </summary>
    [Test]
    public void FindConservedClusters_AddingGenome_CanOnlyRemoveClusters()
    {
        var twoGenomes = ClusterGenomeLabels.Take(2).ToList();   // A B C D E ; D A B C E
        var threeGenomes = ClusterGenomeLabels.Take(3).ToList(); // + E D A B C

        var (g2, m2) = BuildGenomes(twoGenomes);
        var (g3, m3) = BuildGenomes(threeGenomes);

        var two = ComparativeGenomics.FindConservedClusters(g2, m2, 2)
            .Select(c => string.Join(",", c.OrderBy(x => x, StringComparer.Ordinal))).ToHashSet(StringComparer.Ordinal);
        var three = ComparativeGenomics.FindConservedClusters(g3, m3, 2)
            .Select(c => string.Join(",", c.OrderBy(x => x, StringComparer.Ordinal))).ToHashSet(StringComparer.Ordinal);

        three.Should().BeSubsetOf(two, "more genomes can only remove common intervals");
    }

    /// <summary>
    /// Interaction witness (INV-04): with fewer than two genomes the conserved-cluster question is vacuous
    /// and the result is empty. Source: common-interval family definition (K ≥ 2).
    /// </summary>
    [Test]
    public void FindConservedClusters_SingleGenome_IsEmpty()
    {
        var (genomes, map) = BuildGenomes(ClusterGenomeLabels.Take(1).ToList());

        ComparativeGenomics.FindConservedClusters(genomes, map, 2).Should().BeEmpty("fewer than 2 genomes → empty");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: COMPGEN-DOTPLOT-001 — Word-match dot plot (Comparative)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 134.
    // Spec: tests/TestSpecs/COMPGEN-DOTPLOT-001.md (ComparativeGenomics.GenerateDotPlot).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: EMBOSS dottup; Gibbs & McIntyre (1970) — a dot at (x,y) iff sequence1[x..x+w] equals
    // sequence2[y..y+w] (case-insensitive), x sampled by stepSize.
    //
    // Checklist axes wordSize(3) × seqLen(3) × strand(2) map onto the real knobs:
    //   • wordSize → the word (k-tuple) size ∈ {3, 5, 8}.
    //   • seqLen   → sequence length ∈ {10, 20, 40}.
    //   • strand   → the second sequence's orientation {Self (seq2 = seq1, forward), RevComp (seq2 =
    //     reverse complement of seq1)} — forward self-comparison fills the main diagonal; the
    //     reverse-complement plot reveals inverted (palindromic) word matches.
    // Grid = 3 × 3 × 2 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: the dot set is a JOINT function of word size, length and orientation. Each
    // cell re-derives the exact word-match set by brute force and checks production returns exactly it.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// For every (word size, length, orientation) the production dot set equals the exact word-match set
    /// re-derived by brute force; a forward self-comparison additionally contains the full main diagonal.
    /// </summary>
    [Test, Combinatorial]
    public void GenerateDotPlot_WordSizeLengthOrientationGrid_MatchesExactWordMatches(
        [Values(3, 5, 8)] int wordSize,
        [Values(10, 20, 40)] int seqLen,
        [Values(true, false)] bool selfOrientation)
    {
        string seq1 = BuildVariedDna(seqLen);
        string seq2 = selfOrientation ? seq1 : ReverseComplement(seq1);

        var expected = BruteForceDotPlot(seq1, seq2, wordSize);

        var actual = ComparativeGenomics.GenerateDotPlot(seq1, seq2, wordSize).ToHashSet();

        actual.Should().BeEquivalentTo(expected, "[INV-1] dot ⟺ exact word match");
        if (selfOrientation)
            for (int i = 0; i <= seqLen - wordSize; i++)
                actual.Should().Contain((i, i), "[INV-2] self-comparison fills the main diagonal");
    }

    /// <summary>
    /// Interaction witness (length × word, empty): when a sequence is shorter than the word size no dots are
    /// produced. Source: dottup undefined window. INV-3.
    /// </summary>
    [Test]
    public void GenerateDotPlot_SequenceShorterThanWord_IsEmpty()
    {
        ComparativeGenomics.GenerateDotPlot("ACGT", "ACGTACGT", wordSize: 8).Should().BeEmpty("seq1 length 4 < word 8");
    }

    /// <summary>
    /// Witness (INV-5): a non-positive word or step size throws. Source: sibling validation convention.
    /// </summary>
    [Test]
    public void GenerateDotPlot_NonPositiveParameters_Throw()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => ComparativeGenomics.GenerateDotPlot("ACGT", "ACGT", 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => ComparativeGenomics.GenerateDotPlot("ACGT", "ACGT", 2, 0));
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: COMPGEN-ORTHO-001 — Ortholog detection (reciprocal best hits) (Comparative)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 135.
    // Spec: tests/TestSpecs/COMPGEN-ORTHO-001.md (ComparativeGenomics.FindOrthologs).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Tatusov et al. (1997); Moreno-Hagelsieb & Latimer (2008) — orthologs are reciprocal best
    // hits (RBH) above an identity and coverage gate; similarity is alignment-free k-mer Jaccard.
    //
    // Checklist axes nGenes(3) × identityThreshold(3) × eValue(2) map onto the real knobs:
    //   • nGenes           → genes per genome ∈ {2, 4, 6} (each with a unique sequence, ortholog i ↔ i).
    //   • identityThreshold→ minIdentity ∈ {0.1, 0.3, 0.5}.
    //   • eValue           → there is no e-value (k-mer model); the analogous quality gate is minCoverage ∈
    //     {0.0, 0.5}.
    // Grid = 3 × 3 × 2 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: across every (gene count, identity, coverage) combination RBH yields a valid
    // ortholog MATCHING — every reported pair is reciprocal (i ↔ i), no gene is reused, and every pair
    // clears both thresholds. Identity-based gating is shown by a diverged-pair witness.
    // ═══════════════════════════════════════════════════════════════════════

    // Distinct dinucleotide patterns with disjoint 5-mer sets, so different genes never cross-match.
    private static readonly string[] OrthoPatterns = { "AC", "AG", "AT", "CG", "CT", "GT" };

    /// <summary>
    /// For every (gene count, identity threshold, coverage threshold) with identical orthologs (i ↔ i), RBH
    /// returns exactly that reciprocal matching: every pair clears both thresholds and no gene is reused.
    /// </summary>
    [Test, Combinatorial]
    public void FindOrthologs_GeneCountIdentityCoverageGrid_YieldsReciprocalMatching(
        [Values(2, 4, 6)] int nGenes,
        [Values(0.1, 0.3, 0.5)] double minIdentity,
        [Values(0.0, 0.5)] double minCoverage)
    {
        var genome1 = new List<ComparativeGenomics.Gene>();
        var genome2 = new List<ComparativeGenomics.Gene>();
        for (int i = 0; i < nGenes; i++)
        {
            string seq = Repeat(OrthoPatterns[i], 12);
            genome1.Add(new ComparativeGenomics.Gene($"G1_{i}", "G1", i * 100, i * 100 + 12, '+', seq));
            genome2.Add(new ComparativeGenomics.Gene($"G2_{i}", "G2", i * 100, i * 100 + 12, '+', seq));
        }

        var orthologs = ComparativeGenomics.FindOrthologs(genome1, genome2, minIdentity, minCoverage).ToList();

        orthologs.Should().HaveCount(nGenes, "each identical ortholog forms one reciprocal pair");
        orthologs.Should().OnlyContain(o => o.Gene1Id.Replace("G1_", "") == o.Gene2Id.Replace("G2_", ""), "[INV-1] reciprocal i ↔ i");
        orthologs.Should().OnlyContain(o => o.Identity >= minIdentity && o.Coverage >= minCoverage, "[INV-5] pairs clear both thresholds");
        orthologs.Select(o => o.Gene1Id).Should().OnlyHaveUniqueItems("[INV-2] genome-1 genes are matched at most once");
        orthologs.Select(o => o.Gene2Id).Should().OnlyHaveUniqueItems("[INV-2] genome-2 genes are matched at most once");
    }

    /// <summary>
    /// Interaction witness (identity-threshold gating): a diverged ortholog pair (k-mer Jaccard ≈ 1/6) is a
    /// reciprocal best hit at minIdentity 0.1 but excluded at 0.3 — the pair flips on the identity axis.
    /// Source: Moreno-Hagelsieb & Latimer (2008) identity gate.
    /// </summary>
    [Test]
    public void FindOrthologs_IdentityThreshold_GatesDivergedPair()
    {
        var genome1 = new[] { new ComparativeGenomics.Gene("G1_0", "G1", 0, 10, '+', "AAAAAAAAAA") };
        var genome2 = new[] { new ComparativeGenomics.Gene("G2_0", "G2", 0, 10, '+', "AAAAACCCCC") };

        ComparativeGenomics.FindOrthologs(genome1, genome2, minIdentity: 0.1, minCoverage: 0.0)
            .Should().ContainSingle("Jaccard ≈ 0.167 ≥ 0.1");
        ComparativeGenomics.FindOrthologs(genome1, genome2, minIdentity: 0.3, minCoverage: 0.0)
            .Should().BeEmpty("Jaccard ≈ 0.167 < 0.3");
    }

    /// <summary>
    /// Witness: genes without a sequence carry no similarity evidence and never form ortholog pairs. Source:
    /// RBH requires a comparable sequence.
    /// </summary>
    [Test]
    public void FindOrthologs_GenesWithoutSequence_ProduceNoPairs()
    {
        var genome1 = new[] { new ComparativeGenomics.Gene("G1_0", "G1", 0, 10, '+', null) };
        var genome2 = new[] { new ComparativeGenomics.Gene("G2_0", "G2", 0, 10, '+', null) };

        ComparativeGenomics.FindOrthologs(genome1, genome2).Should().BeEmpty("no sequence → no similarity");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Helpers — engineered constructs + independent ANIb ground truth
    // ───────────────────────────────────────────────────────────────────────

    private static string Repeat(string pattern, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++) chars[i] = pattern[i % pattern.Length];
        return new string(chars);
    }


    private static string ReverseComplement(string dna)
    {
        var chars = new char[dna.Length];
        for (int i = 0; i < dna.Length; i++)
        {
            char c = dna[dna.Length - 1 - i];
            chars[i] = c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c };
        }
        return new string(chars);
    }

    private static HashSet<(int, int)> BruteForceDotPlot(string s1, string s2, int wordSize)
    {
        var set = new HashSet<(int, int)>();
        string a = s1.ToUpperInvariant(), b = s2.ToUpperInvariant();
        for (int i = 0; i + wordSize <= a.Length; i++)
        {
            string word = a.Substring(i, wordSize);
            for (int j = 0; j + wordSize <= b.Length; j++)
                if (b.Substring(j, wordSize) == word) set.Add((i, j));
        }
        return set;
    }


    /// <summary>Builds genomes (Gene lists) and the shared gene→group map from per-genome ortholog-group label arrays.</summary>
    private static (List<IReadOnlyList<ComparativeGenomics.Gene>> genomes, Dictionary<string, string> map)
        BuildGenomes(IReadOnlyList<string[]> labelSets)
    {
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>>();
        var map = new Dictionary<string, string>();
        for (int g = 0; g < labelSets.Count; g++)
        {
            string genomeId = $"G{g}";
            var genes = new List<ComparativeGenomics.Gene>();
            for (int i = 0; i < labelSets[g].Length; i++)
            {
                string geneId = $"{genomeId}_{i}";
                genes.Add(new ComparativeGenomics.Gene(geneId, genomeId, i * 100, i * 100 + 50, '+'));
                map[geneId] = labelSets[g][i];
            }
            genomes.Add(genes);
        }
        return (genomes, map);
    }

    /// <summary>
    /// Brute-force common intervals (Bui-Xuan et al. 2013, Def. 1): every label set that occupies a
    /// contiguous, foreign-free window in EVERY genome, filtered to size ≥ max(minSize, 2). Returns the
    /// canonical sorted-joined keys.
    /// </summary>
    private static HashSet<string> BruteForceCommonIntervals(IReadOnlyList<string[]> genomes, int minSize)
    {
        int effectiveMin = Math.Max(minSize, 2);
        var result = new HashSet<string>(StringComparer.Ordinal);
        string[] g0 = genomes[0];

        for (int start = 0; start < g0.Length; start++)
        {
            var set = new SortedSet<string>(StringComparer.Ordinal);
            for (int end = start; end < g0.Length; end++)
            {
                set.Add(g0[end]);
                if (set.Count < effectiveMin) continue;
                if (IsCommonInterval(set, genomes))
                    result.Add(string.Join(",", set));
            }
        }
        return result;
    }

    private static bool IsCommonInterval(SortedSet<string> set, IReadOnlyList<string[]> genomes)
    {
        foreach (string[] genome in genomes)
        {
            int min = int.MaxValue, max = int.MinValue, found = 0;
            for (int i = 0; i < genome.Length; i++)
            {
                if (set.Contains(genome[i]))
                {
                    min = Math.Min(min, i);
                    max = Math.Max(max, i);
                    found++;
                }
            }
            if (found != set.Count) return false;           // a member is missing from this genome
            if (max - min + 1 != set.Count) return false;    // a foreign group sits inside the window
        }
        return true;
    }


    /// <summary>Builds a low-repeat DNA sequence so that the best ungapped placement is at the aligned offset.</summary>
    private static string BuildVariedDna(int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = Bases[(i * 7 + (i / 4) * 3 + 1) % 4];
        return new string(chars);
    }

    /// <summary>Mutates every (1/divergence)-th base to a different base; divergence 0 returns the input unchanged.</summary>
    private static string MutateAtDensity(string sequence, double divergence)
    {
        if (divergence <= 0.0) return sequence;
        int stride = (int)Math.Round(1.0 / divergence);
        var chars = sequence.ToCharArray();
        for (int i = 0; i < chars.Length; i += stride)
            chars[i] = Bases[(Array.IndexOf(Bases.ToCharArray(), chars[i]) + 1) % 4]; // flip to next base
        return new string(chars);
    }

    /// <summary>Independent ANIb ground truth (Goris et al. 2007): mean identity of qualifying fragments.</summary>
    private static double GroundTruthAni(string query, string reference, int fragLen, double minId, double minFrac)
    {
        if (string.IsNullOrEmpty(query) || string.IsNullOrEmpty(reference)) return 0.0;
        string q = query.ToUpperInvariant();
        string r = reference.ToUpperInvariant();

        var qualifying = new List<double>();
        for (int start = 0; start + fragLen <= q.Length; start += fragLen)
        {
            (double id, double frac) = BestUngappedMatch(q.Substring(start, fragLen), r);
            if (id > minId && frac >= minFrac)
                qualifying.Add(id);
        }
        return qualifying.Count > 0 ? qualifying.Average() : 0.0;
    }

    private static (double identity, double alignableFraction) BestUngappedMatch(string fragment, string reference)
    {
        int fragLen = fragment.Length;
        if (reference.Length < fragLen) return (0.0, 0.0);

        int best = 0;
        for (int offset = 0; offset + fragLen <= reference.Length; offset++)
        {
            int matches = 0;
            for (int k = 0; k < fragLen; k++)
                if (fragment[k] == reference[offset + k]) matches++;
            if (matches > best) best = matches;
        }
        return (Math.Min((double)best / fragLen, 1.0), 1.0);
    }
}
