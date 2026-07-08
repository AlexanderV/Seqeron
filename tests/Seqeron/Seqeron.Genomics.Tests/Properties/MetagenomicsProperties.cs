using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for metagenomics diversity analysis:
/// alpha diversity indices and beta diversity metrics.
///
/// Test Units: META-ALPHA-001, META-BETA-001, META-FUNC-001, META-PATHWAY-001, META-RESIST-001, META-TAXA-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Metagenomics")]
public class MetagenomicsProperties
{
    #region Generators

    /// <summary>
    /// Generates a random abundance dictionary with N species, values normalized to sum to 1.0.
    /// </summary>
    private static Arbitrary<IReadOnlyDictionary<string, double>> AbundanceArbitrary(int minSpecies = 2) =>
        Gen.Choose(minSpecies, 20)
            .SelectMany(n =>
                Gen.Choose(1, 100).ArrayOf(n)
                    .Select(counts =>
                    {
                        double total = counts.Sum();
                        var dict = new Dictionary<string, double>();
                        for (int i = 0; i < n; i++)
                            dict[$"species_{i + 1}"] = counts[i] / total;
                        return (IReadOnlyDictionary<string, double>)dict;
                    }))
            .ToArbitrary();

    /// <summary>
    /// Creates a uniform abundance for N species (equal proportions).
    /// </summary>
    private static IReadOnlyDictionary<string, double> UniformAbundance(int n)
    {
        var dict = new Dictionary<string, double>();
        for (int i = 0; i < n; i++)
            dict[$"species_{i + 1}"] = 1.0 / n;
        return dict;
    }

    /// <summary>
    /// Creates a single-species dominance abundance.
    /// </summary>
    private static IReadOnlyDictionary<string, double> SingleSpeciesAbundance()
    {
        return new Dictionary<string, double> { ["species_1"] = 1.0 };
    }

    #endregion

    #region META-ALPHA-001: R: Shannon ≥ 0; R: Simpson ∈ [0,1]; M: more species → higher diversity; D: deterministic

    /// <summary>
    /// INV-1: Shannon index H ≥ 0 for any abundance distribution.
    /// Evidence: H = −∑ pᵢ ln(pᵢ) where pᵢ ∈ (0,1] → each term is ≥ 0, so H ≥ 0.
    /// Source: Shannon (1948); Wikipedia "Diversity index".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_ShannonIndex_NonNegative()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            return (result.ShannonIndex >= -1e-9)
                .Label($"Shannon={result.ShannonIndex:F6} must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-2: Simpson index λ ∈ [0, 1] for any abundance distribution.
    /// Evidence: λ = ∑ pᵢ² where pᵢ ∈ [0,1] and ∑pᵢ = 1; 
    /// by Cauchy-Schwarz: 1/S ≤ λ ≤ 1.
    /// Source: Simpson (1949); Wikipedia "Diversity index".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_SimpsonIndex_InRange()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            return (result.SimpsonIndex >= -1e-9 && result.SimpsonIndex <= 1.0 + 1e-9)
                .Label($"Simpson={result.SimpsonIndex:F6} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-3: Inverse Simpson ≥ 1.0 for any distribution with at least one species.
    /// Evidence: D = 1/λ where λ ≤ 1, so D ≥ 1.
    /// Source: Hill (1973); Wikipedia "Diversity index".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_InverseSimpson_AtLeast1()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            return (result.InverseSimpson >= 1.0 - 1e-9)
                .Label($"InverseSimpson={result.InverseSimpson:F4} must be ≥ 1");
        });
    }

    /// <summary>
    /// INV-4: Shannon entropy is maximized (H = ln S) when all species are equally abundant.
    /// Evidence: By Gibbs' inequality, uniform distribution maximizes entropy.
    /// Source: Shannon (1948); Wikipedia "Diversity index".
    /// </summary>
    [TestCase(2)]
    [TestCase(5)]
    [TestCase(10)]
    [Category("Property")]
    public void AlphaDiversity_UniformDistribution_MaxShannon(int n)
    {
        var abundances = UniformAbundance(n);
        var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        double expectedMaxH = Math.Log(n);
        Assert.That(result.ShannonIndex, Is.EqualTo(expectedMaxH).Within(0.01),
            $"Uniform distribution of {n} species should have H = ln({n}) = {expectedMaxH:F4}, got {result.ShannonIndex:F4}");
    }

    /// <summary>
    /// INV-5: Single species → Shannon H = 0 (no uncertainty).
    /// Evidence: H = −(1 × ln 1) = 0.
    /// Source: Shannon (1948).
    /// </summary>
    [Test]
    [Category("Property")]
    public void AlphaDiversity_SingleSpecies_ShannonIsZero()
    {
        var abundances = SingleSpeciesAbundance();
        var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);

        Assert.That(result.ShannonIndex, Is.EqualTo(0.0).Within(0.001),
            $"Single species should have Shannon H = 0, got {result.ShannonIndex:F4}");
    }

    /// <summary>
    /// INV-6: Pielou evenness J ∈ [0, 1] for distributions with &gt; 1 species.
    /// Evidence: J = H / ln(S) where H ≤ ln(S), so J ∈ [0, 1].
    /// Source: Pielou (1966); Wikipedia "Species evenness".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_PielouEvenness_InRange()
    {
        return Prop.ForAll(AbundanceArbitrary(2), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            return (result.PielouEvenness >= -1e-9 && result.PielouEvenness <= 1.0 + 1e-9)
                .Label($"Pielou J={result.PielouEvenness:F4} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-7: Alpha diversity is deterministic.
    /// Evidence: CalculateAlphaDiversity is a pure mathematical function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_IsDeterministic()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var r1 = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            var r2 = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            bool same = Math.Abs(r1.ShannonIndex - r2.ShannonIndex) < 1e-10 &&
                        Math.Abs(r1.SimpsonIndex - r2.SimpsonIndex) < 1e-10 &&
                        r1.ObservedSpecies == r2.ObservedSpecies;
            return same.Label("CalculateAlphaDiversity must be deterministic");
        });
    }

    /// <summary>
    /// INV-8: ObservedSpecies equals the number of taxa with non-zero abundance.
    /// Evidence: By definition, observed species count = |{i : pᵢ &gt; 0}|.
    /// Source: Wikipedia "Species richness".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaDiversity_ObservedSpecies_EqualsNonZeroTaxa()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
            int expected = abundances.Values.Count(v => v > 0);
            return (result.ObservedSpecies == expected)
                .Label($"ObservedSpecies={result.ObservedSpecies}, expected={expected}");
        });
    }

    #endregion

    #region META-BETA-001: S: dist(a,b)=dist(b,a); R: dist ∈ [0,1]; I: dist(x,x)=0; D: deterministic

    /// <summary>
    /// INV-1: Bray-Curtis ∈ [0, 1].
    /// Evidence: BC = 1 − 2C/(S₁+S₂) where 0 ≤ C ≤ min(S₁,S₂).
    /// Source: Bray &amp; Curtis (1957).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_BrayCurtis_InRange()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var result = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            return (result.BrayCurtis >= -1e-9 && result.BrayCurtis <= 1.0 + 1e-9)
                .Label($"BrayCurtis={result.BrayCurtis:F6} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-2: Jaccard distance ∈ [0, 1].
    /// Evidence: J = 1 − |A∩B|/|A∪B| where |A∩B| ≤ |A∪B|.
    /// Source: Jaccard (1901).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_JaccardDistance_InRange()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var result = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            return (result.JaccardDistance >= -1e-9 && result.JaccardDistance <= 1.0 + 1e-9)
                .Label($"Jaccard={result.JaccardDistance:F6} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-3: Symmetry — BrayCurtis(a,b) = BrayCurtis(b,a).
    /// Evidence: The formula is symmetric in its inputs.
    /// Source: Bray &amp; Curtis (1957).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_BrayCurtis_IsSymmetric()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var r1 = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            var r2 = MetagenomicsAnalyzer.CalculateBetaDiversity("s2", a2, "s1", a1);
            return (Math.Abs(r1.BrayCurtis - r2.BrayCurtis) < 1e-9)
                .Label($"BrayCurtis not symmetric: {r1.BrayCurtis:F6} vs {r2.BrayCurtis:F6}");
        });
    }

    /// <summary>
    /// INV-4: Symmetry — Jaccard(a,b) = Jaccard(b,a).
    /// Evidence: Jaccard distance is symmetric by definition.
    /// Source: Jaccard (1901).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_JaccardDistance_IsSymmetric()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var r1 = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            var r2 = MetagenomicsAnalyzer.CalculateBetaDiversity("s2", a2, "s1", a1);
            return (Math.Abs(r1.JaccardDistance - r2.JaccardDistance) < 1e-9)
                .Label($"Jaccard not symmetric: {r1.JaccardDistance:F6} vs {r2.JaccardDistance:F6}");
        });
    }

    /// <summary>
    /// INV-5: Identity — distance of a sample to itself is 0.
    /// Evidence: BC = 1 − 2C/(S+S) = 1 − 1 = 0 when both samples are identical.
    /// Jaccard = 0 when all species are shared.
    /// Source: Bray &amp; Curtis (1957); Jaccard (1901).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_SelfDistance_IsZero()
    {
        return Prop.ForAll(AbundanceArbitrary(), abundances =>
        {
            var result = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", abundances, "s2", abundances);
            return (Math.Abs(result.BrayCurtis) < 1e-9 && Math.Abs(result.JaccardDistance) < 1e-9)
                .Label($"Self-distance must be 0: BC={result.BrayCurtis:F6}, Jaccard={result.JaccardDistance:F6}");
        });
    }

    /// <summary>
    /// INV-6: Disjoint samples → maximum Jaccard distance = 1.0.
    /// Evidence: When A ∩ B = ∅, Jaccard = 1 − 0/|A∪B| = 1.
    /// Source: Jaccard (1901).
    /// </summary>
    [Test]
    [Category("Property")]
    public void BetaDiversity_DisjointSamples_JaccardIsOne()
    {
        var s1 = new Dictionary<string, double> { ["A"] = 0.5, ["B"] = 0.5 };
        var s2 = new Dictionary<string, double> { ["C"] = 0.5, ["D"] = 0.5 };

        var result = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", s1, "s2", s2);
        Assert.That(result.JaccardDistance, Is.EqualTo(1.0).Within(0.001),
            $"Disjoint samples should have Jaccard = 1.0, got {result.JaccardDistance:F4}");
    }

    /// <summary>
    /// INV-7: Beta diversity is deterministic.
    /// Evidence: CalculateBetaDiversity is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_IsDeterministic()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var r1 = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            var r2 = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            bool same = Math.Abs(r1.BrayCurtis - r2.BrayCurtis) < 1e-10 &&
                        Math.Abs(r1.JaccardDistance - r2.JaccardDistance) < 1e-10 &&
                        r1.SharedSpecies == r2.SharedSpecies;
            return same.Label("CalculateBetaDiversity must be deterministic");
        });
    }

    /// <summary>
    /// INV-8: SharedSpecies + UniqueToSample1 + UniqueToSample2 = total unique species.
    /// Evidence: Every species across both samples is classified into exactly one category.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BetaDiversity_SpeciesCounts_SumToTotal()
    {
        return Prop.ForAll(AbundanceArbitrary(), AbundanceArbitrary(), (a1, a2) =>
        {
            var result = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", a1, "s2", a2);
            int total = a1.Keys.Union(a2.Keys)
                .Count(k => (a1.GetValueOrDefault(k) > 0) || (a2.GetValueOrDefault(k) > 0));
            int sum = result.SharedSpecies + result.UniqueToSample1 + result.UniqueToSample2;
            return (sum == total)
                .Label($"Species count: shared+unique1+unique2={sum}, total={total}");
        });
    }

    #endregion

    #region META-CLASS-001: R: confidence ∈ [0,1]; P: assigned taxon in tree; D: deterministic (Kraken k-mer/LCA/RTL)

    // META-CLASS-001 — Kraken (Wood & Salzberg 2014) k-mer / LCA / RTL taxonomic classification.
    //
    // The tests below treat the Kraken model as the specification and verify the implementation
    // against INDEPENDENT oracles derived purely from the doc + theory:
    //   * canonical k-mer        = min(kmer, revcomp(kmer))   — re-implemented here, not reused.
    //   * database value         = LCA of the taxa owning a shared k-mer (INV-06).
    //   * per-read assignment     = leaf of the max-weight root-to-leaf (RTL) path, with the LCA of
    //                              the tied maximal leaves on a tie (the key Kraken behavior).
    //   * confidence             = C/Q with C = clade-rooted k-mers, Q = non-ambiguous k-mers.
    // A small k = 4 is used throughout so reads can be short and generation is fast.

    private const int ClassK = 4;

    /// <summary>
    /// Fixed taxonomy fixture reused across META-CLASS-001:
    /// <code>
    /// root(1) ── kingdom(2) ─┬─ genus(10) ─┬─ speciesA(100)
    ///                        │             └─ speciesB(101)
    ///                        └─ genus(20) ──── speciesC(200)
    /// </code>
    /// Lca(100,101)=10 (sibling species → genus); Lca(100,200)=2 (cross-genus → kingdom).
    /// </summary>
    private static TaxonomyTree BuildClassFixture() => new(new[]
    {
        new TaxonNode(1, "root", "root", 1),
        new TaxonNode(2, "Bacteria", "kingdom", 1),
        new TaxonNode(10, "Genus_X", "genus", 2),
        new TaxonNode(100, "Species_A", "species", 10),
        new TaxonNode(101, "Species_B", "species", 10),
        new TaxonNode(20, "Genus_Y", "genus", 2),
        new TaxonNode(200, "Species_C", "species", 20),
    });

    /// <summary>Independent reverse-complement oracle over the A/C/G/T alphabet.</summary>
    private static string ReverseComplement(string s)
    {
        var chars = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
        {
            char c = s[s.Length - 1 - i];
            chars[i] = c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c };
        }
        return new string(chars);
    }

    /// <summary>Independent canonical-k-mer oracle: min(kmer, revcomp(kmer)) by ordinal order.</summary>
    private static string Canonical(string kmer)
    {
        string rc = ReverseComplement(kmer);
        return string.CompareOrdinal(kmer, rc) <= 0 ? kmer : rc;
    }

    /// <summary>Distinct non-ambiguous canonical k-mers of a sequence (independent Q oracle).</summary>
    private static List<string> CanonicalKmers(string seq, int k)
    {
        var result = new List<string>();
        string s = (seq ?? "").ToUpperInvariant();
        for (int i = 0; i + k <= s.Length; i++)
        {
            string kmer = s.Substring(i, k);
            if (kmer.All(c => "ACGT".Contains(c)))
                result.Add(Canonical(kmer));
        }
        return result;
    }

    /// <summary>Generates random ACGT sequences of length ≥ k (default ≥ ClassK).</summary>
    private static Gen<string> AcgtReadGen(int minLen = ClassK) =>
        Gen.Choose(minLen, 24)
            .SelectMany(len =>
                Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
                    .Select(cs => new string(cs)));

    private static Arbitrary<string> AcgtReadArbitrary(int minLen = ClassK) =>
        AcgtReadGen(minLen).ToArbitrary();

    private static Arbitrary<string[]> AcgtReadArrayArbitrary() =>
        AcgtReadGen().ArrayOf().ToArbitrary();

    // Reference sequences with controlled canonical k-mer content (k = 4). Non-repetitive so that
    // every length-4 window is a distinct canonical k-mer.
    //   refA   — unique to species A; shares no canonical k-mer with refB or refC.
    //   refB   — unique to species B.
    //   refC   — unique to species C (other genus).
    //   refSharedAB — given once to A and once to B during DB build so its k-mers collapse to Lca=10.
    private const string RefA = "ACGTACAGTC";          // species A
    private const string RefB = "GGATCATGCA";          // species B
    private const string RefC = "TTAGCTGACT";          // species C
    private const string RefShared = "ACGGTACC";       // shared between A and B at build time

    /// <summary>
    /// INV-01: <c>ClassifyReads</c> yields exactly one result per input read, preserving
    /// <c>ReadId</c> and input order. Source: doc §2.4 INV-01.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClassifyReads_OutputCount_EqualsInputCount()
    {
        var taxonomy = BuildClassFixture();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, RefA), (101, RefB), (200, RefC) }, taxonomy, ClassK);

        return Prop.ForAll(AcgtReadArrayArbitrary(), seqs =>
        {
            var reads = seqs.Select((s, i) => ($"read_{i}", s)).ToList();
            var results = MetagenomicsAnalyzer.ClassifyReads(reads, db, taxonomy, ClassK).ToList();
            bool countOk = results.Count == reads.Count;
            bool idsOk = results.Select(r => r.ReadId).SequenceEqual(reads.Select(r => r.Item1));
            return (countOk && idsOk)
                .Label($"expected {reads.Count} results in order, got {results.Count}");
        });
    }

    /// <summary>
    /// INV-02 (R: confidence ∈ [0,1]): every <c>Confidence</c> lies in [0,1] for arbitrary reads,
    /// because it is C/Q with 0 ≤ C ≤ Q (else 0). Source: doc §2.4 INV-02.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClassifyReads_Confidence_InUnitInterval()
    {
        var taxonomy = BuildClassFixture();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, RefA), (101, RefB), (200, RefC) }, taxonomy, ClassK);

        return Prop.ForAll(AcgtReadArbitrary(), seq =>
        {
            var r = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", seq) }, db, taxonomy, ClassK).Single();
            return (r.Confidence >= -1e-12 && r.Confidence <= 1.0 + 1e-12)
                .Label($"Confidence={r.Confidence:F6} must be in [0,1]");
        });
    }

    /// <summary>
    /// INV-03: <c>MatchedKmers</c> (C) ≤ <c>TotalKmers</c> (Q), both ≥ 0 — the assigned clade's
    /// k-mers are a subset of the queried non-ambiguous k-mers. Source: doc §2.4 INV-03.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClassifyReads_MatchedKmers_AtMostTotalKmers()
    {
        var taxonomy = BuildClassFixture();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, RefA), (101, RefB), (200, RefC) }, taxonomy, ClassK);

        return Prop.ForAll(AcgtReadArbitrary(), seq =>
        {
            var r = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", seq) }, db, taxonomy, ClassK).Single();
            return (r.MatchedKmers >= 0 && r.TotalKmers >= 0 && r.MatchedKmers <= r.TotalKmers)
                .Label($"0 ≤ C({r.MatchedKmers}) ≤ Q({r.TotalKmers})");
        });
    }

    /// <summary>
    /// INV-03 (oracle): Q equals the count of non-ambiguous length-k windows (here all windows,
    /// since reads are pure ACGT), computed independently. Source: doc §3.3 (Q = non-ambiguous
    /// k-mers queried).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClassifyReads_TotalKmers_EqualsNonAmbiguousWindowCount()
    {
        var taxonomy = BuildClassFixture();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, RefA), (101, RefB), (200, RefC) }, taxonomy, ClassK);

        return Prop.ForAll(AcgtReadArbitrary(), seq =>
        {
            var r = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", seq) }, db, taxonomy, ClassK).Single();
            int expectedQ = CanonicalKmers(seq, ClassK).Count; // pure ACGT ⇒ every window counts
            return (r.TotalKmers == expectedQ)
                .Label($"Q={r.TotalKmers}, expected {expectedQ}");
        });
    }

    /// <summary>
    /// INV-05 (P: assigned taxon in tree): every assigned <c>TaxonId</c> satisfies
    /// <c>taxonomy.Contains</c>. Source: doc §2.4 INV-05.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClassifyReads_AssignedTaxon_IsInTree()
    {
        var taxonomy = BuildClassFixture();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, RefA), (101, RefB), (200, RefC) }, taxonomy, ClassK);

        return Prop.ForAll(AcgtReadArbitrary(), seq =>
        {
            var r = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", seq) }, db, taxonomy, ClassK).Single();
            return taxonomy.Contains(r.TaxonId)
                .Label($"assigned TaxonId {r.TaxonId} is not in the taxonomy tree");
        });
    }

    /// <summary>
    /// D (determinism): identical inputs produce identical classification sequences (all fields).
    /// Source: doc — exact-k-mer Kraken is a pure deterministic function of (reads, db, tree, k).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClassifyReads_IsDeterministic()
    {
        var taxonomy = BuildClassFixture();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, RefA), (101, RefB), (200, RefC) }, taxonomy, ClassK);

        return Prop.ForAll(AcgtReadArrayArbitrary(), seqs =>
        {
            var reads = seqs.Select((s, i) => ($"read_{i}", s)).ToList();
            var r1 = MetagenomicsAnalyzer.ClassifyReads(reads, db, taxonomy, ClassK).ToList();
            var r2 = MetagenomicsAnalyzer.ClassifyReads(reads, db, taxonomy, ClassK).ToList();
            return r1.SequenceEqual(r2).Label("ClassifyReads must be deterministic");
        });
    }

    /// <summary>
    /// Canonical strand symmetry (I / D-flavored): a read and its reverse complement classify to the
    /// same <c>TaxonId</c> and <c>Confidence</c>, because the DB and lookup index canonical k-mers
    /// min(kmer, revcomp(kmer)) and revcomp reverses the window order without changing the canonical
    /// multiset. Source: doc §2.2 (canonical k-mer for double-stranded DNA).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClassifyReads_ReverseComplement_SameAssignment()
    {
        var taxonomy = BuildClassFixture();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, RefA), (101, RefB), (200, RefC) }, taxonomy, ClassK);

        return Prop.ForAll(AcgtReadArbitrary(), seq =>
        {
            var fwd = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", seq) }, db, taxonomy, ClassK).Single();
            var rev = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", ReverseComplement(seq)) }, db, taxonomy, ClassK).Single();
            return (fwd.TaxonId == rev.TaxonId && Math.Abs(fwd.Confidence - rev.Confidence) < 1e-12)
                .Label($"fwd→{fwd.TaxonId}/{fwd.Confidence:F6} vs rev→{rev.TaxonId}/{rev.Confidence:F6}");
        });
    }

    /// <summary>
    /// INV-04 (no-hit, full-length read): a read whose canonical k-mers are all absent from the DB
    /// is unclassified ⇒ <c>TaxonId == RootId</c>, <c>MatchedKmers == 0</c>; Q is still recorded.
    /// Here the read is classified against an EMPTY database, so no k-mer can hit.
    /// Source: doc §2.2 / §6.1 (no database matches → root, Q recorded).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClassifyReads_NoHits_IsRootWithZeroMatches()
    {
        var taxonomy = BuildClassFixture();
        var emptyDb = new Dictionary<string, int>();

        return Prop.ForAll(AcgtReadArbitrary(), seq =>
        {
            var r = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", seq) }, emptyDb, taxonomy, ClassK).Single();
            int expectedQ = CanonicalKmers(seq, ClassK).Count;
            return (r.TaxonId == TaxonomyTree.RootId && r.MatchedKmers == 0 && r.TotalKmers == expectedQ)
                .Label($"no-hit read: TaxonId={r.TaxonId}, C={r.MatchedKmers}, Q={r.TotalKmers} (expected Q={expectedQ})");
        });
    }

    /// <summary>
    /// INV-04 (edge): a read shorter than k (or empty) is unclassified with Q == 0 — there are no
    /// length-k windows to query. Source: doc §6.1 (empty / shorter-than-k read → root, Q = 0).
    /// </summary>
    [TestCase("")]
    [TestCase("A")]
    [TestCase("ACG")]
    [Category("Property")]
    public void ClassifyReads_ShorterThanK_IsRootWithZeroQ(string seq)
    {
        var taxonomy = BuildClassFixture();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(new[] { (100, RefA) }, taxonomy, ClassK);

        var r = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", seq) }, db, taxonomy, ClassK).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(TaxonomyTree.RootId), "shorter-than-k read must be root");
            Assert.That(r.TotalKmers, Is.EqualTo(0), "shorter-than-k read must have Q = 0");
            Assert.That(r.MatchedKmers, Is.EqualTo(0), "shorter-than-k read must have C = 0");
            Assert.That(r.Confidence, Is.EqualTo(0.0), "Q = 0 ⇒ Confidence = 0");
        });
    }

    /// <summary>
    /// Pure single-taxon read: DB built only from species A's non-repetitive reference; a read equal
    /// to that reference has every canonical k-mer unique to A ⇒ <c>TaxonId == 100</c> and
    /// <c>Confidence == 1.0</c> (C == Q). Source: doc §2.2 (RTL leaf = sole hit taxon; C/Q).
    /// </summary>
    [Test]
    [Category("Property")]
    public void ClassifyReads_PureSpeciesARead_AssignsSpeciesAWithFullConfidence()
    {
        var taxonomy = BuildClassFixture();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(new[] { (100, RefA) }, taxonomy, ClassK);

        var r = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", RefA) }, db, taxonomy, ClassK).Single();

        int expectedQ = CanonicalKmers(RefA, ClassK).Count;
        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(100), "exact species-A read must classify to species A");
            Assert.That(r.TotalKmers, Is.EqualTo(expectedQ), "Q must equal A's window count");
            Assert.That(r.MatchedKmers, Is.EqualTo(expectedQ), "all A k-mers lie in A's clade ⇒ C == Q");
            Assert.That(r.Confidence, Is.EqualTo(1.0).Within(1e-9), "C == Q ⇒ Confidence == 1.0");
            Assert.That(r.Species, Is.EqualTo("Species_A"));
            Assert.That(r.Genus, Is.EqualTo("Genus_X"));
        });
    }

    /// <summary>
    /// INV-06 (DB build, LCA): when two references — species A(100) and species B(101) — share a
    /// k-mer, that shared canonical k-mer maps to <c>Lca(100,101) == 10</c> (the genus); a k-mer
    /// unique to A maps to 100. Verified against the independent canonical/LCA oracles.
    /// Source: doc §2.2 / §2.4 INV-06 (k-mer → LCA of owning taxa).
    /// </summary>
    [Test]
    [Category("Property")]
    public void BuildKmerDatabase_SharedKmer_MapsToLcaOfOwners()
    {
        var taxonomy = BuildClassFixture();
        // RefShared is handed to BOTH A and B; RefA's own k-mers stay unique to A.
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, RefA), (100, RefShared), (101, RefShared) }, taxonomy, ClassK);

        int expectedLca = taxonomy.Lca(100, 101); // = 10 (genus)
        Assert.That(expectedLca, Is.EqualTo(10), "fixture sanity: Lca(A,B) = genus 10");

        // Every canonical k-mer of the shared reference must collapse to the genus LCA.
        foreach (string shared in CanonicalKmers(RefShared, ClassK).Distinct())
            Assert.That(db[shared], Is.EqualTo(expectedLca),
                $"shared k-mer {shared} must map to Lca(100,101)=10");

        // A k-mer unique to A (not present in RefShared) must remain species A.
        var sharedSet = CanonicalKmers(RefShared, ClassK).ToHashSet();
        string uniqueToA = CanonicalKmers(RefA, ClassK).First(km => !sharedSet.Contains(km));
        Assert.That(db[uniqueToA], Is.EqualTo(100), $"k-mer {uniqueToA} unique to A must map to 100");
    }

    /// <summary>
    /// LCA tie-break — sibling species (the key Kraken behavior): a DB in which one canonical k-mer
    /// is owned only by species A(100) and a DISJOINT canonical k-mer only by species B(101) yields a
    /// read with one hit to each. The two leaves tie at score 1, so the assignment is the LCA of the
    /// tied leaves = <c>Lca(100,101) == 10</c> (genus). Built deterministically; exact id asserted.
    /// Source: doc §2.2 / §6.1 (read split equally between siblings → genus).
    /// </summary>
    [Test]
    [Category("Property")]
    public void ClassifyReads_EqualSiblingHits_AssignsGenusLca()
    {
        var taxonomy = BuildClassFixture();
        // Hand-built DB: one k-mer → A, a different k-mer → B. No shared/LCA folding here.
        string kmerA = Canonical("ACGT");
        string kmerB = Canonical("GGCC");
        Assert.That(kmerA, Is.Not.EqualTo(kmerB), "fixture sanity: the two k-mers must differ");
        var db = new Dictionary<string, int> { [kmerA] = 100, [kmerB] = 101 };

        // Read containing exactly one occurrence of each k-mer (overlap-free, k = 4): "ACGTGGCC".
        string read = "ACGTGGCC";
        // Independent oracle: canonical k-mers of the read, count hits per taxon.
        var kmers = CanonicalKmers(read, ClassK);
        int hitsA = kmers.Count(km => km == kmerA);
        int hitsB = kmers.Count(km => km == kmerB);
        Assert.That(hitsA, Is.EqualTo(1).And.EqualTo(hitsB), "fixture sanity: exactly one hit to each sibling");

        var r = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", read) }, db, taxonomy, ClassK).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(taxonomy.Lca(100, 101)), "equal sibling hits ⇒ genus LCA");
            Assert.That(r.TaxonId, Is.EqualTo(10), "genus LCA id must be 10");
            // C = k-mers in the assigned clade (genus 10) = both hits; Q = 2 non-ambiguous windows ('GTGG' is the third).
            Assert.That(r.MatchedKmers, Is.EqualTo(2), "both hits lie in the genus clade ⇒ C = 2");
            Assert.That(r.TotalKmers, Is.EqualTo(CanonicalKmers(read, ClassK).Count));
            Assert.That(r.Confidence, Is.EqualTo(2.0 / CanonicalKmers(read, ClassK).Count).Within(1e-9));
        });
    }

    /// <summary>
    /// LCA tie-break — cross-genus split: one k-mer owned by species A(100) (genus 10) and another by
    /// species C(200) (genus 20) tie at score 1, so the assignment is <c>Lca(100,200) == 2</c>
    /// (the kingdom, the family-level LCA in this fixture). Source: doc §6.1 (split across genera →
    /// higher LCA).
    /// </summary>
    [Test]
    [Category("Property")]
    public void ClassifyReads_EqualCrossGenusHits_AssignsHigherLca()
    {
        var taxonomy = BuildClassFixture();
        string kmerA = Canonical("ACGT");
        string kmerC = Canonical("GGCC");
        var db = new Dictionary<string, int> { [kmerA] = 100, [kmerC] = 200 };

        string read = "ACGTGGCC";
        var r = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", read) }, db, taxonomy, ClassK).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(taxonomy.Lca(100, 200)), "equal cross-genus hits ⇒ kingdom LCA");
            Assert.That(r.TaxonId, Is.EqualTo(2), "cross-genus LCA id must be kingdom 2");
            Assert.That(r.MatchedKmers, Is.EqualTo(2), "both hits lie in the kingdom clade ⇒ C = 2");
        });
    }

    /// <summary>
    /// RTL maximum-weight path: unequal hits do NOT trigger the LCA tie-break — the deeper, higher-
    /// scoring leaf wins. Two A k-mers vs one B k-mer ⇒ species A(100). The C/Q confidence is
    /// computed by hand: C = 2 (A's clade), Q = read window count. Source: doc §2.2 (max-scoring RTL).
    /// </summary>
    [Test]
    [Category("Property")]
    public void ClassifyReads_UnequalSiblingHits_AssignsMajorityLeaf()
    {
        var taxonomy = BuildClassFixture();
        string kmerA = Canonical("ACGT");
        string kmerA2 = Canonical("TACA");
        string kmerB = Canonical("GGCC");
        Assert.That(new[] { kmerA, kmerA2, kmerB }.Distinct().Count(), Is.EqualTo(3), "fixture sanity: distinct k-mers");
        var db = new Dictionary<string, int> { [kmerA] = 100, [kmerA2] = 100, [kmerB] = 101 };

        // Read with two A windows and one B window. "ACGT" (→A), then "TACA" needs a window: build
        // explicitly and verify via the oracle rather than trusting the layout.
        string read = "ACGTACAGGCC"; // windows: ACGT, CGTA, GTAC, TACA, ACAG, CAGG, AGGC, GGCC
        var kmers = CanonicalKmers(read, ClassK);
        int hitsA = kmers.Count(km => km == kmerA) + kmers.Count(km => km == kmerA2);
        int hitsB = kmers.Count(km => km == kmerB);
        Assert.That(hitsA, Is.EqualTo(2), "fixture sanity: two A hits");
        Assert.That(hitsB, Is.EqualTo(1), "fixture sanity: one B hit");

        var r = MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", read) }, db, taxonomy, ClassK).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(100), "A scores 2 > B scores 1 ⇒ species A, no tie-break");
            Assert.That(r.MatchedKmers, Is.EqualTo(2), "C = A's clade k-mers = 2");
            Assert.That(r.TotalKmers, Is.EqualTo(kmers.Count));
            Assert.That(r.Confidence, Is.EqualTo(2.0 / kmers.Count).Within(1e-9));
        });
    }

    /// <summary>
    /// Exceptions (doc §3.3): null reads / db / taxonomy ⇒ <see cref="ArgumentNullException"/>;
    /// non-positive k ⇒ <see cref="ArgumentOutOfRangeException"/>. <c>ClassifyReads</c> validates
    /// eagerly (its body runs before the lazy iterator).
    /// </summary>
    [Test]
    [Category("Property")]
    public void ClassifyReads_InvalidArguments_Throw()
    {
        var taxonomy = BuildClassFixture();
        var db = new Dictionary<string, int>();
        var reads = new[] { ("r", "ACGT") };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() =>
                MetagenomicsAnalyzer.ClassifyReads(null!, db, taxonomy, ClassK));
            Assert.Throws<ArgumentNullException>(() =>
                MetagenomicsAnalyzer.ClassifyReads(reads, null!, taxonomy, ClassK));
            Assert.Throws<ArgumentNullException>(() =>
                MetagenomicsAnalyzer.ClassifyReads(reads, db, null!, ClassK));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                MetagenomicsAnalyzer.ClassifyReads(reads, db, taxonomy, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                MetagenomicsAnalyzer.ClassifyReads(reads, db, taxonomy, -1));
        });
    }

    /// <summary>
    /// Exceptions (doc §3.3): <c>BuildKmerDatabase</c> rejects null args, non-positive k, and a
    /// reference taxon absent from the tree (<see cref="KeyNotFoundException"/>).
    /// </summary>
    [Test]
    [Category("Property")]
    public void BuildKmerDatabase_InvalidArguments_Throw()
    {
        var taxonomy = BuildClassFixture();
        var refs = new[] { (100, RefA) };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() =>
                MetagenomicsAnalyzer.BuildKmerDatabase(null!, taxonomy, ClassK));
            Assert.Throws<ArgumentNullException>(() =>
                MetagenomicsAnalyzer.BuildKmerDatabase(refs, null!, ClassK));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                MetagenomicsAnalyzer.BuildKmerDatabase(refs, taxonomy, 0));
            // Reference taxon 999 is not in the tree.
            Assert.Throws<KeyNotFoundException>(() =>
                MetagenomicsAnalyzer.BuildKmerDatabase(new[] { (999, RefA) }, taxonomy, ClassK));
        });
    }

    #endregion

    #region META-PROF-001: R: abundance ∈ [0,1]; Kingdom sums to 1 (lower ranks ≤ 1, exact shortfall); INV-01..04; D: deterministic

    // META-PROF-001 — Taxonomic profile generation (doc: docs/algorithms/Metagenomics/Taxonomic_Profile.md).
    //
    // GenerateTaxonomicProfile aggregates per-read TaxonomicClassification records into a community
    // profile. The tests below treat the doc §2.2 / §2.4 contract as the specification and verify the
    // implementation against INDEPENDENT oracles recomputed from the generated input:
    //   * INV-01  TotalReads     = input.Count.
    //   * INV-02  ClassifiedReads = #reads whose Kingdom is NEITHER "Unclassified" NOR empty/whitespace.
    //   * R       every abundance value in all four maps ∈ [0,1].
    //   * sum     KingdomAbundance sums to exactly 1 when ClassifiedReads>0 (else empty);
    //             Phylum/Genus/Species each sum to ≤ 1, and == 1 IFF every retained read has a
    //             non-empty value at that rank — with the exact shortfall == emptyCount/ClassifiedReads.
    //   * oracle  abundance[label] = (#retained reads with that non-empty label) / ClassifiedReads.
    //   * INV-04  Shannon = −Σ p·ln p, Simpson = Σ p² over the SPECIES distribution (natural log).
    //   * D       identical input ⇒ identical profile.
    //
    // The generator draws Kingdom/Phylum/Genus/Species from small label sets that INCLUDE the special
    // values "Unclassified" and "" so both the unclassified filter (Kingdom) and the empty-rank skip
    // (lower ranks) are exercised. Non-profiling fields are harmless defaults.

    /// <summary>Label pool for the kingdom rank — includes the two filtered-out specials.</summary>
    private static readonly string[] ProfKingdoms = { "Bacteria", "Archaea", "Unclassified", "" };

    /// <summary>Label pool for the lower ranks — includes "" so empty-rank skipping is exercised.</summary>
    private static readonly string[] ProfPhyla = { "Firmicutes", "Proteobacteria", "" };
    private static readonly string[] ProfGenera = { "Escherichia", "Bacillus", "Clostridium", "" };
    private static readonly string[] ProfSpecies = { "E_coli", "B_subtilis", "C_difficile", "S_aureus", "" };

    /// <summary>
    /// Builds a <see cref="MetagenomicsAnalyzer.TaxonomicClassification"/> from rank labels alone,
    /// filling the (irrelevant-to-profiling) classifier fields with harmless deterministic defaults.
    /// </summary>
    private static MetagenomicsAnalyzer.TaxonomicClassification ProfRead(
        int idx, string kingdom, string phylum, string genus, string species) =>
        new(
            ReadId: $"read_{idx}",
            TaxonId: 0,
            TaxonName: "",
            Rank: "",
            RtlScore: 0,
            Confidence: 0.0,
            MatchedKmers: 0,
            TotalKmers: 0,
            Kingdom: kingdom,
            Phylum: phylum,
            Class: "",
            Order: "",
            Family: "",
            Genus: genus,
            Species: species);

    /// <summary>Generates a single read's four rank labels, each drawn from its pool (which includes specials).</summary>
    private static Gen<(string K, string P, string G, string S)> ProfReadGen() =>
        Gen.Elements(ProfKingdoms).SelectMany(k =>
            Gen.Elements(ProfPhyla).SelectMany(ph =>
                Gen.Elements(ProfGenera).SelectMany(g =>
                    Gen.Elements(ProfSpecies).Select(sp => (k, ph, g, sp)))));

    /// <summary>Generates a (possibly empty) list of classifications drawing each rank from the pools above.</summary>
    private static Arbitrary<List<MetagenomicsAnalyzer.TaxonomicClassification>> ProfClassificationsArbitrary() =>
        Gen.Choose(0, 30)
            .SelectMany(n =>
                ProfReadGen().ArrayOf(n)
                    .Select(tuples =>
                    {
                        var list = new List<MetagenomicsAnalyzer.TaxonomicClassification>(tuples.Length);
                        for (int i = 0; i < tuples.Length; i++)
                            list.Add(ProfRead(i, tuples[i].K, tuples[i].P, tuples[i].G, tuples[i].S));
                        return list;
                    }))
            .ToArbitrary();

    /// <summary>Independent oracle: retained (classified) reads = Kingdom not "Unclassified" and not blank.</summary>
    private static List<MetagenomicsAnalyzer.TaxonomicClassification> Retained(
        IEnumerable<MetagenomicsAnalyzer.TaxonomicClassification> reads) =>
        reads.Where(c => c.Kingdom != "Unclassified" && !string.IsNullOrWhiteSpace(c.Kingdom)).ToList();

    /// <summary>
    /// Independent oracle for a single rank's abundance map: count non-empty labels among retained
    /// reads, divide by ClassifiedReads. Empty/blank labels are excluded from the map.
    /// </summary>
    private static Dictionary<string, double> ExpectedAbundance(
        IReadOnlyList<MetagenomicsAnalyzer.TaxonomicClassification> retained,
        Func<MetagenomicsAnalyzer.TaxonomicClassification, string> rank)
    {
        int classified = retained.Count;
        var map = new Dictionary<string, double>();
        if (classified == 0) return map;
        foreach (var label in retained.Select(rank).Where(s => !string.IsNullOrEmpty(s)))
            map[label] = map.GetValueOrDefault(label) + 1.0;
        foreach (var key in map.Keys.ToList())
            map[key] /= classified;
        return map;
    }

    private static bool MapsApproxEqual(IReadOnlyDictionary<string, double> a, IReadOnlyDictionary<string, double> b)
    {
        if (a.Count != b.Count) return false;
        foreach (var kv in a)
            if (!b.TryGetValue(kv.Key, out double v) || Math.Abs(v - kv.Value) > 1e-9) return false;
        return true;
    }

    /// <summary>
    /// INV-01: <c>TotalReads</c> equals the number of input classification records (counted before
    /// any filtering). Source: doc §2.4 INV-01.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateTaxonomicProfile_TotalReads_EqualsInputCount()
    {
        return Prop.ForAll(ProfClassificationsArbitrary(), reads =>
        {
            var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
            return (p.TotalReads == reads.Count)
                .Label($"TotalReads={p.TotalReads}, expected {reads.Count}");
        });
    }

    /// <summary>
    /// INV-02: <c>ClassifiedReads</c> equals the number of reads whose Kingdom is neither
    /// "Unclassified" nor empty/whitespace — the unclassified filter is on Kingdom only.
    /// Source: doc §2.4 INV-02 / §3.3.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateTaxonomicProfile_ClassifiedReads_EqualsKingdomFilteredCount()
    {
        return Prop.ForAll(ProfClassificationsArbitrary(), reads =>
        {
            var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
            int expected = Retained(reads).Count;
            return (p.ClassifiedReads == expected)
                .Label($"ClassifiedReads={p.ClassifiedReads}, expected {expected}");
        });
    }

    /// <summary>
    /// R (abundance ∈ [0,1]): every value in all four abundance maps lies in [0,1].
    /// Evidence: each is a non-empty-label count divided by ClassifiedReads ≥ that count.
    /// Source: doc §2.2 (abundance = count/Σcount).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateTaxonomicProfile_AllAbundances_InUnitInterval()
    {
        return Prop.ForAll(ProfClassificationsArbitrary(), reads =>
        {
            var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
            bool ok = new[] { p.KingdomAbundance, p.PhylumAbundance, p.GenusAbundance, p.SpeciesAbundance }
                .SelectMany(m => m.Values)
                .All(v => v >= -1e-12 && v <= 1.0 + 1e-12);
            return ok.Label("every abundance value must be in [0,1]");
        });
    }

    /// <summary>
    /// Abundance value oracle: for every rank map, <c>map[label]</c> equals
    /// (#retained reads with that non-empty label) / ClassifiedReads, recomputed independently.
    /// Source: doc §2.2 (abundance = count / Σcount) with denominator = ClassifiedReads.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateTaxonomicProfile_AbundanceMaps_MatchCountOverClassifiedOracle()
    {
        return Prop.ForAll(ProfClassificationsArbitrary(), reads =>
        {
            var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
            var ret = Retained(reads);
            bool kOk = MapsApproxEqual(p.KingdomAbundance, ExpectedAbundance(ret, c => c.Kingdom));
            bool phOk = MapsApproxEqual(p.PhylumAbundance, ExpectedAbundance(ret, c => c.Phylum));
            bool gOk = MapsApproxEqual(p.GenusAbundance, ExpectedAbundance(ret, c => c.Genus));
            bool spOk = MapsApproxEqual(p.SpeciesAbundance, ExpectedAbundance(ret, c => c.Species));
            return (kOk && phOk && gOk && spOk)
                .Label("abundance maps must equal count/ClassifiedReads oracle for every rank");
        });
    }

    /// <summary>
    /// INV-03 (Kingdom): <c>KingdomAbundance</c> sums to exactly 1 when ClassifiedReads &gt; 0, and is
    /// empty (sum 0) when ClassifiedReads == 0. Kingdom labels of retained reads are never empty (the
    /// filter removes blank kingdoms), so no kingdom count is skipped. Source: doc §2.4 INV-03.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateTaxonomicProfile_KingdomAbundance_SumsToExactlyOneWhenClassified()
    {
        return Prop.ForAll(ProfClassificationsArbitrary(), reads =>
        {
            var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
            double sum = p.KingdomAbundance.Values.Sum();
            bool ok = p.ClassifiedReads > 0
                ? Math.Abs(sum - 1.0) < 1e-9
                : p.KingdomAbundance.Count == 0;
            return ok.Label($"ClassifiedReads={p.ClassifiedReads}, kingdom sum={sum:F12}");
        });
    }

    /// <summary>
    /// INV-03 (lower ranks, exact shortfall): each of Phylum/Genus/Species sums to ≤ 1, and the exact
    /// shortfall from 1 equals (#retained reads with an EMPTY value at that rank) / ClassifiedReads —
    /// because empty rank strings are skipped while the denominator stays ClassifiedReads. Hence the
    /// sum is exactly 1 IFF every retained read has a non-empty value at that rank. Both directions are
    /// covered: random inputs include reads with empty and non-empty lower ranks. Source: doc §2.4
    /// INV-03 / §6.1 (missing rank ⇒ sum &lt; 1).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateTaxonomicProfile_LowerRankSums_HaveExactEmptyShortfall()
    {
        return Prop.ForAll(ProfClassificationsArbitrary(), reads =>
        {
            var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
            var ret = Retained(reads);
            int classified = ret.Count;

            bool CheckRank(IReadOnlyDictionary<string, double> map,
                Func<MetagenomicsAnalyzer.TaxonomicClassification, string> rank)
            {
                double sum = map.Values.Sum();
                if (sum > 1.0 + 1e-9) return false;
                if (classified == 0) return map.Count == 0 && Math.Abs(sum) < 1e-9;
                int empties = ret.Count(c => string.IsNullOrEmpty(rank(c)));
                double expectedSum = (classified - empties) / (double)classified;
                bool shortfallOk = Math.Abs(sum - expectedSum) < 1e-9;
                // IFF direction: sum == 1 exactly when there are no empty values at this rank.
                bool iffOk = (Math.Abs(sum - 1.0) < 1e-9) == (empties == 0);
                return shortfallOk && iffOk;
            }

            bool ok = CheckRank(p.PhylumAbundance, c => c.Phylum)
                      && CheckRank(p.GenusAbundance, c => c.Genus)
                      && CheckRank(p.SpeciesAbundance, c => c.Species);
            return ok.Label("lower-rank sum must be 1 minus emptyCount/ClassifiedReads");
        });
    }

    /// <summary>
    /// INV-04 + diversity oracle: <c>ShannonDiversity</c> = −Σ pᵢ ln pᵢ and <c>SimpsonDiversity</c> =
    /// Σ pᵢ² over the SPECIES distribution, where pᵢ = (species count) / (Σ species counts), recomputed
    /// independently with the natural log. Source: doc §2.2 (H, D) / §2.4 INV-04.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateTaxonomicProfile_Diversity_MatchesSpeciesOracle()
    {
        return Prop.ForAll(ProfClassificationsArbitrary(), reads =>
        {
            var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
            var ret = Retained(reads);

            // Independent species-count distribution (non-empty species among retained reads).
            var counts = ret.Select(c => c.Species)
                .Where(s => !string.IsNullOrEmpty(s))
                .GroupBy(s => s)
                .Select(g => (double)g.Count())
                .ToList();
            double sum = counts.Sum();

            double expectedShannon = 0.0, expectedSimpson = 0.0;
            if (sum > 0)
            {
                foreach (double c in counts)
                {
                    double pi = c / sum;
                    expectedShannon -= pi * Math.Log(pi);
                    expectedSimpson += pi * pi;
                }
            }

            bool ok = Math.Abs(p.ShannonDiversity - expectedShannon) < 1e-9
                      && Math.Abs(p.SimpsonDiversity - expectedSimpson) < 1e-9;
            return ok.Label(
                $"Shannon={p.ShannonDiversity:F9} (exp {expectedShannon:F9}), " +
                $"Simpson={p.SimpsonDiversity:F9} (exp {expectedSimpson:F9})");
        });
    }

    /// <summary>
    /// D (determinism): identical input produces an identical profile — all four maps and every scalar.
    /// Evidence: <c>GenerateTaxonomicProfile</c> is a pure function of its input. Source: doc INV-D.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateTaxonomicProfile_IsDeterministic()
    {
        return Prop.ForAll(ProfClassificationsArbitrary(), reads =>
        {
            var a = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
            var b = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
            bool ok = a.TotalReads == b.TotalReads
                      && a.ClassifiedReads == b.ClassifiedReads
                      && Math.Abs(a.ShannonDiversity - b.ShannonDiversity) < 1e-12
                      && Math.Abs(a.SimpsonDiversity - b.SimpsonDiversity) < 1e-12
                      && MapsApproxEqual(a.KingdomAbundance, b.KingdomAbundance)
                      && MapsApproxEqual(a.PhylumAbundance, b.PhylumAbundance)
                      && MapsApproxEqual(a.GenusAbundance, b.GenusAbundance)
                      && MapsApproxEqual(a.SpeciesAbundance, b.SpeciesAbundance);
            return ok.Label("GenerateTaxonomicProfile must be deterministic");
        });
    }

    /// <summary>
    /// Edge (doc §6.1): empty input ⇒ TotalReads=0, ClassifiedReads=0, all maps empty, zero diversity.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GenerateTaxonomicProfile_EmptyInput_IsZeroProfile()
    {
        var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(
            Array.Empty<MetagenomicsAnalyzer.TaxonomicClassification>());

        Assert.Multiple(() =>
        {
            Assert.That(p.TotalReads, Is.EqualTo(0));
            Assert.That(p.ClassifiedReads, Is.EqualTo(0));
            Assert.That(p.KingdomAbundance, Is.Empty);
            Assert.That(p.PhylumAbundance, Is.Empty);
            Assert.That(p.GenusAbundance, Is.Empty);
            Assert.That(p.SpeciesAbundance, Is.Empty);
            Assert.That(p.ShannonDiversity, Is.EqualTo(0.0));
            Assert.That(p.SimpsonDiversity, Is.EqualTo(0.0));
        });
    }

    /// <summary>
    /// Edge (doc §6.1): all-unclassified input (Kingdom "Unclassified" or empty) ⇒ ClassifiedReads=0,
    /// all maps empty, zero diversity, while TotalReads still counts every input record.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GenerateTaxonomicProfile_AllUnclassified_HasNoClassifiedReads()
    {
        var reads = new[]
        {
            ProfRead(0, "Unclassified", "Firmicutes", "Bacillus", "B_subtilis"),
            ProfRead(1, "", "Proteobacteria", "Escherichia", "E_coli"),
            ProfRead(2, "Unclassified", "", "", ""),
        };

        var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);

        Assert.Multiple(() =>
        {
            Assert.That(p.TotalReads, Is.EqualTo(3));
            Assert.That(p.ClassifiedReads, Is.EqualTo(0));
            Assert.That(p.KingdomAbundance, Is.Empty);
            Assert.That(p.PhylumAbundance, Is.Empty);
            Assert.That(p.GenusAbundance, Is.Empty);
            Assert.That(p.SpeciesAbundance, Is.Empty);
            Assert.That(p.ShannonDiversity, Is.EqualTo(0.0));
            Assert.That(p.SimpsonDiversity, Is.EqualTo(0.0));
        });
    }

    /// <summary>
    /// Anchor (doc §6.1): a single distinct species after filtering ⇒ ShannonDiversity = 0,
    /// SimpsonDiversity = 1 (one-category distribution). Three reads of the same species; the kingdom
    /// map sums to 1 and the species map has a single entry equal to 1.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GenerateTaxonomicProfile_SingleSpecies_ShannonZeroSimpsonOne()
    {
        var reads = new[]
        {
            ProfRead(0, "Bacteria", "Firmicutes", "Bacillus", "B_subtilis"),
            ProfRead(1, "Bacteria", "Firmicutes", "Bacillus", "B_subtilis"),
            ProfRead(2, "Bacteria", "Firmicutes", "Bacillus", "B_subtilis"),
        };

        var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);

        Assert.Multiple(() =>
        {
            Assert.That(p.ClassifiedReads, Is.EqualTo(3));
            Assert.That(p.ShannonDiversity, Is.EqualTo(0.0).Within(1e-9), "single species ⇒ H = 0");
            Assert.That(p.SimpsonDiversity, Is.EqualTo(1.0).Within(1e-9), "single species ⇒ D = 1");
            Assert.That(p.KingdomAbundance.Values.Sum(), Is.EqualTo(1.0).Within(1e-9));
            Assert.That(p.SpeciesAbundance.Count, Is.EqualTo(1));
            Assert.That(p.SpeciesAbundance["B_subtilis"], Is.EqualTo(1.0).Within(1e-9));
        });
    }

    /// <summary>
    /// Lower-rank shortfall (concrete, doc §6.1): four classified reads, two with empty species. The
    /// species map must sum to exactly 1 − 2/4 = 0.5 while the kingdom map still sums to 1.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GenerateTaxonomicProfile_SomeEmptySpecies_SpeciesSumHasExactShortfall()
    {
        var reads = new[]
        {
            ProfRead(0, "Bacteria", "Firmicutes", "Bacillus", "B_subtilis"),
            ProfRead(1, "Bacteria", "Firmicutes", "Clostridium", "C_difficile"),
            ProfRead(2, "Bacteria", "Proteobacteria", "Escherichia", ""),
            ProfRead(3, "Bacteria", "Proteobacteria", "Escherichia", ""),
        };

        var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);

        Assert.Multiple(() =>
        {
            Assert.That(p.ClassifiedReads, Is.EqualTo(4));
            Assert.That(p.KingdomAbundance.Values.Sum(), Is.EqualTo(1.0).Within(1e-9),
                "kingdom labels are never empty for retained reads ⇒ sum 1");
            Assert.That(p.SpeciesAbundance.Values.Sum(), Is.EqualTo(0.5).Within(1e-9),
                "2 of 4 retained reads have empty species ⇒ shortfall 2/4 ⇒ sum 0.5");
            Assert.That(p.GenusAbundance.Values.Sum(), Is.EqualTo(1.0).Within(1e-9),
                "all retained reads have a non-empty genus ⇒ genus sum 1");
        });
    }

    #endregion

    #region META-BIN-001: R: each contig in ≤ 1 bin; P: bin GC% consistent within bin; D: deterministic

    // META-BIN-001 — Genome binning (doc: docs/algorithms/Metagenomics/Genome_Binning.md).
    //
    // BinContigs groups contigs into bins (MAGs) via deterministic k-means over composition/coverage.
    // k-means cluster membership cannot be predicted easily, so these tests are HEURISTIC-INDEPENDENT:
    // for WHATEVER bins are produced, every reported field must be CONSISTENT with that bin's reported
    // ContigIds, recomputed independently from the input contig map. The oracles below are derived from
    // the doc §2.2/§4.2 formulas and the source contract, NOT from the code's output:
    //   * TotalLength   = Σ member sequence lengths.
    //   * GcContent     = mean of per-member GC FRACTION = (#G+#C)/len over uppercase ACGT.
    //   * Coverage      = mean of raw member coverage.
    //   * Completeness  = min(TotalLength / expectedGenomeSize · 100, 100).
    //   * Contamination = 0 if <2 members, else min( popStdDev(member GC fractions) / 0.5 · 100, 100 ),
    //                     where popStdDev uses POPULATION variance (divide by count, not count-1).
    //   * PredictedTaxonomy == "".
    // INV-01 (R): contig-disjointness — the multiset of all ContigIds across bins has no duplicate, and
    //             every id is an input id (some inputs may appear in 0 bins; >1 is forbidden).
    // INV-02:     Completeness, Contamination ∈ [0,100]. INV-03: every emitted bin TotalLength ≥ minBinSize.
    // D:          identical input ⇒ identical bin sequence (GC-sorted centroid init ⇒ deterministic).
    //
    // Defaults (minBinSize=500000, expectedGenomeSize=4_000_000) are too large for fast tests, so the
    // tests pass SMALL parameters and short uppercase A/C/G/T contigs so bins are actually emitted.

    /// <summary>Test-scale binning parameters: small so bins are emitted for short test contigs.</summary>
    private const double BinMinBinSize = 50.0;
    private const double BinExpectedGenomeSize = 1000.0;

    /// <summary>Independent GC-fraction oracle = (#G+#C)/len over uppercase ACGT (mirrors CalculateGcFractionFast).</summary>
    private static double GcFractionOracle(string seq)
    {
        if (string.IsNullOrEmpty(seq)) return 0.0;
        int gc = 0, valid = 0;
        foreach (char c in seq)
        {
            switch (c)
            {
                case 'G' or 'C': gc++; valid++; break;
                case 'A' or 'T': valid++; break;
            }
        }
        return valid == 0 ? 0.0 : (double)gc / valid;
    }

    /// <summary>Independent contamination oracle: population GC std-dev / 0.5 · 100, clamped to 100; 0 if &lt;2 members.</summary>
    private static double ContaminationOracle(IReadOnlyList<double> gcFractions)
    {
        if (gcFractions.Count < 2) return 0.0;
        double mean = gcFractions.Average();
        double variance = gcFractions.Sum(g => (g - mean) * (g - mean)) / gcFractions.Count; // population variance
        double stdDev = Math.Sqrt(variance);
        return Math.Min(stdDev / 0.5 * 100.0, 100.0);
    }

    /// <summary>
    /// Generates a small list of contigs with UNIQUE ids, uppercase ACGT sequences (length 20–120),
    /// and finite non-negative coverage. Pure ACGT so the GC-fraction oracle matches the implementation.
    /// </summary>
    private static Arbitrary<List<(string ContigId, string Sequence, double Coverage)>> BinContigsArbitrary() =>
        Gen.Choose(1, 12)
            .SelectMany(n =>
                BinContigPartGen().ArrayOf(n)
                    .Select(parts =>
                    {
                        var list = new List<(string, string, double)>(n);
                        for (int i = 0; i < n; i++)
                            list.Add(($"contig_{i}", parts[i].Seq, parts[i].Cov));
                        return list;
                    }))
            .ToArbitrary();

    /// <summary>Generates one contig's (sequence, coverage): uppercase ACGT of length 20–120, coverage ≥ 0.</summary>
    private static Gen<(string Seq, double Cov)> BinContigPartGen() =>
        Gen.Choose(20, 120).SelectMany(len =>
            Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len).SelectMany(cs =>
                Gen.Choose(0, 1000).Select(cov => (new string(cs), cov / 10.0))));

    /// <summary>
    /// P (GC consistent) + derived-field oracle: for EVERY output bin, all reported numeric fields are
    /// consistent with that bin's reported ContigIds recomputed from the input map — GcContent = mean
    /// member GC fraction, TotalLength = Σ member lengths, Coverage = mean member coverage,
    /// Completeness = min(TotalLength/expectedGenomeSize·100,100), Contamination = population GC-stddev
    /// oracle. Source: doc §2.2 / §3.2 / §4.2.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BinContigs_ReportedFields_MatchMemberOracle()
    {
        return Prop.ForAll(BinContigsArbitrary(), contigs =>
        {
            var map = contigs.ToDictionary(c => c.ContigId, c => c);
            var bins = MetagenomicsAnalyzer.BinContigs(
                contigs, numBins: 10, minBinSize: BinMinBinSize, expectedGenomeSize: BinExpectedGenomeSize).ToList();

            foreach (var bin in bins)
            {
                var members = bin.ContigIds.Select(id => map[id]).ToList();
                var gcFractions = members.Select(m => GcFractionOracle(m.Sequence)).ToList();

                double expTotalLength = members.Sum(m => (double)m.Sequence.Length);
                double expGc = gcFractions.Average();
                double expCoverage = members.Average(m => m.Coverage);
                double expCompleteness = Math.Min(expTotalLength / BinExpectedGenomeSize * 100.0, 100.0);
                double expContamination = ContaminationOracle(gcFractions);

                if (Math.Abs(bin.TotalLength - expTotalLength) > 1e-9 ||
                    Math.Abs(bin.GcContent - expGc) > 1e-9 ||
                    Math.Abs(bin.Coverage - expCoverage) > 1e-9 ||
                    Math.Abs(bin.Completeness - expCompleteness) > 1e-9 ||
                    Math.Abs(bin.Contamination - expContamination) > 1e-9)
                {
                    return false.Label(
                        $"bin {bin.BinId}: len {bin.TotalLength} vs {expTotalLength}, gc {bin.GcContent:F12} vs {expGc:F12}, " +
                        $"cov {bin.Coverage:F12} vs {expCoverage:F12}, comp {bin.Completeness:F12} vs {expCompleteness:F12}, " +
                        $"contam {bin.Contamination:F12} vs {expContamination:F12}");
                }
            }
            return true.Label("all reported bin fields must match the member-derived oracle");
        });
    }

    /// <summary>
    /// INV-01 (R: each contig in ≤ 1 bin): the multiset of all ContigIds across all output bins has NO
    /// duplicate contig id, and every reported id is a member of the input set. Source: doc §2.4 INV-01
    /// (k-means assigns each contig to exactly one cluster ⇒ disjoint bins).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BinContigs_ContigIds_AreDisjointAndFromInput()
    {
        return Prop.ForAll(BinContigsArbitrary(), contigs =>
        {
            var inputIds = contigs.Select(c => c.ContigId).ToHashSet();
            var bins = MetagenomicsAnalyzer.BinContigs(
                contigs, numBins: 10, minBinSize: BinMinBinSize, expectedGenomeSize: BinExpectedGenomeSize).ToList();

            var all = bins.SelectMany(b => b.ContigIds).ToList();
            bool noDuplicate = all.Count == all.Distinct().Count();
            bool allFromInput = all.All(inputIds.Contains);
            return (noDuplicate && allFromInput)
                .Label($"contig ids must be disjoint across bins and ⊆ input: total={all.Count}, distinct={all.Distinct().Count()}");
        });
    }

    /// <summary>
    /// INV-02: every output bin has Completeness ∈ [0,100] and Contamination ∈ [0,100], because both
    /// proxy formulas clamp via Math.Min(..., 100) and are non-negative. Source: doc §2.4 INV-02.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BinContigs_CompletenessAndContamination_InRange()
    {
        return Prop.ForAll(BinContigsArbitrary(), contigs =>
        {
            var bins = MetagenomicsAnalyzer.BinContigs(
                contigs, numBins: 10, minBinSize: BinMinBinSize, expectedGenomeSize: BinExpectedGenomeSize).ToList();
            bool ok = bins.All(b =>
                b.Completeness >= -1e-9 && b.Completeness <= 100.0 + 1e-9 &&
                b.Contamination >= -1e-9 && b.Contamination <= 100.0 + 1e-9);
            return ok.Label("Completeness and Contamination must lie in [0,100] for every bin");
        });
    }

    /// <summary>
    /// INV-03: every emitted bin has TotalLength ≥ minBinSize — sub-threshold clusters are filtered out.
    /// Source: doc §2.4 INV-03 / §4.1 step 7.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BinContigs_EmittedBins_MeetMinBinSize()
    {
        return Prop.ForAll(BinContigsArbitrary(), contigs =>
        {
            var bins = MetagenomicsAnalyzer.BinContigs(
                contigs, numBins: 10, minBinSize: BinMinBinSize, expectedGenomeSize: BinExpectedGenomeSize).ToList();
            return bins.All(b => b.TotalLength >= BinMinBinSize - 1e-9)
                .Label($"every emitted bin must have TotalLength ≥ {BinMinBinSize}");
        });
    }

    /// <summary>
    /// PredictedTaxonomy: the current implementation emits an empty string for every bin (no taxonomy
    /// assignment). Source: doc §3.2 / §5.2.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BinContigs_PredictedTaxonomy_IsEmpty()
    {
        return Prop.ForAll(BinContigsArbitrary(), contigs =>
        {
            var bins = MetagenomicsAnalyzer.BinContigs(
                contigs, numBins: 10, minBinSize: BinMinBinSize, expectedGenomeSize: BinExpectedGenomeSize).ToList();
            return bins.All(b => b.PredictedTaxonomy == "")
                .Label("every bin's PredictedTaxonomy must be the empty string");
        });
    }

    /// <summary>
    /// D (determinism): identical input ⇒ identical bin sequence — BinId order, ContigIds, and all
    /// numeric fields. k-means here uses deterministic GC-sorted centroid initialization. Source: doc
    /// §5.2 (deterministic centroid init).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BinContigs_IsDeterministic()
    {
        return Prop.ForAll(BinContigsArbitrary(), contigs =>
        {
            var r1 = MetagenomicsAnalyzer.BinContigs(
                contigs, numBins: 10, minBinSize: BinMinBinSize, expectedGenomeSize: BinExpectedGenomeSize).ToList();
            var r2 = MetagenomicsAnalyzer.BinContigs(
                contigs, numBins: 10, minBinSize: BinMinBinSize, expectedGenomeSize: BinExpectedGenomeSize).ToList();

            if (r1.Count != r2.Count) return false.Label($"bin count differs: {r1.Count} vs {r2.Count}");
            for (int i = 0; i < r1.Count; i++)
            {
                bool same = r1[i].BinId == r2[i].BinId
                            && r1[i].ContigIds.SequenceEqual(r2[i].ContigIds)
                            && Math.Abs(r1[i].TotalLength - r2[i].TotalLength) < 1e-12
                            && Math.Abs(r1[i].GcContent - r2[i].GcContent) < 1e-12
                            && Math.Abs(r1[i].Coverage - r2[i].Coverage) < 1e-12
                            && Math.Abs(r1[i].Completeness - r2[i].Completeness) < 1e-12
                            && Math.Abs(r1[i].Contamination - r2[i].Contamination) < 1e-12
                            && r1[i].PredictedTaxonomy == r2[i].PredictedTaxonomy;
                if (!same) return false.Label($"bin {i} differs between identical runs");
            }
            return true.Label("BinContigs must be deterministic for identical input");
        });
    }

    /// <summary>
    /// Edge (doc §6.1): empty input collection ⇒ no bins are returned.
    /// </summary>
    [Test]
    [Category("Property")]
    public void BinContigs_EmptyInput_ReturnsNoBins()
    {
        var bins = MetagenomicsAnalyzer.BinContigs(
            Array.Empty<(string, string, double)>(),
            numBins: 10, minBinSize: BinMinBinSize, expectedGenomeSize: BinExpectedGenomeSize).ToList();

        Assert.That(bins, Is.Empty, "empty input must yield no bins");
    }

    /// <summary>
    /// Edge anchor: a single contig of length ≥ minBinSize ⇒ exactly one bin containing it, with
    /// Contamination == 0 (only one member ⇒ &lt;2 GC values) and GcContent == its GC fraction.
    /// Source: doc §6.1 (single-contig bin ⇒ Contamination = 0).
    /// </summary>
    [Test]
    [Category("Property")]
    public void BinContigs_SingleContig_OneBinZeroContaminationGcFraction()
    {
        // 60 bp, 30 G/C ⇒ GC fraction = 0.5, length 60 ≥ minBinSize 50.
        const string seq = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATAT";
        var contigs = new[] { ("contig_solo", seq, 12.5) };

        var bins = MetagenomicsAnalyzer.BinContigs(
            contigs, numBins: 10, minBinSize: BinMinBinSize, expectedGenomeSize: BinExpectedGenomeSize).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(bins, Has.Count.EqualTo(1), "one above-threshold contig ⇒ exactly one bin");
            Assert.That(bins[0].ContigIds, Is.EquivalentTo(new[] { "contig_solo" }));
            Assert.That(bins[0].Contamination, Is.EqualTo(0.0), "single member ⇒ Contamination = 0");
            Assert.That(bins[0].GcContent, Is.EqualTo(GcFractionOracle(seq)).Within(1e-9));
            Assert.That(bins[0].GcContent, Is.EqualTo(0.5).Within(1e-9), "30 G/C of 60 bp ⇒ GC fraction 0.5");
            Assert.That(bins[0].TotalLength, Is.EqualTo(seq.Length));
            Assert.That(bins[0].Coverage, Is.EqualTo(12.5).Within(1e-9));
            Assert.That(bins[0].PredictedTaxonomy, Is.EqualTo(""));
        });
    }

    /// <summary>
    /// Edge (doc §6.1): numBins larger than the contig count is accepted — k is clamped to the contig
    /// count, so two distinct above-threshold contigs still produce at most two bins covering both ids.
    /// </summary>
    [Test]
    [Category("Property")]
    public void BinContigs_NumBinsLargerThanContigCount_IsAccepted()
    {
        const string seqA = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"; // high GC
        const string seqB = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT"; // low GC
        var contigs = new[] { ("contig_a", seqA, 5.0), ("contig_b", seqB, 50.0) };

        var bins = MetagenomicsAnalyzer.BinContigs(
            contigs, numBins: 100, minBinSize: BinMinBinSize, expectedGenomeSize: BinExpectedGenomeSize).ToList();

        var coveredIds = bins.SelectMany(b => b.ContigIds).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(bins.Count, Is.GreaterThanOrEqualTo(1).And.LessThanOrEqualTo(2),
                "k clamped to contig count ⇒ at most two bins");
            Assert.That(coveredIds, Is.Unique, "no contig appears in more than one bin");
            Assert.That(coveredIds, Is.EquivalentTo(new[] { "contig_a", "contig_b" }),
                "both above-threshold contigs are binned");
        });
    }

    #endregion

    #region META-FUNC-001: R: function scores ≥ 0; P: assigned function in DB (signature in gene); D: deterministic

    // PredictFunctions assigns each gene its best-hit (lowest E-value) function whose DB signature is
    // a substring of the gene's protein sequence (BLAST-style ranking).

    private static readonly IReadOnlyDictionary<string, (string Function, string Pathway, string Ko)> FunctionDb =
        new Dictionary<string, (string, string, string)>
        {
            ["GGG"] = ("Glycine-rich", "P1", "K001"),
            ["WWY"] = ("Aromatic", "P2", "K002"),
            ["MKMK"] = ("Initiator-like", "P3", "K003"),
        };

    /// <summary>Genes whose sequences embed a random DB signature inside random protein flanks.</summary>
    private static Arbitrary<(string GeneId, string Seq)[]> FunctionGenesArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            const string aa = "ACDEFGHIKLMNPQRSTVWY";
            string Rand(int n)
            {
                var c = new char[n];
                for (int i = 0; i < n; i++) c[i] = aa[rng.Next(aa.Length)];
                return new string(c);
            }
            var sigs = FunctionDb.Keys.ToList();
            int n = 1 + rng.Next(4);
            var genes = new (string, string)[n];
            for (int i = 0; i < n; i++)
            {
                string sig = sigs[rng.Next(sigs.Count)];
                genes[i] = ($"g{i}", Rand(rng.Next(5)) + sig + Rand(rng.Next(5)));
            }
            return genes;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R + P): every annotation has non-negative scores and an assigned function that comes
    /// from a DB entry whose signature is actually contained in that gene's sequence.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Functions_AreInDbAndScoresNonNegative()
    {
        return Prop.ForAll(FunctionGenesArbitrary(), genes =>
        {
            var byId = genes.ToDictionary(g => g.GeneId, g => g.Seq);
            var annotations = MetagenomicsAnalyzer.PredictFunctions(
                genes.Select(g => (g.GeneId, g.Seq)), FunctionDb).ToList();
            bool ok = annotations.All(a =>
                a.BitScore >= 0 && a.EValue >= 0 &&
                FunctionDb.Any(kv => kv.Value.Function == a.Function && byId[a.GeneId].Contains(kv.Key, StringComparison.Ordinal)));
            return ok.Label("an annotation had a negative score or a function not backed by a contained signature");
        });
    }

    /// <summary>
    /// INV-2 (P, positive control): a gene containing the GGG signature is annotated Glycine-rich.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Functions_EmbeddedSignature_IsAnnotated()
    {
        var genes = new[] { ("gene1", "AKAGGGAK"), ("gene2", "ACDEFHIK") };
        var annotations = MetagenomicsAnalyzer.PredictFunctions(genes, FunctionDb).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(annotations.Any(a => a.GeneId == "gene1" && a.Function == "Glycine-rich"), Is.True);
            Assert.That(annotations.Any(a => a.GeneId == "gene2"), Is.False, "gene with no signature is unannotated");
        });
    }

    /// <summary>
    /// INV-3 (D): Functional prediction is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Functions_AreDeterministic()
    {
        return Prop.ForAll(FunctionGenesArbitrary(), genes =>
        {
            var a = MetagenomicsAnalyzer.PredictFunctions(genes.Select(g => (g.GeneId, g.Seq)), FunctionDb).Select(x => (x.GeneId, x.Function)).ToList();
            var b = MetagenomicsAnalyzer.PredictFunctions(genes.Select(g => (g.GeneId, g.Seq)), FunctionDb).Select(x => (x.GeneId, x.Function)).ToList();
            return a.SequenceEqual(b).Label("PredictFunctions must be deterministic");
        });
    }

    #endregion

    #region META-PATHWAY-001: R: p-value ∈ [0,1]; M: more pathway genes → higher enrichment; D: deterministic

    // FindPathwayEnrichment is a hypergeometric over-representation test: each pathway's p-value is
    // the upper-tail probability of seeing at least the observed query/pathway overlap.

    /// <summary>Generates a pathway DB (3 pathways over a 20-gene universe) and a random query subset.</summary>
    private static Arbitrary<(string[] query, Dictionary<string, IReadOnlyCollection<string>> db)> EnrichmentArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            var universe = Enumerable.Range(0, 20).Select(i => $"gene{i}").ToArray();
            var db = new Dictionary<string, IReadOnlyCollection<string>>();
            for (int p = 0; p < 3; p++)
                db[$"P{p}"] = universe.Where(_ => rng.Next(2) == 0).ToArray();
            var query = universe.Where(_ => rng.Next(3) == 0).ToArray();
            return (query, db);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R): every pathway p-value is a probability in [0,1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Pathway_PValues_InUnitInterval()
    {
        return Prop.ForAll(EnrichmentArbitrary(), input =>
        {
            var (query, db) = input;
            var results = MetagenomicsAnalyzer.FindPathwayEnrichment(query, db);
            return results.All(r => r.PValue is >= 0.0 and <= 1.0)
                .Label("a pathway p-value fell outside [0,1]");
        });
    }

    /// <summary>
    /// INV-2 (M): for a fixed pathway/background and equal query size, a larger query∩pathway overlap
    /// gives a smaller (more significant) p-value.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Pathway_MoreOverlap_LowerPValue()
    {
        var pathwayMembers = Enumerable.Range(0, 10).Select(i => $"p{i}").ToArray();
        var nonMembers = Enumerable.Range(0, 10).Select(i => $"n{i}").ToArray();
        var db = new Dictionary<string, IReadOnlyCollection<string>> { ["P"] = pathwayMembers };
        var background = pathwayMembers.Concat(nonMembers).ToArray();

        // Same query size (5), different overlap with the pathway.
        var lowOverlap = new[] { "p0", "n0", "n1", "n2", "n3" };   // x = 1
        var highOverlap = new[] { "p0", "p1", "p2", "p3", "p4" };  // x = 5

        double pLow = MetagenomicsAnalyzer.FindPathwayEnrichment(lowOverlap, db, background).Single().PValue;
        double pHigh = MetagenomicsAnalyzer.FindPathwayEnrichment(highOverlap, db, background).Single().PValue;

        Assert.That(pHigh, Is.LessThanOrEqualTo(pLow + 1e-12), "more overlap must not increase the p-value");
    }

    /// <summary>
    /// INV-3 (boundary + D): zero overlap gives p = 1; the test is deterministic.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Pathway_NoOverlap_PIsOne_AndDeterministic()
    {
        var db = new Dictionary<string, IReadOnlyCollection<string>> { ["P"] = new[] { "p0", "p1", "p2" } };
        var query = new[] { "q0", "q1" }; // disjoint from the pathway
        double p1 = MetagenomicsAnalyzer.FindPathwayEnrichment(query, db).Single().PValue;
        double p2 = MetagenomicsAnalyzer.FindPathwayEnrichment(query, db).Single().PValue;
        Assert.Multiple(() =>
        {
            Assert.That(p1, Is.EqualTo(1.0).Within(1e-12), "no overlap → p = 1");
            Assert.That(p2, Is.EqualTo(p1), "deterministic");
        });
    }

    #endregion

    #region META-RESIST-001: P: hit matches resistance DB; R: identity ∈ (0,1]; D: deterministic

    // FindResistanceGenes reports a hit when a resistance-marker motif is contained in a gene; the
    // reported name/class come from the matching DB entry and identity = motifLen / geneLen.

    private static readonly IReadOnlyDictionary<string, (string Name, string AntibioticClass)> ResistanceDb =
        new Dictionary<string, (string, string)>
        {
            ["BLA"] = ("blaTEM", "beta-lactam"),
            ["TET"] = ("tetA", "tetracycline"),
            ["VAN"] = ("vanA", "glycopeptide"),
        };

    private static Arbitrary<(string GeneId, string Seq)[]> ResistanceGenesArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            const string aa = "ACDEFGHIKLMNPQRSTVWY";
            string Rand(int n)
            {
                var c = new char[n];
                for (int i = 0; i < n; i++) c[i] = aa[rng.Next(aa.Length)];
                return new string(c);
            }
            var motifs = ResistanceDb.Keys.ToList();
            int n = 1 + rng.Next(4);
            var genes = new (string, string)[n];
            for (int i = 0; i < n; i++)
            {
                // Embed a motif in some genes; leave others to chance.
                string body = rng.Next(2) == 0 ? motifs[rng.Next(motifs.Count)] : "";
                genes[i] = ($"g{i}", Rand(rng.Next(6)) + body + Rand(rng.Next(6)));
            }
            return genes;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (P + R): every hit's name/class come from a DB entry whose motif is contained in the
    /// gene, with identity in (0,1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Resistance_HitsMatchDatabase()
    {
        return Prop.ForAll(ResistanceGenesArbitrary(), genes =>
        {
            var byId = genes.ToDictionary(g => g.GeneId, g => g.Seq);
            var hits = MetagenomicsAnalyzer.FindResistanceGenes(genes.Select(g => (g.GeneId, g.Seq)), ResistanceDb).ToList();
            bool ok = hits.All(h =>
                h.Identity is > 0.0 and <= 1.0 &&
                ResistanceDb.Any(kv => kv.Value.Name == h.ResistanceGene && kv.Value.AntibioticClass == h.AntibioticClass
                                       && byId[h.GeneId].Contains(kv.Key, StringComparison.Ordinal)));
            return ok.Label("a resistance hit did not match a contained DB motif");
        });
    }

    /// <summary>
    /// INV-2 (P, positive/negative control): a gene carrying the BLA motif is reported as blaTEM; a
    /// gene with no marker yields no hit.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Resistance_EmbeddedMotif_IsReported()
    {
        var genes = new[] { ("g1", "AKBLAQR"), ("g2", "AKDEFGH") };
        var hits = MetagenomicsAnalyzer.FindResistanceGenes(genes, ResistanceDb).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(hits.Any(h => h.GeneId == "g1" && h.ResistanceGene == "blaTEM" && h.AntibioticClass == "beta-lactam"), Is.True);
            Assert.That(hits.Any(h => h.GeneId == "g2"), Is.False, "gene without a marker yields no hit");
        });
    }

    /// <summary>
    /// INV-3 (D): Resistance detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Resistance_IsDeterministic()
    {
        return Prop.ForAll(ResistanceGenesArbitrary(), genes =>
        {
            var a = MetagenomicsAnalyzer.FindResistanceGenes(genes.Select(g => (g.GeneId, g.Seq)), ResistanceDb).ToList();
            var b = MetagenomicsAnalyzer.FindResistanceGenes(genes.Select(g => (g.GeneId, g.Seq)), ResistanceDb).ToList();
            return a.SequenceEqual(b).Label("FindResistanceGenes must be deterministic");
        });
    }

    #endregion

    #region META-TAXA-001: R: p-value ∈ [0,1]; P: significant ⟺ p < threshold; D: deterministic

    // FindSignificantTaxa runs a Mann–Whitney U test per taxon between two groups; a taxon is
    // significant when its p-value is below the threshold.

    private static Arbitrary<(IReadOnlyList<IReadOnlyDictionary<string, double>> profiles, IReadOnlyList<int> groups)>
        TaxaProfilesArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int n = 4 + rng.Next(5);
            var profiles = new List<IReadOnlyDictionary<string, double>>();
            var groups = new List<int>();
            for (int i = 0; i < n; i++)
            {
                var prof = new Dictionary<string, double>
                {
                    ["t0"] = rng.NextDouble() * 10,
                    ["t1"] = rng.NextDouble() * 10,
                    ["t2"] = rng.NextDouble() * 10,
                };
                profiles.Add(prof);
                groups.Add((i % 2) + 1); // alternate 1,2 → both groups present
            }
            return ((IReadOnlyList<IReadOnlyDictionary<string, double>>)profiles, (IReadOnlyList<int>)groups);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R + P): each taxon's p-value is in [0,1] and its significance flag equals p &lt; threshold.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Taxa_PValuesAndSignificanceFlag()
    {
        return Prop.ForAll(TaxaProfilesArbitrary(), input =>
        {
            const double threshold = 0.05;
            var (profiles, groups) = input;
            var results = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, threshold);
            bool ok = results.All(r => r.PValue is >= 0.0 and <= 1.0 && r.Significant == (r.PValue < threshold));
            return ok.Label("a taxon p-value was out of range or the significance flag disagreed with p<α");
        });
    }

    /// <summary>
    /// INV-2 (P, positive/negative control): a fully separated taxon is significant; a constant taxon
    /// is not.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Taxa_SeparatedIsSignificant_ConstantIsNot()
    {
        var profiles = new List<IReadOnlyDictionary<string, double>>();
        var groups = new List<int>();
        for (int i = 0; i < 5; i++) { profiles.Add(new Dictionary<string, double> { ["sep"] = 1, ["same"] = 5 }); groups.Add(1); }
        for (int i = 0; i < 5; i++) { profiles.Add(new Dictionary<string, double> { ["sep"] = 100, ["same"] = 5 }); groups.Add(2); }

        var results = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, 0.05);
        Assert.Multiple(() =>
        {
            Assert.That(results.Single(r => r.Taxon == "sep").Significant, Is.True, "fully separated taxon is significant");
            Assert.That(results.Single(r => r.Taxon == "same").Significant, Is.False, "constant taxon is not significant");
        });
    }

    /// <summary>
    /// INV-3 (D): Significant-taxa detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Taxa_IsDeterministic()
    {
        return Prop.ForAll(TaxaProfilesArbitrary(), input =>
        {
            var (profiles, groups) = input;
            var a = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups).Select(r => (r.Taxon, r.PValue, r.Significant)).ToList();
            var b = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups).Select(r => (r.Taxon, r.PValue, r.Significant)).ToList();
            return a.SequenceEqual(b).Label("FindSignificantTaxa must be deterministic");
        });
    }

    #endregion

    #region META-CHECKM-001: R: 0 ≤ completeness ≤ 100; R: contamination ≥ 0; D: deterministic

    // EstimateBinQualityFromMarkerCounts — domain-level CheckM (Parks et al. 2015): completeness is the
    // mean fraction of each collocated marker set present (0–100%); contamination is the mean excess of
    // duplicated markers (≥ 0%).

    private static Arbitrary<(IReadOnlyList<MetagenomicsAnalyzer.MarkerSet> sets, IReadOnlyDictionary<string, int> counts)> CheckmArbitrary() =>
        (from numSets in Gen.Choose(1, 4)
         from sets in (from k in Gen.Choose(1, 5)
                       from ids in Gen.Choose(0, 9).ArrayOf(k)
                       select new MetagenomicsAnalyzer.MarkerSet(ids.Distinct().Select(i => $"M{i}").ToList())).ArrayOf(numSets)
         from countVals in Gen.Choose(0, 3).ArrayOf(10)
         let counts = (IReadOnlyDictionary<string, int>)Enumerable.Range(0, 10)
             .Where(i => countVals[i] > 0).ToDictionary(i => $"M{i}", i => countVals[i])
         select ((IReadOnlyList<MetagenomicsAnalyzer.MarkerSet>)sets, counts)).ToArbitrary();

    /// <summary>
    /// INV-1 (R): completeness lies in [0,100]% and contamination is ≥ 0% for any marker sets / counts.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Checkm_CompletenessInRange_ContaminationNonNegative()
    {
        return Prop.ForAll(CheckmArbitrary(), t =>
        {
            var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(t.sets, t.counts);
            bool ok = q.Completeness >= -1e-9 && q.Completeness <= 100.0 + 1e-9 && q.Contamination >= -1e-9;
            return ok.Label($"completeness={q.Completeness}, contamination={q.Contamination}");
        });
    }

    /// <summary>INV-2 (D): CheckM quality estimation is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property Checkm_IsDeterministic()
    {
        return Prop.ForAll(CheckmArbitrary(), t =>
        {
            var a = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(t.sets, t.counts);
            var b = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(t.sets, t.counts);
            return (a == b).Label("EstimateBinQualityFromMarkerCounts must be deterministic");
        });
    }

    #endregion

    #region META-TETRA-001: R: correlation ∈ [-1,1]; INV: z(ACGT)=√5 on reference; D: deterministic

    // CalculateTetranucleotideZScores — Teeling/Schbath TETRA z-scores (256 tetranucleotides), strand-
    // symmetric via reverse-complement extension; TetranucleotideZScoreCorrelation is their Pearson r.

    private static Arbitrary<string> TetraDnaArbitrary() =>
        (from cs in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= 8) select new string(cs)).ToArbitrary();

    /// <summary>
    /// INV-Z (hand-derived anchor): for the reference sequence "ACGTACGTGGCC", z(ACGT) = √5 (E=3.2, var=0.128)
    /// per the Teeling/Schbath equations. Source: META-TETRA-001 TestSpec M-Z1.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Tetra_AcgtZScore_OnReference_IsSqrt5()
    {
        var z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores("ACGTACGTGGCC");
        Assert.That(z["ACGT"], Is.EqualTo(Math.Sqrt(5.0)).Within(1e-9),
            "z(ACGT) must equal √5 on the reference sequence (E=3.2, var=0.128).");
    }

    /// <summary>INV-1 (R): the tetranucleotide correlation is a Pearson r in [−1,1] for any two sequences.</summary>
    [FsCheck.NUnit.Property]
    public Property Tetra_Correlation_InUnitInterval()
    {
        var gen = (from a in TetraDnaArbitrary().Generator from b in TetraDnaArbitrary().Generator select (a, b)).ToArbitrary();
        return Prop.ForAll(gen, t =>
        {
            double r = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(t.a, t.b);
            return (r >= -1.0 - 1e-9 && r <= 1.0 + 1e-9).Label($"correlation {r} outside [-1,1]");
        });
    }

    /// <summary>
    /// INV-2 (I): a non-degenerate sequence correlates perfectly with itself (r = 1); an all-zero signature
    /// (no over/under-representation) yields the documented r = 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tetra_SelfCorrelation_IsOne()
    {
        return Prop.ForAll(TetraDnaArbitrary(), seq =>
        {
            double self = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(seq, seq);
            bool allZero = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq).Values.All(v => v == 0.0);
            return (allZero ? Math.Abs(self) < 1e-9 : Math.Abs(self - 1.0) < 1e-9)
                .Label($"self-correlation {self} (allZero={allZero})");
        });
    }

    /// <summary>INV-3 (D): the TETRA correlation is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property Tetra_IsDeterministic()
    {
        return Prop.ForAll(TetraDnaArbitrary(), seq =>
            (MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(seq, seq)
                == MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(seq, seq))
                .Label("TetranucleotideZScoreCorrelation must be deterministic"));
    }

    #endregion
}
