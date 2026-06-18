using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for metagenomics diversity analysis:
/// alpha diversity indices and beta diversity metrics.
///
/// Test Units: META-ALPHA-001, META-BETA-001
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
            if (kmer.All(c => "ACGT".IndexOf(c) >= 0))
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
}
