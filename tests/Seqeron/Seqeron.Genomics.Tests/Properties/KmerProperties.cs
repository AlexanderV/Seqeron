using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for K-mer analysis.
/// Verifies counting invariants that must hold for ALL valid DNA sequences.
///
/// Test Units: KMER-COUNT-001, KMER-FREQ-001, KMER-FIND-001 (Property Extensions), KMER-ASYNC-001, KMER-BOTH-001, KMER-DIST-001, KMER-GENERATE-001, KMER-POSITIONS-001, KMER-STATS-001, KMER-UNIQUE-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class KmerProperties
{
    private static Arbitrary<string> DnaArbitrary(int minLen = 5) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    #region KMER-COUNT-001: R: count > 0; P: sum(counts) = seqLen - k + 1; M: larger k → ≤ distinct k-mers

    /// <summary>
    /// INV-1: Total k-mer count == sequence_length - k + 1 for any valid sequence.
    /// Evidence: Sliding window of size k produces exactly (n - k + 1) k-mers.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property TotalKmerCount_EqualsExpected()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var counts = KmerAnalyzer.CountKmers(seq, k);
            int totalCount = counts.Values.Sum();
            int expected = seq.Length - k + 1;
            return (totalCount == expected)
                .Label($"Total={totalCount}, expected={expected}, k={k}, len={seq.Length}");
        });
    }

    /// <summary>
    /// INV-2: Each k-mer count is strictly positive.
    /// Evidence: CountKmers only includes k-mers that appear at least once.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EachKmerCount_IsPositive()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var counts = KmerAnalyzer.CountKmers(seq, k);
            return counts.Values.All(c => c > 0)
                .Label("All k-mer counts must be > 0");
        });
    }

    /// <summary>
    /// INV-3: All k-mer keys have length exactly k.
    /// Evidence: The sliding window extracts substrings of fixed length k.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AllKmers_HaveCorrectLength()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var counts = KmerAnalyzer.CountKmers(seq, k);
            return counts.Keys.All(kmer => kmer.Length == k)
                .Label($"All k-mers must have length {k}");
        });
    }

    /// <summary>
    /// INV-4: The number of distinct k-mers is bounded by min(4^k, n - k + 1).
    /// Evidence: At most 4^k possible k-mers exist for DNA alphabet; at most (n-k+1) windows.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DistinctKmers_BoundedByTheoreticalMax()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            int k = Math.Min(3, seq.Length);
            int distinct = KmerAnalyzer.CountKmers(seq, k).Count;
            int alphabetBound = (int)Math.Pow(4, k);
            int windowBound = seq.Length - k + 1;
            int upperBound = Math.Min(alphabetBound, windowBound);
            return (distinct <= upperBound)
                .Label($"Distinct={distinct} must be ≤ min(4^{k}={alphabetBound}, n-k+1={windowBound})");
        });
    }

    #endregion

    #region KMER-FREQ-001: R: freq ∈ [0,1]; P: sum(freqs) = 1.0; D: deterministic

    /// <summary>
    /// INV-1: K-mer frequencies sum to 1.0 (within floating point tolerance).
    /// Evidence: Frequency = count / totalCount, and Σ count = totalCount.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Frequencies_SumToOne()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var freqs = KmerAnalyzer.GetKmerFrequencies(seq, k);
            double sum = freqs.Values.Sum();
            return (Math.Abs(sum - 1.0) < 0.0001)
                .Label($"Sum of frequencies={sum:F6}, expected=1.0");
        });
    }

    /// <summary>
    /// INV-2: Each frequency is in [0, 1].
    /// Evidence: frequency = count / totalCount where both are positive.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EachFrequency_InRange()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var freqs = KmerAnalyzer.GetKmerFrequencies(seq, k);
            return freqs.Values.All(f => f >= 0.0 && f <= 1.0)
                .Label("All frequencies must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-3: Frequencies are deterministic — same input always yields same result.
    /// Evidence: GetKmerFrequencies is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Frequencies_AreDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var freqs1 = KmerAnalyzer.GetKmerFrequencies(seq, k);
            var freqs2 = KmerAnalyzer.GetKmerFrequencies(seq, k);
            bool equal = freqs1.Count == freqs2.Count &&
                         freqs1.All(kvp => freqs2.ContainsKey(kvp.Key) &&
                                           Math.Abs(freqs2[kvp.Key] - kvp.Value) < 0.0001);
            return equal.Label("GetKmerFrequencies must be deterministic");
        });
    }

    #endregion

    #region KMER-FIND-001: R: positions valid; M: lower minFreq → ≥ k-mers returned; D: deterministic

    /// <summary>
    /// INV-1: FindMostFrequentKmers returns k-mers with the maximum count.
    /// Evidence: "Most frequent" means count equals the global maximum count.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MostFrequent_HasMaxCount()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var counts = KmerAnalyzer.CountKmers(seq, k);
            var mostFrequent = KmerAnalyzer.FindMostFrequentKmers(seq, k).ToList();
            if (mostFrequent.Count == 0 || counts.Count == 0) return true.ToProperty();

            int maxCount = counts.Values.Max();
            return mostFrequent.All(kmer => counts.ContainsKey(kmer) && counts[kmer] == maxCount)
                .Label("Most frequent k-mers must all have maximum count");
        });
    }

    /// <summary>
    /// INV-2: FindKmersWithMinCount returns only k-mers with count ≥ minCount.
    /// Evidence: The method filters by count threshold.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindWithMinCount_AllResultsHaveMinCount()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            int k = Math.Min(3, seq.Length);
            int minCount = 2;
            var results = KmerAnalyzer.FindKmersWithMinCount(seq, k, minCount).ToList();
            var counts = KmerAnalyzer.CountKmers(seq, k);

            return results.All(r => r.Count >= minCount && counts.ContainsKey(r.Kmer))
                .Label("All results must have count ≥ minCount");
        });
    }

    /// <summary>
    /// INV-3: Lower minCount yields more or equal results (monotonicity).
    /// Evidence: Lowering the threshold expands the result set.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindWithMinCount_LowerThreshold_MoreOrEqualResults()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var resultsHigh = KmerAnalyzer.FindKmersWithMinCount(seq, k, 3).ToList();
            var resultsLow = KmerAnalyzer.FindKmersWithMinCount(seq, k, 2).ToList();

            return (resultsLow.Count >= resultsHigh.Count)
                .Label($"minCount=2 → {resultsLow.Count}, minCount=3 → {resultsHigh.Count}");
        });
    }

    /// <summary>
    /// INV-4: FindKmersWithMinCount is deterministic.
    /// Evidence: Pure function with no side effects.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindWithMinCount_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var results1 = KmerAnalyzer.FindKmersWithMinCount(seq, k, 2).ToList();
            var results2 = KmerAnalyzer.FindKmersWithMinCount(seq, k, 2).ToList();

            return results1.SequenceEqual(results2)
                .Label("FindKmersWithMinCount must be deterministic");
        });
    }

    /// <summary>
    /// INV-5: K-mer entropy is non-negative.
    /// Evidence: Shannon entropy ≥ 0 by definition (H = -Σ p·log₂(p)).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Entropy_IsNonNegative()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            double entropy = KmerAnalyzer.CalculateKmerEntropy(seq, k);
            return (entropy >= -0.0001)
                .Label($"Entropy={entropy:F4}, must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-6: Homopolymer has zero k-mer entropy for k=1.
    /// Evidence: Single symbol → p=1 → H = -1·log₂(1) = 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Homopolymer_HasZeroEntropy()
    {
        var baseGen = Gen.Elements('A', 'C', 'G', 'T').ToArbitrary();
        return Prop.ForAll(baseGen, b =>
        {
            string homo = new(b, 20);
            double entropy = KmerAnalyzer.CalculateKmerEntropy(homo, 1);
            return (Math.Abs(entropy) < 0.0001)
                .Label($"Homopolymer '{b}' entropy={entropy:F4}, expected=0");
        });
    }

    #endregion

    #region KMER-ASYNC-001: P: async result = sync result; D: deterministic

    // CountKmersAsync wraps the synchronous counter in Task.Run, so it must return exactly the same
    // k-mer multiplicities as CountKmers for any input — concurrency must not change the result.

    /// <summary>
    /// INV-1 (P): The async counter returns the same k-mer → count map as the synchronous counter.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CountKmersAsync_EqualsSync()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var sync = KmerAnalyzer.CountKmers(seq, k);
            var asyncResult = KmerAnalyzer.CountKmersAsync(seq, k).GetAwaiter().GetResult();
            bool same = sync.Count == asyncResult.Count
                        && sync.All(kv => asyncResult.TryGetValue(kv.Key, out int v) && v == kv.Value);
            return same.Label($"async result differs from sync (k={k}, len={seq.Length})");
        });
    }

    /// <summary>
    /// INV-2 (D): Repeated async invocations on the same input yield identical results.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CountKmersAsync_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var a = KmerAnalyzer.CountKmersAsync(seq, k).GetAwaiter().GetResult();
            var b = KmerAnalyzer.CountKmersAsync(seq, k).GetAwaiter().GetResult();
            bool same = a.Count == b.Count && a.All(kv => b.TryGetValue(kv.Key, out int v) && v == kv.Value);
            return same.Label("CountKmersAsync must be deterministic");
        });
    }

    #endregion

    #region KMER-BOTH-001: P: count = forward + reverse-complement; S: strand-symmetric; D: deterministic

    // CountKmersBothStrands sums each k-mer's forward count and its count on the reverse-complement
    // strand (kPAL; Chargaff/inversion symmetry). Total = 2·(L−k+1); count(w) = count(RC(w)).

    /// <summary>
    /// INV-1 (P): The both-strands count of each k-mer equals its forward count plus the forward
    /// count of the reverse-complement sequence, and the grand total is 2·(L−k+1).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CountBothStrands_EqualsForwardPlusReverseComplement()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var both = KmerAnalyzer.CountKmersBothStrands(seq, k);
            var fwd = KmerAnalyzer.CountKmers(seq, k);
            var rc = KmerAnalyzer.CountKmers(DnaSequence.GetReverseComplementString(seq), k);

            bool decomposes = both.All(kv =>
                kv.Value == fwd.GetValueOrDefault(kv.Key) + rc.GetValueOrDefault(kv.Key));
            bool totalOk = both.Values.Sum() == 2 * (seq.Length - k + 1);
            return (decomposes && totalOk)
                .Label($"decompose={decomposes}, total={both.Values.Sum()} vs {2 * (seq.Length - k + 1)}");
        });
    }

    /// <summary>
    /// INV-2 (S): Strand symmetry — a k-mer and its reverse complement carry equal both-strands counts.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CountBothStrands_IsStrandSymmetric()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var both = KmerAnalyzer.CountKmersBothStrands(seq, k);
            bool symmetric = both.All(kv =>
            {
                string rc = DnaSequence.GetReverseComplementString(kv.Key);
                return both.TryGetValue(rc, out int v) && v == kv.Value;
            });
            return symmetric.Label("count(w) ≠ count(RC(w)) — strand symmetry violated");
        });
    }

    /// <summary>
    /// INV-3 (D): Both-strands counting is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CountBothStrands_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var a = KmerAnalyzer.CountKmersBothStrands(seq, k);
            var b = KmerAnalyzer.CountKmersBothStrands(seq, k);
            return (a.Count == b.Count && a.All(kv => b.TryGetValue(kv.Key, out int v) && v == kv.Value))
                .Label("CountKmersBothStrands must be deterministic");
        });
    }

    #endregion

    #region KMER-DIST-001: R: distance ≥ 0; S: d(a,b)=d(b,a); I: d(x,x)=0; D: deterministic

    // KmerDistance is the Euclidean distance between k-mer relative-frequency vectors (Vinga &
    // Almeida 2003; Zielezinski et al. 2017) — a metric: non-negative, symmetric, zero on equal
    // inputs, and satisfying the triangle inequality.

    private const int DistK = 3;

    /// <summary>INV-1 (R): k-mer distance is non-negative.</summary>
    [FsCheck.NUnit.Property]
    public Property KmerDistance_IsNonNegative()
    {
        return Prop.ForAll(DnaArbitrary(5), DnaArbitrary(5), (a, b) =>
            (KmerAnalyzer.KmerDistance(a, b, DistK) >= 0.0).Label("distance must be ≥ 0"));
    }

    /// <summary>INV-2 (S): distance is symmetric.</summary>
    [FsCheck.NUnit.Property]
    public Property KmerDistance_IsSymmetric()
    {
        return Prop.ForAll(DnaArbitrary(5), DnaArbitrary(5), (a, b) =>
        {
            double ab = KmerAnalyzer.KmerDistance(a, b, DistK);
            double ba = KmerAnalyzer.KmerDistance(b, a, DistK);
            return (Math.Abs(ab - ba) < 1e-12).Label($"d(a,b)={ab} ≠ d(b,a)={ba}");
        });
    }

    /// <summary>INV-3 (I): self-distance is zero.</summary>
    [FsCheck.NUnit.Property]
    public Property KmerDistance_SelfIsZero()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
            (KmerAnalyzer.KmerDistance(seq, seq, DistK) < 1e-12).Label("d(x,x) must be 0"));
    }

    /// <summary>INV-4 (metric): the triangle inequality holds, d(a,c) ≤ d(a,b)+d(b,c).</summary>
    [FsCheck.NUnit.Property]
    public Property KmerDistance_TriangleInequality()
    {
        return Prop.ForAll(DnaArbitrary(5), DnaArbitrary(5), DnaArbitrary(5), (a, b, c) =>
        {
            double ac = KmerAnalyzer.KmerDistance(a, c, DistK);
            double ab = KmerAnalyzer.KmerDistance(a, b, DistK);
            double bc = KmerAnalyzer.KmerDistance(b, c, DistK);
            return (ac <= ab + bc + 1e-9).Label($"triangle violated: {ac} > {ab}+{bc}");
        });
    }

    /// <summary>INV-5 (D): distance is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property KmerDistance_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(5), DnaArbitrary(5), (a, b) =>
            (KmerAnalyzer.KmerDistance(a, b, DistK) == KmerAnalyzer.KmerDistance(a, b, DistK))
                .Label("KmerDistance must be deterministic"));
    }

    #endregion

    #region KMER-GENERATE-001: R: count = |alphabet|^k; P: all distinct, length k, over alphabet; D: deterministic

    // GenerateAllKmers enumerates the k-fold Cartesian product of the alphabet — |alphabet|^k distinct
    // strings (4^k for DNA), in lexicographic order when the alphabet is sorted.

    private static Arbitrary<int> SmallKArbitrary() => Gen.Choose(1, 6).ToArbitrary();

    /// <summary>
    /// INV-1 (R): The number of generated k-mers equals |alphabet|^k for the DNA alphabet (4^k) and
    /// for a 5-letter alphabet (5^k).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateAllKmers_Count_IsAlphabetPowK()
    {
        return Prop.ForAll(SmallKArbitrary(), k =>
        {
            int dna = KmerAnalyzer.GenerateAllKmers(k).Count();
            int five = KmerAnalyzer.GenerateAllKmers(k, "ACGTN").Count();
            return (dna == (int)Math.Pow(4, k) && five == (int)Math.Pow(5, k))
                .Label($"k={k}: DNA {dna} vs {(int)Math.Pow(4, k)}, 5-letter {five} vs {(int)Math.Pow(5, k)}");
        });
    }

    /// <summary>
    /// INV-2 (P): Every generated k-mer is distinct, has length k, and uses only alphabet letters.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateAllKmers_AreDistinctWellFormed()
    {
        return Prop.ForAll(SmallKArbitrary(), k =>
        {
            var all = KmerAnalyzer.GenerateAllKmers(k).ToList();
            bool distinct = all.Distinct().Count() == all.Count;
            bool wellFormed = all.All(m => m.Length == k && m.All(c => "ACGT".Contains(c)));
            return (distinct && wellFormed).Label($"k={k}: distinct={distinct}, wellFormed={wellFormed}");
        });
    }

    /// <summary>
    /// INV-3 (D): Generation is deterministic and, for the sorted DNA alphabet, lexicographically ordered.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenerateAllKmers_IsDeterministicAndSorted()
    {
        return Prop.ForAll(SmallKArbitrary(), k =>
        {
            var a = KmerAnalyzer.GenerateAllKmers(k).ToList();
            var b = KmerAnalyzer.GenerateAllKmers(k).ToList();
            bool deterministic = a.SequenceEqual(b);
            bool sorted = a.SequenceEqual(a.OrderBy(x => x, StringComparer.Ordinal));
            return (deterministic && sorted).Label($"k={k}: deterministic={deterministic}, sorted={sorted}");
        });
    }

    /// <summary>
    /// INV-4 (boundary): non-positive k and an empty alphabet are rejected.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GenerateAllKmers_Boundaries()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.GenerateAllKmers(0).ToList());
            Assert.Throws<ArgumentException>(() => KmerAnalyzer.GenerateAllKmers(2, "").ToList());
        });
    }

    #endregion

    #region KMER-POSITIONS-001: R: positions ∈ [0, len−k]; P: seq[pos..pos+k] = kmer; D: deterministic

    // FindKmerPositions reports every overlapping start position of a k-mer in the sequence
    // (Rosalind BA1D). Each position p satisfies 0 ≤ p ≤ len−k and seq[p..p+k] == kmer.

    /// <summary>Generates a sequence together with one of its own substrings as the query k-mer.</summary>
    private static Arbitrary<(string seq, string kmer)> SeqAndSubstringArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            const string bases = "ACGT";
            int len = 15 + rng.Next(15);
            var chars = new char[len];
            for (int i = 0; i < len; i++) chars[i] = bases[rng.Next(4)];
            string seq = new string(chars);
            int k = 2 + rng.Next(4);
            int start = rng.Next(len - k + 1);
            return (seq, seq.Substring(start, k));
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R + P): every reported position is in [0, len−k], the substring there equals the query
    /// k-mer, positions are ascending, and at least one occurrence is found (the source position).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KmerPositions_AreValidOccurrences()
    {
        return Prop.ForAll(SeqAndSubstringArbitrary(), input =>
        {
            var (seq, kmer) = input;
            var pos = KmerAnalyzer.FindKmerPositions(seq, kmer).ToList();
            bool valid = pos.All(p => p >= 0 && p <= seq.Length - kmer.Length
                                      && seq.Substring(p, kmer.Length) == kmer);
            bool ascending = pos.SequenceEqual(pos.OrderBy(x => x));
            return (pos.Count >= 1 && valid && ascending)
                .Label($"positions invalid/empty/unsorted for kmer '{kmer}'");
        });
    }

    /// <summary>
    /// INV-2 (completeness): the reported positions are exactly all overlapping occurrences, matched
    /// against an independent scan.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KmerPositions_AreComplete()
    {
        return Prop.ForAll(SeqAndSubstringArbitrary(), input =>
        {
            var (seq, kmer) = input;
            var pos = KmerAnalyzer.FindKmerPositions(seq, kmer).ToList();
            var expected = Enumerable.Range(0, seq.Length - kmer.Length + 1)
                .Where(i => seq.Substring(i, kmer.Length) == kmer).ToList();
            return pos.SequenceEqual(expected).Label($"positions ≠ independent scan for '{kmer}'");
        });
    }

    /// <summary>
    /// INV-3 (D): Position finding is deterministic; empty/oversized queries yield no positions.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KmerPositions_IsDeterministic()
    {
        return Prop.ForAll(SeqAndSubstringArbitrary(), input =>
        {
            var (seq, kmer) = input;
            return KmerAnalyzer.FindKmerPositions(seq, kmer)
                .SequenceEqual(KmerAnalyzer.FindKmerPositions(seq, kmer))
                .Label("FindKmerPositions must be deterministic");
        });
    }

    /// <summary>
    /// INV-4 (boundary): empty sequence/k-mer and over-long k-mers return no positions.
    /// </summary>
    [Test]
    [Category("Property")]
    public void KmerPositions_Boundaries()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.FindKmerPositions("", "AC"), Is.Empty);
            Assert.That(KmerAnalyzer.FindKmerPositions("ACGT", ""), Is.Empty);
            Assert.That(KmerAnalyzer.FindKmerPositions("AC", "ACGT"), Is.Empty, "k-mer longer than sequence");
        });
    }

    #endregion

    #region KMER-STATS-001: R: counts ≥ 0; P: total k-mers = len−k+1; D: deterministic

    // AnalyzeKmers summarises k-mer composition: TotalKmers is the number of overlapping windows
    // L−k+1, UniqueKmers the distinct count, with consistent max/min/average multiplicity and entropy.

    /// <summary>
    /// INV-1 (P): The reported total k-mer count equals len − k + 1.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeKmers_TotalEqualsWindowCount()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var stats = KmerAnalyzer.AnalyzeKmers(seq, k);
            return (stats.TotalKmers == seq.Length - k + 1)
                .Label($"TotalKmers={stats.TotalKmers}, expected {seq.Length - k + 1}");
        });
    }

    /// <summary>
    /// INV-2 (R): Statistics are internally consistent — counts non-negative, every distinct k-mer
    /// observed at least once, max ≥ min, distinct ≤ total, entropy ≥ 0, average = total/distinct.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeKmers_StatsAreConsistent()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var s = KmerAnalyzer.AnalyzeKmers(seq, k);
            bool ok = s.MinCount >= 1 && s.MaxCount >= s.MinCount
                      && s.UniqueKmers >= 1 && s.UniqueKmers <= s.TotalKmers
                      && s.Entropy >= -1e-9
                      && Math.Abs(s.AverageCount - Math.Round((double)s.TotalKmers / s.UniqueKmers, 2)) < 1e-9;
            return ok.Label($"inconsistent stats: {s}");
        });
    }

    /// <summary>
    /// INV-3 (D): Statistics are deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeKmers_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(5), seq =>
        {
            int k = Math.Min(3, seq.Length);
            return (KmerAnalyzer.AnalyzeKmers(seq, k) == KmerAnalyzer.AnalyzeKmers(seq, k))
                .Label("AnalyzeKmers must be deterministic");
        });
    }

    /// <summary>
    /// INV-4 (boundary): when k exceeds the length there are no k-mers and all stats are zero.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AnalyzeKmers_NoKmers_AllZero()
    {
        var s = KmerAnalyzer.AnalyzeKmers("ACG", 5);
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalKmers, Is.Zero);
            Assert.That(s.UniqueKmers, Is.Zero);
            Assert.That(s.Entropy, Is.Zero);
        });
    }

    #endregion

    #region KMER-UNIQUE-001: P: unique k-mers have count 1; P: min-count filter is exact and monotone; D: deterministic

    // FindUniqueKmers returns the k-mers occurring exactly once; FindKmersWithMinCount returns those
    // occurring ≥ minCount, ordered by count descending (Compeau & Pevzner recurrent k-mers).

    private const int UniqueK = 2; // small k so repeats (count > 1) actually arise

    /// <summary>
    /// INV-1 (P): FindUniqueKmers is exactly the set of k-mers with occurrence count 1.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property UniqueKmers_AreExactlyCountOne()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            var counts = KmerAnalyzer.CountKmers(seq, UniqueK);
            var unique = KmerAnalyzer.FindUniqueKmers(seq, UniqueK).ToHashSet();
            var expected = counts.Where(c => c.Value == 1).Select(c => c.Key).ToHashSet();
            return unique.SetEquals(expected).Label("FindUniqueKmers ≠ {k-mers with count 1}");
        });
    }

    /// <summary>
    /// INV-2 (P): FindKmersWithMinCount returns exactly the k-mers with count ≥ minCount, ordered by
    /// count descending.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MinCountKmers_AreExactAndOrdered()
    {
        var gen = Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            const string bases = "ACGT";
            int len = 10 + rng.Next(15);
            var chars = new char[len];
            for (int i = 0; i < len; i++) chars[i] = bases[rng.Next(4)];
            return (new string(chars), 1 + rng.Next(4)); // minCount 1..4
        }).ToArbitrary();

        return Prop.ForAll(gen, input =>
        {
            var (seq, minCount) = input;
            var res = KmerAnalyzer.FindKmersWithMinCount(seq, UniqueK, minCount).ToList();
            var counts = KmerAnalyzer.CountKmers(seq, UniqueK);
            var expected = counts.Where(c => c.Value >= minCount).Select(c => c.Key).ToHashSet();

            bool exact = res.Select(r => r.Kmer).ToHashSet().SetEquals(expected);
            bool meetsMin = res.All(r => r.Count >= minCount);
            bool ordered = res.Select(r => r.Count).SequenceEqual(res.Select(r => r.Count).OrderByDescending(x => x));
            return (exact && meetsMin && ordered)
                .Label($"minCount={minCount}: exact={exact}, meetsMin={meetsMin}, ordered={ordered}");
        });
    }

    /// <summary>
    /// INV-3 (P, monotone): a lower minCount returns a superset — min-count 1 contains all min-count 2 k-mers.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MinCountKmers_LowerThreshold_IsSuperset()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            var loose = KmerAnalyzer.FindKmersWithMinCount(seq, UniqueK, 1).Select(r => r.Kmer).ToHashSet();
            var strict = KmerAnalyzer.FindKmersWithMinCount(seq, UniqueK, 2).Select(r => r.Kmer).ToHashSet();
            return strict.IsSubsetOf(loose).Label("min-count 2 result not ⊆ min-count 1 result");
        });
    }

    /// <summary>
    /// INV-4 (D): Both finders are deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property UniqueKmers_AreDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            var a = KmerAnalyzer.FindUniqueKmers(seq, UniqueK).ToHashSet();
            var b = KmerAnalyzer.FindUniqueKmers(seq, UniqueK).ToHashSet();
            return a.SetEquals(b).Label("FindUniqueKmers must be deterministic");
        });
    }

    #endregion
}
