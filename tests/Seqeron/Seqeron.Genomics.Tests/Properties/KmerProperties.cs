using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for K-mer analysis.
/// Verifies counting invariants that must hold for ALL valid DNA sequences.
///
/// Test Units: KMER-COUNT-001, KMER-FREQ-001, KMER-FIND-001 (Property Extensions), KMER-ASYNC-001, KMER-BOTH-001
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
            var async = KmerAnalyzer.CountKmersAsync(seq, k).GetAwaiter().GetResult();
            bool same = sync.Count == async.Count
                        && sync.All(kv => async.TryGetValue(kv.Key, out int v) && v == kv.Value);
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
}
