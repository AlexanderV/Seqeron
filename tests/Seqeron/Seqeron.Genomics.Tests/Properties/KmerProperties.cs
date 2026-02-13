using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for K-mer analysis.
/// Verifies counting invariants that must hold for ALL valid DNA sequences.
///
/// Test Units: KMER-COUNT-001, KMER-FREQ-001, KMER-FIND-001 (Property Extensions)
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

    /// <summary>
    /// INV-1: Total k-mer count == sequence_length - k + 1 for any valid sequence.
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
    /// INV-2: K-mer frequencies sum to 1.0 (within floating point tolerance).
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
    /// INV-3: Each frequency is in [0, 1].
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
    /// INV-4: All k-mer keys have length exactly k.
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
    /// INV-5: FindMostFrequentKmers returns k-mers with the maximum count.
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
    /// INV-6: K-mer entropy is non-negative.
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
    /// INV-7: Homopolymer has zero k-mer entropy for k=1.
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
}
