using System.Linq;
using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the K-mer area (counting, frequencies).
///
/// Algebraic testing pins the conservation identities of the sliding-window
/// k-mer decomposition: the counts of the distinct k-mers partition the
/// (n−k+1) windows exactly, and the normalized frequencies form a probability
/// distribution summing to one.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 32, 33.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("K-mer")]
public class KmerAlgebraicTests
{
    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: KMER-COUNT-001 — K-mer counting (K-mer)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 32.
    //
    // Model: a length-n sequence has exactly n−k+1 overlapping k-mer windows
    //        (1 ≤ k ≤ n). CountKmers groups these windows by content; the counts
    //        therefore partition the windows.
    //   — docs/algorithms/K-mer; KmerAnalyzer.CountKmers.
    //
    // Laws under test (checklist row 32):
    //   • ID   — Σ counts = n − k + 1 (the window-count identity).
    //   • DIST — counting identity / conservation: grouping the windows by k-mer
    //            loses none — Σ over distinct k-mers of count(kmer) equals the
    //            total number of windows, and each count equals the independent
    //            occurrence count of that k-mer.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: Σ counts = n − k + 1 for every 1 ≤ k ≤ n.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KmerCount_Identity_SumEqualsWindowCount()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 1), seq =>
        {
            for (int k = 1; k <= seq.Length; k++)
            {
                int sum = KmerAnalyzer.CountKmers(seq, k).Values.Sum();
                int expected = seq.Length - k + 1;
                if (sum != expected)
                    return false.Label($"k={k}: Σcounts={sum} != n-k+1={expected} for \"{seq}\"");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// DIST: counting identity — the per-k-mer count equals the number of windows
    /// that actually equal that k-mer (an independent recount), so the grouping
    /// conserves every window.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KmerCount_Distributive_CountsMatchOccurrences()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 2), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var counts = KmerAnalyzer.CountKmers(seq, k);
            var upper = seq.ToUpperInvariant();

            bool ok = counts.All(kvp =>
            {
                int occurrences = 0;
                for (int i = 0; i <= upper.Length - k; i++)
                    if (upper.Substring(i, k) == kvp.Key)
                        occurrences++;
                return occurrences == kvp.Value;
            });
            return ok.Label($"a k-mer count did not match its occurrences in \"{seq}\"");
        });
    }
}
