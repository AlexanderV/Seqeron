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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: KMER-FREQ-001 — K-mer frequencies (K-mer)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 33.
    //
    // Model: the relative frequency of a k-mer is its count divided by the total
    //        number of k-mer windows; the frequencies thus form a probability
    //        distribution over the observed k-mers.
    //   — docs/algorithms/K-mer; KmerAnalyzer.GetKmerFrequencies / CountKmers.
    //
    // Laws under test (checklist row 33):
    //   • ID   — Σ frequencies = 1.0 (normalization / probability-simplex identity).
    //   • DIST — freq(kmer) = count(kmer) / total_kmers for every k-mer.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: the k-mer frequencies sum to 1 for any non-empty k-mer set
    /// (1 ≤ k ≤ n).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KmerFreq_Identity_FrequenciesSumToOne()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 1), seq =>
        {
            for (int k = 1; k <= seq.Length; k++)
            {
                var freqs = KmerAnalyzer.GetKmerFrequencies(seq, k);
                double sum = freqs.Values.Sum();
                if (Math.Abs(sum - 1.0) > 1e-9)
                    return false.Label($"k={k}: Σfreq={sum} != 1 for \"{seq}\"");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// DIST: freq(kmer) = count(kmer) / total — each frequency is exactly the
    /// k-mer's share of the n−k+1 windows.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KmerFreq_Distributive_FreqIsCountOverTotal()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 2), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var counts = KmerAnalyzer.CountKmers(seq, k);
            var freqs = KmerAnalyzer.GetKmerFrequencies(seq, k);
            double total = counts.Values.Sum();

            bool ok = counts.All(kvp =>
                Math.Abs(freqs[kvp.Key] - kvp.Value / total) < 1e-12);
            return ok.Label($"a frequency != count/total for \"{seq}\"");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: KMER-BOTH-001 — Both-strand k-mer counting (K-mer)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 157.
    //
    // Model: counting k-mers over both strands sums the forward-strand count and
    //        the reverse-complement-strand count. Because a k-mer on the reverse
    //        strand reads as its reverse complement on the forward strand, the
    //        both-strand histogram is strand-symmetric: count(w) = count(RC(w)).
    //   — docs/algorithms/K-mer; KmerAnalyzer.CountKmersBothStrands.
    //
    // Laws under test (checklist row 157):
    //   • ADD  — both-strand count(w) = forward(w) + revcomp(w); the grand total
    //            equals 2·(n−k+1).
    //   • COMM — strand symmetry: count(w) = count(RC(w)).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>ADD: each both-strand count is the sum of the forward and the
    /// reverse-complement-strand counts; the total is 2·(n−k+1).</summary>
    [FsCheck.NUnit.Property]
    public Property KmerBoth_Additive_ForwardPlusReverse()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 4), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var both = KmerAnalyzer.CountKmersBothStrands(seq, k);
            var fwd = KmerAnalyzer.CountKmers(seq, k);
            var rev = KmerAnalyzer.CountKmers(DnaSequence.GetReverseComplementString(seq), k);

            bool perKmer = both.All(kvp =>
                kvp.Value == fwd.GetValueOrDefault(kvp.Key) + rev.GetValueOrDefault(kvp.Key));
            int expectedTotal = 2 * (seq.Length - k + 1);
            return (perKmer && both.Values.Sum() == expectedTotal)
                .Label($"both-strand counts not additive for \"{seq}\"");
        });
    }

    /// <summary>COMM: the both-strand histogram is strand-symmetric —
    /// count(w) = count(reverseComplement(w)).</summary>
    [FsCheck.NUnit.Property]
    public Property KmerBoth_Commutative_StrandSymmetric()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 4), seq =>
        {
            int k = Math.Min(3, seq.Length);
            var both = KmerAnalyzer.CountKmersBothStrands(seq, k);
            bool ok = both.All(kvp =>
                kvp.Value == both.GetValueOrDefault(DnaSequence.GetReverseComplementString(kvp.Key)));
            return ok.Label($"both-strand histogram not strand-symmetric for \"{seq}\"");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: KMER-DIST-001 — Alignment-free k-mer distance (K-mer)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 158.
    //
    // Model: the k-mer distance is the Euclidean distance between the two
    //        sequences' k-mer frequency vectors — a genuine metric.
    //   — docs/algorithms/K-mer; KmerAnalyzer.KmerDistance.
    //
    // Laws under test (checklist row 158):
    //   • ID       — d(x, x) = 0.
    //   • COMM     — d(a, b) = d(b, a).
    //   • TRIANGLE — d(a, c) ≤ d(a, b) + d(b, c).
    // ═══════════════════════════════════════════════════════════════════════

    private static Arbitrary<(string A, string B, string C)> DnaTriple() =>
        (from a in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(x => x.Length >= 3)
         from b in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(x => x.Length >= 3)
         from c in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(x => x.Length >= 3)
         select (new string(a), new string(b), new string(c)))
        .ToArbitrary();

    /// <summary>ID, COMM and TRIANGLE for the k-mer Euclidean distance.</summary>
    [FsCheck.NUnit.Property]
    public Property KmerDistance_MetricAxioms()
    {
        return Prop.ForAll(DnaTriple(), t =>
        {
            const int k = 3;
            double aa = KmerAnalyzer.KmerDistance(t.A, t.A, k);
            double ab = KmerAnalyzer.KmerDistance(t.A, t.B, k);
            double ba = KmerAnalyzer.KmerDistance(t.B, t.A, k);
            double ac = KmerAnalyzer.KmerDistance(t.A, t.C, k);
            double bc = KmerAnalyzer.KmerDistance(t.B, t.C, k);
            bool id = Math.Abs(aa) < 1e-12;
            bool comm = Math.Abs(ab - ba) < 1e-12;
            bool tri = ac <= ab + bc + 1e-9;
            return (id && comm && tri)
                .Label($"id={id}, comm={comm}, tri={tri} (ac={ac}, ab={ab}, bc={bc})");
        });
    }
}
