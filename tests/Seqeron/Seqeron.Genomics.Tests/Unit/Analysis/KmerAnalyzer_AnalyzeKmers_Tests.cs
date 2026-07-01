// KMER-STATS-001 — K-mer Statistics
// Evidence: docs/Evidence/KMER-STATS-001-Evidence.md
// TestSpec: tests/TestSpecs/KMER-STATS-001.md
// Source: Wikipedia — K-mer (https://en.wikipedia.org/wiki/K-mer);
//         Clavijo B (2018), BioInfoLogics — k-mer counting part I
//         (https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/);
//         Manca V et al. (2021), Spectral concepts in genome informational analysis,
//         arXiv:2106.15351 (k-entropy E_k = -Σ p log2 p, p = mult/(L-k+1)).

using NUnit.Framework;
using System;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Tests for KMER-STATS-001: comprehensive k-mer statistics (KmerAnalyzer.AnalyzeKmers).
///
/// Expected values are derived from the source definitions, not from the implementation:
/// TotalKmers = L-k+1, UniqueKmers = distinct count, Max/Min/Average over the multiplicity
/// table, and Shannon entropy E_k = -Σ p log2 p with p = count/(L-k+1). The GTAGAGCTGT and
/// AGAT worked examples come from the Wikipedia K-mer example tables; ATCGATCAC from
/// BioInfoLogics. Entropy constants are computed independently from the frequency formula.
/// </summary>
[TestFixture]
public class KmerAnalyzer_AnalyzeKmers_Tests
{
    // Wikipedia K-mer example sequence (L=10).
    private const string Gtagagctgt = "GTAGAGCTGT";
    // Wikipedia AGAT example.
    private const string Agat = "AGAT";
    // BioInfoLogics example (L=9).
    private const string Atcgatcac = "ATCGATCAC";

    #region AnalyzeKmers — MUST (worked examples)

    // M1 — GTAGAGCTGT k=1: G4 T3 A2 C1 => total 10, distinct 4, max 4, min 1, avg 2.5.
    // Entropy = -(0.4 log2 0.4 + 0.3 log2 0.3 + 0.2 log2 0.2 + 0.1 log2 0.1) = 1.846439344671 bits.
    [Test]
    public void AnalyzeKmers_GtagagctgtK1_MatchesWikipediaMonomerTable()
    {
        var s = KmerAnalyzer.AnalyzeKmers(Gtagagctgt, 1);

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalKmers, Is.EqualTo(10), "L-k+1 = 10-1+1 = 10 monomers.");
            Assert.That(s.UniqueKmers, Is.EqualTo(4), "Distinct monomers G,T,A,C = 4 (Wikipedia table).");
            Assert.That(s.MaxCount, Is.EqualTo(4), "G occurs 4 times — the maximum.");
            Assert.That(s.MinCount, Is.EqualTo(1), "C occurs once — the minimum.");
            Assert.That(s.AverageCount, Is.EqualTo(2.5).Within(1e-10), "Average = 10/4 = 2.5.");
            Assert.That(s.Entropy, Is.EqualTo(1.846439344671).Within(1e-10),
                "Shannon entropy of {0.4,0.3,0.2,0.1} = 1.84643934... bits.");
        });
    }

    // M2 — GTAGAGCTGT k=2: GT,AG each x2; total 9, distinct 7, max 2, min 1, avg 9/7=1.29(rounded).
    // Entropy over {2,2,1,1,1,1,1}/9 = 2.725480556998 bits.
    [Test]
    public void AnalyzeKmers_GtagagctgtK2_MatchesWikipediaDimerTable()
    {
        var s = KmerAnalyzer.AnalyzeKmers(Gtagagctgt, 2);

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalKmers, Is.EqualTo(9), "L-k+1 = 10-2+1 = 9 dimers.");
            Assert.That(s.UniqueKmers, Is.EqualTo(7), "7 distinct dimers (GT and AG each appear twice).");
            Assert.That(s.MaxCount, Is.EqualTo(2), "GT (and AG) occur twice — the maximum.");
            Assert.That(s.MinCount, Is.EqualTo(1), "The other five dimers occur once.");
            Assert.That(s.AverageCount, Is.EqualTo(1.29).Within(1e-10), "Average = 9/7 ≈ 1.29 (rounded to 2 dp).");
            Assert.That(s.Entropy, Is.EqualTo(2.725480556998).Within(1e-10),
                "Shannon entropy over {2,2,1,1,1,1,1}/9 = 2.72548055... bits.");
        });
    }

    // M3 — GTAGAGCTGT k=3: all 8 windows distinct => entropy = log2(8) = 3 exactly.
    [Test]
    public void AnalyzeKmers_GtagagctgtK3_AllDistinct_EntropyIsLog2Eight()
    {
        var s = KmerAnalyzer.AnalyzeKmers(Gtagagctgt, 3);

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalKmers, Is.EqualTo(8), "L-k+1 = 10-3+1 = 8 trimers.");
            Assert.That(s.UniqueKmers, Is.EqualTo(8), "All 8 trimers are distinct (Wikipedia table).");
            Assert.That(s.MaxCount, Is.EqualTo(1), "Every trimer occurs once.");
            Assert.That(s.MinCount, Is.EqualTo(1), "Every trimer occurs once.");
            Assert.That(s.AverageCount, Is.EqualTo(1.0).Within(1e-10), "Average = 8/8 = 1.0.");
            Assert.That(s.Entropy, Is.EqualTo(3.0).Within(1e-10),
                "8 equiprobable k-mers => H = log2(8) = 3 bits exactly.");
        });
    }

    // M4 — ATCGATCAC k=3: ATC=2, rest=1 => total 7, distinct 6, max 2, min 1, avg 7/6=1.17.
    // Entropy = -(2/7 log2 2/7 + 5*(1/7 log2 1/7)) = 2.521640636343 bits.
    [Test]
    public void AnalyzeKmers_AtcgatcacK3_MatchesBioInfoLogicsTable()
    {
        var s = KmerAnalyzer.AnalyzeKmers(Atcgatcac, 3);

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalKmers, Is.EqualTo(7), "L-k+1 = 9-3+1 = 7 trimers.");
            Assert.That(s.UniqueKmers, Is.EqualTo(6), "6 distinct trimers (ATC appears twice).");
            Assert.That(s.MaxCount, Is.EqualTo(2), "ATC occurs twice — the maximum.");
            Assert.That(s.MinCount, Is.EqualTo(1), "The other five trimers occur once.");
            Assert.That(s.AverageCount, Is.EqualTo(1.17).Within(1e-10), "Average = 7/6 ≈ 1.17 (rounded).");
            Assert.That(s.Entropy, Is.EqualTo(2.521640636343).Within(1e-10),
                "Shannon entropy over {2,1,1,1,1,1}/7 = 2.52164063... bits.");
        });
    }

    // M5 — AGAT k=2: AG,GA,AT each once => total 3, distinct 3, uniform, entropy = log2(3).
    [Test]
    public void AnalyzeKmers_AgatK2_UniformThreeDimers()
    {
        var s = KmerAnalyzer.AnalyzeKmers(Agat, 2);

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalKmers, Is.EqualTo(3), "L-k+1 = 4-2+1 = 3 dimers (Wikipedia AGAT example).");
            Assert.That(s.UniqueKmers, Is.EqualTo(3), "AG, GA, AT are all distinct.");
            Assert.That(s.MaxCount, Is.EqualTo(1), "Each dimer occurs once.");
            Assert.That(s.MinCount, Is.EqualTo(1), "Each dimer occurs once.");
            Assert.That(s.AverageCount, Is.EqualTo(1.0).Within(1e-10), "Average = 3/3 = 1.0.");
            Assert.That(s.Entropy, Is.EqualTo(1.584962500721).Within(1e-10),
                "3 equiprobable k-mers => H = log2(3) = 1.58496250... bits.");
        });
    }

    #endregion

    #region AnalyzeKmers — MUST (invariant cross-checks)

    // M6 — INV-1 & INV-2: TotalKmers equals L-k+1 AND equals the sum of all k-mer counts,
    // cross-checked independently against CountKmers.
    [Test]
    public void AnalyzeKmers_TotalKmers_EqualsLMinusKPlus1AndSumOfCounts()
    {
        const string seq = "ATCGATCACGATCG"; // L=14
        const int k = 3;
        var counts = KmerAnalyzer.CountKmers(seq, k);

        var s = KmerAnalyzer.AnalyzeKmers(seq, k);

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalKmers, Is.EqualTo(seq.Length - k + 1),
                "INV-1: TotalKmers must equal L-k+1.");
            Assert.That(s.TotalKmers, Is.EqualTo(counts.Values.Sum()),
                "INV-2: TotalKmers must equal the sum of all k-mer multiplicities.");
        });
    }

    // M7 — INV-3: UniqueKmers equals the distinct k-mer count from CountKmers.
    [Test]
    public void AnalyzeKmers_UniqueKmers_EqualsDistinctCount()
    {
        const string seq = "ATCGATCACGATCG";
        const int k = 3;
        var distinct = KmerAnalyzer.CountKmers(seq, k).Count;

        var s = KmerAnalyzer.AnalyzeKmers(seq, k);

        Assert.That(s.UniqueKmers, Is.EqualTo(distinct),
            "INV-3: UniqueKmers must equal the number of distinct k-mers.");
    }

    #endregion

    #region AnalyzeKmers — SHOULD (corner cases / invariants)

    // S1 — INV-5 lower bound: a single distinct k-mer (homopolymer) gives entropy 0.
    [Test]
    public void AnalyzeKmers_HomopolymerK2_EntropyZeroAndExtremesEqualTotal()
    {
        var s = KmerAnalyzer.AnalyzeKmers("AAAA", 2); // AA x3

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalKmers, Is.EqualTo(3), "L-k+1 = 4-2+1 = 3.");
            Assert.That(s.UniqueKmers, Is.EqualTo(1), "Only one distinct dimer AA.");
            Assert.That(s.MaxCount, Is.EqualTo(3), "AA occurs 3 times.");
            Assert.That(s.MinCount, Is.EqualTo(3), "Single distinct k-mer => max == min == total.");
            Assert.That(s.AverageCount, Is.EqualTo(3.0).Within(1e-10), "Average = 3/1 = 3.");
            Assert.That(s.Entropy, Is.EqualTo(0.0).Within(1e-10),
                "INV-5: one-component distribution (p=1) => H = -1*log2(1) = 0.");
        });
    }

    // S2 — INV-4: MinCount <= AverageCount <= MaxCount and AverageCount = Total/Unique.
    [Test]
    public void AnalyzeKmers_MinAvgMax_RespectOrderingInvariant()
    {
        var s = KmerAnalyzer.AnalyzeKmers(Gtagagctgt, 2); // total 9, distinct 7

        Assert.Multiple(() =>
        {
            Assert.That(s.AverageCount, Is.GreaterThanOrEqualTo(s.MinCount),
                "INV-4: average must be >= min count.");
            Assert.That(s.AverageCount, Is.LessThanOrEqualTo(s.MaxCount),
                "INV-4: average must be <= max count.");
            Assert.That(s.AverageCount,
                Is.EqualTo(Math.Round((double)s.TotalKmers / s.UniqueKmers, 2)).Within(1e-10),
                "INV-4: AverageCount = round(Total/Unique, 2).");
        });
    }

    // S3 — INV-5 upper bound: entropy never exceeds log2(distinct), checked on several inputs.
    [Test]
    public void AnalyzeKmers_Entropy_WithinLog2DistinctBound()
    {
        foreach (var (seq, k) in new[] { (Gtagagctgt, 1), (Gtagagctgt, 2), (Atcgatcac, 3) })
        {
            var s = KmerAnalyzer.AnalyzeKmers(seq, k);
            double bound = Math.Log2(s.UniqueKmers);
            Assert.That(s.Entropy, Is.LessThanOrEqualTo(bound + 1e-10),
                $"INV-5: entropy of {seq} k={k} must not exceed log2(distinct)={bound}.");
        }
    }

    #endregion

    #region AnalyzeKmers — edge / error inputs

    // M8 — empty sequence => all-zero statistics (no k-mers).
    [Test]
    public void AnalyzeKmers_EmptySequence_ReturnsAllZero()
    {
        var s = KmerAnalyzer.AnalyzeKmers("", 3);

        Assert.That(s, Is.EqualTo(new KmerStatistics(0, 0, 0, 0, 0, 0)),
            "Empty sequence has no k-mers => all-zero KmerStatistics (INV-6).");
    }

    // C2 — null sequence => all-zero (CountKmers treats null as empty).
    [Test]
    public void AnalyzeKmers_NullSequence_ReturnsAllZero()
    {
        var s = KmerAnalyzer.AnalyzeKmers(null!, 3);

        Assert.That(s, Is.EqualTo(new KmerStatistics(0, 0, 0, 0, 0, 0)),
            "Null sequence is treated as empty => all-zero KmerStatistics.");
    }

    // M9 — k > L => L-k+1 < 0 => no k-mers => all-zero.
    [Test]
    public void AnalyzeKmers_KExceedsLength_ReturnsAllZero()
    {
        var s = KmerAnalyzer.AnalyzeKmers("ACG", 5);

        Assert.That(s, Is.EqualTo(new KmerStatistics(0, 0, 0, 0, 0, 0)),
            "k=5 > length 3 gives L-k+1 < 0, so no k-mers (INV-6).");
    }

    // M10 — k <= 0 is invalid (k-mer length must be positive).
    [Test]
    public void AnalyzeKmers_NonPositiveK_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => KmerAnalyzer.AnalyzeKmers("ACGT", 0),
            "k must be positive; k=0 is invalid.");
    }

    // C1 — case-insensitivity: lower-case input yields identical statistics.
    [Test]
    public void AnalyzeKmers_LowerCaseInput_MatchesUpperCase()
    {
        var lower = KmerAnalyzer.AnalyzeKmers("gtagagctgt", 1);
        var upper = KmerAnalyzer.AnalyzeKmers(Gtagagctgt, 1);

        Assert.That(lower, Is.EqualTo(upper),
            "Input is upper-cased internally; lower-case GTAGAGCTGT gives identical statistics.");
    }

    #endregion
}
