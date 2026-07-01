using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Chromosome;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Assembly area — Assembly Statistics (ASSEMBLY-STATS-001), the unit that
/// summarizes assembly contiguity/composition from contig lengths: N50/L50, Nx/Lx, total length,
/// contig count, mean/largest/smallest length and auN, via
/// <see cref="GenomeAssemblyAnalyzer.CalculateStatistics(IEnumerable{ValueTuple{string,string}})"/>,
/// <see cref="GenomeAssemblyAnalyzer.CalculateNx(IReadOnlyList{int}, long, int)"/>,
/// <see cref="GenomeAssemblyAnalyzer.CalculateNx(IEnumerable{int}, int)"/>,
/// <see cref="GenomeAssemblyAnalyzer.CalculateN50(IEnumerable{int})"/> and
/// <see cref="GenomeAssemblyAnalyzer.CalculateAuN(IEnumerable{int})"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate and boundary inputs (an EMPTY assembly, a SINGLE contig,
/// many EQUAL-LENGTH contigs, contigs whose total is odd so the 50% cutoff falls
/// between contigs, contigs sitting exactly on the ≥50% boundary, and large random
/// length sets) and asserts the unit NEVER fails in an undisciplined way:
///   • no DivideByZero / crash on an empty assembly — CalculateStatistics returns an
///     all-zero record, CalculateNx returns Nx=Lx=0, CalculateAuN returns 0 (§3.3, §6.1);
///   • no off-by-one in the canonical N50 boundary — the cumulative threshold is the
///     INCLUSIVE integer-exact test `cumulative*100 ≥ total*x` (INV-05, §4.2 [3]);
///   • no wrong L50 rank on ties — equal-length contigs select the documented rank k
///     where the cumulative FIRST reaches ≥50% (INV-02);
///   • the result is always WELL-FORMED — all stats non-negative, TotalLength = Σ lengths,
///     N50 ≤ LargestContig, SmallestContig ≤ N50 ≤ LargestContig, L50 in [1, count].
/// Every input resolves to a well-formed statistics record (the empty assembly is the
/// guarded all-zero record). A DivideByZero, an off-by-one N50, a wrong L50 on a tie, or
/// a TotalLength mismatch is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ASSEMBLY-STATS-001 — Assembly Statistics (N50 / L50 / Nx / Lx / auN)
/// Checklist: docs/checklists/03_FUZZING.md, row 147.
/// Algorithm doc: docs/algorithms/Assembly/Assembly_Statistics.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row "empty assembly, single contig, equal-length contigs":
///          – EMPTY ASSEMBLY: no contigs → all-zero record (N50=L50=0, total=0,
///            count=0, mean=0); CalculateNx → Nx=Lx=0; CalculateAuN → 0. No
///            DivideByZero on mean = total/count (§6.1, deviation #1).
///          – SINGLE CONTIG: one contig of length L → N50=L, L50=1, total=L,
///            largest=smallest=L, mean=L; the lone contig is the whole assembly
///            and reaches 100% ≥ 50% at rank 1 (§6.1).
///          – EQUAL-LENGTH CONTIGS: k contigs all of length L → N50=L and L50 is
///            the EXACT documented rank where the cumulative first reaches ≥50% of
///            the total k·L, i.e. ⌈k/2⌉ — the canonical tie-handling (INV-02, INV-05).
///   • Also swept: the inclusive ≥50% boundary on odd totals (off-by-one guard) and
///     contigs landing EXACTLY on 50% (INV-05 inclusive selection).
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Assembly_Statistics.md §2.2, §2.4, §3, §4, §6)
/// ───────────────────────────────────────────────────────────────────────────
/// Sort lengths descending L₁≥…≥Lₙ, T=ΣLᵢ. For threshold x, take the smallest prefix
/// L₁…Lₖ with Σᵢ₌₁ᵏ Lᵢ ≥ (x/100)·T; then Nx=Lₖ and Lx=k (INV-01/INV-02). The cutoff is
/// INCLUSIVE and evaluated with integer arithmetic `cumulative*100 ≥ T*x` (INV-05, §4.2).
/// N50/L50 is x=50. auN = ΣLᵢ²/T (INV-04). Empty input → all-zero record / Nx=Lx=0 /
/// auN=0, no exception (§3.3, §6.1). N50 is a length; L50 is a count (§2.2).
///   GenomeAssemblyAnalyzer.CalculateStatistics(IEnumerable&lt;(string Id, string Sequence)&gt;)
///       → AssemblyStatistics (TotalSequences, TotalLength, N50, L50, N90, L90,
///                             LargestContig, SmallestContig, MeanLength, …)
///   GenomeAssemblyAnalyzer.CalculateNx(IReadOnlyList&lt;int&gt; sortedLengths, long total, int x)
///       → NxStatistics (Threshold, Nx, Lx, CumulativeLength)
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyStatsFuzzTests
{
    private const int N50ThresholdPercent = 50;

    #region Helpers

    private static readonly char[] DnaAlphabet = { 'A', 'C', 'G', 'T' };

    private static string RandomSeq(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = DnaAlphabet[rng.Next(DnaAlphabet.Length)];
        return new string(chars);
    }

    /// <summary>Builds (Id, Sequence) contigs of the given lengths over the DNA alphabet.</summary>
    private static (string Id, string Sequence)[] ContigsOfLengths(Random rng, IReadOnlyList<int> lengths)
    {
        var contigs = new (string, string)[lengths.Count];
        for (int i = 0; i < lengths.Count; i++)
            contigs[i] = ($"c{i}", RandomSeq(rng, lengths[i]));
        return contigs;
    }

    /// <summary>
    /// Independent oracle for Nx/Lx straight from the §2.2 definition: sort descending, accumulate,
    /// and pick the FIRST contig whose cumulative reaches ≥ (x/100)·T using the same integer-exact
    /// inclusive test (INV-05, §4.2). Does not call the unit. Returns (Nx, Lx, cumulative).
    /// </summary>
    private static (int Nx, int Lx, long Cumulative) OracleNx(IReadOnlyList<int> lengths, int threshold)
    {
        var sorted = lengths.OrderByDescending(l => l).ToList();
        long total = sorted.Sum(l => (long)l);
        if (sorted.Count == 0 || total == 0)
            return (0, 0, 0);

        long cumulative = 0;
        int count = 0;
        foreach (int len in sorted)
        {
            cumulative += len;
            count++;
            if (cumulative * 100 >= total * threshold)
                return (len, count, cumulative);
        }
        return (sorted[^1], count, cumulative);
    }

    /// <summary>
    /// Asserts an <see cref="GenomeAssemblyAnalyzer.NxStatistics"/> is WELL-FORMED for a non-empty
    /// length set: all fields non-negative; Lx in [1, count]; Nx is one of the actual contig lengths
    /// and lies in [smallest, largest]; the cumulative is ≥ the inclusive cutoff and ≤ total.
    /// </summary>
    private static void AssertWellFormedNx(
        GenomeAssemblyAnalyzer.NxStatistics nx, IReadOnlyList<int> lengths, int threshold)
    {
        long total = lengths.Sum(l => (long)l);
        nx.Threshold.Should().Be(threshold);
        nx.Nx.Should().BeGreaterThanOrEqualTo(0, "a length is never negative");
        nx.Lx.Should().BeInRange(1, lengths.Count, "Lx is a contig count in [1, n] (INV-02)");
        nx.CumulativeLength.Should().BeInRange(1, total, "cumulative is bounded by the assembly total");

        if (total > 0)
        {
            nx.Nx.Should().BeInRange(lengths.Min(), lengths.Max(),
                "Nx is the length of an actual contig (INV-01)");
            lengths.Should().Contain(nx.Nx, "Nx equals one of the input contig lengths");
            (nx.CumulativeLength * 100).Should().BeGreaterThanOrEqualTo(total * threshold,
                "INV-05: the selected prefix's cumulative reaches the inclusive ≥x% cutoff");
        }
    }

    /// <summary>
    /// Asserts a full <see cref="GenomeAssemblyAnalyzer.AssemblyStatistics"/> record is WELL-FORMED for
    /// the given gap-free contig lengths: all stats non-negative; TotalLength = Σ lengths; counts and
    /// largest/smallest match; N50 ≤ LargestContig and SmallestContig ≤ N50; L50 in [1, count];
    /// N90 ≤ N50, L90 ≥ L50 (INV-03); mean = total/count.
    /// </summary>
    private static void AssertWellFormedStats(
        GenomeAssemblyAnalyzer.AssemblyStatistics stats, IReadOnlyList<int> lengths)
    {
        long total = lengths.Sum(l => (long)l);

        stats.TotalSequences.Should().Be(lengths.Count);
        stats.TotalLength.Should().Be(total, "TotalLength = Σ contig lengths (the assembly size)");
        stats.N50.Should().BeGreaterThanOrEqualTo(0);
        stats.L50.Should().BeInRange(1, lengths.Count, "L50 is a contig count in [1, n] (INV-02)");
        stats.LargestContig.Should().Be(lengths.Max());
        stats.SmallestContig.Should().Be(lengths.Min());
        stats.N50.Should().BeLessThanOrEqualTo(stats.LargestContig,
            "N50 ≤ the largest contig (it is one of the contig lengths)");
        stats.N50.Should().BeGreaterThanOrEqualTo(stats.SmallestContig,
            "N50 ≥ the smallest contig");

        // INV-03: Nx non-increasing, Lx non-decreasing in x.
        stats.N90.Should().BeLessThanOrEqualTo(stats.N50, "INV-03: N90 ≤ N50");
        stats.L90.Should().BeGreaterThanOrEqualTo(stats.L50, "INV-03: L90 ≥ L50");

        stats.MeanLength.Should().BeApproximately(total / (double)lengths.Count, 1e-9,
            "mean = total / count (no DivideByZero for a non-empty assembly)");
    }

    #endregion

    #region ASSEMBLY-STATS-001 — Assembly Statistics (BE: empty, single contig, equal-length contigs)

    #region Positive sanity — hand-computed documented N50 / L50 / total / mean

    // Doc §7.1 worked example: lengths {80,70,50,40,30,20}, T=290, 50%·T=145.
    // Cumulative 80 (<145), 80+70=150 (≥145) → N50=70, L50=2.
    [Test]
    public void Nx_DocWorkedExample_N50Is70_L50Is2()
    {
        var lengths = new[] { 80, 70, 50, 40, 30, 20 };

        var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, N50ThresholdPercent);

        nx.Nx.Should().Be(70, "doc §7.1: cumulative first reaches ≥145 (50% of 290) at the 70-contig");
        nx.Lx.Should().Be(2, "L50 = the rank (2) where the cumulative first reaches ≥50% (INV-02)");
        GenomeAssemblyAnalyzer.CalculateN50(lengths).Should().Be(70, "CalculateN50 == CalculateNx(·,50).Nx");
        AssertWellFormedNx(nx, lengths, N50ThresholdPercent);
    }

    // Prompt-specified hand worked set: {2,2,2,3,4,8,8} → total 29, half 14.5, N50 = 8.
    // Sorted desc {8,8,4,3,2,2,2}: cum 8 (<14.5), 16 (≥14.5) → N50=8, L50=2.
    [Test]
    public void Nx_HandComputedSet_N50Is8_L50Is2()
    {
        var lengths = new[] { 2, 2, 2, 3, 4, 8, 8 }; // total 29

        var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, N50ThresholdPercent);

        nx.Nx.Should().Be(8, "cumulative 8 (<14.5) then 16 (≥14.5) → N50 = 8");
        nx.Lx.Should().Be(2, "the cumulative first reaches ≥50% at rank 2");
        nx.CumulativeLength.Should().Be(16);
        AssertWellFormedNx(nx, lengths, N50ThresholdPercent);
    }

    // auN for the doc example: Σl²/Σl = 16700/290 ≈ 57.5862 (INV-04, §7.1).
    [Test]
    public void AuN_DocWorkedExample_MatchesSumOfSquaresOverTotal()
    {
        var lengths = new[] { 80, 70, 50, 40, 30, 20 };

        double aun = GenomeAssemblyAnalyzer.CalculateAuN(lengths);

        aun.Should().BeApproximately(16700.0 / 290.0, 1e-9, "auN = Σl²/Σl (INV-04, §7.1)");
    }

    // Full-record sanity over the doc set: total 290, count 6, mean 290/6, largest 80, smallest 20.
    [Test]
    public void Statistics_DocSet_TotalCountMeanN50AllDocumented()
    {
        var lengths = new[] { 80, 70, 50, 40, 30, 20 };
        var contigs = lengths.Select((l, i) => ($"c{i}", new string('A', l))).ToArray();

        var stats = GenomeAssemblyAnalyzer.CalculateStatistics(contigs);

        stats.TotalSequences.Should().Be(6);
        stats.TotalLength.Should().Be(290);
        stats.N50.Should().Be(70);
        stats.L50.Should().Be(2);
        stats.LargestContig.Should().Be(80);
        stats.SmallestContig.Should().Be(20);
        stats.MeanLength.Should().BeApproximately(290.0 / 6.0, 1e-9);
        AssertWellFormedStats(stats, lengths);
    }

    #endregion

    #region BE — EMPTY ASSEMBLY (all-zero record, no DivideByZero)

    // Empty assembly → all-zero record; mean = total/count must NOT throw DivideByZero (§6.1, dev #1).
    [Test]
    public void Statistics_EmptyAssembly_AllZeroRecord_NoDivideByZero()
    {
        Action act = () => GenomeAssemblyAnalyzer.CalculateStatistics(Array.Empty<(string, string)>());
        act.Should().NotThrow("an empty assembly is handled without throwing (§3.3, §6.1)");

        var stats = GenomeAssemblyAnalyzer.CalculateStatistics(Array.Empty<(string, string)>());

        stats.TotalSequences.Should().Be(0);
        stats.TotalLength.Should().Be(0);
        stats.N50.Should().Be(0, "no defined N50 for an empty assembly → 0 (§6.1)");
        stats.L50.Should().Be(0);
        stats.N90.Should().Be(0);
        stats.L90.Should().Be(0);
        stats.LargestContig.Should().Be(0);
        stats.SmallestContig.Should().Be(0);
        stats.MeanLength.Should().Be(0, "mean is 0, not NaN/∞: no division by a zero count");
        double.IsNaN(stats.MeanLength).Should().BeFalse("DivideByZero would surface as NaN");
    }

    // Empty length set → CalculateNx returns Nx=Lx=0 (no enumeration, no crash) (§3.3).
    [Test]
    public void Nx_EmptyLengths_ReturnsZeroNxZeroLx()
    {
        var nx = GenomeAssemblyAnalyzer.CalculateNx(Array.Empty<int>(), N50ThresholdPercent);

        nx.Nx.Should().Be(0);
        nx.Lx.Should().Be(0);
        nx.CumulativeLength.Should().Be(0);
        GenomeAssemblyAnalyzer.CalculateN50(Array.Empty<int>()).Should().Be(0);
    }

    // All-zero-length contigs make total = 0; the total==0 guard must avoid divide-by-zero in the
    // ≥x% comparison and return Nx=Lx=0 (§6.1).
    [Test]
    public void Nx_AllZeroLengthContigs_TotalZero_ReturnsZero()
    {
        var lengths = new[] { 0, 0, 0 }; // total 0

        var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, N50ThresholdPercent);

        nx.Nx.Should().Be(0, "total 0 → guarded Nx=Lx=0, no division by zero (§6.1)");
        nx.Lx.Should().Be(0);
        GenomeAssemblyAnalyzer.CalculateAuN(lengths).Should().Be(0, "auN guards total==0 → 0");
    }

    // Empty assembly → auN = 0 (§6.1).
    [Test]
    public void AuN_EmptyAssembly_ReturnsZero()
    {
        GenomeAssemblyAnalyzer.CalculateAuN(Array.Empty<int>()).Should().Be(0);
    }

    #endregion

    #region BE — SINGLE CONTIG (N50 = its length, L50 = 1, the whole assembly)

    // A single contig of length L → N50=L, L50=1: it is the whole assembly, reaching 100% ≥ 50%
    // at rank 1 (§6.1 "Single contig").
    [Test]
    public void Nx_SingleContig_N50IsItsLength_L50IsOne()
    {
        var lengths = new[] { 137 };

        var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, N50ThresholdPercent);

        nx.Nx.Should().Be(137, "the lone contig IS the assembly (§6.1)");
        nx.Lx.Should().Be(1, "rank 1 already covers 100% ≥ 50%");
        nx.CumulativeLength.Should().Be(137);
        AssertWellFormedNx(nx, lengths, N50ThresholdPercent);
    }

    // Single contig full record: total = mean = largest = smallest = N50 = L; L50 = 1; no crash.
    [Test]
    public void Statistics_SingleContig_AllStatsEqualItsLength()
    {
        var contigs = new[] { ("only", new string('G', 250)) };

        var stats = GenomeAssemblyAnalyzer.CalculateStatistics(contigs);

        stats.TotalSequences.Should().Be(1);
        stats.TotalLength.Should().Be(250);
        stats.N50.Should().Be(250);
        stats.L50.Should().Be(1);
        stats.N90.Should().Be(250, "with one contig, every Nx is that contig's length");
        stats.L90.Should().Be(1);
        stats.LargestContig.Should().Be(250);
        stats.SmallestContig.Should().Be(250);
        stats.MeanLength.Should().Be(250, "mean of one contig = its length, no DivideByZero");
        AssertWellFormedStats(stats, new[] { 250 });
    }

    // Fuzz: a single random contig over many trials always yields N50 = its length, L50 = 1.
    [Test]
    [CancelAfter(10_000)]
    public void Nx_RandomSingleContig_AlwaysItsLengthRankOne()
    {
        var rng = new Random(147_001);
        for (int trial = 0; trial < 1000; trial++)
        {
            int len = rng.Next(1, 100_000);
            var lengths = new[] { len };

            var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, N50ThresholdPercent);

            nx.Nx.Should().Be(len, "a single contig's N50 is its own length");
            nx.Lx.Should().Be(1, "a single contig's L50 is rank 1");
            AssertWellFormedNx(nx, lengths, N50ThresholdPercent);
        }
    }

    #endregion

    #region BE — EQUAL-LENGTH CONTIGS (N50 = L, L50 = exact tie rank ⌈k/2⌉)

    // Two equal contigs {L, L}: T = 2L, 50%·T = L. Cumulative L (= L, ≥L inclusive) at rank 1
    // → N50=L, L50=1 (INV-05 inclusive boundary: exactly 50% selects that contig).
    [Test]
    public void Nx_TwoEqualContigs_FirstAlreadyReachesExactly50Percent()
    {
        var lengths = new[] { 10, 10 }; // T=20, 50%=10, cum 10 ≥ 10 at rank 1

        var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, N50ThresholdPercent);

        nx.Nx.Should().Be(10);
        nx.Lx.Should().Be(1, "INV-05: cumulative exactly = 50% selects that contig (inclusive)");
        nx.CumulativeLength.Should().Be(10);
        AssertWellFormedNx(nx, lengths, N50ThresholdPercent);
    }

    // Three equal contigs {L,L,L}: T=3L, 50%·T=1.5L. Cumulative L (<1.5L), 2L (≥1.5L) at rank 2
    // → N50=L, L50=2 = ⌈3/2⌉ (the documented tie rank).
    [Test]
    public void Nx_ThreeEqualContigs_L50IsCeilHalf()
    {
        var lengths = new[] { 7, 7, 7 }; // T=21, half=10.5, cum 7 (<10.5), 14 (≥10.5) rank 2

        var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, N50ThresholdPercent);

        nx.Nx.Should().Be(7, "all contigs equal → N50 = L");
        nx.Lx.Should().Be(2, "L50 = ⌈3/2⌉ = 2 where cumulative first reaches ≥50% (INV-02)");
        AssertWellFormedNx(nx, lengths, N50ThresholdPercent);
    }

    // Equal-length sweep: for k equal contigs of length L, N50 = L and L50 = ⌈k/2⌉ exactly.
    // ⌈k/2⌉ is the smallest k' with k'·L ≥ ⌈(k·L)/2⌉, i.e. k'*2 ≥ k (integer-exact INV-05).
    [Test]
    [CancelAfter(15_000)]
    public void Nx_EqualLengthContigs_N50IsL_L50IsCeilHalf_AllK()
    {
        var rng = new Random(147_002);
        foreach (int L in new[] { 1, 5, 100, 9999 })
        {
            for (int k = 1; k <= 64; k++)
            {
                var lengths = Enumerable.Repeat(L, k).ToArray();

                var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, N50ThresholdPercent);

                int expectedL50 = (k + 1) / 2; // ⌈k/2⌉
                nx.Nx.Should().Be(L, $"k={k} equal contigs of length {L} → N50 = L");
                nx.Lx.Should().Be(expectedL50,
                    $"k={k}: L50 = ⌈k/2⌉ = {expectedL50} (exact tie rank, INV-02/INV-05)");

                // Cross-check the full statistics record with shuffled contig order (sort-independence).
                var shuffled = lengths.OrderBy(_ => rng.Next()).ToArray();
                var contigs = ContigsOfLengths(rng, shuffled);
                var stats = GenomeAssemblyAnalyzer.CalculateStatistics(contigs);
                stats.N50.Should().Be(L, "N50 is order-independent (sorted internally)");
                stats.L50.Should().Be(expectedL50, "L50 is order-independent");
                stats.TotalLength.Should().Be((long)L * k);
                stats.MeanLength.Should().BeApproximately(L, 1e-9, "equal contigs → mean = L");
                AssertWellFormedStats(stats, lengths);
            }
        }
    }

    #endregion

    #region BE — inclusive ≥50% boundary (off-by-one guard on odd totals & exact-50% landings)

    // Odd total: {3,3,3}, T=9, 50%·T=4.5. Cumulative 3 (<4.5), 6 (≥4.5) at rank 2 → N50=3, L50=2.
    // Integer-exact test `cum*100 ≥ T*50`: 3*100=300 < 9*50=450; 6*100=600 ≥ 450. No off-by-one.
    [Test]
    public void Nx_OddTotal_InclusiveIntegerExactBoundary()
    {
        var lengths = new[] { 3, 3, 3 };

        var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, N50ThresholdPercent);

        nx.Nx.Should().Be(3);
        nx.Lx.Should().Be(2, "integer-exact ≥50% test selects rank 2, not 1 (no off-by-one, §4.2)");
        AssertWellFormedNx(nx, lengths, N50ThresholdPercent);
    }

    // A contig landing EXACTLY on the ≥50% cutoff must be selected (INV-05 inclusive).
    // {6,4} T=10, 50%=5: cum 6 ≥ 5 at rank 1 → N50=6, L50=1.
    // {5,5} T=10, 50%=5: cum 5 = 5 (inclusive) at rank 1 → N50=5, L50=1.
    [Test]
    public void Nx_ExactlyAtFiftyPercent_InclusiveSelectsThatContig()
    {
        var a = GenomeAssemblyAnalyzer.CalculateNx(new[] { 6, 4 }, N50ThresholdPercent);
        a.Nx.Should().Be(6);
        a.Lx.Should().Be(1, "cum 6 ≥ 5 at rank 1");

        var b = GenomeAssemblyAnalyzer.CalculateNx(new[] { 5, 5 }, N50ThresholdPercent);
        b.Nx.Should().Be(5);
        b.Lx.Should().Be(1, "INV-05: cumulative EXACTLY 50% is inclusive → selected at rank 1");
    }

    // INV-03 boundary on a skewed set: one big + many tiny contigs. Big contig alone exceeds 50%,
    // so N50 = big, L50 = 1; N90 ≤ N50 and L90 ≥ L50.
    [Test]
    public void Nx_OneDominantContig_N50IsDominant_L50One()
    {
        var lengths = new[] { 1000, 1, 1, 1, 1 }; // T=1004, 50%=502, big alone 1000 ≥ 502

        var n50 = GenomeAssemblyAnalyzer.CalculateNx(lengths, 50);
        var n90 = GenomeAssemblyAnalyzer.CalculateNx(lengths, 90);

        n50.Nx.Should().Be(1000);
        n50.Lx.Should().Be(1, "the dominant contig alone covers ≥50%");
        n90.Nx.Should().BeLessThanOrEqualTo(n50.Nx, "INV-03: N90 ≤ N50");
        n90.Lx.Should().BeGreaterThanOrEqualTo(n50.Lx, "INV-03: L90 ≥ L50");
    }

    #endregion

    #region BE — large random sweep (well-formed, matches the definition oracle, deterministic)

    // Fuzz: arbitrary random length sets NEVER throw, are DETERMINISTIC, match the §2.2 definition
    // oracle for both N50 and N90, and produce a well-formed record (total = Σ, N50 ≤ max, L50 ≤ n).
    [Test]
    [CancelAfter(60_000)]
    public void Statistics_RandomLengthSets_WellFormed_MatchOracle_Deterministic()
    {
        var rng = new Random(147_003);
        for (int trial = 0; trial < 3000; trial++)
        {
            int n = rng.Next(1, 40);
            var lengths = new int[n];
            for (int i = 0; i < n; i++)
            {
                // Mix tiny boundary lengths (0/1) with larger spreads to hit odd totals and ties.
                lengths[i] = rng.Next(6) switch
                {
                    0 => 0,
                    1 => 1,
                    2 => rng.Next(1, 5),       // many small → frequent ties
                    _ => rng.Next(1, 100_000), // large spread
                };
            }

            // CalculateNx(lengths,·) must match the independent definition oracle exactly.
            foreach (int x in new[] { 50, 90 })
            {
                var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, x);
                var (oNx, oLx, oCum) = OracleNx(lengths, x);
                nx.Nx.Should().Be(oNx, "Nx matches the §2.2 definition oracle");
                nx.Lx.Should().Be(oLx, "Lx matches the §2.2 definition oracle");
                nx.CumulativeLength.Should().Be(oCum);

                // Determinism: a repeat call is identical.
                GenomeAssemblyAnalyzer.CalculateNx(lengths, x).Should().Be(nx,
                    "Nx is a deterministic pure function");
            }

            long total = lengths.Sum(l => (long)l);
            if (total == 0) continue; // all-zero handled by the empty/zero-total tests

            var contigs = ContigsOfLengths(rng, lengths);
            var stats = GenomeAssemblyAnalyzer.CalculateStatistics(contigs);
            AssertWellFormedStats(stats, lengths);
            AssertWellFormedNx(
                GenomeAssemblyAnalyzer.CalculateNx(lengths, 50), lengths, 50);
        }
    }

    #endregion

    #endregion
}
