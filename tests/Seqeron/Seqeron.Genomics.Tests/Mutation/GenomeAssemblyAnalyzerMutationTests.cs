using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// ASSEMBLY-STATS-001 mutation killers: exact-value tests for the assembly-statistics helper methods
/// that the canonical AssemblyStatistics fixture left uncovered — N50/Nx/Lx/auN contiguity metrics,
/// gap finding/classification/distribution, contig extraction, and the length filter/sort/distribution
/// utilities. Each assertion reproduces the published definition (Miller, Koren &amp; Sutton 2010 §1.2;
/// QUAST N50; auN = Σℓ²/Σℓ) so an operator/threshold mutation diverges.
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzerMutationTests
{
    private const double Tol = 1e-9;
    private static readonly int[] Lengths = { 100, 200, 300, 400, 500 }; // total 1500

    #region N50 / Nx / Lx / auN

    [Test]
    public void CalculateN50_AndNx_ShortestContigReaching50Percent()
    {
        // sorted desc 500,400,…; cumulative reaches ≥50%·1500=750 at the 400-contig (cum 900).
        Assert.That(CalculateN50(Lengths), Is.EqualTo(400));

        var n50 = CalculateNx(Lengths, 50);
        Assert.That(n50.Nx, Is.EqualTo(400));
        Assert.That(n50.Lx, Is.EqualTo(2));
        Assert.That(n50.CumulativeLength, Is.EqualTo(900));

        // N90 needs ≥90%·1500=1350 ⇒ reached at the 200-contig (cum 1400).
        var n90 = CalculateNx(Lengths, 90);
        Assert.That(n90.Nx, Is.EqualTo(200));
        Assert.That(n90.Lx, Is.EqualTo(4));

        Assert.That(CalculateN50(Array.Empty<int>()), Is.EqualTo(0));
    }

    [Test]
    public void CalculateNxCurve_OrdersByThresholdAscending()
    {
        var curve = CalculateNxCurve(Lengths, 90, 50).ToList(); // passed out of order
        Assert.That(curve.Select(c => c.Threshold), Is.EqualTo(new[] { 50, 90 }));
        Assert.That(curve[0].Nx, Is.EqualTo(400));
        Assert.That(curve[1].Nx, Is.EqualTo(200));
    }

    [Test]
    public void CalculateAuN_SumOfSquaresOverTotal()
    {
        // auN = Σ ℓ² / Σ ℓ = (500²+400²+300²+200²+100²)/1500 = 550000/1500.
        Assert.That(CalculateAuN(Lengths), Is.EqualTo(550000.0 / 1500.0).Within(Tol));
        Assert.That(CalculateAuN(Array.Empty<int>()), Is.EqualTo(0.0).Within(Tol));
    }

    #endregion

    #region Gaps

    [Test]
    public void FindGaps_LocatesNRunsAndClassifiesByLength()
    {
        var gaps = FindGaps(new[] { ("s1", "AAANNNAAA") }).ToList();
        Assert.That(gaps, Has.Count.EqualTo(1));
        Assert.That(gaps[0].SequenceId, Is.EqualTo("s1"));
        Assert.That(gaps[0].Start, Is.EqualTo(3));
        Assert.That(gaps[0].End, Is.EqualTo(5));
        Assert.That(gaps[0].Length, Is.EqualTo(3));
        Assert.That(gaps[0].GapType, Is.EqualTo("Short")); // <10

        // Length-class boundaries: 10→Medium, 100→Long, 1000→Scaffold.
        Assert.That(FindGaps(new[] { ("s", "A" + new string('N', 10) + "A") }).Single().GapType, Is.EqualTo("Medium"));
        Assert.That(FindGaps(new[] { ("s", "A" + new string('N', 100) + "A") }).Single().GapType, Is.EqualTo("Long"));
        Assert.That(FindGaps(new[] { ("s", "A" + new string('N', 1000) + "A") }).Single().GapType, Is.EqualTo("Scaffold"));
    }

    [Test]
    public void FindGaps_RespectsMinGapLength()
    {
        Assert.That(FindGaps(new[] { ("s", "ANA") }, minGapLength: 2).ToList(), Is.Empty); // gap len 1 < 2
        Assert.That(FindGaps(new[] { ("s", "ANNA") }, minGapLength: 2).ToList(), Has.Count.EqualTo(1));
    }

    [Test]
    public void AnalyzeGapDistribution_ExactSummary()
    {
        var gaps = FindGaps(new[] { ("s", "NNAANNNN") }).ToList(); // gaps of length 2 and 4
        var dist = AnalyzeGapDistribution(gaps);
        Assert.That(dist.Count, Is.EqualTo(2));
        Assert.That(dist.MeanLength, Is.EqualTo(3.0).Within(Tol));
        Assert.That(dist.MedianLength, Is.EqualTo(4.0).Within(Tol)); // upper-mid at index count/2
        Assert.That(dist.MaxLength, Is.EqualTo(4));
        Assert.That(dist.TypeCounts["Short"], Is.EqualTo(2));
    }

    #endregion

    #region Contig extraction & length utilities

    [Test]
    public void ExtractContigs_SplitsOnNRunsAndFiltersByLength()
    {
        var contigs = ExtractContigs(new[] { ("s", "AAANNAAAA") }, minContigLength: 3).ToList();
        Assert.That(contigs, Is.EqualTo(new[] { ("s_contig1", "AAA"), ("s_contig2", "AAAA") }));

        // Raising the threshold drops the 3-bp contig; the 4-bp contig is renumbered from 1.
        var filtered = ExtractContigs(new[] { ("s", "AAANNAAAA") }, minContigLength: 4).ToList();
        Assert.That(filtered, Is.EqualTo(new[] { ("s_contig1", "AAAA") }));
    }

    [Test]
    public void FilterByLength_InclusiveBounds()
    {
        var seqs = new[] { ("a", "AAA"), ("b", "AAAAA"), ("c", "A") };
        Assert.That(FilterByLength(seqs, minLength: 3).Select(s => s.Id), Is.EqualTo(new[] { "a", "b" }));
        Assert.That(FilterByLength(seqs, minLength: 3, maxLength: 4).Select(s => s.Id), Is.EqualTo(new[] { "a" }));
    }

    [Test]
    public void SortByLength_DescendingAndAscending()
    {
        var seqs = new[] { ("a", "AAA"), ("b", "AAAAA"), ("c", "A") };
        Assert.That(SortByLength(seqs).Select(s => s.Id), Is.EqualTo(new[] { "b", "a", "c" }));
        Assert.That(SortByLength(seqs, descending: false).Select(s => s.Id), Is.EqualTo(new[] { "c", "a", "b" }));
    }

    [Test]
    public void CalculateLengthDistribution_BinsByFirstUpperBound()
    {
        var dist = CalculateLengthDistribution(new[] { 50, 150, 600 }, 100, 500, 1000);
        Assert.That(dist["<100"], Is.EqualTo(1));
        Assert.That(dist["<500"], Is.EqualTo(1));
        Assert.That(dist["<1000"], Is.EqualTo(1));
        Assert.That(dist[">=1000"], Is.EqualTo(0));
    }

    #endregion
}
