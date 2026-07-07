using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// ASSEMBLY-STATS-001 mutation killers (batch 6): pins the duplication-percent ratio with more than
/// one complete marker, the multi-window local-quality scan (WindowSize at a non-zero offset and the
/// exact window count), and the assembly-comparison disjoint case.
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_MutationKillers6_Tests
{
    private const double Tol = 1e-9;

    private static string RandomDna(int length, int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT";
        return string.Concat(Enumerable.Range(0, length).Select(_ => bases[rng.Next(4)]));
    }

    [Test]
    public void AssessCompleteness_SingleAndDuplicatedMarkers_ExactDuplicationPercent()
    {
        // g1 present once (single-copy), g2 present in two assembly sequences (duplicated):
        // complete = 2, duplicated = 1 ⇒ DuplicationPercent = 1·100/2 = 50; completeness = 100.
        string g1 = RandomDna(100, 21);
        string g2 = RandomDna(100, 22);
        var r = AssessCompleteness(
            new[] { ("a", g1), ("b", g2), ("c", g2) },
            new[] { ("g1", g1), ("g2", g2) });

        Assert.That(r.Complete, Is.EqualTo(2));
        Assert.That(r.CompleteSingleCopy, Is.EqualTo(1));
        Assert.That(r.CompleteDuplicated, Is.EqualTo(1));
        Assert.That(r.CompletenessPercent, Is.EqualTo(100.0).Within(Tol));
        Assert.That(r.DuplicationPercent, Is.EqualTo(50.0).Within(Tol)); // kills duplicated*100*complete
    }

    [Test]
    public void CalculateLocalQuality_MultipleWindows_ExactCountAndOffsetSize()
    {
        // 1500-bp sequence, default window 1000, step 500 ⇒ windows start at 0, 500, 1000 (exactly 3).
        // The window at offset 500 has size end−start = 1500−500 = 1000 (kills an 'end + i' mutant);
        // a 'i <= length' loop-bound mutant would emit a spurious 4th empty window at 1500.
        var windows = CalculateLocalQuality(new[] { ("s", RandomDna(1500, 9)) }).ToList();
        Assert.That(windows, Has.Count.EqualTo(3));
        Assert.That(windows[1].Position, Is.EqualTo(500));
        Assert.That(windows[1].WindowSize, Is.EqualTo(1000));
        Assert.That(windows[2].Position, Is.EqualTo(1000));
        Assert.That(windows[2].WindowSize, Is.EqualTo(500)); // 1500 − 1000
    }

    [Test]
    public void CompareAssemblies_DisjointHasZeroAlignedFractions()
    {
        var cmp = CompareAssemblies(
            new[] { ("s", RandomDna(60, 1)) },
            new[] { ("s", RandomDna(60, 2)) },
            kmerSize: 21);
        Assert.That(cmp.AlignedFraction1, Is.EqualTo(0.0).Within(Tol));
        Assert.That(cmp.AlignedFraction2, Is.EqualTo(0.0).Within(Tol));
        Assert.That(cmp.SequenceIdentity, Is.EqualTo(0.0).Within(Tol));
    }
}
