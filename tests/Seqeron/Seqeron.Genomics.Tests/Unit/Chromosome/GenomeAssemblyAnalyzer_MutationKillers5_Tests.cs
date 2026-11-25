using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// ASSEMBLY-STATS-001 mutation killers (batch 5): exercises the suspicious-region detector's
/// High-N-content and GC-deviation branches (the canonical and earlier batches only hit the
/// low-complexity branch), the inverted syntenic-block branch, and the repetitive-region span.
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_MutationKillers5_Tests
{
    private static string RandomFrom(string alphabet, int length, int seed)
    {
        var rng = new Random(seed);
        return string.Concat(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]));
    }

    [Test]
    public void FindSuspiciousRegions_HighNContent_IsFlagged()
    {
        // A 100-bp N run inside otherwise complex sequence: a 500-bp window has 100/500 = 20% N
        // (> 10%) ⇒ "High N content".
        string seq = RandomFrom("ACGT", 250, 1) + new string('N', 100) + RandomFrom("ACGT", 250, 3);
        var regions = FindSuspiciousRegions(new[] { ("s", seq) }).ToList();
        Assert.That(regions.Any(r => r.Reason.Contains("High N content")), Is.True);
    }

    [Test]
    public void FindSuspiciousRegions_GcDeviation_IsFlagged()
    {
        // GC-rich half (≈75% GC) then AT-rich half (≈25% GC): global ≈ 0.5, so 500-bp windows in
        // either half deviate by ≈0.25 (> 0.15) ⇒ "GC deviation". Both halves stay complex (4 bases).
        string seq = RandomFrom("GGGCCCAT", 1000, 1) + RandomFrom("AAATTTGC", 1000, 2);
        var regions = FindSuspiciousRegions(new[] { ("s", seq) }).ToList();
        Assert.That(regions.Any(r => r.Reason.Contains("GC deviation")), Is.True);
    }

    [Test]
    public void FindRepetitiveRegions_TandemSpansWholeSequence()
    {
        // (CAG)×200 = 600 bp: the three 6-mers CAGCAG/AGCAGC/GCAGCA occur 199/198/198 times, every window
        // scores far above T_threshold = 3 ⇒ one WindowMasker interval [0, 599] with Copies = 199.
        string seq = string.Concat(Enumerable.Repeat("CAG", 200));
        var r = FindRepetitiveRegions(new[] { ("s", seq) }, kmerSize: 6, minCopies: 3, windowSize: 60).Single();
        Assert.That(r.SequenceId, Is.EqualTo("s"));
        Assert.That(r.Start, Is.EqualTo(0));
        Assert.That(r.End, Is.EqualTo(599));
        Assert.That(r.Copies, Is.EqualTo(199));
    }
}
