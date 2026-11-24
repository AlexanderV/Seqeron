using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// B09 WP78 R14-19 (F361): the result records of the ORF, repeat and common-region searches keep the text and the coordinates and cut their string out on first access
/// (<see cref="OrfInfo.Sequence"/>, <see cref="RepeatInfo.Sequence"/>, <see cref="CommonRegion.Sequence"/>; <c>SequenceMemory</c> is the view without a copy). The values, the equality and the order of the
/// results are those of the eager records; the nested-ORF input of the worklist (ATG x 8000 + TAA, minLength 3: 8 000 ORFs of 3 ... 24 003 nucleotides, 96 036 000 characters) is enumerated without building them.
/// </summary>
[TestFixture]
[Category("GENOMIC-ORF-001")]
[Category("GENOMIC-REPEAT-001")]
public class GenomicAnalyzer_LazyResults_Tests
{
    private static string Nested(int copies) => string.Concat(Enumerable.Repeat("ATG", copies)) + "TAA";

    [Test]
    public void NestedOrfs_AreEnumeratedWithoutBuildingTheirStrings()
    {
        var dna = new DnaSequence(Nested(8000));
        GC.Collect();
        long before = GC.GetAllocatedBytesForCurrentThread();
        int count = 0;
        long totalLength = 0;
        foreach (var orf in GenomicAnalyzer.FindOpenReadingFrames(dna, 3))
        {
            count++;
            totalLength += orf.Length;
        }

        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(count, Is.EqualTo(8000));
        Assert.That(totalLength, Is.EqualTo(96_036_000L), "the sum of the lengths of the ORFs (3 * (k + 1) for k = 1 ... 8000)");
        Assert.That(allocated, Is.LessThan(40L * 1024 * 1024), "the eager strings were 183 MB (96 036 000 characters)");
    }

    [Test]
    public void OrfSequence_IsTheSubstringOfTheScannedStrand_AndEqualsTheEagerRecord()
    {
        var dna = new DnaSequence("CCATGAAATGTAACATGGCCTGAGG");
        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, 3).ToList();
        Assert.That(orfs, Is.Not.Empty);
        foreach (var orf in orfs)
        {
            string strand = orf.IsReverseComplement ? DnaSequence.GetReverseComplementString(dna.Sequence) : dna.Sequence;
            Assert.That(orf.Sequence, Is.EqualTo(strand.Substring(orf.Position, orf.Length)));
            Assert.That(orf.SequenceMemory.ToString(), Is.EqualTo(orf.Sequence));
            string built = orf.Sequence;
            Assert.That(orf.Sequence, Is.SameAs(built), "built once");
            Assert.That(orf.CodonCount, Is.EqualTo(orf.Length / 3));
            var eager = new OrfInfo(orf.Sequence, orf.Position, orf.Frame, orf.IsReverseComplement, orf.IsTerminated, orf.RegionStart, orf.RegionLength);
            Assert.That(orf, Is.EqualTo(eager));
            Assert.That(orf.GetHashCode(), Is.EqualTo(eager.GetHashCode()));
            Assert.That(orf.GetForwardStrandRange(), Is.EqualTo(eager.GetForwardStrandRange()));
        }

        var again = GenomicAnalyzer.FindOpenReadingFrames(dna, 3).ToList();
        Assert.That(again, Is.EqualTo(orfs), "two enumerations give equal records");
    }

    [Test]
    public void OrfEquality_DistinguishesTheSequenceAndTheCoordinates()
    {
        var a = new OrfInfo("ATGTAA", 0, 1, false);
        Assert.That(a, Is.EqualTo(new OrfInfo("ATGTAA", 0, 1, false)));
        Assert.That(a, Is.Not.EqualTo(new OrfInfo("ATGTGA", 0, 1, false)));
        Assert.That(a, Is.Not.EqualTo(new OrfInfo("ATGTAA", 3, 1, false)));
        Assert.That(a, Is.Not.EqualTo(new OrfInfo("ATGTAA", 0, 1, true)));
        Assert.That(a.Length, Is.EqualTo(6));
    }

    [Test]
    public void RepeatSequences_AreTheLabelPrefixes_AndEqualTheEagerRecords()
    {
        var dna = new DnaSequence("ACGTACGTACGTTTACGTACGT");
        var repeats = GenomicAnalyzer.FindRepeats(dna, 2).ToList();
        Assert.That(repeats, Is.Not.Empty);
        foreach (var repeat in repeats)
        {
            Assert.That(repeat.Length, Is.EqualTo(repeat.Sequence.Length));
            Assert.That(repeat.SequenceMemory.ToString(), Is.EqualTo(repeat.Sequence));
            Assert.That(repeat.IsEmpty, Is.False);
            foreach (int p in repeat.Positions)
                Assert.That(dna.Sequence.AsSpan(p, repeat.Length).SequenceEqual(repeat.Sequence), Is.True, repeat.Sequence + " at " + p);
            var eager = new RepeatInfo(repeat.Sequence, repeat.Positions);
            Assert.That(repeat, Is.EqualTo(eager));
            Assert.That(repeat.GetHashCode(), Is.EqualTo(eager.GetHashCode()));
        }

        Assert.That(RepeatInfo.None.IsEmpty, Is.True);
        Assert.That(RepeatInfo.None.Length, Is.EqualTo(0));
    }

    [Test]
    public void CommonRegionSequences_AreTheQuerySubstrings_AndEqualTheEagerRecords()
    {
        var a = new DnaSequence("ACGTACGTTTGACCAGTA");
        var b = new DnaSequence("TTGACCAGTAACGTACGG");
        var regions = GenomicAnalyzer.FindCommonRegions(a, b, 4).ToList();
        Assert.That(regions, Is.Not.Empty);
        foreach (var region in regions)
        {
            Assert.That(a.Sequence.Substring(region.PositionInFirst, region.Length), Is.EqualTo(region.Sequence));
            Assert.That(region.SequenceMemory.ToString(), Is.EqualTo(region.Sequence));
            var eager = new CommonRegion(region.Sequence, region.PositionInFirst, region.PositionInSecond, region.IsReverseComplement);
            Assert.That(region, Is.EqualTo(eager));
            Assert.That(region.GetHashCode(), Is.EqualTo(eager.GetHashCode()));
        }

        Assert.That(CommonRegion.None.IsEmpty, Is.True);
        Assert.That(CommonRegion.None.Equals(new CommonRegion(string.Empty, -1, -1)), Is.True);
    }
}
