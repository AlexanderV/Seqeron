using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// GENOMIC-MOTIFS-001 (audit R14-3, B09 F290): <see cref="GenomicAnalyzer.FindMotif(DnaSequence, string)"/> / <see cref="MotifFinder.FindExactMotif"/> no longer build the suffix tree of the sequence for one query
/// (4.37 s / 542 MB at 2 Mbp against 92 ms / 7 MB for the linear scan). A sequence whose tree is already built (<see cref="DnaSequence.HasSuffixTree"/>) keeps the O(m + occ) route; both give the positions of
/// a brute-force scan, ascending, overlapping occurrences included. <see cref="DnaSequence.ReleaseSuffixTree"/> frees the cached tree.
/// </summary>
[TestFixture]
[Category("GENOMIC-MOTIFS-001")]
public class GenomicAnalyzer_FindMotifIndex_Tests
{
    private static List<int> BruteForce(string text, string motif)
    {
        var list = new List<int>();
        for (int i = 0; i + motif.Length <= text.Length; i++)
        {
            if (string.CompareOrdinal(text, i, motif, 0, motif.Length) == 0)
                list.Add(i);
        }

        return list;
    }

    [Test]
    public void FindMotif_WithoutATree_DoesNotBuildOne_AndEqualsBruteForce_AndTheTreeRoute()
    {
        var rnd = new Random(14003);
        for (int t = 0; t < 120; t++)
        {
            int n = rnd.Next(1, 300);
            string text = string.Concat(Enumerable.Range(0, n).Select(_ => "ACGT"[rnd.Next(rnd.Next(2) == 0 ? 2 : 4)]));
            string motif = rnd.Next(3) == 0
                ? string.Concat(Enumerable.Range(0, rnd.Next(1, 8)).Select(_ => "ACGT"[rnd.Next(4)]))
                : text.Substring(rnd.Next(n), 1).PadRight(rnd.Next(1, Math.Min(8, n) + 1), text[rnd.Next(n)]);
            var dna = new DnaSequence(text);
            var scan = GenomicAnalyzer.FindMotif(dna, motif);
            Assert.That(dna.HasSuffixTree, Is.False, "one query must not build the index");
            Assert.That(scan, Is.EqualTo(BruteForce(text, motif.ToUpperInvariant())), $"{text} / {motif}");
            _ = dna.SuffixTree;
            Assert.That(dna.HasSuffixTree, Is.True);
            Assert.That(GenomicAnalyzer.FindMotif(dna, motif), Is.EqualTo(scan), "the tree route gives the same ascending positions");
        }
    }

    [Test]
    public void FindMotif_LowerCaseMotif_IsCaseInsensitive_AndEmptyOrLongerMotifIsEmpty()
    {
        var dna = new DnaSequence("ACGTACGTAC");
        Assert.That(GenomicAnalyzer.FindMotif(dna, "acgt"), Is.EqualTo(new[] { 0, 4 }));
        Assert.That(GenomicAnalyzer.FindMotif(dna, "AC"), Is.EqualTo(new[] { 0, 4, 8 }));
        Assert.That(GenomicAnalyzer.FindMotif(dna, ""), Is.Empty);
        Assert.That(GenomicAnalyzer.FindMotif(dna, "ACGTACGTACG"), Is.Empty);
        Assert.That(GenomicAnalyzer.FindMotif(dna, "NNN"), Is.Empty, "an IUPAC code is a literal symbol");
        Assert.That(dna.HasSuffixTree, Is.False);
    }

    [Test]
    public void ReleaseSuffixTree_FreesTheTree_AndTheNextReadRebuildsIt()
    {
        var dna = new DnaSequence("ACGTACGTAC");
        dna.ReleaseSuffixTree();       // nothing built: no effect
        Assert.That(dna.HasSuffixTree, Is.False);
        var tree = dna.SuffixTree;
        Assert.That(dna.HasSuffixTree, Is.True);
        Assert.That(dna.SuffixTree, Is.SameAs(tree));
        dna.ReleaseSuffixTree();
        Assert.That(dna.HasSuffixTree, Is.False);
        Assert.That(GenomicAnalyzer.FindMotif(dna, "ACG"), Is.EqualTo(new[] { 0, 4 }));
        Assert.That(dna.SuffixTree, Is.Not.SameAs(tree));
    }

    [Test]
    public void FindMotif_OnOneMegabase_AllocatesFarLessThanTheTree()
    {
        // R14-3 measured the tree at about 170 MB per Mbp (542 MB at 2 Mbp); the scan allocates the upper-cased copy of the text and the result
        var rnd = new Random(14004);
        var chars = new char[1_000_000];
        for (int i = 0; i < chars.Length; i++)
            chars[i] = "ACGT"[rnd.Next(4)];
        var dna = new DnaSequence(new string(chars));
        long before = GC.GetAllocatedBytesForCurrentThread();
        var hits = GenomicAnalyzer.FindMotif(dna, "ACGTACGTAC");
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(allocated, Is.LessThan(30L * 1024 * 1024), $"allocated {allocated / 1048576} MB for one query on 1 Mbp");
        Assert.That(hits, Is.EqualTo(BruteForce(dna.Sequence, "ACGTACGTAC")));
        Assert.That(dna.HasSuffixTree, Is.False);
    }
}
