using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class AssembleDeBruijnTests
{
    [Test]
    public void AssembleDeBruijn_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.AssembleDeBruijn(new[] { "AAACCC", "AACCCG" }, 20, 0.9, 3, 3));
        // Empty read set yields an empty assembly (no throw).
        var empty = AlignmentTools.AssembleDeBruijn(Array.Empty<string>(), 20, 0.9, 3, 3);
        Assert.That(empty.Contigs, Is.Empty);
        Assert.That(empty.TotalLength, Is.EqualTo(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.AssembleDeBruijn(new[] { "ACGT" }, 20, 0.9, 1, 1));
    }

    [Test]
    public void AssembleDeBruijn_Binding_InvokesSuccessfully()
    {
        // Overlapping reads share the 4-mers GGCG/GCGT, counted once: one non-branching path (single-strand model).
        var r = AlignmentTools.AssembleDeBruijn(new[] { "ATGGCGT", "GGCGTGCA" }, 20, 0.9, 4, 4, includeReverseComplement: false, clipTips: false);
        Assert.Multiple(() =>
        {
            Assert.That(r.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA" }));
            Assert.That(r.TotalReads, Is.EqualTo(2));
            Assert.That(r.LongestContig, Is.EqualTo(10));
            Assert.That(r.TotalLength, Is.EqualTo(10));
            Assert.That(r.N50, Is.EqualTo(10));
        });
    }

    [Test]
    public void AssembleDeBruijn_BranchingGraph_EmitsUnitigs()
    {
        // k=3: AA and CC carry self-loops (branch nodes), so the contigs are the maximal
        // non-branching paths AAA, AACC, CCC, CCG (Compeau & Pevzner BA3M; Python oracle, single strand).
        var r = AlignmentTools.AssembleDeBruijn(new[] { "AAACCC", "AACCCG" }, 20, 0.9, 3, 3, includeReverseComplement: false, clipTips: false);
        Assert.Multiple(() =>
        {
            Assert.That(r.Contigs, Is.EqualTo(new[] { "AAA", "AACC", "CCC", "CCG" }));
            Assert.That(r.TotalLength, Is.EqualTo(13));
            Assert.That(r.LongestContig, Is.EqualTo(4));
        });
    }

    [Test]
    public void AssembleDeBruijn_ReverseComplement_JoinsStrands()
    {
        var r = AlignmentTools.AssembleDeBruijn(new[] { "GATTACAGGT", "TAAGGACCTG" }, 20, 0.9, 5, 1, includeReverseComplement: true);
        Assert.That(r.Contigs, Is.EqualTo(new[] { "GATTACAGGTCCTTA" }));
    }

    [Test]
    public void AssembleDeBruijn_DefaultIsDoubleStranded_AndFalseIsReachable()
    {
        // A2-4: Velvet records every k-mer with its reverse complement; the default of includeReverseComplement is true.
        var reads = new[] { "GATTACAGGT", "TAAGGACCTG" };
        var byDefault = AlignmentTools.AssembleDeBruijn(reads, 20, 0.9, 5, 1);
        var single = AlignmentTools.AssembleDeBruijn(reads, 20, 0.9, 5, 1, includeReverseComplement: false, clipTips: false);
        Assert.Multiple(() =>
        {
            Assert.That(byDefault.Contigs, Is.EqualTo(new[] { "GATTACAGGTCCTTA" }));
            Assert.That(byDefault.AssembledReads, Is.EqualTo(2));
            Assert.That(single.Contigs, Is.EqualTo(new[] { "GATTACAGGT", "TAAGGACCTG" }));
        });
    }

    [Test]
    public void AssembleDeBruijn_NonAcgtSymbols_AreReadAsA()
    {
        // A2-5: Velvet velvetifySequence + writeNucleotide (N and every other non-ACGT symbol -> A).
        var r = AlignmentTools.AssembleDeBruijn(new[] { "ACGTNGGTCA" }, 20, 0.9, 4, 1, includeReverseComplement: false, clipTips: false);
        Assert.That(r.Contigs, Is.EqualTo(new[] { "ACGTAGGTCA" }));
    }

    // A2-2 (Velvet clipTips_pg + clipTipsHard): the unit-test genome G (50 nt, k = 7) is read at 2x coverage by two
    // reads plus one erroneous read whose last base differs (a one-k-mer tip). wp5_oracle.py (literal Velvet transcription).
    private const string TipGenome = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCA";

    private static readonly string[] TipReads =
    {
        "GCTAAAGACAATTACATAACATACACGTCAGCACGA", "GCTAAAGACAATTACATAACATACACGTCAGCACGA",
        "CATAACATACACGTCAGCACGAAACTTGTTGGCCCA", "CATAACATACACGTCAGCACGAAACTTGTTGGCCCA",
        "GCTAAAGACAATTACATAACATACACGTCC",
    };

    [Test]
    public void AssembleDeBruijn_TipClipping_IsOnByDefault_AndFalseKeepsTheTip()
    {
        var clipped = AlignmentTools.AssembleDeBruijn(TipReads, 20, 0.9, 7, 1);
        var unclipped = AlignmentTools.AssembleDeBruijn(TipReads, 20, 0.9, 7, 1, clipTips: false);
        Assert.Multiple(() =>
        {
            Assert.That(clipped.Contigs, Is.EqualTo(new[] { TipGenome }));
            Assert.That(unclipped.Contigs, Is.EqualTo(new[]
            {
                "CACGTCAGCACGAAACTTGTTGGCCCA", "CACGTCC", "GACGTGTATGTTATGTAATTGTCTTTAGC",
            }));
        });
    }

    [Test]
    public void AssembleDeBruijn_MinTipLength_IsInKmers_AndZeroMeansTwoK()
    {
        // The erroneous tail is a single k-mer: minTipLength 1 keeps it (1 < 1 is false), 2 clips it (Velvet: length < cutoff).
        var keep = AlignmentTools.AssembleDeBruijn(TipReads, 20, 0.9, 7, 1, minTipLength: 1);
        var clip = AlignmentTools.AssembleDeBruijn(TipReads, 20, 0.9, 7, 1, minTipLength: 2);
        Assert.Multiple(() =>
        {
            Assert.That(keep.Contigs.Count, Is.EqualTo(3));
            Assert.That(clip.Contigs, Is.EqualTo(new[] { TipGenome }));
        });
    }

    [Test]
    public void AssembleDeBruijn_CoverageCutoff_UsesIntegerQuotientOfCoverageAndLength()
    {
        // Velvet graphStats.c: getTotalCoverage(node) / getNodeLength(node) is an integer division (7 / 2 = 3), so the
        // unitig ACGT (k = 3, coverage 4 + 3, length 2) goes at a cut-off of 3.5 but stays at 3.
        var reads = new[] { "ACGT", "ACGT", "ACGT", "ACG", "GGCATTAGC", "GGCATTAGC", "GGCATTAGC", "GGCATTAGC", "GGCATTAGC" };
        var removed = AlignmentTools.AssembleDeBruijn(reads, 20, 0.9, 3, 1, includeReverseComplement: false, clipTips: false, coverageCutoff: 3.5);
        var kept = AlignmentTools.AssembleDeBruijn(reads, 20, 0.9, 3, 1, includeReverseComplement: false, clipTips: false, coverageCutoff: 3);
        Assert.Multiple(() =>
        {
            Assert.That(removed.Contigs, Is.EqualTo(new[] { "GCATTAGC", "GGC" }));
            Assert.That(kept.Contigs, Is.EqualTo(new[] { "ACGT", "GCATTAGC", "GGC" }));
        });
    }

    [Test]
    public void AssembleDeBruijn_AutoCoverageCutoff_IsHalfTheLengthWeightedMedianCoverage()
    {
        // Velvet estimated_cov: unitigs A (length 11, coverage 3 per k-mer), B and C (length 2, coverage 1): weighted
        // median 3 -> cut-off 1.5; B and C go. An unweighted median (1 -> cut-off 0.5) would keep them.
        var reads = new[] { "GAATTCAGGTCCAC", "GAATTCAGGTCCAC", "GAATTCAGGTCCAC", "ACATG", "GCTGG" };
        var auto = AlignmentTools.AssembleDeBruijn(reads, 20, 0.9, 4, 1, includeReverseComplement: false, clipTips: false, autoCoverageCutoff: true);
        var off = AlignmentTools.AssembleDeBruijn(reads, 20, 0.9, 4, 1, includeReverseComplement: false, clipTips: false);
        Assert.Multiple(() =>
        {
            Assert.That(auto.Contigs, Is.EqualTo(new[] { "GAATTCAGGTCCAC" }));
            Assert.That(off.Contigs, Is.EqualTo(new[] { "ACATG", "GAATTCAGGTCCAC", "GCTGG" }));
        });
    }

    [Test]
    public void AssembleDeBruijn_InvalidCleaningParameters_Throw()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.AssembleDeBruijn(TipReads, 20, 0.9, 7, 1, minTipLength: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.AssembleDeBruijn(TipReads, 20, 0.9, 7, 1, coverageCutoff: -0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.AssembleDeBruijn(TipReads, 20, 0.9, 7, 1, coverageCutoff: double.NaN));
    }
}
