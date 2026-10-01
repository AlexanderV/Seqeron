using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>kmer_jaccard</c> MCP tool. Expected values: Python set Jaccard, sourmash 4.9.4
/// MinHash(scaled=1) (Jaccard, <c>contained_by</c>) and the Mash 2.3 binary (<c>mash dist -s 100000 [-n]</c>, <c>-k 4 -s 10</c>), not the wrapper output.
/// </summary>
[TestFixture]
public class KmerJaccardTests
{
    [Test]
    public void KmerJaccard_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.KmerJaccard("ACGT", "ACGA", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard("", "ACGT", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard("ACGT", null!, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard("ACGT", "ACGT", 0));
    }

    [Test]
    public void KmerJaccard_Binding_MatchesMashAndSourmash()
    {
        // ACGTTGCAACGGT / ACGTAGCATCGGTA k=2: literal J = 6/13; canonical J = 0.5 (mash 5/10, sourmash 0.5),
        // Mash distance 0.20273255405408222 (mash dist 0.202733); mash -n (acgtOnly) 0.229766.
        var literal = AnalysisTools.KmerJaccard("ACGTTGCAACGGT", "ACGTAGCATCGGTA", 2);
        var canon = AnalysisTools.KmerJaccard("ACGTTGCAACGGT", "ACGTAGCATCGGTA", 2, canonical: true);
        var acgt = AnalysisTools.KmerJaccard("ACGTTGCAACGGT", "ACGTAGCATCGGTA", 2, acgtOnly: true);
        Assert.Multiple(() =>
        {
            Assert.That(literal.Jaccard, Is.EqualTo(6.0 / 13).Within(1e-12));
            Assert.That(canon.Jaccard, Is.EqualTo(0.5));
            Assert.That(canon.MashDistance, Is.EqualTo(0.20273255405408222).Within(1e-12));
            Assert.That(acgt.MashDistance, Is.EqualTo(0.2297661646892201).Within(1e-12));
        });
    }

    [Test]
    public void KmerJaccard_Containment_MatchesSourmashScaled1()
    {
        // sourmash 4.9.4 MinHash(scaled=1).contained_by: S1 in S2 0.35, S2 in S1 3/7 (k=4, canonical);
        // literal ATGTGTG/CATGTG k=3: {ATG,GTG,TGT} ⊂ {ATG,CAT,GTG,TGT} → 1 and 0.75.
        var canon = AnalysisTools.KmerJaccard(S1, S2, 4, canonical: true);
        var literal = AnalysisTools.KmerJaccard("ATGTGTG", "CATGTG", 3);
        Assert.Multiple(() =>
        {
            Assert.That(canon.ContainmentSeq1InSeq2, Is.EqualTo(0.35).Within(1e-15));
            Assert.That(canon.ContainmentSeq2InSeq1, Is.EqualTo(0.42857142857142855).Within(1e-15));
            Assert.That(canon.Jaccard, Is.EqualTo(0.23863636363636365).Within(1e-15));
            Assert.That(canon.SharedHashes, Is.Null);
            Assert.That(canon.PValue, Is.Null);
            Assert.That(literal.ContainmentSeq1InSeq2, Is.EqualTo(1.0));
            Assert.That(literal.ContainmentSeq2InSeq1, Is.EqualTo(0.75));
        });
    }

    [Test]
    public void KmerJaccard_Sketch_MatchesMashDist()
    {
        // mash dist -k 4 -s 10 S1.fa S2.fa → 0.139904  0.0291702  4/10; containment stays exact.
        var r = AnalysisTools.KmerJaccard(S1, S2, 4, canonical: true, sketchSize: 10);
        Assert.Multiple(() =>
        {
            Assert.That(r.SharedHashes, Is.EqualTo(4));
            Assert.That(r.SketchDenominator, Is.EqualTo(10));
            Assert.That(r.Jaccard, Is.EqualTo(0.4));
            Assert.That(r.MashDistance, Is.EqualTo(0.139904).Within(5e-7));
            Assert.That(r.PValue!.Value, Is.EqualTo(0.0291702).Within(5e-8));
            Assert.That(r.ContainmentSeq1InSeq2, Is.EqualTo(0.35).Within(1e-15));
            Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard(S1, S2, 4, sketchSize: -1));
            Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard(S1, S2, 33, sketchSize: 10));
        });
    }

    private const string S1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
    private const string S2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";
}
