using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>kmer_jaccard</c> MCP tool. Expected values: Python set Jaccard, sourmash 4.9.4
/// MinHash(scaled=1|3[, track_abundance=True]) (Jaccard, <c>contained_by</c>, <c>angular_similarity</c>) and the Mash 2.3 binary (<c>mash dist -s 100000 [-n]</c>, <c>-k 4 -s 10</c>), not the wrapper output.
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

    [Test]
    public void KmerJaccard_Scaled_MatchesSourmashFracMinHash()
    {
        // sourmash 4.9.4 MinHash(n=0, ksize=4, scaled=3) on S1/S2: len 21 / 17, count_common 9, jaccard 0.3103448275862069,
        // S1.contained_by(S2) 0.4285714285748822, S2.contained_by(S1) 0.52941176525941, max_containment 0.52941176525941;
        // Mash formula -ln(2J/(1+J))/4 = 0.18680360045755526. scaled=1: jaccard 0.23863636363636365 (exact).
        var r = AnalysisTools.KmerJaccard(S1, S2, 4, canonical: true, scaled: 3);
        Assert.Multiple(() =>
        {
            Assert.That(r.SharedHashes, Is.EqualTo(9));
            Assert.That(r.SketchDenominator, Is.EqualTo(29));
            Assert.That(r.Jaccard, Is.EqualTo(0.3103448275862069).Within(1e-15));
            Assert.That(r.ContainmentSeq1InSeq2, Is.EqualTo(0.4285714285748822).Within(1e-15));
            Assert.That(r.ContainmentSeq2InSeq1, Is.EqualTo(0.52941176525941).Within(1e-15));
            Assert.That(r.MaxContainment!.Value, Is.EqualTo(0.52941176525941).Within(1e-15));
            Assert.That(r.MashDistance, Is.EqualTo(0.18680360045755526).Within(1e-15));
            Assert.That(r.PValue, Is.Null);
            Assert.That(AnalysisTools.KmerJaccard(S1, S2, 4, canonical: true, scaled: 1).Jaccard, Is.EqualTo(0.23863636363636365).Within(1e-15));
            Assert.That(AnalysisTools.KmerJaccard(S1, S2, 4).MaxContainment, Is.Null);
            Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard(S1, S2, 4, scaled: -1));
            Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard(S1, S2, 4, canonical: true, sketchSize: 10, scaled: 10));
            Assert.That(r.AngularSimilarity, Is.Null);
        });
    }

    [Test]
    public void KmerJaccard_Scaled_CanonicalRequiredAndU32Range()
    {
        Assert.Multiple(() =>
        {
            // WP10 (B2): canonical=false used to be ignored silently with scaled > 0; sourmash DNA hashing is always canonical
            // (signature.rs SeqToHashes: _hash_murmur(min(kmer, krc))), so it is rejected.
            Assert.That(Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard(S1, S2, 4, scaled: 3))!.ParamName,
                Is.EqualTo("canonical"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard(S1, S2, 4, canonical: false, acgtOnly: true, scaled: 3));
            Assert.That(Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard(S1, S2, 4, trackAbundance: true))!.ParamName,
                Is.EqualTo("trackAbundance"));
            // sourmash ScaledType is u32: 4294967295 accepted (max_hash 4294967297), 4294967296 rejected.
            Assert.DoesNotThrow(() => AnalysisTools.KmerJaccard(S1, S2, 4, canonical: true, scaled: 4294967295L));
            Assert.That(Assert.Throws<ArgumentException>(() => AnalysisTools.KmerJaccard(S1, S2, 4, canonical: true, scaled: 4294967296L))!.ParamName,
                Is.EqualTo("scaled"));
        });
    }

    [Test]
    public void KmerJaccard_Scaled_NonAcgtSkippedAsSourmashForce()
    {
        // sourmash 4.9.4 MinHash(0, 4, scaled=1[, track_abundance=True]).add_sequence(seq, force=True) on N1/N2 (N, R, Y,
        // lower case): len 53 / 44, count_common 18, jaccard 0.22784810126582278, contained_by 0.33962264150943394 /
        // 0.4090909090909091, angular_similarity 0.23696183358819045 (without force=True sourmash raises
        // "invalid DNA character in input k-mer: GTGN"). acgtOnly does not change the result (implied by canonical).
        const string n1 = "AGGTAAGGTGNGTTGAGATctggacTTTTGACGCCTRGAGCCCGCAGTGCTCCTCGAAAAGTAGCNNATGCCTTGGGCTGCT";
        const string n2 = "CAAAGGCCCTACCTTCTTATAGTCCTTYCAACATACAAGTAtagttgGAAGTTCTAAGTTCAGNTTAATC";
        var r = AnalysisTools.KmerJaccard(n1, n2, 4, canonical: true, scaled: 1, trackAbundance: true);
        Assert.Multiple(() =>
        {
            Assert.That(r.SharedHashes, Is.EqualTo(18));
            Assert.That(r.SketchDenominator, Is.EqualTo(53 + 44 - 18));
            Assert.That(r.Jaccard, Is.EqualTo(0.22784810126582278).Within(1e-15));
            Assert.That(r.ContainmentSeq1InSeq2, Is.EqualTo(0.33962264150943394).Within(1e-15));
            Assert.That(r.ContainmentSeq2InSeq1, Is.EqualTo(0.4090909090909091).Within(1e-15));
            Assert.That(r.AngularSimilarity!.Value, Is.EqualTo(0.23696183358819045).Within(1e-15));
            Assert.That(AnalysisTools.KmerJaccard(n1, n2, 4, canonical: true, acgtOnly: true, scaled: 1, trackAbundance: true), Is.EqualTo(r));
        });
    }

    [Test]
    public void KmerJaccard_TrackAbundance_MatchesSourmashAngularSimilarity()
    {
        // sourmash 4.9.4 MinHash(0, 4, scaled=S, track_abundance=True) on S1/S2: angular_similarity 0.2363801370444173 (S = 1),
        // 0.3123095603640216 (S = 3); Jaccard / containments unchanged by abundance tracking.
        var r1 = AnalysisTools.KmerJaccard(S1, S2, 4, canonical: true, scaled: 1, trackAbundance: true);
        var r3 = AnalysisTools.KmerJaccard(S1, S2, 4, canonical: true, scaled: 3, trackAbundance: true);
        Assert.Multiple(() =>
        {
            Assert.That(r1.AngularSimilarity!.Value, Is.EqualTo(0.2363801370444173).Within(1e-15));
            Assert.That(r1.Jaccard, Is.EqualTo(0.23863636363636365).Within(1e-15));
            Assert.That(r3.AngularSimilarity!.Value, Is.EqualTo(0.3123095603640216).Within(1e-15));
            Assert.That(r3.Jaccard, Is.EqualTo(0.3103448275862069).Within(1e-15));
            Assert.That(r3.ContainmentSeq1InSeq2, Is.EqualTo(0.4285714285748822).Within(1e-15));
        });
    }

    private const string S1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
    private const string S2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";
}
