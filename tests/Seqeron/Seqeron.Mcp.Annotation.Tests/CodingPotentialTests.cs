using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CodingPotentialTests
{
    private static readonly Dictionary<string, double> Coding = new() { ["ATGAAA"] = 8, ["AAACCC"] = 2 };
    private static readonly Dictionary<string, double> Noncoding = new() { ["ATGAAA"] = 2, ["AAACCC"] = 4 };

    [Test]
    public void CodingPotential_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CodingPotential("ATGAAACCC", Coding, Noncoding));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CodingPotential("", Coding, Noncoding));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CodingPotential(null!, Coding, Noncoding));
    }

    [Test]
    public void CodingPotential_TwoInFrameHexamers_ReturnsMeanLogRatio()
    {
        // Mirrors GenomeAnnotator_CalculateCodingPotential_Tests M1.
        // ATGAAA: ln(8/2)=1.3862943611...; AAACCC: ln(2/4)=-0.6931471805...; mean=0.34657359027997264.
        var result = AnnotationTools.CodingPotential("ATGAAACCC", Coding, Noncoding);
        Assert.That(result.Score, Is.EqualTo(0.34657359027997264).Within(1e-10));
    }

    [Test]
    public void CodingPotential_SingleHexamerBothTables_ReturnsNaturalLogRatio()
    {
        // Mirrors M2: ln(4/1) = ln 4 = 1.3862943611198906.
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 4 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 1 };
        var result = AnnotationTools.CodingPotential("ATGAAA", coding, noncoding);
        Assert.That(result.Score, Is.EqualTo(1.3862943611198906).Within(1e-10));
    }

    [Test]
    public void CodingPotential_CodingAndNoncodingOnly_ContributePlusMinusOne()
    {
        // Mirrors M3/M4: coding-only -> +1, noncoding-only -> -1 for a single hexamer.
        var codingOnly = AnnotationTools.CodingPotential(
            "ATGAAA",
            new Dictionary<string, double> { ["ATGAAA"] = 5 },
            new Dictionary<string, double> { ["ATGAAA"] = 0 });
        var noncodingOnly = AnnotationTools.CodingPotential(
            "ATGAAA",
            new Dictionary<string, double> { ["ATGAAA"] = 0 },
            new Dictionary<string, double> { ["ATGAAA"] = 5 });

        Assert.Multiple(() =>
        {
            Assert.That(codingOnly.Score, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(noncodingOnly.Score, Is.EqualTo(-1.0).Within(1e-10));
        });
    }

    [Test]
    public void CodingPotential_SequenceShorterThanHexamer_ReturnsZero()
    {
        // Mirrors M7: len < wordSize -> 0.
        var result = AnnotationTools.CodingPotential(
            "ATGAA",
            new Dictionary<string, double> { ["ATGAA"] = 5 },
            new Dictionary<string, double> { ["ATGAA"] = 5 });
        Assert.That(result.Score, Is.EqualTo(0.0));
    }

    [Test]
    public void CodingPotential_NoScorableHexamer_ReturnsMinusOne()
    {
        // CPAT FrameKmer.kmer_ratio returns -1 when no in-frame hexamer is scored (F11).
        var result = AnnotationTools.CodingPotential("GGGGGGCCC", Coding, Noncoding);
        Assert.That(result.Score, Is.EqualTo(-1.0));
    }

    [Test]
    public void FickettScore_DelegatesToGenomeAnnotator()
    {
        const string seq = "TTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGG";
        var fickett1982 = AnnotationTools.FickettScore(seq);
        var cpat = AnnotationTools.FickettScore(seq, "cpat");
        Assert.Multiple(() =>
        {
            Assert.That(fickett1982.Score, Is.EqualTo(0.8787), "exact Fickett 1982 TESTCODE");
            Assert.That(fickett1982.Call, Is.EqualTo("NoOpinion"));
            Assert.That(cpat.Score, Is.EqualTo(0.8847), "CPAT fickett_value");
        });
        Assert.Throws<ArgumentException>(() => AnnotationTools.FickettScore(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FickettScore(seq, "other"));
    }

    [Test]
    public void FickettScan_MatchesEmbossWindowLayout()
    {
        string seq = string.Concat(Enumerable.Repeat("ATGGCCAAGTTT", 20)); // 240 nt
        var result = AnnotationTools.FickettScan(seq);
        // EMBOSS loop pos + 200 < 240 → pos 0, 3, ..., 39 → 14 windows.
        Assert.That(result.Windows, Has.Count.EqualTo(14));
        Assert.That(result.Windows[0].Start, Is.EqualTo(0));
        Assert.That(result.Windows[13].End, Is.EqualTo(239));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.FickettScan(seq, windowSize: 100));
    }

    [Test]
    public void CpatCodingProbability_DelegatesAndReportsBest()
    {
        // Synthetic table: coding/non-coding = 2 for every hexamer → hexamer score ln 2 for any ORF.
        var coding = new Dictionary<string, double>();
        var noncoding = new Dictionary<string, double>();
        const string b = "ACGT";
        for (int i = 0; i < 4096; i++)
        {
            var c = new char[6];
            for (int p = 0, v = i; p < 6; p++, v /= 4) c[5 - p] = b[v % 4];
            coding[new string(c)] = 2; noncoding[new string(c)] = 1;
        }
        string orf = "ATG" + string.Concat(Enumerable.Repeat("GCC", 40)) + "TAA"; // 126 nt
        string seq = "CC" + orf + "GG";
        var r = AnnotationTools.CpatCodingProbability(seq, coding, noncoding,
            -7.8284881840360567, -5.4831277888113823e-05, 1.1187815405519514e-02, 3.4471695641643225, 6.0682759241106110, 0.364);
        Assert.That(r.Candidates, Has.Count.EqualTo(1));
        var o = r.Candidates[0];
        double fickett = Seqeron.Genomics.Annotation.GenomeAnnotator.CalculateFickettScore(orf,
            Seqeron.Genomics.Annotation.GenomeAnnotator.FickettTableSet.Cpat);
        double eta = -7.8284881840360567 - 5.4831277888113823e-05 * 130 + 1.1187815405519514e-02 * 126
            + 3.4471695641643225 * fickett + 6.0682759241106110 * Math.Log(2);
        Assert.Multiple(() =>
        {
            Assert.That((o.Strand, o.Frame, o.Start, o.End, o.OrfLength, o.MrnaLength), Is.EqualTo(("+", 3, 3, 128, 126, 130)));
            Assert.That(o.HexamerScore, Is.EqualTo(Math.Log(2)).Within(1e-15));
            Assert.That(o.CodingProbability, Is.EqualTo(1 / (1 + Math.Exp(-eta))).Within(1e-15));
            Assert.That(r.BestIndex, Is.EqualTo(0));
        });
        Assert.That(AnnotationTools.CpatCodingProbability("ACGT", coding, noncoding, 0, 0, 0, 0, 0, 0.5).BestIndex, Is.EqualTo(-1));
    }
}
