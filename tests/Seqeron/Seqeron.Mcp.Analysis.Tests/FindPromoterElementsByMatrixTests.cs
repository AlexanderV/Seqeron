using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// find_promoter_elements_by_matrix — delegation to <c>MotifFinder.FindPromoterElementsByMatrix</c>.
/// Expected values: Biopython pssm.search at threshold_fpr(1e-3) with the JASPAR Bucher matrices,
/// locked in MotifFinder_PwmStrandsAndThresholds_Tests.
/// </summary>
[TestFixture]
public class FindPromoterElementsByMatrixTests
{
    private const string Promoter =
        "GGGGCTATAAAAGGGGGTGGGGGCGCGTTCGTCCTCACTCTCTTCCGCATCGCTGTCTGCGAGGGCCAGCCAATCAGCGCCCCGCCCATTGGCTGGGCGGAGCC";

    [Test]
    public void FindPromoterElementsByMatrix_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindPromoterElementsByMatrix(Promoter));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindPromoterElementsByMatrix(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindPromoterElementsByMatrix("ACGU"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindPromoterElementsByMatrix(Promoter, 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindPromoterElementsByMatrix(Promoter, double.NaN));
    }

    [Test]
    public void FindPromoterElementsByMatrix_Fpr1e3_EqualsBiopythonSearch()
    {
        var items = AnalysisTools.FindPromoterElementsByMatrix(Promoter, 1e-3).Items;
        int n = Promoter.Length;
        var expected = new (string Id, int Pos, string Strand, double Score)[]
        {
            ("POL012.1", 4, "+", 14.674918174743652),
            ("POL002.1", 34, "+", 6.793403148651123),
            ("POL004.1", 64, "+", 12.077281951904297),
            ("POL004.1", n - 19, "-", 12.315014839172363),
            ("POL003.1", 11, "+", 13.369436264038086),
            ("POL003.1", 17, "+", 8.659268379211426),
            ("POL003.1", n - 28, "-", 13.07020378112793),
        };
        Assert.That(items.Select(h => (h.MatrixId, h.Position, h.Strand)), Is.EqualTo(expected.Select(e => (e.Id, e.Pos, e.Strand))));
        Assert.That(items.Select(h => h.Score), Is.EqualTo(expected.Select(e => e.Score)).Within(1e-5));
        Assert.That(items.First(h => h.MatrixId == "POL012.1").Threshold, Is.EqualTo(7.042753822501211).Within(1e-9));
    }

    [Test]
    public void FindPromoterElementsByMatrix_SingleStrand_EqualsLibrary()
    {
        var items = AnalysisTools.FindPromoterElementsByMatrix(Promoter, 1e-3, bothStrands: false).Items;
        var lib = global::Seqeron.Genomics.Analysis.MotifFinder.FindPromoterElementsByMatrix(
            new global::Seqeron.Genomics.Core.DnaSequence(Promoter), 1e-3, bothStrands: false).ToArray();
        Assert.That(items, Has.Length.EqualTo(5));
        Assert.That(items.Select(h => (h.Name, h.MatrixId, h.Position, h.Strand, h.Sequence, h.Score, h.Threshold)),
            Is.EqualTo(lib.Select(h => (h.Name, h.MatrixId, h.Position, h.Strand.ToString(), h.Sequence, h.Score, h.Threshold))));
    }
}
