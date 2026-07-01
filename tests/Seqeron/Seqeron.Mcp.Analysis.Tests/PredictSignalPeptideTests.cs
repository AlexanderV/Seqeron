using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>predict_signal_peptide</c> MCP tool.
/// Expected values from ProteinMotifFinder's own unit test
/// (ProteinMotifFinder_PredictSignalPeptide_Tests, ACH2_DROME vs EMBOSS sigcleave:
/// cleavage 42, score 13.7390400704164, window "LLVLLLLCETVQANP"), NOT the wrapper output.
/// </summary>
[TestFixture]
public class PredictSignalPeptideTests
{
    private const string Ach2Drome =
        "MAPGCCTTRPRPIALLAHIWRHCKPLCLLLVLLLLCETVQANPDAKRLYDDLLSNYNRLI" +
        "RPVSNNTDTVLVKLGLRLSQLIDLNLKDQILTTNVWLEHEWQDHKFKWDPSEYGGVTELY" +
        "VPSEHIWLPDIVLYNNADGEYVVTTMTKAILHYTGKVVWTPPAIFKSSCEIDVRYFPFDQ" +
        "QTCFMKFGSWTYDGDQIDLKHISQKNDKDNKVEIGIDLREYYPSVEWDILGVPAERHEKY" +
        "YPCCAEPYPDIFFNITLRRKTLFYTVNLIIPCVGISYLSVLVFYLPADSGEKIALCISIL" +
        "LSQTMFFLLISEIIPSTSLALPLLGKYLLFTMLLVGLSVVITIIILNIHYRKPSTHKMRP" +
        "WIRSFFIKRLPKLLLMRVPKDLLRDLAANKINYGLKFSKTKFGQALMDEMQMNSGGSSPD" +
        "SLRRMQGRVGAGGCNGMHVTTATNRFSGLVGALGGGLSTLSGYNGLPSVLSGLDDSLSDV" +
        "AARKKYPFELEKAIHNVMFIQHHMQRQDEFNAEDQDWGFVAMVMDRLFLWLFMIASLVGT" +
        "FVILGEAPSLYDDTKAIDVQLSDVAKQIYNLTEKKN";

    [Test]
    public void PredictSignalPeptide_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PredictSignalPeptide(Ach2Drome));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictSignalPeptide(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictSignalPeptide(null!));
    }

    [Test]
    public void PredictSignalPeptide_Binding_InvokesSuccessfully()
    {
        var sp = AnalysisTools.PredictSignalPeptide(Ach2Drome);
        Assert.Multiple(() =>
        {
            Assert.That(sp.Found, Is.True);
            Assert.That(sp.CleavagePosition, Is.EqualTo(42));
            Assert.That(sp.Score, Is.EqualTo(13.7390400704164).Within(1e-3));
            Assert.That(sp.WindowSequence, Is.EqualTo("LLVLLLLCETVQANP"));
            Assert.That(sp.IsLikelySignalPeptide, Is.True);
        });
    }
}
