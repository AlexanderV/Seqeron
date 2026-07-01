using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>mask_low_complexity</c> MCP tool.
/// Expected values from SequenceComplexity's own unit test
/// (SequenceComplexityTests.MaskLowComplexity: A*100 window 64 -> fully masked;
/// high-complexity sequence at threshold 10.0 -> unmasked), NOT the wrapper output.
/// </summary>
[TestFixture]
public class MaskLowComplexityTests
{
    [Test]
    public void MaskLowComplexity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.MaskLowComplexity(new string('A', 100), 64, 1.0, 'X'));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MaskLowComplexity("", 64, 1.0, 'N'));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MaskLowComplexity(null!, 64, 1.0, 'N'));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MaskLowComplexity("XYZ", 64, 1.0, 'N'));
    }

    [Test]
    public void MaskLowComplexity_Binding_InvokesSuccessfully()
    {
        // 100 A, window 64, threshold 1.0, mask 'X' -> every position masked.
        var masked = AnalysisTools.MaskLowComplexity(new string('A', 100), 64, 1.0, 'X').Masked;
        Assert.Multiple(() =>
        {
            Assert.That(masked, Has.Length.EqualTo(100));
            Assert.That(masked.Count(c => c == 'X'), Is.EqualTo(100));
        });

        // High-complexity 78bp sequence at threshold 10.0 -> nothing masked.
        const string highComplexity = "ATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCA";
        var preserved = AnalysisTools.MaskLowComplexity(highComplexity, 64, 10.0, 'N').Masked;
        Assert.That(preserved, Does.Not.Contain("N"));
    }
}
