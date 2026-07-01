using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>entropy_profile</c> MCP tool.
/// Expected values computed by hand from SequenceStatistics.CalculateEntropyProfile:
/// per window H = -sum p*log2(p) over its letter frequencies. NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class EntropyProfileTests
{
    [Test]
    public void EntropyProfile_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.EntropyProfile("ACGT", 2, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.EntropyProfile("", 2, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.EntropyProfile(null!, 2, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.EntropyProfile("ACGT", 0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.EntropyProfile("ACGT", 2, 0));
    }

    [Test]
    public void EntropyProfile_Binding_InvokesSuccessfully()
    {
        // "AAAA", window 2, step 1: three "AA" windows, each H = 0.
        var homo = AnalysisTools.EntropyProfile("AAAA", 2, 1).Values;
        Assert.That(homo, Is.EqualTo(new[] { 0.0, 0.0, 0.0 }));

        // "ATGC", window 4: one window, all four bases equal -> H = 2.0 bits.
        var uniform = AnalysisTools.EntropyProfile("ATGC", 4, 1).Values;
        Assert.That(uniform, Has.Length.EqualTo(1));
        Assert.That(uniform[0], Is.EqualTo(2.0).Within(1e-12));

        // "AT", window 2: p=0.5 each -> H = 1.0 bit.
        var two = AnalysisTools.EntropyProfile("AT", 2, 1).Values;
        Assert.That(two, Has.Length.EqualTo(1));
        Assert.That(two[0], Is.EqualTo(1.0).Within(1e-12));
    }
}
