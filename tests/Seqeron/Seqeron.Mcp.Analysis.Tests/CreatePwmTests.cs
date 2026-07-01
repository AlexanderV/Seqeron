using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>create_pwm</c> MCP tool.
/// Expected values computed by hand from the log-odds formula in MotifFinder.CreatePwm:
/// freq = (count + pc) / (N + 4*pc), score = log2(freq / 0.25). NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class CreatePwmTests
{
    [Test]
    public void CreatePwm_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.CreatePwm(new[] { "ACGT", "ACGT" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CreatePwm(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CreatePwm(Array.Empty<string>()));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.CreatePwm(new[] { "ACGT" }, -0.5));
        // Unequal lengths / invalid bases propagate from the algorithm.
        Assert.Throws<ArgumentException>(() => AnalysisTools.CreatePwm(new[] { "ACGT", "AC" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CreatePwm(new[] { "ACXT" }));
    }

    [Test]
    public void CreatePwm_Binding_InvokesSuccessfully()
    {
        // Two identical "ACGT" sequences, pseudocount 0.25.
        // Per position: present base count 2 -> freq = 2.25/3 = 0.75 -> log2(0.75/0.25) = log2(3).
        //               absent base count 0 -> freq = 0.25/3 -> log2((1/12)/0.25) = log2(1/3) = -log2(3).
        double log2Three = Math.Log2(3.0);
        var r = AnalysisTools.CreatePwm(new[] { "ACGT", "ACGT" }, 0.25);

        Assert.Multiple(() =>
        {
            Assert.That(r.Length, Is.EqualTo(4));
            Assert.That(r.Consensus, Is.EqualTo("ACGT"));
            // Rows are A,C,G,T. Diagonal (present base) = +log2(3); off-diagonal = -log2(3).
            Assert.That(r.Matrix[0][0], Is.EqualTo(log2Three).Within(1e-9));   // A at pos0
            Assert.That(r.Matrix[1][0], Is.EqualTo(-log2Three).Within(1e-9));  // C at pos0 (absent)
            Assert.That(r.Matrix[1][1], Is.EqualTo(log2Three).Within(1e-9));   // C at pos1
            Assert.That(r.MaxScore, Is.EqualTo(4 * log2Three).Within(1e-9));
            Assert.That(r.MinScore, Is.EqualTo(-4 * log2Three).Within(1e-9));
        });
    }
}
