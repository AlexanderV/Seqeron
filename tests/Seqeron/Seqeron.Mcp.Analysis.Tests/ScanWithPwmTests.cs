using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>scan_with_pwm</c> MCP tool.
/// Expected scores follow the documented PWM scoring rule
/// (score = sum over positions of Matrix[baseIndex, j], rows A,C,G,T), applied to a
/// hand-built matrix that scores "ATGC" as +4 — NOT the wrapper output.
/// </summary>
[TestFixture]
public class ScanWithPwmTests
{
    // Rows A, C, G, T; columns = positions. +1 for the ATGC-matching base, -1 otherwise.
    // A at pos0, T at pos1, G at pos2, C at pos3.
    private static PwmInput AtgcPwm() => new(
        new[]
        {
            new[] {  1.0, -1.0, -1.0, -1.0 }, // A
            new[] { -1.0, -1.0, -1.0,  1.0 }, // C
            new[] { -1.0, -1.0,  1.0, -1.0 }, // G
            new[] { -1.0,  1.0, -1.0, -1.0 }  // T
        },
        4);

    [Test]
    public void ScanWithPwm_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.ScanWithPwm("ATGC", AtgcPwm()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ScanWithPwm("", AtgcPwm()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ScanWithPwm(null!, AtgcPwm()));
        Assert.Throws<ArgumentNullException>(() => AnalysisTools.ScanWithPwm("ATGC", null!));
        // Wrong row count -> ArgumentException from the jagged->matrix conversion.
        Assert.Throws<ArgumentException>(() =>
            AnalysisTools.ScanWithPwm("ATGC", new PwmInput(new[] { new[] { 1.0, 1.0, 1.0, 1.0 } }, 4)));
    }

    [Test]
    public void ScanWithPwm_Binding_InvokesSuccessfully()
    {
        // "ATGCATGC" with threshold 0: ATGC scores +4 at positions 0 and 4.
        var matches = AnalysisTools.ScanWithPwm("ATGCATGC", AtgcPwm(), threshold: 0.0).Items;
        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Length.EqualTo(2));
            Assert.That(matches.Select(m => m.Position), Is.EqualTo(new[] { 0, 4 }));
            Assert.That(matches.Select(m => m.MatchedSequence), Is.All.EqualTo("ATGC"));
            Assert.That(matches.Select(m => m.Score), Is.All.EqualTo(4.0).Within(1e-10));
        });

        // "AAAA" scores +1-1-1-1 = -2 < 0 -> no match at threshold 0.
        var none = AnalysisTools.ScanWithPwm("AAAA", AtgcPwm(), threshold: 0.0).Items;
        Assert.That(none, Is.Empty);
    }
}
