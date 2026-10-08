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
    [Test]
    [Description("Biopython 1.88: counts.normalize(pseudocounts={A:.1,C:.4,G:.2,T:.3}).log_odds({A:.3,C:.2,G:.2,T:.3}) and jaspar.calculate_pseudocounts.")]
    public void CreatePwm_PseudocountOptions_EqualBiopython()
    {
        var seqs = new[] { "TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC" };
        var bg = new[] { 0.3, 0.2, 0.2, 0.3 };
        var perBase = AnalysisTools.CreatePwm(seqs, pseudocounts: new[] { 0.1, 0.4, 0.2, 0.3 }, background: bg);
        Assert.That(perBase.Matrix[0][0], Is.EqualTo(0.3692338096657191).Within(1e-12));
        Assert.That(perBase.Matrix[0][2], Is.EqualTo(-4.584962500721156).Within(1e-12));
        Assert.That(perBase.Matrix[3][0], Is.EqualTo(0.8413022539809418).Within(1e-12));

        var jaspar = AnalysisTools.CreatePwm(seqs, background: bg, jasparPseudocounts: true);
        Assert.That(jaspar.Matrix[0][0], Is.EqualTo(0.39068723333471306).Within(1e-12));
        Assert.That(jaspar.Matrix[2][3], Is.EqualTo(0.8713553378958487).Within(1e-12));

        var jasparUniform = AnalysisTools.CreatePwm(seqs, jasparPseudocounts: true);
        Assert.That(jasparUniform.Matrix[0][0], Is.EqualTo(0.6025166839942642).Within(1e-12));

        var scalarBg = AnalysisTools.CreatePwm(seqs, 0.25, background: bg);
        var expected = global::Seqeron.Genomics.Analysis.MotifFinder.CreatePwm(seqs, 0.25, bg);
        Assert.That(scalarBg.Matrix[1][2], Is.EqualTo(expected.Matrix[1, 2]));

        Assert.Throws<ArgumentException>(() => AnalysisTools.CreatePwm(seqs, pseudocounts: new[] { 1.0, 1.0 }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CreatePwm(seqs, pseudocounts: new[] { 1.0, 1, 1, 1 }, jasparPseudocounts: true));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CreatePwm(seqs, background: new[] { 1.0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.CreatePwm(seqs, pseudocounts: new[] { 1.0, -1, 1, 1 }));
    }
}
