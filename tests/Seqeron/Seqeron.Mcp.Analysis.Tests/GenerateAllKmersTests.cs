using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>generate_all_kmers</c> MCP tool.
/// Expected values from the k-fold Cartesian product / odometer-order definition
/// (size = alphabet^k), NOT the wrapper output.
/// </summary>
[TestFixture]
public class GenerateAllKmersTests
{
    [Test]
    public void GenerateAllKmers_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.GenerateAllKmers(1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateAllKmers(0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateAllKmers(-2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateAllKmers(2, ""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateAllKmers(2, null!));
    }

    [Test]
    public void GenerateAllKmers_Binding_InvokesSuccessfully()
    {
        // DNA k=1 -> A,C,G,T in odometer order.
        var mono = AnalysisTools.GenerateAllKmers(1).Kmers;
        Assert.That(mono, Is.EqualTo(new[] { "A", "C", "G", "T" }));

        // Binary alphabet "AT" k=2 -> AA,AT,TA,TT.
        var bin = AnalysisTools.GenerateAllKmers(2, "AT").Kmers;
        Assert.That(bin, Is.EqualTo(new[] { "AA", "AT", "TA", "TT" }));

        // Size invariant: |alphabet|^k = 4^3 = 64.
        var trimers = AnalysisTools.GenerateAllKmers(3).Kmers;
        Assert.Multiple(() =>
        {
            Assert.That(trimers, Has.Length.EqualTo(64));
            Assert.That(trimers[0], Is.EqualTo("AAA"));
            Assert.That(trimers[^1], Is.EqualTo("TTT"));
        });
    }

    /// <summary>
    /// Audit round 1 WP4: the |Σ|^k result is materialised as one array, so sizes above the documented cap
    /// (1,048,576 = 4^10) are refused before enumeration instead of overflowing the maximum array length (DNA k ≥ 16).
    /// </summary>
    [Test]
    public void GenerateAllKmers_AboveCap_ThrowsBeforeEnumerating_AtCapSucceeds()
    {
        Assert.That(AnalysisTools.MaxGeneratedKmers, Is.EqualTo(1_048_576));
        Assert.That(AnalysisTools.GenerateAllKmers(10).Kmers, Has.Length.EqualTo(1_048_576)); // 4^10, the cap
        var ex = Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateAllKmers(11));
        Assert.That(ex!.ParamName, Is.EqualTo("k"));
        Assert.That(ex.Message, Does.Contain("1,048,576"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateAllKmers(16));            // > max array length
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateAllKmers(int.MaxValue));  // no overflow, no hang
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateAllKmers(5, "ACDEFGHIKLMNPQRSTVWY")); // 20^5 = 3.2e6
        Assert.That(AnalysisTools.GenerateAllKmers(20, "A").Kmers, Is.EqualTo(new[] { new string('A', 20) }));
    }
}
