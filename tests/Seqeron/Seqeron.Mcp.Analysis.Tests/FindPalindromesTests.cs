using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_palindromes</c> MCP tool.
/// Expected values from the "identical to reverse complement" definition on the
/// EcoRI site GAATTC, NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindPalindromesTests
{
    [Test]
    public void FindPalindromes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindPalindromes("GAATTC", 4, 12));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindPalindromes("", 4, 12));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindPalindromes(null!, 4, 12));
        // minLength must be even and >= 4 (enforced by the algorithm).
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindPalindromes("GAATTC", 3, 12));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindPalindromes("GAATTC", 5, 12));
    }

    [Test]
    public void FindPalindromes_Binding_InvokesSuccessfully()
    {
        // GAATTC -> AATT@1 (len 4) and GAATTC@0 (len 6).
        var eco = AnalysisTools.FindPalindromes("GAATTC", 4, 12).Items;
        Assert.Multiple(() =>
        {
            Assert.That(eco, Has.Length.EqualTo(2));
            var full = eco.Single(p => p.Length == 6);
            Assert.That(full.Sequence, Is.EqualTo("GAATTC"));
            Assert.That(full.Position, Is.EqualTo(0));
            var inner = eco.Single(p => p.Length == 4);
            Assert.That(inner.Sequence, Is.EqualTo("AATT"));
            Assert.That(inner.Position, Is.EqualTo(1));
        });

        // AAAA has no palindrome (revcomp TTTT).
        var none = AnalysisTools.FindPalindromes("AAAA", 4, 12).Items;
        Assert.That(none, Is.Empty);
    }
}
