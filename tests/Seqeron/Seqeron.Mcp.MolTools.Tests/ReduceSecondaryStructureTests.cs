using NUnit.Framework;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class ReduceSecondaryStructureTests
{
    private static CodonUsageTableInput EColi() => new(Preset: "EColiK12");

    [Test]
    public void ReduceSecondaryStructure_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.reduce_secondary_structure("ATGATG", EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.reduce_secondary_structure("", EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.reduce_secondary_structure(null!, EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.reduce_secondary_structure("ATGATG", EColi(), 0));
    }

    [Test]
    public void ReduceSecondaryStructure_Binding_InvokesSuccessfully()
    {
        // A sequence shorter than the window is returned unchanged (documented behaviour).
        var shortSeq = MolToolsTools.reduce_secondary_structure("ATGATG", EColi(), window_size: 40);
        Assert.That(shortSeq.OptimizedSequence, Is.EqualTo("ATGATG"));

        // A sequence at/above the window is processed; synonymous swaps preserve the codon count.
        // 15 codons (45 nt) > 40-nt window.
        string coding = string.Concat(System.Linq.Enumerable.Repeat("GCG", 15)); // 15x Ala
        var result = MolToolsTools.reduce_secondary_structure(coding, EColi(), window_size: 40);
        Assert.That(result.OptimizedSequence.Length, Is.EqualTo(coding.Length));
    }
}
