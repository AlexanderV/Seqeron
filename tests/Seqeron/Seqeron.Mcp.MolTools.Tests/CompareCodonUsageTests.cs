using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class CompareCodonUsageTests
{
    [Test]
    public void CompareCodonUsage_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.compare_codon_usage("ATGATG", "ATGATG"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.compare_codon_usage(null!, "ATG"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.compare_codon_usage("ATG", null!));
        // Empty is a defined input (contributes 0 similarity), not an error.
        Assert.DoesNotThrow(() => MolToolsTools.compare_codon_usage("", "ATGATG"));
    }

    [Test]
    public void CompareCodonUsage_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Identical distributions -> 1.0.
            Assert.That(MolToolsTools.compare_codon_usage("ATGATG", "ATGATG").Similarity,
                Is.EqualTo(1.0).Within(1e-9));

            // Disjoint codon sets -> 0.0. (AUG,AUG) vs (UUU,UUU): |1-0|+|0-1| = 2, sim = 1 - 1 = 0.
            Assert.That(MolToolsTools.compare_codon_usage("ATGATG", "TTTTTT").Similarity,
                Is.EqualTo(0.0).Within(1e-9));

            // Half overlap -> 0.5. (AUG:0.5,UUU:0.5) vs (AUG:1): |0.5-1|+|0.5-0| = 1, sim = 1 - 0.5 = 0.5.
            Assert.That(MolToolsTools.compare_codon_usage("ATGTTT", "ATGATG").Similarity,
                Is.EqualTo(0.5).Within(1e-9));

            // Empty input -> 0.
            Assert.That(MolToolsTools.compare_codon_usage("", "ATGATG").Similarity,
                Is.EqualTo(0.0).Within(1e-9));
        });
    }
}
