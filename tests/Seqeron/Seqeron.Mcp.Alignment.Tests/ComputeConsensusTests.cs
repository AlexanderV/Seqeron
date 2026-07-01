using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class ComputeConsensusTests
{
    [Test]
    public void ComputeConsensus_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.ComputeConsensus(new[] { "ACGT", "ACGT" }));
        // Ragged lengths are rejected by the binding when more than one read is supplied.
        Assert.Throws<ArgumentException>(() => AlignmentTools.ComputeConsensus(new[] { "ACGT", "ACG" }));
        // Empty input yields empty consensus.
        Assert.That(AlignmentTools.ComputeConsensus(Array.Empty<string>()).Consensus, Is.EqualTo(""));
    }

    [Test]
    public void ComputeConsensus_Binding_InvokesSuccessfully()
    {
        // Majority vote per column (Biopython dumb_consensus, threshold 0.5):
        // col 4: A,A,A -> A. col 4 of "ACGA" differs but A/G tie broken by majority A.
        // "ACGT","ACGT","ACGA": columns A,C,G,{T,T,A} -> T wins (2/3 >= 0.5).
        var r = AlignmentTools.ComputeConsensus(new[] { "ACGT", "ACGT", "ACGA" });
        Assert.That(r.Consensus, Is.EqualTo("ACGT"));

        // No majority at a column (tie) -> ambiguous 'N'.
        var tie = AlignmentTools.ComputeConsensus(new[] { "A", "C" });
        Assert.That(tie.Consensus, Is.EqualTo("N"));
    }
}
