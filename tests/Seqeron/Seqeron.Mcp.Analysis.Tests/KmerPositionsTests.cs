using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>kmer_positions</c> MCP tool.
/// Expected values from the Pattern Matching Problem (Rosalind BA1D, overlapping,
/// 0-based, ascending), NOT the wrapper output.
/// </summary>
[TestFixture]
public class KmerPositionsTests
{
    [Test]
    public void KmerPositions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.KmerPositions("AAAA", "AA"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerPositions("", "AA"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerPositions(null!, "AA"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerPositions("AAAA", ""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerPositions("AAAA", null!));
    }

    [Test]
    public void KmerPositions_Binding_InvokesSuccessfully()
    {
        // Overlapping AA in AAAA -> 0,1,2.
        var overlap = AnalysisTools.KmerPositions("AAAA", "AA").Positions;
        Assert.That(overlap, Is.EqualTo(new[] { 0, 1, 2 }));

        // ATG in ATGATG -> 0,3.
        var repeat = AnalysisTools.KmerPositions("ATGATG", "ATG").Positions;
        Assert.That(repeat, Is.EqualTo(new[] { 0, 3 }));

        // Absent k-mer -> empty.
        var absent = AnalysisTools.KmerPositions("ATGATG", "CCC").Positions;
        Assert.That(absent, Is.Empty);
    }
}
