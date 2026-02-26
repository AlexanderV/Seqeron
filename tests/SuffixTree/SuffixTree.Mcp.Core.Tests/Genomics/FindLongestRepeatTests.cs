using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class FindLongestRepeatTests
{
    [Test]
    public void FindLongestRepeat_InvalidArguments_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestRepeat(""));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestRepeat(null!));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.FindLongestRepeat("ATGN")); // invalid DNA symbol
    }

    [Test]
    public void FindLongestRepeat_ReturnsExpectedForRepresentativeInputs()
    {
        var repeated = SuffixTreeTools.FindLongestRepeat("ATGATGATG");
        Assert.That(repeated.Repeat, Is.EqualTo("ATGATG"));
        Assert.That(repeated.Length, Is.EqualTo(6));
        Assert.That(repeated.Positions, Is.EqualTo(new[] { 0, 3 }));

        var none = SuffixTreeTools.FindLongestRepeat("ACGT");
        Assert.That(none.Repeat, Is.Empty);
        Assert.That(none.Length, Is.EqualTo(0));
        Assert.That(none.Positions, Is.Empty);
    }

    [Test]
    public void FindLongestRepeat_MatchesGenomicAnalyzer()
    {
        string[] sequences = { "ATGATGATG", "ACGT", "AACAACTAAC", "TTTTGGGGTTTT" };

        foreach (string sequence in sequences)
        {
            var expected = GenomicAnalyzer.FindLongestRepeat(new DnaSequence(sequence));
            var actual = SuffixTreeTools.FindLongestRepeat(sequence);

            Assert.That(actual.Repeat, Is.EqualTo(expected.Sequence), $"sequence={sequence}: repeat");
            Assert.That(actual.Length, Is.EqualTo(expected.Length), $"sequence={sequence}: length");
            Assert.That(actual.Positions.OrderBy(x => x).ToArray(),
                Is.EqualTo(expected.Positions.OrderBy(x => x).ToArray()),
                $"sequence={sequence}: positions");
        }
    }
}

