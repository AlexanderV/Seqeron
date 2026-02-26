using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
[Category("McpCore")]
public class SuffixTreeCoreConsistencyTests
{
    private static readonly string[] SampleTexts =
    {
        "banana",
        "abracadabra",
        "mississippi",
        "aaaaa",
        "AbCdAbCd"
    };

    [TestCaseSource(nameof(SampleTexts))]
    public void SuffixTreeCoreTools_AreMutuallyConsistent_AndMatchBruteForceOracle(string text)
    {
        foreach (string pattern in BuildPatterns(text))
        {
            bool found = SuffixTreeTools.SuffixTreeContains(text, pattern).Found;
            int count = SuffixTreeTools.SuffixTreeCount(text, pattern).Count;
            int[] positions = SuffixTreeTools.SuffixTreeFindAll(text, pattern).Positions.OrderBy(x => x).ToArray();
            int[] expectedPositions = BruteForcePositions(text, pattern);

            Assert.Multiple(() =>
            {
                Assert.That(found, Is.EqualTo(count > 0), $"contains/count mismatch for pattern='{pattern}'");
                Assert.That(count, Is.EqualTo(positions.Length), $"count/find_all mismatch for pattern='{pattern}'");
                Assert.That(positions, Is.EqualTo(expectedPositions), $"oracle mismatch for pattern='{pattern}'");
            });
        }
    }

    private static IEnumerable<string> BuildPatterns(string text)
    {
        var patterns = new HashSet<string>(StringComparer.Ordinal)
        {
            string.Empty,
            text,
            text + "x",
            text.Substring(0, 1),
            text.Substring(text.Length - 1, 1)
        };

        int midStart = text.Length / 3;
        int midLength = Math.Max(1, text.Length / 3);
        if (midStart + midLength <= text.Length)
            patterns.Add(text.Substring(midStart, midLength));

        // Add several short windows to hit overlaps and repeated symbols.
        for (int i = 0; i < text.Length; i += Math.Max(1, text.Length / 4))
        {
            int len = Math.Min(3, text.Length - i);
            patterns.Add(text.Substring(i, len));
        }

        patterns.Add("xyz");
        return patterns;
    }

    private static int[] BruteForcePositions(string text, string pattern)
    {
        if (pattern.Length == 0)
            return Enumerable.Range(0, text.Length).ToArray();

        var result = new List<int>();
        for (int i = 0; i <= text.Length - pattern.Length; i++)
        {
            if (text.AsSpan(i, pattern.Length).SequenceEqual(pattern.AsSpan()))
                result.Add(i);
        }

        return result.ToArray();
    }
}
