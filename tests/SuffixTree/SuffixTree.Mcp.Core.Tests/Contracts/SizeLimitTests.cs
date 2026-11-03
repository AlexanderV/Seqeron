using System.Reflection;
using Microsoft.Extensions.Configuration;
using ModelContextProtocol.Server;
using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

/// <summary>
/// ST-S0 WP 0.6: every MCP tool that builds a suffix tree from its string arguments must reject an oversized call
/// with a clear <see cref="ArgumentException"/> BEFORE the tree is built, and the limit must be configurable.
/// The tests lower the limit to a few dozen characters, so nothing large is ever built.
/// </summary>
[TestFixture]
[Category("McpCore")]
[NonParallelizable]
public class SizeLimitTests
{
    private const int TestLimit = 50;
    private static readonly string Big = new('A', TestLimit + 1);     // alone one character above the limit
    private static readonly string Half = new('A', TestLimit / 2 + 1); // two of these exceed the limit together

    [SetUp]
    public void SetUp() => SuffixTreeMcpLimits.SetMaxTextChars(TestLimit);

    [TearDown]
    public void TearDown() => SuffixTreeMcpLimits.SetMaxTextChars(SuffixTreeMcpLimits.DefaultMaxTextChars);

    /// <summary>Every tool that builds a tree, called with arguments whose total size exceeds the limit.</summary>
    private static IEnumerable<TestCaseData> OversizedCalls()
    {
        static TestCaseData Case(string tool, Action call) => new TestCaseData(call).SetName(tool);

        yield return Case("suffix_tree_contains", () => SuffixTreeCoreTools.SuffixTreeContains(Big, "A"));
        yield return Case("suffix_tree_contains/pattern", () => SuffixTreeCoreTools.SuffixTreeContains("ACGT", Big));
        yield return Case("suffix_tree_count", () => SuffixTreeCoreTools.SuffixTreeCount(Big, "A"));
        yield return Case("suffix_tree_find_all", () => SuffixTreeCoreTools.SuffixTreeFindAll(Big, "A"));
        yield return Case("suffix_tree_lrs", () => SuffixTreeCoreTools.SuffixTreeLrs(Big));
        yield return Case("suffix_tree_lcs", () => SuffixTreeCoreTools.SuffixTreeLcs(Half, Half));
        yield return Case("suffix_tree_stats", () => SuffixTreeCoreTools.SuffixTreeStats(Big));
        yield return Case("suffix_tree_all_lrs", () => SuffixTreeCoreTools.SuffixTreeAllLrs(Big));
        yield return Case("suffix_tree_find_mems", () => SuffixTreeCoreTools.SuffixTreeFindMems(Half, Half, 3));
        yield return Case("suffix_tree_find_mums", () => SuffixTreeCoreTools.SuffixTreeFindMums(Half, Half, 3));
        yield return Case("suffix_tree_maximal_repeats", () => SuffixTreeCoreTools.SuffixTreeMaximalRepeats(Big, 3));
        yield return Case("suffix_tree_all_lcs", () => SuffixTreeCoreTools.SuffixTreeAllLcs(Half, Half));
        yield return Case("suffix_tree_k_common_substrings", () => SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(new[] { Half, Half }));
        yield return Case("find_longest_repeat", () => SuffixTreeGenomicsTools.FindLongestRepeat(Big));
        yield return Case("find_longest_common_region", () => SuffixTreeGenomicsTools.FindLongestCommonRegion(Half, Half));
    }

    [TestCaseSource(nameof(OversizedCalls))]
    public void OversizedCall_IsRejectedWithAClearArgumentException(Action call)
    {
        var ex = Assert.Throws<ArgumentException>(() => call());

        Assert.Multiple(() =>
        {
            Assert.That(ex!.Message, Does.StartWith("Input too large"));
            Assert.That(ex.Message, Does.Contain("limit of 50 characters"));
            Assert.That(ex.Message, Does.Contain("SuffixTree:MaxTextChars"));
            Assert.That(ex.Message, Does.Contain("SuffixTree__MaxTextChars"));
            Assert.That(ex.ParamName, Is.Not.Null.And.Not.Empty);
        });
    }

    [Test]
    public void EveryCoreTool_IsCoveredByTheOversizedCallTable()
    {
        // A new suffix_tree_* tool must be added to OversizedCalls (and guarded), otherwise this test fails.
        var coreTools = typeof(SuffixTreeCoreTools)
            .GetMethods(BindingFlags.Public | BindingFlags.Static)
            .Select(m => m.GetCustomAttribute<McpServerToolAttribute>()?.Name)
            .Where(n => n is not null)
            .ToHashSet();
        var covered = OversizedCalls().Select(c => c.TestName!.Split('/')[0]).ToHashSet();

        Assert.That(covered, Is.SupersetOf(coreTools!));
    }

    [Test]
    public void InputsAtTheLimit_AreAccepted_AndTheLimitCountsTheSumOfAllTexts()
    {
        string atLimit = new('A', TestLimit);

        Assert.DoesNotThrow(() => SuffixTreeCoreTools.SuffixTreeStats(atLimit));
        Assert.DoesNotThrow(() => SuffixTreeCoreTools.SuffixTreeLcs(new string('A', 25), new string('A', 25)));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeLcs(new string('A', 25), new string('A', 26)));
        // k texts: every text counts its characters plus one separator.
        Assert.DoesNotThrow(() => SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(new[] { new string('A', 24), new string('A', 24) }));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeKCommonSubstrings(new[] { new string('A', 25), new string('A', 25) }));
    }

    [Test]
    public void ArgumentValidation_StillRunsBeforeTheSizeCheck()
    {
        // Existing contract: empty/null inputs are reported as such, not as "too large".
        var ex = Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeStats(""));
        Assert.That(ex!.Message, Does.Not.Contain("too large"));
    }

    [Test]
    public void NonTreeTools_AreNotAffected()
    {
        // Hamming distance and k-mer similarity do not build a suffix tree; they keep working above the limit.
        string a = new('A', 200), b = new('C', 200);

        Assert.That(SuffixTreeGenomicsTools.HammingDistance(a, b).Distance, Is.EqualTo(200));
        Assert.DoesNotThrow(() => SuffixTreeGenomicsTools.CalculateSimilarity(a, b));
    }

    // ───────────────────────── default and derivation ─────────────────────────

    [Test]
    public void DefaultLimit_IsFiveMillionCharacters_DerivedFromTheMeasuredMemory()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SuffixTreeMcpLimits.DefaultMaxTextChars, Is.EqualTo(5_000_000));
            // Design constant >= the worst measured value (193 B/char, repetitive 10 Mbp DNA; 169 for random DNA).
            Assert.That(SuffixTreeMcpLimits.DesignBytesPerChar, Is.GreaterThanOrEqualTo(193));
            // One call stays within 1 GiB at the design constant.
            Assert.That((long)SuffixTreeMcpLimits.DefaultMaxTextChars * SuffixTreeMcpLimits.DesignBytesPerChar, Is.LessThanOrEqualTo(1L << 30));
            // ... and the default is the round number just below the 1 GiB / 200 B limit of 5,368,709 characters.
            Assert.That(SuffixTreeMcpLimits.DefaultMaxTextChars, Is.GreaterThan((1L << 30) / SuffixTreeMcpLimits.DesignBytesPerChar / 2));
        });
    }

    [Test]
    public void DefaultLimit_RejectsAFiveMillionAndOneCharacterTextWithoutBuildingATree()
    {
        SuffixTreeMcpLimits.SetMaxTextChars(SuffixTreeMcpLimits.DefaultMaxTextChars);
        string tooBig = new('A', SuffixTreeMcpLimits.DefaultMaxTextChars + 1); // 10 MB string; a tree would need ~1 GB and ~7 s

        var watch = System.Diagnostics.Stopwatch.StartNew();
        var ex = Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeStats(tooBig));
        watch.Stop();

        Assert.That(ex!.Message, Does.Contain("5,000,001").And.Contain("5,000,000"));
        Assert.That(watch.ElapsedMilliseconds, Is.LessThan(1000), "rejection must not build the tree");
    }

    // ───────────────────────── configuration ─────────────────────────

    private static IConfiguration Config(params (string Key, string? Value)[] values)
        => new ConfigurationBuilder()
            .AddInMemoryCollection(values.ToDictionary(v => v.Key, v => v.Value))
            .Build();

    [Test]
    public void Configure_NoValue_KeepsTheDefault()
    {
        SuffixTreeMcpLimits.Configure(Config());

        Assert.That(SuffixTreeMcpLimits.MaxTextChars, Is.EqualTo(SuffixTreeMcpLimits.DefaultMaxTextChars));
    }

    [TestCase("1")]
    [TestCase("123")]
    [TestCase("20000000")]
    [TestCase("2147483647")]
    public void Configure_PositiveInteger_SetsTheLimit(string value)
    {
        SuffixTreeMcpLimits.Configure(Config((SuffixTreeMcpLimits.ConfigurationKey, value)));

        Assert.That(SuffixTreeMcpLimits.MaxTextChars, Is.EqualTo(int.Parse(value)));
    }

    [TestCase("abc")]
    [TestCase("0")]
    [TestCase("-5")]
    [TestCase("1.5")]
    [TestCase("1e6")]
    [TestCase("2147483648")]
    [TestCase("1 000")]
    public void Configure_InvalidValue_FailsFastWithTheKeyInTheMessage(string value)
    {
        var ex = Assert.Throws<InvalidOperationException>(() => SuffixTreeMcpLimits.Configure(Config((SuffixTreeMcpLimits.ConfigurationKey, value))));

        Assert.That(ex!.Message, Does.Contain("SuffixTree:MaxTextChars").And.Contain(value));
    }

    [Test]
    public void ConfiguredLimit_IsUsedByTheTools()
    {
        SuffixTreeMcpLimits.Configure(Config((SuffixTreeMcpLimits.ConfigurationKey, "10")));

        Assert.DoesNotThrow(() => SuffixTreeCoreTools.SuffixTreeStats("ACGTACGTAC"));
        var ex = Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeStats("ACGTACGTACG"));
        Assert.That(ex!.Message, Does.Contain("limit of 10 characters"));
    }
}
