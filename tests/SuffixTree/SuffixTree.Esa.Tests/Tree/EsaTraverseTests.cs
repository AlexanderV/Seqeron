using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Tree;

/// <summary>
/// <c>Traverse</c> of the ESA (work package 1b.1): the event stream (<c>VisitNode(start, end, leafCount, childCount, depth)</c>, <c>EnterBranch</c>
/// keys, <c>ExitBranch</c>) is identical to the Ukkonen tree's and to the persistent v6 tree's, hence so is
/// <see cref="SuffixTreeSerializer.CalculateLogicalHash"/> (hypothesis H-OCC, ADR-001 section 4). The independent L0 model is
/// <see cref="NaiveEsaModel"/>.
/// </summary>
[TestFixture]
public class EsaTraverseTests
{
    private static List<string> Stream(Action<ISuffixTreeVisitor> traverse)
    {
        var recorder = new TraverseRecorder();
        traverse(recorder);
        return recorder.Events;
    }

    private static void AssertSameStream(List<string> expected, List<string> actual, string what, string context)
    {
        if (expected.SequenceEqual(actual)) return;
        int i = 0;
        while (i < expected.Count && i < actual.Count && expected[i] == actual[i]) i++;
        Assert.Fail($"{context}: {what}: first difference at event {i}: expected '{(i < expected.Count ? expected[i] : "<end>")}' vs actual '{(i < actual.Count ? actual[i] : "<end>")}' ({expected.Count} vs {actual.Count} events)");
    }

    private static void AssertEqualToUkkonen(string text, ChildLayoutKind layout, string context)
    {
        EsaIndex index = EsaIndex.Build(text, layout);
        List<string> actual = Stream(index.Traverse);
        List<string> expected = Stream(UkkonenTree.Build(text).Traverse);
        AssertSameStream(expected, actual, $"ESA({layout}) vs Ukkonen", context);
    }

    private static void AssertEqualToPersistent(string text, string context)
    {
        using PersistentSuffixTree v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
        EsaIndex index = EsaIndex.Build(text);
        AssertSameStream(Stream(v6.Traverse), Stream(index.Traverse), "ESA vs persistent v6", context);
    }

    private static void AssertHashesEqual(string text, string context)
    {
        byte[] ukkonen = SuffixTreeSerializer.CalculateLogicalHash(UkkonenTree.Build(text));
        byte[] esa = SuffixTreeSerializer.CalculateLogicalHash(new EsaTreeShim(EsaIndex.Build(text)));
        using PersistentSuffixTree v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
        byte[] persistent = SuffixTreeSerializer.CalculateLogicalHash(v6);
        Assert.That(Convert.ToHexString(esa), Is.EqualTo(Convert.ToHexString(ukkonen)), $"{context}: ESA vs Ukkonen hash");
        Assert.That(Convert.ToHexString(persistent), Is.EqualTo(Convert.ToHexString(ukkonen)), $"{context}: v6 vs Ukkonen hash");
    }

    [Test]
    public void Banana_StreamIsTheUkkonenStream_ForEveryLayout()
    {
        foreach (ChildLayoutKind layout in Enum.GetValues<ChildLayoutKind>())
            foreach (string s in new[] { "banana", "mississippi", "abracadabra", "a", "aa", "ab", "aba" })
                AssertEqualToUkkonen(s, layout, $"'{s}'");
    }

    [Test]
    public void EmptyText_StreamIsTheSingleRootEvent()
    {
        List<string> stream = Stream(EsaIndex.Build("").Traverse);
        Assert.That(stream, Is.EqualTo(new[] { "N 0 0 0 0 0" }));
        AssertSameStream(Stream(UkkonenTree.Build("").Traverse), stream, "ESA vs Ukkonen", "''");
        AssertEqualToPersistent("", "''");
        AssertHashesEqual("", "''");
    }

    [Test]
    public void Banana_StreamBeginsWithTheRootAndItsTerminatorLeaf()
    {
        List<string> stream = Stream(EsaIndex.Build("banana").Traverse);
        // root: (0, 0, 7 leaves, 4 children, depth 0); first child key -1 = the empty-suffix leaf (start = n = 6, boundless end)
        Assert.That(stream.Take(4), Is.EqualTo(new[] { "N 0 0 7 4 0", "E -1", "N 6 -1 1 0 0", "X" }));
    }

    [Test]
    public void AllSmallStrings_StreamEqualsUkkonen()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 10), (3, 6), (4, 4) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
                AssertEqualToUkkonen(s, ChildLayoutKind.Packed, $"'{s}'");
    }

    [Test]
    public void AllSmallStrings_AllLayoutsGiveTheSameStream()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 8), (3, 5) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
            {
                List<string> reference = Stream(EsaIndex.Build(s, ChildLayoutKind.Reference3).Traverse);
                foreach (ChildLayoutKind layout in new[] { ChildLayoutKind.Packed, ChildLayoutKind.Tagged })
                    AssertSameStream(reference, Stream(EsaIndex.Build(s, layout).Traverse), $"{layout} vs Reference3", $"'{s}'");
            }
    }

    [Test]
    public void SmallStrings_StreamEqualsThePersistentV6Tree()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 8), (3, 5), (4, 4) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
                AssertEqualToPersistent(s, $"'{s}'");
    }

    [Test]
    public void SmallStrings_StreamEqualsTheNaiveModel()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 9), (3, 5), (4, 4) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
                AssertSameStream(Stream(v => NaiveEsaModel.EmitTraverse(s, v)), Stream(EsaIndex.Build(s).Traverse), "ESA vs naive model", $"'{s}'");
    }

    [Test]
    public void AdversarialFamilies_StreamEqualsUkkonenAndV6()
    {
        var rng = new SplitMix64(0x1B06);
        for (int i = 0; i < 20; i++)
        {
            int n = rng.NextInt(250);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 4 == 0))
            {
                AssertEqualToUkkonen(text, ChildLayoutKind.Packed, $"{name} n {n}");
                AssertEqualToPersistent(text, $"{name} n {n}");
            }
        }
    }

    [Test]
    public void AdversarialFamilies_AtFixedSizes()
    {
        foreach (int n in new[] { 0, 1, 2, 3, 4, 5, 64, 300 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 11))
                AssertEqualToUkkonen(text, ChildLayoutKind.Packed, $"{name} n {n}");
    }

    [Test]
    public void Fibonacci_ThueMorse_AllSame_RunsOfN_StreamEqualsUkkonen()
    {
        foreach (int n in new[] { 1000, 3000 })
        {
            AssertEqualToUkkonen(StringFamilies.Fibonacci(n), ChildLayoutKind.Packed, $"fibonacci {n}");
            AssertEqualToUkkonen(StringFamilies.ThueMorse(n), ChildLayoutKind.Packed, $"thue-morse {n}");
            AssertEqualToUkkonen(StringFamilies.AllSame(n), ChildLayoutKind.Packed, $"a^{n}");
            AssertEqualToUkkonen(StringFamilies.AbRepeat(n), ChildLayoutKind.Packed, $"abab {n}");
            AssertEqualToUkkonen(new string('N', n / 2) + StringFamilies.RandomDna(n / 4, 3) + new string('N', n / 4), ChildLayoutKind.Packed, $"N runs {n}");
        }
    }

    [Test]
    public void LogicalHash_IsEqualAcrossEsaUkkonenAndV6()
    {
        foreach (string s in new[] { "", "a", "banana", "mississippi", "abababab", "aaaaaaa", "abracadabraabracadabra", "\uD800a\uD800\0￿" })
            AssertHashesEqual(s, $"'{s}'");
        foreach (string s in ExhaustiveStrings.All(2, 7).Concat(ExhaustiveStrings.All(3, 4)))
            AssertHashesEqual(s, $"'{s}'");
        var rng = new SplitMix64(0x1B07);
        for (int i = 0; i < 10; i++)
        {
            int n = 1 + rng.NextInt(200);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 6 == 0))
                AssertHashesEqual(text, $"{name} n {n}");
        }
    }

    [Test]
    public void Traverse_NullVisitorThrows()
    {
        Assert.That(() => EsaIndex.Build("abc").Traverse(null!), Throws.TypeOf<ArgumentNullException>());
    }

    [Test]
    public void Traverse_DeepChainDoesNotOverflowTheStack()
    {
        // a^n has a chain of n internal nodes: the walk must not recurse
        string text = StringFamilies.AllSame(50_000);
        var counter = new CountingVisitor();
        EsaIndex index = EsaIndex.Build(text);
        index.Traverse(counter);
        // root + (n - 1) internal nodes + (n + 1) leaves
        Assert.That(counter.Nodes, Is.EqualTo(2 * 50_000 + 1));
        Assert.That(counter.Nodes, Is.EqualTo(index.Statistics.NodeCount));
        Assert.That(counter.Enters, Is.EqualTo(counter.Nodes - 1));
        Assert.That(counter.Exits, Is.EqualTo(counter.Enters));
    }

    private sealed class CountingVisitor : ISuffixTreeVisitor
    {
        internal int Nodes;
        internal int Enters;
        internal int Exits;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Nodes++;

        public void EnterBranch(int key) => Enters++;

        public void ExitBranch() => Exits++;
    }
}
