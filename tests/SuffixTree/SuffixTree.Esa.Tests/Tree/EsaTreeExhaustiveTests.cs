using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Tree;

/// <summary>
/// Heavy checks of work package 1b.1. Category <c>Exhaustive</c>, <c>[Explicit]</c>:
/// <c>dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter Category=Exhaustive</c>.
/// Every string over sigma 2 up to length 14, sigma 3 up to 9, sigma 4 up to 7 (the lists of PROGRESS-S1.md, "Test strategy"): the <c>Traverse</c> stream of the ESA (both
/// single-array layouts) equals the Ukkonen stream, the node statistics are equal; the persistent v6 tree is
/// compared on the strings up to length 11 / 7 / 5; random and adversarial texts up to 2000 symbols against the naive model.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log 1b.1)")]
public class EsaTreeExhaustiveTests
{
    private static List<string> Stream(Action<ISuffixTreeVisitor> traverse)
    {
        var recorder = new TraverseRecorder();
        traverse(recorder);
        return recorder.Events;
    }

    private static void CheckString(string s, bool withV6)
    {
        UkkonenTree ukkonen = UkkonenTree.Build(s);
        List<string> expected = Stream(ukkonen.Traverse);
        foreach (ChildLayoutKind layout in new[] { ChildLayoutKind.Tagged, ChildLayoutKind.Packed })
        {
            EsaIndex esa = EsaIndex.Build(s, layout);
            if (!expected.SequenceEqual(Stream(esa.Traverse))) Assert.Fail($"Traverse stream of ESA({layout}) differs from Ukkonen on '{s}'");
            Assert.That(esa.Statistics.NodeCount, Is.EqualTo(ukkonen.NodeCount), $"NodeCount '{s}'");
            Assert.That(esa.Statistics.LeafCount, Is.EqualTo(ukkonen.LeafCount), $"LeafCount '{s}'");
            Assert.That(esa.Statistics.MaxDepth, Is.EqualTo(ukkonen.MaxDepth), $"MaxDepth '{s}'");
        }
        if (withV6)
        {
            using PersistentSuffixTree v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(s));
            if (!Stream(v6.Traverse).SequenceEqual(expected)) Assert.Fail($"v6 stream differs from Ukkonen on '{s}' (oracle disagreement)");
        }
    }

    [Test]
    public void AllBinaryStrings_UpTo14()
    {
        foreach (string s in ExhaustiveStrings.All(2, 14)) CheckString(s, withV6: s.Length <= 11);
    }

    [Test]
    public void AllTernaryStrings_UpTo9()
    {
        foreach (string s in ExhaustiveStrings.All(3, 9)) CheckString(s, withV6: s.Length <= 7);
    }

    [Test]
    public void AllFourLetterStrings_UpTo7()
    {
        foreach (string s in ExhaustiveStrings.All(4, 7)) CheckString(s, withV6: s.Length <= 5);
    }

    [Test]
    public void RandomAndAdversarialTexts_UpTo2000_EqualTheNaiveModelAndUkkonen()
    {
        var rng = new SplitMix64(0x1B30);
        for (int i = 0; i < 60; i++)
        {
            int n = rng.NextInt(2001);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 3 == 0))
            {
                List<string> naive = Stream(v => NaiveEsaModel.EmitTraverse(text, v));
                foreach (ChildLayoutKind layout in new[] { ChildLayoutKind.Tagged, ChildLayoutKind.Packed, ChildLayoutKind.Reference3 })
                {
                    EsaIndex esa = EsaIndex.Build(text, layout);
                    if (!naive.SequenceEqual(Stream(esa.Traverse))) Assert.Fail($"{name} n {n} seed {seed}: ESA({layout}) differs from the naive model");
                }
                List<string> ukkonen = Stream(UkkonenTree.Build(text).Traverse);
                if (!naive.SequenceEqual(ukkonen)) Assert.Fail($"{name} n {n} seed {seed}: naive model differs from Ukkonen (oracle disagreement)");
            }
        }
    }
}
