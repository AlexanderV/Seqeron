using System.Collections.Concurrent;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// Several readers over one facade (ADR-001 rule 4: engines are safe for concurrent readers, lazily built auxiliary structures are published once): the answers of every thread equal
/// the sequential answers, and the navigation pack is built exactly once however many threads are the first to need it.
/// </summary>
[TestFixtureSource(typeof(SuffixTreeBackends), nameof(SuffixTreeBackends.EsaFixtures))]
public class FacadeConcurrencyTests(SuffixTreeBackend backend) : PortedBackendTestBase(backend)
{
    private const int Threads = 8;

    private static void RunThreads(int count, Action<int> body)
    {
        var errors = new ConcurrentBag<Exception>();
        var gate = new Barrier(count);
        var threads = Enumerable.Range(0, count).Select(i => new Thread(() =>
        {
            try
            {
                gate.SignalAndWait();
                body(i);
            }
            catch (Exception ex)
            {
                errors.Add(ex);
            }
        })).ToList();
        foreach (Thread t in threads) t.Start();
        foreach (Thread t in threads) t.Join();
        if (!errors.IsEmpty) throw new AggregateException(errors);
    }

    [Test]
    public void ManyReaders_ContainsCountFindAll_MatchTheSequentialAnswers()
    {
        string text = StringFamilies.RandomDna(4000, 77) + "ACGTACGTAC";
        ISuffixTree tree = Tree(text);
        var patterns = Enumerable.Range(0, 200).Select(i => text.Substring((i * 37) % (text.Length - 12), 3 + i % 9)).ToArray();
        int[] counts = patterns.Select(tree.CountOccurrences).ToArray();
        int[][] positions = patterns.Select(p => tree.FindAllOccurrences(p).Order().ToArray()).ToArray();
        Assert.DoesNotThrow(() => RunThreads(Threads, t =>
        {
            for (int round = 0; round < 5; round++)
                for (int i = t; i < patterns.Length; i += Threads)
                {
                    if (tree.CountOccurrences(patterns[i]) != counts[i]) throw new InvalidOperationException($"count of '{patterns[i]}'");
                    if (!tree.FindAllOccurrences(patterns[i]).Order().SequenceEqual(positions[i])) throw new InvalidOperationException($"positions of '{patterns[i]}'");
                    if (!tree.Contains(patterns[i])) throw new InvalidOperationException("contains");
                }
        }));
    }

    [Test]
    public void ManyFirstUsers_OfThePack_BuildItOnce_AndAgreeWithTheSequentialAnswer()
    {
        string text = StringFamilies.RandomDna(6000, 5);
        string query = string.Concat(text.AsSpan(1000, 800), "TTTT", text.AsSpan(3000, 500));
        EsaSuffixTree tree = (EsaSuffixTree)Tree(text);
        Assert.That(tree.Index.NavigationBuilt, Is.False);
        var results = new ConcurrentBag<string>();
        Assert.DoesNotThrow(() => RunThreads(Threads, _ =>
        {
            var mems = tree.FindMaximalExactMatches(query, 12);
            results.Add(string.Join(";", mems));
        }));
        Assert.Multiple(() =>
        {
            Assert.That(tree.Index.NavigationBuildCount, Is.EqualTo(1), "the pack is built once");
            Assert.That(results.Distinct().Count(), Is.EqualTo(1), "every thread got the same answer");
            Assert.That(results.First(), Is.EqualTo(string.Join(";", global::SuffixTree.SuffixTree.Build(text).FindMaximalExactMatches(query, 12))));
        });
    }

    [Test]
    public void MixedOperations_AcrossTheWholeSurface_ThreadSafe()
    {
        string text = StringFamilies.Fibonacci(2500);
        EsaSuffixTree tree = (EsaSuffixTree)Tree(text);
        global::SuffixTree.SuffixTree reference = global::SuffixTree.SuffixTree.Build(text);
        string query = string.Concat(text.AsSpan(100, 300), "ab", text.AsSpan(900, 200));
        string lrs = reference.LongestRepeatedSubstring();
        long distinct = reference.CountDistinctSubstrings();
        var mem = reference.FindMaximalExactMatches(query, 6);
        var mum = reference.FindMaximalUniqueMatches(query, 6);
        var lcs = reference.LongestCommonSubstring(query);
        Assert.DoesNotThrow(() => RunThreads(Threads, t =>
        {
            for (int round = 0; round < 6; round++)
            {
                switch ((t + round) % 8)
                {
                    case 0:
                        if (tree.LongestRepeatedSubstring().Length != lrs.Length) throw new InvalidOperationException("lrs");
                        break;
                    case 1:
                        if (tree.CountDistinctSubstrings() != distinct) throw new InvalidOperationException("distinct");
                        break;
                    case 2:
                        if (!tree.FindMaximalExactMatches(query, 6).SequenceEqual(mem)) throw new InvalidOperationException("mem");
                        break;
                    case 3:
                        if (!tree.FindMaximalUniqueMatches(query, 6).SequenceEqual(mum)) throw new InvalidOperationException("mum");
                        break;
                    case 4:
                        if (tree.LongestCommonSubstring(query) != lcs) throw new InvalidOperationException("lcs");
                        break;
                    case 5:
                        if (tree.GetAllSuffixes().Count != text.Length) throw new InvalidOperationException("suffixes");
                        break;
                    case 6:
                        if (tree.Traverse_Count() != reference.Traverse_Count()) throw new InvalidOperationException("traverse");
                        break;
                    default:
                        {
                            long rank = 1 + (t * 31 + round * 7) % text.Length;
                            if (tree.RankOf(tree.SuffixAt(rank)) != rank) throw new InvalidOperationException("rank");
                            if (!tree.TryExtendRight(tree.Root, text[(int)tree.SuffixAt(rank)], out IndexLocus locus) || tree.Count(locus) < 1) throw new InvalidOperationException("extend");
                            break;
                        }
                }
            }
        }));
        Assert.That(tree.Index.NavigationBuildCount, Is.EqualTo(1));
    }

    [Test]
    public void ManyTrees_BuiltAndDisposedConcurrently_DoNotInterfere()
    {
        Assert.DoesNotThrow(() => RunThreads(Threads, t =>
        {
            for (int i = 0; i < 6; i++)
            {
                string text = StringFamilies.Random(300 + t * 10 + i, 3, (ulong)(t * 100 + i));
                ISuffixTree tree = Backend.FromString(text);
                try
                {
                    string query = text.Length > 40 ? text.Substring(20, 20) : text;
                    if (tree.CountOccurrences(query) < 1) throw new InvalidOperationException("missing own substring");
                    if (tree.FindMaximalExactMatches(query + "z", 5).Count < 1) throw new InvalidOperationException("no MEM");
                }
                finally
                {
                    (tree as IDisposable)?.Dispose();
                }
            }
        }));
    }
}
