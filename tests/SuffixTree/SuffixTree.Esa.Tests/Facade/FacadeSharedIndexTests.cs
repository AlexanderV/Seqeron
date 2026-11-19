using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Fuzz;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The shared-index test of work package 1b.8 (gate G3): one tree, four threads, every read operation of the contract over a fixed script of queries (the members of <see cref="ISuffixTree"/> through
/// <see cref="FuzzDigest"/>, the navigation and suffix-array members through <see cref="FuzzContractDigest"/>), the checksum of every thread's answers equal to the checksum of the single-thread run, in at least 20
/// repetitions on a tree that is new for every repetition, on the in-memory tree and on the file-backed tree at both ends of the validation levels (Header and Checked, Structure and Trusted); the once-only build of the
/// navigation pack is checked through <c>NavigationBuildCount</c> in every repetition. The same code runs in the budgeted harness (<see cref="FullFuzzConcurrency"/>) on generated texts.
/// </summary>
[TestFixture]
public class FacadeSharedIndexTests
{
    private const int Repetitions = 20;
    private const int Threads = 4;

    private static void Check(string text, ulong seed, string[]? kinds = null)
    {
        using var dir = new TempDir();
        var ctx = new FullFuzzContext { Counters = new MemberCounters(FuzzSubjects.EngineNames), WorkDir = dir.Path };
        var stats = new ConcurrencyStats();
        string? failure = FullFuzzConcurrency.Run(text, seed, ctx, stats, Repetitions, Threads, kinds);
        Assert.That(failure, Is.Null);
        int kindCount = kinds?.Length ?? FullFuzzConcurrency.Kinds.Length;
        Assert.Multiple(() =>
        {
            Assert.That(stats.Repetitions, Is.EqualTo(Repetitions * kindCount), "repetitions");
            Assert.That(stats.ThreadRuns, Is.EqualTo(Repetitions * kindCount * Threads), "thread runs");
            Assert.That(stats.Checksums, Is.EqualTo(Repetitions * kindCount * Threads), "checksums compared with the single-thread run");
            Assert.That(stats.PackBuildsChecked, Is.EqualTo(Repetitions * kindCount), "repetitions in which the pack was built exactly once");
        });
    }

    [Test]
    public void RandomDna_EveryMember_FourThreads_20Repetitions_InMemoryAndFileBacked() => Check(StringFamilies.RandomDna(1000, 41) + "ACGTACGTAC", 11);

    [Test]
    public void FibonacciWord_EveryMember_FourThreads_20Repetitions_InMemoryAndFileBacked() => Check(StringFamilies.Fibonacci(800), 12);

    [Test]
    public void RunsAndRepeats_EveryMember_FourThreads_20Repetitions_InMemoryAndFileBacked() => Check(new string('a', 300) + "b" + new string('a', 200) + StringFamilies.Random(150, 3, 5), 13);

    [Test]
    public void ATextAboveByte_InMemoryOnly_FourThreads_20Repetitions() => Check("abĀāabĀĂ" + StringFamilies.Random(600, 5, 7) + "Āāab", 14, ["memory-default"]);

    [Test]
    public void ThePackBuild_IsDoneOnce_WhenFourThreadsAreTheFirstToNeedIt_OnTheFileBackedTree()
    {
        string text = StringFamilies.RandomDna(5000, 3);
        string query = string.Concat(text.AsSpan(500, 400), "TTTT", text.AsSpan(3000, 300));
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, text);
        for (int round = 0; round < Repetitions; round++)
        {
            using EsaSuffixTree tree = SesaTestFile.Open(path, round % 2 == 0 ? EsaValidationLevel.Header : EsaValidationLevel.Structure, round % 2 == 0 ? EsaAccessMode.Checked : EsaAccessMode.Trusted);
            Assert.That(tree.Index.NavigationBuilt, Is.False);
            var results = new string[Threads];
            using var gate = new Barrier(Threads);
            Thread[] workers = [.. Enumerable.Range(0, Threads).Select(t => new Thread(() =>
            {
                gate.SignalAndWait();
                results[t] = string.Join(";", tree.FindMaximalExactMatches(query, 12));
            }))];
            foreach (Thread worker in workers) worker.Start();
            foreach (Thread worker in workers) worker.Join();
            Assert.Multiple(() =>
            {
                Assert.That(tree.Index.NavigationBuildCount, Is.EqualTo(1), "the pack is built once");
                Assert.That(results.Distinct().Count(), Is.EqualTo(1), "every thread got the same answer");
            });
        }
    }
}
