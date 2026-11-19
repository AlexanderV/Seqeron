using System.Collections.Concurrent;
using System.Text;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>What the concurrency test counted: repetitions, trees, thread runs and the checksums that were compared.</summary>
internal sealed class ConcurrencyStats
{
    private long _repetitions, _threadRuns, _checksums, _packBuildsChecked;

    internal void Repetition(int threads, int checksums, bool packChecked)
    {
        Interlocked.Increment(ref _repetitions);
        Interlocked.Add(ref _threadRuns, threads);
        Interlocked.Add(ref _checksums, checksums);
        if (packChecked) Interlocked.Increment(ref _packBuildsChecked);
    }

    internal long Repetitions => Interlocked.Read(ref _repetitions);

    internal long ThreadRuns => Interlocked.Read(ref _threadRuns);

    internal long Checksums => Interlocked.Read(ref _checksums);

    internal long PackBuildsChecked => Interlocked.Read(ref _packBuildsChecked);
}

/// <summary>
/// The shared-index test (work package 1b.8, gate G3): one tree, several threads, every read operation of the contract over a fixed script of queries (<see cref="FuzzDigest"/> for the members of
/// <see cref="ISuffixTree"/> and <see cref="FuzzContractDigest"/> for the navigation and suffix-array members), and the checksum of every thread's answers equal to the checksum of the single-thread run, on a tree
/// that is new for every repetition so that the once-only build of the navigation pack is raced by the threads that are the first to need it and is checked through <c>NavigationBuildCount</c>. The kinds: the
/// in-memory tree, the file-backed tree at the Header level in the Checked mode and at the Structure level in the Trusted mode.
/// </summary>
internal static class FullFuzzConcurrency
{
    internal static readonly string[] Kinds = ["memory-default", "file-header-checked", "file-structure-trusted"];

    /// <summary>Runs <paramref name="repetitions"/> repetitions of <paramref name="threads"/> threads on each kind of tree over <paramref name="text"/>; null when every checksum was equal.</summary>
    internal static string? Run(string text, ulong seed, FullFuzzContext ctx, ConcurrencyStats stats, int repetitions, int threads, IReadOnlyList<string>? kinds = null)
    {
        FuzzScript script = FuzzScript.Create(text, seed);
        DigestLimits limits = DigestLimits.Default;
        bool fileCapable = text.Length > 0 && !text.Any(c => c > 255);
        string input = Path.Combine(ctx.WorkDir, "cc-in.txt");
        string container = Path.Combine(ctx.WorkDir, "cc.sesa");
        try
        {
            if (fileCapable)
            {
                File.WriteAllBytes(input, Encoding.Latin1.GetBytes(text));
                SesaBuilder.Build(input, container, FuzzSubjects.ContainerOptions(text.Length, FullFuzzConfig.Of(seed)));
            }
            // the single-thread answer: the in-memory tree, new, one thread
            (List<string> expected, ulong expectedChecksum, int expectedPackBuilds) = Single(Open("memory-default", text, container, ctx), text, seed, script, limits, ctx);
            foreach (string kind in kinds ?? Kinds)
            {
                if (kind.StartsWith("file", StringComparison.Ordinal) && !fileCapable) continue;
                (List<string> single, ulong singleChecksum, _) = Single(Open(kind, text, container, ctx), text, seed, script, limits, ctx);
                if (singleChecksum != expectedChecksum)
                    return $"[{kind}] single thread: {FuzzDigest.Diff(single, expected, kind, "memory-default") ?? "the checksum differs"}";
                for (int r = 0; r < repetitions; r++)
                {
                    string? failure = Repeat(kind, text, container, seed, script, limits, expected, expectedChecksum, expectedPackBuilds, threads, r, ctx, stats);
                    if (failure is not null) return failure;
                }
            }
            return null;
        }
        finally
        {
            File.Delete(input);
            File.Delete(container);
        }
    }

    private static CountingTree Open(string kind, string text, string container, FullFuzzContext ctx)
    {
        int e = ctx.Counters.EngineIndex(kind);
        EsaSuffixTree tree = kind switch
        {
            "memory-default" => EsaSuffixTree.Build(text),
            "file-header-checked" => (EsaSuffixTree)SuffixIndexFactory.Open(container, new EsaOpenOptions { Validation = EsaValidationLevel.Header, Access = EsaAccessMode.Checked }),
            "file-structure-trusted" => (EsaSuffixTree)SuffixIndexFactory.Open(container, new EsaOpenOptions { Validation = EsaValidationLevel.Structure, Access = EsaAccessMode.Trusted }),
            _ => throw new ArgumentException($"unknown kind '{kind}'", nameof(kind)),
        };
        return new CountingTree(tree, ctx.Counters, e);
    }

    private static List<string> Script(CountingTree tree, string text, ulong seed, FuzzScript script, DigestLimits limits, bool contractFirst) =>
        FuzzSteps.Run(FuzzSteps.Concurrency, () =>
        {
            if (contractFirst)
            {
                List<string> c = FuzzContractDigest.Compute(tree, text, seed);
                List<string> d = FuzzDigest.Compute(tree, text, script, limits);
                return [.. d, .. c];
            }
            List<string> digest = FuzzDigest.Compute(tree, text, script, limits);
            digest.AddRange(FuzzContractDigest.Compute(tree, text, seed));
            return digest;
        });

    private static (List<string> Lines, ulong Checksum, int PackBuilds) Single(CountingTree tree, string text, ulong seed, FuzzScript script, DigestLimits limits, FullFuzzContext ctx)
    {
        using (tree)
        {
            List<string> lines = Script(tree, text, seed, script, limits, contractFirst: false);
            return (lines, FuzzDigest.Checksum(lines), tree.Inner.Index.NavigationBuildCount);
        }
    }

    private static string? Repeat(
        string kind, string text, string container, ulong seed, FuzzScript script, DigestLimits limits, List<string> expected, ulong expectedChecksum,
        int expectedPackBuilds, int threads, int repetition, FullFuzzContext ctx, ConcurrencyStats stats)
    {
        using CountingTree tree = Open(kind, text, container, ctx);
        var results = new (List<string> Lines, ulong Checksum)[threads];
        var errors = new ConcurrentBag<string>();
        using var gate = new Barrier(threads);
        var workers = Enumerable.Range(0, threads).Select(t => new Thread(() =>
        {
            try
            {
                gate.SignalAndWait();
                // odd repetitions: every thread runs the digest first; even repetitions: half of the threads start with the contract surface (the lines are the same, only the order of the calls differs), so that different lazily built structures are first used concurrently
                List<string> lines = Script(tree, text, seed, script, limits, contractFirst: repetition % 2 == 0 && t % 2 == 1);
                results[t] = (lines, FuzzDigest.Checksum(lines));
            }
            catch (Exception ex)
            {
                errors.Add($"thread {t}: {ex.GetType().Name}: {ex.Message}");
            }
        })).ToArray();
        foreach (Thread worker in workers) worker.Start();
        foreach (Thread worker in workers)
            if (!worker.Join(TimeSpan.FromMinutes(5))) return $"[{kind}] repetition {repetition}: a thread did not finish (a hang)";
        if (!errors.IsEmpty) return $"[{kind}] repetition {repetition}: {string.Join(" | ", errors.Take(3))}";
        for (int t = 0; t < threads; t++)
        {
            if (!results[t].Lines.SequenceEqual(expected))
                return $"[{kind}] repetition {repetition}, thread {t}: {FuzzDigest.Diff(results[t].Lines, expected, "the thread", "the single-thread run") ?? "different lines"}";
            if (results[t].Checksum != expectedChecksum)
                return $"[{kind}] repetition {repetition}, thread {t}: checksum {results[t].Checksum:x16} against the single-thread {expectedChecksum:x16}";
        }
        int builds = tree.Inner.Index.NavigationBuildCount;
        if (builds != expectedPackBuilds) return $"[{kind}] repetition {repetition}: the navigation pack was built {builds} times, once is expected ({expectedPackBuilds} in the single-thread run)";
        stats.Repetition(threads, threads, packChecked: expectedPackBuilds == 1);
        return null;
    }
}
