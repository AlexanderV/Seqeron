using System.Diagnostics;
using System.Text.RegularExpressions;
using SuffixTree.Esa.Tests.Contract;
using SuffixTree.Esa.Tests.Facade;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>Per-engine time and check counts of a run (build time of the engine, time of its checks, number of checks), thread safe.</summary>
internal sealed class EngineStats
{
    private readonly string[] _names = [.. FuzzSubjects.EngineNames, .. FuzzSubjects.ReferenceNames];
    private readonly long[] _checks;
    private readonly long[] _checkTicks;
    private readonly long[] _buildTicks;

    internal EngineStats()
    {
        _checks = new long[_names.Length];
        _checkTicks = new long[_names.Length];
        _buildTicks = new long[_names.Length];
    }

    /// <summary>The names, engines under test first, then the references.</summary>
    internal IReadOnlyList<string> Names => _names;

    internal void AddBuild(string engine, long ticks) => Interlocked.Add(ref _buildTicks[Array.IndexOf(_names, engine)], ticks);

    internal void AddCheck(string engine, long ticks)
    {
        int i = Array.IndexOf(_names, engine);
        Interlocked.Increment(ref _checks[i]);
        Interlocked.Add(ref _checkTicks[i], ticks);
    }

    internal long Checks(string engine) => Interlocked.Read(ref _checks[Array.IndexOf(_names, engine)]);

    internal double CheckSeconds(string engine) => Interlocked.Read(ref _checkTicks[Array.IndexOf(_names, engine)]) / (double)Stopwatch.Frequency;

    internal double BuildSeconds(string engine) => Interlocked.Read(ref _buildTicks[Array.IndexOf(_names, engine)]) / (double)Stopwatch.Frequency;
}

/// <summary>What the check of a text needs besides the text and its configuration.</summary>
internal sealed class FullFuzzContext
{
    /// <summary>The call counts (a throw-away instance during minimization).</summary>
    internal required MemberCounters Counters { get; init; }

    /// <summary>The directory of the container builder (one per worker).</summary>
    internal required string WorkDir { get; init; }

    /// <summary>Per-engine time.</summary>
    internal EngineStats Stats { get; init; } = new();

    /// <summary>Optional fault injection of the harness tests: the count that an engine reports for a pattern, from the engine, the pattern and the true count.</summary>
    internal Func<string, string, int, int>? CountFault { get; init; }

    /// <summary>Texts above this length get no v6 tree (the harness bands end at 10^5 symbols, where v6 builds in about 70 ms, so the default is every case).</summary>
    internal int V6UpTo { get; init; } = int.MaxValue;

    /// <summary>Incremented for every case whose text has a code unit above 255 (no file-backed engines).</summary>
    internal long SkippedForFile;
}

/// <summary>
/// The check of one text on every engine (work package 1b.8): the in-memory trees and the file-backed trees each run the definition-level checks of <see cref="FacadeChecks"/> where the text is short, the search against
/// the definitions otherwise, the digest of every member against the Ukkonen tree and the v6 tree, the contract checker against its oracle, the contract digest against the first in-memory tree, and the suffix-array and
/// LCP arrays through <see cref="ISuffixArrayAccess"/> against the naive or the verified arrays. A failure is a <see cref="FuzzMismatch"/> whose operation names the engine and the member group, so that the minimizer keeps
/// the same failure while it shrinks.
/// </summary>
internal static class FullFuzzCheck
{
    /// <summary>Texts up to this length are checked against the definitions with <see cref="FacadeChecks.CheckEverything"/>.</summary>
    internal const int DefinitionsUpTo = 120;

    /// <summary>Texts up to this length are checked by the contract checker.</summary>
    internal const int ContractUpTo = 150;

    private static readonly string[] Members =
    [
        "FindAllLongestRepeatedSubstrings", "LongestRepeatedSubstringMemory", "LongestRepeatedSubstring", "CountDistinctSubstringsByLength", "CountDistinctSubstrings", "GetAllSuffixes", "EnumerateSuffixes",
        "FindMaximalRepeatedPairs", "FindAllDistinctLongestCommonSubstrings", "FindAllLongestCommonSubstrings", "LongestCommonSubstring", "MEM", "MUM", "anchors", "FindAllOccurrences", "CountOccurrences",
        "Contains", "Traverse", "PrintTree", "NodeCount", "LeafCount", "MaxDepth", "IsEmpty", "Text is not", "info",
    ];

    /// <summary>The member group named in a failure message (the first known name), the stable part of the operation.</summary>
    internal static string OpOf(string message)
    {
        string body = message;
        int cut = message.IndexOf("' query", StringComparison.Ordinal);
        if (cut < 0) cut = message.IndexOf("text '", StringComparison.Ordinal);
        if (cut >= 0) body = message[cut..];
        string best = "assert";
        int at = int.MaxValue;
        foreach (string member in Members)
        {
            int i = body.IndexOf(member, StringComparison.Ordinal);
            if (i >= 0 && i < at)
            {
                at = i;
                best = member;
            }
        }
        return best;
    }

    private static string StepOf(string? step) => step is null ? "?" : step[(step.IndexOf('.', StringComparison.Ordinal) + 1)..];

    private static string DigestOp(string diff)
    {
        string label = diff[..diff.IndexOf(':', StringComparison.Ordinal)];
        return "digest:" + Regex.Replace(label, @"\[[^\]]*\]", string.Empty);
    }

    /// <summary>Checks <paramref name="text"/> on every engine; null when every check passed.</summary>
    internal static FuzzMismatch? CheckText(string text, FullFuzzConfig cfg, FullFuzzContext ctx)
    {
        int n = text.Length;
        FuzzSubjects set;
        try
        {
            set = FuzzSubjects.Build(text, cfg, ctx.Counters, ctx.WorkDir, ctx.V6UpTo);
        }
        catch (Exception ex)
        {
            return new FuzzMismatch("build:" + ex.GetType().Name, ex.Message);
        }

        using (set)
        {
            if (set.FileSkipped) Interlocked.Increment(ref ctx.SkippedForFile);
            ctx.Stats.AddBuild("ukkonen", set.ReferenceBuildTicks);
            FuzzScript script = FuzzScript.Create(text, cfg.Seed);
            DigestLimits limits = DigestLimits.Default;
            List<string> ukkonenDigest;
            List<string>? v6Digest = null;
            try
            {
                ukkonenDigest = Timed(ctx, "ukkonen", () => FuzzSteps.Run(FuzzSteps.Digest, () => FuzzDigest.Compute(set.Ukkonen, text, script, limits)));
                if (set.V6 is not null) v6Digest = Timed(ctx, "v6", () => FuzzSteps.Run(FuzzSteps.Digest, () => FuzzDigest.Compute(set.V6, text, script, limits)));
            }
            catch (Exception ex)
            {
                return new FuzzMismatch("references:exception:" + ex.GetType().Name, ex.Message);
            }
            if (v6Digest is not null && FuzzDigest.Diff(ukkonenDigest, v6Digest, "ukkonen", "v6") is { } referenceDiff)
                return new FuzzMismatch("references:" + DigestOp(referenceDiff), referenceDiff);

            var oracles = new CaseOracles(text);
            List<string>? firstContract = null;
            foreach (FuzzSubject subject in set.Subjects)
            {
                ctx.Stats.AddBuild(subject.Engine, subject.BuildTicks);
                long t = Stopwatch.GetTimestamp();
                FuzzMismatch? mismatch = CheckSubject(subject, set, text, cfg, script, limits, ukkonenDigest, v6Digest, oracles, ctx, ref firstContract);
                ctx.Stats.AddCheck(subject.Engine, Stopwatch.GetTimestamp() - t);
                if (mismatch is not null) return mismatch;
            }
        }
        return null;
    }

    private static T Timed<T>(FullFuzzContext ctx, string engine, Func<T> body)
    {
        long t = Stopwatch.GetTimestamp();
        try
        {
            return body();
        }
        finally
        {
            ctx.Stats.AddCheck(engine, Stopwatch.GetTimestamp() - t);
        }
    }

    private static FuzzMismatch? CheckSubject(
        FuzzSubject subject, FuzzSubjects set, string text, FullFuzzConfig cfg, FuzzScript script, DigestLimits limits,
        List<string> ukkonenDigest, List<string>? v6Digest, CaseOracles oracles, FullFuzzContext ctx, ref List<string>? firstContract)
    {
        int n = text.Length;
        string engine = subject.Engine;
        string where = $"[{engine}] n={n}";
        FuzzSteps.TakeFailedStep();
        try
        {
            if (ctx.CountFault is { } fault) subject.Tree.CountFault = (pattern, count) => fault(engine, pattern, count);
            ISuffixTree view = subject.Tree;
            if (n <= DefinitionsUpTo)
                FuzzSteps.Run(FuzzSteps.Definitions, () => FacadeChecks.CheckEverything(view, set.Ukkonen, set.V6, text, cfg.Seed, where, maxPatternLength: 4, printTree: n <= 64));
            else
                FuzzSteps.Run(FuzzSteps.Search, () => FacadeChecks.CheckSearch(view, text, script.Patterns, where));

            List<string> digest = FuzzSteps.Run(FuzzSteps.Digest, () => FuzzDigest.Compute(view, text, script, limits));
            if (FuzzDigest.Diff(digest, ukkonenDigest, engine, "ukkonen") is { } d1) return new FuzzMismatch($"{engine}:{DigestOp(d1)}", d1);
            if (v6Digest is not null && FuzzDigest.Diff(digest, v6Digest, engine, "v6") is { } d2) return new FuzzMismatch($"{engine}:{DigestOp(d2)}", d2);

            if (n <= ContractUpTo)
                FuzzSteps.Run(FuzzSteps.Contract, () => EsaContractChecker.CheckAll(subject.Tree, oracles.Contract, where, n <= 24 ? CheckOptions.Thorough : CheckOptions.Sampled));
            List<string> contract = FuzzSteps.Run(FuzzSteps.ContractDigest, () => FuzzContractDigest.Compute(subject.Tree, text, cfg.Seed));
            firstContract ??= contract;
            if (!ReferenceEquals(firstContract, contract) && FuzzDigest.Diff(contract, firstContract, engine, "the first in-memory tree") is { } d3)
                return new FuzzMismatch($"{engine}:contract-digest:{FuzzDigest.LabelOf(contract[0])}", d3);

            string? arrays = FuzzSteps.Run(FuzzSteps.Arrays, () => SuffixArrayLevel(subject.Tree, text, oracles));
            if (arrays is not null) return new FuzzMismatch($"{engine}:arrays", arrays);
            return null;
        }
        catch (AssertionException ex)
        {
            return new FuzzMismatch($"{engine}:{StepOf(FuzzSteps.TakeFailedStep())}:{OpOf(ex.Message)}", ex.Message);
        }
        catch (Exception ex)
        {
            return new FuzzMismatch($"{engine}:{StepOf(FuzzSteps.TakeFailedStep())}:exception:{ex.GetType().Name}", ex.Message);
        }
    }

    /// <summary>The suffix array and the LCP array through <see cref="ISuffixArrayAccess"/> against the naive arrays (up to 2000 symbols) or the verifier and Kasai's LCP (above).</summary>
    private static string? SuffixArrayLevel(ISuffixArrayAccess tree, string text, CaseOracles oracles)
    {
        int n = text.Length;
        var sa = new int[n + 1];
        var lcp = new int[n + 1];
        for (int r = 0; r <= n; r++)
        {
            sa[r] = (int)tree.SuffixAt(r);
            lcp[r] = (int)tree.LcpAt(r);
        }
        if (n <= FuzzOracle.NaiveLimit)
        {
            int[] expected = oracles.NaiveSa;
            for (int r = 0; r <= n; r++)
                if (sa[r] != expected[r]) return $"SuffixAt({r}) = {sa[r]}, the naive suffix array says {expected[r]}";
        }
        else if (SuffixArrayVerifier.Check(text, sa) is { } verdict)
        {
            return "SuffixAt: " + verdict;
        }
        int[] expectedLcp = oracles.ExpectedLcp(sa);
        for (int r = 0; r <= n; r++)
            if (lcp[r] != expectedLcp[r]) return $"LcpAt({r}) = {lcp[r]}, expected {expectedLcp[r]}";
        var rng = new SplitMix64((ulong)n);
        for (int k = 0; k < Math.Min(64, n + 1); k++)
        {
            int r = rng.NextInt(n + 1);
            if (tree.RankOf(sa[r]) != r) return $"RankOf(SuffixAt({r})) = {tree.RankOf(sa[r])}";
        }
        return null;
    }

    /// <summary>The oracles that are the same for every engine of a case, computed on first use.</summary>
    private sealed class CaseOracles(string text)
    {
        private int[]? _naiveSa;
        private int[]? _lcp;
        private IndexContractOracle? _contract;

        internal int[] NaiveSa => _naiveSa ??= NaiveSuffixArray.Build(text);

        internal IndexContractOracle Contract => _contract ??= new IndexContractOracle(text);

        internal int[] ExpectedLcp(int[] verifiedSa) => _lcp ??= text.Length <= FuzzOracle.NaiveLimit ? NaiveLcp.Build(text, NaiveSa) : KasaiLcp.Build(text, verifiedSa);
    }
}
