using System.Collections.Concurrent;
using System.Diagnostics;
using System.Globalization;
using System.Text;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>One lane of the budgeted run: a name and its share of the time budget.</summary>
internal sealed record FullFuzzLane(string Name, double Share);

/// <summary>Parameters of one run of the full differential fuzz (work package 1b.8).</summary>
internal sealed class FullFuzzOptions
{
    /// <summary>The lanes of the matrix with the shares of the time budget: the three size bands of the differential fuzz over every engine, the container corruption fuzz and the shared-index test.</summary>
    internal static readonly FullFuzzLane[] DefaultLanes =
    [
        new("tiny", 0.24), new("small", 0.24), new("large", 0.16), new("xl", 0.06), new("corruption", 0.18), new("concurrency", 0.12),
    ];

    /// <summary>The largest band of the differential fuzz: 10^5 to 5 x 10^5 symbols (a case costs seconds, almost all of it in the references: the v6 tree builds in about 1.5 s per 10^6 symbols, the digest of the Ukkonen and the v6 tree is the slow part).</summary>
    internal static readonly FuzzTier XlTier = new("n<=5e5", 100_001, 500_000, true, 0, 6);

    /// <summary>The size band of the concurrency lane (its texts are short: every case runs 3 kinds x repetitions x threads).</summary>
    internal static readonly FuzzTier ConcurrencyTier = new("n<=1500", 20, 1500, true, 0, 4);

    internal ulong MasterSeed { get; init; } = 1;

    /// <summary>The whole run; every lane gets its share.</summary>
    internal TimeSpan Budget { get; init; } = TimeSpan.FromSeconds(30);

    /// <summary>Worker threads of the differential and corruption lanes (the concurrency lane runs one worker, which starts four threads per repetition).</summary>
    internal int Threads { get; init; } = Environment.ProcessorCount;

    internal IReadOnlyList<FullFuzzLane> Lanes { get; init; } = DefaultLanes;

    /// <summary>Stops a lane after this many cases (tests), whatever the budget.</summary>
    internal long? MaxCasesPerLane { get; init; }

    /// <summary>The run stops after this many failures.</summary>
    internal int MaxFailures { get; init; } = 3;

    /// <summary>Time the minimizer may spend on one failure.</summary>
    internal TimeSpan MinimizeLimit { get; init; } = TimeSpan.FromSeconds(30);

    /// <summary>Repetitions per kind of tree in the concurrency lane (G3 asks for at least 20).</summary>
    internal int ConcurrencyRepetitions { get; init; } = 20;

    internal int ConcurrencyThreads { get; init; } = 4;

    /// <summary>A worker that stays in one case longer than this is a hang: the run is ended with the case named (null: no watchdog).</summary>
    internal TimeSpan? HangLimit { get; init; }

    /// <summary>When set, every worker writes the case it is in to a file here before it runs it (the case that was running when the process died).</summary>
    internal string? JournalDirectory { get; init; }

    /// <summary>Fault injection of the harness tests: the count that an engine reports for a pattern, from the engine, the pattern and the true count.</summary>
    internal Func<string, string, int, int>? CountFault { get; init; }

    internal Action<string>? Log { get; init; }
}

/// <summary>The result of one lane.</summary>
internal sealed class FullFuzzLaneResult(string name, TimeSpan budget)
{
    internal string Name { get; } = name;

    internal TimeSpan Budget { get; } = budget;

    internal TimeSpan Elapsed { get; set; }

    internal long Cases;

    internal long Units;

    internal long Skipped;

    internal string FirstCaseId { get; set; } = "";

    internal string LastCaseId { get; set; } = "";

    internal ConcurrentDictionary<string, long> PerFamily { get; } = new();

    internal ConcurrentBag<FuzzFailure> Failures { get; } = [];

    internal string Extra { get; set; } = "";
}

/// <summary>The outcome of a run.</summary>
internal sealed class FullFuzzSummary
{
    internal required ulong MasterSeed { get; init; }

    internal required MemberCounters Counters { get; init; }

    internal required EngineStats Stats { get; init; }

    internal List<FullFuzzLaneResult> Lanes { get; } = [];

    internal CorruptionStats Corruption { get; } = new();

    internal ConcurrencyStats Concurrency { get; } = new();

    internal long SkippedForFile;

    internal TimeSpan Elapsed { get; set; }

    internal IEnumerable<FuzzFailure> Failures => Lanes.SelectMany(l => l.Failures);

    /// <summary>The report: lanes, engines, members.</summary>
    internal string Report()
    {
        var ci = CultureInfo.InvariantCulture;
        var sb = new StringBuilder();
        sb.AppendLine(ci, $"full differential fuzz: master seed {MasterSeed}, {Elapsed.TotalSeconds:F1} s, {Lanes.Sum(l => l.Cases)} cases, {Failures.Count()} failure(s)");
        foreach (FullFuzzLaneResult l in Lanes)
            sb.AppendLine(ci, $"lane {l.Name}: budget {l.Budget.TotalSeconds:F0} s, ran {l.Elapsed.TotalSeconds:F1} s, {l.Cases} cases ({l.Units} code units), {l.Skipped} skipped, {l.Failures.Count} failure(s); first {l.FirstCaseId}, last {l.LastCaseId}{(l.Extra.Length > 0 ? "; " + l.Extra : "")}");
        foreach (FullFuzzLaneResult l in Lanes.Where(l => l.PerFamily.Count > 0))
            sb.Append("families of lane ").Append(l.Name).Append(": ").AppendLine(string.Join(", ", l.PerFamily.OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => $"{kv.Key}={kv.Value}")));
        sb.AppendLine(ci, $"texts with a code unit above 255 (no file-backed engines): {SkippedForFile}");
        sb.AppendLine("engine\tchecks\tcheck_s\tchecks_per_s\tbuild_s\tmember_calls");
        for (int i = 0; i < Stats.Names.Count; i++)
        {
            string name = Stats.Names[i];
            long checks = Stats.Checks(name);
            double sec = Stats.CheckSeconds(name);
            long calls = i < Counters.Engines.Count ? Counters.Total(i) : 0;
            sb.AppendLine(ci, $"{name}\t{checks}\t{sec:F1}\t{(sec > 0 ? checks / sec : 0):F1}\t{Stats.BuildSeconds(name):F1}\t{calls}");
        }
        sb.AppendLine(ci, $"corruption: {Corruption.Cases} cases ({Corruption.Skipped} skipped), batteries {Corruption.BatteryRuns} ({Corruption.GuardedRuns} between guard pages, {Corruption.BatteryOk} operations answered, {Corruption.BatteryInvalid} refused as InvalidDataException); {Corruption}");
        sb.AppendLine(ci, $"concurrency: {Concurrency.Repetitions} repetitions, {Concurrency.ThreadRuns} thread runs, {Concurrency.Checksums} checksums equal to the single-thread run, {Concurrency.PackBuildsChecked} repetitions with the once-only pack build checked");
        sb.AppendLine("member\t" + string.Join("\t", Counters.Engines) + "\ttotal\tsteps");
        foreach (string key in FuzzMemberKeys.All)
            sb.AppendLine(ci, $"{key}\t{string.Join("\t", Enumerable.Range(0, Counters.Engines.Count).Select(e => Counters.Calls(e, key)))}\t{Counters.Calls(key)}\t{string.Join(";", Counters.StepsOf(key))}");
        return sb.ToString();
    }
}

/// <summary>
/// The budgeted run of the full differential fuzz (work package 1b.8): the matrix of generator families x alphabet sizes x size bands, every engine (in-memory trees, file-backed trees at the Header and Structure levels in the
/// Checked and Trusted modes), every member of <see cref="ISuffixTree"/> and of the contract interfaces, against the definitions, the Ukkonen tree, the v6 tree, then the container corruption fuzz and the shared-index test.
/// The time budget is split among the lanes by share; a case is a pure function of (master seed, lane, index), so a failure is reproduced by its id whatever the budget and the number of threads were.
/// </summary>
internal static class FullFuzzRunner
{
    /// <summary>Runs the lanes in sequence, each with its share of the budget.</summary>
    internal static FullFuzzSummary Run(FullFuzzOptions o)
    {
        var summary = new FullFuzzSummary
        {
            MasterSeed = o.MasterSeed,
            Counters = new MemberCounters(FuzzSubjects.EngineNames),
            Stats = new EngineStats(),
        };
        long start = Stopwatch.GetTimestamp();
        double total = o.Lanes.Sum(l => l.Share);
        var stop = new StopFlag(0);
        using var watchdog = new Watchdog(o.HangLimit, o.Log);
        foreach (FullFuzzLane lane in o.Lanes)
        {
            if (stop.Value != 0) break;
            TimeSpan budget = o.Budget * (lane.Share / total);
            var result = new FullFuzzLaneResult(lane.Name, budget);
            summary.Lanes.Add(result);
            long laneStart = Stopwatch.GetTimestamp();
            o.Log?.Invoke($"lane {lane.Name}: budget {budget.TotalSeconds:F0} s");
            switch (lane.Name)
            {
                case "tiny": Differential(o, FuzzTier.Tiny, result, summary, budget, watchdog, stop); break;
                case "small": Differential(o, FuzzTier.Small, result, summary, budget, watchdog, stop); break;
                case "large": Differential(o, FuzzTier.Large, result, summary, budget, watchdog, stop); break;
                case "xl": Differential(o, FullFuzzOptions.XlTier, result, summary, budget, watchdog, stop); break;
                case "corruption": Corruption(o, result, summary, budget, watchdog, stop); break;
                case "concurrency": Concurrency(o, result, summary, budget, watchdog, stop); break;
                default: throw new ArgumentException($"unknown lane '{lane.Name}'", nameof(o));
            }
            result.Elapsed = Stopwatch.GetElapsedTime(laneStart);
            summary.SkippedForFile += result.Name is "tiny" or "small" or "large" or "xl" ? result.Skipped : 0;
            o.Log?.Invoke($"lane {lane.Name} done: {result.Cases} cases in {result.Elapsed.TotalSeconds:F1} s, {result.Failures.Count} failure(s)");
        }
        summary.Elapsed = Stopwatch.GetElapsedTime(start);
        return summary;
    }

    private static string NewWorkDir()
    {
        string dir = Path.Combine(Path.GetTempPath(), "esa-fuzz-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        return dir;
    }

    private static void RunWorkers(int count, Action<int> body)
    {
        var errors = new ConcurrentBag<Exception>();
        Thread[] threads = [.. Enumerable.Range(0, count).Select(w => new Thread(() =>
        {
            try
            {
                body(w);
            }
            catch (Exception e)
            {
                errors.Add(e);
            }
        }))];
        // the workers run without the execution context of the caller: the checks report with Assert.Fail, which NUnit records in the result of the test that is current on the thread (and throws), so a failure that the
        // harness catches and minimizes as data would still fail the test that started the run; on a thread without the context the assertion is only an exception
        foreach (Thread t in threads)
        {
            using (ExecutionContext.SuppressFlow()) t.Start();
        }
        foreach (Thread t in threads) t.Join();
        if (!errors.IsEmpty) throw new AggregateException("a worker of the fuzz run failed outside a case", errors);
    }

    private static void Journal(FullFuzzOptions o, string lane, int worker, string text)
    {
        if (o.JournalDirectory is null) return;
        try
        {
            File.WriteAllText(Path.Combine(o.JournalDirectory, $"current-{lane}-{worker}.txt"), text);
        }
        catch (IOException)
        {
            // the journal is a help for the forensics of a crash, never a reason to fail a case
        }
    }

    private static void Differential(FullFuzzOptions o, FuzzTier tier, FullFuzzLaneResult result, FullFuzzSummary summary, TimeSpan budget, Watchdog watchdog, StopFlag stop)
    {
        long next = -1;
        long ticksPerCase = 0;
        long deadline = Stopwatch.GetTimestamp() + (long)(budget.TotalSeconds * Stopwatch.Frequency);
        long maxCases = o.MaxCasesPerLane ?? long.MaxValue;
        RunWorkers(o.Threads, w =>
        {
            string workDir = NewWorkDir();
            try
            {
                var ctx = new FullFuzzContext { Counters = summary.Counters, WorkDir = workDir, Stats = summary.Stats, CountFault = o.CountFault };
                while (stop.Value == 0)
                {
                    if (!MayStart(result, deadline, ticksPerCase, o.Threads)) break;
                    long i = Interlocked.Increment(ref next);
                    if (i >= maxCases) break;
                    FuzzCase c = FuzzRunner.CaseAt(o.MasterSeed, tier, i);
                    string text = c.Text();
                    if (i == 0) result.FirstCaseId = c.Id;
                    result.LastCaseId = c.Id;
                    Journal(o, tier.Name, w, c.Id);
                    watchdog.Enter(w, $"{tier.Name} {c.Id}");
                    long caseStart = Stopwatch.GetTimestamp();
                    long before = Interlocked.Read(ref ctx.SkippedForFile);
                    FuzzMismatch? mismatch;
                    try
                    {
                        mismatch = FullFuzzCheck.CheckText(text, FullFuzzConfig.Of(c), ctx);
                    }
                    catch (Exception e)
                    {
                        mismatch = new FuzzMismatch("harness:" + e.GetType().Name, e.ToString());
                    }
                    watchdog.Leave(w);
                    Interlocked.Add(ref ticksPerCase, Stopwatch.GetTimestamp() - caseStart);
                    Interlocked.Increment(ref result.Cases);
                    Interlocked.Add(ref result.Units, text.Length);
                    if (Interlocked.Read(ref ctx.SkippedForFile) > before) Interlocked.Increment(ref result.Skipped);
                    result.PerFamily.AddOrUpdate(c.Family, 1, (_, v) => v + 1);
                    if (mismatch is { } m) Fail(o, c, m, text, FullFuzzConfig.Of(c), workDir, result, stop);
                }
            }
            finally
            {
                Directory.Delete(workDir, true);
            }
        });
    }

    /// <summary>A worker starts a case only if the lane's time is not up and the case, at the mean cost of the cases so far, is likely to end before it (a lane of long cases does not run past its budget by a whole case).</summary>
    private static bool MayStart(FullFuzzLaneResult result, long deadline, long totalTicks, int workers)
    {
        long now = Stopwatch.GetTimestamp();
        if (now >= deadline) return false;
        long cases = Interlocked.Read(ref result.Cases);
        if (cases < workers) return true;
        long mean = totalTicks / cases;
        return now + (mean * 7 / 10) < deadline;
    }

    internal sealed class StopFlag(int initial)
    {
        private int _value = initial;

        internal int Value => Volatile.Read(ref _value);

        internal void Set() => Interlocked.Exchange(ref _value, 1);
    }

    private static void Fail(FullFuzzOptions o, FuzzCase c, FuzzMismatch m, string text, FullFuzzConfig cfg, string workDir, FullFuzzLaneResult result, StopFlag stop)
    {
        long t0 = Stopwatch.GetTimestamp();
        var scratch = new FullFuzzContext { Counters = new MemberCounters(FuzzSubjects.EngineNames), WorkDir = workDir, CountFault = o.CountFault };
        string minimized = FuzzMinimizer.Minimize(
            text,
            t => Stopwatch.GetElapsedTime(t0) < o.MinimizeLimit && FullFuzzCheck.CheckText(t, cfg, scratch) is { } x && x.Op == m.Op,
            maxChecks: 5_000);
        var failure = new FuzzFailure(c, m.Op, m.Detail, minimized);
        result.Failures.Add(failure);
        o.Log?.Invoke($"FAIL {c.Id} op {m.Op}: {Short(m.Detail)}; minimized to {minimized.Length} unit(s): {FuzzRegressions.EncodeUnits(minimized)}");
        if (result.Failures.Count >= o.MaxFailures) stop.Set();
    }

    /// <summary>True when the corruption case still fails on the shrunk text; an exception of the harness itself on a shrunk text is not the failure that is being minimized.</summary>
    private static bool FailsToo(Func<string?> run)
    {
        try
        {
            return run() is not null;
        }
        catch (Exception)
        {
            return false;
        }
    }

    private static string Short(string s) => s.Length <= 400 ? s : s[..400] + "...";

    private static void Corruption(FullFuzzOptions o, FullFuzzLaneResult result, FullFuzzSummary summary, TimeSpan budget, Watchdog watchdog, StopFlag stop)
    {
        long next = -1;
        long deadline = Stopwatch.GetTimestamp() + (long)(budget.TotalSeconds * Stopwatch.Frequency);
        long maxCases = o.MaxCasesPerLane ?? long.MaxValue;
        RunWorkers(o.Threads, w =>
        {
            string workDir = NewWorkDir();
            try
            {
                var ctx = new CorruptionContext { Counters = summary.Counters, WorkDir = workDir, Stats = summary.Corruption };
                while (stop.Value == 0)
                {
                    long i = Interlocked.Increment(ref next);
                    if (i >= maxCases || Stopwatch.GetTimestamp() >= deadline) break;
                    (FuzzCase c, DamageKind kind, ulong damageSeed) = CorruptionFuzz.Draw(o.MasterSeed, i);
                    string text = c.Text();
                    if (i == 0) result.FirstCaseId = c.Id + "/" + kind;
                    result.LastCaseId = c.Id + "/" + kind;
                    Journal(o, "corruption", w, $"{i} {c.Id} {kind} {damageSeed:x16}");
                    watchdog.Enter(w, $"corruption {i} {c.Id} {kind} {damageSeed:x16}");
                    string? failure;
                    try
                    {
                        failure = CorruptionFuzz.Run(c.Id, text, kind, damageSeed, ctx);
                    }
                    catch (Exception e)
                    {
                        failure = $"{c.Id} {kind}: harness exception {e.GetType().Name}: {e}";
                    }
                    watchdog.Leave(w);
                    Interlocked.Increment(ref result.Cases);
                    Interlocked.Add(ref result.Units, text.Length);
                    result.PerFamily.AddOrUpdate(c.Family, 1, (_, v) => v + 1);
                    if (failure is not null)
                    {
                        long t0 = Stopwatch.GetTimestamp();
                        var scratch = new CorruptionContext { Counters = new MemberCounters(FuzzSubjects.EngineNames), WorkDir = workDir };
                        string minimized = FuzzMinimizer.Minimize(
                            text,
                            t => Stopwatch.GetElapsedTime(t0) < o.MinimizeLimit && FailsToo(() => CorruptionFuzz.Run(c.Id, t, kind, damageSeed, scratch)),
                            maxChecks: 2_000);
                        result.Failures.Add(new FuzzFailure(c, CorruptionFuzz.Op(kind, damageSeed), failure, minimized));
                        o.Log?.Invoke($"FAIL corruption {c.Id} {kind}: {Short(failure)}; minimized to {minimized.Length} unit(s): {FuzzRegressions.EncodeUnits(minimized)}");
                        if (result.Failures.Count >= o.MaxFailures) stop.Set();
                    }
                }
            }
            finally
            {
                Directory.Delete(workDir, true);
            }
        });
        result.Skipped = summary.Corruption.Skipped;
    }

    private static void Concurrency(FullFuzzOptions o, FullFuzzLaneResult result, FullFuzzSummary summary, TimeSpan budget, Watchdog watchdog, StopFlag stop)
    {
        long next = -1;
        long deadline = Stopwatch.GetTimestamp() + (long)(budget.TotalSeconds * Stopwatch.Frequency);
        long maxCases = o.MaxCasesPerLane ?? long.MaxValue;
        RunWorkers(1, w =>
        {
            string workDir = NewWorkDir();
            try
            {
                var ctx = new FullFuzzContext { Counters = summary.Counters, WorkDir = workDir, Stats = summary.Stats };
                while (stop.Value == 0)
                {
                    long i = Interlocked.Increment(ref next);
                    if (i >= maxCases || Stopwatch.GetTimestamp() >= deadline) break;
                    FuzzCase c = FuzzRunner.CaseAt(o.MasterSeed, FullFuzzOptions.ConcurrencyTier, i);
                    string text = c.Text();
                    if (i == 0) result.FirstCaseId = c.Id;
                    result.LastCaseId = c.Id;
                    Journal(o, "concurrency", w, c.Id);
                    watchdog.Enter(w, $"concurrency {c.Id}");
                    string? failure;
                    try
                    {
                        failure = FullFuzzConcurrency.Run(text, c.Seed, ctx, summary.Concurrency, o.ConcurrencyRepetitions, o.ConcurrencyThreads);
                    }
                    catch (Exception e)
                    {
                        failure = $"harness exception {e.GetType().Name}: {e}";
                    }
                    watchdog.Leave(w);
                    Interlocked.Increment(ref result.Cases);
                    Interlocked.Add(ref result.Units, text.Length);
                    result.PerFamily.AddOrUpdate(c.Family, 1, (_, v) => v + 1);
                    if (failure is not null)
                    {
                        result.Failures.Add(new FuzzFailure(c, "concurrency", failure, text));
                        o.Log?.Invoke($"FAIL concurrency {c.Id}: {Short(failure)}");
                        if (result.Failures.Count >= o.MaxFailures) stop.Set();
                    }
                }
            }
            finally
            {
                Directory.Delete(workDir, true);
            }
        });
        result.Extra = $"{summary.Concurrency.Repetitions} repetitions of {o.ConcurrencyThreads} threads";
    }
}

/// <summary>Ends the process when a worker has been in one case longer than the limit (a hang cannot be recovered from inside the process): the case is written to the log and to the standard error first.</summary>
internal sealed class Watchdog : IDisposable
{
    private readonly TimeSpan? _limit;
    private readonly Action<string>? _log;
    private readonly ConcurrentDictionary<int, (long Since, string Case)> _active = new();
    private readonly Timer? _timer;

    internal Watchdog(TimeSpan? limit, Action<string>? log)
    {
        _limit = limit;
        _log = log;
        if (limit is not null) _timer = new Timer(_ => Check(), null, TimeSpan.FromSeconds(2), TimeSpan.FromSeconds(2));
    }

    internal void Enter(int worker, string caseText)
    {
        if (_limit is not null) _active[worker] = (Stopwatch.GetTimestamp(), caseText);
    }

    internal void Leave(int worker)
    {
        if (_limit is not null) _active.TryRemove(worker, out _);
    }

    private void Check()
    {
        foreach (var (worker, (since, caseText)) in _active)
        {
            if (Stopwatch.GetElapsedTime(since) <= _limit) continue;
            string message = $"HANG: worker {worker} has been in the case '{caseText}' for more than {_limit!.Value.TotalSeconds:F0} s";
            _log?.Invoke(message);
            Console.Error.WriteLine(message);
            Environment.FailFast(message);
        }
    }

    public void Dispose() => _timer?.Dispose();
}
