using System.Globalization;
using SuffixTree.Esa.Tests.Facade;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// The budgeted run of the full differential fuzz (work package 1b.8, gate G2): the whole matrix (the families x alphabets x size bands over the in-memory and the file-backed engines, the container corruption fuzz, the
/// shared-index test) for a time budget from the environment. Category <c>Fuzz</c>, <c>[Explicit]</c>. Environment: <c>ESA_FUZZ_SEED</c> (master seed, default 1), <c>ESA_FUZZ_BUDGET_S</c> (seconds, default 60),
/// <c>ESA_FUZZ_THREADS</c> (workers per lane, default the number of cores), <c>ESA_FUZZ_LANES</c> (a comma separated subset of tiny, small, large, xl, corruption, concurrency; default all with the shares of
/// <see cref="FullFuzzOptions.DefaultLanes"/>), <c>ESA_FUZZ_REPORT</c> (a file for the report), <c>ESA_FUZZ_JOURNAL</c> (a directory for the journal of the case each worker is in), <c>ESA_FUZZ_HANG_S</c> (seconds a worker may stay in one case, default 300), <c>ESA_FUZZ_MAXFAIL</c> (stop after this many failures, default 3), <c>ESA_FUZZ_MINIMIZE_S</c> (seconds for the minimizer per failure, default 30).
/// A failure prints <c>family:seed:n:sigma:param</c>, the failed operation and the minimized text, writes the regression lines to <c>fuzz-failures.txt</c> in the test work directory, and the lines are appended to
/// <c>Fuzz/regressions.txt</c> once the defect is fixed. <c>analysis/em/run_1b8_fuzz.sh</c> drives the run.
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Full differential fuzz of work package 1b.8 with a time budget; run with --filter FullyQualifiedName~FullDifferentialFuzzTests (see analysis/em/run_1b8_fuzz.sh)")]
public class FullDifferentialFuzzTests
{
    private static string? Env(string name) => Environment.GetEnvironmentVariable(name);

    [Test]
    public void BudgetedRun_AllEngines_AllMembers_Corruption_Concurrency_HasNoFailure()
    {
        ulong seed = ulong.Parse(Env("ESA_FUZZ_SEED") ?? "1", CultureInfo.InvariantCulture);
        double budget = double.Parse(Env("ESA_FUZZ_BUDGET_S") ?? "60", CultureInfo.InvariantCulture);
        int threads = int.Parse(Env("ESA_FUZZ_THREADS") ?? Environment.ProcessorCount.ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
        double hang = double.Parse(Env("ESA_FUZZ_HANG_S") ?? "300", CultureInfo.InvariantCulture);
        string[]? names = Env("ESA_FUZZ_LANES")?.Split(',', StringSplitOptions.RemoveEmptyEntries);
        FullFuzzLane[] lanes = names is null ? [.. FullFuzzOptions.DefaultLanes] : [.. FullFuzzOptions.DefaultLanes.Where(l => names.Contains(l.Name))];
        int maxFailures = int.Parse(Env("ESA_FUZZ_MAXFAIL") ?? "3", CultureInfo.InvariantCulture);
        double minimize = double.Parse(Env("ESA_FUZZ_MINIMIZE_S") ?? "30", CultureInfo.InvariantCulture);
        string? journal = Env("ESA_FUZZ_JOURNAL");
        if (journal is not null) Directory.CreateDirectory(journal);
        string? reportFile = Env("ESA_FUZZ_REPORT");
        void Log(string line)
        {
            TestContext.Progress.WriteLine(line);
            if (reportFile is not null) File.AppendAllText(reportFile + ".log", line + "\n");
        }

        FullFuzzSummary summary = FullFuzzRunner.Run(new FullFuzzOptions
        {
            MasterSeed = seed,
            Budget = TimeSpan.FromSeconds(budget),
            Threads = threads,
            Lanes = lanes,
            HangLimit = TimeSpan.FromSeconds(hang),
            MaxFailures = maxFailures,
            MinimizeLimit = TimeSpan.FromSeconds(minimize),
            JournalDirectory = journal,
            Log = Log,
        });
        string report = summary.Report();
        TestContext.Out.WriteLine(report);
        if (reportFile is not null) File.WriteAllText(reportFile, report);
        Log($"done: {summary.Lanes.Sum(l => l.Cases)} cases in {summary.Elapsed.TotalSeconds:F0} s, {summary.Failures.Count()} failure(s)");

        FuzzFailure[] failures = [.. summary.Failures];
        if (failures.Length > 0)
        {
            string lines = string.Join("\n", failures.Select(FuzzRegressions.Format));
            File.WriteAllText(Path.Combine(TestContext.CurrentContext.WorkDirectory, "fuzz-failures.txt"), lines + "\n");
            Assert.Fail($"{failures.Length} failure(s); regression lines:\n{lines}\nfirst: {failures[0].Op}: {failures[0].Detail}");
        }
    }
}

/// <summary>The harness of work package 1b.8 on a small slice (the default run): the coverage of the members, the fault injection on the file-backed path, the corruption lane, the regression corpus.</summary>
[TestFixture]
public class FullFuzzHarnessTests
{
    /// <summary>A slice of every lane: enough cases that every member of every interface is called on every kind of engine.</summary>
    internal static FullFuzzSummary Smoke(ulong seed = 1, bool large = false)
    {
        FullFuzzLane[] lanes = large
            ? [new FullFuzzLane("tiny", 1), new FullFuzzLane("small", 1), new FullFuzzLane("large", 1), new FullFuzzLane("corruption", 1), new FullFuzzLane("concurrency", 1)]
            : [new FullFuzzLane("tiny", 1), new FullFuzzLane("small", 1), new FullFuzzLane("corruption", 1), new FullFuzzLane("concurrency", 1)];
        return FullFuzzRunner.Run(new FullFuzzOptions
        {
            MasterSeed = seed,
            Budget = TimeSpan.FromMinutes(10),
            Threads = 2,
            Lanes = lanes,
            MaxCasesPerLane = 40,
            ConcurrencyRepetitions = 2,
        });
    }

    private static readonly Lazy<FullFuzzSummary> SmokeRun = new(() => Smoke());

    [Test]
    public void TheSmokeSlice_HasNoFailure_AndRunsEveryLane()
    {
        FullFuzzSummary s = SmokeRun.Value;
        Assert.Multiple(() =>
        {
            Assert.That(s.Failures.Select(f => $"{f.Case.Id} {f.Op}: {f.Detail}"), Is.Empty);
            Assert.That(s.Lanes.Select(l => l.Name), Is.EqualTo(new[] { "tiny", "small", "corruption", "concurrency" }));
            Assert.That(s.Lanes.All(l => l.Cases > 0), Is.True, "every lane ran cases");
            Assert.That(s.Corruption.Cases, Is.GreaterThan(0));
            Assert.That(s.Concurrency.Repetitions, Is.GreaterThan(0));
            Assert.That(s.Stats.Names.All(name => s.Stats.Checks(name) > 0), Is.True, "every engine, the Ukkonen and the v6 tree included, was checked");
        });
    }

    [Test]
    public void TheHarness_CallsEveryMember_OnTheInMemoryAndOnTheFileBackedTree()
    {
        MemberCounters c = SmokeRun.Value.Counters;
        var uncalled = new List<string>();
        foreach (string key in FuzzMemberKeys.All)
        {
            bool memoryOnly = key.StartsWith("EsaSuffixTreeFactory.Create", StringComparison.Ordinal);
            bool fileOnly = key.Contains(".Open(", StringComparison.Ordinal) || key.StartsWith("EsaOpenOptions.", StringComparison.Ordinal);
            long memory = 0, file = 0;
            for (int e = 0; e < c.Engines.Count; e++)
            {
                if (c.Engines[e].StartsWith("memory", StringComparison.Ordinal)) memory += c.Calls(e, key);
                else file += c.Calls(e, key);
            }
            if ((!fileOnly && memory == 0) || (!memoryOnly && file == 0)) uncalled.Add($"{key}: in-memory {memory}, file-backed {file}");
        }
        Assert.That(uncalled, Is.Empty, "members the harness did not call on every kind of engine");
    }

    [Test]
    public void EveryKindOfDamage_IsDrawn_AndEveryOutcomeWasAllowed()
    {
        CorruptionStats stats = SmokeRun.Value.Corruption;
        Assert.Multiple(() =>
        {
            foreach (DamageKind kind in Enum.GetValues<DamageKind>()) Assert.That(stats.Total(kind), Is.GreaterThan(0), $"{kind} was drawn");
            Assert.That(stats.BatteryRuns, Is.GreaterThan(0), "a damaged file that passed the Header level was read");
            Assert.That(stats.BatteryOk, Is.GreaterThan(0));
        });
    }

    [Test]
    public void AFaultOfTheFileBackedPath_IsFound_AndMinimizedToAShortText()
    {
        // the file-backed tree at the Structure level reports one occurrence too many for a pattern that contains "aba"
        FullFuzzSummary s = FullFuzzRunner.Run(new FullFuzzOptions
        {
            MasterSeed = 5,
            Budget = TimeSpan.FromMinutes(5),
            Threads = 1,
            Lanes = [new FullFuzzLane("tiny", 1)],
            MaxCasesPerLane = 400,
            MaxFailures = 1,
            CountFault = (engine, pattern, count) => engine == "file-structure-trusted" && count > 0 && pattern.Contains("aba", StringComparison.Ordinal) ? count + 1 : count,
        });
        FuzzFailure[] failures = [.. s.Failures];
        Assert.That(failures, Has.Length.EqualTo(1), "the fault was found");
        Assert.Multiple(() =>
        {
            Assert.That(failures[0].Op, Does.StartWith("file-structure-trusted:"), failures[0].Detail);
            Assert.That(failures[0].MinimizedText, Is.EqualTo("aba"), "the minimizer shrinks the text of the case to the three symbols of the pattern");
        });
    }

    [Test]
    public void AFaultOfTheInMemoryPath_IsFound_AndMinimizedToAShortText()
    {
        FullFuzzSummary s = FullFuzzRunner.Run(new FullFuzzOptions
        {
            MasterSeed = 6,
            Budget = TimeSpan.FromMinutes(5),
            Threads = 1,
            Lanes = [new FullFuzzLane("tiny", 1)],
            MaxCasesPerLane = 400,
            MaxFailures = 1,
            CountFault = (engine, pattern, count) => engine == "memory-default" && count > 1 && pattern.Length == 2 ? count - 1 : count,
        });
        FuzzFailure[] failures = [.. s.Failures];
        Assert.That(failures, Has.Length.EqualTo(1), "the fault was found");
        Assert.Multiple(() =>
        {
            Assert.That(failures[0].Op, Does.StartWith("memory-default:"), failures[0].Detail);
            Assert.That(failures[0].MinimizedText.Length, Is.LessThanOrEqualTo(3), "a pattern of two symbols that occurs twice needs at most three symbols");
        });
    }

    [Test]
    public void AFaultOfASingleCase_IsReproducedFromItsId_Alone()
    {
        FuzzCase c = FuzzRunner.CaseAt(5, FuzzTier.Tiny, 3);
        using var dir = new External.TempDir();
        var ctx = new FullFuzzContext { Counters = new MemberCounters(FuzzSubjects.EngineNames), WorkDir = dir.Path };
        Assert.That(FullFuzzCheck.CheckText(c.Text(), FullFuzzConfig.Of(c), ctx), Is.Null, c.Id);
    }

    [Test]
    public void TheDigest_OfEveryEngine_IsTheSame_OnNamedTexts()
    {
        using var dir = new External.TempDir();
        var ctx = new FullFuzzContext { Counters = new MemberCounters(FuzzSubjects.EngineNames), WorkDir = dir.Path };
        string[] texts = ["", "a", "aa", "ab", "banana", "mississippi", "abracadabra", new string('a', 130), StringFamilies.Fibonacci(200), StringFamilies.RandomDna(400, 3), "ab\u0001ab\0￿Ā"];
        foreach (string text in texts)
        {
            FuzzMismatch? m = FullFuzzCheck.CheckText(text, FullFuzzConfig.Of((ulong)text.Length * 7919), ctx);
            Assert.That(m, Is.Null, $"'{FacadeChecks.Short(text)}': {m?.Op}: {m?.Detail}");
        }
    }

    [Test]
    public void TheCommittedRegressionCorpus_ReplaysThroughTheFullHarness()
    {
        using var dir = new External.TempDir();
        var ctx = new FullFuzzContext { Counters = new MemberCounters(FuzzSubjects.EngineNames), WorkDir = dir.Path };
        var cctx = new CorruptionContext { Counters = ctx.Counters, WorkDir = dir.Path };
        foreach (FuzzRegression r in FuzzRegressions.Load())
        {
            foreach (string text in new[] { r.Case.Text(), r.MinimizedText })
            {
                if (CorruptionFuzz.ParseOp(r.Op) is { } damage)
                {
                    Assert.That(CorruptionFuzz.Run(r.Case.Id, text, damage.Kind, damage.DamageSeed, cctx), Is.Null, $"{r.Case.Id} {r.Op}");
                }
                else
                {
                    FuzzMismatch? m = FullFuzzCheck.CheckText(text, FullFuzzConfig.Of(r.Case), ctx);
                    Assert.That(m, Is.Null, $"{r.Case.Id} {r.Op}: {m?.Op}: {m?.Detail}");
                }
            }
        }
    }

    [Test]
    public void OpOf_NamesTheMemberOfAFailureMessage()
    {
        Assert.Multiple(() =>
        {
            Assert.That(FullFuzzCheck.OpOf("[file-header-checked] n=9 text 'abababab': pattern 'ab': CountOccurrences 3 against 2"), Is.EqualTo("CountOccurrences"));
            Assert.That(FullFuzzCheck.OpOf("[x] n=9 text 'q': query 'abc' L 3: MEM 1:2:3 against the definition"), Is.EqualTo("MEM"));
            Assert.That(FullFuzzCheck.OpOf("something else"), Is.EqualTo("assert"));
        });
    }

    [Test]
    public void TheMemberKeys_AreUnique_AndTheCountingTreeKnowsEveryKeyItUses()
    {
        Assert.That(FuzzMemberKeys.All.Distinct().Count(), Is.EqualTo(FuzzMemberKeys.All.Count));
        Assert.That(() => FuzzMemberKeys.IndexOf("ISuffixTreeSearch.NoSuchMember"), Throws.ArgumentException);
        Assert.That(FuzzMemberKeys.All, Does.Contain("ISuffixTreeSearch.Contains(ReadOnlySpan<Char>)").And.Contain("ISuffixIndex.Locate(IndexLocus,Int64,Span<Int64>)").And.Contain("SuffixIndexFactory.Open(String,EsaOpenOptions)"));
    }
}

/// <summary>The committed operation coverage list: one row per member (a reflection test), and the explicit generator.</summary>
[TestFixture]
public class FullFuzzCoverageTableTests
{
    [Test]
    public void EveryMemberOfTheAuditedInterfaces_HasARowInTheCommittedTable_AndEveryRowNamesAMember()
    {
        string table = File.ReadAllText(FuzzCoverageTable.DefaultPath());
        Assert.That(FuzzCoverageTable.Problems(table), Is.Empty);
    }

    [Test]
    public void ARemovedRow_IsReported_ASurplusRowToo()
    {
        string table = File.ReadAllText(FuzzCoverageTable.DefaultPath());
        string[] lines = table.Split('\n');
        int victim = Array.FindIndex(lines, l => l.StartsWith("ISuffixTreeAnalysis.FindMaximalUniqueMatches(", StringComparison.Ordinal));
        Assert.That(victim, Is.GreaterThanOrEqualTo(0));
        string without = string.Join("\n", lines.Where((_, i) => i != victim));
        Assert.That(FuzzCoverageTable.Problems(without), Is.EqualTo(new[] { "no row for the member ISuffixTreeAnalysis.FindMaximalUniqueMatches(String,Int32,MumUniqueness)" }));
        Assert.That(FuzzCoverageTable.Problems(table + "ISuffixTreeSearch.Removed\t0\t0\tFuzzDigest.Compute\n"), Is.EqualTo(new[] { "the row ISuffixTreeSearch.Removed names no member of the audited types" }));
        Assert.That(FuzzCoverageTable.Problems(table + lines[victim] + "\n"), Is.EqualTo(new[] { "duplicate row for ISuffixTreeAnalysis.FindMaximalUniqueMatches(String,Int32,MumUniqueness)" }));
    }

    [Test]
    public void EveryInterfaceMember_IsInTheKeyListOfTheTypes_IncludingTheInheritedDefaults()
    {
        // the four interfaces of the brief and their members by name: a rename or an addition changes the list, which must then be reflected in the table
        string[] names = ["Text", "NodeCount", "LeafCount", "MaxDepth", "IsEmpty", "Contains", "FindAllOccurrences", "CountOccurrences", "LongestRepeatedSubstring", "LongestRepeatedSubstringMemory", "GetAllSuffixes", "EnumerateSuffixes",
            "LongestCommonSubstring", "LongestCommonSubstringInfo", "FindAllLongestCommonSubstrings", "FindAllDistinctLongestCommonSubstrings", "FindExactMatchAnchors", "FindMaximalExactMatches", "FindMaximalUniqueMatches",
            "PrintTree", "Traverse", "CountDistinctSubstringsByLength", "CountDistinctSubstrings", "FindAllLongestRepeatedSubstrings", "FindMaximalRepeatedPairs"];
        foreach (string name in names)
            Assert.That(FuzzMemberKeys.All.Any(k => k.Split('.')[1].StartsWith(name, StringComparison.Ordinal)), Is.True, name);
    }

    [Test]
    [Explicit("regenerates analysis/s1b/harness_member_coverage.tsv from a smoke slice of the harness; run with --filter FullyQualifiedName~RegenerateTheCommittedTable")]
    public void RegenerateTheCommittedTable()
    {
        FullFuzzSummary s = FullFuzzHarnessTests.Smoke(seed: 3);
        Assert.That(s.Failures, Is.Empty);
        File.WriteAllText(FuzzCoverageTable.DefaultPath(), FuzzCoverageTable.Render(s.Counters));
        Assert.That(FuzzCoverageTable.Problems(File.ReadAllText(FuzzCoverageTable.DefaultPath())), Is.Empty);
    }
}
