using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-12b, B09 F311): the structure methods of T-Coffee 13.41 that call SAP / DaliLite, against the runs 1, 3, 4, 5 and 6 of the user's Docker package <c>handoff/A2-12b-dali-sap-stderr</c>
/// (TestData/Alignment/tcoffee_sap_dali_runs_oracle.tsv.gz, scripts/alignment/convert_tcoffee_sap_dali_returned.py): the same 130 <c>sap_real</c> sets, 8 <c>-mode 3dcoffee</c> sets and 60 <c>dali_pair</c> sets in every run, with t_coffee's standard error.
/// SAP and the structure aligners are NOT deterministic in the user's environment (the same sets and files, pairs of structures without a result in one run and with an identical one in another, a process killed by the kernel in some runs), so each row is
/// checked against the port FED WITH ITS OWN RUN'S SAP OUTPUT: t_coffee's side of the interface is a function of what the external program wrote.
/// </summary>
[TestFixture]
public class SequenceAligner_TCoffeeSapDaliRuns_Tests
{
    private sealed record Row(string Run, string Section, int Index, string Flags, string[] Seqs, string Rows, JsonElement Extra)
    {
        public bool Failed => Rows.StartsWith("FAIL:", StringComparison.Ordinal);
        public string?[] Pdb => Extra.GetProperty("pdb").EnumerateArray().Select(e => e.GetString() is { Length: > 0 } s ? s : null).ToArray();
        public string[] Err => Extra.GetProperty("err").EnumerateArray().Select(e => e.GetString()!).ToArray();
        public bool Killed => Rows == "FAIL:-9" || Err.Any(l => l.Contains("Killed", StringComparison.Ordinal));
    }

    private static readonly Lazy<List<Row>> All = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.tcoffee_sap_dali_runs_oracle.tsv";
        using var stream = TestDataResource.Open(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], int.Parse(f[2]), f[3], f[4].Split('|'), f[5], JsonDocument.Parse(f[6]).RootElement.Clone()));
        }

        return list;
    });

    private static IEnumerable<Row> Section(string section) => All.Value.Where(r => r.Section == section);

    private static Dictionary<(int, int), string> SapFiles(Row r)
    {
        var map = new Dictionary<(int, int), string>();
        foreach (var f in r.Extra.GetProperty("sap").EnumerateObject())
        {
            var m = System.Text.RegularExpressions.Regex.Match(f.Name, @"^p(\d+)\.pdb_p(\d+)\.pdb\.sap_results$");
            Assert.That(m.Success, Is.True, f.Name);
            map[(int.Parse(m.Groups[1].Value), int.Parse(m.Groups[2].Value))] = f.Value.GetString()!;
        }

        return map;
    }

    private static TCoffeeProgramOptions Options(Row r, TCoffeePairMethod? installed, IReadOnlyDictionary<(int First, int Second), string> outputs)
    {
        var map = new Dictionary<TCoffeePairMethod, IReadOnlyDictionary<(int First, int Second), string>>();
        if (installed is { } m) map[m] = outputs;
        if (r.Extra.TryGetProperty("tm", out var tm))
        {
            var byPair = new Dictionary<(int, int), string>();
            foreach (var p in tm.EnumerateObject())
            {
                var ij = p.Name.Split(',');
                byPair[(int.Parse(ij[0]), int.Parse(ij[1]))] = p.Value.GetString()!;
            }

            map[TCoffeePairMethod.TMalignPair] = byPair;
        }

        return TCoffeeProgramOptions.FromCommandLine(r.Flags) with { StructureTemplates = r.Pdb, StructureAlignerOutputs = map.Count > 0 ? map : null };
    }

    private static string Run(Row r, TCoffeeProgramOptions o) => string.Join('|', SequenceAligner.MultipleAlignTCoffeeProgram(r.Seqs, o).AlignedSequences);

    private static string Describe(Row r) => $"run {r.Run} {r.Section} #{r.Index} '{r.Flags}' {string.Join('|', r.Seqs)}";

    // ---- s1: the real SAP -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// <summary>
    /// 5 x 130 sets (runs 1, 3, 4, 5, 6): every alignment t_coffee wrote equals the port fed with the SAP files of that run; every stop (exit status 1) is the port's <see cref="TCoffeeAbortedException"/>, and t_coffee's standard error of a stop is
    /// the one message "You have not provided any method or enough Sequences [FATAL]" (t_coffee.c L5005: the library has no entry; the warning "SAP failed to align" of util_dp_drivers.c L2146 is suppressed by -no_warning).
    /// </summary>
    [Test]
    public void RealSap_EveryRowOfEveryRunEqualsThePortFedWithItsOwnSapOutput()
    {
        var rows = Section("s1").ToList();
        Assert.That(rows, Has.Count.EqualTo(650));
        var failures = new List<string>();
        int aborted = 0, compared = 0, killed = 0;
        foreach (var r in rows)
        {
            if (r.Killed)
            {
                killed++;
                continue;
            }

            var options = Options(r, TCoffeePairMethod.SapPair, SapFiles(r));
            if (r.Failed)
            {
                aborted++;
                Assert.That(r.Rows, Is.EqualTo("FAIL:1"), Describe(r));
                Assert.That(r.Err, Does.Contain("You have not provided any method or enough Sequences[FATAL]"), Describe(r));
                var ex = Assert.Throws<TCoffeeAbortedException>(() => Run(r, options), Describe(r));
                Assert.That(ex!.Message, Does.Contain("You have not provided any method or enough Sequences"), Describe(r));
                continue;
            }

            compared++;
            string got = Run(r, options);
            if (got != r.Rows) failures.Add($"{Describe(r)}: got {got} expected {r.Rows}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count} sets differ; first: {failures.FirstOrDefault()}");
        Assert.That(aborted, Is.GreaterThan(40));
        Assert.That(compared + aborted + killed, Is.EqualTo(650));
        TestContext.Out.WriteLine($"s1: {rows.Count} sets, {compared} alignments equal, {aborted} stops reproduced, {killed} refused (killed process or not reproducible)");
    }

    /// <summary>
    /// SAP's result for a pair is NOT deterministic in the user's environment: with identical sets and identical files a pair has a result file in one run and none in another (the files that exist are identical after the process id of the
    /// temporary names is removed). The evidence of the nondeterminism, not a defect of the port.
    /// </summary>
    [Test]
    public void RealSap_ResultFilesAreIdenticalWhenPresentButTheirPresenceVariesBetweenRuns()
    {
        var byKey = new Dictionary<(int, string), List<(string Run, string Text)>>();
        foreach (var r in Section("s1"))
            foreach (var f in SapFiles(r))
            {
                var key = (r.Index, $"{f.Key.Item1},{f.Key.Item2}");
                if (!byKey.TryGetValue(key, out var l)) byKey[key] = l = [];
                l.Add((r.Run, System.Text.RegularExpressions.Regex.Replace(f.Value, @"\d+_[12]\.\d+\.sap_tmp", "TMP")));
            }

        Assert.That(byKey.Values.Where(l => l.Select(x => x.Text).Distinct().Count() > 1), Is.Empty, "a SAP file that differs between runs");
        int partial = byKey.Values.Count(l => l.Count < 5);
        Assert.That(partial, Is.GreaterThan(50), "pairs with a SAP result in some runs only");
        TestContext.Out.WriteLine($"s1: {byKey.Count} pairs with a result, {partial} of them in fewer than the 5 runs; the files that exist are identical");
    }

    // ---- s2: -mode 3dcoffee with the real SAP ----------------------------------------------------------------------------------------------------------------------------------------------------

    [Test]
    public void ThreeDCoffee_RowsAreTheSameInEveryRunAndNoStructureIsFound()
    {
        var rows = Section("s2").ToList();
        Assert.That(rows, Has.Count.EqualTo(40));
        foreach (var g in rows.GroupBy(r => r.Index))
            Assert.That(g.Select(r => r.Rows).Distinct().Count(), Is.EqualTo(1), $"-mode 3dcoffee set {g.Key}: the five runs differ");
        // SELF_P_: the template of a sequence is looked up by the sequence NAME (q000, ...), which is no structure: the program warns and every pair is aligned by proba_pair
        foreach (var r in rows)
        {
            Assert.That(r.Failed, Is.False, Describe(r));
            Assert.That(r.Err.Any(l => l.Contains("_P_ Template | q000 | Could Not Be Found", StringComparison.Ordinal)), Is.True, Describe(r));
        }
    }

    /// <summary>
    /// The sets of the 60-set dali_pair plan (the same sets in s3 and s4 of every run) whose outcome is not reproducible in the user's environment, with the evidence: the TM-align pairs that wrote an alignment differ between runs for identical files
    /// (the generator's own TM-align failed in one run and not in another, so the TMalign_pair part of the port's input is not t_coffee's), or t_coffee itself stopped in some runs and ran through in others (the glue's lock files between the processes of the method
    /// are timing dependent). Such a set is refused in every run.
    /// </summary>
    private static readonly Lazy<Dictionary<int, string>> Unreproducible = new(() =>
    {
        var refused = new Dictionary<int, string>();
        foreach (var g in All.Value.Where(r => r.Section is "s3" or "s4").GroupBy(r => r.Index))
        {
            Assert.That(g.Select(r => string.Join('|', r.Seqs)).Distinct().Count(), Is.EqualTo(1), $"set {g.Key} is the same in every run");
            Assert.That(g.Select(r => string.Join(',', r.Extra.GetProperty("pdb").EnumerateArray().Select(e => e.GetString()))).Distinct().Count(), Is.EqualTo(1), $"the structures of set {g.Key} are the same in every run");
            var tm = g.Where(r => r.Extra.TryGetProperty("tm", out _)).Select(r => string.Join(',', r.Extra.GetProperty("tm").EnumerateObject().Select(p => p.Name).Order())).Distinct().Count();
            if (tm > 1) refused[g.Key] = "TM-align wrote different pairs in different runs";
            else if (g.Where(r => r.Section == "s3" && !r.Killed).Select(r => r.Failed).Distinct().Count() > 1) refused[g.Key] = "t_coffee stopped in some runs and ran through in others";
            else if (g.Where(r => r.Section == "s4").Select(r => r.Rows).Distinct().Count() > 1) refused[g.Key] = "runs 5 and 6 of the same files and the same DaliLite behaviour (no alignment) ended differently";
        }

        return refused;
    });

    [Test]
    public void ThreeDCoffee_EqualsTheSapPairLibraryWithoutAnyStructure_OnEverySetOfEveryRun()
    {
        var failures = new List<string>();
        foreach (var r in Section("s2"))
        {
            // -mode 3dcoffee = -in Msap_pair -template_file SELF_P_: no structure is found (the names q000, ... are no PDB entries), so every pair is aligned by proba_pair (once: the default method proba_pair is NOT added to Msap_pair)
            string got = string.Join('|', SequenceAligner.MultipleAlignTCoffeeProgram(r.Seqs, TCoffeeProgramOptions.FromCommandLine("-mode 3dcoffee")).AlignedSequences);
            if (got != r.Rows) failures.Add($"{Describe(r)}: got {got} expected {r.Rows}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of 40; first: {failures.FirstOrDefault()}");
    }

    // ---- s3: DaliLite v5 on the PATH (no `DaliLite` command) -------------------------------------------------------------------------------------------------------------------------------------

    [Test]
    public void DaliPair_NoDaliLiteCommand_EveryRowOfEveryRunEqualsThePort()
    {
        var rows = Section("s3").ToList();
        Assert.That(rows, Has.Count.EqualTo(300));
        var failures = new List<string>();
        int stops = 0, compared = 0, refused = 0, withGlueError = 0;
        foreach (var r in rows)
        {
            if (r.Killed || Unreproducible.Value.ContainsKey(r.Index))
            {
                refused++;
                continue;
            }

            var options = Options(r, null, new Dictionary<(int, int), string>());
            string? got = null;
            string? message = null;
            try { got = Run(r, options); } catch (TCoffeeAbortedException ex) { message = ex.Message; }
            if (r.Failed)
            {
                stops++;
                if (message is null) failures.Add($"{Describe(r)}: the port does not stop, t_coffee ended with {r.Rows}");
                else if (!message.Contains(SequenceAligner.DaliLiteMissingMessage, StringComparison.Ordinal) && !message.Contains("You have not provided any method", StringComparison.Ordinal))
                    failures.Add($"{Describe(r)}: unexpected stop {message}");
                else if (message.Contains(SequenceAligner.DaliLiteMissingMessage, StringComparison.Ordinal) && r.Err.Contains("N -- ERROR: UNSPECIFIED UNSPECIFIED"))
                {
                    // t_coffee's own text: the failed glue command line and the status
                    withGlueError++;
                    Assert.That(r.Err.Any(l => l.StartsWith("N -- COM: <tmp>/N.pl -mode=pdb_pair -method=DaliLite -pdbfile1=<tmp>/N -pdbfile2=<tmp>/N -outfile=<tmp>/N -tmpdir=<tmp>/", StringComparison.Ordinal)), Is.True, Describe(r));
                    Assert.That(r.Err.Any(l => l.StartsWith("# TERMINATION STATUS: FAILURE", StringComparison.Ordinal)), Is.True, Describe(r));
                }

                continue;
            }

            compared++;
            if (got != r.Rows) failures.Add($"{Describe(r)}: got {got ?? message} expected {r.Rows}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} sets differ; first: {failures.FirstOrDefault()}");
        Assert.That(stops, Is.GreaterThan(150));
        Assert.That(withGlueError, Is.GreaterThan(140), "stops whose stderr shows the failed glue command line");
        TestContext.Out.WriteLine($"s3: {rows.Count} sets, {compared} alignments equal, {stops} stops reproduced ({withGlueError} with the glue's error text in stderr), {refused} refused (killed process or not reproducible)");
    }

    // ---- s4: a DaliLite command exists but (run 5) cannot run: every call ends without an alignment -----------------------------------------------------------------------------------

    private static int UsablePairs(Row r)
    {
        var pdb = r.Pdb;
        var usable = Enumerable.Range(0, r.Seqs.Length).Select(i => pdb[i] is { } p && SequenceAligner.TcPdbTemplateUsable(r.Seqs[i], p)).ToArray();
        int n = 0;
        for (int i = 0; i < usable.Length; i++)
            for (int j = i + 1; j < usable.Length; j++)
                if (usable[i] && usable[j]) n++;
        return n;
    }

    [Test]
    public void DaliPair_DaliLiteThatWritesNoAlignment_AddsNoEntry()
    {
        var rows = Section("s4").ToList();
        Assert.That(rows, Has.Count.EqualTo(120));
        var callCounts = new Dictionary<string, (int Calls, int Rc0, int Rc2)>();
        int equal = 0, stops = 0, died = 0, refused = 0;
        var failures = new List<string>();
        foreach (var r in rows)
        {
            int rowCalls = 0;
            callCounts.TryAdd(r.Run, default);
            foreach (var c in r.Extra.GetProperty("calls").EnumerateArray())
            {
                rowCalls++;
                var (n, z, t) = callCounts[r.Run];
                int rc = c.GetProperty("rc").GetInt32();
                callCounts[r.Run] = (n + 1, z + (rc == 0 ? 1 : 0), t + (rc == 2 ? 1 : 0));
                Assert.That(c.GetProperty("args").GetString(), Is.EqualTo("-pairwise pdb1.pdb pdb2.pdb"));
                string output = c.GetProperty("out").GetString()!;
                string error = c.GetProperty("err").GetString()!;
                Assert.That(SequenceAligner.TCoffeeDaliLiteAlignment(output), Is.Null, "the glue finds no Query / Sbjct line: no result.aln");
                if (r.Run == "5")
                    Assert.That(error, Does.Contain("forwarder.pl").Or.Contain("getcwd() failed"));     // run 5: the install path compiled into DaliLite 3.3 did not exist
                else if (rc == 0)
                    Assert.That(output, Does.Contain("Analyse () vs. ()"), "DaliLite 3.3 accepts the CA-only structures and analyses nothing");
                else
                    Assert.That(error, Does.Contain("getcwd() failed"), "the working directory was removed while the program ran");
            }

            if (r.Killed || Unreproducible.Value.ContainsKey(r.Index))
            {
                refused++;
                continue;
            }

            int pairs = UsablePairs(r);
            Assert.That(rowCalls, Is.LessThanOrEqualTo(pairs), Describe(r));
            if (rowCalls < pairs)
            {
                // the process that runs the method ended before it had run every pair job: pdb_pair prints "\n\tpdb_pair/DaliLite failed:\n\t%s\n" with a buffer it never wrote (util_dp_drivers.c L1757-1803: `char command[10000]` is read uninitialised),
                // and the program then stops with an empty library; which pair job ends the process is undefined behaviour of the program: refused
                died++;
                Assert.That(r.Rows, Is.EqualTo("FAIL:1"), Describe(r));
                Assert.That(r.Err, Does.Contain("You have not provided any method or enough Sequences[FATAL]"), Describe(r));
                continue;
            }

            var options = Options(r, TCoffeePairMethod.DaliPair, new Dictionary<(int, int), string>());
            string? got = null;
            string? message = null;
            try { got = Run(r, options); } catch (TCoffeeAbortedException ex) { message = ex.Message; }
            if (r.Failed)
            {
                stops++;
                Assert.That(r.Err, Does.Contain("You have not provided any method or enough Sequences[FATAL]"), Describe(r));
                if (message is null) failures.Add($"{Describe(r)}: every pair job ran, the port has entries, t_coffee stopped");
                continue;
            }

            equal++;
            if (got != r.Rows) failures.Add($"{Describe(r)}: got {got ?? message} expected {r.Rows}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} sets differ; first: {failures.FirstOrDefault()}");
        Assert.That(callCounts["5"], Is.EqualTo((181, 175, 6)));
        Assert.That(callCounts["6"], Is.EqualTo((163, 158, 5)));
        TestContext.Out.WriteLine($"s4: {rows.Count} sets of runs 5 and 6, {callCounts["5"].Calls + callCounts["6"].Calls} DaliLite calls; {equal} alignments equal, {stops} stops (empty library) reproduced, {died} refused (process ended mid-method: uninitialised read), {refused} refused (killed or not reproducible)");
    }

    /// <summary>
    /// The proof of the nondeterminism of the dali_pair sets: runs 5 and 6 ran the same 60 sets (the same files, the same DaliLite that finds no alignment); t_coffee ended differently on some of them (stopped in one run, aligned in the other, or
    /// aligned differently): those sets are the refused ones of <see cref="Unreproducible"/>.
    /// </summary>
    [Test]
    public void DaliPair_RunsFiveAndSixOfTheSameFilesEndDifferentlyOnSomeSets()
    {
        var r5 = Section("s4").Where(r => r.Run == "5").ToDictionary(r => r.Index);
        var r6 = Section("s4").Where(r => r.Run == "6").ToDictionary(r => r.Index);
        var differing = r5.Keys.Where(i => r5[i].Rows != r6[i].Rows).ToList();
        Assert.That(differing, Has.Count.GreaterThan(5), "sets whose result differs between two runs of the same files");
        foreach (int i in differing) Assert.That(Unreproducible.Value.ContainsKey(i), Is.True);
        Assert.That(r5.Keys.Count(i => r5[i].Rows == r6[i].Rows), Is.GreaterThan(30), "the other sets agree");
        TestContext.Out.WriteLine($"s4: {differing.Count} of 60 sets end differently in run 5 and run 6; the files, the pdb sequences and the DaliLite behaviour (no alignment) are the same");
    }

    // ---- the DaliLite v3 glue -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    [Test]
    public void DaliLiteAlignment_IsTheQueryAndSbjctTokensWithDotsAsGaps()
    {
        const string output = "# DaliLite\nQuery         AC..DEF\n              ::  ::\nSbjct         ACGG.EF\nQuery   GHI\nSbjct   G-I\n";
        Assert.That(SequenceAligner.TCoffeeDaliLiteAlignment(output), Is.EqualTo(("AC--DEFGHI", "ACGG-EFG-I")));
        Assert.That(SequenceAligner.TCoffeeDaliLiteAlignment("Query AC\n"), Is.Null);
        Assert.That(SequenceAligner.TCoffeeDaliLiteAlignment(""), Is.Null);
    }
}
