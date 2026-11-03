using System.Diagnostics;
using System.Globalization;
using System.Text;
using SuffixTree.PerfHarness;

// Repeatable macro performance harness for the suffix-tree library (ST-S0 WP 0.7). See README.md.
//
//   dotnet run -c Release --project apps/SuffixTree.PerfHarness -- run --data <dir> [options]
//
// Not part of any test run: nothing here is discovered by `dotnet test`.
return Cli.Run(args);

namespace SuffixTree.PerfHarness
{
    internal static class Cli
    {
        private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

        private const string Usage = """
            SuffixTree.PerfHarness - reproduces the 10 Mbp row of docs/refactoring/suffix-tree-upgrade/BENCHMARK.md

              run --data <dir> [--input dna_10m] [--repeats 3] [--tolerance 0.15] [--work <dir>]
                  [--targets inmem,persistent,load] [--out <file.tsv>] [--report-only] [--keep-trees]
                  <dir> holds data/<input>.txt and queries/<input>.* made by analysis/bench/gen/gen_inputs.py.
                  Runs every measurement in a fresh child process <repeats> times, takes the median and compares it with the
                  baseline of BENCHMARK.md. Exit code 0 = all metrics inside the band (or better), 1 = regression or wrong
                  query result, 2 = usage error. --report-only always exits 0.

              (internal) child inmem|pers-build|pers-query|pers-load <args...>
            """;

        public static int Run(string[] args)
        {
            if (args.Length == 0 || args[0] is "-h" or "--help" or "help")
            {
                Console.WriteLine(Usage);
                return args.Length == 0 ? 2 : 0;
            }

            try
            {
                return args[0] switch
                {
                    "child" => RunChild(args[1..]),
                    "run" => RunParent(args[1..]),
                    _ => UsageError($"unknown command '{args[0]}'"),
                };
            }
            catch (Exception ex) when (ex is ArgumentException or FileNotFoundException or DirectoryNotFoundException or InvalidOperationException)
            {
                return UsageError(ex.Message);
            }
        }

        private static int UsageError(string message)
        {
            Console.Error.WriteLine("error: " + message);
            Console.Error.WriteLine();
            Console.Error.WriteLine(Usage);
            return 2;
        }

        private static int RunChild(string[] a) => a switch
        {
            ["inmem", var text, var q] => Measurements.InMemory(text, q),
            ["pers-build", var text, var tree] => Measurements.PersistentBuild(text, tree),
            ["pers-query", var tree, var q] => Measurements.PersistentQuery(tree, q),
            ["pers-load", var tree] => Measurements.PersistentLoad(tree),
            _ => throw new ArgumentException("bad child arguments"),
        };

        // ───────────────────────── parent ─────────────────────────

        private sealed record Options(string Data, string Input, int Repeats, double Tolerance, string Work,
            HashSet<string> Targets, string? Out, bool ReportOnly, bool KeepTrees);

        private static Options Parse(string[] a)
        {
            string? data = null, input = "dna_10m", work = null, outFile = null;
            int repeats = 3;
            double tolerance = Baseline.DefaultTolerance;
            var targets = new HashSet<string>(StringComparer.Ordinal) { "inmem", "persistent", "load" };
            bool reportOnly = false, keep = false;

            var queue = new Queue<string>(a);
            string Next(string option) => queue.Count > 0 ? queue.Dequeue() : throw new ArgumentException($"missing value after {option}");
            while (queue.Count > 0)
            {
                string option = queue.Dequeue();
                switch (option)
                {
                    case "--data": data = Next(option); break;
                    case "--input": input = Next(option); break;
                    case "--work": work = Next(option); break;
                    case "--out": outFile = Next(option); break;
                    case "--repeats": repeats = int.Parse(Next(option), Inv); break;
                    case "--tolerance": tolerance = double.Parse(Next(option), Inv); break;
                    case "--targets": targets = new HashSet<string>(Next(option).Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries), StringComparer.Ordinal); break;
                    case "--report-only": reportOnly = true; break;
                    case "--keep-trees": keep = true; break;
                    default: throw new ArgumentException($"unknown option '{option}'");
                }
            }
            if (data is null)
                throw new ArgumentException("--data <dir> is required");
            if (repeats < 1)
                throw new ArgumentException("--repeats must be >= 1");
            if (tolerance <= 0 || tolerance >= 1)
                throw new ArgumentException("--tolerance must be in (0, 1)");
            if (targets.Except(new[] { "inmem", "persistent", "load" }).Any())
                throw new ArgumentException("--targets accepts inmem, persistent, load");
            return new Options(Path.GetFullPath(data), input, repeats, tolerance,
                Path.GetFullPath(work ?? Path.Combine(data, "trees")), targets, outFile, reportOnly, keep);
        }

        private static int RunParent(string[] a)
        {
            Options o = Parse(a);
            string text = Path.Combine(o.Data, "data", o.Input + ".txt");
            string queries = Path.Combine(o.Data, "queries", o.Input);
            if (!File.Exists(text))
                throw new FileNotFoundException($"input text not found: {text} (generate it with analysis/bench/gen/gen_inputs.py)");
            Directory.CreateDirectory(o.Work);
            string treeFile = Path.Combine(o.Work, o.Input + ".st");

            // values[(target, op, metric)] = one value per repeat
            var values = new Dictionary<string, List<double>>(StringComparer.Ordinal);
            var started = Stopwatch.StartNew();

            Console.WriteLine($"harness: input={o.Input} repeats={o.Repeats} tolerance=+-{o.Tolerance:P0} targets={string.Join(",", o.Targets.Order())}");
            Console.WriteLine($"runtime: {Environment.Version} {(Environment.Is64BitProcess ? "x64" : "x86")} cores={Environment.ProcessorCount} load-avg-before={LoadAverage()}");
            for (int r = 1; r <= o.Repeats; r++)
            {
                Console.WriteLine($"-- repeat {r}/{o.Repeats} ({started.Elapsed.TotalSeconds:F0} s elapsed)");
                if (o.Targets.Contains("inmem"))
                    Collect(values, "inmem", Child("inmem", text, queries));
                if (o.Targets.Contains("persistent"))
                {
                    Collect(values, "persistent", Child("pers-build", text, treeFile));
                    Collect(values, "persistent", Child("pers-query", treeFile, queries));
                }
                if (o.Targets.Contains("load") && File.Exists(treeFile))
                    Collect(values, "persistent", Child("pers-load", treeFile));
            }

            if (!o.KeepTrees)
            {
                foreach (string f in new[] { treeFile, treeFile + ".children.tmp", treeFile + ".depth.tmp" })
                {
                    if (File.Exists(f))
                        File.Delete(f);
                }
            }

            return Report(o, values);
        }

        private static string LoadAverage()
        {
            try { return File.ReadAllText("/proc/loadavg").Split(' ')[0]; }
            catch (IOException) { return "n/a"; }
        }

        private static List<(string Op, string Metric, double Value)> Child(params string[] childArgs)
        {
            var psi = new ProcessStartInfo
            {
                FileName = Environment.ProcessPath ?? throw new InvalidOperationException("cannot determine the host process"),
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                UseShellExecute = false,
            };
            // `dotnet SuffixTree.PerfHarness.dll ...` needs the dll as first argument; the apphost does not.
            if (Path.GetFileNameWithoutExtension(psi.FileName) == "dotnet")
                psi.ArgumentList.Add(typeof(Cli).Assembly.Location);
            psi.ArgumentList.Add("child");
            foreach (string c in childArgs)
                psi.ArgumentList.Add(c);

            var sw = Stopwatch.StartNew();
            using var p = Process.Start(psi) ?? throw new InvalidOperationException("child process did not start");
            string err = string.Empty;
            p.ErrorDataReceived += (_, e) => err += e.Data + Environment.NewLine;
            p.BeginErrorReadLine();
            var rows = new List<(string, string, double)>();
            string? line;
            while ((line = p.StandardOutput.ReadLine()) is not null)
            {
                string[] f = line.Split('\t');
                if (f.Length >= 5 && f[0] == "R")
                    rows.Add((f[1], f[2], double.Parse(f[3], Inv)));
            }
            p.WaitForExit();
            Console.WriteLine($"   {childArgs[0],-11} {sw.Elapsed.TotalSeconds,6:F1} s  rows={rows.Count}  exit={p.ExitCode}");
            if (p.ExitCode != 0)
                throw new InvalidOperationException($"child '{string.Join(' ', childArgs)}' failed with exit code {p.ExitCode}: {err}");
            return rows;
        }

        private static void Collect(Dictionary<string, List<double>> values, string target, List<(string Op, string Metric, double Value)> rows)
        {
            foreach (var (op, metric, value) in rows)
            {
                string key = $"{target}/{op}/{metric}";
                if (!values.TryGetValue(key, out var list))
                    values[key] = list = new List<double>();
                list.Add(value);
            }
        }

        private static double Median(List<double> v)
        {
            var s = v.Order().ToList();
            return s.Count % 2 == 1 ? s[s.Count / 2] : (s[(s.Count / 2) - 1] + s[s.Count / 2]) / 2;
        }

        private static int Report(Options o, Dictionary<string, List<double>> values)
        {
            var expectations = Baseline.For(o.Input, o.Tolerance);
            var sb = new StringBuilder();
            sb.AppendLine("key\tbaseline\tmedian\tmin\tmax\tratio\tverdict");

            int bad = 0;
            var covered = new HashSet<string>(StringComparer.Ordinal);
            Console.WriteLine();
            Console.WriteLine($"{"metric",-52} {"baseline",16} {"median",16} {"min..max",-30} {"ratio",7}  verdict");
            foreach (var e in expectations)
            {
                if (!o.Targets.Contains(e.Target) && !(e.Target == "persistent" && o.Targets.Contains("load") && e.Op == "warm_load"))
                    continue;
                values.TryGetValue(e.Key, out var list);
                covered.Add(e.Key);
                double? median = list is { Count: > 0 } ? Median(list) : null;
                var verdict = Baseline.Judge(e, median);
                if (verdict is Baseline.Verdict.Regression or Baseline.Verdict.Mismatch or Baseline.Verdict.Missing)
                    bad++;

                string range = list is { Count: > 0 } ? $"{Baseline.Format(list.Min())}..{Baseline.Format(list.Max())}" : "-";
                double ratio = median is { } m && e.Value != 0 ? m / e.Value : double.NaN;
                string ratioText = double.IsNaN(ratio) ? "-" : ratio.ToString("0.000", Inv);
                string name = e.Direction == Direction.Exact ? e.Key + " (exact)" : $"{e.Key} [{e.Unit}]";
                Console.WriteLine($"{name,-52} {Baseline.Format(e.Value),16} {(median is { } mm ? Baseline.Format(mm) : "-"),16} {range,-30} {ratioText,7}  {verdict.ToString().ToUpperInvariant()}");
                sb.AppendLine(string.Join('\t', e.Key, e.Value.ToString("R", Inv), median?.ToString("R", Inv) ?? "", list?.Min().ToString("R", Inv) ?? "", list?.Max().ToString("R", Inv) ?? "", ratioText, verdict));
            }

            // informational rows without a baseline (e.g. strict validation time, p99 latency, other inputs)
            Console.WriteLine();
            foreach (var kv in values.Where(k => !covered.Contains(k.Key) && !k.Key.EndsWith("/checksum", StringComparison.Ordinal)).OrderBy(k => k.Key, StringComparer.Ordinal))
            {
                Console.WriteLine($"{"info " + kv.Key,-52} {"",16} {Baseline.Format(Median(kv.Value)),16}");
                sb.AppendLine(string.Join('\t', kv.Key, "", Median(kv.Value).ToString("R", Inv), kv.Value.Min().ToString("R", Inv), kv.Value.Max().ToString("R", Inv), "", "Info"));
            }

            if (expectations.Count == 0)
                Console.WriteLine($"note: no recorded baseline for input '{o.Input}'; measurements are reported only.");

            if (o.Out is not null)
                File.WriteAllText(o.Out, sb.ToString());

            Console.WriteLine();
            if (covered.Count == 0)
                Console.WriteLine("RESULT: nothing to compare (no baseline for this input).");
            else
                Console.WriteLine(bad == 0
                    ? $"RESULT: all {covered.Count} checked metrics are inside the +-{o.Tolerance:P0} band (or better) and every oracle matches."
                    : $"RESULT: {bad} metric(s) regressed, mismatched or missing.");
            return bad == 0 || o.ReportOnly ? 0 : 1;
        }
    }
}
