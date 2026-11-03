using System.Diagnostics;
using System.Globalization;
using ChildProcess = System.Diagnostics.Process;
using System.Text;

namespace SuffixTree.Esa.Bench;

/// <summary>Parameters of one A/B session (see <c>Cli</c> usage).</summary>
internal sealed class AbOptions
{
    internal required string Input { get; init; }

    internal required string Name { get; init; }

    internal string? Verify { get; init; }

    internal int Pairs { get; init; } = 9;

    internal int CalibrationPairs { get; init; } = 7;

    internal string[] ArgsA { get; init; } = [];

    internal string[] ArgsB { get; init; } = [];

    internal Dictionary<string, string> EnvA { get; init; } = [];

    internal Dictionary<string, string> EnvB { get; init; } = [];

    internal string? DllA { get; init; }

    internal string? DllB { get; init; }

    internal string? Libsais { get; init; }

    internal int LibsaisRuns { get; init; } = 5;

    internal double MaxLoad { get; init; } = 0.30;

    internal string? RawOut { get; init; }
}

/// <summary>
/// The paired interleaved A/B runner with A/A calibration (work package 1o.0, PROCESS.md amendment 2026-11-03). Every measurement is a FRESH process of the
/// <c>build</c> command (cold JIT, own peak RSS). A session is: quiet-machine check, libsais reference runs, A/A calibration (the configuration A
/// against itself, N pairs), the A/B comparison (N pairs), libsais reference runs again, report and verdict. Inside a pair the order alternates
/// (A B, B A, A B, ...) so that a drift or an order effect cancels in the ratio. A child that does not reproduce the libsais vectors aborts the session.
/// </summary>
internal static class AbRunner
{
    private static readonly string[] Metrics = ["build_total_s", "build_sa_s", "build_lcp_s", "peak_rss_hwm_after_build"];

    internal static int Run(AbOptions o)
    {
        var raw = new StringBuilder();
        Console.WriteLine($"session {DateTime.UtcNow:yyyy-MM-dd HH:mm:ss}Z input {o.Name} pairs {o.Pairs} calibration pairs {o.CalibrationPairs}");
        Console.WriteLine($"A: dll {o.DllA ?? "(this)"} args [{string.Join(' ', o.ArgsA)}] env [{Env(o.EnvA)}]");
        Console.WriteLine($"B: dll {o.DllB ?? "(this)"} args [{string.Join(' ', o.ArgsB)}] env [{Env(o.EnvB)}]");
        double load = LoadAverage();
        Console.WriteLine($"load average at start (1 min): {load:F2} (limit {o.MaxLoad:F2})");
        if (load > o.MaxLoad)
        {
            Console.WriteLine("the machine is not quiet: refusing to time anything (check uptime / ps, wait, retry)");
            return 4;
        }

        List<double> libsaisPre = LibsaisRuns(o, "pre");
        var aaA = new List<Dictionary<string, double>>();
        var aaB = new List<Dictionary<string, double>>();
        if (!Pairs(o, "AA", o.CalibrationPairs, o.ArgsA, o.EnvA, o.DllA, o.ArgsA, o.EnvA, o.DllA, aaA, aaB, raw)) return 3;
        var abA = new List<Dictionary<string, double>>();
        var abB = new List<Dictionary<string, double>>();
        if (!Pairs(o, "AB", o.Pairs, o.ArgsA, o.EnvA, o.DllA, o.ArgsB, o.EnvB, o.DllB, abA, abB, raw)) return 3;
        List<double> libsaisPost = LibsaisRuns(o, "post");
        Console.WriteLine($"load average at end (1 min): {LoadAverage():F2}");

        Report(o, aaA, aaB, abA, abB, libsaisPre, libsaisPost);
        if (o.RawOut is not null) File.WriteAllText(o.RawOut, "stage\tpair\tside\tmetric\tvalue\n" + raw);
        return 0;
    }

    private static bool Pairs(
        AbOptions o,
        string stage,
        int pairs,
        string[] argsA,
        Dictionary<string, string> envA,
        string? dllA,
        string[] argsB,
        Dictionary<string, string> envB,
        string? dllB,
        List<Dictionary<string, double>> resultsA,
        List<Dictionary<string, double>> resultsB,
        StringBuilder raw)
    {
        for (int p = 0; p < pairs; p++)
        {
            bool aFirst = p % 2 == 0;
            Dictionary<string, double>? first = null, second = null;
            for (int s = 0; s < 2; s++)
            {
                bool isA = (s == 0) == aFirst;
                Dictionary<string, double>? r = Child(o, isA ? argsA : argsB, isA ? envA : envB, isA ? dllA : dllB);
                if (r is null)
                {
                    Console.WriteLine($"{stage} pair {p + 1}: child failed or did not reproduce the libsais vectors: session aborted");
                    return false;
                }
                if (s == 0) first = r;
                else second = r;
                foreach (string m in Metrics)
                    if (r.TryGetValue(m, out double v)) raw.Append(CultureInfo.InvariantCulture, $"{stage}\t{p + 1}\t{(isA ? "A" : "B")}\t{m}\t{v:R}\n");
            }
            resultsA.Add(aFirst ? first! : second!);
            resultsB.Add(aFirst ? second! : first!);
            Console.WriteLine($"{stage} pair {p + 1,2} ({(aFirst ? "A B" : "B A")}): total A {resultsA[^1]["build_total_s"]:F3} s, B {resultsB[^1]["build_total_s"]:F3} s, B/A {resultsB[^1]["build_total_s"] / resultsA[^1]["build_total_s"]:F4}");
        }
        return true;
    }

    private static Dictionary<string, double>? Child(AbOptions o, string[] extra, Dictionary<string, string> env, string? dll)
    {
        string entry = dll ?? Assembly();
        var psi = new ProcessStartInfo { RedirectStandardOutput = true, UseShellExecute = false };
        string host = Environment.ProcessPath ?? "dotnet";
        psi.FileName = host;
        if (Path.GetFileNameWithoutExtension(host).Equals("dotnet", StringComparison.OrdinalIgnoreCase)) psi.ArgumentList.Add(entry);
        foreach (string a in new[] { "build", "--input", o.Input, "--name", o.Name }) psi.ArgumentList.Add(a);
        if (o.Verify is not null)
        {
            psi.ArgumentList.Add("--verify");
            psi.ArgumentList.Add(o.Verify);
        }
        foreach (string a in extra) psi.ArgumentList.Add(a);
        psi.Environment["DOTNET_gcServer"] = "0";
        psi.Environment["DOTNET_CLI_TELEMETRY_OPTOUT"] = "1";
        foreach (var (k, v) in env) psi.Environment[k] = v;

        using ChildProcess process = ChildProcess.Start(psi)!;
        string output = process.StandardOutput.ReadToEnd();
        process.WaitForExit();
        if (process.ExitCode != 0) return null;
        var values = new Dictionary<string, double>();
        foreach (string line in output.Split('\n'))
        {
            string[] f = line.TrimEnd('\r').Split('\t');
            if (f.Length >= 4 && f[0] == "R" && f[1] == "build" && double.TryParse(f[3], NumberStyles.Float, CultureInfo.InvariantCulture, out double v))
                values[f[2]] = v;
            else if (f.Length >= 4 && f[0] == "R" && f[2] == "vector_match" && f[3] != "1") return null;
        }
        return values.ContainsKey("build_total_s") ? values : null;
    }

    private static List<double> LibsaisRuns(AbOptions o, string when)
    {
        var totals = new List<double>();
        if (o.Libsais is null) return totals;
        for (int i = 0; i < o.LibsaisRuns; i++)
        {
            var psi = new ProcessStartInfo(o.Libsais) { RedirectStandardOutput = true, UseShellExecute = false };
            psi.ArgumentList.Add(o.Input);
            using ChildProcess process = ChildProcess.Start(psi)!;
            string output = process.StandardOutput.ReadToEnd();
            process.WaitForExit();
            foreach (string token in output.Split([' ', '\n'], StringSplitOptions.RemoveEmptyEntries))
                if (token.StartsWith("total=", StringComparison.Ordinal))
                    totals.Add(double.Parse(token[6..], CultureInfo.InvariantCulture));
        }
        Console.WriteLine($"libsais SA+PLCP+LCP ({when}): {string.Join(' ', totals.Select(t => t.ToString("F3", CultureInfo.InvariantCulture)))}; median {(totals.Count > 0 ? PairedStats.Median(totals) : double.NaN):F3} s");
        return totals;
    }

    private static void Report(
        AbOptions o,
        List<Dictionary<string, double>> aaA,
        List<Dictionary<string, double>> aaB,
        List<Dictionary<string, double>> abA,
        List<Dictionary<string, double>> abB,
        List<double> libsaisPre,
        List<double> libsaisPost)
    {
        Console.WriteLine();
        Console.WriteLine("== report ==");
        foreach (string m in Metrics)
        {
            double[] a = [.. abA.Select(r => r[m])], b = [.. abB.Select(r => r[m])];
            double[] ra = [.. aaA.Select(r => r[m])], rb = [.. aaB.Select(r => r[m])];
            List<double> ratios = PairedStats.Ratios(a, b), aaRatios = PairedStats.Ratios(ra, rb);
            var band = PairedStats.NoiseBand(aaRatios);
            var ci = PairedStats.MedianInterval(ratios);
            string unit = m.StartsWith("peak", StringComparison.Ordinal) ? "B" : "s";
            Console.WriteLine($"{m} [{unit}]: A median {PairedStats.Median(a):G5} ({a.Min():G5}-{a.Max():G5}), B median {PairedStats.Median(b):G5} ({b.Min():G5}-{b.Max():G5}); paired B/A median {PairedStats.Median(ratios):F4} (min {ratios.Min():F4}, max {ratios.Max():F4}, 95 % interval of the median {ci.Low:F4}-{ci.High:F4}); A/A median {PairedStats.Median(aaRatios):F4} (min {aaRatios.Min():F4}, max {aaRatios.Max():F4}), noise band [{band.Low:F4}, {band.High:F4}]");
        }
        List<double> total = PairedStats.Ratios([.. abA.Select(r => r["build_total_s"])], [.. abB.Select(r => r["build_total_s"])]);
        var totalBand = PairedStats.NoiseBand(PairedStats.Ratios([.. aaA.Select(r => r["build_total_s"])], [.. aaB.Select(r => r["build_total_s"])]));
        var (accept, _, reason) = PairedStats.Decide(PairedStats.Median(total), totalBand, PairedStats.MedianInterval(total));
        Console.WriteLine($"VERDICT (SA+LCP, build_total_s): {(accept ? "ACCEPT" : "REJECT")} - {reason}");
        List<double> libsais = [.. libsaisPre, .. libsaisPost];
        if (libsais.Count > 0)
        {
            double lib = PairedStats.Median(libsais);
            double a = PairedStats.Median([.. abA.Select(r => r["build_total_s"])]), b = PairedStats.Median([.. abB.Select(r => r["build_total_s"])]);
            Console.WriteLine($"same-session libsais SA+PLCP+LCP median {lib:F3} s (pre {PairedStats.Median(libsaisPre):F3}, post {PairedStats.Median(libsaisPost):F3}, drift {(PairedStats.Median(libsaisPost) / PairedStats.Median(libsaisPre)) - 1:P1}); A/libsais {a / lib:F2}x, B/libsais {b / lib:F2}x");
        }
        _ = o;
    }

    private static string Env(Dictionary<string, string> env) => string.Join(' ', env.Select(kv => $"{kv.Key}={kv.Value}"));

    private static string Assembly() => System.Reflection.Assembly.GetEntryAssembly()!.Location;

    private static double LoadAverage()
    {
        try
        {
            string[] f = File.ReadAllText("/proc/loadavg").Split(' ');
            return double.Parse(f[0], CultureInfo.InvariantCulture);
        }
        catch (IOException)
        {
            return 0;
        }
    }
}
