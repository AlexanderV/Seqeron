using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa;

namespace SuffixTree.Esa.Bench;

internal static class Cli
{
    private const string Usage = """
        SuffixTree.Esa.Bench - SA + LCP build measurements of the ESA backend (work package 1a.8)

          build   --input <file> [--name <label>] [--verify <sa_vectors.tsv>] [--sa-kernel reference|pointer]
          build   --input <file> [--name <label>] [--verify <sa_vectors.tsv>] [--lcp-kernel reference|pointer|direct]
                  one cold build in a fresh process: build_prepare_s / build_sais_s / build_sa_s / build_lcp_s / build_total_s, peak RSS
                  (VmHWM), GC allocation and collections, index bytes of SA + LCP + text, SHA-256 of SA and LCP (vector_match with --verify)
          steady  --input <file> [--builds 5] [--name <label>]
                  <builds> builds in one process; the last three are reported (steady state, JIT warm)
          profile --input <file> [--builds 3] [--name <label>]
                  <builds> instrumented builds in one process; phase times of the first (cold) and the last (warm) build
          query   --input <file> [--name <label>] [--queries <bench dir>] [--layout reference3|packed|tagged|all; default: the engine default for the text length] [--trace true] [--min-seconds 1.0]
          query   --synthetic <a^n|abab|fibonacci|thue-morse|n-runs>[:n] ...
                  work packages 1b.1 / 1b.2 (prototype P3 and the throughput gates): SA + LCP once, then per child-table layout the scan time, bytes, write trace
                  (--trace), 20-mer count (present / absent) and 8-mer locate throughput of the in-memory ESA
          analysis --input <file> [--name <label>] [--check true] [--generic true] [--min-length 30]
          analysis --synthetic <a^n|abab|fibonacci|thue-morse|n-runs>[:n] ...
                  work package 1b.3: one full EsaIndex build, then the time of the native analysis operations (longest repeated substring, all ties, distinct-substring
                  counts, suffix enumeration; first call cold and the median of three more); --generic adds the operations that run on Traverse through the generic algorithms
                  (maximal repeated pairs of at least --min-length, k-common substrings of the two halves of the text); --check compares the native results with the
                  Traverse-default path of the same index (exit code 3 on a difference)
          links   --input <file>|--synthetic <a^n|fibonacci|thue-morse|abab|n-runs>[:n] [--name <label>] [--sample 200000] [--seed 1] [--worst 1000] [--fanout 16]
                  [--kernels all|scan,descent,rmq,psv] [--sources all|isa,psi] [--rotate k] [--link-budget 3000000] [--total-budget 400000000] [--time-cap 4] [--search-check 300] [--sweeps true]
                  work package 1b.4a, prototype P1: suffix-link kernels (a) LCP scan, (b) child-table descent, (c) RMQ + PSV/NSV, (c2) PSV/NSV with ISA or Psi endpoints: ns per link
                  (median, p99, max; cold first call), steps, build time and bytes per character of ISA / Psi / the LCP pyramid, definition checks, worst-node sweeps
          mem     --input <file> --queries <bench dir> [--name <label>] [--kernel scan|descent|rmq|psv] [--source isa|psi] [--loci 1|2] [--calls 10] [--min-length 20] [--fanout 16]
                  [--ukkonen true] [--time-links true]
                  work package 1b.4a, prototype P2: MEM of the benchmark's 100 kbp query (1 % substitutions) with the throwaway ESA matching-statistics driver (cold first call, warm calls,
                  links followed, extend steps) and, with --ukkonen, the in-memory Ukkonen tree on the same query in the same process (reference time, MEM-set equality)
          ab      --input <file> [--name <label>] [--verify <sa_vectors.tsv>] [--pairs 9] [--calibration-pairs 7]
                  [--a-args "<build options>"] [--b-args "<build options>"] [--a-env KEY=VALUE]... [--b-env KEY=VALUE]...
                  [--a-dll <path>] [--b-dll <path>] [--libsais <driver>] [--libsais-runs 5] [--max-load 0.30] [--out <raw.tsv>]
                  paired interleaved A/B session: libsais reference, A/A calibration, A/B pairs (fresh process each, order alternates inside the pairs),
                  libsais again, report with median paired ratio, A/A noise band and the verdict of PROCESS.md (>= 3 % and outside the band)
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
            if (args[0] == "ab") return AbRunner.Run(ParseAb(args[1..]));
            if (args[0] == "links")
            {
                LinksCommands.Run(args[1..]);
                return Environment.ExitCode;
            }
            if (args[0] == "mem")
            {
                MemCommands.Run(args[1..]);
                return Environment.ExitCode;
            }
            var options = Options.Parse(args[1..]);
            switch (args[0])
            {
                case "build": Commands.Build(options); return 0;
                case "steady": Commands.Steady(options); return 0;
                case "profile": Commands.Profile(options); return 0;
                case "query": QueryCommands.Query(options); return 0;
                case "analysis": AnalysisCommands.Run(options); return Environment.ExitCode;
                default: return UsageError($"unknown command '{args[0]}'");
            }
        }
        catch (Exception ex) when (ex is ArgumentException or FileNotFoundException or DirectoryNotFoundException or FormatException)
        {
            return UsageError(ex.Message);
        }
    }

    private static AbOptions ParseAb(string[] args)
    {
        string? input = null, name = null, verify = null, dllA = null, dllB = null, libsais = null, raw = null;
        int pairs = 9, calibration = 7, libsaisRuns = 5;
        double maxLoad = 0.30;
        string[] argsA = [], argsB = [];
        var envA = new Dictionary<string, string>();
        var envB = new Dictionary<string, string>();
        for (int i = 0; i < args.Length; i += 2)
        {
            string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
            switch (args[i])
            {
                case "--input": input = value; break;
                case "--name": name = value; break;
                case "--verify": verify = value; break;
                case "--pairs": pairs = int.Parse(value, CultureInfo.InvariantCulture); break;
                case "--calibration-pairs": calibration = int.Parse(value, CultureInfo.InvariantCulture); break;
                case "--a-args": argsA = value.Split(' ', StringSplitOptions.RemoveEmptyEntries); break;
                case "--b-args": argsB = value.Split(' ', StringSplitOptions.RemoveEmptyEntries); break;
                case "--a-env": AddEnv(envA, value); break;
                case "--b-env": AddEnv(envB, value); break;
                case "--a-dll": dllA = value; break;
                case "--b-dll": dllB = value; break;
                case "--libsais": libsais = value; break;
                case "--libsais-runs": libsaisRuns = int.Parse(value, CultureInfo.InvariantCulture); break;
                case "--max-load": maxLoad = double.Parse(value, CultureInfo.InvariantCulture); break;
                case "--out": raw = value; break;
                default: throw new ArgumentException($"unknown option '{args[i]}'");
            }
        }
        if (input is null) throw new ArgumentException("--input is required");
        if (pairs < 7) throw new ArgumentException("--pairs must be at least 7 (PROCESS.md)");
        if (calibration < 5) throw new ArgumentException("--calibration-pairs must be at least 5");
        return new AbOptions
        {
            Input = input,
            Name = name ?? Path.GetFileNameWithoutExtension(input),
            Verify = verify,
            Pairs = pairs,
            CalibrationPairs = calibration,
            ArgsA = argsA,
            ArgsB = argsB,
            EnvA = envA,
            EnvB = envB,
            DllA = dllA,
            DllB = dllB,
            Libsais = libsais,
            LibsaisRuns = libsaisRuns,
            MaxLoad = maxLoad,
            RawOut = raw,
        };
    }

    private static void AddEnv(Dictionary<string, string> env, string assignment)
    {
        int eq = assignment.IndexOf('=', StringComparison.Ordinal);
        if (eq <= 0) throw new ArgumentException($"environment setting '{assignment}' must be KEY=VALUE");
        env[assignment[..eq]] = assignment[(eq + 1)..];
    }

    private static int UsageError(string message)
    {
        Console.Error.WriteLine("error: " + message);
        Console.Error.WriteLine();
        Console.Error.WriteLine(Usage);
        return 2;
    }
}

internal sealed class Options
{
    internal string Input { get; private init; } = "";

    internal string Name { get; private init; } = "";

    internal string? Verify { get; private init; }

    internal int Builds { get; private init; }

    internal SaisKernelKind SaisKernel { get; private init; }
    internal LcpKernelKind LcpKernel { get; private init; }

    internal string? Queries { get; private init; }

    internal string? Synthetic { get; private init; }

    internal ChildLayoutKind[]? Layouts { get; private init; }

    internal bool Trace { get; private init; }

    internal double MinSeconds { get; private init; } = 1.0;

    internal bool Check { get; private init; }

    internal bool Generic { get; private init; }

    internal int MinLength { get; private init; } = 30;

    internal static Options Parse(string[] args)
    {
        string? input = null, name = null, verify = null, queries = null, synthetic = null;
        int? builds = null;
        var saisKernel = SaisKernelKind.Reference;
        var lcpKernel = LcpKernelKind.Reference;
        ChildLayoutKind[]? layouts = null;
        bool trace = false;
        double minSeconds = 1.0;
        bool check = false, generic = false;
        int minLength = 30;
        for (int i = 0; i < args.Length; i += 2)
        {
            string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
            switch (args[i])
            {
                case "--input": input = value; break;
                case "--name": name = value; break;
                case "--verify": verify = value; break;
                case "--check": check = bool.Parse(value); break;
                case "--generic": generic = bool.Parse(value); break;
                case "--min-length": minLength = int.Parse(value, CultureInfo.InvariantCulture); break;
                case "--builds": builds = int.Parse(value, CultureInfo.InvariantCulture); break;
                case "--sa-kernel": saisKernel = Enum.Parse<SaisKernelKind>(value, ignoreCase: true); break;
                case "--lcp-kernel": lcpKernel = Enum.Parse<LcpKernelKind>(value, ignoreCase: true); break;
                case "--queries": queries = value; break;
                case "--synthetic": synthetic = value; break;
                case "--layout": layouts = value.Equals("all", StringComparison.OrdinalIgnoreCase) ? Enum.GetValues<ChildLayoutKind>() : [Enum.Parse<ChildLayoutKind>(value, ignoreCase: true)]; break;
                case "--trace": trace = bool.Parse(value); break;
                case "--min-seconds": minSeconds = double.Parse(value, CultureInfo.InvariantCulture); break;
                default: throw new ArgumentException($"unknown option '{args[i]}'");
            }
        }
        if (input is null && synthetic is null) throw new ArgumentException("--input is required");
        return new Options
        {
            Input = input ?? "",
            Queries = queries,
            Synthetic = synthetic,
            Layouts = layouts,
            Trace = trace,
            MinSeconds = minSeconds,
            Check = check,
            Generic = generic,
            MinLength = minLength,
            Name = name ?? (synthetic is not null ? synthetic.Split(':')[0] : Path.GetFileNameWithoutExtension(input!)),
            Verify = verify,
            Builds = builds ?? 0,
            SaisKernel = saisKernel,
            LcpKernel = lcpKernel,
        };
    }
}
