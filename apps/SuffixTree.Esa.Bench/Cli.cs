using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa;

namespace SuffixTree.Esa.Bench;

internal static class Cli
{
    private const string Usage = """
        SuffixTree.Esa.Bench - SA + LCP build measurements of the ESA backend (work package 1a.8)

          build   --input <file> [--name <label>] [--verify <sa_vectors.tsv>] [--sa-kernel reference|pointer]
                  one cold build in a fresh process: build_prepare_s / build_sais_s / build_sa_s / build_lcp_s / build_total_s, peak RSS
                  (VmHWM), GC allocation and collections, index bytes of SA + LCP + text, SHA-256 of SA and LCP (vector_match with --verify)
          steady  --input <file> [--builds 5] [--name <label>]
                  <builds> builds in one process; the last three are reported (steady state, JIT warm)
          profile --input <file> [--builds 3] [--name <label>]
                  <builds> instrumented builds in one process; phase times of the first (cold) and the last (warm) build
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
            var options = Options.Parse(args[1..]);
            switch (args[0])
            {
                case "build": Commands.Build(options); return 0;
                case "steady": Commands.Steady(options); return 0;
                case "profile": Commands.Profile(options); return 0;
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

    internal static Options Parse(string[] args)
    {
        string? input = null, name = null, verify = null;
        int? builds = null;
        var saisKernel = SaisKernelKind.Reference;
        for (int i = 0; i < args.Length; i += 2)
        {
            string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
            switch (args[i])
            {
                case "--input": input = value; break;
                case "--name": name = value; break;
                case "--verify": verify = value; break;
                case "--builds": builds = int.Parse(value, CultureInfo.InvariantCulture); break;
                case "--sa-kernel": saisKernel = Enum.Parse<SaisKernelKind>(value, ignoreCase: true); break;
                default: throw new ArgumentException($"unknown option '{args[i]}'");
            }
        }
        if (input is null) throw new ArgumentException("--input is required");
        return new Options
        {
            Input = input,
            Name = name ?? Path.GetFileNameWithoutExtension(input),
            Verify = verify,
            Builds = builds ?? 0,
            SaisKernel = saisKernel,
        };
    }
}
