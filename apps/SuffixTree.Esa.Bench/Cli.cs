using System.Diagnostics;
using System.Globalization;

namespace SuffixTree.Esa.Bench;

internal static class Cli
{
    private const string Usage = """
        SuffixTree.Esa.Bench - SA + LCP build measurements of the ESA backend (work package 1a.8)

          build   --input <file> [--name <label>] [--verify <sa_vectors.tsv>]
                  one cold build in a fresh process: build_prepare_s / build_sais_s / build_sa_s / build_lcp_s / build_total_s, peak RSS
                  (VmHWM), GC allocation and collections, index bytes of SA + LCP + text, SHA-256 of SA and LCP (vector_match with --verify)
          steady  --input <file> [--builds 5] [--name <label>]
                  <builds> builds in one process; the last three are reported (steady state, JIT warm)
          profile --input <file> [--builds 3] [--name <label>]
                  <builds> instrumented builds in one process; phase times of the first (cold) and the last (warm) build
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

    internal static Options Parse(string[] args)
    {
        string? input = null, name = null, verify = null;
        int? builds = null;
        for (int i = 0; i < args.Length; i += 2)
        {
            string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
            switch (args[i])
            {
                case "--input": input = value; break;
                case "--name": name = value; break;
                case "--verify": verify = value; break;
                case "--builds": builds = int.Parse(value, CultureInfo.InvariantCulture); break;
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
        };
    }
}
