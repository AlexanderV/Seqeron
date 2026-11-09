using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto;

internal static class Cli
{
    private const string Usage = """
        SuffixTree.EmProto - prototypes E0-E3 of ADR-004 (external-memory construction of the persistent ESA)

          e0   --dir <scratch directory> [--pairs 15] [--io-mb 256] [--out raw.tsv]
               harness gates: RssAnon sampler overhead, I/O counters against /proc/self/io, arena budget
        """;

    internal static int Run(string[] args)
    {
        if (args.Length == 0 || args[0] is "-h" or "--help" or "help")
        {
            Console.WriteLine(Usage);
            return args.Length == 0 ? 2 : 0;
        }
        try
        {
            Options o = Options.Parse(args[1..]);
            var report = new Report(o.Get("out"));
            switch (args[0])
            {
                case "e0": return E0.Run(o, report);
                default:
                    Console.Error.WriteLine($"unknown command '{args[0]}'");
                    return 2;
            }
        }
        catch (ArgumentException ex)
        {
            Console.Error.WriteLine(ex.Message);
            return 2;
        }
    }
}
