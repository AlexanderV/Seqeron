using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto;

internal static class Cli
{
    private const string Usage = """
        SuffixTree.EmProto - prototypes E0-E3 of ADR-004 (external-memory construction of the persistent ESA)

          e1   --family <spec> [--data <dir>] [--vectors <sa_vectors.tsv>] [--threads 1] [--k-div 1000 | --k <K>] [--fan 64] [--limit 600] [--tie-sort merge|intro] [--run <id>] [--rows raw.tsv]
          e2   the same with --dcs <X> [--exact] [--dcs-base 16000000]
               in-RAM simulation of the bucket partition and sort, suffix array only; spec: dna_10m, rep_10m, prot_1m, a:N, fib:N, tm:N, tandem:171x10000, nrun:N:RUN, genome:N, dna:N
          pack     --family <spec> --out <packed file>                      text to the b-bit packed file (+ .meta)
          mkranks  --text <packed> --out <rank file> [--period 1024 [--exact]]  difference-cover sample ranks to a file
          e3       --text <packed> --sa <out> --scratch <dir> [--ranks <file>] [--k-div 3000 | --k K] [--threads 4] [--budget-mb 160] [--cgroup-usage <file>]
          hash-sa  --file <sa file>      SHA-256 of entries 1..n (libsais convention)
          oracle   --family <spec>       same hash from SaisKernel in RAM
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
            // --rows names the raw-row file; for e0 / e1 / e2 the older spelling --out is accepted too (in the E3 commands --out is a data file)
            var report = new Report(o.Get("rows") ?? (args[0] is "e0" or "e1" or "e2" ? o.Get("out") : null), o.Get("run"));
            switch (args[0])
            {
                case "e0": return E0.Run(o, report);
                case "e1": return E1.Run(o, report, "E1");
                case "e2": return E1.Run(o, report, "E2");
                case "pack": return E3.Pack(o, report);
                case "mkranks": return E3.MakeRanks(o, report);
                case "e3": return E3.Build(o, report);
                case "hash-sa": return E3.HashSa(o, report);
                case "oracle": return E3.OracleHash(o, report);
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
