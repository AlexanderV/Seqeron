using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto;

internal static class Cli
{
    private const string Usage = """
        SuffixTree.EmProto - prototypes E0-E3, E3b, E4 and E5 of ADR-004 (external-memory construction of the persistent ESA)

          e1   --family <spec> [--data <dir>] [--vectors <sa_vectors.tsv>] [--threads 1] [--k-div 1000 | --k <K>] [--fan 64] [--limit 600] [--tie-sort merge|intro] [--run <id>] [--rows raw.tsv]
          e2   the same with --dcs <X> [--exact] [--dcs-base 16000000]
               in-RAM simulation of the bucket partition and sort, suffix array only; spec: dna_10m, rep_10m, prot_1m, a:N, fib:N, tm:N, tandem:171x10000, nrun:N:RUN, genome:N, dna:N
          pack     --family <spec> --out <packed file>                      text to the b-bit packed file (+ .meta)
          mkranks  --text <packed> --out <rank file> [--period 1024 [--exact]]  difference-cover sample ranks to a file
          e3       --text <packed> --sa <out> --scratch <dir> [--ranks <file>] [--k-div 3000 | --k K] [--threads 4] [--budget-mb 160] [--cgroup-usage <file>]
          e3b      --text <packed> --out <rank file> --scratch <dir> [--period 1024 [--exact]] [--threads 4] [--budget-mb 160] [--scale-n 3000000000] [--verify <oracle rank file>]
                   [--k-s K] [--nbase N] [--chunk R] [--sort-k K] [--fan 1024] [--reservoir 16384] [--run-spill-kb 64] [--io-kb 64] [--merge-budget-mb 16] [--window-kb 1024] [--seed 1] [--touch-buffers] [--case name] [--series]
                   external-memory PH1 (E3b): the difference-cover rank file built with bounded RAM; --scale-n scales K_s, N_base, R, K_sort by n / scale-n from their 3e9 reference values
          e4       --text <packed> --sa <sa file> --out <prefix> --scratch <dir> [--threads 4] [--budget-mb 160] [--scale-n 3000000000] [--chunk R] [--mode distributed|bounded] [--rounds r] [--b2-workers w]
                   [--io-kb 1024] [--tuple-kb 32] [--window-kb 1024] [--touch-buffers] [--verify <oracle prefix>] [--case name] [--series]
                   external-memory PH4 (E4): the compact LCP (<prefix>.lcp base array, <prefix>.ovf overflow pairs) from a packed text and its suffix array file; --scale-n scales R = 8,000,000 by n / scale-n
          e5       --text <packed> --lcp <prefix of .lcp/.ovf> | --fused --sa <sa file>  --out <prefix> --scratch <dir> [--layout a|f|both|stats] [--codec tagged|packed] [--block-entries 262144]
                   [--window-kb 1024] [--slots 8] [--no-flush] [--spill-flush-mb 16] [--spill-dir <dir>] [--ph5-budget-mb 64] [--io-kb 1024] [--verify <oracle prefix>] [--case name] [--series <file>]
                   [--threads 4] [--budget-mb 160] [--scale-n 3000000000] [--chunk R]  (fused: the PH4 options of e4, base mode)
                   external-memory PH5 (E5): <prefix>.cld child table (layout A), <prefix>.psv and <prefix>.nsv (layout F), NODE_COUNT, LRS; the LCP from a compact LCP file or through the PH4 hook
          e5-oracle --text <packed> --lcp <prefix> --out <prefix>   the reference files (.cldT .cldP .psv .nsv .stats) from the in-RAM structures of the product (n < 2^31)
          lcp-oracle --text <packed> --sa <sa file> --out <prefix>   the same two files from the in-RAM LcpKernel (n < 2^31)
          check-lcp  --text <packed> --sa <sa file> --lcp <prefix> [--samples 1000] [--seed s]   structure of the overflow list and sampled ranks against a direct comparison of the text
          packrand --n <symbols> --out <packed file> [--seed s]        packed 2-bit random DNA of any length (streamed; for runs above the 2^31 symbols of pack)
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
                case "e3b": return E3b.Build(o, report);
                case "e4": return E4.Build(o, report);
                case "e5": return E5.Build(o, report);
                case "e5-oracle": return E5.Oracle(o, report);
                case "lcp-oracle": return E4.Oracle(o, report);
                case "check-lcp": return E4.CheckSample(o, report);
                case "packrand": return E3b.PackRandom(o, report);
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
