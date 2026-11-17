using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto;

internal static class Cli
{
    private const string Usage = """
        SuffixTree.EmProto - prototypes E0-E3, E3b, E4, E5, E6, E8 (hand-off package: generator, verifier) and E9 of ADR-004 (external-memory construction of the persistent ESA)

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
                   [--window-kb 1024] [--slots 8] [--no-flush] [--spill-flush-mb 16] [--spill-dir <dir>] [--ph5-budget-mb 64] [--io-kb 1024] [--parallel-nsv] [--nextsym --sa <sa file>]
                   [--verify <oracle prefix>] [--case name] [--series <file>] [--touch-buffers]
                   [--threads 4] [--budget-mb 160] [--scale-n 3000000000] [--chunk R]  (fused: the PH4 options of e4, base mode)
                   external-memory PH5 (E5): <prefix>.cld child table (layout A), <prefix>.psv and <prefix>.nsv (layout F), NODE_COUNT, LRS; the LCP from a compact LCP file or through the PH4 hook
          e5-oracle --text <packed> --lcp <prefix> --out <prefix> [--sa <sa file>]   the reference files (.cldT .cldP .psv .nsv .stats, .nsym with --sa) from the in-RAM structures of the product (n < 2^31)
          lcp-oracle --text <packed> --sa <sa file> --out <prefix>   the same two files from the in-RAM LcpKernel (n < 2^31)
          check-lcp  --text <packed> --sa <sa file> --lcp <prefix> [--samples 1000] [--seed s]   structure of the overflow list and sampled ranks against a direct comparison of the text
          e6       --input <text file, one byte per symbol> --out <prefix> --scratch <dir> [--threads 4] [--scale-n 3000000000 [--rule-n 100000000]] [--period 1024 [--exact]] [--codec packed|tagged] [--budget-mb 160] [--ph5-budget-mb 64]
                   [--no-sync] [--touch-buffers] [--case name] [--series <file>] [--rss-cap-mb 256]
                   [--k K] [--fan F] [--k-s K] [--nbase N] [--dcs-chunk R] [--sort-k K] [--chunk R] [--b2-workers w] [--block-entries b] [--slots s] [--window-kb w] [--prepare-block symbols] ...  (explicit values win over the scale rule; --rule-n evaluates the rule at a fixed n for every text: E7 series A)
                   [--verify-pk <packed file>] [--verify-sa <kernel SA file> | --verify-sa-hash <sha256>] [--verify-lcp <prefix of .lcp/.ovf>] [--verify-child <.cldP>] [--verify-stats <.stats>]
                   end to end (E6): PH0 prepare, PH1 DCS rank file, PH2/PH3 suffix array, PH4 LCP (base mode) fused with PH5 (child table, layout A) in one process; outputs <prefix>.pk (+ .pk.meta) .sa .lcp .ovf .cld .stats; the scratch is deleted
          e9-bwt   --text <packed> --sa <sa file> --out <bwt file> [--read-kb 1024] [--window-kb 1024] [--budget-mb 160] [--no-sampler] [--case name]
                   E9, route (ii): BWT[j] = T[SA[j] - 1] as a separate pass over the SA file (one byte per row; the row of the suffix 0 holds a dummy 0; <bwt>.meta: n, primary row, symbol counts)
          e9-psi   --bwt <bwt file> --out <psi file> [--window-kb 64] [--symbols-per-pass 64] [--read-kb 1024] [--cyclic] [--budget-mb 160] [--pk-meta <.pk.meta>] [--verify <ref psi file>] [--no-sampler] [--case name]
                   E9 / PH6: PSI[i] = ISA[SA[i] + 1] by sigma-way distribution of the BWT (one output cursor per symbol, one mapped window each; ceil(sigma / 64) passes); PSI[0] = 0, or the row of the suffix 0 with --cyclic
          e9-isa   --sa <sa file> --out <isa file> --scratch <dir> [--scale-n 3000000000 | --chunk R] [--io-kb 64] [--read-kb 1024] [--window-kb 1024] [--verify <ref isa file>] [--case name]
                   E9: ISA by DIST (one-level radix distribution of (SA[i], i) into n / R chunk files, then scatter), 24 B per entry
          e9-psi-isa --sa <sa file> --isa <isa file> --out <psi file> [--verify <ref psi file>]    E9: PSI[i] = ISA[SA[i] + 1] from the ISA file (SA sequential, ISA at random)
          e9-oracle  --sa <sa file> --out <prefix>      E9: the in-RAM ISA and Psi of P1 (LinkTables) written as <prefix>.isa and <prefix>.psi (PSI[0] = -1 of P1 stored as 0)
          e9-smoke   --input <text file> --sa <sa file> --isa <isa file> --psi <psi file> [--psi2 <psi from isa>] [--sample 200000] [--seed 1] [--fanout 16] [--link-budget 3000000]
                   E9: every kernel x source of the suffix-link prototype P1 with the in-RAM tables and with the streamed ones, on the P1 node sample; results must be identical
          e8-gen   --n <symbols, decimal suffix k/m/g allowed> --out <text file> [--seed 20240002] [--sums <checksum file>] [--buffer-kb 1024] [--force] [--case name]
                   E8: the genome-like text of E3 to E7 (Families.GenomeLike, same bytes) streamed with a bounded buffer for any n; <out>.sums (default) holds n, seed, SHA-256 and the letter counts
          e8-env   [--case name]    E8: what this process sees: architecture, cores, GC memory, page size, vector instructions (rows env_*)
          e8-verify --prefix <e6 output prefix> [--input <text file>] [--threads 4] [--lce-cap 65536] [--stack-cap-m 64] [--scratch <dir>] [--io-kb 4096] [--no-hash] [--no-child] [--sums <file>] [--quiet] [--case name]
                   E8: checks the outputs .pk .sa .lcp .ovf .cld .stats of e6 without an in-RAM reference (any n below 2^32 - 2): SA is a permutation, adjacent suffixes strictly ordered, stored LCP equal to the common prefix,
                   statistics recomputed from the LCP, the child table cell by cell against the lcp-interval definition, SHA-256 of every output; exit code 0 only when every check passed
          e6-mktext --family <spec> --out <text file> [--data dir]   the ASCII text of a family (setup, outside the pipeline)
          packrand --n <symbols> --out <packed file> [--seed s]        packed 2-bit random DNA of any length (streamed; for runs above the 2^31 symbols of pack)
          hash-sa  --file <sa file>      SHA-256 of entries 1..n (libsais convention)
          oracle   --family <spec> [--out <sa file>]      same hash from SaisKernel in RAM (and the kernel's suffix array as a file of n + 1 little-endian int32)
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
                case "e6": return E6.Run(o, report);
                case "e6-mktext": return E6.MakeText(o, report);
                case "e8-gen": return E8.Gen(o, report);
                case "e8-verify": return E8.Verify(o, report);
                case "e8-env": return E8.Env(o, report);
                case "e9-bwt": return E9.Bwt(o, report);
                case "e9-psi": return E9.Psi(o, report);
                case "e9-isa": return E9.Isa(o, report);
                case "e9-psi-isa": return E9.PsiFromIsaCommand(o, report);
                case "e9-oracle": return E9.Oracle(o, report);
                case "e9-smoke": return E9.Smoke(o, report);
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
