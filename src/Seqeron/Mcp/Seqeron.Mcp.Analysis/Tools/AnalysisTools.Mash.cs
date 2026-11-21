using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Analysis.Tools;

// B09 A5-10: the text commands of Mash 2.3 (bounds, dist, triangle) on top of the MinHash sketching of KmerAnalyzer; every tool delegates to MashCommands.
public partial class AnalysisTools
{
    private static global::Seqeron.Genomics.Analysis.MashInput[] ToMashInputs(string[] files, string[]? names, string prefix, string paramName)
    {
        if (files is null || files.Length == 0)
            throw new ArgumentException("At least one input is required", paramName);
        if (names is not null && names.Length != files.Length)
            throw new ArgumentException("names needs one name per input", nameof(names));
        var inputs = new global::Seqeron.Genomics.Analysis.MashInput[files.Length];
        for (int i = 0; i < files.Length; i++)
        {
            string text = files[i] ?? throw new ArgumentException($"Input {i} is null", paramName);
            var records = new List<global::Seqeron.Genomics.Analysis.MashRecord>();
            string id = "", comment = "";
            var seq = new System.Text.StringBuilder();
            bool open = false;
            void Close()
            {
                if (open)
                    records.Add(new global::Seqeron.Genomics.Analysis.MashRecord(id, comment, seq.ToString()));
                seq.Clear();
            }

            foreach (string raw in text.Split('\n'))
            {
                string line = raw.TrimEnd('\r');
                if (line.StartsWith('>'))
                {
                    Close();
                    open = true;
                    string header = line[1..].Trim();
                    int sp = header.IndexOfAny([' ', '\t']);
                    id = sp < 0 ? header : header[..sp];
                    comment = sp < 0 ? "" : header[(sp + 1)..].Trim();
                }
                else if (open)
                {
                    foreach (char c in line)
                        if (!char.IsWhiteSpace(c))
                            seq.Append(c);
                }
            }

            Close();
            if (records.Count == 0)
                throw new ArgumentException($"Input {i} has no FASTA record", paramName);
            inputs[i] = new global::Seqeron.Genomics.Analysis.MashInput(names is null ? prefix + (i + 1).ToString(System.Globalization.CultureInfo.InvariantCulture) : names[i], records);
        }

        return inputs;
    }

    private static uint RequireMashSeed(long seed)
    {
        if (seed is < 0 or > uint.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(seed), seed, "seed must be in 0..4294967295 (mash -S).");
        return (uint)seed;
    }

    [McpServerTool(Name = "mash_bounds", Title = "Mash — Error Bounds Table (mash bounds)", ReadOnly = true)]
    [Description("Mash 2.3 'mash bounds -k K -p P', byte for byte: the table of error bounds of the Mash distance and of the screen distance for the sketch sizes 100 ... 1000000 and the distances 0.05 ... 0.4 (CommandBounds.cpp): for the binomial model of the number of shared hashes the smallest x whose cumulative probability exceeds (1 - p) / 2 gives the lower Jaccard bound, converted back to a distance and printed as the error (bound minus distance; inf when no hash is certain to be shared). The numbers are printed as C++ cout does (6 significant digits). VERIFICATION: 220 / 220 random (k, p) equal the installed mash 2.3.")]
    public static MashTextResult MashBounds(
        [Description("k-mer size (1 - 32, default 21).")] int k = 21,
        [Description("Confidence of the bounds (0 - 1 exclusive, default 0.99).")] double probability = 0.99)
        => new(global::Seqeron.Genomics.Analysis.MashCommands.Bounds(k, probability));

    [McpServerTool(Name = "mash_dist", Title = "Mash — Distances of Queries to References (mash dist)", ReadOnly = true)]
    [Description("Mash 2.3 'mash dist': MinHash (bottom-s) sketches of the reference and query files (canonical k-mers, MurmurHash3, the sketch of Ondov et al. 2016: KmerAnalyzer.CreateMinHashSketch) and, for every query in order and every reference, the line 'reference TAB query TAB distance TAB p-value TAB shared/total' of the pairs whose p-value is at most maxPValue (-v) and distance at most maxDistance (-d); table = true (-t) writes the '#query' header and a row of distances per query (a pair that fails a filter leaves its field blank); individual = true (-i) sketches every sequence on its own, named by its id (sequences shorter than k are skipped). A file without a record of at least k symbols is refused as the program refuses it (exit status 1: 'All fasta records in input files were shorter than the k-mer size'). Numbers are printed as C++ cout does (6 significant digits). VERIFICATION: 220 / 220 random sets (1 - 4 references x 1 - 3 queries of 1 - 3 records, k 4 - 32, sketch size 5 - 2000, seeds, -n, -i, -t, -v, -d) equal the installed mash 2.3 byte for byte.")]
    public static MashTextResult MashDist(
        [Description("The reference files as FASTA texts (records start with '>': the first word of the header is the sequence id, the rest the comment).")] string[] references,
        [Description("The query files (FASTA texts).")] string[] queries,
        [Description("Names of the reference files as the program prints them (default reference1, ...).")] string[]? referenceNames = null,
        [Description("Names of the query files (default query1, ...).")] string[]? queryNames = null,
        [Description("k-mer size (-k, 1 - 32, default 21).")] int k = 21,
        [Description("Sketch size (-s, default 1000).")] int sketchSize = 1000,
        [Description("Hash seed (-S, 0 - 4294967295, default 42).")] long seed = 42,
        [Description("-n: preserve the strand (no canonical k-mers).")] bool noncanonical = false,
        [Description("-i: sketch individual sequences instead of whole files.")] bool individual = false,
        [Description("-t: table output.")] bool table = false,
        [Description("-v: maximum p-value to report (0 - 1, default 1).")] double maxPValue = 1.0,
        [Description("-d: maximum distance to report (0 - 1, default 1).")] double maxDistance = 1.0,
        CancellationToken cancellationToken = default)
        => new(global::Seqeron.Genomics.Analysis.MashCommands.Dist(ToMashInputs(references, referenceNames, "reference", nameof(references)), ToMashInputs(queries, queryNames, "query", nameof(queries)), k, sketchSize, RequireMashSeed(seed),
            !noncanonical, individual, table, maxPValue, maxDistance, cancellationToken));

    [McpServerTool(Name = "mash_triangle", Title = "Mash — Lower-Triangular Distance Matrix (mash triangle)", ReadOnly = true)]
    [Description("Mash 2.3 'mash triangle': the distance of each input to every other input as a lower-triangular matrix in relaxed Phylip format (a first line with a TAB and the number of inputs, then one line per input: its name and the distances to the inputs before it), or with edgeList = true (-E; implied by a p-value or distance filter) the edge list 'name TAB name TAB distance TAB p-value TAB shared/total' of the pairs that pass maxPValue (-v) and maxDistance (-d). Sketches as mash_dist (a file without a record of at least k symbols is refused). VERIFICATION: 220 / 220 random sets (2 - 6 files, k 4 - 32, sketch size 5 - 2000, -n, -i, -E, -v, -d) equal the installed mash 2.3 byte for byte.")]
    public static MashTextResult MashTriangle(
        [Description("The input files (FASTA texts), at least two for a matrix.")] string[] sequences,
        [Description("Names of the inputs as the program prints them (default input1, ...).")] string[]? names = null,
        [Description("k-mer size (-k, 1 - 32, default 21).")] int k = 21,
        [Description("Sketch size (-s, default 1000).")] int sketchSize = 1000,
        [Description("Hash seed (-S, 0 - 4294967295, default 42).")] long seed = 42,
        [Description("-n: preserve the strand (no canonical k-mers).")] bool noncanonical = false,
        [Description("-i: sketch individual sequences instead of whole files.")] bool individual = false,
        [Description("-E: output the edge list instead of the matrix.")] bool edgeList = false,
        [Description("-v: maximum p-value to report in the edge list (0 - 1, default 1).")] double maxPValue = 1.0,
        [Description("-d: maximum distance to report in the edge list (0 - 1, default 1).")] double maxDistance = 1.0,
        CancellationToken cancellationToken = default)
        => new(global::Seqeron.Genomics.Analysis.MashCommands.Triangle(ToMashInputs(sequences, names, "input", nameof(sequences)), k, sketchSize, RequireMashSeed(seed), !noncanonical, individual, edgeList, maxPValue, maxDistance, cancellationToken));

    [McpServerTool(Name = "mash_screen", Title = "Mash — Containment of Queries in a Mixture (mash screen)", ReadOnly = true)]
    [Description("Mash 2.3 'mash screen' for nucleotide queries, byte for byte (CommandScreen.cpp): the MinHash sketches of the query files (as mash sketch -k -s -S -n -i would write them into the query .msh) are screened against a mixture (reads or contigs, any number of files): the hashes of the queries are counted in every window of the mixture and for each query one line identity TAB shared/hashes TAB median-multiplicity TAB p-value TAB query-ID TAB query-comment is written, identity = (shared / hashes)^(1/k), the median multiplicity over the shared hashes, the p-value the binomial tail with r = estimated distinct k-mers of the mixture / 4^k. winnerTakesAll = -w (a hash held by several queries counts for the query with the best identity, ties the longer one), minIdentity = -i (-1 writes every query, 0 those with a shared hash), maxPValue = -v. The query-comment is 'name comment' of the first record of a file ('[N seqs] ... [...]' for several records; the record's own comment with individual). Amino-acid queries (6-frame translation of the mixture) are not part of the port. VERIFICATION: 240 / 240 random screens equal the installed mash 2.3.")]
    public static MashTextResult MashScreen(
        [Description("The query files as FASTA texts (records start with '>': the first word of the header is the sequence id, the rest the comment).")] string[] queries,
        [Description("The mixture files as FASTA texts (reads or contigs).")] string[] mixtures,
        [Description("Names of the query files as the program prints them (default query1, ...).")] string[]? queryNames = null,
        [Description("k-mer size of the query sketch (-k of mash sketch, 1 - 32, default 21).")] int k = 21,
        [Description("Sketch size of the query sketch (-s, default 1000).")] int sketchSize = 1000,
        [Description("Hash seed of the query sketch (-S, 0 - 4294967295, default 42).")] long seed = 42,
        [Description("The query sketch was made with -n (strand preserved).")] bool noncanonical = false,
        [Description("The query sketch was made with -i (one reference per sequence).")] bool individual = false,
        [Description("-w: winner-takes-all strategy for identity estimates.")] bool winnerTakesAll = false,
        [Description("-i: minimum identity to report (-1 - 1, default 0).")] double minIdentity = 0.0,
        [Description("-v: maximum p-value to report (0 - 1, default 1).")] double maxPValue = 1.0,
        CancellationToken cancellationToken = default)
        => new(global::Seqeron.Genomics.Analysis.MashCommands.Screen(ToMashInputs(queries, queryNames, "query", nameof(queries)), ToMashInputs(mixtures, null, "mixture", nameof(mixtures)), k, sketchSize, RequireMashSeed(seed),
            !noncanonical, individual, winnerTakesAll, minIdentity, maxPValue, cancellationToken));
}
