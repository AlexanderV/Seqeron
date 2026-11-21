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
            if (text.StartsWith(MshPrefix, StringComparison.Ordinal))
            {
                // a sketch file (x.msh): "msh:" + the base64 of the file (the output of mash_sketch / mash_paste)
                inputs[i] = new global::Seqeron.Genomics.Analysis.MashInput(names is null ? prefix + (i + 1).ToString(System.Globalization.CultureInfo.InvariantCulture) : names[i], [], ReadMshText(text, paramName));
                continue;
            }

            var parsed = global::Seqeron.Genomics.Analysis.MashCommands.ParseFasta(names is null ? prefix + (i + 1).ToString(System.Globalization.CultureInfo.InvariantCulture) : names[i], text);
            if (parsed.Records.Count == 0)
                throw new ArgumentException($"Input {i} has no FASTA record", paramName);
            inputs[i] = parsed;
        }

        return inputs;
    }

    private const string MshPrefix = "msh:";

    private static global::Seqeron.Genomics.Analysis.MashSketchFileData ReadMshText(string text, string paramName)
    {
        try
        {
            return global::Seqeron.Genomics.Analysis.MashSketchFile.Read(Convert.FromBase64String(text[MshPrefix.Length..].Trim()));
        }
        catch (Exception e) when (e is FormatException or InvalidDataException)
        {
            throw new ArgumentException("The sketch file is not valid base64 of a Mash .msh file: " + e.Message, paramName);
        }
    }

    private static global::Seqeron.Genomics.Analysis.MashCommands.MashReadsOptions? ReadsOptions(bool reads, int minCopies, double targetCoverage, long genomeSize, string? bloomFilter = null)
    {
        long bloomBytes = string.IsNullOrEmpty(bloomFilter) ? 0 : global::Seqeron.Genomics.Analysis.MashCommands.ParseSize(bloomFilter);
        if (!reads && bloomBytes == 0)
            return null;
        if (bloomBytes > 0 && minCopies != 1)
            throw new ArgumentException("minCopies cannot be used with bloomFilter (mash: the option -m cannot be used with -b).", nameof(minCopies));
        if (minCopies < 1)
            throw new ArgumentOutOfRangeException(nameof(minCopies), minCopies, "minCopies must be at least 1 (mash -m).");
        if (targetCoverage < 0 || double.IsNaN(targetCoverage))
            throw new ArgumentOutOfRangeException(nameof(targetCoverage), targetCoverage, "targetCoverage must be 0 or more (mash -c).");
        if (genomeSize < 0)
            throw new ArgumentOutOfRangeException(nameof(genomeSize), genomeSize, "genomeSize must be 0 or more (mash -g).");
        return new global::Seqeron.Genomics.Analysis.MashCommands.MashReadsOptions(minCopies, targetCoverage, genomeSize, bloomBytes);
    }

    private static int MashK(int? k, bool protein) => global::Seqeron.Genomics.Analysis.MashAlphabet.ResolveKmerSize(k, protein);

    private static global::Seqeron.Genomics.Analysis.MashAlphabet MashAlphabetOf(bool protein, string? alphabet, bool preserveCase)
        => global::Seqeron.Genomics.Analysis.MashAlphabet.FromOptions(protein, alphabet, preserveCase);

    private static uint RequireMashSeed(long seed)
    {
        if (seed is < 0 or > uint.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(seed), seed, "seed must be in 0..4294967295 (mash -S).");
        return global::Seqeron.Genomics.Analysis.MashCommands.SeedFromOption(seed);   // the program stores 4294967295 as 0 (measured)
    }

    [McpServerTool(Name = "mash_bounds", Title = "Mash — Error Bounds Table (mash bounds)", ReadOnly = true)]
    [Description("Mash 2.3 'mash bounds -k K -p P', byte for byte: the table of error bounds of the Mash distance and of the screen distance for the sketch sizes 100 ... 1000000 and the distances 0.05 ... 0.4 (CommandBounds.cpp): for the binomial model of the number of shared hashes the smallest x whose cumulative probability exceeds (1 - p) / 2 gives the lower Jaccard bound, converted back to a distance and printed as the error (bound minus distance; inf when no hash is certain to be shared). The numbers are printed as C++ cout does (6 significant digits). VERIFICATION: 220 / 220 random (k, p) equal the installed mash 2.3.")]
    public static MashTextResult MashBounds(
        [Description("k-mer size (1 - 32, default 21).")] int k = 21,
        [Description("Confidence of the bounds (0 - 1 exclusive, default 0.99).")] double probability = 0.99)
        => new(global::Seqeron.Genomics.Analysis.MashCommands.Bounds(k, probability));

    [McpServerTool(Name = "mash_dist", Title = "Mash — Distances of Queries to References (mash dist)", ReadOnly = true)]
    [Description("Mash 2.3 'mash dist': MinHash (bottom-s) sketches of the reference and query files (canonical k-mers, MurmurHash3, the sketch of Ondov et al. 2016: KmerAnalyzer.CreateMinHashSketch) and, for every query in order and every reference, the line 'reference TAB query TAB distance TAB p-value TAB shared/total' of the pairs whose p-value is at most maxPValue (-v) and distance at most maxDistance (-d); table = true (-t) writes the '#query' header and a row of distances per query (a pair that fails a filter leaves its field blank); individual = true (-i) sketches every sequence on its own, named by its id (sequences shorter than k are skipped). A file without a record of at least k symbols is refused as the program refuses it (exit status 1: 'All fasta records in input files were shorter than the k-mer size'). Numbers are printed as C++ cout does (6 significant digits). The alphabet options -a (protein: amino acids, implies -n -k 9), -z (own alphabet, implies -n) and -Z (case kept) as in the program (k-mer space alphabetSize^k in the p-value, 64-bit hashes when alphabetSize^k > 2^32). VERIFICATION: 220 / 220 random sets (1 - 4 references x 1 - 3 queries of 1 - 3 records, k 4 - 32, sketch size 5 - 2000, seeds, -n, -i, -t, -v, -d) and 300 / 300 sets with -a / -z / -Z equal the installed mash 2.3 byte for byte.")]
    public static MashTextResult MashDist(
        [Description("The reference files as FASTA texts (records start with '>': the first word of the header is the sequence id, the rest the comment).")] string[] references,
        [Description("The query files (FASTA texts).")] string[] queries,
        [Description("Names of the reference files as the program prints them (default reference1, ...).")] string[]? referenceNames = null,
        [Description("Names of the query files (default query1, ...).")] string[]? queryNames = null,
        [Description("k-mer size (-k, 1 - 32; default 21, or 9 with protein).")] int? k = null,
        [Description("Sketch size (-s, default 1000).")] int sketchSize = 1000,
        [Description("Hash seed (-S, 0 - 4294967295, default 42).")] long seed = 42,
        [Description("-n: preserve the strand (no canonical k-mers).")] bool noncanonical = false,
        [Description("-i: sketch individual sequences instead of whole files.")] bool individual = false,
        [Description("-t: table output.")] bool table = false,
        [Description("-v: maximum p-value to report (0 - 1, default 1).")] double maxPValue = 1.0,
        [Description("-d: maximum distance to report (0 - 1, default 1).")] double maxDistance = 1.0,
        CancellationToken cancellationToken = default,
        [Description("-r: the inputs are read sets (noise filter -m, target coverage -c, genome size -g); incompatible with individual.")] bool reads = false,
        [Description("-m: minimum copies of each k-mer required to pass the noise filter of a read set (default 1; only with reads).")] int minCopies = 1,
        [Description("-c: target coverage of a read set: sketching concludes when it is reached (default 0 = read everything; only with reads).")] double targetCoverage = 0.0,
        [Description("-g: genome size used as the sketch length of a read set (default 0 = estimated; only with reads).")] long genomeSize = 0,
        [Description("-C: show the comment fields with the names (name:comment for dist, the comment instead of the id for triangle; the comment of a file sketch is the first record's name and comment, '[N seqs] ... [...]' for several).")] bool showComments = false,
        [Description("-b: a Bloom filter of this size (raw bytes or with K/M/G/T, e.g. 500M) filters out the k-mers seen once, instead of the exact filter -m (cannot be combined with minCopies > 1; implies reads; copies cannot be counted beyond 2; some unique k-mers pass erroneously; the filter has one hash function over bytes x 8 bits as the program builds it).")] string? bloomFilter = null,
        [Description("-a: the amino-acid alphabet (A-Z except BJOUXZ; implies -n and, without k, -k 9): k-mers of the 20 letters, 64-bit hashes when 20^k > 2^32, p-values over the k-mer space 20^k.")] bool protein = false,
        [Description("-z: the alphabet the hashes are based on (case ignored unless preserveCase; k-mers with other characters are ignored; implies -n; ASCII only: the program has no defined result for a byte of 0x80 or more). protein takes precedence.")] string? alphabet = null,
        [Description("-Z: preserve the case of the letters in the k-mers and the alphabet (a letter whose case is not in the alphabet is skipped).")] bool preserveCase = false)
        => new(global::Seqeron.Genomics.Analysis.MashCommands.Dist(ToMashInputs(references, referenceNames, "reference", nameof(references)), ToMashInputs(queries, queryNames, "query", nameof(queries)), MashK(k, protein), sketchSize, RequireMashSeed(seed),
            !noncanonical, individual, table, maxPValue, maxDistance, cancellationToken, ReadsOptions(reads, minCopies, targetCoverage, genomeSize, bloomFilter), showComments, MashAlphabetOf(protein, alphabet, preserveCase)));

    [McpServerTool(Name = "mash_triangle", Title = "Mash — Lower-Triangular Distance Matrix (mash triangle)", ReadOnly = true)]
    [Description("Mash 2.3 'mash triangle': the distance of each input to every other input as a lower-triangular matrix in relaxed Phylip format (a first line with a TAB and the number of inputs, then one line per input: its name and the distances to the inputs before it), or with edgeList = true (-E; implied by a p-value or distance filter) the edge list 'name TAB name TAB distance TAB p-value TAB shared/total' of the pairs that pass maxPValue (-v) and maxDistance (-d). Sketches as mash_dist (a file without a record of at least k symbols is refused). The alphabet options -a / -z / -Z as mash_dist. VERIFICATION: 220 / 220 random sets (2 - 6 files, k 4 - 32, sketch size 5 - 2000, -n, -i, -E, -v, -d) and the 300 / 300 sets with -a / -z / -Z (shared with mash_dist) equal the installed mash 2.3 byte for byte.")]
    public static MashTextResult MashTriangle(
        [Description("The input files (FASTA texts), at least two for a matrix.")] string[] sequences,
        [Description("Names of the inputs as the program prints them (default input1, ...).")] string[]? names = null,
        [Description("k-mer size (-k, 1 - 32; default 21, or 9 with protein).")] int? k = null,
        [Description("Sketch size (-s, default 1000).")] int sketchSize = 1000,
        [Description("Hash seed (-S, 0 - 4294967295, default 42).")] long seed = 42,
        [Description("-n: preserve the strand (no canonical k-mers).")] bool noncanonical = false,
        [Description("-i: sketch individual sequences instead of whole files.")] bool individual = false,
        [Description("-E: output the edge list instead of the matrix.")] bool edgeList = false,
        [Description("-v: maximum p-value to report in the edge list (0 - 1, default 1).")] double maxPValue = 1.0,
        [Description("-d: maximum distance to report in the edge list (0 - 1, default 1).")] double maxDistance = 1.0,
        CancellationToken cancellationToken = default,
        [Description("-r: the inputs are read sets (noise filter -m, target coverage -c, genome size -g); incompatible with individual.")] bool reads = false,
        [Description("-m: minimum copies of each k-mer required to pass the noise filter of a read set (default 1; only with reads).")] int minCopies = 1,
        [Description("-c: target coverage of a read set: sketching concludes when it is reached (default 0 = read everything; only with reads).")] double targetCoverage = 0.0,
        [Description("-g: genome size used as the sketch length of a read set (default 0 = estimated; only with reads).")] long genomeSize = 0,
        [Description("-C: show the comment fields with the names (name:comment for dist, the comment instead of the id for triangle; the comment of a file sketch is the first record's name and comment, '[N seqs] ... [...]' for several).")] bool showComments = false,
        [Description("-b: a Bloom filter of this size (raw bytes or with K/M/G/T, e.g. 500M) filters out the k-mers seen once, instead of the exact filter -m (cannot be combined with minCopies > 1; implies reads; copies cannot be counted beyond 2; some unique k-mers pass erroneously; the filter has one hash function over bytes x 8 bits as the program builds it).")] string? bloomFilter = null,
        [Description("-a: the amino-acid alphabet (A-Z except BJOUXZ; implies -n and, without k, -k 9): k-mers of the 20 letters, 64-bit hashes when 20^k > 2^32, p-values over the k-mer space 20^k.")] bool protein = false,
        [Description("-z: the alphabet the hashes are based on (case ignored unless preserveCase; k-mers with other characters are ignored; implies -n; ASCII only: the program has no defined result for a byte of 0x80 or more). protein takes precedence.")] string? alphabet = null,
        [Description("-Z: preserve the case of the letters in the k-mers and the alphabet (a letter whose case is not in the alphabet is skipped).")] bool preserveCase = false)
        => new(global::Seqeron.Genomics.Analysis.MashCommands.Triangle(ToMashInputs(sequences, names, "input", nameof(sequences)), MashK(k, protein), sketchSize, RequireMashSeed(seed), !noncanonical, individual, edgeList, maxPValue, maxDistance, cancellationToken, ReadsOptions(reads, minCopies, targetCoverage, genomeSize, bloomFilter), showComments, MashAlphabetOf(protein, alphabet, preserveCase)));

    [McpServerTool(Name = "mash_screen", Title = "Mash — Containment of Queries in a Mixture (mash screen)", ReadOnly = true)]
    [Description("Mash 2.3 'mash screen', byte for byte (CommandScreen.cpp): the MinHash sketches of the query files (as mash sketch -k -s -S -n -i would write them into the query .msh) are screened against a mixture (reads or contigs, any number of files): the hashes of the queries are counted in every window of the mixture and for each query one line identity TAB shared/hashes TAB median-multiplicity TAB p-value TAB query-ID TAB query-comment is written, identity = (shared / hashes)^(1/k), the median multiplicity over the shared hashes, the p-value the binomial tail with r = estimated distinct k-mers of the mixture / 4^k. winnerTakesAll = -w (a hash held by several queries counts for the query with the best identity, ties the longer one), minIdentity = -i (-1 writes every query, 0 those with a shared hash), maxPValue = -v. The query-comment is 'name comment' of the first record of a file ('[N seqs] ... [...]' for several records; the record's own comment with individual). Amino-acid queries (protein = -a) translate the mixture in six frames as the program does (the frames count from the start of the 2^20-character chunk the records of all files are appended to, one '*' in front of each record), other alphabets (-z, -Z) screen the mixture as nucleotides with their own alphabet; the table of the mixture upper-cases the alphabet whatever -Z says (a -z abc -Z sketch screens nothing). VERIFICATION: 240 / 240 random nucleotide screens and 300 / 300 screens of -a / -z / -Z query sketches equal the installed mash 2.3.")]
    public static MashTextResult MashScreen(
        [Description("The query files as FASTA texts (records start with '>': the first word of the header is the sequence id, the rest the comment).")] string[] queries,
        [Description("The mixture files as FASTA texts (reads or contigs).")] string[] mixtures,
        [Description("Names of the query files as the program prints them (default query1, ...).")] string[]? queryNames = null,
        [Description("k-mer size of the query sketch (-k of mash sketch, 1 - 32; default 21, or 9 with protein).")] int? k = null,
        [Description("Sketch size of the query sketch (-s, default 1000).")] int sketchSize = 1000,
        [Description("Hash seed of the query sketch (-S, 0 - 4294967295, default 42).")] long seed = 42,
        [Description("The query sketch was made with -n (strand preserved).")] bool noncanonical = false,
        [Description("The query sketch was made with -i (one reference per sequence).")] bool individual = false,
        [Description("-w: winner-takes-all strategy for identity estimates.")] bool winnerTakesAll = false,
        [Description("-i: minimum identity to report (-1 - 1, default 0).")] double minIdentity = 0.0,
        [Description("-v: maximum p-value to report (0 - 1, default 1).")] double maxPValue = 1.0,
        CancellationToken cancellationToken = default,
        [Description("-a: the amino-acid alphabet (A-Z except BJOUXZ; implies -n and, without k, -k 9): k-mers of the 20 letters, 64-bit hashes when 20^k > 2^32, p-values over the k-mer space 20^k.")] bool protein = false,
        [Description("-z: the alphabet the hashes are based on (case ignored unless preserveCase; k-mers with other characters are ignored; implies -n; ASCII only: the program has no defined result for a byte of 0x80 or more). protein takes precedence.")] string? alphabet = null,
        [Description("-Z: preserve the case of the letters in the k-mers and the alphabet (a letter whose case is not in the alphabet is skipped).")] bool preserveCase = false)
        => new(global::Seqeron.Genomics.Analysis.MashCommands.Screen(ToMashInputs(queries, queryNames, "query", nameof(queries)), ToMashInputs(mixtures, null, "mixture", nameof(mixtures)), MashK(k, protein), sketchSize, RequireMashSeed(seed),
            !noncanonical, individual, winnerTakesAll, minIdentity, maxPValue, cancellationToken, MashAlphabetOf(protein, alphabet, preserveCase)));

    [McpServerTool(Name = "mash_taxscreen", Title = "Mash — Taxonomic Screen (mash taxscreen)", ReadOnly = true)]
    [Description("Mash 2.3 'mash taxscreen', byte for byte (CommandTaxScreen.cpp + taxdb.hpp): the Kraken-style taxonomic report of how well query references are contained in a pool of sequences. The pass over the pool is that of mash_screen (reads or contigs; 6-frame translation for amino-acid queries); the taxonomy ID of every query reference is the one of the mapping file (lines 'taxid SEP reference-name') or the number after the word 'taxid' in its comment, the hashes of references that share a hash are assigned to their lowest common ancestor in the taxonomy (names.dmp / nodes.dmp, the NCBI taxonomy dump), and the report lists from the root: percent of the hashes, number of contained hashes in the clade, in the taxon, number of hashes of the references in the clade and the taxon, rank, taxonomy ID and the name indented by two blanks per level, children by decreasing clade count. The report is the program's, including its order of equal counts (std::sort of libstdc++) and the effect of its unordered_map growing while the clade counts are summed. Reads the dump from the texts given (no download: a taxonomy of your own in the dump format works). A node whose parent is missing from nodes.dmp is refused (the program uses an uninitialised pointer there: valgrind-verified). VERIFICATION: 300 / 300 random runs (synthetic taxonomies of 4 - 300 taxa in the NCBI dump format: random, deep and star-shaped trees with up to 45 equal children, nucleotide and amino-acid queries, mapping files, taxonomy IDs in comments, IDs missing from the dump) equal the installed mash 2.3.")]
    public static MashTextResult MashTaxScreen(
        [Description("The query files as FASTA texts, or one sketch file 'msh:' + base64 (the output of mash_sketch); the reference names are the file names (or the sequence ids with individual).")] string[] queries,
        [Description("The pool files as FASTA texts (reads or contigs, nucleotides).")] string[] pools,
        [Description("The content of names.dmp (tax_id | name | unique name | name class |, one taxon per line).")] string namesDump,
        [Description("The content of nodes.dmp (tax_id | parent tax_id | rank | ...).")] string nodesDump,
        [Description("-m: the mapping file ('taxonomy ID, one separator character, reference name' per line); null = the IDs come from the comments of the references ('taxid N').")] string? mappingFile = null,
        [Description("Names of the query files (default query1, ...).")] string[]? queryNames = null,
        [Description("k-mer size of the query sketch (-k of mash sketch, 1 - 32; default 21, or 9 with protein).")] int? k = null,
        [Description("Sketch size of the query sketch (-s, default 1000).")] int sketchSize = 1000,
        [Description("Hash seed of the query sketch (-S, 0 - 4294967295, default 42).")] long seed = 42,
        [Description("The query sketch was made with -n (strand preserved).")] bool noncanonical = false,
        [Description("The query sketch was made with -i (one reference per sequence).")] bool individual = false,
        CancellationToken cancellationToken = default,
        [Description("-a: the amino-acid alphabet of the query sketch (implies -n and, without k, -k 9; the pool is translated in six frames).")] bool protein = false,
        [Description("-z: the alphabet of the query sketch (implies -n).")] string? alphabet = null,
        [Description("-Z: the case of the letters is preserved in the query sketch.")] bool preserveCase = false)
        => new(global::Seqeron.Genomics.Analysis.MashTaxScreen.Run(ToMashInputs(queries, queryNames, "query", nameof(queries)), ToMashInputs(pools, null, "pool", nameof(pools)),
            new global::Seqeron.Genomics.Analysis.MashTaxonomy(namesDump ?? throw new ArgumentException("namesDump is required", nameof(namesDump)), nodesDump ?? throw new ArgumentException("nodesDump is required", nameof(nodesDump)), mappingFile),
            MashK(k, protein), sketchSize, RequireMashSeed(seed), !noncanonical, individual, MashAlphabetOf(protein, alphabet, preserveCase), cancellationToken));

    [McpServerTool(Name = "mash_sketch", Title = "Mash — Sketch File (mash sketch -o)", ReadOnly = true)]
    [Description("Mash 2.3 'mash sketch -o x.msh': the sketch file of the inputs, byte for byte as the program writes it (120 random runs equal, and 300 more with the alphabet options -a -z -Z: genomes and read sets, k 4 - 32, sketch sizes 1 - 3000, seeds, -n, -i, -r -m -c -g with all the files of a read set in ONE sketch, -I / -C, -l). The file is a Cap'n Proto message of MinHash.capnp built as capnp::MallocMessageBuilder builds it (segments of 1024, 1024, 2048, ... words, far pointers with landing pads); it is returned as base64 and, prefixed with 'msh:', it is accepted by mash_dist / mash_triangle / mash_screen / mash_info / mash_paste as an input. One sketch per file (its name is the file name, its comment the first record's name and comment) or, with individual, one per sequence. -l (a list of file names) is the caller's list of inputs.")]
    public static MashSketchFileResult MashSketch(
        [Description("The input files as FASTA texts (records start with '>': the first word of the header is the sequence id, the rest the comment).")] string[] sequences,
        [Description("Names of the inputs as the program stores them (the file names; default input1, ...).")] string[]? names = null,
        [Description("k-mer size (-k, 1 - 32; default 21, or 9 with protein).")] int? k = null,
        [Description("Sketch size (-s, default 1000).")] int sketchSize = 1000,
        [Description("Hash seed (-S, 0 - 4294967295, default 42; 42 is stored in referenceListOld, another seed in referenceList).")] long seed = 42,
        [Description("-n: preserve the strand (no canonical k-mers).")] bool noncanonical = false,
        [Description("-i: sketch individual sequences instead of whole files.")] bool individual = false,
        [Description("-r: the inputs are read sets (noise filter -m, target coverage -c, genome size -g); the hash counts are stored.")] bool reads = false,
        [Description("-m: minimum copies of each k-mer required to pass the noise filter of a read set (default 1; only with reads).")] int minCopies = 1,
        [Description("-c: target coverage of a read set (default 0 = read everything; only with reads).")] double targetCoverage = 0.0,
        [Description("-g: genome size used as the sketch length of a read set (default 0 = estimated; only with reads).")] long genomeSize = 0,
        [Description("-I: the id (name) of the first sketch.")] string? id = null,
        [Description("-C: the comment of the first sketch.")] string? comment = null,
        CancellationToken cancellationToken = default,
        [Description("-b: a Bloom filter of this size (raw bytes or with K/M/G/T, e.g. 500M) filters out the k-mers seen once, instead of the exact filter -m (cannot be combined with minCopies > 1; implies reads; copies cannot be counted beyond 2; some unique k-mers pass erroneously; the filter has one hash function over bytes x 8 bits as the program builds it).")] string? bloomFilter = null,
        [Description("-a: the amino-acid alphabet (A-Z except BJOUXZ; implies -n and, without k, -k 9): k-mers of the 20 letters, 64-bit hashes when 20^k > 2^32, p-values over the k-mer space 20^k.")] bool protein = false,
        [Description("-z: the alphabet the hashes are based on (case ignored unless preserveCase; k-mers with other characters are ignored; implies -n; ASCII only: the program has no defined result for a byte of 0x80 or more). protein takes precedence.")] string? alphabet = null,
        [Description("-Z: preserve the case of the letters in the k-mers and the alphabet (a letter whose case is not in the alphabet is skipped).")] bool preserveCase = false)
    {
        var data = global::Seqeron.Genomics.Analysis.MashCommands.Sketch(ToMashInputs(sequences, names, "input", nameof(sequences)), MashK(k, protein), sketchSize, RequireMashSeed(seed), !noncanonical, individual, cancellationToken,
            ReadsOptions(reads, minCopies, targetCoverage, genomeSize, bloomFilter), id, comment, MashAlphabetOf(protein, alphabet, preserveCase));
        return MshResult(data);
    }

    private static MashSketchFileResult MshResult(global::Seqeron.Genomics.Analysis.MashSketchFileData data)
    {
        var bytes = global::Seqeron.Genomics.Analysis.MashSketchFile.Write(data);
        return new MashSketchFileResult(MshPrefix + Convert.ToBase64String(bytes), bytes.Length, data.References.Count);
    }

    [McpServerTool(Name = "mash_info", Title = "Mash — Sketch File Information (mash info)", ReadOnly = true)]
    [Description("Mash 2.3 'mash info', byte for byte (120 sketch files x 5 output modes equal): the header (hash function and seed, k-mer size with the hash width, the alphabet with '(canonical)' / '(case-sensitive)', the target sketch size, the number of sketches) and the table of the sketches (hashes, length, id, comment; an empty field is '-'), or with headerOnly only the header, with tabular the tab separated table, with counts the hash count histogram of each read-set sketch, with dump the JSON dump of the sketches and hashes. The options headerOnly / tabular / counts / dump exclude each other as in the program (the errors are the program's messages).")]
    public static MashTextResult MashInfo(
        [Description("The sketch file: 'msh:' + base64 (the output of mash_sketch / mash_paste).")] string sketch,
        [Description("-H: only the header.")] bool headerOnly = false,
        [Description("-t: tabular output (no header).")] bool tabular = false,
        [Description("-c: hash count histograms (only for read-set sketches).")] bool counts = false,
        [Description("-d: dump the sketches as JSON.")] bool dump = false)
    {
        ArgumentException.ThrowIfNullOrEmpty(sketch);
        var data = ReadMshText(sketch.StartsWith(MshPrefix, StringComparison.Ordinal) ? sketch : MshPrefix + sketch, nameof(sketch));
        return new(global::Seqeron.Genomics.Analysis.MashCommands.Info(data, new global::Seqeron.Genomics.Analysis.MashCommands.MashInfoOptions(headerOnly, tabular, counts, dump)));
    }

    [McpServerTool(Name = "mash_paste", Title = "Mash — Paste Sketch Files (mash paste)", ReadOnly = true)]
    [Description("Mash 2.3 'mash paste', byte for byte (35 random runs equal): one sketch file of the sketches of several. The parameters (k, sketch size, seed, alphabet, strand mode, case) are those of the first file; a later file with another alphabet, seed, k or strand mode or a smaller sketch size is skipped (the program warns on the error stream), a larger sketch is reduced to the first file's size; hash counts of read-set sketches are not kept.")]
    public static MashSketchFileResult MashPaste(
        [Description("The sketch files: 'msh:' + base64 each.")] string[] sketches)
    {
        if (sketches is null || sketches.Length == 0)
            throw new ArgumentException("At least one sketch file is required", nameof(sketches));
        var files = sketches.Select(t => ReadMshText(t.StartsWith(MshPrefix, StringComparison.Ordinal) ? t : MshPrefix + t, nameof(sketches))).ToArray();
        return MshResult(global::Seqeron.Genomics.Analysis.MashCommands.Paste(files));
    }
}
