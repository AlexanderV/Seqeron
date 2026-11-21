using System.ComponentModel;
using ModelContextProtocol.Server;
using global::Seqeron.Genomics.Alignment;
using global::Seqeron.Genomics.Analysis;
using Seqeron.Mcp.Analysis.Models;

namespace Seqeron.Mcp.Analysis.Tools;

// MUMmer 3.23 programs on the Analysis server (B09 A3-9b); every tool delegates to GenomicAnalyzer / the ports in Seqeron.Genomics.Alignment.
public partial class AnalysisTools
{
    [McpServerTool(Name = "nucmer", Title = "Genome — MUMmer nucmer (delta alignments)", ReadOnly = true)]
    [Description("MUMmer 3.23 nucmer, the whole pipeline: prenuc (the reference records joined with an x after each) | mummer (anchors: matchMode 'mum' = unique in reference and query, 'mumreference' = unique in the reference (default), 'maxmatch' = all; minMatch = -l; only A/C/G/T match; both strands unless forwardOnly / reverseOnly) | mgaps (clusters: maxGap = -g, minCluster = -c, diagonalDifference = -D, diagonalFactor = -d) | postnuc (the banded Smith-Waterman extension of the clusters: gaps between the matches of a cluster, then outward off both ends until the score has not improved for breakLength = -b bases; clusters that meet fuse; with simplify (default) clusters shadowed by an earlier alignment are dropped; extend / optimize / banded as --[no]extend / --[no]optimize (nooptimize forces the extension to the sequence ends) / --banded (hard diagonal band = diagonalDifference)). Returns the text of the .delta file - the alignments between every reference record and every query record as the distances between insertions and deletions (first line = the two file names, then NUCMER, then per record pair '>refId qryId refLen qryLen' and per alignment 'sR eR sQ eQ errors simErrors nonAlpha' followed by the delta values and 0; a reverse-strand alignment has sQ > eQ) - byte for byte what the installed nucmer 3.23 writes (434 runs: multi-record inputs, both strands, repeats, N runs, every option family, the abort paths). With delta = false (--nodelta) the .cluster text is returned instead. A break length outside 1 - 10000 is ignored with a warning (messages), as the program does. The inputs on which the programs abort (empty file, no header, a file whose only record is empty, an unparsable header) are refused with an error.")]
    public static NucmerToolResult Nucmer(
        [Description("Reference FASTA text (multi-record).")] string reference,
        [Description("Query FASTA text (multi-record).")] string query,
        [Description("Anchor matches: 'mumreference' (default, --mumreference / --mumcand), 'mum' (--mum) or 'maxmatch' (--maxmatch).")] string matchMode = "mumreference",
        [Description("Minimum length of a single match (-l, default 20).")] int minMatch = 20,
        [Description("Maximum gap between two adjacent matches of a cluster (-g, default 90).")] long maxGap = 90,
        [Description("Minimum length of a cluster of matches (-c, default 65).")] long minCluster = 65,
        [Description("Maximum diagonal difference of two adjacent anchors (-D, default 5).")] int diagonalDifference = 5,
        [Description("Diagonal difference as a fraction of the gap length (-d, default 0.12).")] double diagonalFactor = 0.12,
        [Description("Distance an extension tries to extend poor scoring regions before giving up (-b, default 200; outside 1 - 10000 ignored with a warning).")] int breakLength = 200,
        [Description("Enforce absolute banding of the dynamic programming matrix with diagonalDifference (--banded, default false).")] bool banded = false,
        [Description("Use only the forward strand of the query (-f, default false).")] bool forwardOnly = false,
        [Description("Use only the reverse complement of the query (-r, default false); both flags = both strands.")] bool reverseOnly = false,
        [Description("Extend the clusters (--extend, default true).")] bool extend = true,
        [Description("Create the delta (--delta, default true); false returns the cluster file (--nodelta).")] bool delta = true,
        [Description("Alignment score optimisation (--optimize, default true); false = --nooptimize: an extension that reaches the end of a sequence is not backtracked.")] bool optimize = true,
        [Description("Remove shadowed alignments (--simplify, default true); false = --nosimplify (align a sequence to itself to find repeats).")] bool simplify = true,
        [Description("Reference file name printed in the first line of the output (nucmer prints the absolute path).")] string referenceName = "reference.fa",
        [Description("Query file name printed in the first line of the output.")] string queryName = "query.fa",
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(query));
        var mode = (matchMode ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "mumreference" or "mumcand" => CommonRegionMode.MaximalUniqueReference,
            "mum" => CommonRegionMode.MaximalUnique,
            "maxmatch" => CommonRegionMode.MaximalExactMatches,
            _ => throw new ArgumentException("matchMode must be 'mumreference', 'mum' or 'maxmatch'", nameof(matchMode)),
        };
        var r = GenomicAnalyzer.Nucmer(reference, query, referenceName, queryName, new NucmerOptions
        {
            MatchMode = mode, MinMatch = minMatch, MaxGap = maxGap, MinCluster = minCluster, DiagonalDifference = diagonalDifference, DiagonalFactor = diagonalFactor, BreakLength = breakLength,
            Banded = banded, ForwardOnly = forwardOnly, ReverseOnly = reverseOnly, Extend = extend, Delta = delta, Optimize = optimize, Simplify = simplify,
        }, cancellationToken);
        return new NucmerToolResult(r.Delta, r.Cluster, r.Messages);
    }

    [McpServerTool(Name = "delta_filter", Title = "Genome — MUMmer delta-filter", ReadOnly = true)]
    [Description("MUMmer 3.23 delta-filter: reads a delta file (nucmer / promer) and keeps the alignments that survive the switches, in the order -i -l -u -q -r -g -m -1; the output is a delta file again (grouped by query record). minIdentity (-i, 0 - 100) and minLength (-l) drop weak alignments, minUnique (-u) those whose reference AND query part is not unique; queryLis (-q) maps each query position to its best hit in the reference (allowing reference overlaps), referenceLis (-r) each reference position to its best hit, globalLis (-g) the longest set of mutually consistent alignments (no rearrangements), manyToMany (-m) the union and oneToOne (-1) the intersection of the -r and -q chains; maxOverlap (-o) is the largest overlap of two chained alignments as a percent of their length; epsilon (-e, default -1) is the negligible score loss: a negative value breaks ties among equally good chains with glibc rand() seeded with 1, exactly as the program (reproducible). Byte for byte equal to the installed delta-filter on 300 + runs of random option sets over delta files of rearranged, duplicated and inverted multi-contig inputs. A percentage outside [0, 100] or a negative length is an error (the program's range errors).")]
    public static MummerTextResult DeltaFilter(
        [Description("Delta file text (the .delta output of nucmer).")] string delta,
        [Description("-q: query-based longest increasing subset (default false).")] bool queryLis = false,
        [Description("-r: reference-based longest increasing subset (default false).")] bool referenceLis = false,
        [Description("-g: 1-to-1 global alignment not allowing rearrangements (default false).")] bool globalLis = false,
        [Description("-1: 1-to-1 alignment allowing rearrangements = intersection of -r and -q (default false).")] bool oneToOne = false,
        [Description("-m: many-to-many alignment allowing rearrangements = union of -r and -q (default false).")] bool manyToMany = false,
        [Description("-l: minimum alignment length (default 0).")] long minLength = 0,
        [Description("-i: minimum alignment identity, 0 - 100 (default 0).")] float minIdentity = 0,
        [Description("-u: minimum alignment uniqueness, 0 - 100 (default 0).")] float minUnique = 0,
        [Description("-o: maximum alignment overlap for -r and -q as a percent of the alignment length, 0 - 100 (default 100).")] float maxOverlap = 100,
        [Description("-e: negligible alignment score; negative (default -1) picks randomly among equal chains with the program's rand().")] float epsilon = -1)
    {
        if (string.IsNullOrEmpty(delta))
            throw new ArgumentException("Delta file text cannot be null or empty", nameof(delta));
        return new MummerTextResult(MummerDeltaFilter.Run(delta, new MummerDeltaFilterOptions
        {
            QueryLis = queryLis, ReferenceLis = referenceLis, GlobalLis = globalLis, OneToOne = oneToOne, ManyToMany = manyToMany, MinLength = minLength, MinIdentity = minIdentity, MinUnique = minUnique, MaxOverlap = maxOverlap, Epsilon = epsilon,
        }));
    }

    [McpServerTool(Name = "show_coords", Title = "Genome — MUMmer show-coords", ReadOnly = true)]
    [Description("MUMmer 3.23 show-coords: the coordinates of the alignments of a delta file as the human-readable table (header, [S1] [E1] | [S2] [E2] | [LEN 1] [LEN 2] | [% IDY] | tags), the tab-delimited table (tabular, -T) or btab (btab, -B). header = false is -H; brief (-b) merges overlapping alignments regardless of direction and drops the identity columns; coverage (-c) adds the percent coverage columns, lengths (-l) the sequence lengths, direction (-d) the FRM columns; minIdentity (-I) and minLength (-L) filter; sortByReference (-r) / sortByQuery (-q) sort; annotateOverlaps (-o) tags alignments that cover a whole sequence end ([IDENTITY] [CONTAINED] [CONTAINS] [END] [BEGIN]); the hidden warnings (-w) tag [DUPLICATE] [CONTAINS] [SHADOWED] [OVERLAPS]; las (-g) / weightedLas (-G) keep the longest ascending subset (ties broken with rand() seeded by randomSeed, the program seeds it with the current time); knockout (-k) is for promer files. Percentages are printed from the float values with %.2f. Byte for byte equal to the installed show-coords on 450 + runs; PROmer delta files are accepted.")]
    public static MummerProgramResult ShowCoords(
        [Description("Delta file text.")] string delta,
        [Description("Print the header (default true; false = -H).")] bool header = true,
        [Description("-b: merge overlapping alignments, no identity columns (default false).")] bool brief = false,
        [Description("-B: btab output (default false).")] bool btab = false,
        [Description("-c: percent coverage columns (default false).")] bool coverage = false,
        [Description("-d: alignment direction in the FRM columns (default false; always for promer).")] bool direction = false,
        [Description("-I: minimum alignment identity, 0 - 100 (default 0).")] float minIdentity = 0,
        [Description("-k: knockout (promer data only, default false).")] bool knockout = false,
        [Description("-l: sequence length columns (default false).")] bool lengths = false,
        [Description("-L: minimum alignment length (default 0).")] long minLength = 0,
        [Description("-o: annotate overlaps (default false).")] bool annotateOverlaps = false,
        [Description("-q: sort by query (default false).")] bool sortByQuery = false,
        [Description("-r: sort by reference (default false).")] bool sortByReference = false,
        [Description("-g: longest ascending subset, deprecated (default false).")] bool las = false,
        [Description("-G: weighted longest ascending subset, hidden (default false).")] bool weightedLas = false,
        [Description("-T: tab-delimited output (default false).")] bool tabular = false,
        [Description("-w: warnings column, hidden (default false).")] bool warnings = false,
        [Description("Seed of the rand() that breaks ties of -g / -G (the program uses the current time); default = the current time.")] uint? randomSeed = null)
    {
        if (string.IsNullOrEmpty(delta))
            throw new ArgumentException("Delta file text cannot be null or empty", nameof(delta));
        var r = MummerShowCoords.Run(delta, new MummerShowCoordsOptions
        {
            PrintHeader = header, Brief = brief, Btab = btab, ShowCoverage = coverage, ShowDirection = direction, MinIdentity = minIdentity, Knockout = knockout, ShowSeqLengths = lengths, MinLength = minLength,
            AnnotateOverlaps = annotateOverlaps, SortByQuery = sortByQuery, SortByReference = sortByReference, Las = las, WeightedLas = weightedLas, Tabular = tabular, ShowWarnings = warnings, RandomSeed = randomSeed,
        });
        return new MummerProgramResult(r.Text, r.Messages);
    }

    [McpServerTool(Name = "show_aligns", Title = "Genome — MUMmer show-aligns", ReadOnly = true)]
    [Description("MUMmer 3.23 show-aligns: the alignments between one reference record and one query record of a nucmer delta file drawn as three rows (reference, query, a ^ under every difference; an indel is a . in the other row) in blocks of the screen width (-w, default 60; 11 or less is replaced with a warning), each row prefixed with the 1-based forward-strand position of its first base; '-- BEGIN alignment [ +1 sR - eR | +1 sQ - eQ ]' / END lines frame every alignment (a reverse-strand alignment shows -1). The alignments are in file order unless sortByReference (-r) / sortByQuery (-q). The program reads the sequences from the FASTA files named in the delta header: pass their text. Byte for byte equal to the installed show-aligns on 300 + runs; a promer (amino-acid) delta file is refused with NotSupportedException (registered A3-9d).")]
    public static MummerProgramResult ShowAligns(
        [Description("Delta file text (nucmer).")] string delta,
        [Description("Reference FASTA text (the file named in the delta header).")] string referenceFasta,
        [Description("Query FASTA text (the file named in the delta header).")] string queryFasta,
        [Description("Reference record id.")] string referenceId,
        [Description("Query record id.")] string queryId,
        [Description("-q: sort the alignments by the query start (default false).")] bool sortByQuery = false,
        [Description("-r: sort the alignments by the reference start (default false).")] bool sortByReference = false,
        [Description("-w: screen width (default 60).")] int screenWidth = 60)
    {
        if (string.IsNullOrEmpty(delta))
            throw new ArgumentException("Delta file text cannot be null or empty", nameof(delta));
        var r = MummerShowAligns.Run(delta, referenceFasta, queryFasta, referenceId, queryId, new MummerShowAlignsOptions { SortByQuery = sortByQuery, SortByReference = sortByReference, ScreenWidth = screenWidth });
        return new MummerProgramResult(r.Text, r.Messages);
    }

    [McpServerTool(Name = "show_snps", Title = "Genome — MUMmer show-snps", ReadOnly = true)]
    [Description("MUMmer 3.23 show-snps: the single-base substitutions and indels of every alignment of a nucmer delta file: [P1] reference position, [SUB] the two characters ('.' = indel), [P2] query position, [BUFF] the distance to the nearest other difference in the same alignment, [DIST] the distance to the nearest sequence end, [R] [Q] the number of other alignments over the position (showConflict = false is -C: only unambiguous SNPs and no [R] [Q] columns), optional sequence lengths (-l), context (-x characters each side), frames and ids; a difference between compatible IUPAC codes (N, R/A ...) is not a SNP (the program's asymmetric CompareIUPAC table). sortByReference (-r, default) / sortByQuery (-q), tabular (-T), header = false (-H), showIndels = false (-I), alignments (-S): only the alignments of these show-coords lines. Byte for byte equal to the installed show-snps on 450 + runs; a promer delta file is refused with NotSupportedException (registered A3-9d).")]
    public static MummerProgramResult ShowSnps(
        [Description("Delta file text (nucmer).")] string delta,
        [Description("Reference FASTA text (the file named in the delta header).")] string referenceFasta,
        [Description("Query FASTA text (the file named in the delta header).")] string queryFasta,
        [Description("Report SNPs of ambiguous alignments and show the [R] [Q] columns (default true; false = -C).")] bool showConflict = true,
        [Description("Report indels (default true; false = -I).")] bool showIndels = true,
        [Description("-l: sequence length columns (default false).")] bool lengths = false,
        [Description("-q: sort by query ids and positions (default false; -r wins when both are set).")] bool sortByQuery = false,
        [Description("-r: sort by reference ids and positions (the default when neither is set).")] bool sortByReference = false,
        [Description("-T: tab-delimited output (default false).")] bool tabular = false,
        [Description("Print the header (default true; false = -H).")] bool header = true,
        [Description("-x: characters of SNP context on each side (default 0).")] int context = 0,
        [Description("-S: report only the alignments of these show-coords lines (null = all).")] string[]? alignments = null)
    {
        if (string.IsNullOrEmpty(delta))
            throw new ArgumentException("Delta file text cannot be null or empty", nameof(delta));
        var r = MummerShowSnps.Run(delta, referenceFasta, queryFasta, new MummerShowSnpsOptions
        {
            ShowConflict = showConflict, ShowIndels = showIndels, ShowLength = lengths, SortQuery = sortByQuery, SortReference = sortByReference, Tabular = tabular, PrintHeader = header, Context = context, SelectedAlignments = alignments,
        });
        return new MummerProgramResult(r.Text, r.Messages);
    }

    [McpServerTool(Name = "show_diff", Title = "Genome — MUMmer show-diff", ReadOnly = true)]
    [Description("MUMmer 3.23 show-diff: the structural differences between the sequences of a delta file read off the many-to-many mapping (delta-filter -m) along every reference record (referenceDiff, -r, the default) or query record (queryDiff, -q; not both): BRK breaks (unaligned stretches), GAP (gap sizes in both sequences and their difference), JMP (a jump within a sequence), INV (an inversion), DUP (a duplication) and SEQ (a jump to another sequence, with the two neighbour ids); columns [SEQ] [TYPE] [S1] [E1] [LEN 1] (+ extra fields per type); amos (-f) writes the AMOS layout with the coordinates repeated. header = false is -H. The ties of the longest-increasing-subset choice use the program's rand() (seed 1). Byte for byte equal to the installed show-diff on 300 + runs.")]
    public static MummerProgramResult ShowDiff(
        [Description("Delta file text.")] string delta,
        [Description("-f: AMOS output (default false).")] bool amos = false,
        [Description("-r: differences of the reference (default; the default when neither -r nor -q is set).")] bool referenceDiff = false,
        [Description("-q: differences of the query (cannot be combined with -r).")] bool queryDiff = false,
        [Description("Print the header (default true; false = -H).")] bool header = true)
    {
        if (string.IsNullOrEmpty(delta))
            throw new ArgumentException("Delta file text cannot be null or empty", nameof(delta));
        var r = MummerShowDiff.Run(delta, new MummerShowDiffOptions { Amos = amos, ReferenceDiff = referenceDiff, QueryDiff = queryDiff, PrintHeader = header });
        return new MummerProgramResult(r.Text, r.Messages);
    }

    [McpServerTool(Name = "show_tiling", Title = "Genome — MUMmer show-tiling", ReadOnly = true)]
    [Description("MUMmer 3.23 show-tiling: places the aligned query contigs of a delta file on the reference. For every contig and reference the longest consistent subset of its alignments (maxGap -g, default 1000 for nucmer files / -1 for promer) is taken; a contig is tiled when its coverage is at least minCoverage (-v, default 95) and better than the second best by minCoverageDifference (-V, default 10) and its identity at least minIdentity (-i, default 90), or, when it hangs over the reference ends, by the same rule on the overhang-corrected coverage; contigs shorter than minContigLength (-l) are ignored. The tiled contigs are linked along each reference (circular = -c lets them span the origin) and printed as '>refId N bases' followed by 'start end gapToNext length coverage identity strand id' (the tiling path), or with printAlignments (-a) as the tab-delimited alignment rows, or with xml (-x) as the XML link evidence. contigList (-t), unusable (-u) and pseudoMolecule (-p, needs the query FASTA) return the program's files. randomRepeats (-R) places repetitive contigs at a random copy (rand() seeded by randomSeed; the program uses the current time; implies minCoverageDifference 0). Byte for byte equal to the installed show-tiling on 350 + runs.")]
    public static ShowTilingToolResult ShowTiling(
        [Description("Delta file text.")] string delta,
        [Description("Query FASTA text (only needed for pseudoMolecule).")] string? queryFasta = null,
        [Description("-a: print the alignment rows of the tiling path (default false).")] bool printAlignments = false,
        [Description("-c: circular reference sequences (default false).")] bool circular = false,
        [Description("-g: maximum gap between clustered alignments, -1 = infinity (nucmer default 1000, promer -1).")] long? maxGap = null,
        [Description("-i: minimum percent identity to tile, 0 - 100 (nucmer default 90, promer 55).")] float? minIdentity = null,
        [Description("-l: minimum contig length to report (default 1).")] long minContigLength = 1,
        [Description("-R: place repetitive contigs at a random one of their copies (default false).")] bool randomRepeats = false,
        [Description("-v: minimum contig coverage to tile, 0 - 100 (nucmer default 95, promer 50).")] float? minCoverage = null,
        [Description("-V: minimum contig coverage difference, 0 - 100 (nucmer default 10, promer 30).")] float? minCoverageDifference = null,
        [Description("-x: print the XML contig linking information (default false).")] bool xml = false,
        [Description("-t: also return the TIGR-style contig list (default false).")] bool contigList = false,
        [Description("-u: also return the alignment rows of the unusable contigs (default false).")] bool unusable = false,
        [Description("-p: also return the pseudo molecule of the query contigs (default false).")] bool pseudoMolecule = false,
        [Description("Seed of rand() for -R (the program uses the current time); default = the current time.")] uint? randomSeed = null)
    {
        if (string.IsNullOrEmpty(delta))
            throw new ArgumentException("Delta file text cannot be null or empty", nameof(delta));
        var r = MummerShowTiling.Run(delta, queryFasta, new MummerShowTilingOptions
        {
            PrintAlignments = printAlignments, Circular = circular, MaxGap = maxGap, MinIdentity = minIdentity, MinContigLength = minContigLength, RandomRepeats = randomRepeats, MinCoverage = minCoverage,
            MinCoverageDifference = minCoverageDifference, Xml = xml, ContigList = contigList, Unusable = unusable, PseudoMolecule = pseudoMolecule, RandomSeed = randomSeed,
        });
        return new ShowTilingToolResult(r.Text, r.ContigList, r.Unusable, r.PseudoMolecule);
    }

    [McpServerTool(Name = "dnadiff", Title = "Genome — MUMmer dnadiff", ReadOnly = true)]
    [Description("MUMmer 3.23 dnadiff: the comparative analysis of two sequence sets, the perl script over nucmer --maxmatch, delta-filter -1 / -m, show-coords -rclTH, show-snps -rlTHC and show-diff -rH / -qH, with its report: sequences, aligned / unaligned sequences and bases, the 1-to-1 and M-to-M alignments (count, total and average length, average identity), breakpoints, relocations, translocations, inversions, insertions, tandem insertions (count, sum, average) and the SNPs / indels with the 12 substitution classes, all also for the 'good' ones (buffer >= 20). Returns every file of the program: report, delta, 1delta, mdelta, 1coords, mcoords, snps, rdiff, qdiff and (when something is unaligned) unref / unqry. Pass delta (dnadiff -d) to analyse a precomputed nucmer delta file instead of running nucmer. The order of the per-substitution lines inside the SNP / indel blocks of the report and of the unref / unqry lines is the order of perl's hash iteration, which differs from run to run in the program (6 runs of the installed dnadiff on one input gave 6 different reports); here it is fixed (. A C G T, then other symbols). Everything else is byte for byte the program's (300 + runs, the unordered lines compared as sets). A FASTA line with white space (carriage return) is the program's 'Whitespace found in FastA' error.")]
    public static DnadiffToolResult Dnadiff(
        [Description("Reference FASTA text.")] string reference,
        [Description("Query FASTA text.")] string query,
        [Description("Precomputed delta file text (dnadiff -d); null = run nucmer --maxmatch.")] string? delta = null,
        [Description("Reference file name printed in the delta header (default reference.fa).")] string referenceName = "reference.fa",
        [Description("Query file name printed in the delta header (default query.fa).")] string queryName = "query.fa",
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(query));
        var r = GenomicAnalyzer.Dnadiff(reference, query, referenceName, queryName, delta, cancellationToken);
        return new DnadiffToolResult(r.Report, r.Delta, r.Delta1, r.DeltaM, r.Coords1, r.CoordsM, r.Snps, r.ReferenceDiff, r.QueryDiff, r.UnalignedReference, r.UnalignedQuery);
    }
}
