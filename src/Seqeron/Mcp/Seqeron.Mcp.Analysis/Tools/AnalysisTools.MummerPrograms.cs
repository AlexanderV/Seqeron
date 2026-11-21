using System.ComponentModel;
using ModelContextProtocol.Server;
using global::Seqeron.Genomics.Alignment;
using Seqeron.Mcp.Analysis.Models;

namespace Seqeron.Mcp.Analysis.Tools;

// The other programs of the MUMmer 3.23 package on the Analysis server (B09 A3-9d / A3-9e); every tool delegates to the ports in Seqeron.Genomics.Alignment.
public partial class AnalysisTools
{
    [McpServerTool(Name = "mummer", Title = "Genome — MUMmer mummer (maximal matches)", ReadOnly = true)]
    [Description("MUMmer 3.23 mummer (Kurtz' maxmat3): the maximal exact matches between a reference file and one or more query files. matchMode 'mumreference' (default, -mumreference / -mumcand) reports matches unique in the reference, 'mum' (-mum) matches unique in both, 'maxmatch' (-maxmatch) every maximal match; minLength = -l (default 20); nucleotidesOnly (-n) matches only a, c, g, t (any case) - without it every character is a letter of the alphabet (also protein text), both (-b) adds the reverse-complement strand of the query, reverseOnly (-r) only that strand. Output per query record: a '> id' line ('> id Reverse' for the reverse strand, '  Len = n' appended with showSequenceLengths = -L) and one line per match: reference position, query position, length (%8lu columns; with several reference records or fourColumn = -F the reference id is the first column; for reverse matches the query position is a position in the reverse-complemented query, or with showReversePositions = -c relative to the original query), followed by the matching text with showString (-s). Matches come in the order the suffix tree of the program enumerates them (per query position, the reference occurrences in lexicographic order of their suffixes; -mum by reference position). The reference records are joined by a separator no match can span; text before the first '>' belongs to the first record; the residues are lower-cased. The file read errors of the program (empty file, no header, only empty records) and its option conflicts (-b with -r, -c without -b / -r, a non-positive -l, no query file, more than 32 query files) are errors. Byte for byte equal to the installed mummer on 3 200 random runs (every mode, every option, multi-record files, N / IUPAC / odd characters, CR LF).")]
    public static MummerProgramResult Mummer(
        [Description("Reference FASTA text (multi-record).")] string reference,
        [Description("Query FASTA texts (at most 32 files); each is matched in turn.")] string[] queries,
        [Description("Matching algorithm: 'mumreference' (default, -mumreference / -mumcand), 'mum' (-mum) or 'maxmatch' (-maxmatch).")] string matchMode = "mumreference",
        [Description("-n: match only the characters a, c, g, t (default false).")] bool nucleotidesOnly = false,
        [Description("-l: minimum length of a match (default 20).")] long minLength = 20,
        [Description("-b: compute forward and reverse-complement matches (default false).")] bool both = false,
        [Description("-r: compute only reverse-complement matches (default false).")] bool reverseOnly = false,
        [Description("-s: show the matching substrings (default false).")] bool showString = false,
        [Description("-c: report the query position of a reverse match relative to the original query (needs both or reverseOnly; default false).")] bool showReversePositions = false,
        [Description("-F: force the 4 column format also for a single reference record (default false).")] bool fourColumn = false,
        [Description("-L: show the length of the query sequence on the header line (default false).")] bool showSequenceLengths = false,
        [Description("The reference file name shown in the messages (default reference.fa).")] string referenceName = "reference.fa",
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(reference));
        if (queries is null || queries.Length == 0)
            throw new ArgumentException("At least one query file is required", nameof(queries));
        var mode = (matchMode ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "mumreference" or "mumcand" => MummerMatchMode.MumReference,
            "mum" => MummerMatchMode.Mum,
            "maxmatch" => MummerMatchMode.MaxMatch,
            _ => throw new ArgumentException("matchMode must be 'mumreference', 'mum' or 'maxmatch'", nameof(matchMode)),
        };
        var r = MummerProgram.Run(reference, queries, new MummerOptions
        {
            Mode = mode, MatchNucleotidesOnly = nucleotidesOnly, MinLength = minLength, Both = both, OnlyReverse = reverseOnly, ShowString = showString, ShowReversePositions = showReversePositions,
            FourColumn = fourColumn, ShowSequenceLengths = showSequenceLengths, ReferenceName = referenceName, QueryNames = queries.Select((_, i) => queries.Length == 1 ? "query.fa" : $"query{i + 1}.fa").ToList(),
        }, cancellationToken);
        return new MummerProgramResult(r.Text, r.Messages);
    }

    [McpServerTool(Name = "promer", Title = "Genome — MUMmer promer (amino-acid alignments)", ReadOnly = true)]
    [Description("MUMmer 3.23 promer (PROtein MUMmer 3.07): amino-acid alignments between two multi-FASTA DNA files, the whole pipeline: prepro (both files translated in all six reading frames; a stretch of at most maskLength amino acids between stop codons is masked) | mummer (anchors in amino acids: matchMode 'mumreference' (default), 'mum', 'maxmatch'; minMatch = -l, default 6) | mgaps (clusters: maxGap = -g (30), minCluster = -c (20), diagonalFactor = -d (0.11)) | postpro (the affine Smith-Waterman extension of the clusters in amino-acid space with the BLOSUM matrix of matrix (-x: 1 = BLOSUM 45, 2 = BLOSUM 62 default, 3 = BLOSUM 80; 0 and other values are refused with a warning), give-up distance breakLength = -b (60 amino acids); extend / optimize as --[no]extend / --[no]optimize). Returns the .delta file: first line = the two file names, then PROMER, per record pair '>refId qryId refLen qryLen' and per alignment 'sR eR sQ eQ errors similarityErrors stopCodons' in nucleotide coordinates of the original DNA (a reverse-strand alignment has sR > eR / sQ > eQ) followed by the delta values - amino-acid distances, one delta = three nucleotides - and 0. With delta = false (--nodelta) the .cluster text is returned instead; coords (-o) adds the PROmer 1.1 coords file (show-coords -r of the delta). The prefix (-p) names the files of the result and the temporary files in the messages. Messages are the progress lines of the script, the warnings of the programs ('Forcing unrecognized DNA char to N', invalid break length / matrix) and mummer's reading lines. Byte for byte equal to the installed promer on 320 recorded and 1 920 further random runs (delta, cluster, coords, the temporary files and the messages; every option). The inputs on which the pipeline aborts (a file without sequence, a non-positive minMatch) are refused with an error.")]
    public static PromerToolResult Promer(
        [Description("Reference FASTA text (multi-record DNA).")] string reference,
        [Description("Query FASTA text (multi-record DNA).")] string query,
        [Description("Anchor matches: 'mumreference' (default, --mumreference / --mumcand), 'mum' (--mum) or 'maxmatch' (--maxmatch).")] string matchMode = "mumreference",
        [Description("Minimum length of a single match in amino acids (-l, default 6).")] int minMatch = 6,
        [Description("Maximum gap between two adjacent matches of a cluster in amino acids (-g, default 30).")] long maxGap = 30,
        [Description("Minimum length of a cluster of matches in amino acids (-c, default 20).")] long minCluster = 20,
        [Description("Diagonal difference as a fraction of the gap length (-d, default 0.11).")] double diagonalFactor = 0.11,
        [Description("Distance an extension tries to extend poor scoring regions before giving up, in amino acids (-b, default 60; outside 1 - 10000 ignored with a warning).")] int breakLength = 60,
        [Description("Alignment matrix (-x): 1 = BLOSUM 45, 2 = BLOSUM 62 (default), 3 = BLOSUM 80.")] int matrix = 2,
        [Description("Maximum bookend masking length in amino acids (-m, default 8; negative = invalid, ignored with a warning).")] int maskLength = 8,
        [Description("Extend the clusters (--extend, default true).")] bool extend = true,
        [Description("Create the delta (--delta, default true); false returns the cluster file (--nodelta).")] bool delta = true,
        [Description("Alignment score optimisation (--optimize, default true); false = --nooptimize.")] bool optimize = true,
        [Description("Also return the PROmer 1.1 coords file (-o / --coords: show-coords -r of the delta; default false; not with delta = false).")] bool coords = false,
        [Description("Output prefix (-p, default out): names the files in the messages.")] string prefix = "out",
        [Description("Reference file name printed in the first line of the delta (promer prints the absolute path).")] string referenceName = "reference.fa",
        [Description("Query file name printed in the first line of the delta.")] string queryName = "query.fa",
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(query));
        var mode = (matchMode ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "mumreference" or "mumcand" => MummerMatchMode.MumReference,
            "mum" => MummerMatchMode.Mum,
            "maxmatch" => MummerMatchMode.MaxMatch,
            _ => throw new ArgumentException("matchMode must be 'mumreference', 'mum' or 'maxmatch'", nameof(matchMode)),
        };
        var r = PromerProgram.Run(reference, query, referenceName, queryName, new PromerOptions
        {
            MatchMode = mode, MinMatch = minMatch, MaxGap = maxGap, MinCluster = minCluster, DiagonalFactor = diagonalFactor, BreakLength = breakLength, Matrix = matrix, MaskLength = maskLength,
            Extend = extend, Delta = delta, Optimize = optimize, Coords = coords, Prefix = prefix,
        }, cancellationToken);
        return new PromerToolResult(r.Delta, r.Cluster, r.Messages, r.Coords, r.AminoReference, r.AminoQuery, r.Mgaps);
    }

    [McpServerTool(Name = "mummer_gaps", Title = "Genome — MUMmer gaps (consistent matches)", ReadOnly = true)]
    [Description("MUMmer 3.23 gaps: reads a list of matches between two strings (start1 start2 length, as mummer prints them without the header lines; text that is not a number is refused - the program loops forever on it - and an empty list too - the program crashes) and prints first the longest consistent set of matches (ordered along both strings; overlaps are cut off: the overlap column shows the cut) after '> fileName [reverse] Consistent matches', then all the other matches after '> fileName [reverse] Other matches', in the format of mgaps: start1 start2 length overlap gap1 gap2 (the first line of a block has none - -; the other matches show a gap only after a match they chain to). Ties in the second start keep their input order. reverseComplement (-r) writes 'reverse' into the header lines, which mummer-annotate reads. Byte for byte equal to the installed gaps on 259 random match lists (+ the empty-list crash).")]
    public static MummerTextResult MummerGaps(
        [Description("The matches: white-space separated numbers, three per match.")] string matches,
        [Description("The first argument of the program, printed in the header lines (the file name of the first sequence).")] string fileName = "reference.fa",
        [Description("The second argument -r: the header lines read 'reverse' (default false).")] bool reverseComplement = false)
    {
        if (string.IsNullOrWhiteSpace(matches))
            throw new ArgumentException("Matches cannot be null or empty", nameof(matches));
        return new MummerTextResult(global::Seqeron.Genomics.Alignment.MummerGaps.Run(fileName, reverseComplement, matches));
    }

    [McpServerTool(Name = "delta2blocks", Title = "Genome — MUMmer delta2blocks (CIGAR blocks)", ReadOnly = true)]
    [Description("delta2blocks of the Debian MUMmer 3.23 package: the alignments of a nucmer delta file between one reference record and one query record as one line each, 'refId qryId sR eR sQ eQ CIGAR', the CIGAR written with M (matches), X (mismatches), I (an extra base in the reference, a positive delta) and D (an extra base in the query, a negative delta); the coordinates are those of the delta file (a reverse-strand alignment has sR > eR / sQ > eQ and is compared on the reverse complement). sortByReference (-r) / sortByQuery (-q) order the alignments. The sequences are read from the FASTA texts named in the delta header. A promer delta file is refused (the program reads the sequences through a null pointer and crashes). Byte for byte equal to the installed delta2blocks on 226 + runs.")]
    public static MummerProgramResult Delta2Blocks(
        [Description("Delta file text (nucmer).")] string delta,
        [Description("Reference FASTA text.")] string referenceFasta,
        [Description("Query FASTA text.")] string queryFasta,
        [Description("Reference record id.")] string referenceId,
        [Description("Query record id.")] string queryId,
        [Description("-q: sort the alignments by the query start (default false).")] bool sortByQuery = false,
        [Description("-r: sort the alignments by the reference start (default false).")] bool sortByReference = false)
    {
        if (string.IsNullOrEmpty(delta))
            throw new ArgumentException("Delta file text cannot be null or empty", nameof(delta));
        try
        {
            var r = MummerDelta2Blocks.Run(delta, referenceFasta, queryFasta, referenceId, queryId, new MummerShowAlignsOptions { SortByQuery = sortByQuery, SortByReference = sortByReference });
            return new MummerProgramResult(r.Text, r.Messages);
        }
        catch (InvalidOperationException e)
        {
            throw new NotSupportedException(e.Message, e);
        }
    }

    [McpServerTool(Name = "mummer_annotate", Title = "Genome — MUMmer mummer-annotate (gap alignments)", ReadOnly = true)]
    [Description("MUMmer 3.23 mummer-annotate (annotate): reads the output of gaps (the first sequence is named on its '>' lines - supply it in files by that name - and the word 'reverse' there reverse-complements the query), repeats the gap file and, under every gap line, prints the edit-distance alignment of the stretch between the two matches (10 bases of each flanking match included): '    Errors = n' and the rows 'T:' (first sequence) / 'S:' (query) / mismatch marks, 60 columns per block; a gap longer than 10 000 prints '*** Too long ***'. withErrors is the file witherrors.gaps the program writes: the gap lines with the number of errors appended. A window that reaches before the first base of a sequence reads memory in front of the program's buffer (valgrind: invalid read) and is refused with an error; beyond the end the buffer holds zeros. Byte for byte equal to the installed program on 230 + runs.")]
    public static MummerAnnotateToolResult MummerAnnotate(
        [Description("The gap file text (the output of mummer_gaps).")] string gaps,
        [Description("The data file (the query): its first record is used.")] string dataFasta,
        [Description("The FASTA texts of the files named on the '>' lines of the gap file, by name.")] Dictionary<string, string> files)
    {
        if (string.IsNullOrEmpty(gaps))
            throw new ArgumentException("Gap file text cannot be null or empty", nameof(gaps));
        if (string.IsNullOrEmpty(dataFasta))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dataFasta));
        try
        {
            var r = global::Seqeron.Genomics.Alignment.MummerAnnotate.Run(gaps, dataFasta, files ?? []);
            return new MummerAnnotateToolResult(r.Text, r.WithErrors);
        }
        catch (InvalidOperationException e)
        {
            throw new NotSupportedException(e.Message, e);
        }
    }

    [McpServerTool(Name = "run_mummer1", Title = "Genome — MUMmer run-mummer1 (mummer 1 suite)", ReadOnly = true)]
    [Description("MUMmer 3.23 run-mummer1, the 'severely antiquated' script that runs the mummer 1 suite: mummer -mum -l 20 -r (reverse-complement matches only; the first output line is dropped) | gaps | mummer-annotate. reference and query are single-sequence FASTA texts; referenceName is the file name printed in the gaps file. Returns the files of the script: out (matches), gaps, align (the gaps with their alignments) and errorsGaps (the gaps with the number of errors). Without reverse (the script's fourth argument -r) the script passes an empty argument to mummer, which fails, and the tool refuses; a query with a second record makes gaps loop forever (refused), no match makes gaps crash (refused). Byte for byte equal to the installed script on 55 + runs.")]
    public static RunMummer1ToolResult RunMummer1(
        [Description("Reference FASTA text (one sequence).")] string reference,
        [Description("Query FASTA text (one sequence).")] string query,
        [Description("The reference file name as the script receives it (printed in the gaps file).")] string referenceName = "reference.fa",
        [Description("The prefix of the output files (names the files in the result).")] string prefix = "out",
        [Description("The fourth argument -r (default true; without it the script fails).")] bool reverse = true)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(query));
        try
        {
            var r = MummerRunPrograms.RunMummer1(reference, query, referenceName, prefix, reverse);
            return new RunMummer1ToolResult(r.Out, r.Gaps, r.Align, r.ErrorsGaps);
        }
        catch (InvalidOperationException e)
        {
            throw new NotSupportedException(e.Message, e);
        }
    }

    [McpServerTool(Name = "delta2maf", Title = "Genome — MUMmer delta2maf (MAF alignments)", ReadOnly = true)]
    [Description("delta2maf of the Debian MUMmer 3.23 package: the alignments of a nucmer delta file as a MAF file: '##maf version=1 scoring=single_cov2', per alignment 'a score=n' (n = |end - start| in the reference) and two 's' lines (source id, 0-based start, aligned length, strand, source length, the gapped sequence text in lower case; a reverse-strand alignment is written on the reverse complement with strand '-'), and '##eof maf' (no final newline). referenceId / queryId select the alignments of one reference record, one query record or one pair (a query id alone needs an empty referenceId); all alignments when both are empty. sortByReference (-r) / sortByQuery (-q) order the alignments. The reference and query files named in the delta header are read into one map by id (the first record of an id wins). An id without alignments is an error (the program prints the maf header and exits with 1); a promer delta file is refused (the program reads its amino-acid deltas as nucleotide deltas and runs past the sequences). Byte for byte equal to the installed delta2maf on 315 + runs.")]
    public static MummerProgramResult Delta2Maf(
        [Description("Delta file text (nucmer).")] string delta,
        [Description("Reference FASTA text.")] string referenceFasta,
        [Description("Query FASTA text.")] string queryFasta,
        [Description("Only the alignments of this reference record (default: any).")] string referenceId = "",
        [Description("Only the alignments of this query record (default: any; both ids must match when both are given).")] string queryId = "",
        [Description("-q: sort the alignments by the query start (default false).")] bool sortByQuery = false,
        [Description("-r: sort the alignments by the reference start (default false).")] bool sortByReference = false)
    {
        if (string.IsNullOrEmpty(delta))
            throw new ArgumentException("Delta file text cannot be null or empty", nameof(delta));
        var r = MummerDelta2Maf.Run(delta, referenceFasta, queryFasta, referenceId ?? string.Empty, queryId ?? string.Empty, new MummerDelta2MafOptions { SortByQuery = sortByQuery, SortByReference = sortByReference });
        return new MummerProgramResult(r.Text, r.Messages);
    }

    [McpServerTool(Name = "nucmer2xfig", Title = "Genome — MUMmer nucmer2xfig (xfig plot)", ReadOnly = true)]
    [Description("nucmer2xfig of MUMmer 3.23 (perl script, Steven Salzberg 2001): an xfig file (FIG 3.2) plotting a coords file of 'show-coords -c -l' (19 white-space separated fields per alignment: s1 e1 | s2 e2 | len1 len2 | %idy | lenR lenQ | covR covQ | idR idQ; no header lines, sorted by query id and start) of the contigs of one genome against a reference chromosome: the reference as a red line (colour 4) with its id, x-axis tick marks every 250 000 bases and y-axis tick marks for 50 - 100 percent identity, every alignment as a green line (colour 2, x scaled by 0.005, y = 20 x (identity - 50)), every contig - the stretch its alignments span, split where they are more than twice the query length apart - as a blue line (colour 1) at 100 percent with a vertical label (4.7124 rad) when the line is wider than the label. Perl's numeric conversion of the fields is reproduced (text is 0) and printf %.0f rounds half to even. Byte for byte equal to the installed script on 260 + runs including malformed input (the output of an empty file is just the reference line with empty fields).")]
    public static MummerTextResult Nucmer2Xfig(
        [Description("The text of the coords file (show-coords -c -l without its header).")] string coords)
    {
        if (coords is null)
            throw new ArgumentException("Coords text cannot be null", nameof(coords));
        return new MummerTextResult(MummerNucmer2Xfig.Run(coords));
    }

    [McpServerTool(Name = "combine_mums", Title = "Genome — MUMmer combineMUMs (extend and combine clusters)", ReadOnly = true)]
    [Description("MUMmer 3.23 combineMUMs: takes the clusters of exact matches of mgaps (gaps = its output: '> id' / '> id Reverse' lines and 'start1 start2 length ...' lines; a '#' line ends a cluster) and extends every cluster off both ends and fills the gaps between its matches by edit-distance alignments (Myers / Landau-Vishkin diagonals; branch point score 0.272 per match and -0.728 per error in single precision, or 1 - (1 - errorRate) with errorRate = -e; give up 200 bases past the best score, at most 500 errors, extensions in steps of 10 000 bases). reference = a FASTA text (only its first record is read), query = the multi-FASTA texts of the query records; the n-th forward block of the gaps file belongs to the n-th query record and the reverse block after it to its reverse complement. Output (text) = the gaps file with, under every gap line, '     Errors = n' and the alignment (rows A: reference, B: query and a mark row with the 10 matching bases of each flanking match, 60 columns per block) or '   No extension', and a 'Region:  refStart .. refEnd  qryStart .. qryEnd  errors / length  percent' line after every combined region (a region on the reverse strand shows its query coordinates relative to the original query). onlyDifferences (-D) writes only the positions and characters of the differences instead; showDifferences (-S) every difference of the extensions and gaps; verbose (-v 1) every extension with its alignment; nucleotidesOnly (-n) lets only a, c, g, t match (others become '2' in the reference and '1' in the query); consecutiveNonAcgt (-N) breaks matches at that many consecutive non-ACGT characters; tagFromFasta (-t) tags the matches with the first word of the FASTA headers; queryTag (-q, default 'Query') / referenceTag (-r, default 'Ref') are suffixes of the cover file names. errorsFile = the file -W names (default witherrors.gaps): the gap lines with the number of errors appended. Without noCoverFiles (-x) the files '<reference base name><referenceTag>.cover' and '<query base name><queryTag>.cover' (regions and gaps with their N counts) are returned too; messages = standard error (a summary of the matches and the coverage of both sequences). A failed assertion of the program (no reference record, fewer query records than forward blocks, a line of the gaps file that is not three numbers) and a coordinate beyond a sequence (the program reads outside its buffers) are errors; the coverage report of a reference whose length is 7 modulo 16 reads the byte after its realloc'd buffer (valgrind: invalid read; the value printed depends on the allocator) and is refused unless noCoverFiles is set. Byte for byte equal to the installed combineMUMs on 240 + runs (every option, both strands, N / IUPAC characters, extensions of 24 kb) and 3 000 more random runs.")]
    public static CombineMumsToolResult CombineMums(
        [Description("Reference FASTA text (only the first record is read).")] string reference,
        [Description("Multi-FASTA text of the query records.")] string query,
        [Description("The gaps file (output of mgaps).")] string gaps,
        [Description("-D: write only the positions and characters of the differences (default false).")] bool onlyDifferences = false,
        [Description("-e: the error rate cutoff of the extensions, e.g. 0.02 (default: the built-in branch point score).")] double? errorRate = null,
        [Description("-n: only a, c, g, t match (default false).")] bool nucleotidesOnly = false,
        [Description("-N: break matches at this many consecutive non-ACGT characters (default 0 = never).")] int consecutiveNonAcgt = 0,
        [Description("-q: suffix of the query tag (default Query).")] string queryTag = "Query",
        [Description("-r: suffix of the reference tag (default Ref).")] string referenceTag = "Ref",
        [Description("-S: show every difference of the extensions and gaps (default false).")] bool showDifferences = false,
        [Description("-t: tag the matches with the first word of the FASTA headers (default false).")] bool tagFromFasta = false,
        [Description("-v: verbose level, above 0 writes every extension with its alignment (default 0).")] int verbose = 0,
        [Description("-x: do not write the cover files (default false).")] bool noCoverFiles = false,
        [Description("-W: name of the errors file (default witherrors.gaps).")] string errorsFileName = "witherrors.gaps",
        [Description("The reference file name as given on the command line (its base name is the reference tag; default reference.fa).")] string referenceName = "reference.fa",
        [Description("The query file name as given on the command line (its base name is the query tag; default query.fa).")] string queryName = "query.fa")
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(query));
        if (gaps is null)
            throw new ArgumentException("Gaps file text cannot be null", nameof(gaps));
        var r = MummerCombineMums.Run(reference, query, gaps, new MummerCombineMumsOptions
        {
            OnlyDifferencePositions = onlyDifferences, ErrorRate = errorRate, NucleotidesOnly = nucleotidesOnly, ConsecutiveNonAcgt = consecutiveNonAcgt, QuerySuffix = queryTag ?? "Query", ReferenceSuffix = referenceTag ?? "Ref",
            ShowDifferences = showDifferences, TagFromFastaLine = tagFromFasta, Verbose = verbose, NoCoverFiles = noCoverFiles, ErrorsFileName = errorsFileName ?? "witherrors.gaps", ReferenceName = referenceName ?? "reference.fa", QueryName = queryName ?? "query.fa",
        });
        return new CombineMumsToolResult(r.Text, r.ErrorsFile, r.Log, r.Files.ToDictionary(f => f.Key, f => f.Value));
    }

    [McpServerTool(Name = "run_mummer3", Title = "Genome — MUMmer run-mummer3 (mummer 3 suite)", ReadOnly = true)]
    [Description("MUMmer 3.23 run-mummer3, the script that runs the basic mummer 3 suite: mummer -mumreference -b -l 20 (maximal matches unique in the reference, both strands) | mgaps -l 100 -f .12 -s 600 (clusters) | combineMUMs -x -e .10 -W prefix.errorsgaps (extension of the clusters; no cover files). reference is a FASTA text (only its first record is aligned), query the multi-FASTA text; referenceName / queryName are the file names as the script receives them (printed in the messages of mummer). Returns the files of the script - out (the matches; the reverse matches are relative to the reverse complement of the query), gaps (the clusters), align (the alignments of the gaps: see combine_mums) and errorsGaps - the script's standard output (the three progress lines) and standard error (messages: the lines of mummer and the summary of combineMUMs). An empty input file or a query without a sequence stops the script (sh -e): refused. Byte for byte equal to the installed script on 216 + runs (files and messages) and 1 600 more random runs.")]
    public static RunMummer3ToolResult RunMummer3(
        [Description("Reference FASTA text (only the first record is aligned).")] string reference,
        [Description("Multi-FASTA text of the query records.")] string query,
        [Description("The reference file name as the script receives it.")] string referenceName = "reference.fa",
        [Description("The query file name as the script receives it.")] string queryName = "query.fa",
        [Description("The prefix of the output files (names the files in the result).")] string prefix = "out")
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(query));
        var r = MummerRunPrograms.RunMummer3(reference, query, referenceName, queryName, prefix);
        return new RunMummer3ToolResult(r.Out, r.Gaps, r.Align, r.ErrorsGaps, r.Stdout, r.Log);
    }

    [McpServerTool(Name = "repeat_match", Title = "Genome — MUMmer repeat-match (exact repeats)", ReadOnly = true)]
    [Description("MUMmer 3.23 repeat-match (Delcher): the maximal exact repeats of the first sequence of a FASTA file that are at least minLength long (-n, default 20), within the forward strand and, unless forwardOnly (-f), between the forward and the reverse-complement strand. Output: 'Long Exact Matches:' and a line per match 'start1 start2[r] length' (a line with the letter r after the second position matches the reverse-complement strand: start2 is read on the original sequence). tandemOnly (-t, implies forwardOnly) lists only the repeats whose two copies touch (start1 + length >= start2). The matches come in the order of the program's own suffix tree (a suffix tree of the sequence and its reverse complement in which a new child is put in front of its siblings), reproduced by the same construction. exhaustive (-E) uses the program's quadratic search instead (a different order); verbose (-V 1) prints the sizes of the program's structures, -V 2 the construction and the traversal of the tree. messages = standard error ('Genome Length = n   Used m internal nodes', the warnings of the FASTA reader). A genome that is an exact reverse-complement palindrome stops the program (exitCode 255, '*** Genome is exact palindrome ***' on standard output). An empty file, a file without a '>' line (the program works on an uninitialised buffer: valgrind use of uninitialised values; the output changes under MALLOC_PERTURB_) and a header whose first word has 500 or more characters (buffer overflow, abort) are refused with an error. Byte for byte equal to the installed repeat-match on 370 + runs and 3 300 more random runs (every option, both strands, N / IUPAC characters, genomes of 30 kb).")]
    public static RepeatMatchToolResult RepeatMatch(
        [Description("Genome FASTA text (only the first record is used).")] string genome,
        [Description("-n: minimum length of an exact match (default 20; at least 1).")] int minLength = 20,
        [Description("-f: forward strand only (default false).")] bool forwardOnly = false,
        [Description("-t: only tandem repeats; implies forwardOnly (default false).")] bool tandemOnly = false,
        [Description("-E: exhaustive (quadratic) search instead of the suffix tree (default false).")] bool exhaustive = false,
        [Description("-V: verbose level (default 0; 1 prints the structure sizes, 2 the tree construction and traversal).")] int verbose = 0)
    {
        if (string.IsNullOrEmpty(genome))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(genome));
        if (minLength < 1)
            throw new ArgumentException("ERROR:  Illegal min match length \"" + minLength.ToString(System.Globalization.CultureInfo.InvariantCulture) + "\"", nameof(minLength));
        var r = MummerRepeatMatch.Run(genome, new MummerRepeatMatchOptions { MinMatchLength = minLength, ForwardOnly = forwardOnly, TandemOnly = tandemOnly, Exhaustive = exhaustive, Verbose = verbose });
        return new RepeatMatchToolResult(r.Text, r.Messages, r.ExitCode);
    }

    [McpServerTool(Name = "exact_tandems", Title = "Genome — MUMmer exact-tandems (tandem repeats)", ReadOnly = true)]
    [Description("MUMmer 3.23 exact-tandems: the exact tandem repeats of the first sequence of a FASTA file with a duplicated unit of at least minLength bases. The script runs repeat-match -t -n N, drops its two header lines, sorts the matches by the first and second column (sort -k1n -k2n) and feeds them to tandem-repeat.awk, which writes the table 'Start Extent UnitLen Copies' with one line per repeat region (the first match of every start whose copies touch; extent = start2 + length - start1, unit = start2 - start1, copies = extent / unit with one decimal). By default (installedScript = false) the result is that intended pipeline. The script as installed by the Debian package is broken - its last command is split over two lines (awk -f $ and, on the next line, $scriptdir/tandem-repeat.awk) - so awk looks for a file named $ and fails with status 2 ('awk: fatal: cannot open source file `$' for reading: No such file or directory'); sh -e then stops the script without a table and without removing its temporary file: installedScript = true reproduces that (text = the two progress lines 'Finding matches' and 'Tandem repeats', exitCode 2, filesLeft = $$.tmp.matches). messages = standard error (the lines of repeat-match, then the awk error with installedScript). A minLength below 1 makes repeat-match fail (the table is just the header). An empty file or a file without a '>' line is refused (see repeat_match). Byte for byte equal to the installed script (installedScript = true) and to the pipeline run by hand on 215 + runs each and 1 200 more random runs each.")]
    public static ExactTandemsToolResult ExactTandems(
        [Description("Genome FASTA text (only the first record is used).")] string genome,
        [Description("Minimum length of the duplicated unit (repeat-match -n; default 20).")] int minLength = 20,
        [Description("Reproduce the script as installed (broken: awk -f $), instead of the intended pipeline (default false).")] bool installedScript = false)
    {
        if (string.IsNullOrEmpty(genome))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(genome));
        var r = installedScript ? MummerExactTandems.RunInstalled(genome, minLength) : MummerExactTandems.Run(genome, minLength);
        return new ExactTandemsToolResult(r.Text, r.Log, r.ExitCode, [.. r.LeftBehind]);
    }

    [McpServerTool(Name = "mapview", Title = "Genome — MUMmer mapview (FIG alignment maps)", ReadOnly = true)]
    [Description("MUMmer 3.23 mapview (perl script): turns the output of show-coords -rl (nucmer or promer; or of mgaps with mgapsReference = -g) into FIG 3.2 files that show the alignments along the reference - the percent identity (50 - 100) on the vertical axis, a red match line per alignment, a coloured link line between the alignments of one query contig (one of three colours; no link when the next match of the contig is more than maxDistance = -d (default 50 000) bases away) - partitioned into numberOfFiles = -n (default 10) files of an equal share of the reference (a multi-reference input gets one file per group of references that fit a share), named <prefix>_<n>.fig (prefix = -p; default 'NUCMER_graph' / 'PROMER_graph': the program name in the second line of the coords file), with a kilobase scale, the percent identity ticks and a legend. With a UTR file and a CDS file (GFF: sequence id, source, exon type single-exon / initial-exon / internal-exon / last-exon, start, end, score, strand, frame, gene name) the genes are drawn above the map (5' UTR, CDS, 3' UTR, introns; overlapping genes on levels), x1 / x2 (-x1 / -x2) bound the displayed range, printQueryNames (-I) / printGeneNames (-Ir) write the names. format (-f) 'pdf' or 'ps' makes the script convert every FIG file with fig2dev: the result then lists the tmp<prefix>_<n>.fig texts the script feeds to it (figs), the fig2dev command lines (commands) and the names of the converted files (outputs) - the conversion itself is done by fig2dev, which is not part of this library. verbose (-v) lists the files and their ranges on standard output (stdout). A UTR file without a CDS file, or GFF files whose first column names no reference of the coords file, stop the script (exitCode 1, the message on messages). An mgaps input (-g) also writes the file <coordsName>coords (side): the clusters in the format of show-coords, sorted by position. The messages of a wrong option (x1 > x2, a negative bound) are printed and the script goes on, as the program does. Byte for byte equal to the installed mapview on 263 + runs (nucmer and promer coords, one or several references, mgaps input, GFF genes on both strands, every option) and 800 more random runs.")]
    public static MapviewToolResult Mapview(
        [Description("The coords file (show-coords -rl output, or mgaps output with mgapsReference).")] string coords,
        [Description("UTR coordinates in GFF format (needs the CDS file too).")] string? utrGff = null,
        [Description("CDS coordinates in GFF format.")] string? cdsGff = null,
        [Description("-d: the largest distance in base pairs between linked matches (default 50 000).")] int maxDistance = 0,
        [Description("-f: output format 'fig' (default), 'pdf' or 'ps' (the last two through fig2dev: see the description).")] string format = "fig",
        [Description("-m: magnification handed to fig2dev (default none).")] double magnification = 0,
        [Description("-n: number of output files the reference is partitioned into (default 10).")] int numberOfFiles = 0,
        [Description("-p: prefix of the output files (default 'NUCMER_graph' / 'PROMER_graph' by the coords file).")] string prefix = "",
        [Description("-v: list the output files and their ranges (default false).")] bool verbose = false,
        [Description("-x1: lower coordinate bound of the display (default none).")] int x1 = 0,
        [Description("-x2: upper coordinate bound of the display (default none).")] int x2 = 0,
        [Description("-g: the input is mgaps output; the reference sequence id (as in the first column of the GFF files).")] string mgapsReference = "",
        [Description("-I: display the names of the query sequences (default false).")] bool printQueryNames = false,
        [Description("-Ir: display the names of the reference genes (default false).")] bool printGeneNames = false,
        [Description("The file name of the coords file as the script receives it (the mgaps side file is <name>coords; default coords).")] string coordsName = "coords",
        [Description("The file name of the UTR file; when it equals cdsName the exon lists are not split (default utr.gff).")] string utrName = "utr.gff",
        [Description("The file name of the CDS file (default cds.gff).")] string cdsName = "cds.gff")
    {
        if (string.IsNullOrEmpty(coords))
            throw new ArgumentException("Coords text cannot be null or empty", nameof(coords));
        if (format is not ("fig" or "pdf" or "ps"))
            throw new ArgumentException("format must be 'fig', 'pdf' or 'ps'", nameof(format));
        var r = MummerMapview.Run(coords, utrGff, cdsGff, new MummerMapviewOptions
        {
            MaxDistance = maxDistance, Format = format, Magnification = magnification, FileCount = numberOfFiles, Prefix = prefix ?? string.Empty, Verbose = verbose, X1 = x1, X2 = x2,
            MgapsReference = mgapsReference ?? string.Empty, PrintQueryNames = printQueryNames, PrintGeneNames = printGeneNames, CoordsName = coordsName ?? "coords", UtrName = utrName ?? "utr.gff", CdsName = cdsName ?? "cds.gff",
        });
        return new MapviewToolResult(r.Figs.ToDictionary(f => f.Key, f => f.Value), [.. r.Commands], [.. r.Outputs], r.Stdout, r.Stderr, r.ExitCode, r.Side.ToDictionary(f => f.Key, f => f.Value));
    }
}
