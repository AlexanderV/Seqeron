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
}
