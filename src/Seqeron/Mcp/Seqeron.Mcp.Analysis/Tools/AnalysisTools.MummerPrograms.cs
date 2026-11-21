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
}
