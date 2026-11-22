using System.Text.Json.Serialization;

namespace Seqeron.Mcp.Analysis.Tools;

// B09 A3-8 / A3-12 / A3-13: results of the EMBOSS program ports (etandem / equicktandem, fuzzpro / fuzztran, seqmatchall, palindrome / einverted).

/// <summary>A repeat of <c>find_tandem_repeats_emboss</c>: 0-based start, length, the program's score, unit size, copies, and (etandem only) the percent identity and the consensus unit.</summary>
public record EmbossTandemItem(
    int Position,
    int Length,
    int Score,
    int RepeatSize,
    int RepeatCount,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] double? Identity = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] string? Consensus = null);

/// <summary>Result of <c>find_tandem_repeats_emboss</c>: the repeats in the program's report order.</summary>
public record FindEmbossTandemRepeatsResult(string Program, EmbossTandemItem[] Repeats, [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] string? Report = null, [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] string? OriginalFile = null);

/// <summary>A match of <c>seqmatchall</c>: sequence indices (0-based, first &lt; second), 0-based starts and the length.</summary>
public record WordMatchItem(int Sequence1, int Sequence2, int Position1, int Position2, int Length);

/// <summary>Result of <c>seqmatchall</c>: every maximal exact word match of every pair.</summary>
public record SeqMatchAllResult(
    WordMatchItem[] Matches,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] ReportedWordMatchItem[]? ReportedMatches = null);

/// <summary>A seqmatchall match as the program's report prints it with <c>sequenceQualifiers</c>: 1-based positions in the original sequences, start &gt; end on a reversed one (A5-7).</summary>
public record ReportedWordMatchItem(int Sequence1, int Sequence2, int Start1, int End1, int Start2, int End2, int Length, bool Reversed);

/// <summary>A hit of <c>fuzzpro</c>: 0-based start and length (<c>mismatches</c>: substituted residues of a search with <c>mismatch</c> &gt; 0, omitted for an exact search).</summary>
public record ProteinPatternHitItem(int Position, int Length, [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)] int Mismatches = 0);

/// <summary>Result of <c>fuzzpro</c>: the normalised pattern and every matching window.</summary>
public record FuzzproResult(string Pattern, ProteinPatternHitItem[] Hits, [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] string? Report = null);

/// <summary>A hit of <c>fuzztran</c>: frame, protein start / end (1-based), the matched residues, the forward-strand nucleotide location (0-based) and the coordinates the program prints.</summary>
public record TranslatedPatternHitItem(
    int Frame, int ProteinStart, int ProteinEnd, string Translation, int Position, int Length,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] int? EmbossStart = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] int? EmbossEnd = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)] int Mismatches = 0);

/// <summary>Result of <c>fuzztran</c>: the normalised pattern and the hits per frame.</summary>
public record FuzztranResult(string Pattern, TranslatedPatternHitItem[] Hits, [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] string? Report = null);

/// <summary>Result of <c>showalign</c>: the text the program writes (blocks of number line, ruler, reference, rows and consensus).</summary>
public record ShowAlignResult(string Text);

/// <summary>Result of <c>megamerger</c>: the merged sequence, the program's report and the minimal set of anchoring matches (0-based starts; empty merge = nothing to merge).</summary>
public record MegaMergerResult(string MergedSequence, string Report, WordMatchItem[] Matches);

/// <summary>One region of difference of <c>diffseq</c>: 1-based coordinates in both sequences (a gap has length 0: its start is after the gap, its end before it).</summary>
public record DiffseqDifferenceItem(int Start1, int End1, int Length1, int Start2, int End2, int Length2);

/// <summary>Result of <c>seqret</c>: the text of the sequence output file, the text of the feature output file (empty when no feature was written), the number of sequences written, the messages the library printed and, with outputSingle, the file of every sequence by path (directory/name.extension).</summary>
public record SeqretToolResult(string SequenceText, string FeatureText, int Sequences, string[] Messages, IReadOnlyDictionary<string, string>? SingleFiles = null);

/// <summary>Result of <c>union</c>: the text of the sequence output file, the overlap of every sequence (empty without overlapFile) and the messages the library printed.</summary>
public record UnionToolResult(string SequenceText, string OverlapText, string[] Messages);

/// <summary>Result of <c>seqretsplit</c>: the text of the sequence output file (empty when every sequence has a file of its own), the file of every sequence by path, the number of sequences written and the messages the library printed.</summary>
public record SeqretsplitToolResult(string SequenceText, IReadOnlyDictionary<string, string> SingleFiles, int Sequences, string[] Messages);

/// <summary>Result of <c>maskfeat</c>: the text of the sequence output file (FASTA), the number of sequences written and the messages the library printed.</summary>
public record MaskfeatToolResult(string SequenceText, int Sequences, string[] Messages);

/// <summary>Result of <c>descseq</c>: the text of the sequence output file (FASTA), the number of sequences written and the messages the library printed.</summary>
public record DescseqToolResult(string SequenceText, int Sequences, string[] Messages);

/// <summary>Result of <c>feattext</c>: the text of every entry the feature readers read (only the EMBL / GenBank and the GFF3 readers keep it), the number of tables and the messages the library printed.</summary>
public record FeattextToolResult(string Text, int Tables, string[] Messages);

/// <summary>Result of <c>diffseq</c>: the differences, the overlap of the two sequences (-1 when no match), the substitution counts of a nucleotide pair, the report text (default format <c>diffseq</c>) and the two feature files (GFF3 unless <c>-offormat4</c> / <c>-offormat5</c> name another format; <c>FeatureFileName1</c> / <c>FeatureFileName2</c> = the file each is written to).</summary>
public record WordmatchHitItem(string Pattern, int PatternStart, string Query, int QueryStart, int Length);

/// <summary>Result of <c>wordmatch</c>: the matches (0-based starts; every match is extended to the right beyond its word), the statistics and the program's files.</summary>
public record WordmatchToolResult(WordmatchHitItem[] Matches, int PatternCount, long AllMatches, long SumLength, string Alignments, string Log, string FeatureFile1, string FeatureFile2, string? FeatureFileName1 = null, string? FeatureFileName2 = null);

/// <summary>Result of <c>featcopy</c>: the output file name, its text (every feature table of the input in the output format), the format, the number of tables written and the messages the library printed.</summary>
public record FeatcopyToolResult(string FileName, string Text, string Format, int Tables, string[] Messages);

/// <summary>Result of <c>diffseq</c>: see <see cref="DiffseqToolResult"/>.</summary>
public record DiffseqToolResult(DiffseqDifferenceItem[] Differences, int OverlapStart1, int OverlapEnd1, int OverlapStart2, int OverlapEnd2, int Snps, int Transitions, int Transversions, string Report, string FeatureFile1, string FeatureFile2, string? FeatureFileName1 = null, string? FeatureFileName2 = null);

/// <summary>A match of a <c>nucmer_clusters</c> cluster: 1-based starts (the second one in the reverse-complemented query on the reverse strand), the length left after the overlap cut and the overlap.</summary>
public record NucmerClusterMatchItem(long Start1, long Start2, long Length, long Overlap);

/// <summary>A cluster of <c>nucmer_clusters</c>: its strand and matches.</summary>
public record NucmerClusterItem(string Strand, NucmerClusterMatchItem[] Matches);

/// <summary>Result of <c>nucmer_clusters</c>: the clusters of the forward strand, then of the reverse-complement strand.</summary>
public record NucmerClustersResult(NucmerClusterItem[] Clusters);

/// <summary>The text a Mash command writes (<c>mash_bounds</c>, <c>mash_dist</c>, <c>mash_triangle</c>), byte for byte as the program prints it.</summary>
public record MashTextResult(string Text);

/// <summary>A Mash sketch file (<c>mash_sketch</c>, <c>mash_paste</c>): <c>msh:</c> + the base64 of the <c>.msh</c> bytes (accepted as an input by the other Mash tools), their number and the number of sketches in it.</summary>
public record MashSketchFileResult(string Msh, int Bytes, int Sketches);

/// <summary>Result of <c>emboss_report_file</c>: the <c>-outfile</c> report text of the program and, for etandem / equicktandem, the <c>-origfile</c> text.</summary>
public record EmbossReportFileResult(string Report, [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] string? OriginalFile = null);
