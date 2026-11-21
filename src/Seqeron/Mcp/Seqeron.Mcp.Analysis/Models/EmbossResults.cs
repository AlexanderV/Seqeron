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

/// <summary>A hit of <c>fuzzpro</c>: 0-based start and length.</summary>
public record ProteinPatternHitItem(int Position, int Length);

/// <summary>Result of <c>fuzzpro</c>: the normalised pattern and every matching window.</summary>
public record FuzzproResult(string Pattern, ProteinPatternHitItem[] Hits, [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] string? Report = null);

/// <summary>A hit of <c>fuzztran</c>: frame, protein start / end (1-based), the matched residues, the forward-strand nucleotide location (0-based) and the coordinates the program prints.</summary>
public record TranslatedPatternHitItem(
    int Frame, int ProteinStart, int ProteinEnd, string Translation, int Position, int Length,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] int? EmbossStart = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] int? EmbossEnd = null);

/// <summary>Result of <c>fuzztran</c>: the normalised pattern and the hits per frame.</summary>
public record FuzztranResult(string Pattern, TranslatedPatternHitItem[] Hits, [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] string? Report = null);

/// <summary>Result of <c>showalign</c>: the text the program writes (blocks of number line, ruler, reference, rows and consensus).</summary>
public record ShowAlignResult(string Text);

/// <summary>Result of <c>megamerger</c>: the merged sequence, the program's report and the minimal set of anchoring matches (0-based starts; empty merge = nothing to merge).</summary>
public record MegaMergerResult(string MergedSequence, string Report, WordMatchItem[] Matches);

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
