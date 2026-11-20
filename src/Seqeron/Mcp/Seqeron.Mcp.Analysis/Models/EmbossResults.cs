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
public record FindEmbossTandemRepeatsResult(string Program, EmbossTandemItem[] Repeats);

/// <summary>A match of <c>seqmatchall</c>: sequence indices (0-based, first &lt; second), 0-based starts and the length.</summary>
public record WordMatchItem(int Sequence1, int Sequence2, int Position1, int Position2, int Length);

/// <summary>Result of <c>seqmatchall</c>: every maximal exact word match of every pair.</summary>
public record SeqMatchAllResult(WordMatchItem[] Matches);

/// <summary>A hit of <c>fuzzpro</c>: 0-based start and length.</summary>
public record ProteinPatternHitItem(int Position, int Length);

/// <summary>Result of <c>fuzzpro</c>: the normalised pattern and every matching window.</summary>
public record FuzzproResult(string Pattern, ProteinPatternHitItem[] Hits);

/// <summary>A hit of <c>fuzztran</c>: frame, protein start / end (1-based), the matched residues, the forward-strand nucleotide location (0-based) and the coordinates the program prints.</summary>
public record TranslatedPatternHitItem(
    int Frame, int ProteinStart, int ProteinEnd, string Translation, int Position, int Length,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] int? EmbossStart = null,
    [property: JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)] int? EmbossEnd = null);

/// <summary>Result of <c>fuzztran</c>: the normalised pattern and the hits per frame.</summary>
public record FuzztranResult(string Pattern, TranslatedPatternHitItem[] Hits);

/// <summary>Result of <c>showalign</c>: the text the program writes (blocks of number line, ruler, reference, rows and consensus).</summary>
public record ShowAlignResult(string Text);

/// <summary>Result of <c>megamerger</c>: the merged sequence, the program's report and the minimal set of anchoring matches (0-based starts; empty merge = nothing to merge).</summary>
public record MegaMergerResult(string MergedSequence, string Report, WordMatchItem[] Matches);
