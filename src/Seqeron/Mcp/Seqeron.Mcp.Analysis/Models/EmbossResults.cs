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
