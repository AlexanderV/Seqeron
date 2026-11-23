namespace SuffixTree.Mcp.Core.Tools;

// ================================
// Suffix Tree Core Results
// ================================

/// <summary>
/// Result of suffix_tree_contains operation.
/// </summary>
/// <param name="Found">Whether the pattern was found in the text.</param>
public record SuffixTreeContainsResult(bool Found);

/// <summary>
/// Result of suffix_tree_count operation.
/// </summary>
/// <param name="Count">Number of pattern occurrences in the text.</param>
public record SuffixTreeCountResult(int Count);

/// <summary>
/// Result of suffix_tree_find_all operation.
/// </summary>
/// <param name="Positions">Array of positions where pattern was found.</param>
public record SuffixTreeFindAllResult(int[] Positions);

/// <summary>
/// Result of suffix_tree_lrs operation.
/// </summary>
public record SuffixTreeLrsResult(string Substring, int Length);

/// <summary>
/// Result of suffix_tree_lcs operation.
/// </summary>
public record SuffixTreeLcsResult(string Substring, int Length);

/// <summary>
/// Result of suffix_tree_stats operation.
/// </summary>
public record SuffixTreeStatsResult(int NodeCount, int LeafCount, int MaxDepth, int TextLength);

/// <summary>
/// One repeated substring with all its 0-based start positions (ascending).
/// </summary>
public record RepeatedSubstringItem(string Substring, int[] Positions);

/// <summary>
/// Result of suffix_tree_all_lrs operation.
/// </summary>
/// <param name="Substrings">Every longest repeated substring, ordered by first occurrence.</param>
/// <param name="Length">Common length of the substrings (0 when none).</param>
public record SuffixTreeAllLrsResult(RepeatedSubstringItem[] Substrings, int Length);

/// <summary>
/// One maximal match between the reference text and the query (0-based, forward strand).
/// </summary>
public record MaximalMatchItem(int PositionInText, int PositionInQuery, int Length);

/// <summary>
/// Result of suffix_tree_find_mems / suffix_tree_find_mums operations.
/// </summary>
/// <param name="Matches">Matches sorted by query position, then text position.</param>
public record SuffixTreeMaximalMatchesResult(MaximalMatchItem[] Matches);

/// <summary>
/// One maximal repeated pair (0-based, forward strand, firstPosition &lt; secondPosition).
/// </summary>
public record MaximalRepeatItem(int FirstPosition, int SecondPosition, int Length);

/// <summary>
/// Result of suffix_tree_maximal_repeats operation.
/// </summary>
/// <param name="Pairs">Maximal repeated pairs sorted by first, then second position.</param>
public record SuffixTreeMaximalRepeatsResult(MaximalRepeatItem[] Pairs);

/// <summary>
/// One longest common substring with all its 0-based start positions in both texts (ascending).
/// </summary>
public record CommonSubstringItem(string Substring, int[] PositionsInText1, int[] PositionsInText2);

/// <summary>
/// Result of suffix_tree_all_lcs operation.
/// </summary>
/// <param name="Substrings">Every longest common substring, ordered by first occurrence in text2.</param>
/// <param name="Length">Common length of the substrings (0 when none).</param>
public record SuffixTreeAllLcsResult(CommonSubstringItem[] Substrings, int Length);

/// <summary>
/// Result of suffix_tree_k_common_substrings operation.
/// </summary>
/// <param name="Substrings">Every longest substring present in at least MinSupport texts, sorted ordinally.</param>
/// <param name="Length">Their length (0 when none).</param>
/// <param name="MinSupport">The support threshold used.</param>
/// <param name="LengthsBySupport">Element q-1 = longest length present in at least q texts (q = 1..k).</param>
public record SuffixTreeKCommonSubstringsResult(string[] Substrings, int Length, int MinSupport, int[] LengthsBySupport);

// ================================
// Genomics Results
// ================================

/// <summary>
/// Result of find_longest_repeat operation.
/// </summary>
/// <param name="Repeat">The repeated string.</param>
/// <param name="Positions">0-based forward-strand starts of the string.</param>
/// <param name="Length">Length of the repeat (0 when nothing repeats).</param>
/// <param name="ReversePositions">With bothStrands: the 0-based forward-strand starts of the windows that spell the reverse complement of the repeat (omitted otherwise).</param>
public record FindLongestRepeatResult(
    string Repeat,
    int[] Positions,
    int Length,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    int[]? ReversePositions = null);

/// <summary>
/// Result of find_longest_common_region operation.
/// </summary>
/// <param name="Region">The common string as it reads in the first sequence (with maxMismatches: the window of the first sequence).</param>
/// <param name="Position1">0-based start in the first sequence.</param>
/// <param name="Position2">0-based forward-strand start in the second sequence (for a reverse-complement region: the window that spells the reverse complement).</param>
/// <param name="Length">Length of the region (0 when there is none).</param>
/// <param name="IsReverseComplement">True when the region is shared with the reverse complement of the second sequence (omitted when false).</param>
/// <param name="Mismatches">With maxMismatches > 0: the number of mismatching positions of the window pair (omitted otherwise).</param>
/// <param name="Region2">With maxMismatches > 0: the window of the second sequence (the reverse complement strand for a reverse region; omitted otherwise).</param>
public record FindLongestCommonRegionResult(
    string Region,
    int Position1,
    int Position2,
    int Length,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingDefault)]
    bool IsReverseComplement = false,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    int? Mismatches = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    string? Region2 = null);

/// <summary>
/// Result of calculate_similarity operation.
/// </summary>
public record CalculateSimilarityResult(
    double Similarity,
    string? Estimator = null,
    double? JaccardDistance = null,
    double? MashDistance = null,
    double? JaccardAni = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? ContainmentAInB = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? ContainmentBInA = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? MaxContainment = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? ContainmentAniAInB = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? ContainmentAniBInA = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    int? SharedHashes = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    int? DenominatorHashes = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? MashPValue = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    AniEstimatesDto? AniEstimates = null);

/// <summary>One sourmash ANI result of calculate_similarity (scaled): distance, interval, ANI and the quality flags.</summary>
public record AniEstimateDto(
    double Distance,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? DistanceLow,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? DistanceHigh,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? Ani,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? AniLow,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? AniHigh,
    double ProbabilityNothingInCommon,
    bool ProbabilityExceedsThreshold,
    bool SizeIsInaccurate,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    double? JaccardError = null,
    [property: System.Text.Json.Serialization.JsonIgnore(Condition = System.Text.Json.Serialization.JsonIgnoreCondition.WhenWritingNull)]
    bool? JaccardErrorExceedsThreshold = null);

/// <summary>The four sourmash ANI results of calculate_similarity with scaled: containment_ani both ways, max_containment_ani, jaccard_ani.</summary>
public record AniEstimatesDto(AniEstimateDto ContainmentAInB, AniEstimateDto ContainmentBInA, AniEstimateDto MaxContainment, AniEstimateDto Jaccard);

/// <summary>
/// Result of hamming_distance operation.
/// </summary>
public record HammingDistanceResult(int Distance);

/// <summary>
/// Result of edit_distance operation.
/// </summary>
public record EditDistanceResult(int Distance);

/// <summary>
/// Result of count_approximate_occurrences operation.
/// </summary>
public record CountApproximateOccurrencesResult(int Count);
