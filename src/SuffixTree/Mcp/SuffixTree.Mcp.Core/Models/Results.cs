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

// ================================
// Genomics Results
// ================================

/// <summary>
/// Result of find_longest_repeat operation.
/// </summary>
public record FindLongestRepeatResult(string Repeat, int[] Positions, int Length);

/// <summary>
/// Result of find_longest_common_region operation.
/// </summary>
public record FindLongestCommonRegionResult(string Region, int Position1, int Position2, int Length);

/// <summary>
/// Result of calculate_similarity operation.
/// </summary>
public record CalculateSimilarityResult(double Similarity);

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
