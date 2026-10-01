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
