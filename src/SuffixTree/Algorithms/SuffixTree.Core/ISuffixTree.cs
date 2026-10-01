using System;
using System.Collections.Generic;

namespace SuffixTree;
/// <summary>
/// Interface for suffix tree operations.
/// <para>
/// Provides efficient substring search and pattern matching capabilities with O(m) time complexity
/// where m is the pattern length.
/// </para>
/// </summary>
/// <remarks>
/// <para>
/// A suffix tree is a compressed trie of all suffixes of a string. This data structure enables
/// solving many string problems in optimal time:
/// </para>
/// <list type="table">
/// <listheader>
/// <term>Operation</term>
/// <description>Time Complexity</description>
/// </listheader>
/// <item>
/// <term>Substring search</term>
/// <description>O(m)</description>
/// </item>
/// <item>
/// <term>Count occurrences</term>
/// <description>O(m)</description>
/// </item>
/// <item>
/// <term>Find all occurrences</term>
/// <description>O(m + k) where k = number of matches</description>
/// </item>
/// <item>
/// <term>Longest repeated substring</term>
/// <description>O(1) cached</description>
/// </item>
/// <item>
/// <term>Longest common substring</term>
/// <description>O(m)</description>
/// </item>
/// </list>
/// </remarks>
public interface ISuffixTree : ISuffixTreeSearch, ISuffixTreeAnalysis, ISuffixTreeDiagnostics
{
    /// <summary>
    /// Counts distinct substrings by length: element <c>i</c> is the number of distinct
    /// substrings of length <c>i</c> (<c>i = 1..min(maxLength, n)</c>); element 0 is 1.
    /// The default implementation traverses the tree; implementations may override.
    /// </summary>
    long[] CountDistinctSubstringsByLength(int maxLength)
        => SuffixTreeAlgorithms.CountDistinctSubstringsByLength(this, maxLength);

    /// <summary>
    /// Total number of distinct non-empty substrings of the text.
    /// </summary>
    long CountDistinctSubstrings() => SuffixTreeAlgorithms.CountDistinctSubstrings(this);

    /// <summary>
    /// Every distinct longest repeated substring (all length ties, unlike the single representative
    /// of <see cref="ISuffixTreeAnalysis.LongestRepeatedSubstring"/>), each with all 0-based start
    /// positions ascending (occurrences may overlap), ordered by first occurrence; empty when no
    /// character repeats. Identical for every implementation
    /// (see <see cref="SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings"/>).
    /// </summary>
    IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> FindAllLongestRepeatedSubstrings()
        => SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings(this);

    /// <summary>
    /// Every maximal repeated pair (FirstPosition i &lt; SecondPosition j, Length L ≥ <paramref name="minLength"/>)
    /// of the text: text[i..i+L) = text[j..j+L), not extendable to the left nor to the right (Gusfield 1997
    /// §7.12; MUMmer 3 <c>repeat-match -f -n minLength</c>), 0-based, forward strand, copies may overlap.
    /// Characters for which <paramref name="isUniqueSymbol"/> returns true never match (not even
    /// themselves). Ordered by FirstPosition, then SecondPosition; identical for every implementation
    /// (see <see cref="SuffixTreeAlgorithms.FindMaximalRepeatedPairs"/>).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minLength"/> &lt; 1.</exception>
    IReadOnlyList<(int FirstPosition, int SecondPosition, int Length)> FindMaximalRepeatedPairs(
        int minLength, Func<char, bool>? isUniqueSymbol = null)
        => SuffixTreeAlgorithms.FindMaximalRepeatedPairs(this, minLength, isUniqueSymbol);
}

/// <summary>
/// Visitor interface for deterministic tree traversal.
/// </summary>
public interface ISuffixTreeVisitor
{
    /// <summary>
    /// Called when entering a node.
    /// </summary>
    void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth);

    /// <summary>
    /// Called before visiting a child branch.
    /// </summary>
    void EnterBranch(int key);

    /// <summary>
    /// Called after visiting a child branch.
    /// </summary>
    void ExitBranch();
}
