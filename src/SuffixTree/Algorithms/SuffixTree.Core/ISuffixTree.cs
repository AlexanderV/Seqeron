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
