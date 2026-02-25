namespace SuffixTree;

/// <summary>
/// Diagnostics and deterministic traversal contract for suffix trees.
/// </summary>
public interface ISuffixTreeDiagnostics
{
    /// <summary>
    /// Creates a detailed string representation of the tree structure.
    /// Useful for debugging and visualization.
    /// </summary>
    string PrintTree();

    /// <summary>
    /// Performs a deterministic traversal of the tree nodes.
    /// </summary>
    void Traverse(ISuffixTreeVisitor visitor);
}
