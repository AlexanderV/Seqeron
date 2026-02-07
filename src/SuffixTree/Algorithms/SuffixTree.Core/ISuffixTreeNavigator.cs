using System.Collections.Generic;

namespace SuffixTree
{
    /// <summary>
    /// Abstraction for navigating a suffix tree's node structure.
    /// Used by <see cref="SuffixTreeAlgorithms"/> to implement shared algorithms
    /// (LCS, FindExactMatchAnchors) that work with any node representation.
    /// <para>
    /// Implementations should be <c>struct</c> to enable JIT specialization
    /// and zero-overhead generic dispatch.
    /// </para>
    /// </summary>
    /// <typeparam name="TNode">The node type (e.g., SuffixTreeNode or PersistentSuffixTreeNode).</typeparam>
    public interface ISuffixTreeNavigator<TNode>
    {
        /// <summary>Gets the text source backing this tree.</summary>
        ITextSource Text { get; }

        /// <summary>Gets the root node of the tree.</summary>
        TNode Root { get; }

        /// <summary>Gets a sentinel value representing "no node".</summary>
        TNode NullNode { get; }

        /// <summary>Returns true if the node is the null sentinel.</summary>
        bool IsNull(TNode node);

        /// <summary>Returns true if the node is the root.</summary>
        bool IsRoot(TNode node);

        /// <summary>
        /// Gets the character/symbol at the given position on this node's edge.
        /// Returns -1 for the implicit terminator.
        /// </summary>
        int GetEdgeSymbol(TNode node, int offset);

        /// <summary>Gets the length of this node's edge label.</summary>
        int LengthOf(TNode node);

        /// <summary>
        /// Gets the total depth from root to the END of this node's edge
        /// (DepthFromRoot + LengthOf).
        /// </summary>
        int GetNodeDepth(TNode node);

        /// <summary>
        /// Gets the depth from root to the START of this node's edge (DepthFromRoot).
        /// </summary>
        int GetDepthFromRoot(TNode node);

        /// <summary>
        /// Gets the suffix link target. Returns Root if no suffix link exists.
        /// </summary>
        TNode GetSuffixLink(TNode node);

        /// <summary>
        /// Tries to get a child node by edge key (character as int, -1 for terminator).
        /// </summary>
        bool TryGetChild(TNode node, int key, out TNode child);

        /// <summary>
        /// Collects all leaf positions in the subtree rooted at <paramref name="node"/>.
        /// </summary>
        void CollectLeaves(TNode node, int depth, List<int> results);

        /// <summary>
        /// Finds any single leaf position in the subtree rooted at <paramref name="node"/>.
        /// Returns -1 if not found.
        /// </summary>
        int FindAnyLeafPosition(TNode node);
    }
}
