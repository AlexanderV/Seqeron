using System.Collections.Generic;

namespace SuffixTree
{
    /// <summary>
    /// Internal node representation for suffix tree.
    /// Each edge is implicitly stored as (Start, End) indices into the character array.
    /// </summary>
    internal class SuffixTreeNode
    {
        /// <summary>
        /// Sentinel value indicating an open-ended (growing) edge.
        /// Leaf edges use this to automatically extend as new characters are added.
        /// </summary>
        internal const int BOUNDLESS = -1;

        /// <summary>Start index of edge label in character array.</summary>
        public int Start { get; set; }

        /// <summary>
        /// End index of edge label (exclusive). 
        /// BOUNDLESS (-1) means this is a leaf edge that grows automatically.
        /// </summary>
        public int End { get; set; }

        /// <summary>Children edges, keyed by first character.</summary>
        public Dictionary<char, SuffixTreeNode> Children { get; } = new Dictionary<char, SuffixTreeNode>();

        /// <summary>
        /// Suffix link: connects node for "xα" to node for "α".
        /// Used for O(1) jumps between suffixes during construction.
        /// </summary>
        public SuffixTreeNode SuffixLink { get; set; }

        /// <summary>True if this is a leaf node (edge grows with string).</summary>
        public bool IsLeaf => End == BOUNDLESS;
    }
}
