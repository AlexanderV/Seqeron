using System;

namespace SuffixTree
{
    /// <summary>
    /// Implementation of Ukkonen's online suffix tree construction algorithm.
    /// Thread-safe for multiple readers once constructed.
    /// </summary>
    public partial class SuffixTree : ISuffixTree
    {
        private const int TERMINATOR_KEY = -1;
        private const int MAX_TOSTRING_CONTENT_LENGTH = 50;

        private readonly string _text;
        private readonly SuffixTreeNode _root;
        private readonly int _cachedNodeCount;
        private readonly int _cachedLeafCount;
        private readonly int _cachedMaxDepth;
        private readonly SuffixTreeNode? _deepestInternalNode;
        private readonly int _maxInternalDepth;

        // Construction state (transient)
        private int _remainder;
        private int _position = -1;
        private SuffixTreeNode? _lastCreatedInternalNode;
        private SuffixTreeNode? _activeNode;
        private int _activeEdgeIndex = -1;
        private int _activeLength;

        /// <summary>
        /// Creates and returns a suffix tree for the specified string.
        /// </summary>
        public static SuffixTree Build(string value)
        {
            return new SuffixTree(value);
        }

        private SuffixTree(string value)
        {
            ArgumentNullException.ThrowIfNull(value);
            _text = value;
            _root = new SuffixTreeNode(0, 0, 0);
            _root.SuffixLink = _root;
            _activeNode = _root;

            if (value.Length > 0)
            {
                // Run Ukkonen's algorithm
                foreach (var c in value)
                    ExtendTree(c);

                ExtendTree(TERMINATOR_KEY);
            }

            // Calculate statistics once and store in readonly fields
            var (nodeCount, leafCount, maxDepth, deepestInternalNode, maxInternalDepth) = CalculateTreeStatisticsInternal();
            _cachedNodeCount = nodeCount;
            _cachedLeafCount = leafCount;
            _cachedMaxDepth = maxDepth;
            _deepestInternalNode = deepestInternalNode;
            _maxInternalDepth = maxInternalDepth;

#if DEBUG
            ValidateSuffixLinks();
#endif

            // Clear construction state to help GC
            _lastCreatedInternalNode = null;
            _activeNode = null;
            _activeEdgeIndex = -1;
            _activeLength = 0;
            _remainder = 0;
        }

        public string Text => _text;
        public int NodeCount => _cachedNodeCount;
        public int LeafCount => _cachedLeafCount;
        public int MaxDepth => _cachedMaxDepth;

        /// <summary>
        /// Gets the total depth from root to the END of this node's edge.
        /// </summary>
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        private int GetNodeDepth(SuffixTreeNode? node)
        {
            if (node == null || node == _root) return 0;
            return node.DepthFromRoot + LengthOf(node);
        }

        public override string ToString()
        {
            if (_text.Length == 0)
                return "SuffixTree (empty)";

            var content = _text.Length <= MAX_TOSTRING_CONTENT_LENGTH
                ? _text
                : string.Concat(_text.AsSpan(0, MAX_TOSTRING_CONTENT_LENGTH), "...");
            return $"SuffixTree (Nodes: {NodeCount}, Leaves: {LeafCount}, Text: \"{content}\")";
        }
    }
}
