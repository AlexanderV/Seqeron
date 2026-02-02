using System;
using System.Diagnostics.CodeAnalysis;

namespace SuffixTree
{
    /// <summary>
    /// High-performance suffix tree implementation using Ukkonen's online algorithm.
    /// <para>
    /// A suffix tree is a compressed trie of all suffixes of a string, enabling O(m) substring
    /// search where m is the pattern length, regardless of the text size.
    /// </para>
    /// <para>
    /// <b>Thread Safety:</b> The tree is immutable after construction and safe for concurrent reads.
    /// Construction is not thread-safe.
    /// </para>
    /// <para>
    /// <b>Performance:</b>
    /// <list type="bullet">
    /// <item>Construction: O(n) time, O(n) space</item>
    /// <item>Contains: O(m) time, zero allocations</item>
    /// <item>CountOccurrences: O(m) time using precomputed leaf counts</item>
    /// <item>FindAllOccurrences: O(m + k) time where k is number of occurrences</item>
    /// <item>LongestRepeatedSubstring: O(1) after first call (cached)</item>
    /// <item>LongestCommonSubstring: O(m) using suffix link traversal</item>
    /// </list>
    /// </para>
    /// </summary>
    /// <example>
    /// <code>
    /// // Build a suffix tree
    /// var tree = SuffixTree.Build("banana");
    /// 
    /// // Search for substrings
    /// bool found = tree.Contains("nan");  // true
    /// int count = tree.CountOccurrences("ana");  // 2
    /// var positions = tree.FindAllOccurrences("ana");  // [1, 3]
    /// 
    /// // Find patterns
    /// string lrs = tree.LongestRepeatedSubstring();  // "ana"
    /// string lcs = tree.LongestCommonSubstring("bandana");  // "ana"
    /// </code>
    /// </example>
    public partial class SuffixTree : ISuffixTree
    {
        private const int TERMINATOR_KEY = -1;
        private const int MAX_TOSTRING_CONTENT_LENGTH = 50;

        private static readonly SuffixTree _empty = new SuffixTree(new StringTextSource(string.Empty));

        private readonly ITextSource _text;
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
        /// <param name="value">The text to build the suffix tree from. Cannot be null.</param>
        /// <returns>A new suffix tree instance.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="value"/> is null.</exception>
        /// <remarks>
        /// Construction time is O(n) where n is the length of the input string.
        /// The tree uses approximately 20 bytes per character for typical texts.
        /// </remarks>
        /// <example>
        /// <code>
        /// var tree = SuffixTree.Build("mississippi");
        /// Console.WriteLine(tree.NodeCount);  // Number of nodes in tree
        /// </code>
        /// </example>
        public static SuffixTree Build(string value)
        {
            return new SuffixTree(new StringTextSource(value));
        }

        /// <summary>
        /// Creates and returns a suffix tree for the specified text source.
        /// </summary>
        /// <param name="value">The text source to build the suffix tree from. Cannot be null.</param>
        /// <returns>A new suffix tree instance.</returns>
        public static SuffixTree Build(ITextSource value)
        {
            return new SuffixTree(value);
        }

        /// <summary>
        /// Attempts to create a suffix tree for the specified string.
        /// </summary>
        /// <param name="value">The text to build the suffix tree from.</param>
        /// <param name="tree">When this method returns, contains the suffix tree if successful; otherwise, null.</param>
        /// <returns>True if the suffix tree was created successfully; otherwise, false.</returns>
        /// <remarks>
        /// This method never throws exceptions. Use this when you want to handle
        /// null or invalid input without exception handling overhead.
        /// </remarks>
        /// <example>
        /// <code>
        /// if (SuffixTree.TryBuild(userInput, out var tree))
        /// {
        ///     bool found = tree.Contains("pattern");
        /// }
        /// else
        /// {
        ///     Console.WriteLine("Invalid input");
        /// }
        /// </code>
        /// </example>
        public static bool TryBuild(string? value, [NotNullWhen(true)] out SuffixTree? tree)
        {
            if (value == null)
            {
                tree = null;
                return false;
            }

            tree = new SuffixTree(new StringTextSource(value));
            return true;
        }

        /// <summary>
        /// Attempts to create a suffix tree for the specified text source.
        /// </summary>
        /// <param name="value">The text source to build the suffix tree from.</param>
        /// <param name="tree">When this method returns, contains the suffix tree if successful; otherwise, null.</param>
        /// <returns>True if the suffix tree was created successfully; otherwise, false.</returns>
        public static bool TryBuild(ITextSource? value, [NotNullWhen(true)] out SuffixTree? tree)
        {
            if (value == null)
            {
                tree = null;
                return false;
            }

            tree = new SuffixTree(value);
            return true;
        }

        /// <summary>
        /// Gets an empty suffix tree (singleton instance).
        /// </summary>
        /// <value>A cached empty suffix tree instance.</value>
        /// <remarks>
        /// An empty tree has no leaves and Contains() returns true only for empty string.
        /// This property returns the same instance on every access.
        /// </remarks>
        public static SuffixTree Empty => _empty;

        /// <summary>
        /// Creates a suffix tree from a memory region containing characters.
        /// </summary>
        /// <param name="value">The memory region containing the text to build the suffix tree from.</param>
        /// <returns>A new suffix tree instance.</returns>
        /// <remarks>
        /// This overload is useful when working with buffers or memory pools.
        /// The content is copied to a string internally.
        /// </remarks>
        /// <example>
        /// <code>
        /// Memory&lt;char&gt; buffer = GetBufferFromPool();
        /// var tree = SuffixTree.Build(buffer);
        /// </code>
        /// </example>
        public static SuffixTree Build(ReadOnlyMemory<char> value)
        {
            return new SuffixTree(new StringTextSource(new string(value.Span)));
        }

        /// <summary>
        /// Creates a suffix tree from a character span.
        /// </summary>
        /// <param name="value">The character span to build the suffix tree from.</param>
        /// <returns>A new suffix tree instance.</returns>
        /// <remarks>
        /// This overload is useful for zero-allocation scenarios where you have a span.
        /// The content is copied to a string internally for storage.
        /// </remarks>
        public static SuffixTree Build(ReadOnlySpan<char> value)
        {
            return new SuffixTree(new StringTextSource(new string(value)));
        }

        private SuffixTree(ITextSource value)
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

        /// <inheritdoc/>
        public ITextSource Text => _text;

        /// <inheritdoc/>
        public int NodeCount => _cachedNodeCount;

        /// <inheritdoc/>
        public int LeafCount => _cachedLeafCount;

        /// <inheritdoc/>
        public int MaxDepth => _cachedMaxDepth;

        /// <summary>
        /// Gets a value indicating whether this tree is empty (built from an empty string).
        /// </summary>
        public bool IsEmpty => _text.Length == 0;

        /// <summary>
        /// Gets the total depth from root to the END of this node's edge.
        /// </summary>
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        private int GetNodeDepth(SuffixTreeNode? node)
        {
            if (node == null || node == _root) return 0;
            return node.DepthFromRoot + LengthOf(node);
        }

        /// <summary>
        /// Returns a string representation of this suffix tree.
        /// </summary>
        /// <returns>A string containing node count, leaf count, and a preview of the text.</returns>
        public override string ToString()
        {
            if (_text.Length == 0)
                return "SuffixTree (empty)";

            var content = _text.Length <= MAX_TOSTRING_CONTENT_LENGTH
                ? _text.ToString()
                : string.Concat(_text.Slice(0, MAX_TOSTRING_CONTENT_LENGTH), "...");
            return $"SuffixTree (Nodes: {NodeCount}, Leaves: {LeafCount}, Text: \"{content}\")";
        }
    }
}
