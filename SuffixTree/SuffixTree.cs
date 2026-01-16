using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;

namespace SuffixTree
{
    /// <summary>
    /// Implementation of Ukkonen's online suffix tree construction algorithm.
    ///
    /// A suffix tree is a compressed trie containing all suffixes of a given string.
    /// This implementation builds the tree in O(n) time and allows O(m) substring search,
    /// where n is the text length and m is the pattern length.
    ///
    /// Key concepts:
    /// - Active Point: tracks current position (activeNode, activeEdge, activeLength)
    /// - Suffix Links: enable O(1) jumps between internal nodes
    /// - Remainder: counts pending suffixes to be inserted
    /// - Terminator: ensures all suffixes are explicit (end at leaves)
    /// </summary>
    public class SuffixTree : ISuffixTree
    {
        /// <summary>
        /// Unique terminator character appended to ensure all suffixes are explicit.
        /// Without this, some suffixes may remain implicit (e.g., "a" in "aa").
        /// </summary>
        private const char TERMINATOR = '\0';

        /// <summary>
        /// Maximum content length displayed in ToString() before truncation.
        /// </summary>
        private const int MAX_TOSTRING_CONTENT_LENGTH = 50;

        /// <summary>Number of suffixes that still need explicit insertion.</summary>
        private int _remainder;

        /// <summary>Current position in _chars (0-based index of last added character).</summary>
        private int _position = -1;

        private readonly SuffixTreeNode _root;

        /// <summary>
        /// Last internal node created in current phase.
        /// Used to set up suffix link chain: when we create a new internal node,
        /// we link the previous one to it.
        /// </summary>
        private SuffixTreeNode? _lastCreatedInternalNode;

        /// <summary>Deepest internal node found during construction (for O(1) LRS).</summary>
        private SuffixTreeNode? _deepestInternalNode;

        /// <summary>Depth of the deepest internal node (sum of edge lengths from root).</summary>
        private int _maxInternalDepth;

        /// <summary>Cached total node count (calculated once after construction).</summary>
        private int _cachedNodeCount;

        /// <summary>Cached total leaf count (calculated once after construction).</summary>
        private int _cachedLeafCount;

        /// <summary>Cached max tree depth (calculated once after construction).</summary>
        private int _cachedMaxDepth;

        /// <summary>The string content stored as a character array.</summary>
        private char[] _chars = [];

        // ============================================================
        // Active Point - tracks where we are in the tree during construction
        // ============================================================
        // When Rule 3 (showstopper) fires, we don't insert but remember WHERE
        // the suffix exists implicitly. The active point tracks this location.
        // ============================================================

        /// <summary>The node we're currently at (or descended from).</summary>
        private SuffixTreeNode? _activeNode;

        /// <summary>Index in _chars of the first character of the active edge (-1 = none).</summary>
        private int _activeEdgeIndex = -1;

        /// <summary>How many characters along the active edge we've matched.</summary>
        private int _activeLength;

        /// <summary>Cached original text (without terminator).</summary>
        private string _text = string.Empty;

        /// <summary>
        /// Gets the original text that this suffix tree was built from (without terminator).
        /// </summary>
        public string Text => _text;

        /// <summary>
        /// Gets the total number of nodes in the tree (including root).
        /// Time complexity: O(1) - calculated once after construction.
        /// </summary>
        public int NodeCount => _cachedNodeCount;

        /// <summary>
        /// Gets the number of leaf nodes in the tree.
        /// Equal to the length of the original text (one leaf per suffix).
        /// Time complexity: O(1) - calculated once after construction.
        /// </summary>
        public int LeafCount => _cachedLeafCount;

        /// <summary>
        /// Gets the maximum depth of the tree (longest path from root to leaf in characters).
        /// Time complexity: O(1) - calculated once after construction.
        /// </summary>
        public int MaxDepth => _cachedMaxDepth;

        /// <summary>
        /// Calculates all tree statistics in a single post-order DFS pass.
        /// Sets LeafCount for each node and caches NodeCount, LeafCount, MaxDepth.
        /// Time complexity: O(n) where n is the number of nodes.
        /// </summary>
        private void CalculateTreeStatistics()
        {
            if (_chars == null || _chars.Length == 0)
            {
                _cachedNodeCount = 1; // Just root
                _cachedLeafCount = 0;
                _cachedMaxDepth = 0;
                return;
            }

            // Post-order DFS to calculate leaf counts bottom-up
            // We need to process children before parents
            var stack = new Stack<(SuffixTreeNode Node, int Depth, bool Visited)>();
            stack.Push((_root, 0, false));

            int nodeCount = 0;
            int leafCount = 0;
            int maxDepth = 0;

            while (stack.Count > 0)
            {
                var (node, depth, visited) = stack.Pop();

                if (visited)
                {
                    // Post-order: all children processed, calculate this node's leaf count
                    nodeCount++;

                    if (node.IsLeaf)
                    {
                        node.LeafCount = 1;
                        leafCount++;
                        if (depth > maxDepth) maxDepth = depth;
                    }
                    else
                    {
                        // Sum of all children's leaf counts
                        int sum = 0;
                        foreach (var child in node.GetChildren())
                            sum += child.LeafCount;
                        node.LeafCount = sum;
                    }
                }
                else
                {
                    // Pre-order: push self for post-processing, then push children
                    int currentDepth = depth + (node == _root ? 0 : LengthOf(node));
                    stack.Push((node, currentDepth, true));

                    if (node.HasChildren)
                    {
                        foreach (var child in node.GetChildren())
                            stack.Push((child, currentDepth, false));
                    }
                }
            }

            _cachedNodeCount = nodeCount;
            _cachedLeafCount = leafCount;
            // Subtract 1 for terminator character (not counted in MaxDepth)
            _cachedMaxDepth = maxDepth > 0 ? maxDepth - 1 : 0;
        }

        private SuffixTree()
        {
            _root = new SuffixTreeNode { Start = 0, End = 0 };
            _root.SuffixLink = _root; // Root's suffix link points to itself
            _activeNode = _root;
        }

        /// <summary>
        /// Creates and returns a suffix tree for the specified string.
        /// </summary>
        /// <param name="value">The string to build the tree from.</param>
        /// <returns>A new SuffixTree instance.</returns>
        /// <exception cref="ArgumentNullException">If value is null.</exception>
        public static SuffixTree Build(string value)
        {
            var t = new SuffixTree();
            t.BuildInternal(value);
            return t;
        }

        /// <summary>
        /// Internal method to construct the suffix tree from a string.
        /// A unique terminator character is automatically appended to ensure
        /// all suffixes are explicit (end at leaf nodes).
        /// </summary>
        /// <param name="value">The string to build tree from.</param>
        /// <exception cref="ArgumentNullException">If value is null.</exception>
        /// <exception cref="ArgumentException">If value contains the null character '\0'.</exception>
        private void BuildInternal(string value)
        {
            ArgumentNullException.ThrowIfNull(value);

            if (value.Contains(TERMINATOR))
                throw new ArgumentException(
                    $"Input string cannot contain the null character '\\0' as it is used as internal terminator.",
                    nameof(value));

            if (value.Length == 0)
            {
                _chars = Array.Empty<char>();
                _text = string.Empty;
                CalculateTreeStatistics(); // Sets _cachedNodeCount = 1 for root
                return;
            }

            // Pre-allocate array for string + terminator
            _chars = new char[value.Length + 1];

            // Process each character - Ukkonen's online construction
            foreach (var c in value)
                ExtendTree(c);

            // Add terminator to convert implicit suffixes to explicit leaves
            ExtendTree(TERMINATOR);

            // Cache the original text (without terminator)
            _text = value;

            // Calculate leaf counts and statistics (one-time O(n) pass)
            CalculateTreeStatistics();

            // Clear construction state (not needed anymore, helps GC)
            _lastCreatedInternalNode = null;
            _activeNode = null;
            _activeEdgeIndex = -1;
            _activeLength = 0;
            _remainder = 0;
        }

        /// <summary>
        /// Extends the suffix tree with one character.
        /// This is the core of Ukkonen's algorithm - each call processes one "phase".
        ///
        /// In each phase, we must extend ALL current suffixes with the new character.
        /// The key insight is that we can do this in amortized O(1) per suffix using:
        /// - Rule 3 (showstopper): if suffix already exists, stop immediately
        /// - Suffix links: jump between suffixes in O(1) instead of walking from root
        /// - Open-ended leaves: automatically extend without explicit updates
        /// </summary>
        private void ExtendTree(char c)
        {
            _position++;
            _chars[_position] = c;
            _remainder++;
            _lastCreatedInternalNode = null;

            // Process all pending suffixes (remainder count)
            while (_remainder > 0)
            {
                // If activeLength is 0, we're exactly at a node
                // Set active edge to start with current character
                if (_activeLength == 0)
                    _activeEdgeIndex = _position;

                char activeEdgeChar = _chars[_activeEdgeIndex];

                if (!_activeNode!.TryGetChild(activeEdgeChar, out SuffixTreeNode? activeEdge))
                {
                    // ============================================================
                    // RULE 1: No edge starts with this character - create new leaf
                    // ============================================================
                    // This happens when we're at a node and need to add a new suffix
                    // that doesn't share any prefix with existing edges.
                    //
                    //      [node]                    [node]
                    //         |           →            |  \
                    //       (existing)              (existing) c
                    // ============================================================
                    var leaf = new SuffixTreeNode { Start = _position, End = SuffixTreeNode.BOUNDLESS, Parent = _activeNode };
                    _activeNode.SetChild(c, leaf);

                    // Set suffix link for previously created internal node
                    AddSuffixLink(_activeNode);
                }
                else
                {
                    // Edge exists - need to check if we should walk down first
                    if (WalkDownIfNeeded(activeEdge!))
                        continue; // Active point moved, restart this iteration

                    // Check if next character on edge matches
                    if (_chars[activeEdge!.Start + _activeLength] == c)
                    {
                        // ============================================================
                        // RULE 3: Character matches - suffix exists implicitly (SHOWSTOPPER)
                        // ============================================================
                        // The suffix we're trying to insert already exists in the tree!
                        // This is the key optimization: we stop processing immediately.
                        //
                        // Why? All remaining suffixes are SHORTER, so if this one exists,
                        // they must also exist (as prefixes of this one).
                        //
                        // We don't decrement remainder - these suffixes are still "pending"
                        // and will be handled when a mismatch occurs in a future phase.
                        // ============================================================
                        AddSuffixLink(_activeNode);
                        _activeLength++;
                        break; // STOP - this is the "showstopper"
                    }

                    // ============================================================
                    // RULE 2: Character doesn't match - split the edge
                    // ============================================================
                    // We're in the middle of an edge and found a mismatch.
                    // Need to create a new internal node at this position.
                    //
                    //      [node]                    [node]
                    //         |                         |
                    //       "abc"          →          "ab"     (new internal node)
                    //         |                      /    \
                    //       [leaf]               "c"      "x"  (new leaf + old suffix)
                    //                          [new]     [old]
                    // ============================================================

                    // Calculate depth for the new split node
                    int splitNodeDepth = GetNodeDepth(_activeNode) + _activeLength;

                    var splitNode = new SuffixTreeNode
                    {
                        Start = activeEdge.Start,
                        End = activeEdge.Start + _activeLength,
                        Parent = _activeNode,
                        DepthFromRoot = GetNodeDepth(_activeNode)
                    };

                    // Replace old edge with new internal node
                    _activeNode.SetChild(activeEdgeChar, splitNode);

                    // Create new leaf for current character
                    var newLeaf = new SuffixTreeNode { Start = _position, End = SuffixTreeNode.BOUNDLESS, Parent = splitNode };
                    splitNode.SetChild(c, newLeaf);

                    // Move old edge to be child of split node
                    activeEdge.Start += _activeLength;
                    activeEdge.Parent = splitNode;
                    splitNode.SetChild(_chars[activeEdge.Start], activeEdge);

                    // Update deepest internal node tracking for O(1) LRS
                    if (splitNodeDepth > _maxInternalDepth)
                    {
                        _maxInternalDepth = splitNodeDepth;
                        _deepestInternalNode = splitNode;
                    }

                    // Track this internal node for suffix link setup
                    SetLastCreatedInternalNode(splitNode);
                }

                // Successfully inserted a suffix - decrement counter
                _remainder--;

                // Move active point to next suffix
                if (_activeNode == _root && _activeLength > 0)
                {
                    // ============================================================
                    // At root: adjust active point for next shorter suffix
                    // ============================================================
                    // When at root, we can't follow a suffix link.
                    // Instead, we manually adjust: the next suffix to process
                    // is one character shorter, starting one position later.
                    // ============================================================
                    _activeLength--;
                    _activeEdgeIndex = _position - _remainder + 1;
                }
                else if (_activeNode != _root)
                {
                    // ============================================================
                    // Not at root: follow suffix link
                    // ============================================================
                    // Suffix link connects node for "xα" to node for "α".
                    // This lets us jump to the next suffix in O(1) time!
                    // If no suffix link exists, fall back to root.
                    // ============================================================
                    _activeNode = _activeNode.SuffixLink ?? _root;
                }
            }
        }

        /// <summary>
        /// Walk down the tree if activeLength exceeds current edge length.
        ///
        /// This handles the case where our "position" conceptually extends
        /// past the current edge into a child edge. We need to move activeNode
        /// to the end of this edge and adjust activeLength accordingly.
        ///
        /// Example: If activeLength=5 but edge "abc" has length 3:
        ///   Before: activeNode=parent, activeEdge='a', activeLength=5
        ///   After:  activeNode=node_after_abc, activeLength=2
        /// </summary>
        /// <returns>True if we moved (caller should restart loop), false otherwise.</returns>
        private bool WalkDownIfNeeded(SuffixTreeNode edge)
        {
            int edgeLength = LengthOf(edge);
            if (_activeLength >= edgeLength)
            {
                _activeEdgeIndex += edgeLength;
                _activeLength -= edgeLength;
                _activeNode = edge;
                return true;
            }
            return false;
        }

        /// <summary>
        /// Sets suffix link from previously created internal node to the given node.
        /// Called after each extension to maintain suffix link invariants.
        /// </summary>
        private void AddSuffixLink(SuffixTreeNode node)
        {
            if (_lastCreatedInternalNode != null)
            {
                _lastCreatedInternalNode.SuffixLink = node;
                _lastCreatedInternalNode = null;
            }
        }

        /// <summary>
        /// Tracks a newly created internal node for suffix linking.
        /// The next internal node created (or activeNode if Rule 1/3) will receive
        /// a suffix link FROM this node.
        /// </summary>
        private void SetLastCreatedInternalNode(SuffixTreeNode node)
        {
            if (_lastCreatedInternalNode != null)
                _lastCreatedInternalNode.SuffixLink = node;

            _lastCreatedInternalNode = node;
        }

        /// <summary>
        /// Gets the total depth from root to the END of this node's edge.
        /// For root, returns 0.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int GetNodeDepth(SuffixTreeNode node)
        {
            if (node == _root) return 0;
            return node.DepthFromRoot + LengthOf(node);
        }

        /// <summary>
        /// Calculates the length of an edge.
        /// For leaves (End == BOUNDLESS), the edge extends to current position.
        /// </summary>
        private int LengthOf(SuffixTreeNode edge)
            => (edge.End == SuffixTreeNode.BOUNDLESS ? _position + 1 : edge.End) - edge.Start;

        /// <summary>Returns the first character of an edge's label.</summary>
        private char FirstCharOf(SuffixTreeNode edge)
            => _chars[edge.Start];

        /// <summary>Returns the full label of an edge (for debugging).</summary>
        private string LabelOf(SuffixTreeNode edge)
        {
            var length = LengthOf(edge);
            var sb = new StringBuilder(length);
            for (int i = edge.Start; i < edge.Start + length; i++)
            {
                var c = _chars[i];
                sb.Append(c == TERMINATOR ? '$' : c);
            }
            return sb.ToString();
        }

        /// <summary>
        /// Checks if the specified string is a substring of the tree content.
        ///
        /// Algorithm: Walk down the tree following edges that match the pattern.
        /// If we can match all characters, the substring exists.
        ///
        /// Time complexity: O(m) where m is the length of the search string.
        /// </summary>
        /// <param name="value">The substring to search for.</param>
        /// <returns>True if the substring exists, false otherwise.</returns>
        /// <exception cref="ArgumentNullException">If value is null.</exception>
        public bool Contains(string value)
        {
            ArgumentNullException.ThrowIfNull(value);

            return Contains(value.AsSpan());
        }

        /// <summary>
        /// Checks if the specified character span is a substring of the tree content.
        /// This overload avoids string allocation when searching with spans or slices.
        ///
        /// Time complexity: O(m) where m is the length of the search span.
        /// </summary>
        /// <param name="value">The character span to search for.</param>
        /// <returns>True if the substring exists, false otherwise.</returns>
        /// <exception cref="ArgumentException">If value contains the null character '\0'.</exception>
        public bool Contains(ReadOnlySpan<char> value)
        {
            if (value.Length == 0)
                return true; // Empty string is always a substring

            // Check for terminator character in input
            if (value.Contains(TERMINATOR))
                throw new ArgumentException(
                    "Search value cannot contain the null character '\\0' as it is used as internal terminator.",
                    nameof(value));

            var node = _root;
            int i = 0;

            while (i < value.Length)
            {
                // Try to find edge starting with current character
                if (!node.TryGetChild(value[i], out var child) || child is null)
                    return false; // No matching edge - substring doesn't exist

                int edgeLength = LengthOf(child);
                int j = 0;

                // Match characters along this edge
                while (j < edgeLength && i < value.Length)
                {
                    if (_chars[child.Start + j] != value[i])
                        return false; // Mismatch - substring doesn't exist
                    i++;
                    j++;
                }

                // If we consumed the entire edge, move to child node
                if (j == edgeLength)
                    node = child;
                // Otherwise we matched in middle of edge - that's fine, we're done
            }

            return true; // All characters matched
        }

        /// <summary>
        /// Finds all starting positions where the pattern occurs in the original string.
        ///
        /// Algorithm:
        /// 1. Navigate to the node representing the pattern (like Contains)
        /// 2. Collect all leaf positions in the subtree (each leaf = one occurrence)
        ///
        /// Time complexity: O(m + k) where m is pattern length and k is number of occurrences.
        /// </summary>
        /// <param name="pattern">The pattern to search for.</param>
        /// <returns>Collection of 0-based starting positions of all occurrences.</returns>
        /// <exception cref="ArgumentNullException">If pattern is null.</exception>
        public IReadOnlyList<int> FindAllOccurrences(string pattern)
        {
            ArgumentNullException.ThrowIfNull(pattern);

            return FindAllOccurrences(pattern.AsSpan());
        }

        /// <summary>
        /// Finds all starting positions where the pattern occurs in the original string.
        /// Zero-allocation overload for performance-critical scenarios.
        ///
        /// Algorithm:
        /// 1. Navigate down the tree following the pattern
        /// 2. Collect all leaf positions in the subtree (each leaf = one occurrence)
        ///
        /// Time complexity: O(m + k) where m is pattern length and k is number of occurrences.
        /// </summary>
        /// <param name="pattern">The pattern to search for.</param>
        /// <returns>Collection of 0-based starting positions of all occurrences.</returns>
        /// <exception cref="ArgumentException">If pattern contains the null character '\0'.</exception>
        public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
        {
            if (pattern.Contains(TERMINATOR))
                throw new ArgumentException(
                    $"Pattern cannot contain the null character '\\0' as it is used as internal terminator.",
                    nameof(pattern));

            var results = new List<int>();

            if (pattern.Length == 0)
            {
                // Empty pattern matches at every position in original text
                int textLength = _chars.Length - 1; // Exclude terminator
                for (int p = 0; p < textLength; p++)
                    results.Add(p);
                return results;
            }

            // Navigate to the node/edge representing the pattern
            var node = _root;
            int i = 0;
            int depthFromRoot = 0; // Total characters matched from root

            while (i < pattern.Length)
            {
                if (!node.TryGetChild(pattern[i], out var child) || child is null)
                    return results; // Pattern not found

                int edgeLength = LengthOf(child);
                int j = 0;

                while (j < edgeLength && i < pattern.Length)
                {
                    if (_chars[child.Start + j] != pattern[i])
                        return results; // Mismatch
                    i++;
                    j++;
                }

                if (j == edgeLength)
                {
                    // Fully traversed this edge
                    depthFromRoot += edgeLength;
                    node = child;
                }
                else
                {
                    // Pattern ends in middle of edge 'child'
                    // We pass depthFromRoot (depth to parent), CollectLeaves will add child's edge
                    CollectLeaves(child, depthFromRoot, results);
                    return results;
                }
            }

            // Pattern matched exactly to a node - collect all leaves in subtree
            // depthFromRoot includes all edges up to and including node's edge
            // CollectLeaves expects depth EXCLUDING node's own edge
            CollectLeaves(node, depthFromRoot - LengthOf(node), results);
            return results;
        }

        /// <summary>
        /// Counts the number of occurrences of a pattern in the text.
        /// More efficient than FindAllOccurrences when you only need the count.
        ///
        /// Time complexity: O(m + k) where m is pattern length and k is number of occurrences.
        /// Space complexity: O(h) where h is tree height (for traversal stack).
        /// </summary>
        /// <param name="pattern">The pattern to count.</param>
        /// <returns>Number of times the pattern occurs in the text.</returns>
        /// <exception cref="ArgumentNullException">If pattern is null.</exception>
        public int CountOccurrences(string pattern)
        {
            ArgumentNullException.ThrowIfNull(pattern);

            return CountOccurrences(pattern.AsSpan());
        }

        /// <summary>
        /// Counts the number of occurrences of a pattern in the text.
        /// Zero-allocation overload for performance-critical scenarios.
        /// </summary>
        /// <param name="pattern">The pattern to count.</param>
        /// <returns>Number of times the pattern occurs in the text.</returns>
        /// <exception cref="ArgumentException">If pattern contains the null character '\0'.</exception>
        public int CountOccurrences(ReadOnlySpan<char> pattern)
        {
            if (pattern.Contains(TERMINATOR))
                throw new ArgumentException(
                    $"Pattern cannot contain the null character '\\0' as it is used as internal terminator.",
                    nameof(pattern));

            if (pattern.Length == 0)
                return _chars.Length - 1; // Empty pattern matches at every position (excluding terminator)

            // Navigate to the node/edge representing the pattern
            var node = _root;
            int i = 0;

            while (i < pattern.Length)
            {
                if (!node.TryGetChild(pattern[i], out var child) || child is null)
                    return 0; // Pattern not found

                int edgeLength = LengthOf(child);
                int j = 0;

                while (j < edgeLength && i < pattern.Length)
                {
                    if (_chars[child.Start + j] != pattern[i])
                        return 0; // Mismatch
                    i++;
                    j++;
                }

                if (j == edgeLength)
                {
                    node = child;
                }
                else
                {
                    // Pattern ends in middle of edge - return leaf count of this subtree
                    return child.LeafCount;
                }
            }

            // Pattern matched exactly to a node - count all leaves in subtree
            return node.LeafCount;
        }

        /// <summary>
        /// Finds the longest substring that appears at least twice in the text.
        ///
        /// Algorithm: The longest repeated substring corresponds to the deepest
        /// internal node in the suffix tree. Internal nodes represent branching points
        /// where multiple suffixes share a common prefix.
        ///
        /// Optimization: The deepest internal node is tracked during construction,
        /// making this operation O(depth) instead of O(n).
        ///
        /// Time complexity: O(depth) where depth is the length of the result.
        /// </summary>
        /// <returns>
        /// The longest repeated substring, or empty string if no repetition exists
        /// (i.e., all characters are unique).
        /// </returns>
        public string LongestRepeatedSubstring()
        {
            if (_chars == null || _chars.Length <= 1)
                return string.Empty;

            // Use cached deepest internal node (tracked during construction)
            if (_deepestInternalNode == null || _maxInternalDepth == 0)
                return string.Empty;

            // Reconstruct the substring from root to deepestNode
            return ReconstructPath(_deepestInternalNode, _maxInternalDepth);
        }

        /// <summary>
        /// Returns all suffixes of the original string in sorted order.
        ///
        /// Implementation uses iterative DFS with explicit stack to avoid StackOverflow
        /// on deep trees (e.g., strings with many repeated characters).
        ///
        /// Algorithm:
        /// 1. Traverse tree depth-first, visiting children in sorted character order
        /// 2. Track current path in StringBuilder (edge labels from root to current position)
        /// 3. When reaching a leaf, the path represents a complete suffix
        /// 4. Backtrack by removing characters when returning from child nodes
        ///
        /// Memory: O(n²) for result (n suffixes of average length n/2).
        /// For large strings, prefer EnumerateSuffixes() for lazy evaluation.
        /// </summary>
        /// <returns>All suffixes sorted lexicographically.</returns>
        public IReadOnlyList<string> GetAllSuffixes()
        {
            if (_chars == null || _chars.Length <= 1)
                return Array.Empty<string>();

            return EnumerateSuffixesCore().ToList();
        }

        /// <summary>
        /// Enumerates all suffixes of the original string in sorted order lazily.
        /// Use this for large strings to avoid O(n²) memory allocation from GetAllSuffixes.
        ///
        /// Each suffix is yielded one at a time, allowing early termination
        /// and streaming processing without loading all suffixes into memory.
        /// </summary>
        /// <returns>Lazy enumerable of suffixes sorted lexicographically.</returns>
        public IEnumerable<string> EnumerateSuffixes()
        {
            if (_chars == null || _chars.Length <= 1)
                return Enumerable.Empty<string>();

            return EnumerateSuffixesCore();
        }

        /// <summary>
        /// Core implementation for enumerating suffixes via iterative DFS.
        /// </summary>
        private IEnumerable<string> EnumerateSuffixesCore()
        {
            // Iterative DFS traversal yielding suffixes (paths to leaves)
            // Stack stores: (node, childIndex, sortedKeys) where childIndex tracks which children we've visited
            var stack = new Stack<(SuffixTreeNode Node, int ChildIndex, List<char> SortedKeys)>();
            var path = new StringBuilder(_chars!.Length); // Capacity = max possible path length
            var keyBuffer = new List<char>(8); // Reusable buffer for keys

            // Start with root's sorted children
            _root.GetKeys(keyBuffer);
            var rootKeys = new List<char>(keyBuffer);
            rootKeys.Sort();
            stack.Push((_root, 0, rootKeys));

            while (stack.Count > 0)
            {
                var (node, childIndex, sortedKeys) = stack.Pop();

                if (childIndex < sortedKeys.Count)
                {
                    // More children to process - push back with next index
                    stack.Push((node, childIndex + 1, sortedKeys));

                    var childKey = sortedKeys[childIndex];
                    if (!node.TryGetChild(childKey, out var child) || child is null)
                        continue;

                    // Add child's edge label to path
                    int edgeLen = LengthOf(child);
                    int charsAdded = 0;
                    for (int i = 0; i < edgeLen; i++)
                    {
                        char c = _chars[child.Start + i];
                        if (c == TERMINATOR) break;
                        path.Append(c);
                        charsAdded++;
                    }

                    if (child.IsLeaf)
                    {
                        // Yield suffix
                        if (path.Length > 0)
                            yield return path.ToString();
                        // Backtrack
                        path.Length -= charsAdded;
                    }
                    else
                    {
                        // Push child for processing with its sorted children
                        child.GetKeys(keyBuffer);
                        var childKeys = new List<char>(keyBuffer);
                        childKeys.Sort();
                        stack.Push((child, 0, childKeys));
                    }
                }
                else
                {
                    // All children processed - backtrack (remove this node's edge from path)
                    if (node != _root)
                    {
                        int edgeLen = LengthOf(node);
                        int charsToRemove = 0;
                        for (int i = 0; i < edgeLen; i++)
                        {
                            if (_chars[node.Start + i] == TERMINATOR) break;
                            charsToRemove++;
                        }
                        if (charsToRemove > 0 && path.Length >= charsToRemove)
                            path.Length -= charsToRemove;
                    }
                }
            }
        }

        /// <summary>
        /// Finds the longest common substring between this tree's text and another string.
        ///
        /// Algorithm: Uses suffix links for O(m) traversal instead of O(m²) naive approach.
        /// Maintains position in tree while scanning 'other', using suffix links to backtrack
        /// efficiently when mismatches occur.
        ///
        /// Time complexity: O(m) where m is the length of 'other'.
        /// </summary>
        /// <param name="other">The string to compare against.</param>
        /// <returns>The longest common substring, or empty string if none exists.</returns>
        /// <exception cref="ArgumentNullException">If other is null.</exception>
        /// <exception cref="ArgumentException">If other contains the null character '\0'.</exception>
        public string LongestCommonSubstring(string other)
        {
            ArgumentNullException.ThrowIfNull(other);

            return LongestCommonSubstring(other.AsSpan());
        }

        /// <summary>
        /// Finds the longest common substring between this tree's text and another character span.
        /// Zero-allocation overload for performance-critical scenarios.
        ///
        /// Algorithm: O(m) using suffix links.
        /// - Walk through 'other' character by character
        /// - Maintain current position in tree (node + offset on edge)
        /// - On mismatch, use suffix link to efficiently move to next shorter suffix
        /// - Track maximum match length throughout
        ///
        /// Time complexity: O(m) where m is the length of 'other'.
        /// </summary>
        /// <param name="other">The character span to compare against.</param>
        /// <returns>The longest common substring, or empty string if none exists.</returns>
        /// <exception cref="ArgumentException">If other contains the null character '\0'.</exception>
        public string LongestCommonSubstring(ReadOnlySpan<char> other)
        {
            if (other.Contains(TERMINATOR))
                throw new ArgumentException(
                    $"Input string cannot contain the null character '\\0' as it is used as internal terminator.",
                    nameof(other));

            if (other.Length == 0 || _chars.Length <= 1)
                return string.Empty;

            int maxLen = 0;
            int maxEndPos = -1; // End position in 'other' (exclusive)

            // Current position in tree
            var currentNode = _root;
            int edgeOffset = 0; // Position within current edge (0 = at node)
            SuffixTreeNode? currentEdge = null; // The edge we're currently on (null if at node)
            int currentMatchLen = 0; // Current match length

            int i = 0;
            while (i < other.Length)
            {
                char c = other[i];
                bool matched = false;

                if (currentEdge != null)
                {
                    // We're in the middle of an edge - check if next char matches
                    int edgeLen = LengthOf(currentEdge);
                    if (edgeOffset < edgeLen)
                    {
                        char edgeChar = _chars[currentEdge.Start + edgeOffset];
                        if (edgeChar != TERMINATOR && edgeChar == c)
                        {
                            // Match! Continue along edge
                            edgeOffset++;
                            currentMatchLen++;
                            matched = true;

                            // If we've consumed the entire edge, move to the child node
                            if (edgeOffset >= edgeLen)
                            {
                                currentNode = currentEdge;
                                currentEdge = null;
                                edgeOffset = 0;
                            }
                        }
                    }
                }
                else
                {
                    // We're at a node - look for child edge starting with c
                    if (currentNode.TryGetChild(c, out var child))
                    {
                        // Found edge - start matching
                        currentEdge = child;
                        edgeOffset = 1; // Already matched first char
                        currentMatchLen++;
                        matched = true;

                        // Check if edge is just one character
                        int edgeLen = LengthOf(child!);
                        if (edgeOffset >= edgeLen)
                        {
                            currentNode = child!;
                            currentEdge = null;
                            edgeOffset = 0;
                        }
                    }
                }

                if (matched)
                {
                    // Update max if this is the longest match so far
                    if (currentMatchLen > maxLen)
                    {
                        maxLen = currentMatchLen;
                        maxEndPos = i + 1; // End position (exclusive)
                    }
                    i++; // Move to next character in 'other'
                }
                else
                {
                    // Mismatch - use suffix link to try shorter suffix
                    if (currentMatchLen == 0)
                    {
                        // No match at all - just move to next character
                        i++;
                    }
                    else
                    {
                        // Use suffix link to efficiently backtrack
                        // We need to find where to continue from after removing first character of current match
                        if (currentEdge != null)
                        {
                            // In middle of edge - go back to parent and use its suffix link
                            currentNode = currentEdge.Parent ?? _root;
                        }

                        // Follow suffix link
                        if (currentNode != _root && currentNode.SuffixLink != null)
                        {
                            currentNode = currentNode.SuffixLink;
                        }
                        else
                        {
                            currentNode = _root;
                        }

                        currentMatchLen--;
                        currentEdge = null;
                        edgeOffset = 0;

                        // After suffix link, we need to "rescan" to find our position
                        // This is needed because suffix link destination might not have same edge structure
                        if (currentMatchLen > 0)
                        {
                            // Rescan from i - currentMatchLen to find position
                            int scanStart = i - currentMatchLen;
                            currentNode = _root;
                            currentEdge = null;
                            edgeOffset = 0;
                            int newMatchLen = 0;

                            for (int j = scanStart; j < i; j++)
                            {
                                char sc = other[j];
                                if (currentEdge != null)
                                {
                                    int edgeLen = LengthOf(currentEdge);
                                    if (edgeOffset < edgeLen)
                                    {
                                        char edgeChar = _chars[currentEdge.Start + edgeOffset];
                                        if (edgeChar != TERMINATOR && edgeChar == sc)
                                        {
                                            edgeOffset++;
                                            newMatchLen++;
                                            if (edgeOffset >= edgeLen)
                                            {
                                                currentNode = currentEdge;
                                                currentEdge = null;
                                                edgeOffset = 0;
                                            }
                                            continue;
                                        }
                                    }
                                    break; // Mismatch during rescan - shouldn't happen
                                }
                                else
                                {
                                    if (currentNode.TryGetChild(sc, out var child))
                                    {
                                        currentEdge = child;
                                        edgeOffset = 1;
                                        newMatchLen++;
                                        int edgeLen = LengthOf(child!);
                                        if (edgeOffset >= edgeLen)
                                        {
                                            currentNode = child!;
                                            currentEdge = null;
                                            edgeOffset = 0;
                                        }
                                        continue;
                                    }
                                    break; // No edge - shouldn't happen if tree is correct
                                }
                            }
                            currentMatchLen = newMatchLen;
                        }
                        // Don't increment i - try matching current char c again with new position
                    }
                }
            }

            if (maxLen == 0)
                return string.Empty;

            return other.Slice(maxEndPos - maxLen, maxLen).ToString();
        }

        /// <summary>
        /// Finds the longest common substring with position information.
        ///
        /// Time complexity: O(n * m) where n is tree text length, m is other length.
        /// </summary>
        /// <param name="other">The string to compare against.</param>
        /// <returns>
        /// A tuple containing: the substring, position in tree's text, position in other.
        /// Returns (empty string, -1, -1) if no common substring exists.
        /// </returns>
        /// <exception cref="ArgumentNullException">If other is null.</exception>
        /// <exception cref="ArgumentException">If other contains the null character '\0'.</exception>
        public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other)
        {
            ArgumentNullException.ThrowIfNull(other);

            if (other.Contains(TERMINATOR))
                throw new ArgumentException(
                    $"Input string cannot contain the null character '\\0' as it is used as internal terminator.",
                    nameof(other));

            if (other.Length == 0 || _chars.Length <= 1)
                return (string.Empty, -1, -1);

            int maxLen = 0;
            int maxStartInOther = 0;
            int maxStartInText = -1;

            for (int start = 0; start < other.Length; start++)
            {
                var (matchLen, startInText) = MatchFromRootWithPosition(other.AsSpan(), start);
                if (matchLen > maxLen)
                {
                    maxLen = matchLen;
                    maxStartInOther = start;
                    maxStartInText = startInText;
                }
            }

            if (maxLen == 0)
                return (string.Empty, -1, -1);

            string substring = other.Substring(maxStartInOther, maxLen);
            return (substring, maxStartInText, maxStartInOther);
        }

        /// <summary>
        /// Finds all positions where the longest common substring occurs.
        /// Time complexity: O(n * m) where n is tree text length, m is other length.
        /// </summary>
        /// <param name="other">The string to compare against.</param>
        /// <returns>
        /// A tuple containing: the substring, all positions in tree's text, all positions in other.
        /// Returns (empty string, empty list, empty list) if no common substring exists.
        /// </returns>
        /// <exception cref="ArgumentNullException">If other is null.</exception>
        /// <exception cref="ArgumentException">If other contains the null character '\0'.</exception>
        public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
        {
            ArgumentNullException.ThrowIfNull(other);

            if (other.Contains(TERMINATOR))
                throw new ArgumentException(
                    $"Input string cannot contain the null character '\\0' as it is used as internal terminator.",
                    nameof(other));

            if (other.Length == 0 || _chars.Length <= 1)
                return (string.Empty, [], []);

            int maxLen = 0;
            var matchesInOther = new List<(int StartInOther, int StartInText)>();

            for (int start = 0; start < other.Length; start++)
            {
                var (matchLen, startInText) = MatchFromRootWithPosition(other.AsSpan(), start);
                if (matchLen > maxLen)
                {
                    maxLen = matchLen;
                    matchesInOther.Clear();
                    matchesInOther.Add((start, startInText));
                }
                else if (matchLen == maxLen && maxLen > 0)
                {
                    matchesInOther.Add((start, startInText));
                }
            }

            if (maxLen == 0)
                return (string.Empty, Array.Empty<int>(), Array.Empty<int>());

            string substring = other.Substring(matchesInOther[0].StartInOther, maxLen);

            var positionsInText = new List<int>();
            var positionsInOther = new List<int>();

            foreach (var (startInOther, startInText) in matchesInOther)
            {
                positionsInOther.Add(startInOther);
                positionsInText.Add(startInText);
            }

            return (substring, positionsInText, positionsInOther);
        }

        /// <summary>
        /// Matches characters from 'other' starting at 'start' against the tree from root.
        /// Returns the number of characters matched and the starting position in tree's text.
        /// </summary>
        /// <param name="other">The character span to match against.</param>
        /// <param name="start">Starting position in 'other'.</param>
        /// <returns>Tuple of (match length, start position in text). Position is -1 if no match.</returns>
        private (int MatchLength, int StartPositionInText) MatchFromRootWithPosition(ReadOnlySpan<char> other, int start)
        {
            var node = _root;
            int i = start;
            int matched = 0;
            int matchStartInText = -1;

            while (i < other.Length)
            {
                char c = other[i];

                if (!node.TryGetChild(c, out var child) || child is null)
                    break; // No edge starting with c

                // Track the start position in tree's text (first character matched)
                if (matched == 0)
                    matchStartInText = child.Start;

                // Match along this edge
                int edgeLen = LengthOf(child);
                int j = 0;

                while (j < edgeLen && i < other.Length)
                {
                    char edgeChar = _chars[child.Start + j];
                    if (edgeChar == TERMINATOR || edgeChar != other[i])
                        return (matched, matchStartInText);

                    matched++;
                    i++;
                    j++;
                }

                if (j < edgeLen)
                    break; // Didn't finish edge (hit terminator or end of other)

                // Move to child node
                node = child;
            }

            return (matched, matchStartInText);
        }

        /// <summary>
        /// Reconstructs the path label from root to a given node.
        /// Uses Parent references for O(depth) traversal instead of O(n) DFS search.
        /// </summary>
        private string ReconstructPath(SuffixTreeNode targetNode, int pathLength)
        {
            // Walk from targetNode up to root, collecting nodes
            var pathFromTarget = new List<SuffixTreeNode>();
            var current = targetNode;

            while (current != null && current != _root)
            {
                pathFromTarget.Add(current);
                current = current.Parent;
            }

            // Reverse to get root-to-target order
            var sb = new StringBuilder(pathLength);
            for (int i = pathFromTarget.Count - 1; i >= 0; i--)
            {
                var node = pathFromTarget[i];
                int len = LengthOf(node);
                for (int j = 0; j < len; j++)
                {
                    char c = _chars[node.Start + j];
                    if (c == TERMINATOR) break; // Don't include terminator
                    sb.Append(c);
                }
            }

            return sb.ToString();
        }

        /// <summary>
        /// Collects all leaf positions starting from the given node using iterative traversal.
        /// Avoids stack overflow for deep trees.
        /// </summary>
        /// <param name="node">Node to start collecting from.</param>
        /// <param name="depth">Depth from root to this node (sum of edge lengths on path, NOT including this node's edge).</param>
        /// <param name="results">List to collect results into.</param>
        private void CollectLeaves(SuffixTreeNode node, int depth, List<int> results)
        {
            // Stack stores (node, depthBeforeNode) pairs
            var stack = new Stack<(SuffixTreeNode Node, int Depth)>();
            stack.Push((node, depth));

            while (stack.Count > 0)
            {
                var (current, currentDepthBefore) = stack.Pop();
                int currentDepth = currentDepthBefore + LengthOf(current);

                if (current.IsLeaf)
                {
                    // This leaf represents a suffix of length = currentDepth
                    int suffixLength = currentDepth;
                    int startPosition = _chars.Length - suffixLength;
                    results.Add(startPosition);
                }
                else
                {
                    // Push children with current depth as their starting depth
                    foreach (var child in current.GetChildren())
                    {
                        stack.Push((child, currentDepth));
                    }
                }
            }
        }

        /// <summary>
        /// Returns a brief string representation of the tree.
        /// </summary>
        public override string ToString()
        {
            if (_chars.Length == 0)
                return "SuffixTree (empty)";

            // Exclude terminator - it's always the last character
            var content = new string(_chars, 0, _chars.Length - 1);
            if (content.Length > MAX_TOSTRING_CONTENT_LENGTH)
                content = string.Concat(content.AsSpan(0, MAX_TOSTRING_CONTENT_LENGTH - 3), "...");

            return string.Create(System.Globalization.CultureInfo.InvariantCulture,
                $"SuffixTree (length: {_chars.Length}, content: \"{content}\")");
        }

        /// <summary>
        /// Creates a detailed string representation of the tree structure.
        /// Useful for debugging and visualization.
        /// Uses iterative traversal to avoid StackOverflow on deep trees.
        /// </summary>
        public string PrintTree()
        {
            // Estimate capacity: ~50 chars per node, string of n chars has at most 2n-1 nodes
            int estimatedNodes = Math.Max(1, _chars.Length * 2);
            var sb = new StringBuilder(Math.Max(256, estimatedNodes * 50));
            var ci = System.Globalization.CultureInfo.InvariantCulture;
            sb.Append(ci, $"Content length: {_chars.Length}").AppendLine();
            sb.AppendLine();

            // Iterative DFS: stack stores (node, depth, childIndex, sortedChildren)
            var stack = new Stack<(SuffixTreeNode Node, int Depth, int ChildIndex, List<SuffixTreeNode> SortedChildren)>();
            var keyBuffer = new List<char>(8);

            // Print root first
            var rootLabel = LabelOf(_root);
            sb.Append(ci, $"0:{rootLabel}").AppendLine();

            var rootChildren = GetSortedChildren(_root, keyBuffer);
            if (rootChildren.Count > 0)
                stack.Push((_root, 0, 0, rootChildren));

            while (stack.Count > 0)
            {
                var (node, depth, childIndex, sortedChildren) = stack.Pop();

                if (childIndex < sortedChildren.Count)
                {
                    // More children to process
                    stack.Push((node, depth, childIndex + 1, sortedChildren));

                    var child = sortedChildren[childIndex];
                    int childDepth = depth + 1;

                    // Print this child
                    var nodeLabel = LabelOf(child);
                    var leafMark = child.IsLeaf ? "..." : "";
                    var linkMark = child.SuffixLink != null && child.SuffixLink != _root && !child.IsLeaf
                        ? $" -> {FirstCharOf(child.SuffixLink)}"
                        : "";
                    sb.Append(' ', childDepth);
                    sb.Append(ci, $"{childDepth}:{nodeLabel}{leafMark}{linkMark}").AppendLine();

                    // If child has children, push it for processing
                    if (!child.IsLeaf && child.ChildCount > 0)
                    {
                        var grandChildren = GetSortedChildren(child, keyBuffer);
                        stack.Push((child, childDepth, 0, grandChildren));
                    }
                }
            }

            return sb.ToString();
        }

        /// <summary>
        /// Returns children of a node sorted by edge character.
        /// </summary>
        private static List<SuffixTreeNode> GetSortedChildren(SuffixTreeNode node, List<char> keyBuffer)
        {
            node.GetKeys(keyBuffer);
            keyBuffer.Sort();
            var result = new List<SuffixTreeNode>(keyBuffer.Count);
            foreach (var key in keyBuffer)
            {
                node.TryGetChild(key, out var child);
                result.Add(child!);
            }
            return result;
        }
    }
}
