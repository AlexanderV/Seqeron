using System;
using System.Collections.Generic;
using System.Linq;
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
        private SuffixTreeNode _lastCreatedInternalNode;

        /// <summary>The string content stored as a character array.</summary>
        private char[] _chars;

        // ============================================================
        // Active Point - tracks where we are in the tree during construction
        // ============================================================
        // When Rule 3 (showstopper) fires, we don't insert but remember WHERE
        // the suffix exists implicitly. The active point tracks this location.
        // ============================================================

        /// <summary>The node we're currently at (or descended from).</summary>
        private SuffixTreeNode _activeNode;

        /// <summary>Index in _chars of the first character of the active edge (-1 = none).</summary>
        private int _activeEdgeIndex = -1;

        /// <summary>How many characters along the active edge we've matched.</summary>
        private int _activeLength = 0;

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
            if (value == null)
                throw new ArgumentNullException(nameof(value));

            if (value.Contains(TERMINATOR))
                throw new ArgumentException(
                    $"Input string cannot contain the null character '\\0' as it is used as internal terminator.",
                    nameof(value));

            if (value.Length == 0)
            {
                _chars = Array.Empty<char>();
                return;
            }

            // Pre-allocate array for string + terminator
            _chars = new char[value.Length + 1];

            // Process each character - Ukkonen's online construction
            foreach (var c in value)
                ExtendTree(c);

            // Add terminator to convert implicit suffixes to explicit leaves
            ExtendTree(TERMINATOR);

            // Clear construction state (not needed anymore, helps GC)
            _lastCreatedInternalNode = null;
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

                if (!_activeNode.Children.TryGetValue(activeEdgeChar, out SuffixTreeNode activeEdge))
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
                    var leaf = new SuffixTreeNode { Start = _position, End = SuffixTreeNode.BOUNDLESS };
                    _activeNode.Children[c] = leaf;

                    // Set suffix link for previously created internal node
                    AddSuffixLink(_activeNode);
                }
                else
                {
                    // Edge exists - need to check if we should walk down first
                    if (WalkDownIfNeeded(activeEdge))
                        continue; // Active point moved, restart this iteration

                    // Check if next character on edge matches
                    if (_chars[activeEdge.Start + _activeLength] == c)
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
                    var splitNode = new SuffixTreeNode
                    {
                        Start = activeEdge.Start,
                        End = activeEdge.Start + _activeLength
                    };

                    // Replace old edge with new internal node
                    _activeNode.Children[activeEdgeChar] = splitNode;

                    // Create new leaf for current character
                    var newLeaf = new SuffixTreeNode { Start = _position, End = SuffixTreeNode.BOUNDLESS };
                    splitNode.Children[c] = newLeaf;

                    // Move old edge to be child of split node
                    activeEdge.Start += _activeLength;
                    splitNode.Children[_chars[activeEdge.Start]] = activeEdge;

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
            if (value == null)
                throw new ArgumentNullException(nameof(value));

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
                if (!node.Children.TryGetValue(value[i], out var child))
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
            if (pattern == null)
                throw new ArgumentNullException(nameof(pattern));

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
                if (!node.Children.TryGetValue(pattern[i], out var child))
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
            // depthFromRoot already includes path to 'node', so we pass depth BEFORE node's edge
            // But wait - depthFromRoot was updated to include node's edge. Let's trace:
            // After loop: depthFromRoot = sum of all edges INCLUDING node's edge
            // In CollectLeaves we'll add node's edge again - that's wrong!
            // 
            // Actually, 'node' here IS the node we landed on after following edges.
            // So depthFromRoot includes all edges UP TO AND INCLUDING the edge to 'node'.
            // But CollectLeaves expects depth EXCLUDING node's own edge.
            // So we need: CollectLeaves(node, depthFromRoot - LengthOf(node), results)
            // But that's getting complicated. Let me redesign.
            CollectLeaves(node, depthFromRoot - LengthOf(node), results);
            return results;
        }

        /// <summary>
        /// Counts the number of occurrences of a pattern in the text.
        /// More efficient than FindAllOccurrences when you only need the count.
        /// 
        /// Time complexity: O(m + k) where m is pattern length and k is number of occurrences.
        /// Space complexity: O(1) additional space (no list allocation).
        /// </summary>
        /// <param name="pattern">The pattern to count.</param>
        /// <returns>Number of times the pattern occurs in the text.</returns>
        /// <exception cref="ArgumentNullException">If pattern is null.</exception>
        public int CountOccurrences(string pattern)
        {
            if (pattern == null)
                throw new ArgumentNullException(nameof(pattern));

            if (pattern.Length == 0)
                return _chars.Length - 1; // Empty pattern matches at every position (excluding terminator)

            // Navigate to the node/edge representing the pattern
            var node = _root;
            int i = 0;
            int depthFromRoot = 0;

            while (i < pattern.Length)
            {
                if (!node.Children.TryGetValue(pattern[i], out var child))
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
                    depthFromRoot += edgeLength;
                    node = child;
                }
                else
                {
                    // Pattern ends in middle of edge - count leaves from here
                    return CountLeaves(child);
                }
            }

            // Pattern matched exactly to a node - count all leaves in subtree
            return CountLeaves(node);
        }

        /// <summary>
        /// Counts all leaves in a subtree using iterative traversal.
        /// Avoids stack overflow for deep trees.
        /// </summary>
        private int CountLeaves(SuffixTreeNode node)
        {
            int count = 0;
            var stack = new Stack<SuffixTreeNode>();
            stack.Push(node);

            while (stack.Count > 0)
            {
                var current = stack.Pop();

                if (current.IsLeaf)
                {
                    count++;
                }
                else
                {
                    foreach (var child in current.Children.Values)
                    {
                        stack.Push(child);
                    }
                }
            }

            return count;
        }

        /// <summary>
        /// Finds the longest substring that appears at least twice in the text.
        /// 
        /// Algorithm: The longest repeated substring corresponds to the deepest 
        /// internal node in the suffix tree. Internal nodes represent branching points
        /// where multiple suffixes share a common prefix.
        /// 
        /// Time complexity: O(n) where n is the text length.
        /// </summary>
        /// <returns>
        /// The longest repeated substring, or empty string if no repetition exists
        /// (i.e., all characters are unique).
        /// </returns>
        public string LongestRepeatedSubstring()
        {
            if (_chars == null || _chars.Length <= 1)
                return string.Empty;

            int maxDepth = 0;
            SuffixTreeNode deepestNode = null;

            FindDeepestInternalNode(_root, 0, ref maxDepth, ref deepestNode);

            if (deepestNode == null)
                return string.Empty;

            // Reconstruct the substring from root to deepestNode
            return ReconstructPath(deepestNode, maxDepth);
        }

        /// <summary>
        /// Recursively finds the deepest internal node.
        /// </summary>
        private void FindDeepestInternalNode(SuffixTreeNode node, int depth, ref int maxDepth, ref SuffixTreeNode deepestNode)
        {
            int currentDepth = depth + LengthOf(node);

            // Only internal nodes (non-leaves) represent repeated substrings
            if (!node.IsLeaf && node.Children.Count > 0)
            {
                if (currentDepth > maxDepth)
                {
                    maxDepth = currentDepth;
                    deepestNode = node;
                }
            }

            foreach (var child in node.Children.Values)
            {
                FindDeepestInternalNode(child, currentDepth, ref maxDepth, ref deepestNode);
            }
        }

        /// <summary>
        /// Reconstructs the path label from root to a given node.
        /// </summary>
        private string ReconstructPath(SuffixTreeNode targetNode, int pathLength)
        {
            // We need to trace from root to targetNode
            // Since nodes don't store parent references, we'll do a DFS to find the path
            var path = new List<SuffixTreeNode>();
            if (FindPathToNode(_root, targetNode, path))
            {
                var sb = new StringBuilder(pathLength);
                foreach (var node in path)
                {
                    if (node == _root) continue;
                    int len = LengthOf(node);
                    for (int i = 0; i < len; i++)
                    {
                        char c = _chars[node.Start + i];
                        if (c == TERMINATOR) break; // Don't include terminator
                        sb.Append(c);
                    }
                }
                return sb.ToString();
            }
            return string.Empty;
        }

        /// <summary>
        /// Finds the path from source to target node using DFS.
        /// </summary>
        private bool FindPathToNode(SuffixTreeNode current, SuffixTreeNode target, List<SuffixTreeNode> path)
        {
            path.Add(current);

            if (current == target)
                return true;

            foreach (var child in current.Children.Values)
            {
                if (FindPathToNode(child, target, path))
                    return true;
            }

            path.RemoveAt(path.Count - 1);
            return false;
        }

        /// <summary>
        /// Recursively collects all leaf positions starting from the given node.
        /// </summary>
        /// <param name="node">Node to start collecting from.</param>
        /// <param name="depth">Depth from root to this node (sum of edge lengths on path, NOT including this node's edge).</param>
        /// <param name="results">List to collect results into.</param>
        private void CollectLeaves(SuffixTreeNode node, int depth, List<int> results)
        {
            // Add current node's edge length to get total depth
            int currentDepth = depth + LengthOf(node);

            if (node.IsLeaf)
            {
                // This leaf represents a suffix of length = currentDepth
                int suffixLength = currentDepth;
                int startPosition = _chars.Length - suffixLength;
                results.Add(startPosition);
                return;
            }

            // Recurse into children, passing currentDepth as their starting depth
            foreach (var child in node.Children.Values)
            {
                CollectLeaves(child, currentDepth, results);
            }
        }

        /// <summary>
        /// Returns a brief string representation of the tree.
        /// </summary>
        public override string ToString()
        {
            if (_chars == null || _chars.Length == 0)
                return "SuffixTree (empty)";

            var content = new string(_chars.Where(c => c != TERMINATOR).ToArray());
            if (content.Length > 50)
                content = content.Substring(0, 47) + "...";

            return $"SuffixTree (length: {_chars.Length}, content: \"{content}\")";
        }

        /// <summary>
        /// Creates a detailed string representation of the tree structure.
        /// Useful for debugging and visualization.
        /// </summary>
        public string PrintTree()
        {
            var sb = new StringBuilder();
            sb.AppendLine($"Content length: {_chars?.Length ?? 0}");
            sb.AppendLine();
            PrintNode(sb, _root, 0);
            return sb.ToString();
        }

        private void PrintNode(StringBuilder sb, SuffixTreeNode node, int depth)
        {
            var nodeLabel = LabelOf(node);
            var leafMark = node.IsLeaf ? "..." : "";
            var linkMark = node.SuffixLink != null && node.SuffixLink != _root && !node.IsLeaf
                ? $" -> {FirstCharOf(node.SuffixLink)}"
                : "";

            sb.AppendLine($"{new string(' ', depth)}{depth}:{nodeLabel}{leafMark}{linkMark}");

            foreach (var child in node.Children.OrderBy(kvp => kvp.Key))
            {
                PrintNode(sb, child.Value, depth + 1);
            }
        }
    }
}