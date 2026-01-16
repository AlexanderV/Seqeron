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
    public class SuffixTree
    {
        /// <summary>
        /// Sentinel value indicating an open-ended (growing) edge.
        /// Leaf edges use this to automatically extend as new characters are added.
        /// </summary>
        private const int BOUNDLESS = -1;

        /// <summary>
        /// Unique terminator character appended to ensure all suffixes are explicit.
        /// Without this, some suffixes may remain implicit (e.g., "a" in "aa").
        /// </summary>
        private const char TERMINATOR = '\0';

        /// <summary>Number of suffixes that still need explicit insertion.</summary>
        private int _remainder;

        /// <summary>Current position in _chars (0-based index of last added character).</summary>
        private int _position = -1;

        private Node _root;

        /// <summary>
        /// Last internal node created in current phase.
        /// Used to set up suffix link chain: when we create a new internal node,
        /// we link the previous one to it.
        /// </summary>
        private Node _lastCreatedInternalNode;

        /// <summary>The string content stored as a list for efficient appending.</summary>
        private List<char> _chars = new List<char>();

        // ============================================================
        // Active Point - tracks where we are in the tree during construction
        // ============================================================
        // When Rule 3 (showstopper) fires, we don't insert but remember WHERE
        // the suffix exists implicitly. The active point tracks this location.
        // ============================================================

        /// <summary>The node we're currently at (or descended from).</summary>
        private Node _activeNode;

        /// <summary>Index in _chars of the first character of the active edge (-1 = none).</summary>
        private int _activeEdgeIndex = -1;

        /// <summary>How many characters along the active edge we've matched.</summary>
        private int _activeLength = 0;

        /// <summary>
        /// Internal node representation.
        /// Each edge is implicitly stored as (Start, End) indices into _chars.
        /// </summary>
        private class Node
        {
            /// <summary>Start index of edge label in _chars.</summary>
            public int Start { get; set; }

            /// <summary>
            /// End index of edge label (exclusive). 
            /// BOUNDLESS (-1) means this is a leaf edge that grows automatically.
            /// </summary>
            public int End { get; set; }

            /// <summary>Children edges, keyed by first character.</summary>
            public Dictionary<char, Node> Children { get; } = new Dictionary<char, Node>();

            /// <summary>
            /// Suffix link: connects node for "xα" to node for "α".
            /// Used for O(1) jumps between suffixes during construction.
            /// </summary>
            public Node SuffixLink { get; set; }

            /// <summary>True if this is a leaf node (edge grows with string).</summary>
            public bool IsLeaf => End == BOUNDLESS;
        }

        public SuffixTree()
        {
            _root = new Node { Start = 0, End = 0 };
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
            t.AddString(value);
            return t;
        }

        /// <summary>
        /// Extends the suffix tree with the specified string.
        /// A unique terminator character ($) is automatically appended to ensure
        /// all suffixes are explicit (end at leaf nodes).
        /// </summary>
        /// <param name="value">The string to add.</param>
        /// <exception cref="ArgumentNullException">If value is null.</exception>
        public void AddString(string value)
        {
            if (value == null)
                throw new ArgumentNullException(nameof(value));

            if (value.Length == 0)
                return;

            // Process each character - Ukkonen's online construction
            foreach (var c in value)
                ExtendTree(c);

            // Add terminator to convert implicit suffixes to explicit leaves
            ExtendTree(TERMINATOR);

            // Reset state for potential next string
            _remainder = 0;
            _activeNode = _root;
            _activeEdgeIndex = -1;
            _activeLength = 0;
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
            _chars.Add(c);
            _position++;
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

                if (!_activeNode.Children.TryGetValue(activeEdgeChar, out Node activeEdge))
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
                    var leaf = new Node { Start = _position, End = BOUNDLESS };
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
                    var splitNode = new Node
                    {
                        Start = activeEdge.Start,
                        End = activeEdge.Start + _activeLength
                    };

                    // Replace old edge with new internal node
                    _activeNode.Children[activeEdgeChar] = splitNode;

                    // Create new leaf for current character
                    var newLeaf = new Node { Start = _position, End = BOUNDLESS };
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
        private bool WalkDownIfNeeded(Node edge)
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
        private void AddSuffixLink(Node node)
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
        private void SetLastCreatedInternalNode(Node node)
        {
            if (_lastCreatedInternalNode != null)
                _lastCreatedInternalNode.SuffixLink = node;

            _lastCreatedInternalNode = node;
        }

        /// <summary>
        /// Calculates the length of an edge.
        /// For leaves (End == BOUNDLESS), the edge extends to current position.
        /// </summary>
        private int LengthOf(Node edge)
            => (edge.End == BOUNDLESS ? _position + 1 : edge.End) - edge.Start;

        /// <summary>Returns the first character of an edge's label.</summary>
        private char FirstCharOf(Node edge)
            => _chars[edge.Start];

        /// <summary>Returns the full label of an edge (for debugging).</summary>
        private string LabelOf(Node edge)
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

            if (value.Length == 0)
                return true; // Empty string is always a substring

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
        /// Returns a brief string representation of the tree.
        /// </summary>
        public override string ToString()
        {
            if (_chars.Count == 0)
                return "SuffixTree (empty)";

            var content = new string(_chars.Where(c => c != TERMINATOR).ToArray());
            if (content.Length > 50)
                content = content.Substring(0, 47) + "...";

            return $"SuffixTree (length: {_chars.Count}, content: \"{content}\")";
        }

        /// <summary>
        /// Creates a detailed string representation of the tree structure.
        /// Useful for debugging and visualization.
        /// </summary>
        public string PrintTree()
        {
            var sb = new StringBuilder();
            sb.AppendLine($"Content length: {_chars.Count}");
            sb.AppendLine();
            PrintNode(sb, _root, 0);
            return sb.ToString();
        }

        private void PrintNode(StringBuilder sb, Node node, int depth)
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