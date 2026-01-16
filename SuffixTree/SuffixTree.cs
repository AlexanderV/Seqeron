using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SuffixTree
{
    /// <summary>
    /// Implementation of Ukkonen's online suffix tree construction algorithm.
    /// Time complexity: O(n) for construction, O(m) for substring search.
    /// </summary>
    public class SuffixTree
    {
        private const int BOUNDLESS = -1;
        private const char TERMINATOR = '\0';

        private int _remainder, _position = -1;
        private Node _root, _lastCreatedInternalNode;
        private List<char> _chars = new List<char>();

        // Active point - represents where we are in the tree during construction
        private Node _activeNode;
        private int _activeEdgeIndex = -1; // Index in _chars of the first char of active edge (-1 = no edge)
        private int _activeLength = 0;

        private class Node
        {
            public int Start { get; set; }
            public int End { get; set; }
            public Dictionary<char, Node> Children { get; } = new Dictionary<char, Node>();
            public Node SuffixLink { get; set; }
            public bool IsLeaf => End == BOUNDLESS;
        }

        public SuffixTree()
        {
            _root = new Node { Start = 0, End = 0 };
            _root.SuffixLink = _root; // Root links to itself
            _activeNode = _root;
        }

        /// <summary>
        /// Creates and returns a tree with the specified value.
        /// </summary>
        public static SuffixTree Build(string value)
        {
            var t = new SuffixTree();
            t.AddString(value);
            return t;
        }

        /// <summary>
        /// Extends the suffix tree with the specified value.
        /// A unique terminator character is automatically appended.
        /// </summary>
        public void AddString(string value)
        {
            if (value == null)
                throw new ArgumentNullException(nameof(value));

            if (value.Length == 0)
                return;

            foreach (var c in value)
                ExtendTree(c);

            // Add terminating character to ensure all suffixes end at leaves
            ExtendTree(TERMINATOR);

            // Reset state for potential next string
            _remainder = 0;
            _activeNode = _root;
            _activeEdgeIndex = -1;
            _activeLength = 0;
            _lastCreatedInternalNode = null;
        }

        private void ExtendTree(char c)
        {
            _chars.Add(c);
            _position++;
            _remainder++;
            _lastCreatedInternalNode = null;

            while (_remainder > 0)
            {
                // If activeLength is 0, we're at a node - set active edge to current char
                if (_activeLength == 0)
                    _activeEdgeIndex = _position;

                // Check if there's an edge starting with active edge char
                char activeEdgeChar = _chars[_activeEdgeIndex];

                if (!_activeNode.Children.TryGetValue(activeEdgeChar, out Node activeEdge))
                {
                    // No edge with this char - create new leaf directly from activeNode
                    var leaf = new Node { Start = _position, End = BOUNDLESS };
                    _activeNode.Children[c] = leaf;

                    // If there was a previously created internal node, link it to activeNode
                    AddSuffixLink(_activeNode);
                }
                else
                {
                    // Edge exists - check if we can walk down
                    if (WalkDownIfNeeded(activeEdge))
                        continue; // Active point changed, restart this phase

                    // Check if next char on edge matches
                    if (_chars[activeEdge.Start + _activeLength] == c)
                    {
                        // Character matches - increment active length and stop (Rule 3 extension)
                        // This is a "showstopper" - we don't create any new nodes, just extend implicitly
                        // Note: We still need to set suffix link if we created an internal node earlier
                        AddSuffixLink(_activeNode);
                        _activeLength++;
                        break;
                    }

                    // Character doesn't match - need to split the edge
                    var splitNode = new Node
                    {
                        Start = activeEdge.Start,
                        End = activeEdge.Start + _activeLength
                    };

                    _activeNode.Children[activeEdgeChar] = splitNode;

                    // Create new leaf for current char
                    var newLeaf = new Node { Start = _position, End = BOUNDLESS };
                    splitNode.Children[c] = newLeaf;

                    // Update old edge and make it child of split node
                    activeEdge.Start += _activeLength;
                    splitNode.Children[_chars[activeEdge.Start]] = activeEdge;

                    // This is a new internal node - set up suffix link chain
                    SetLastCreatedInternalNode(splitNode);
                }

                _remainder--;

                if (_activeNode == _root && _activeLength > 0)
                {
                    // Rule 1: If active node is root, decrement active length
                    // and set active edge to next suffix start
                    _activeLength--;
                    _activeEdgeIndex = _position - _remainder + 1;
                }
                else if (_activeNode != _root)
                {
                    // Rule 3: Follow suffix link
                    _activeNode = _activeNode.SuffixLink ?? _root;
                }
            }
        }

        /// <summary>
        /// Walk down the tree if activeLength >= edge length.
        /// Returns true if we moved (caller should restart), false otherwise.
        /// </summary>
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

        private void AddSuffixLink(Node node)
        {
            // Link previously created internal node to this node
            if (_lastCreatedInternalNode != null)
            {
                _lastCreatedInternalNode.SuffixLink = node;
                _lastCreatedInternalNode = null;
            }
        }

        private void SetLastCreatedInternalNode(Node node)
        {
            // Track this internal node for suffix linking in next iteration
            if (_lastCreatedInternalNode != null)
                _lastCreatedInternalNode.SuffixLink = node;

            _lastCreatedInternalNode = node;
        }

        private int LengthOf(Node edge)
            => (edge.End == BOUNDLESS ? _position + 1 : edge.End) - edge.Start;

        private char FirstCharOf(Node edge)
            => _chars[edge.Start];

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
        /// Checks if the specified value is a substring of the tree content.
        /// Executes with O(m) time complexity, where m is the length of the value.
        /// </summary>
        public bool Contains(string value)
        {
            if (value == null)
                throw new ArgumentNullException(nameof(value));

            if (value.Length == 0)
                return true;

            var node = _root;
            int i = 0;

            while (i < value.Length)
            {
                if (!node.Children.TryGetValue(value[i], out var child))
                    return false;

                int edgeLength = LengthOf(child);
                int j = 0;

                // Match characters along this edge
                while (j < edgeLength && i < value.Length)
                {
                    if (_chars[child.Start + j] != value[i])
                        return false;
                    i++;
                    j++;
                }

                // If we consumed the entire edge, move to next node
                if (j == edgeLength)
                    node = child;
                // Otherwise we're done (matched in middle of edge)
            }

            return true;
        }

        /// <summary>
        /// Returns a string representation of the tree content.
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
        /// Useful for debugging.
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