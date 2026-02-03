using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace SuffixTree
{
    /// <summary>
    /// Internal node representation for suffix tree.
    /// Each edge is implicitly stored as (Start, End) indices into the character array.
    ///
    /// Memory optimization: Uses hybrid children storage.
    /// - Inline array for ≤4 children (most common case) - avoids Dictionary overhead
    /// - Dictionary for >4 children - efficient lookup for large alphabet
    /// </summary>
    internal sealed class SuffixTreeNode
    {
        /// <summary>
        /// Sentinel value indicating an open-ended (growing) edge.
        /// Leaf edges use this to automatically extend as new characters are added.
        /// </summary>
        internal const int BOUNDLESS = -1;

        /// <summary>
        /// Threshold for switching from inline array to Dictionary.
        /// Most suffix tree nodes have 1-3 children, so 4 is a good cutoff.
        /// </summary>
        private const int INLINE_THRESHOLD = 4;

        /// <summary>Start index of edge label in character array.</summary>
        public int Start { get; set; }

        /// <summary>
        /// End index of edge label (exclusive).
        /// BOUNDLESS (-1) means this is a leaf edge that grows automatically.
        /// </summary>
        public int End { get; set; }

        // Hybrid children storage: inline array for small, Dictionary for large
        private byte _inlineCount;
        private char _key0, _key1, _key2, _key3;
        private SuffixTreeNode? _child0, _child1, _child2, _child3;
        private Dictionary<char, SuffixTreeNode>? _overflow;
        private SuffixTreeNode? _terminatorChild;

        /// <summary>
        /// Returns true if this node has any children.
        /// </summary>
        public bool HasChildren => _inlineCount > 0 || (_overflow != null && _overflow.Count > 0) || _terminatorChild != null;

        /// <summary>
        /// Gets the number of children.
        /// </summary>
        public int ChildCount => (_overflow != null ? _overflow.Count : _inlineCount) + (_terminatorChild != null ? 1 : 0);

        public SuffixTreeNode() { }

        public SuffixTreeNode(int start, int end, int depthFromRoot)
        {
            Start = start;
            End = end;
            DepthFromRoot = depthFromRoot;
        }

        /// <summary>
        /// Suffix link: connects node for "xα" to node for "α".
        /// Used for O(1) jumps between suffixes during construction.
        /// </summary>
        public SuffixTreeNode? SuffixLink { get; set; }

        /// <summary>
        /// Depth from root (sum of edge lengths on path from root to THIS node's parent,
        /// i.e., depth BEFORE this node's edge). Used for O(1) LRS lookups.
        /// </summary>
        public int DepthFromRoot { get; set; }

        /// <summary>
        /// Number of leaves in this node's subtree. For leaf nodes, this is 1.
        /// Calculated once after tree construction for O(1) CountOccurrences.
        /// </summary>
        public int LeafCount { get; set; }

        /// <summary>True if this is a leaf node (edge grows with string).</summary>
        public bool IsLeaf => End == BOUNDLESS;

        /// <summary>
        /// Tries to get a child by key. Key is char cast to int, or -1 for terminator.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool TryGetChild(int key, out SuffixTreeNode? child)
        {
            if (key < 0)
            {
                child = _terminatorChild;
                return child != null;
            }

            char cKey = (char)key;
            if (_overflow != null)
                return _overflow.TryGetValue(cKey, out child);

            if (_inlineCount > 0 && _key0 == cKey) { child = _child0; return true; }
            if (_inlineCount > 1 && _key1 == cKey) { child = _child1; return true; }
            if (_inlineCount > 2 && _key2 == cKey) { child = _child2; return true; }
            if (_inlineCount > 3 && _key3 == cKey) { child = _child3; return true; }

            child = null;
            return false;
        }

        /// <summary>
        /// Sets a child by key. Key is char cast to int, or -1 for terminator.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void SetChild(int key, SuffixTreeNode child)
        {
            if (key < 0)
            {
                _terminatorChild = child;
                return;
            }

            char cKey = (char)key;
            if (_overflow != null)
            {
                _overflow[cKey] = child;
                return;
            }

            // Check if key already exists in inline storage
            if (_inlineCount > 0 && _key0 == cKey) { _child0 = child; return; }
            if (_inlineCount > 1 && _key1 == cKey) { _child1 = child; return; }
            if (_inlineCount > 2 && _key2 == cKey) { _child2 = child; return; }
            if (_inlineCount > 3 && _key3 == cKey) { _child3 = child; return; }

            // Add new child
            if (_inlineCount < INLINE_THRESHOLD)
            {
                switch (_inlineCount)
                {
                    case 0: _key0 = cKey; _child0 = child; break;
                    case 1: _key1 = cKey; _child1 = child; break;
                    case 2: _key2 = cKey; _child2 = child; break;
                    case 3: _key3 = cKey; _child3 = child; break;
                }
                _inlineCount++;
            }
            else
            {
                // Promote to Dictionary
                _overflow = new Dictionary<char, SuffixTreeNode>(8)
                {
                    [_key0] = _child0!,
                    [_key1] = _child1!,
                    [_key2] = _child2!,
                    [_key3] = _child3!,
                    [cKey] = child
                };
                // Clear inline storage (helps GC)
                _child0 = _child1 = _child2 = _child3 = null;
            }
        }

        /// <summary>
        /// Populates the provided list with all children. 
        /// Avoids enumerator allocation.
        /// </summary>
        public void GetChildren(List<SuffixTreeNode> buffer)
        {
            buffer.Clear();
            if (_terminatorChild != null)
                buffer.Add(_terminatorChild);

            if (_overflow != null)
            {
                foreach (var child in _overflow.Values)
                    buffer.Add(child);
            }
            else
            {
                if (_inlineCount > 0) buffer.Add(_child0!);
                if (_inlineCount > 1) buffer.Add(_child1!);
                if (_inlineCount > 2) buffer.Add(_child2!);
                if (_inlineCount > 3) buffer.Add(_child3!);
            }
        }

        /// <summary>
        /// Enumerates all children using yield. Keep for compatibility/linq but use GetChildren(List) for performance.
        /// </summary>
        public IEnumerable<SuffixTreeNode> GetChildren()
        {
            if (_terminatorChild != null)
                yield return _terminatorChild;

            if (_overflow != null)
            {
                foreach (var child in _overflow.Values)
                    yield return child;
            }
            else
            {
                if (_inlineCount > 0) yield return _child0!;
                if (_inlineCount > 1) yield return _child1!;
                if (_inlineCount > 2) yield return _child2!;
                if (_inlineCount > 3) yield return _child3!;
            }
        }

        /// <summary>
        /// Enumerates all (key, child) pairs. Key is int (-1 for terminator).
        /// </summary>
        public IEnumerable<(int Key, SuffixTreeNode Child)> GetChildrenWithKeys()
        {
            if (_terminatorChild != null)
                yield return (-1, _terminatorChild);

            if (_overflow != null)
            {
                foreach (var kvp in _overflow)
                    yield return (kvp.Key, kvp.Value);
            }
            else
            {
                if (_inlineCount > 0) yield return (_key0, _child0!);
                if (_inlineCount > 1) yield return (_key1, _child1!);
                if (_inlineCount > 2) yield return (_key2, _child2!);
                if (_inlineCount > 3) yield return (_key3, _child3!);
            }
        }

        /// <summary>
        /// Gets all child keys for sorting/enumeration.
        /// </summary>
        public void GetKeys(List<int> keys)
        {
            keys.Clear();
            if (_terminatorChild != null)
                keys.Add(-1);

            if (_overflow != null)
            {
                foreach (var key in _overflow.Keys)
                    keys.Add(key);
            }
            else
            {
                if (_inlineCount > 0) keys.Add(_key0);
                if (_inlineCount > 1) keys.Add(_key1);
                if (_inlineCount > 2) keys.Add(_key2);
                if (_inlineCount > 3) keys.Add(_key3);
            }
        }
    }
}
