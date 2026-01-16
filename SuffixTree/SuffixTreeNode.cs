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
    internal class SuffixTreeNode
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
        private SuffixTreeNode _child0, _child1, _child2, _child3;
        private Dictionary<char, SuffixTreeNode> _overflow;

        /// <summary>
        /// Returns true if this node has any children.
        /// </summary>
        public bool HasChildren => _inlineCount > 0 || (_overflow != null && _overflow.Count > 0);

        /// <summary>
        /// Gets the number of children.
        /// </summary>
        public int ChildCount => _overflow != null ? _overflow.Count : _inlineCount;

        /// <summary>
        /// Suffix link: connects node for "xα" to node for "α".
        /// Used for O(1) jumps between suffixes during construction.
        /// </summary>
        public SuffixTreeNode SuffixLink { get; set; }

        /// <summary>
        /// Parent node reference for O(depth) path reconstruction.
        /// Null for root node.
        /// </summary>
        public SuffixTreeNode Parent { get; set; }

        /// <summary>True if this is a leaf node (edge grows with string).</summary>
        public bool IsLeaf => End == BOUNDLESS;

        /// <summary>
        /// Tries to get a child by key character.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool TryGetChild(char key, out SuffixTreeNode child)
        {
            if (_overflow != null)
                return _overflow.TryGetValue(key, out child);

            if (_inlineCount > 0 && _key0 == key) { child = _child0; return true; }
            if (_inlineCount > 1 && _key1 == key) { child = _child1; return true; }
            if (_inlineCount > 2 && _key2 == key) { child = _child2; return true; }
            if (_inlineCount > 3 && _key3 == key) { child = _child3; return true; }

            child = null;
            return false;
        }

        /// <summary>
        /// Sets a child by key character.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void SetChild(char key, SuffixTreeNode child)
        {
            if (_overflow != null)
            {
                _overflow[key] = child;
                return;
            }

            // Check if key already exists in inline storage
            if (_inlineCount > 0 && _key0 == key) { _child0 = child; return; }
            if (_inlineCount > 1 && _key1 == key) { _child1 = child; return; }
            if (_inlineCount > 2 && _key2 == key) { _child2 = child; return; }
            if (_inlineCount > 3 && _key3 == key) { _child3 = child; return; }

            // Add new child
            if (_inlineCount < INLINE_THRESHOLD)
            {
                switch (_inlineCount)
                {
                    case 0: _key0 = key; _child0 = child; break;
                    case 1: _key1 = key; _child1 = child; break;
                    case 2: _key2 = key; _child2 = child; break;
                    case 3: _key3 = key; _child3 = child; break;
                }
                _inlineCount++;
            }
            else
            {
                // Promote to Dictionary
                _overflow = new Dictionary<char, SuffixTreeNode>(8)
                {
                    [_key0] = _child0,
                    [_key1] = _child1,
                    [_key2] = _child2,
                    [_key3] = _child3,
                    [key] = child
                };
                // Clear inline storage (optional, helps GC)
                _child0 = _child1 = _child2 = _child3 = null;
            }
        }

        /// <summary>
        /// Enumerates all children. Allocates only for Dictionary case.
        /// </summary>
        public IEnumerable<SuffixTreeNode> GetChildren()
        {
            if (_overflow != null)
            {
                foreach (var child in _overflow.Values)
                    yield return child;
            }
            else
            {
                if (_inlineCount > 0) yield return _child0;
                if (_inlineCount > 1) yield return _child1;
                if (_inlineCount > 2) yield return _child2;
                if (_inlineCount > 3) yield return _child3;
            }
        }

        /// <summary>
        /// Enumerates all (key, child) pairs.
        /// </summary>
        public IEnumerable<(char Key, SuffixTreeNode Child)> GetChildrenWithKeys()
        {
            if (_overflow != null)
            {
                foreach (var kvp in _overflow)
                    yield return (kvp.Key, kvp.Value);
            }
            else
            {
                if (_inlineCount > 0) yield return (_key0, _child0);
                if (_inlineCount > 1) yield return (_key1, _child1);
                if (_inlineCount > 2) yield return (_key2, _child2);
                if (_inlineCount > 3) yield return (_key3, _child3);
            }
        }

        /// <summary>
        /// Gets all child keys for sorting/enumeration.
        /// </summary>
        public void GetKeys(List<char> keys)
        {
            keys.Clear();
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

        // Legacy property for compatibility - avoid using in new code
        [System.Obsolete("Use TryGetChild/SetChild instead for better performance")]
        public Dictionary<char, SuffixTreeNode> Children
        {
            get
            {
                if (_overflow != null)
                    return _overflow;

                // Create Dictionary on demand (legacy compatibility)
                var dict = new Dictionary<char, SuffixTreeNode>(_inlineCount);
                if (_inlineCount > 0) dict[_key0] = _child0;
                if (_inlineCount > 1) dict[_key1] = _child1;
                if (_inlineCount > 2) dict[_key2] = _child2;
                if (_inlineCount > 3) dict[_key3] = _child3;
                return dict;
            }
        }
    }
}
