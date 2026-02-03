using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;

namespace SuffixTree
{
    public partial class SuffixTree
    {
        [ThreadStatic]
        private static List<SuffixTreeNode>? _sharedSearchBuffer;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static List<SuffixTreeNode> GetSearchBuffer()
        {
            _sharedSearchBuffer ??= new List<SuffixTreeNode>(64);
            _sharedSearchBuffer.Clear();
            return _sharedSearchBuffer;
        }

        public bool Contains(string value)
        {
            ArgumentNullException.ThrowIfNull(value);
            if (value.Length == 0) return true;
            return ContainsCore(value.AsSpan());
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private bool ContainsCore(ReadOnlySpan<char> pattern)
        {
            var node = _root;
            int i = 0;

            while (i < pattern.Length)
            {
                if (!node.TryGetChild(pattern[i], out var child))
                    return false;

                int edgeStart = child!.Start;
                int edgeLength = LengthOf(child);
                int remaining = pattern.Length - i;
                int compareLen = edgeLength < remaining ? edgeLength : remaining;

                if (edgeStart + compareLen > _text.Length)
                    return false;

                var edgeSpan = _text.Slice(edgeStart, compareLen);

                // Use SIMD for longer comparisons (>=8 chars), scalar for short
                if (compareLen >= 8)
                {
                    if (!edgeSpan.SequenceEqual(pattern.Slice(i, compareLen)))
                        return false;
                }
                else
                {
                    for (int j = 0; j < compareLen; j++)
                    {
                        if (edgeSpan[j] != pattern[i + j])
                            return false;
                    }
                }

                i += compareLen;
                node = child;
            }
            return true;
        }

        /// <summary>
        /// Core pattern matching with hybrid SIMD optimization.
        /// Returns the node where pattern ends and whether full match was found.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private (SuffixTreeNode node, bool matched) MatchPatternCore(ReadOnlySpan<char> pattern)
        {
            var node = _root;
            int i = 0;

            while (i < pattern.Length)
            {
                if (!node.TryGetChild(pattern[i], out var child))
                    return (node, false);

                int edgeStart = child!.Start;
                int edgeLength = LengthOf(child);
                int remaining = pattern.Length - i;
                int compareLen = edgeLength < remaining ? edgeLength : remaining;

                if (edgeStart + compareLen > _text.Length)
                    return (node, false);

                var edgeSpan = _text.Slice(edgeStart, compareLen);

                // Use SIMD for longer comparisons (>=8 chars), scalar for short
                if (compareLen >= 8)
                {
                    if (!edgeSpan.SequenceEqual(pattern.Slice(i, compareLen)))
                        return (node, false);
                }
                else
                {
                    for (int j = 0; j < compareLen; j++)
                    {
                        if (edgeSpan[j] != pattern[i + j])
                            return (node, false);
                    }
                }

                i += compareLen;
                node = child;
            }
            return (node, true);
        }

        public IReadOnlyList<int> FindAllOccurrences(string pattern)
        {
            ArgumentNullException.ThrowIfNull(pattern);
            var results = new List<int>();
            if (pattern.Length == 0)
            {
                for (int m = 0; m < _text.Length; m++) results.Add(m);
                return results;
            }

            var (node, matched) = MatchPatternCore(pattern.AsSpan());
            if (!matched) return results;

            CollectLeaves(node, node.DepthFromRoot, results);
            return results;
        }

        public int CountOccurrences(string pattern)
        {
            ArgumentNullException.ThrowIfNull(pattern);
            if (pattern.Length == 0) return _text.Length;

            var (node, matched) = MatchPatternCore(pattern.AsSpan());
            return matched ? node.LeafCount : 0;
        }

        [ThreadStatic]
        private static Stack<(SuffixTreeNode Node, int Depth)>? _sharedStack;

        private void CollectLeaves(SuffixTreeNode node, int depth, List<int> results)
        {
            _sharedStack ??= new Stack<(SuffixTreeNode Node, int Depth)>(64);
            var stack = _sharedStack;
            stack.Clear();
            stack.Push((node, depth));

            var buffer = GetSearchBuffer();

            while (stack.Count > 0)
            {
                var (current, currentDepthBefore) = stack.Pop();
                int currentDepth = currentDepthBefore + LengthOf(current);

                if (current.IsLeaf)
                {
                    int suffixLength = currentDepth;
                    int startPosition = _text!.Length + 1 - suffixLength;
                    if (startPosition < _text.Length)
                        results.Add(startPosition);
                }
                else
                {
                    current.GetChildren(buffer);
                    foreach (var child in buffer)
                    {
                        stack.Push((child, currentDepth));
                    }
                }
            }
        }

        public bool Contains(ReadOnlySpan<char> value)
        {
            if (value.IsEmpty) return true;
            return ContainsCore(value);
        }

        public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
        {
            var results = new List<int>();
            if (pattern.IsEmpty) return results;

            var (node, matched) = MatchPatternCore(pattern);
            if (!matched) return results;

            CollectLeaves(node, node.DepthFromRoot, results);
            return results;
        }

        public int CountOccurrences(ReadOnlySpan<char> pattern)
        {
            if (pattern.IsEmpty) return 0;

            var (node, matched) = MatchPatternCore(pattern);
            return matched ? node.LeafCount : 0;
        }

        public IEnumerable<string> EnumerateSuffixes() => EnumerateSuffixesCore();

        public IReadOnlyList<string> GetAllSuffixes() => EnumerateSuffixesCore().ToList();

        private IEnumerable<string> EnumerateSuffixesCore()
        {
            var stack = new Stack<(SuffixTreeNode Node, int ChildIndex, List<int> SortedKeys)>();
            var path = new StringBuilder(_text!.Length);
            var keyBuffer = new List<int>(8);

            _root.GetKeys(keyBuffer);
            var rootKeys = new List<int>(keyBuffer);
            rootKeys.Sort();
            stack.Push((_root, 0, rootKeys));

            while (stack.Count > 0)
            {
                var (node, childIndex, sortedKeys) = stack.Pop();

                if (childIndex < sortedKeys.Count)
                {
                    stack.Push((node, childIndex + 1, sortedKeys));

                    var childKey = sortedKeys[childIndex];
                    if (!node.TryGetChild(childKey, out var child) || child is null)
                        continue;

                    int edgeLen = LengthOf(child);
                    int charsAdded = 0;
                    for (int i = 0; i < edgeLen; i++)
                    {
                        int s = GetSymbolAt(child.Start + i);
                        if (s == TERMINATOR_KEY) break;
                        path.Append((char)s);
                        charsAdded++;
                    }

                    if (child.IsLeaf)
                    {
                        if (path.Length > 0)
                            yield return path.ToString();
                        path.Length -= charsAdded;
                    }
                    else
                    {
                        child.GetKeys(keyBuffer);
                        var childKeys = new List<int>(keyBuffer);
                        childKeys.Sort();
                        stack.Push((child, 0, childKeys));
                    }
                }
                else
                {
                    if (node != _root)
                    {
                        int edgeLen = LengthOf(node);
                        int charsToRemove = 0;
                        for (int i = 0; i < edgeLen; i++)
                        {
                            int s = GetSymbolAt(node.Start + i);
                            if (s == TERMINATOR_KEY) break;
                            charsToRemove++;
                        }
                        if (charsToRemove > 0 && path.Length >= charsToRemove)
                            path.Length -= charsToRemove;
                    }
                }
            }
        }
    }
}
