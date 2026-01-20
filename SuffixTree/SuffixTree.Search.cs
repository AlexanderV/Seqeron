using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SuffixTree
{
    public partial class SuffixTree
    {
        [ThreadStatic]
        private static List<SuffixTreeNode>? _sharedSearchBuffer;

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
            var node = _root;
            int i = 0;

            while (i < value.Length)
            {
                if (!node.TryGetChild(value[i], out var child) || child is null)
                    return false;

                int edgeLength = LengthOf(child);
                int j = 0;
                while (j < edgeLength && i < value.Length)
                {
                    if (GetSymbolAt(child.Start + j) != value[i])
                        return false;
                    i++;
                    j++;
                }
                node = child;
            }
            return true;
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

            var node = _root;
            int i = 0;

            while (i < pattern.Length)
            {
                if (!node.TryGetChild(pattern[i], out var child) || child is null)
                    return results;

                int edgeLength = LengthOf(child);
                int j = 0;
                while (j < edgeLength && i < pattern.Length)
                {
                    if (GetSymbolAt(child.Start + j) != pattern[i])
                        return results;
                    i++;
                    j++;
                }

                if (j == edgeLength)
                {
                    node = child;
                }
                else
                {
                    CollectLeaves(child, child.DepthFromRoot, results);
                    return results;
                }
            }

            CollectLeaves(node, node.DepthFromRoot, results);
            return results;
        }

        public int CountOccurrences(string pattern)
        {
            ArgumentNullException.ThrowIfNull(pattern);
            if (pattern.Length == 0) return _text.Length;

            var node = _root;
            int i = 0;

            while (i < pattern.Length)
            {
                if (!node.TryGetChild(pattern[i], out var child) || child is null)
                    return 0;

                int edgeLength = LengthOf(child);
                int j = 0;
                while (j < edgeLength && i < pattern.Length)
                {
                    if (GetSymbolAt(child.Start + j) != pattern[i])
                        return 0;
                    i++;
                    j++;
                }
                node = child;
            }

            return node.LeafCount;
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
            var node = _root;
            int i = 0;

            while (i < value.Length)
            {
                if (!node.TryGetChild(value[i], out var child) || child is null)
                    return false;

                int edgeLength = LengthOf(child);
                int j = 0;
                while (j < edgeLength && i < value.Length)
                {
                    if (GetSymbolAt(child.Start + j) != value[i])
                        return false;
                    i++;
                    j++;
                }
                node = child;
            }
            return true;
        }

        public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
        {
            var results = new List<int>();
            if (pattern.IsEmpty) return results;

            var node = _root;
            int i = 0;

            while (i < pattern.Length)
            {
                if (!node.TryGetChild(pattern[i], out var child) || child is null)
                    return results;

                int edgeLength = LengthOf(child);
                int j = 0;
                while (j < edgeLength && i < pattern.Length)
                {
                    if (GetSymbolAt(child.Start + j) != pattern[i])
                        return results;
                    i++;
                    j++;
                }

                if (j == edgeLength)
                {
                    node = child;
                }
                else
                {
                    CollectLeaves(child, child.DepthFromRoot, results);
                    return results;
                }
            }

            CollectLeaves(node, node.DepthFromRoot, results);
            return results;
        }

        public int CountOccurrences(ReadOnlySpan<char> pattern)
        {
            if (pattern.IsEmpty) return 0;

            var node = _root;
            int i = 0;

            while (i < pattern.Length)
            {
                if (!node.TryGetChild(pattern[i], out var child) || child is null)
                    return 0;

                int edgeLength = LengthOf(child);
                int j = 0;
                while (j < edgeLength && i < pattern.Length)
                {
                    if (GetSymbolAt(child.Start + j) != pattern[i])
                        return 0;
                    i++;
                    j++;
                }
                node = child;
            }

            return node.LeafCount;
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
