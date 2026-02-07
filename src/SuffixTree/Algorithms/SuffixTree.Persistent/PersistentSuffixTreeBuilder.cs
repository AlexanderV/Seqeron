using System;
using System.Buffers;
using System.Buffers.Binary;
using System.Collections.Generic;
using System.Text;

namespace SuffixTree.Persistent;

/// <summary>
/// Handles construction of a persistent suffix tree using Ukkonen's algorithm.
/// Writes nodes and child entries directly to the storage provider.
/// </summary>
public class PersistentSuffixTreeBuilder
{
    private readonly IStorageProvider _storage;
    private long _rootOffset;
    private ITextSource _text = new StringTextSource(string.Empty);
    private int _nodeCount = 0;

    // Algorithm state
    private long _activeNodeOffset;
    private int _activeEdgeIndex;
    private int _activeLength;
    private int _remainder;
    private int _position = -1;
    private long _lastCreatedInternalNodeOffset = PersistentConstants.NULL_OFFSET;
    private readonly Dictionary<long, List<(uint Key, long ChildOffset)>> _children = new();
    private bool _built;

    public PersistentSuffixTreeBuilder(IStorageProvider storage)
    {
        _storage = storage;
        _storage.Allocate(PersistentConstants.HEADER_SIZE); // Reserve space for header

        _rootOffset = _storage.Allocate(PersistentConstants.NODE_SIZE);
        _nodeCount = 1;
        var root = new PersistentSuffixTreeNode(_storage, _rootOffset);
        root.Start = 0;
        root.End = 0;
        root.SuffixLink = PersistentConstants.NULL_OFFSET;
        root.ChildrenHead = PersistentConstants.NULL_OFFSET;

        _activeNodeOffset = _rootOffset;
    }

    public long Build(ITextSource text)
    {
        if (_built)
            throw new InvalidOperationException("Build() has already been called. Create a new builder instance to build another tree.");
        _built = true;

        _text = text;
        if (text.Length > 0)
        {
            for (int i = 0; i < text.Length; i++)
            {
                ExtendTree((uint)text[i]);
            }
            ExtendTree(PersistentConstants.TERMINATOR_KEY);
        }

        FinalizeTree();
        return _rootOffset;
    }

    private void ExtendTree(uint key)
    {
        _position++;
        _remainder++;
        _lastCreatedInternalNodeOffset = PersistentConstants.NULL_OFFSET;

        while (_remainder > 0)
        {
            if (_activeLength == 0)
                _activeEdgeIndex = _position;

            var activeNode = new PersistentSuffixTreeNode(_storage, _activeNodeOffset);
            uint activeEdgeKey = GetSymbolAt(_activeEdgeIndex);

            if (!BuilderTryGetChild(_activeNodeOffset, activeEdgeKey, out var nextChildOffset))
            {
                var leafOffset = CreateNode((uint)_position, PersistentConstants.BOUNDLESS, GetNodeDepth(activeNode));
                BuilderSetChild(_activeNodeOffset, activeEdgeKey, leafOffset);
                AddSuffixLink(_activeNodeOffset);
            }
            else
            {
                var nextChild = new PersistentSuffixTreeNode(_storage, nextChildOffset);
                int edgeLen = LengthOf(nextChild);
                if (_activeLength >= edgeLen)
                {
                    _activeEdgeIndex += edgeLen;
                    _activeLength -= edgeLen;
                    _activeNodeOffset = nextChildOffset;
                    continue;
                }

                if (GetSymbolAt((int)(nextChild.Start + (uint)_activeLength)) == key)
                {
                    _activeLength++;
                    AddSuffixLink(_activeNodeOffset);
                    break;
                }

                // Split edge
                long splitOffset = CreateNode(nextChild.Start, nextChild.Start + (uint)_activeLength, nextChild.DepthFromRoot);
                var split = new PersistentSuffixTreeNode(_storage, splitOffset);
                BuilderSetChild(_activeNodeOffset, activeEdgeKey, splitOffset);

                long leafOffset = CreateNode((uint)_position, PersistentConstants.BOUNDLESS, split.DepthFromRoot + (uint)LengthOf(split));
                BuilderSetChild(splitOffset, key, leafOffset);

                nextChild.Start += (uint)_activeLength;
                nextChild.DepthFromRoot = split.DepthFromRoot + (uint)LengthOf(split);
                BuilderSetChild(splitOffset, GetSymbolAt((int)nextChild.Start), nextChildOffset);

                AddSuffixLink(splitOffset);
            }

            _remainder--;
            if (_activeNodeOffset == _rootOffset && _activeLength > 0)
            {
                _activeLength--;
                _activeEdgeIndex = _position - _remainder + 1;
            }
            else if (_activeNodeOffset != _rootOffset)
            {
                var node = new PersistentSuffixTreeNode(_storage, _activeNodeOffset);
                _activeNodeOffset = node.SuffixLink != PersistentConstants.NULL_OFFSET ? node.SuffixLink : _rootOffset;
            }
        }
    }

    private long CreateNode(uint start, uint end, uint depthFromRoot)
    {
        _nodeCount++;
        long offset = _storage.Allocate(PersistentConstants.NODE_SIZE);
        var node = new PersistentSuffixTreeNode(_storage, offset);
        node.Start = start;
        node.End = end;
        node.DepthFromRoot = depthFromRoot;
        node.SuffixLink = PersistentConstants.NULL_OFFSET;
        node.ChildrenHead = PersistentConstants.NULL_OFFSET;
        node.LeafCount = 0;
        return offset;
    }

    private void AddSuffixLink(long nodeOffset)
    {
        if (_lastCreatedInternalNodeOffset != PersistentConstants.NULL_OFFSET)
        {
            var lastNode = new PersistentSuffixTreeNode(_storage, _lastCreatedInternalNodeOffset);
            lastNode.SuffixLink = nodeOffset;
        }
        _lastCreatedInternalNodeOffset = nodeOffset;
    }

    private uint GetSymbolAt(int index)
    {
        if (index > _position) return PersistentConstants.TERMINATOR_KEY;
        return (index < _text.Length) ? (uint)_text[index] : PersistentConstants.TERMINATOR_KEY;
    }

    private int LengthOf(PersistentSuffixTreeNode node)
        => (int)((node.End == PersistentConstants.BOUNDLESS ? (uint)(_position + 1) : node.End) - node.Start);

    private uint GetNodeDepth(PersistentSuffixTreeNode node)
        => node.DepthFromRoot + (uint)LengthOf(node);

    private void FinalizeTree()
    {
        // Calculate leaf counts
        CalculateLeafCount(_rootOffset);

        // Write sorted children arrays to storage
        WriteChildrenArrays();

        // Store text in storage for true persistence (chunked write, no full-string copy)
        long textOffset = _storage.Allocate(_text.Length * 2);
        const int ChunkChars = 4096;
        byte[] chunkBuf = ArrayPool<byte>.Shared.Rent(ChunkChars * 2);
        try
        {
            int written = 0;
            while (written < _text.Length)
            {
                int remaining = _text.Length - written;
                int chunkLen = remaining < ChunkChars ? remaining : ChunkChars;
                int byteCount = Encoding.Unicode.GetBytes(
                    _text.Slice(written, chunkLen).ToString(), chunkBuf);
                _storage.WriteBytes(textOffset + written * 2, chunkBuf, 0, byteCount);
                written += chunkLen;
            }
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(chunkBuf);
        }

        // Write Header
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_MAGIC, PersistentConstants.MAGIC_NUMBER);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_VERSION, PersistentConstants.CURRENT_VERSION);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, _rootOffset);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF, textOffset);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN, _text.Length);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_NODE_COUNT, _nodeCount);
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_SIZE, _storage.Size);
    }

    private void CalculateLeafCount(long rootOffset)
    {
        var workStack = new Stack<long>();
        var resultStack = new Stack<long>();

        workStack.Push(rootOffset);

        while (workStack.Count > 0)
        {
            long offset = workStack.Pop();
            resultStack.Push(offset);

            if (_children.TryGetValue(offset, out var childList))
            {
                foreach (var (_, childOffset) in childList)
                    workStack.Push(childOffset);
            }
        }

        while (resultStack.Count > 0)
        {
            long offset = resultStack.Pop();
            var node = new PersistentSuffixTreeNode(_storage, offset);

            if (node.IsLeaf)
            {
                node.LeafCount = 1;
            }
            else
            {
                uint totalLeaves = 0;
                if (_children.TryGetValue(offset, out var childList))
                {
                    foreach (var (_, childOffset) in childList)
                    {
                        var child = new PersistentSuffixTreeNode(_storage, childOffset);
                        totalLeaves += child.LeafCount;
                    }
                }
                node.LeafCount = totalLeaves;
            }
        }
    }

    private bool BuilderTryGetChild(long nodeOffset, uint key, out long childOffset)
    {
        if (_children.TryGetValue(nodeOffset, out var childList))
        {
            foreach (var entry in childList)
            {
                if (entry.Key == key)
                {
                    childOffset = entry.ChildOffset;
                    return true;
                }
            }
        }
        childOffset = PersistentConstants.NULL_OFFSET;
        return false;
    }

    private void BuilderSetChild(long nodeOffset, uint key, long childOffset)
    {
        if (!_children.TryGetValue(nodeOffset, out var childList))
        {
            childList = new List<(uint Key, long ChildOffset)>();
            _children[nodeOffset] = childList;
        }

        for (int i = 0; i < childList.Count; i++)
        {
            if (childList[i].Key == key)
            {
                childList[i] = (key, childOffset);
                return;
            }
        }

        childList.Add((key, childOffset));
    }

    private void WriteChildrenArrays()
    {
        foreach (var (nodeOffset, childList) in _children)
        {
            // Sort by key using signed comparison (terminator=-1 first)
            childList.Sort((a, b) => ((int)a.Key).CompareTo((int)b.Key));

            int count = childList.Count;
            int totalBytes = count * PersistentConstants.CHILD_ENTRY_SIZE;
            long arrayOffset = _storage.Allocate(totalBytes);

            // Serialize all entries into a single buffer, then batch-write
            byte[] buf = ArrayPool<byte>.Shared.Rent(totalBytes);
            try
            {
                for (int i = 0; i < count; i++)
                {
                    int off = i * PersistentConstants.CHILD_ENTRY_SIZE;
                    BinaryPrimitives.WriteUInt32LittleEndian(buf.AsSpan(off, 4), childList[i].Key);
                    BinaryPrimitives.WriteInt64LittleEndian(buf.AsSpan(off + 4, 8), childList[i].ChildOffset);
                }
                _storage.WriteBytes(arrayOffset, buf, 0, totalBytes);
            }
            finally
            {
                ArrayPool<byte>.Shared.Return(buf);
            }

            var node = new PersistentSuffixTreeNode(_storage, nodeOffset);
            node.ChildrenHead = arrayOffset;
            node.ChildCount = count;
        }
    }
}
