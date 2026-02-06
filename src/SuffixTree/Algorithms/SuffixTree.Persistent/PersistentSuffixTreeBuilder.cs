using System;
using System.Collections.Generic;

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
            if (!activeNode.TryGetChild(GetSymbolAt(_activeEdgeIndex), out var nextChild) || nextChild.IsNull)
            {
                var leafOffset = CreateNode((uint)_position, PersistentConstants.BOUNDLESS, GetNodeDepth(activeNode));
                activeNode.SetChild(GetSymbolAt(_activeEdgeIndex), new PersistentSuffixTreeNode(_storage, leafOffset));
                AddSuffixLink(_activeNodeOffset);
            }
            else
            {
                int edgeLen = LengthOf(nextChild);
                if (_activeLength >= edgeLen)
                {
                    _activeEdgeIndex += edgeLen;
                    _activeLength -= edgeLen;
                    _activeNodeOffset = nextChild.Offset;
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
                activeNode.SetChild(GetSymbolAt(_activeEdgeIndex), split);

                long leafOffset = CreateNode((uint)_position, PersistentConstants.BOUNDLESS, split.DepthFromRoot + (uint)LengthOf(split));
                split.SetChild(key, new PersistentSuffixTreeNode(_storage, leafOffset));

                nextChild.Start += (uint)_activeLength;
                nextChild.DepthFromRoot = split.DepthFromRoot + (uint)LengthOf(split);
                split.SetChild(GetSymbolAt((int)nextChild.Start), nextChild);

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
        // Calculate leaf counts recursively
        CalculateLeafCount(_rootOffset);

        // Store text in storage for true persistence
        long textOffset = _storage.Allocate(_text.Length * 2);
        for (int i = 0; i < _text.Length; i++)
        {
            _storage.WriteChar(textOffset + (i * 2), _text[i]);
        }

        // Write Header
        _storage.WriteInt64(PersistentConstants.HEADER_OFFSET_MAGIC, PersistentConstants.MAGIC_NUMBER);
        _storage.WriteInt32(PersistentConstants.HEADER_OFFSET_VERSION, 1);
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

            var node = new PersistentSuffixTreeNode(_storage, offset);
            long childEntryOffset = node.ChildrenHead;
            while (childEntryOffset != PersistentConstants.NULL_OFFSET)
            {
                var entry = new PersistentChildEntry(_storage, childEntryOffset);
                workStack.Push(entry.ChildNodeOffset);
                childEntryOffset = entry.NextEntryOffset;
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
                long childEntryOffset = node.ChildrenHead;
                while (childEntryOffset != PersistentConstants.NULL_OFFSET)
                {
                    var entry = new PersistentChildEntry(_storage, childEntryOffset);
                    var child = new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset);
                    totalLeaves += child.LeafCount;
                    childEntryOffset = entry.NextEntryOffset;
                }
                node.LeafCount = totalLeaves;
            }
        }
    }
}
