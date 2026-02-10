namespace SuffixTree;

public partial class SuffixTree
{

    private void ExtendTree(int key)
    {
        _position++;
        _remainder++;
        _lastCreatedInternalNode = null;

        while (_remainder > 0)
        {
            if (_activeLength == 0)
                _activeEdgeIndex = _position;

            if (!_activeNode!.TryGetChild(GetSymbolAt(_activeEdgeIndex), out var nextChild) || nextChild == null)
            {
                var leaf = new SuffixTreeNode(_position, SuffixTreeNode.BOUNDLESS, GetNodeDepth(_activeNode));
                _activeNode.SetChild(GetSymbolAt(_activeEdgeIndex), leaf);
                AddSuffixLink(_activeNode);
            }
            else
            {
                int edgeLen = LengthOf(nextChild);
                if (_activeLength >= edgeLen)
                {
                    _activeEdgeIndex += edgeLen;
                    _activeLength -= edgeLen;
                    _activeNode = nextChild;
                    continue;
                }

                if (GetSymbolAt(nextChild.Start + _activeLength) == key)
                {
                    _activeLength++;
                    AddSuffixLink(_activeNode);
                    break;
                }

                var split = new SuffixTreeNode(nextChild.Start, nextChild.Start + _activeLength, nextChild.DepthFromRoot);
                _activeNode.SetChild(GetSymbolAt(_activeEdgeIndex), split);

                var leaf = new SuffixTreeNode(_position, SuffixTreeNode.BOUNDLESS, split.DepthFromRoot + LengthOf(split));
                split.SetChild(key, leaf);

                nextChild.Start += _activeLength;
                nextChild.DepthFromRoot = split.DepthFromRoot + LengthOf(split);
                split.SetChild(GetSymbolAt(nextChild.Start), nextChild);

                AddSuffixLink(split);
            }

            _remainder--;
            if (_activeNode == _root && _activeLength > 0)
            {
                _activeLength--;
                _activeEdgeIndex = _position - _remainder + 1;
            }
            else if (_activeNode != _root)
            {
                _activeNode = _activeNode.SuffixLink ?? _root;
            }
        }
    }

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private void AddSuffixLink(SuffixTreeNode node)
    {
        if (_lastCreatedInternalNode != null)
            _lastCreatedInternalNode.SuffixLink = node;
        _lastCreatedInternalNode = node;
    }

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private int GetSymbolAt(int index)
    {
        if (index > _position) return TERMINATOR_KEY;
        var raw = _rawString;
        if (raw != null)
            return index < raw.Length ? raw[index] : TERMINATOR_KEY;
        return index < _text.Length ? _text[index] : TERMINATOR_KEY;
    }

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private int LengthOf(SuffixTreeNode edge)
        => (edge.End == SuffixTreeNode.BOUNDLESS ? _position + 1 : edge.End) - edge.Start;

    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private int FirstCharOf(SuffixTreeNode edge)
        => GetSymbolAt(edge.Start);
}
