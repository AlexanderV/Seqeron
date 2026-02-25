using System.Text;

namespace SuffixTree.Persistent;

public sealed partial class PersistentSuffixTree
{
    /// <inheritdoc />
    public string PrintTree()
    {
        ThrowIfDisposed();
        var sb = new StringBuilder(Math.Max(256, Math.Min(_textSource.Length, 100_000) * 100));
        var ci = System.Globalization.CultureInfo.InvariantCulture;
        sb.Append(ci, $"Content length: {_textSource.Length}").AppendLine();
        sb.AppendLine();
        sb.Append(ci, $"0:ROOT").AppendLine();

        var root = NodeAt(_rootOffset);
        int rootChildCount;
        long rootArrayBase;
        NodeLayout rootEntryLayout;
        {
            var (ab, el, cc) = ReadChildArrayInfo(root);
            rootChildCount = cc;
            rootArrayBase = ab;
            rootEntryLayout = el;
        }

        // Frame: (parent arrayBase, entryLayout, child count, current index, depth)
        var stack = new Stack<(long ArrayBase, NodeLayout EntryLayout, int ChildCount, int Index, int Depth)>();
        if (rootChildCount > 0)
            stack.Push((rootArrayBase, rootEntryLayout, rootChildCount, 0, 0));

        while (stack.Count > 0)
        {
            var (arrBase, entryLay, childCount, index, depth) = stack.Pop();
            if (index >= childCount) continue;

            // Push continuation for next sibling
            stack.Push((arrBase, entryLay, childCount, index + 1, depth));

            // Read child directly from sorted array
            long entryOffset = arrBase + (long)index * entryLay.ChildEntrySize;
            long childNodeOffset = entryLay.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
            var child = NodeAt(childNodeOffset);

            int childDepth = depth + 1;
            sb.Append(' ', childDepth * 2);
            sb.Append(ci, $"{childDepth}: ");
            AppendLabel(sb, child);

            if (child.IsLeaf)
            {
                sb.AppendLine(" (Leaf)");
            }
            else
            {
                long suffLink = ResolveSuffixLink(child);
                if (suffLink != PersistentConstants.NULL_OFFSET)
                {
                    var linkNode = NodeAt(suffLink);
                    if (linkNode.Offset != _rootOffset)
                    {
                        int firstChar = GetSymbolAt((int)linkNode.Start);
                        if (firstChar >= 0)
                            sb.Append(ci, $" -> [{(char)firstChar}]");
                    }
                }
                sb.AppendLine();

                var (cBase, cLayout, cCount) = ReadChildArrayInfo(child);
                if (cCount > 0)
                    stack.Push((cBase, cLayout, cCount, 0, childDepth));
            }
        }

        return sb.ToString();
    }

    private void AppendLabel(StringBuilder sb, PersistentSuffixTreeNode node)
    {
        int len = LengthOf(node);
        for (int i = 0; i < len; i++)
        {
            int s = GetSymbolAt((int)node.Start + i);
            if (s == -1) { sb.Append('#'); break; }
            sb.Append((char)s);
        }
    }

    /// <inheritdoc/>
    public void Traverse(ISuffixTreeVisitor visitor)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(visitor);
        TraverseCore(NodeAt(_rootOffset), 0, visitor);
    }

    private void TraverseCore(PersistentSuffixTreeNode node, int depth, ISuffixTreeVisitor visitor)
    {
        // Frame: arrayBase, entryLayout, child count, current child index, depth (character-based)
        var stack = new Stack<(long ArrayBase, NodeLayout EntryLayout, int ChildCount, int Index, int Depth)>();

        // Visit root node
        var (rootAB, rootEL, rootCC) = ReadChildArrayInfo(node);
        visitor.VisitNode((int)node.Start, (int)node.End, (int)node.LeafCount, rootCC, depth);
        if (!node.IsLeaf)
        {
            stack.Push((rootAB, rootEL, rootCC, 0, depth));
        }

        while (stack.Count > 0)
        {
            var (arrBase, entryLay, childCount, index, parentDepth) = stack.Pop();

            if (index >= childCount)
            {
                // All children processed — exit branch (unless this is the root frame)
                if (stack.Count > 0)
                    visitor.ExitBranch();
                continue;
            }

            // Push continuation for next sibling
            stack.Push((arrBase, entryLay, childCount, index + 1, parentDepth));

            // Read child directly from sorted array
            long entryOffset = arrBase + (long)index * entryLay.ChildEntrySize;
            uint key = _storage.ReadUInt32(entryOffset + NodeLayout.ChildOffsetKey);
            long childOffset = entryLay.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
            var child = NodeAt(childOffset);

            var (cAB, cEL, cCC) = ReadChildArrayInfo(child);
            visitor.EnterBranch((int)key);
            visitor.VisitNode((int)child.Start, (int)child.End, (int)child.LeafCount, cCC, parentDepth);

            int childDepth = parentDepth + LengthOf(child);

            if (!child.IsLeaf)
            {
                stack.Push((cAB, cEL, cCC, 0, childDepth));
            }
            else
            {
                visitor.ExitBranch();
            }
        }
    }
}
