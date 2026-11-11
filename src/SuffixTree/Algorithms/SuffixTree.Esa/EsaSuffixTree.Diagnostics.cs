using System.Globalization;
using System.Text;

namespace SuffixTree.Esa;

public sealed partial class EsaSuffixTree
{
    private const int PrintTreeInitialCapacityLimit = 100_000;

    /// <inheritdoc />
    /// <remarks>
    /// The layout of the other backends: the header <c>Content length: n</c>, a blank line, <c>0:ROOT</c>, then one line per node, indented by two spaces per depth, children in key order
    /// (the terminator first), <c>(Leaf)</c> after a leaf, and for an internal node whose suffix link is not the root <c>-&gt; [c]</c> with <c>c</c> the first symbol of the incoming edge of the
    /// link target (as the persistent tree prints it). The link marks need the navigation pack, which the first call builds.
    /// </remarks>
    public string PrintTree()
    {
        ThrowIfDisposed();
        int n = _index.Length;
        var ci = CultureInfo.InvariantCulture;
        var sb = new StringBuilder(Math.Max(256, Math.Min(n, PrintTreeInitialCapacityLimit) * 100));
        sb.Append(ci, $"Content length: {n}").AppendLine();
        sb.AppendLine();
        sb.Append(ci, $"0:ROOT").AppendLine();
        if (n == 0) return sb.ToString();

        EsaNavigator nav = _index.CreateNavigator();
        var stack = new Stack<(List<EsaNode> Children, int Index, int Depth)>();
        var rootChildren = new List<EsaNode>();
        nav.GetChildren(nav.Root, rootChildren);
        stack.Push((rootChildren, 0, 0));
        while (stack.Count > 0)
        {
            var (children, index, depth) = stack.Pop();
            if (index >= children.Count) continue;
            stack.Push((children, index + 1, depth));

            EsaNode child = children[index];
            int childDepth = depth + 1;
            sb.Append(' ', childDepth * 2);
            sb.Append(ci, $"{childDepth}: ");
            AppendLabel(sb, in nav, child);
            if (child.IsLeaf)
            {
                sb.AppendLine(" (Leaf)");
                continue;
            }

            EsaNode link = nav.GetSuffixLink(child);
            if (!nav.IsRoot(link)) sb.Append(ci, $" -> [{(char)nav.GetEdgeSymbol(link, 0)}]");
            sb.AppendLine();
            var grandChildren = new List<EsaNode>();
            nav.GetChildren(child, grandChildren);
            stack.Push((grandChildren, 0, childDepth));
        }
        return sb.ToString();
    }

    private static void AppendLabel(StringBuilder sb, in EsaNavigator nav, EsaNode node)
    {
        int length = nav.LengthOf(node);
        for (int i = 0; i < length; i++)
        {
            int symbol = nav.GetEdgeSymbol(node, i);
            if (symbol < 0)
            {
                sb.Append('#');
                return;
            }
            sb.Append((char)symbol);
        }
    }

    /// <inheritdoc />
    public void Traverse(ISuffixTreeVisitor visitor)
    {
        ThrowIfDisposed();
        _index.Traverse(visitor);
    }
}
