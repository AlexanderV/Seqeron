using System.Text;
using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// <see cref="ISuffixTreeDiagnostics.PrintTree"/> of the persistent v6 tree and of the ESA facade: the node lines equal the Ukkonen tree's, and every link mark (<c>-&gt; [c]</c>: the first
/// symbol of the incoming edge of the suffix-link target; none when the target is the root) follows from the definition of the suffix link on the printed tree itself.
/// </summary>
[TestFixtureSource(typeof(FacadePrintTreeTests), nameof(Fixtures))]
public class FacadePrintTreeTests(SuffixTreeBackend backend) : PortedBackendTestBase(backend)
{
    public static IEnumerable<TestFixtureData> Fixtures()
    {
        foreach (SuffixTreeBackend b in new[] { SuffixTreeBackends.PersistentHeap, SuffixTreeBackends.Esa, SuffixTreeBackends.EsaPacked, SuffixTreeBackends.EsaReference3 })
            yield return new TestFixtureData(b).SetArgDisplayNames(b.Name);
    }

    // a label with a line break or a ':' would break the line-oriented parse below
    private static bool IsPrintable(string text) => text.All(c => c is >= ' ' and < (char)0x7f && c != ':');

    private sealed record Line(int Depth, string Label, bool IsLeaf, char? Mark, string Path);

    private static List<Line> Parse(string printed)
    {
        var lines = new List<Line>();
        var paths = new List<string> { string.Empty };
        foreach (string raw in printed.Split('\n', StringSplitOptions.RemoveEmptyEntries).Select(l => l.TrimEnd('\r')).Skip(1))
        {
            int colon = raw.IndexOf(':');
            int depth = int.Parse(raw[..colon].Trim(), System.Globalization.CultureInfo.InvariantCulture);
            string content = raw[(colon + 1)..].TrimStart();
            bool leaf = content.EndsWith(" (Leaf)", StringComparison.Ordinal);
            if (leaf) content = content[..^" (Leaf)".Length];
            char? mark = null;
            int k = content.IndexOf(" -> [", StringComparison.Ordinal);
            if (k >= 0)
            {
                Assert.That(content, Does.EndWith("]"));
                mark = content[k + 5];
                Assert.That(content.Length, Is.EqualTo(k + 7), "the mark is one character");
                content = content[..k];
            }
            if (depth == 0)
            {
                paths = [string.Empty];
                lines.Add(new Line(0, string.Empty, false, null, string.Empty));
                continue;
            }
            string label = content;
            while (paths.Count > depth) paths.RemoveAt(paths.Count - 1);
            paths.Add(paths[depth - 1] + label);
            lines.Add(new Line(depth, label, leaf, mark, paths[^1]));
        }
        return lines;
    }

    [Test]
    public void NodeLines_EqualTheUkkonenTrees_OnceTheLinkMarksAreRemoved()
    {
        foreach (string text in new[] { "", "a", "ab", "banana", "mississippi", "abracadabra", "aaaaaaa", "abababab", "xabxac", "🧬αβγ🧪$", "a#b#a" }
                     .Concat(StringFamilies.Adversarial(60, 4).Select(f => f.Text)))
            FacadeChecks.CheckPrintTree(Tree(text), UkkonenTree.Build(text), text, Backend.Name);
    }

    [Test]
    public void LinkMarks_FollowTheDefinitionOfTheSuffixLink()
    {
        foreach (string text in new[] { "banana", "mississippi", "abracadabra", "aaaaaaa", "abababab", "xabxac", "abcabxabcd", "GATTACAGATTACA" }
                     .Concat(StringFamilies.Adversarial(50, 6).Select(f => f.Text).Where(IsPrintable)))
        {
            List<Line> lines = Parse(Tree(text).PrintTree());
            var byPath = lines.ToDictionary(l => l.Path);
            foreach (Line line in lines.Where(l => !l.IsLeaf))
            {
                string path = line.Path;
                if (path.Length <= 1)
                {
                    Assert.That(line.Mark, Is.Null, $"[{Backend.Name}] '{text}' node '{path}': the link of a depth-1 node is the root");
                    continue;
                }
                string target = path[1..];
                Assert.That(byPath.TryGetValue(target, out Line? link), Is.True, $"[{Backend.Name}] '{text}': the suffix-link target '{target}' of '{path}' is an explicit node");
                char expected = link!.Label[0];
                Assert.That(line.Mark, Is.EqualTo(expected), $"[{Backend.Name}] '{text}' node '{path}' -> '{target}': first symbol of the incoming edge");
            }
        }
    }

    [Test]
    public void Output_UsesTheDocumentedShape()
    {
        string printed = Tree("ab").PrintTree();
        string[] lines = printed.Split('\n', StringSplitOptions.RemoveEmptyEntries).Select(l => l.TrimEnd('\r')).ToArray();
        Assert.That(lines, Is.EqualTo(new[] { "Content length: 2", "0:ROOT", "  1: # (Leaf)", "  1: ab# (Leaf)", "  1: b# (Leaf)" }));
        Assert.That(printed, Does.Contain("\n\n").Or.Contain("\r\n\r\n"), "a blank line separates the header and the root");
        Assert.That(new StringBuilder(printed).ToString(), Does.StartWith("Content length: 2"));
    }
}
