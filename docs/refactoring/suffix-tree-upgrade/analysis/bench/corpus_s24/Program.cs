// The runner of the corruption rows (work package S2.4): corpus <container> <Checked|Trusted>. Opens the container at the weakest level that accepts it, then runs a battery of public calls; prints one line:
//   refused                      no level accepts the file
//   answered <hash>              every call returned (the hash is of the results: equal builds give equal hashes)
//   invalid <call>               a call threw InvalidDataException (the damage was found by the read)
//   other <type> <call>          a call threw another exception (a finding in Checked mode)
// A crash of the process (AccessViolation, exit 134) and a hang are for the driver to see.
using System.Security.Cryptography;
using System.Text;
using SuffixTree;
using SuffixTree.Esa;
using SuffixTree.Indexing;

static class Program
{
    sealed class CountingVisitor : ISuffixTreeVisitor
    {
        public long Nodes, Branches;
        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Nodes += 1 + (leafCount & 1);
        public void EnterBranch(int key) => Branches += key & 7;
        public void ExitBranch() { }
    }

    static int Main(string[] args)
    {
        string path = args[0];
        var access = args[1] == "Trusted" ? EsaAccessMode.Trusted : EsaAccessMode.Checked;
        ISuffixTree? tree = null;
        foreach (EsaValidationLevel level in new[] { EsaValidationLevel.Header, EsaValidationLevel.Sections, EsaValidationLevel.Structure })
        {
            try
            {
                tree = SuffixIndexFactory.Open(path, new EsaOpenOptions { Validation = level, Access = access });
                break;
            }
            catch (Exception e) when (e is InvalidDataException or NotSupportedException)
            {
                // refused at this level: try the next (the weakest level that accepts the file is the interesting one)
            }
        }
        if (tree is null)
        {
            Console.WriteLine("refused");
            return 0;
        }
        var sb = new StringBuilder();
        string call = "open";
        try
        {
            string[] patterns = ["A", "C", "G", "T", "AC", "GT", "ACG", "TTT", "AAAA", "ACGTACGT", "CGCGCG", "L", "LL", "KLM", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "ACGTTGCAACGTTGCA", "Q", "ZZ"];
            foreach (string p in patterns)
            {
                call = "Contains"; sb.Append(tree.Contains(p)).Append(' ');
                call = "CountOccurrences"; sb.Append(tree.CountOccurrences(p)).Append(' ');
                call = "FindAllOccurrences"; foreach (int x in tree.FindAllOccurrences(p).Take(30)) sb.Append(x).Append(',');
            }
            call = "LongestRepeatedSubstring"; sb.Append(tree.LongestRepeatedSubstring().Length).Append(' ');
            call = "CountDistinctSubstrings"; sb.Append(tree.CountDistinctSubstrings()).Append(' ');
            call = "CountDistinctSubstringsByLength"; sb.Append(string.Join(',', tree.CountDistinctSubstringsByLength(8))).Append(' ');
            call = "FindAllLongestRepeatedSubstrings"; sb.Append(tree.FindAllLongestRepeatedSubstrings().Count).Append(' ');
            call = "GetAllSuffixes"; if (tree.Text.Length <= 5000) sb.Append(tree.GetAllSuffixes().Count).Append(' ');
            call = "Traverse"; var v = new CountingVisitor(); tree.Traverse(v); sb.Append(v.Nodes).Append(' ').Append(v.Branches).Append(' ');
            call = "PrintTree"; if (tree.Text.Length <= 3000) sb.Append(tree.PrintTree().Length).Append(' ');
            call = "LongestCommonSubstring"; sb.Append(tree.LongestCommonSubstring("ACGTACGTTGCAACGTACGTAAACCCGGGTTT").Length).Append(' ');
            call = "FindMaximalExactMatches"; sb.Append(tree.FindMaximalExactMatches("ACGTACGTTGCAACGTACGTAAACCCGGGTTT", 6).Count).Append(' ');
            call = "FindExactMatchAnchors"; sb.Append(tree.FindExactMatchAnchors("ACGTACGTTGCAACGTACGTAAACCCGGGTTT", 6).Count).Append(' ');
            if (tree is ISuffixIndex index)
            {
                call = "TextLength"; long n = index.TextLength; sb.Append(n).Append(' ');
                foreach (string p in new[] { "A", "AC", "ACG", "L", "AAAAAA" })
                {
                    call = "TryExtendRight"; bool ok = index.TryExtendRight(index.Root, p.AsSpan(), out IndexLocus locus); sb.Append(ok).Append(' ');
                    if (!ok) continue;
                    call = "Count"; sb.Append(index.Count(locus)).Append(' ');
                    var positions = new long[16];
                    call = "Locate"; sb.Append(index.Locate(locus, 1, positions)).Append(' '); foreach (long x in positions) sb.Append(x).Append(',');
                    var children = new IndexLocus[8];
                    call = "GetChildren"; if (index.IsExplicitNode(locus)) sb.Append(index.GetChildren(locus, children)).Append(' ');
                }
            }
            if (tree is ISuffixArrayAccess sa)
            {
                long n = ((ISuffixIndex)tree).TextLength;
                foreach (long r in new[] { 0L, 1, 2, n / 3, n / 2, n - 1, n })
                {
                    call = "SuffixAt"; sb.Append(sa.SuffixAt(r)).Append(' ');
                    call = "LcpAt"; sb.Append(sa.LcpAt(r)).Append(' ');
                }
            }
            if (tree is ISuffixTreeNavigation nav)
            {
                long n = ((ISuffixIndex)tree).TextLength;
                foreach (long r in new[] { 1L, n / 4, n / 2 })
                {
                    call = "RankOf"; sb.Append(nav.RankOf(r)).Append(' ');
                    call = "Lcp"; sb.Append(nav.Lcp(r, n - r)).Append(' ');
                    call = "Lca"; sb.Append(nav.Lca(r, n - r)).Append(' ');
                }
            }
        }
        catch (InvalidDataException)
        {
            Console.WriteLine("invalid " + call);
            return 0;
        }
        catch (Exception e)
        {
            Console.WriteLine("other " + e.GetType().Name + " " + call);
            return 0;
        }
        Console.WriteLine("answered " + Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(sb.ToString())))[..16]);
        return 0;
    }
}
