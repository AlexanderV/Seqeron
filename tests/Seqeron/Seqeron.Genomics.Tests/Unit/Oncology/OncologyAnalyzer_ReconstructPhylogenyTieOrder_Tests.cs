// ONCO-PHYLO-001 — LICHeE equal-score tie order = Java HashMap iteration order of the presence profiles (B24 F45)
// Evidence: docs/Evidence/ONCO-PHYLO-001-Evidence.md (§ F45)
// TestSpec: tests/TestSpecs/ONCO-PHYLO-001.md (M28–M30)
// Source: github.com/viq854/lichee (master 26c2a70) SNVDataStore.loadSNVFileWithClusters (tag2SNVs HashMap<String, …>,
//         containsKey + put in clusters-file order), LineageEngine.buildLineage step 2 (groups = tag2SNVs.keySet() order),
//         PHYNetwork constructor (node ids in group order), PHYTree.compareTo + Collections.sort (stable ranking);
//         OpenJDK java.util.HashMap (putVal / resize / treeifyBin / TreeNode).
//
// Expected key orders are the output of java.util.HashMap (OpenJDK 21) for the same insertion sequence; expected trees
// are the output of lichee.jar driven by the F4445 harness with the HashMap group order (default JVM, no hash flags —
// String hash codes are deterministic). Not copied from the implementation.

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_ReconstructPhylogenyTieOrder_Tests
{
    private static OncologyAnalyzer.CcfCluster C(int id, params double[] ccf) => new(id, ccf);

    private static int?[] Parents(OncologyAnalyzer.ClonalPhylogeny p, int n) =>
        Enumerable.Range(1, n).Select(id => p.ParentOf(id) == p.RootId ? -1 : p.ParentOf(id)).ToArray();

    #region JavaStringHashMapOrder vs java.util.HashMap

    // Java: printf '<keys>' | HmOrder (HashMap<String,Integer>.put in order, then keySet()).
    [TestCase("001 111", "001 111")]
    [TestCase("10 01 11", "11 01 10")]
    [TestCase("101 110 100", "110 100 101")]
    [TestCase("00110 00111 11001 00001", "11001 00110 00111 00001")]
    [TestCase("1111 1110 1101 1100 1011 1010 1001 1000 0111 0110 0101 0100 0011 0010 0001",
        "1110 1011 1000 1010 0010 0101 1101 0110 0011 1111 1100 1001 0111 0100 0001")]
    [TestCase("00001 00010 00100 01000 10000 00011 00101 01001 10001 11111 11110 11101 11011 10111 01111",
        "00101 00011 00001 00100 01001 01111 01000 00010 10000 11111 11110 11011 10001 11101 10111")]
    public void KeyOrder_MatchesJavaHashMap(string inserted, string javaOrder)
    {
        IReadOnlyList<string> order = OncologyAnalyzer.JavaStringHashMapOrder.KeyOrder(inserted.Split(' '));

        Assert.That(order, Is.EqualTo(javaOrder.Split(' ')));
    }

    // 50 16-character keys (binary of the ints below), the first 12 colliding in the low 6 bits of the spread hash:
    // the bin is resized at 16/32 and treeified at 64 (red-black TreeNode bin, root moved to the front), then the
    // table grows to 128 and the tree bin is split. java.util.HashMap keySet() order as indices into the insertion list.
    [Test]
    public void KeyOrder_TreeifiedAndSplitBin_MatchesJavaHashMap()
    {
        int[] ints =
        {
            33877, 3735, 39229, 54940, 4351, 60470, 21437, 2643, 13065, 63991, 9591, 23649, 2503, 5725, 28859, 27835,
            4660, 16036, 6042, 36669, 28253, 3938, 55032, 37624, 8250, 63055, 14869, 41972, 41759, 38794, 63075, 4126,
            38400, 38966, 26402, 3303, 64979, 14727, 3106, 37041, 57132, 8875, 19286, 27898, 9610, 35979, 35538, 43380,
            5985, 26805,
        };
        int[] javaOrder =
        {
            20, 26, 16, 7, 3, 49, 4, 47, 6, 9, 10, 41, 14, 21, 19, 15, 18, 38, 34, 24, 30, 22, 33, 23, 29, 1, 0, 2, 11,
            46, 5, 8, 48, 31, 35, 39, 12, 17, 36, 13, 28, 45, 40, 44, 27, 42, 37, 43, 25, 32,
        };
        string[] keys = ints.Select(i => Convert.ToString(i, 2).PadLeft(16, '0')).ToArray();

        IReadOnlyList<string> order = OncologyAnalyzer.JavaStringHashMapOrder.KeyOrder(keys);

        Assert.That(order, Is.EqualTo(javaOrder.Select(i => keys[i])));
    }

    [Test]
    public void JavaStringHash_MatchesJavaStringHashCode()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.JavaStringHashMapOrder.JavaStringHash("001"), Is.EqualTo(47665));
            Assert.That(OncologyAnalyzer.JavaStringHashMapOrder.JavaStringHash("111"), Is.EqualTo(48657));
            // "1111111111" overflows int: Java "1111111111".hashCode() = 2075774624 (31-polynomial mod 2^32).
            Assert.That(OncologyAnalyzer.JavaStringHashMapOrder.JavaStringHash("1111111111"), Is.EqualTo(2075774624));
        });
    }

    #endregion

    #region Equal-score trees: LICHeE picks the first in HashMap-ordered enumeration

    // t02042 (ε 0.1, 3 samples): profiles in first appearance 101, 110, 100; Java HashMap order 110, 100, 101.
    // lichee.jar: 3 valid trees, all error 0; top tree root->{2,4}, 4->1, 4->3. With first-appearance order the
    // first enumerated tree would be 2->3 instead.
    [Test]
    public void EqualScoreTie_T02042_FollowsHashMapProfileOrder()
    {
        var clusters = new[] { C(1, 0.25, 0, 0.8), C(2, 0.2, 0.68, 0), C(3, 0.2, 0, 0), C(4, 0.73, 0, 0.9) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters, 0.1);

        Assert.Multiple(() =>
        {
            Assert.That(p.ValidTreeCount, Is.EqualTo(3));
            Assert.That(p.ErrorScore, Is.EqualTo(0.0));
            Assert.That(Parents(p, 4), Is.EqualTo(new int?[] { 4, -1, 4, -1 }));
        });
    }

    // t02235 (ε 0.2): profiles 010, 101, 110, 100; HashMap order 110, 100, 101, 010. lichee.jar: 3 trees tied at
    // error 0.10000000000000009; top root->{1,3,4}, 4->2, 4->5 (first-appearance order: 3->5).
    [Test]
    public void EqualScoreTie_T02235_FollowsHashMapProfileOrder()
    {
        var clusters = new[]
        {
            C(1, 0, 0.74, 0), C(2, 0.4, 0, 0.6), C(3, 0.3, 0.21, 0), C(4, 0.8, 0, 0.8), C(5, 0.3, 0, 0),
        };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters, 0.2);

        Assert.Multiple(() =>
        {
            Assert.That(p.ValidTreeCount, Is.EqualTo(3));
            Assert.That(p.ErrorScore, Is.EqualTo(0.10000000000000009));
            Assert.That(Parents(p, 5), Is.EqualTo(new int?[] { -1, 4, -1, -1, 4 }));
        });
    }

    // t00750 (ε 0, 5 samples): profiles 00110, 00111, 11001, 00001; HashMap order 11001, 00110, 00111, 00001.
    // lichee.jar: 2 trees tied at 0; top root->{2,3,5}, 5->1, 3->4 (first-appearance order: 2->4).
    [Test]
    public void EqualScoreTie_T00750_FollowsHashMapProfileOrder()
    {
        var clusters = new[]
        {
            C(1, 0, 0, 0.23, 0.3, 0), C(2, 0, 0, 0.44, 0.34, 0.9), C(3, 0.1, 0.3, 0, 0, 0.1),
            C(4, 0, 0, 0, 0, 0.07), C(5, 0, 0, 0.3, 0.6, 0),
        };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters, 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(p.ValidTreeCount, Is.EqualTo(2));
            Assert.That(p.ErrorScore, Is.EqualTo(0.0));
            Assert.That(Parents(p, 5), Is.EqualTo(new int?[] { 5, -1, -1, 3, -1 }));
        });
    }

    // t01561 (ε 0): profiles 101, 001, 011; HashMap order 011, 001, 101. lichee.jar: 2 trees tied at 0; top
    // root->{1,4}, 1->2, 2->3 (first-appearance order: 4->3).
    [Test]
    public void EqualScoreTie_T01561_FollowsHashMapProfileOrder()
    {
        var clusters = new[] { C(1, 0.46, 0, 0.6), C(2, 0, 0, 0.4), C(3, 0, 0, 0.3), C(4, 0, 0.45, 0.3) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters, 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(p.ValidTreeCount, Is.EqualTo(2));
            Assert.That(p.ErrorScore, Is.EqualTo(0.0));
            Assert.That(Parents(p, 4), Is.EqualTo(new int?[] { -1, 1, 2, -1 }));
        });
    }

    // t01396 (ε 0.05, 5 samples): a single valid tree, but node ids follow the HashMap order (11010, 00101, 11110,
    // 01101), which changes the summation order of computeErrorScore: lichee.jar error 0.03741657386773935
    // (first-appearance node ids give 0.03741657386773934).
    [Test]
    public void NodeIdOrder_T01396_ErrorScoreBitIdenticalToLichee()
    {
        var clusters = new[]
        {
            C(1, 0, 0, 0.5, 0, 0.1), C(2, 0, 0.19, 0.48, 0, 0.19), C(3, 0.41, 0.42, 0, 0.7, 0), C(4, 0.4, 0.42, 0.04, 0.67, 0),
        };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters, 0.05);

        Assert.Multiple(() =>
        {
            Assert.That(p.ValidTreeCount, Is.EqualTo(1));
            Assert.That(p.ErrorScore, Is.EqualTo(0.03741657386773935));
            Assert.That(Parents(p, 4), Is.EqualTo(new int?[] { 2, -1, 4, -1 }));
        });
    }

    #endregion
}
