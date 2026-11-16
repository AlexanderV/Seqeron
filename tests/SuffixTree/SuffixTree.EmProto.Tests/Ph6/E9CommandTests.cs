using System.Text;
using SuffixTree.EmProto.Ph6;

namespace SuffixTree.EmProto.Tests.Ph6;

/// <summary>The commands of E9 (<c>e9-bwt</c>, <c>e9-psi</c>, <c>e9-isa</c>, <c>e9-psi-isa</c>, <c>e9-oracle</c>, <c>e9-smoke</c>): exit codes, outputs, and the rows the measurement scripts read.</summary>
[TestFixture]
public class E9CommandTests
{
    private static Dictionary<string, string> Rows(string path, string phaseOrCase) => File.ReadAllLines(path)
        .Select(l => l.Split('\t'))
        .Where(f => f.Length >= 6 && f[0] == "R" && f[1] == "E9" && f[2] == phaseOrCase)
        .GroupBy(f => f[3])
        .ToDictionary(g => g.Key, g => g.Last()[4]);

    private static int Run(params string[] args) => Cli.Run(args);

    [Test]
    public void TheChain_ProducesTheOracleFiles_AndTheRowsOfTheProtocol()
    {
        Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(20000, 4, 3));
        using var dir = new Ph6Dir();
        string sa = c.WriteSa(dir), pk = c.WritePacked(dir), rows = dir.File("rows.tsv");
        string input = dir.File("in.txt");
        File.WriteAllBytes(input, Encoding.Latin1.GetBytes(c.Text));
        string bwt = dir.File("x.bwt"), psi = dir.File("x.psi"), isa = dir.File("x.isa"), psi2 = dir.File("x.psi2"), prefix = dir.File("ref");
        long n = c.N;

        Assert.That(Run("e9-oracle", "--sa", sa, "--out", prefix, "--rows", rows, "--case", "t"), Is.EqualTo(0));
        Assert.That(File.ReadAllBytes(prefix + ".psi"), Is.EqualTo(c.ExpectedPsiBytes(false)));
        Assert.That(Run("e9-bwt", "--text", pk, "--sa", sa, "--out", bwt, "--read-kb", "1", "--window-kb", "1", "--rows", rows, "--case", "t"), Is.EqualTo(0));
        Assert.That(File.ReadAllBytes(bwt), Is.EqualTo(c.ExpectedBwt()));
        Assert.That(Run("e9-psi", "--bwt", bwt, "--out", psi, "--window-kb", "1", "--symbols-per-pass", "3", "--read-kb", "1", "--verify", prefix + ".psi", "--rows", rows, "--case", "t"), Is.EqualTo(0));
        Assert.That(Run("e9-isa", "--sa", sa, "--out", isa, "--scratch", dir.File("scr"), "--scale-n", "3000", "--io-kb", "1", "--verify", prefix + ".isa", "--rows", rows, "--case", "t"), Is.EqualTo(0));
        Assert.That(Run("e9-psi-isa", "--sa", sa, "--isa", isa, "--out", psi2, "--verify", prefix + ".psi", "--rows", rows, "--case", "t"), Is.EqualTo(0));
        Assert.That(Run("e9-smoke", "--input", input, "--sa", sa, "--isa", isa, "--psi", psi, "--psi2", psi2, "--sample", "500", "--rows", rows, "--case", "t"), Is.EqualTo(0));

        Dictionary<string, string> r = Rows(rows, "t");
        Assert.That(r["psi_equals_ref"], Is.EqualTo("yes"));
        Assert.That(r["isa_equals_ref"], Is.EqualTo("yes"));
        Assert.That(r["psi_from_isa_equals_ref"], Is.EqualTo("yes"));
        Assert.That(r["smoke_ok"], Is.EqualTo("yes"));
        Assert.That(r["smoke_psi_stream_mismatches"], Is.EqualTo("0"));
        Assert.That(long.Parse(r["psi_io_counted_total"]), Is.EqualTo((2 * (n + 1)) + (4 * (n + 1))), "BWT read passes x (n + 1) + Psi mapped 4 (n + 1)");
        Assert.That(long.Parse(r["psi_io_explicit_read"]), Is.EqualTo(2 * (n + 1)), "sigma 4 in groups of 3: two passes over the BWT");
        Assert.That(long.Parse(r["psi_io_mapped_write"]), Is.EqualTo(4 * (n + 1)));
        Assert.That(long.Parse(r["bwt_io_explicit_read"]), Is.EqualTo(4 * (n + 1)));
        Assert.That(long.Parse(r["bwt_io_mapped_write"]), Is.EqualTo(n + 1));
        Assert.That(long.Parse(r["isa_io_explicit_read"]), Is.GreaterThan(4 * (n + 1)), "SA once, plus the tuples that went through the chunk files (1 KB buffers)");
        Assert.That(long.Parse(r["isa_io_mapped_write"]), Is.EqualTo(4 * (n + 1)));
        Assert.That(File.Exists(BwtMeta.MetaPath(bwt)), Is.True);
    }

    [Test]
    public void AWrongReference_ExitsWithOne_AndSaysNo()
    {
        Ph6Case c = Ph6Case.OfNaive("abracadabra");
        using var dir = new Ph6Dir();
        string sa = c.WriteSa(dir), pk = c.WritePacked(dir), rows = dir.File("rows.tsv");
        string bwt = dir.File("x.bwt");
        Assert.That(Run("e9-bwt", "--text", pk, "--sa", sa, "--out", bwt), Is.EqualTo(0));
        string wrong = dir.File("wrong.psi");
        File.WriteAllBytes(wrong, c.ExpectedPsiBytes(false).Select(b => (byte)(b ^ 1)).ToArray());
        Assert.That(Run("e9-psi", "--bwt", bwt, "--out", dir.File("x.psi"), "--verify", wrong, "--rows", rows, "--case", "w"), Is.EqualTo(1));
        Assert.That(Rows(rows, "w")["psi_equals_ref"], Is.EqualTo("NO"));
    }

    [Test]
    public void CountsOfPh0_AreCrossChecked()
    {
        Ph6Case c = Ph6Case.OfNaive("abracadabra");
        using var dir = new Ph6Dir();
        string sa = c.WriteSa(dir), pk = c.WritePacked(dir), rows = dir.File("rows.tsv");
        string bwt = dir.File("x.bwt");
        Assert.That(Run("e9-bwt", "--text", pk, "--sa", sa, "--out", bwt), Is.EqualTo(0));
        string good = dir.File("good.meta"), bad = dir.File("bad.meta"), none = dir.File("none.meta");
        File.WriteAllText(good, "n=11\ncounts=" + string.Join(',', c.Counts()) + "\n");
        File.WriteAllText(bad, "n=11\ncounts=" + string.Join(',', c.Counts().Reverse()) + "\n");
        File.WriteAllText(none, "n=11\n");
        Assert.That(Run("e9-psi", "--bwt", bwt, "--out", dir.File("a.psi"), "--pk-meta", good, "--rows", rows, "--case", "g"), Is.EqualTo(0));
        Assert.That(Rows(rows, "g")["counts_equal_ph0"], Is.EqualTo("yes"));
        Assert.That(Run("e9-psi", "--bwt", bwt, "--out", dir.File("b.psi"), "--pk-meta", bad, "--rows", rows, "--case", "b"), Is.EqualTo(1));
        Assert.That(Rows(rows, "b")["counts_equal_ph0"], Is.EqualTo("NO"));
        Assert.That(Run("e9-psi", "--bwt", bwt, "--out", dir.File("c.psi"), "--pk-meta", none, "--rows", rows, "--case", "n"), Is.EqualTo(0));
        Assert.That(Rows(rows, "n")["counts_equal_ph0"], Is.EqualTo("n/a"));
    }

    [Test]
    public void ASmokeWithAWrongPsi_ExitsWithOne()
    {
        Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(3000, 4, 8));
        using var dir = new Ph6Dir();
        string sa = c.WriteSa(dir), rows = dir.File("rows.tsv"), input = dir.File("in.txt"), prefix = dir.File("ref");
        File.WriteAllBytes(input, Encoding.Latin1.GetBytes(c.Text));
        Assert.That(Run("e9-oracle", "--sa", sa, "--out", prefix), Is.EqualTo(0));
        Assert.That(Run("e9-smoke", "--input", input, "--sa", sa, "--isa", prefix + ".isa", "--psi", prefix + ".psi", "--rows", rows, "--case", "ok"), Is.EqualTo(0));
        byte[] psi = File.ReadAllBytes(prefix + ".psi");
        (psi[400], psi[404]) = (psi[404], psi[400]);
        (psi[800], psi[808]) = (psi[808], psi[800]);
        File.WriteAllBytes(dir.File("bad.psi"), psi);
        Assert.That(Run("e9-smoke", "--input", input, "--sa", sa, "--isa", prefix + ".isa", "--psi", dir.File("bad.psi"), "--rows", rows, "--case", "bad"), Is.EqualTo(1));
        Assert.That(Rows(rows, "bad")["smoke_ok"], Is.EqualTo("NO"));
    }

    [Test]
    public void MissingArguments_ExitWithTwo()
    {
        Assert.That(Run("e9-psi", "--bwt", "x"), Is.EqualTo(2));
        Assert.That(Run("e9-bwt"), Is.EqualTo(2));
        Assert.That(Run("e9-isa", "--sa", "x", "--out", "y"), Is.EqualTo(2));
        Assert.That(Run("e9-psi-isa", "--sa", "x"), Is.EqualTo(2));
        Assert.That(Run("e9-oracle", "--sa", "x"), Is.EqualTo(2));
        Assert.That(Run("e9-smoke", "--input", "x"), Is.EqualTo(2));
    }
}
