using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Oracles;

[TestFixture]
public class CorpusManifestTests
{
    private const string Sample =
        "75942359851c45df9f7ba5327a83502b2d6c287f3fa3e7d4c4c9f70e05bfc750  data/dna_10m.fa\n" +
        "\n" +
        "# a comment\n" +
        "67e49c4bc074f073c89cc0f6bc48bdfcf6b0bd5bc79321bdd30c8c36d2b6d8ee *data/dna_10m.txt\r\n";

    [Test]
    public void Parse_ReadsHashAndPath_IgnoringBlankAndCommentLines()
    {
        var m = CorpusManifest.Parse(Sample);
        Assert.That(m.Entries, Has.Count.EqualTo(2));
        Assert.That(m.Entries[0].Path, Is.EqualTo("data/dna_10m.fa"));
        Assert.That(m.Entries[0].Sha256, Is.EqualTo("75942359851c45df9f7ba5327a83502b2d6c287f3fa3e7d4c4c9f70e05bfc750"));
        Assert.That(m.Entries[1].Path, Is.EqualTo("data/dna_10m.txt"));
        Assert.That(m.TryGet("data/dna_10m.txt", out var e), Is.True);
        Assert.That(e.Sha256, Does.StartWith("67e49c4b"));
        Assert.That(m.TryGet("data/none", out _), Is.False);
    }

    [TestCase("zz942359851c45df9f7ba5327a83502b2d6c287f3fa3e7d4c4c9f70e05bfc750  x")]
    [TestCase("75942359851c45df9f7ba5327a83502b2d6c287f3fa3e7d4c4c9f70e05bfc75  x")]
    [TestCase("75942359851c45df9f7ba5327a83502b2d6c287f3fa3e7d4c4c9f70e05bfc750x")]
    [TestCase("75942359851c45df9f7ba5327a83502b2d6c287f3fa3e7d4c4c9f70e05bfc750   ")]
    public void Parse_RejectsMalformedLines(string line)
    {
        Assert.Throws<FormatException>(() => CorpusManifest.Parse(line));
    }

    [Test]
    public void Parse_RejectsDuplicatePaths()
    {
        string h = new('a', 64);
        Assert.Throws<FormatException>(() => CorpusManifest.Parse($"{h}  p\n{h}  p\n"));
    }

    [Test]
    public void Verify_ComparesTheFileHash_AndFailsOnMissingFile()
    {
        string dir = Path.Combine(Path.GetTempPath(), "esa-manifest-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        try
        {
            File.WriteAllText(Path.Combine(dir, "abc.txt"), "abc");
            // SHA-256("abc")
            const string abc = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
            var m = CorpusManifest.Parse($"{abc}  abc.txt\n{new string('0', 64)}  bad.txt\n");
            File.WriteAllText(Path.Combine(dir, "bad.txt"), "abc");
            Assert.That(m.Verify("abc.txt", dir), Is.True);
            Assert.That(m.Verify("bad.txt", dir), Is.False);
            Assert.Throws<FileNotFoundException>(() => CorpusManifest.Parse($"{abc}  missing.txt").Verify("missing.txt", dir));
            Assert.Throws<KeyNotFoundException>(() => m.Verify("unlisted.txt", dir));
        }
        finally
        {
            Directory.Delete(dir, true);
        }
    }

    [Test]
    public void TheCommittedBenchManifest_ParsesAndListsTheGateInputs()
    {
        var m = CorpusManifest.LoadBenchManifest();
        Assert.That(m.Entries, Has.Count.GreaterThan(20));
        Assert.That(m.TryGet("data/dna_10m.txt", out var dna10), Is.True);
        Assert.That(dna10.Sha256, Is.EqualTo("67e49c4bc074f073c89cc0f6bc48bdfcf6b0bd5bc79321bdd30c8c36d2b6d8ee"));
        foreach (string name in new[] { "data/dna_1m.txt", "data/rep_10m.txt", "data/prot_1m.txt", "queries/dna_10m.memq100k.txt" })
            Assert.That(m.TryGet(name, out _), Is.True, name);
    }
}
