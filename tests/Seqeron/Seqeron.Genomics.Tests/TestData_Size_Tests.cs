using System.IO.Compression;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Audit R14-20 (B09 F299): the test data is stored compressed. TestData was 274 MB (the test assembly 270 MB, <c>biopython_coord_oracle.tsv</c> alone 68 MB of text); every text oracle above 1 MiB is now a <c>.gz</c>
/// (<c>scripts/alignment/compress_test_data.py</c>) that the tests open through <see cref="TestDataResource"/> by the name of the plain file. This contract keeps it so: a plain text file above 1 MiB is a failure,
/// and the helper opens both forms with the same text.
/// </summary>
[TestFixture]
public class TestData_Size_Tests
{
    private static IEnumerable<string> DataDirectories()
    {
        string tests = Path.GetDirectoryName(Path.GetDirectoryName(DataDirectory()))!;
        return Directory.EnumerateDirectories(tests).Select(d => Path.Combine(d, "TestData")).Where(Directory.Exists);
    }

    private static string DataDirectory()
    {
        var dir = new DirectoryInfo(TestContext.CurrentContext.TestDirectory);
        while (dir is not null && !Directory.Exists(Path.Combine(dir.FullName, "TestData")) && !File.Exists(Path.Combine(dir.FullName, "Seqeron.sln")))
            dir = dir.Parent;
        string? root = dir?.FullName;
        while (root is not null && !File.Exists(Path.Combine(root, "Seqeron.sln")))
            root = Path.GetDirectoryName(root);
        return Path.Combine(root ?? throw new InvalidOperationException("Seqeron.sln not found."), "tests", "Seqeron", "Seqeron.Genomics.Tests", "TestData");
    }

    [Test]
    public void NoPlainTextOracleAboveOneMebibyte_IsCommitted()
    {
        string[] text = [".tsv", ".csv", ".jsonl", ".json", ".txt", ".fa", ".fasta"];
        // every test project (the tests of the separately licensed assemblies keep their own oracle data, B09 WP79)
        var big = DataDirectories().SelectMany(dir => Directory.EnumerateFiles(dir, "*", SearchOption.AllDirectories))
            .Where(f => text.Contains(Path.GetExtension(f).ToLowerInvariant()) && new FileInfo(f).Length > 1 << 20)
            .Select(f => $"{Path.GetRelativePath(Path.GetDirectoryName(Path.GetDirectoryName(DataDirectory()))!, f)} ({new FileInfo(f).Length / 1048576} MB)")
            .ToList();
        Assert.That(big, Is.Empty, "gzip it: python3 scripts/alignment/compress_test_data.py");
    }

    [Test]
    public void TestDataResource_OpensTheGzipForm_ByThePlainName_AndReturnsTheSameText()
    {
        // biopython_coord_oracle.tsv is embedded as .gz only: the plain name opens it, the .gz name opens the raw bytes
        const string plain = "Seqeron.Genomics.Tests.TestData.Alignment.biopython_coord_oracle.tsv";
        using var viaHelper = TestDataResource.Open(plain);
        Assert.That(viaHelper, Is.Not.Null);
        using var reader = new StreamReader(viaHelper!);
        string header = reader.ReadLine()!;
        Assert.That(header, Does.Contain('\t'));
        using var raw = typeof(TestDataResource).Assembly.GetManifestResourceStream(plain + ".gz");
        Assert.That(raw, Is.Not.Null);
        using var gz = new GZipStream(raw!, CompressionMode.Decompress);
        using var direct = new StreamReader(gz);
        Assert.That(direct.ReadLine(), Is.EqualTo(header));
        Assert.That(typeof(TestDataResource).Assembly.GetManifestResourceStream(plain), Is.Null, "the plain form is not embedded");
        Assert.That(TestDataResource.Open("Seqeron.Genomics.Tests.TestData.Alignment.no_such_oracle.tsv"), Is.Null);
    }
}
