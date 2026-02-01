using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class FastaWriteTests
{
    [Test]
    public void FastaWrite_Schema_ValidatesCorrectly()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            Assert.DoesNotThrow(() => ParsersTools.FastaWrite(tempFile, "seq1", "ATGCATGC"));
            Assert.Throws<ArgumentException>(() => ParsersTools.FastaWrite("", "seq1", "ATGC"));
            Assert.Throws<ArgumentException>(() => ParsersTools.FastaWrite(null!, "seq1", "ATGC"));
            Assert.Throws<ArgumentException>(() => ParsersTools.FastaWrite(tempFile, "", "ATGC"));
            Assert.Throws<ArgumentException>(() => ParsersTools.FastaWrite(tempFile, "seq1", ""));
            Assert.Throws<ArgumentException>(() => ParsersTools.FastaWrite(tempFile, "seq1", "ATGC", null, 0));
        }
        finally
        {
            if (File.Exists(tempFile)) File.Delete(tempFile);
        }
    }

    [Test]
    public void FastaWrite_Binding_InvokesSuccessfully()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            var result = ParsersTools.FastaWrite(tempFile, "seq1", "ATGCATGCATGC", "Human gene");

            Assert.That(result.FilePath, Is.EqualTo(tempFile));
            Assert.That(result.EntriesWritten, Is.EqualTo(1));
            Assert.That(result.TotalBases, Is.EqualTo(12));

            // Verify file content
            var content = File.ReadAllText(tempFile);
            Assert.That(content, Does.StartWith(">seq1 Human gene"));
            Assert.That(content, Does.Contain("ATGCATGCATGC"));
        }
        finally
        {
            if (File.Exists(tempFile)) File.Delete(tempFile);
        }
    }
}
