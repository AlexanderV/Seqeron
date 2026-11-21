using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>MCP <c>spoa_command</c> (B09 A2-26c, F232): the output of the installed spoa 4.0.8 on the same command line.</summary>
[TestFixture]
public class SpoaTests
{
    private static readonly Dictionary<string, string> Files = new() { ["u.fa"] = ">r1\nACGTACGTAC\n>r2\nACGTTCGTAC\n>r3\nACGACGTAC\n" };

    [Test]
    public void Msa_Gfa_And_Version_AsTheProgram()
    {
        var msa = AlignmentTools.SpoaCommand(["-l", "1", "-r", "1", "u.fa"], Files);
        var version = AlignmentTools.SpoaCommand(["--version"]);
        var dot = AlignmentTools.SpoaCommand(["-d", "g.dot", "u.fa"], Files);
        Assert.Multiple(() =>
        {
            Assert.That(msa.ExitCode, Is.EqualTo(0));
            Assert.That(msa.Stdout, Is.EqualTo(">r1\nACGTACGTAC\n>r2\nACGTTCGTAC\n>r3\nACG-ACGTAC\n"));
            Assert.That(version.Stdout, Is.EqualTo("4.0.8\n"));
            Assert.That(dot.Stdout, Is.EqualTo(">Consensus LN:i:10\nACGTACGTAC\n"));
            Assert.That(dot.Files["g.dot"], Does.StartWith("digraph 3 {"));
        });
    }

    [Test]
    public void GzipInput_And_Errors()
    {
        using var ms = new MemoryStream();
        using (var gz = new System.IO.Compression.GZipStream(ms, System.IO.Compression.CompressionLevel.Optimal, true))
        {
            var bytes = System.Text.Encoding.ASCII.GetBytes(Files["u.fa"]);
            gz.Write(bytes, 0, bytes.Length);
        }

        var r = AlignmentTools.SpoaCommand(["-r", "0", "u.fa.gz"], null, new Dictionary<string, string> { ["u.fa.gz"] = Convert.ToBase64String(ms.ToArray()) });
        var bad = AlignmentTools.SpoaCommand(["-x"]);
        Assert.Multiple(() =>
        {
            Assert.That(r.Stdout, Is.EqualTo(">Consensus LN:i:10\nACGTACGTAC\n"));
            Assert.That(bad.ExitCode, Is.EqualTo(1));
            Assert.That(bad.Stderr, Is.EqualTo("spoa: invalid option -- 'x'\n"));
        });
    }
}
