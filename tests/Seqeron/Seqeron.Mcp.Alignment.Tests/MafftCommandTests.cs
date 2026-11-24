using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>mafft_command (B09 WP70, F315 ff.): the script mafft of MAFFT 7.505 as a command line. The tool delegates to MafftProgram.Run; the oracle comparison against the installed program lives in the library tests (SequenceAligner_MafftProgram_Tests).</summary>
[TestFixture]
public class MafftCommandTests
{
    private const string Fasta = ">a first\nMKVLAAGIVGLLLA\n>b\nMKVLAGIVGLLA\n>c\nMKVIAAGLVGLLLA\n";

    private static Dictionary<string, string> Files() => new() { ["in.fa"] = Fasta };

    [Test]
    public void MafftCommand_RunsTheCommandLine_FastaOutput_AndStandardError()
    {
        var plain = AlignmentTools.MafftCommand(["--quiet", "in.fa"], Files());
        var verbose = AlignmentTools.MafftCommand(["in.fa"], Files());
        Assert.Multiple(() =>
        {
            Assert.That(plain.ExitCode, Is.Zero);
            Assert.That(plain.Stdout, Is.EqualTo(">a first\nMKVLAAGIVGLLLA\n>b\nMKVL-AGIVGLLA-\n>c\nMKVIAAGLVGLLLA\n"));
            Assert.That(plain.Stderr, Is.Empty);
            Assert.That(verbose.Stdout, Is.EqualTo(plain.Stdout));
            Assert.That(verbose.Stderr, Does.Contain("Strategy:\n FFT-NS-2 (Fast but rough)"));
        });
    }

    [Test]
    public void MafftCommand_OutputFormats_AndGuideTreeFile()
    {
        var clustal = AlignmentTools.MafftCommand(["--quiet", "--clustalout", "in.fa"], Files());
        var phylip = AlignmentTools.MafftCommand(["--quiet", "--phylipout", "in.fa"], Files());
        var tree = AlignmentTools.MafftCommand(["--quiet", "--treeout", "--distout", "--out", "o.fa", "in.fa"], Files());
        Assert.Multiple(() =>
        {
            Assert.That(clustal.Stdout, Does.StartWith("CLUSTAL format alignment by MAFFT FFT-NS-2 (v7.505)\n\n\na               MKVLAAGIVGLLLA"));
            Assert.That(phylip.Stdout, Does.StartWith(" 3 14\na          MKVLAAGIVG LLLA"));
            Assert.That(tree.Files.Keys, Is.EquivalentTo(new[] { "o.fa", "in.fa.tree", "in.fa.hat2" }));
            Assert.That(tree.Stdout, Is.Empty, "--out writes the alignment to the file");
            Assert.That(tree.Files["in.fa.tree"], Does.Contain("1_a").And.EndWith(";\n"));
        });
    }

    [Test]
    public void MafftCommand_RefusalsOfTheScript_AreTheScriptsOwnMessages()
    {
        var unknown = AlignmentTools.MafftCommand(["--bogus", "in.fa"], Files());
        var missing = AlignmentTools.MafftCommand(["nofile.fa"], Files());
        var version = AlignmentTools.MafftCommand(["--version"]);
        Assert.Multiple(() =>
        {
            Assert.That(unknown.ExitCode, Is.EqualTo(1));
            Assert.That(unknown.Stderr, Does.StartWith("Unknown option:  --bogus\n").And.Contain("High speed:"));
            Assert.That(missing.ExitCode, Is.EqualTo(1));
            Assert.That(missing.Stderr, Does.StartWith("/usr/bin/mafft: Cannot open nofile.fa."));
            Assert.That(version.Stderr, Is.EqualTo("v7.505 (2022/Apr/10)\n"));
            Assert.That(AlignmentTools.MafftCommand(["--lastpair", "in.fa"], Files()).Stderr, Does.Contain("== Install LAST"));
            Assert.Throws<NotSupportedException>(() => AlignmentTools.MafftCommand(["--sp", "in.fa"], Files()), "--sp (D__align) is registered, not ported yet");
        });
    }

    [Test]
    public void MafftCommand_Validation()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => AlignmentTools.MafftCommand(null!));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MafftCommand([]));
            Assert.Throws<ArgumentNullException>(() => AlignmentTools.MafftCommand(["in.fa", null!]));
        });
    }
}
