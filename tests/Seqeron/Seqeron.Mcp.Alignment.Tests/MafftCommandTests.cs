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
    public void MafftCommand_GuideTreeOptions_AreAccepted_AndEqualMultipleAlignMafft()
    {
        // B09 WP72 (F326 ff.): --pileup / --randomchain / --large / --stepadd / --youngestlinkage / --treein / --memsave / --skipiterate / --core are ported; the tool multiple_align_mafft has the same parameters
        string[] seqs = ["MKVLAAGIVGLLLA", "MKVLAGIVGLLA", "MKVIAAGLVGLLLA", "MKVLAAGIVGLLLG"];
        string fasta = string.Concat(seqs.Select((x, i) => ">s" + i + "\n" + x + "\n"));
        var files = new Dictionary<string, string> { ["in.fa"] = fasta };
        Assert.Multiple(() =>
        {
            foreach (var (tree, option) in new[] { ("pileup", "--pileup"), ("randomchain", "--randomchain"), ("large", "--large"), ("stepadd", "--stepadd"), ("youngestlinkage", "--youngestlinkage") })
            {
                var viaCommand = AlignmentTools.MafftCommand(["--quiet", option, "in.fa"], files);
                var viaTool = AlignmentTools.MultipleAlignMafft(seqs, guideTree: tree, anySymbol: false, quiet: true, retree: tree is "pileup" or "randomchain" ? 1 : 2);
                var rows = viaCommand.Stdout.Split('>', StringSplitOptions.RemoveEmptyEntries).Select(x => string.Concat(x.Split('\n').Skip(1))).ToArray();
                Assert.That(viaCommand.ExitCode, Is.Zero, option);
                Assert.That(viaTool.AlignedSequences, Is.EqualTo(rows), option);
            }
            var tl = AlignmentTools.MultipleAlignMafft(seqs, guideTree: "treein", guideTreeText: "1 2 0.1 0.2\n1 3 0.1 0.3\n1 4 0.2 0.1\n", anySymbol: false, quiet: true);
            Assert.That(tl.AlignedSequences, Has.Length.EqualTo(4));
            var ms = AlignmentTools.MultipleAlignMafft(seqs, memoryMode: "memsave", anySymbol: false, quiet: true);
            var msCommand = AlignmentTools.MafftCommand(["--quiet", "--memsave", "in.fa"], files);
            Assert.That(ms.AlignedSequences, Is.EqualTo(msCommand.Stdout.Split('>', StringSplitOptions.RemoveEmptyEntries).Select(x => string.Concat(x.Split('\n').Skip(1))).ToArray()));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMafft(seqs, guideTree: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMafft(seqs, memoryMode: "bogus"));
            var skip = AlignmentTools.MafftCommand(["--quiet", "--maxiterate", "2", "--skipiterate", "0.1", "in.fa"], files);
            Assert.That(skip.ExitCode, Is.Zero);
            var core = AlignmentTools.MafftCommand(["--core", "--maxiterate", "2", "--corewin", "5", "in.fa"], files);
            Assert.That(core.Stderr, Does.Contain("Loading 'hat2' ... done."));
            Assert.That(core.Stdout, Does.StartWith(">=s0\n"));
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
