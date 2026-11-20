using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>multiple_align_poa / poa_command: POA 2.0 (audit A2-26, B09 F128). The tools delegate to SequenceAligner.PoaRun / PoaProgram.Run; the oracle comparison lives in the library tests.</summary>
[TestFixture]
public class PoaTests
{
    private static readonly string[] Seqs = ["MKVLAAGIVGLLLAQPSTAEEKV", "MKVLAGIVGLLAQPSTAEEKVLL", "MKILAAGIAGLLLQPSAEEKV", "MRVLAAGLVGLLLAQSTAEEK"];

    [Test]
    public void MultipleAlignPoa_EqualsTheProgram_RowsAndPoFile()
    {
        var r = AlignmentTools.MultipleAlignPoa(Seqs, ["a", "b", "c", "d"]);
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "MKVLAAGIVGLLLAQPSTAEEKV--", "MKVLA-GIVGLL-AQPSTAEEKVLL", "MKILAAGIAGLL-LQPS-AEEKV--", "MRVLAAGLVGLLLAQ-STAEEK---" }), "poa -read_fasta t1.fa");
            Assert.That(r.Names, Is.EqualTo(new[] { "a", "b", "c", "d" }));
            Assert.That(r.Po, Does.StartWith("VERSION=LPO.0.1\nNAME=a\n"));
            Assert.That(r.Pir, Does.StartWith(">a untitled\n"));
            Assert.That(r.Clustal, Does.StartWith("CLUSTAL W (1.74) multiple sequence alignment\n"));
        });
    }

    [Test]
    public void MultipleAlignPoa_OptionsAndValidation()
    {
        var hb = AlignmentTools.MultipleAlignPoa(Seqs, heaviestBundle: true);
        Assert.Multiple(() =>
        {
            Assert.That(hb.Names, Has.Length.EqualTo(6));
            Assert.That(hb.Names[4], Is.EqualTo("CONSENS0"));
            Assert.That(hb.Order, Is.EqualTo(new[] { 0, 1, 2, 3, -1, -1 }));
            Assert.That(AlignmentTools.MultipleAlignPoa(["mkv", "MKV"], @case: "upper").AlignedSequences, Is.EqualTo(new[] { "MKV", "MKV" }));
            Assert.That(AlignmentTools.MultipleAlignPoa(Seqs, global: true, progressive: true, preserveSequenceOrder: true).Stdout, Is.EqualTo("0 2 3 1 \n"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignPoa(Seqs, @case: "mixed"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignPoa(Seqs, matrix: "A C\nA 1\nC 3 4\n"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignPoa([]));
            Assert.Throws<ArgumentNullException>(() => AlignmentTools.MultipleAlignPoa(null!));
        });
    }

    [Test]
    public void PoaCommand_RunsTheCommandLine_FilesAndExitCodes()
    {
        var files = new Dictionary<string, string>
        {
            ["in.fa"] = ">a one\n" + Seqs[0] + "\n>b\n" + Seqs[1] + "\n>c\n" + Seqs[2] + "\n>d\n" + Seqs[3] + "\n",
            ["m.mat"] = "# t\nGAP-PENALTIES=12 6 6\n A R\nA 7 -3\nR -3 9\n",
        };
        var res = AlignmentTools.PoaCommand(["-read_fasta", "in.fa", "-pir", "p.pir", "m.mat"], files);
        Assert.Multiple(() =>
        {
            Assert.That(res.ExitCode, Is.Zero);
            Assert.That(res.Files["p.pir"], Does.StartWith(">a one\n"));
            Assert.That(res.Stderr, Does.Contain("...Read 4 sequences from sequence file in.fa..."));
            Assert.That(AlignmentTools.PoaCommand(["-read_fasta", "in.fa", "missing.mat"], files).ExitCode, Is.EqualTo(1));
            Assert.That(AlignmentTools.PoaCommand([]).ExitCode, Is.EqualTo(255));
            Assert.Throws<ArgumentException>(() => AlignmentTools.PoaCommand(null!));
        });
    }
}
