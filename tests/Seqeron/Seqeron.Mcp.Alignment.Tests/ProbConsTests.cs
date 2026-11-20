using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>multiple_align_probcons / probcons_command: ProbCons 1.12 (audit A2-23, B09 F127). The tools delegate to SequenceAligner.ProbConsRun / ProbConsProgram.Run; the oracle comparison lives in the library tests.</summary>
[TestFixture]
public class ProbConsTests
{
    private static readonly string[] Seqs = ["MKVLAAGIVGLLLAQPSTAEEKV", "MKVLAGIVGLLAQPSTAEEKVLL", "MKILAAGIAGLLLQPSAEEKV", "MRVLAAGLVGLLLAQSTAEEK"];

    [Test]
    public void MultipleAlignProbCons_EqualsTheProgram_TreeAndAnnotation()
    {
        var r = AlignmentTools.MultipleAlignProbCons(Seqs, ["a one", "b", "c", "d"], annotate: true);
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "MKVLAAGIVGLLLAQPSTAEEK--V", "MKVL-AGIVGLLA-QPSTAEEKVLL", "MKILAAGIAGLLL-QPS-AEEK--V", "MRVLAAGLVGLLL-AQSTAEEK---" }), "probcons t1.fa");
            Assert.That(r.GuideTree, Is.EqualTo("((a one (c d)) b)"));
            Assert.That(r.Order, Is.EqualTo(new[] { 0, 1, 2, 3 }));
            Assert.That(r.Annotation, Has.Length.EqualTo(25));
        });
    }

    [Test]
    public void MultipleAlignProbCons_OptionsAndValidation()
    {
        var noRefine = AlignmentTools.MultipleAlignProbCons(Seqs, consistency: 0, iterativeRefinement: 0, annotate: true);
        Assert.Multiple(() =>
        {
            Assert.That(noRefine.Annotation!.Take(5), Is.EqualTo(new[] { 99, 99, 98, 95, 99 }), "probcons -c 0 -ir 0 -annot");
            Assert.That(AlignmentTools.MultipleAlignProbCons(Seqs, alignmentOrder: true).Order, Is.EqualTo(new[] { 0, 2, 3, 1 }));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignProbCons(Seqs, consistency: 9));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignProbCons(Seqs, parameters: "0.5\n"));
            Assert.Throws<ArgumentNullException>(() => AlignmentTools.MultipleAlignProbCons(null!));
        });
    }

    [Test]
    public void ProbConsCommand_RunsTheCommandLine_ClustalOutputAndStderr()
    {
        var files = new Dictionary<string, string> { ["in.fa"] = ">a one\n" + Seqs[0] + "\n>b\n" + Seqs[1] + "\n>c\n" + Seqs[2] + "\n>d\n" + Seqs[3] + "\n" };
        var res = AlignmentTools.ProbConsCommand(["-clustalw", "in.fa"], files);
        Assert.Multiple(() =>
        {
            Assert.That(res.ExitCode, Is.Zero);
            Assert.That(res.Stdout, Does.StartWith("PROBCONS version 1.12 multiple sequence alignment\n\na    MKVLAAGIVGLLLAQPSTAEEK--V\n"));
            Assert.That(res.Stderr, Does.Contain("Alignment tree: ((a one (c d)) b)"));
            Assert.That(AlignmentTools.ProbConsCommand(["-c", "9", "in.fa"], files).ExitCode, Is.EqualTo(1));
            Assert.That(AlignmentTools.ProbConsCommand(["missing.fa"], files).Stderr, Does.Contain("Could not open file 'missing.fa' for reading."));
            Assert.That(AlignmentTools.ProbConsCommand(["-pairs", "in.fa"], files).Files.Keys, Has.Member("a one-b.fasta").And.Count.EqualTo(6));
            Assert.Throws<ArgumentException>(() => AlignmentTools.ProbConsCommand(null!));
        });
    }
}
