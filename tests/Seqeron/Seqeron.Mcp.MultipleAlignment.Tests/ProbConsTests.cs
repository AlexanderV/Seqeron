// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.MultipleAlignment.Tools;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

/// <summary>multiple_align_probcons / probcons_command: ProbCons 1.12 (audit A2-23, B09 F127). The tools delegate to Seqeron.Genomics.TCoffee.ProbConsAligner.ProbConsRun / ProbConsProgram.Run; the oracle comparison lives in the library tests.</summary>
[TestFixture]
public class ProbConsTests
{
    private static readonly string[] Seqs = ["MKVLAAGIVGLLLAQPSTAEEKV", "MKVLAGIVGLLAQPSTAEEKVLL", "MKILAAGIAGLLLQPSAEEKV", "MRVLAAGLVGLLLAQSTAEEK"];

    [Test]
    public void MultipleAlignProbCons_EqualsTheProgram_TreeAndAnnotation()
    {
        var r = MultipleAlignmentTools.MultipleAlignProbCons(Seqs, ["a one", "b", "c", "d"], annotate: true);
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
        var noRefine = MultipleAlignmentTools.MultipleAlignProbCons(Seqs, consistency: 0, iterativeRefinement: 0, annotate: true);
        Assert.Multiple(() =>
        {
            Assert.That(noRefine.Annotation!.Take(5), Is.EqualTo(new[] { 99, 99, 98, 95, 99 }), "probcons -c 0 -ir 0 -annot");
            Assert.That(MultipleAlignmentTools.MultipleAlignProbCons(Seqs, alignmentOrder: true).Order, Is.EqualTo(new[] { 0, 2, 3, 1 }));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignProbCons(Seqs, consistency: 9));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignProbCons(Seqs, parameters: "0.5\n"));
            Assert.Throws<ArgumentNullException>(() => MultipleAlignmentTools.MultipleAlignProbCons(null!));
        });
    }

    [Test]
    public void ProbConsCommand_RunsTheCommandLine_ClustalOutputAndStderr()
    {
        var files = new Dictionary<string, string> { ["in.fa"] = ">a one\n" + Seqs[0] + "\n>b\n" + Seqs[1] + "\n>c\n" + Seqs[2] + "\n>d\n" + Seqs[3] + "\n" };
        var res = MultipleAlignmentTools.ProbConsCommand(["-clustalw", "in.fa"], files);
        Assert.Multiple(() =>
        {
            Assert.That(res.ExitCode, Is.Zero);
            Assert.That(res.Stdout, Does.StartWith("PROBCONS version 1.12 multiple sequence alignment\n\na    MKVLAAGIVGLLLAQPSTAEEK--V\n"));
            Assert.That(res.Stderr, Does.Contain("Alignment tree: ((a one (c d)) b)"));
            Assert.That(MultipleAlignmentTools.ProbConsCommand(["-c", "9", "in.fa"], files).ExitCode, Is.EqualTo(1));
            Assert.That(MultipleAlignmentTools.ProbConsCommand(["missing.fa"], files).Stderr, Does.Contain("Could not open file 'missing.fa' for reading."));
            Assert.That(MultipleAlignmentTools.ProbConsCommand(["-pairs", "in.fa"], files).Files.Keys, Has.Member("a one-b.fasta").And.Count.EqualTo(6));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ProbConsCommand(null!));
        });
    }
}
