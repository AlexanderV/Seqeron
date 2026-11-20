using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-23, B09 F127): ProbCons 1.12 port (<see cref="SequenceAligner.ProbConsRun"/>, <see cref="ProbConsParameters"/>, <see cref="ProbConsProgram"/>).
/// The expected values were produced by the installed Debian probcons 1.12-14 (the committed oracle covers 407 command lines in <see cref="SequenceAligner_ProbConsOracle_Tests"/>).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ProbCons_Tests
{
    private static readonly string[] Names = ["a one", "b", "c", "d"];
    private static readonly string[] Seqs = ["MKVLAAGIVGLLLAQPSTAEEKV", "MKVLAGIVGLLAQPSTAEEKVLL", "MKILAAGIAGLLLQPSAEEKV", "MRVLAAGLVGLLLAQSTAEEK"];

    [Test]
    public void Default_RowsAndTreeAreTheProgramsOutput()
    {
        var r = SequenceAligner.ProbConsRun(Seqs, new ProbConsOptions { Names = Names });
        Assert.Multiple(() =>
        {
            Assert.That(r.Alignment.AlignedSequences, Is.EqualTo(new[]
            {
                "MKVLAAGIVGLLLAQPSTAEEK--V", "MKVL-AGIVGLLA-QPSTAEEKVLL", "MKILAAGIAGLLL-QPS-AEEK--V", "MRVLAAGLVGLLL-AQSTAEEK---",
            }), "probcons t1.fa");
            Assert.That(r.GuideTree, Is.EqualTo("((a one (c d)) b)"));
            Assert.That(r.Order, Is.EqualTo(new[] { 0, 1, 2, 3 }));
            Assert.That(r.Log, Does.Contain("Alignment tree: ((a one (c d)) b)\n\n\n"));
        });
    }

    [Test]
    public void NoConsistencyNoRefinement_AnnotationIsTheProgramsColumnScores()
    {
        var r = SequenceAligner.ProbConsRun(Seqs, new ProbConsOptions { ConsistencyReps = 0, IterativeRefinementReps = 0, Annotate = true });
        Assert.That(r.Annotation!.Take(25), Is.EqualTo(new[] { 99, 99, 98, 95, 99, 75, 96, 99, 99, 99, 92, 86, 79, 68, 65, 0, 83, 70, 99, 62, 64, 65, 0, 0, 59 }), "probcons -c 0 -ir 0 -annot");
    }

    [Test]
    public void AlignmentOrder_ReturnsTheRowsInTheProgramsAlignmentOrder()
    {
        var r = SequenceAligner.ProbConsRun(Seqs, new ProbConsOptions { Names = Names, AlignmentOrder = true });
        Assert.Multiple(() =>
        {
            Assert.That(r.Order, Is.EqualTo(new[] { 0, 2, 3, 1 }), "probcons -a: a, c, d, b");
            Assert.That(r.Alignment.AlignedSequences[0], Is.EqualTo("MKVLAAGIVGLLLAQPSTAEEK--V"));
            Assert.That(r.Alignment.AlignedSequences[3], Is.EqualTo("MKVL-AGIVGLLA-QPSTAEEKVLL"));
        });
    }

    [Test]
    public void Viterbi_PairsAreTheProgramsPairwiseAlignments()
    {
        var res = SequenceAligner.ProbConsPairs(Seqs, viterbi: true, new ProbConsOptions { Names = Names });
        Assert.Multiple(() =>
        {
            Assert.That(res.Pairs, Has.Count.EqualTo(6));
            Assert.That(res.Pairs[0].FirstRow, Is.EqualTo("MKVLAAGIVGLLLAQPSTAEEKV--"), "a one-b.fasta");
            Assert.That(res.Pairs[0].SecondRow, Is.EqualTo("MKVL-AGIVGL-LAQPSTAEEKVLL"));
            Assert.That(res.Pairs[2].SecondRow, Is.EqualTo("MRVLAAGLVGLLLAQ-STAEEK-"), "a one-d.fasta");
        });
    }

    [Test]
    public void GlibcRandom_IsTheCLibrarysRandAtSeedOne()
    {
        var r = new SequenceAligner.GlibcRandom();
        Assert.That(Enumerable.Range(0, 5).Select(_ => r.Next()), Is.EqualTo(new[] { 1804289383, 846930886, 1681692777, 1714636915, 1957747793 }));
    }

    [Test]
    public void CFormatGeneral_IsPrintfG()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CFormat.General(0.01993141696f, 10), Is.EqualTo("0.01993141696"));
            Assert.That(CFormat.General(8.615339902e-05f, 10), Is.EqualTo("8.615339902e-05"));
            Assert.That(CFormat.General(0.6814756989f, 10), Is.EqualTo("0.6814756989"));
            Assert.That(CFormat.General(0.0, 10), Is.EqualTo("0"));
            Assert.That(CFormat.General(1234567890123.0, 10), Is.EqualTo("1.234567890e+12".Replace(".234567890", ".23456789")));
            Assert.That(CFormat.General(100.0, 10), Is.EqualTo("100"));
            Assert.That(CFormat.General(0.0001, 10), Is.EqualTo("0.0001"));
            Assert.That(CFormat.General(0.00001, 10), Is.EqualTo("1e-05"));
            Assert.That(CFormat.General(-2.5, 3), Is.EqualTo("-2.5"));
            Assert.That(CFormat.General(double.NaN, 10), Is.EqualTo("-nan"), "double.NaN has the sign bit set on x86, as the default NaN of 0 / 0");
        });
    }

    [Test]
    public void Parameters_DefaultSetEqualsTheTCoffeeEmbeddedModel_AndTheFileRoundTrips()
    {
        var p = ProbConsParameters.Default;
        var back = ProbConsParameters.Parse(p.ToFileText());
        Assert.Multiple(() =>
        {
            Assert.That(p.InitDistrib, Is.EqualTo(new[] { 0.6814756989f, 8.615339902e-05f, 8.615339902e-05f, 0.1591759622f, 0.1591759622f }));
            Assert.That(p.GapOpen, Is.EqualTo(new[] { 0.0119511066f, 0.0119511066f, 0.008008334786f, 0.008008334786f }));
            Assert.That(p.GapExtend, Is.EqualTo(new[] { 0.3965826333f, 0.3965826333f, 0.8988758326f, 0.8988758326f }));
            Assert.That(p.EmitPairs['a', 'R'], Is.EqualTo(0.00244502f), "emitPairsDefault[1][0] for both cases and both orders");
            Assert.That(p.EmitPairs['r', 'A'], Is.EqualTo(0.00244502f));
            Assert.That(p.EmitPairs['B', 'A'], Is.EqualTo(1e-10f), "a letter outside the alphabet");
            Assert.That(p.EmitSingle['v'], Is.EqualTo(0.07343426f));
            Assert.That(p.EmitSingle['X'], Is.EqualTo(1e-5f));
            Assert.That(back.InitDistrib[3], Is.EqualTo(0.1591759622f).Within(1e-9f));
            Assert.That(back.EmitPairs['W', 'w'], Is.EqualTo(p.EmitPairs['W', 'w']).Within(1e-9f));
            Assert.That(back.Alphabet, Is.EqualTo("ARNDCQEGHILKMFPSTWYV"));
        });
    }

    [Test]
    public void Training_OfOneFileGivesTheProgramsParameterFileLayout()
    {
        var res = ProbConsProgram.Run(["-t", "params.out", "a.fa"], new Dictionary<string, string> { ["a.fa"] = ">a\nMKVLAAGIVGLLLAQPSTAEEKV\n>b\nMKVLAGIVGLLAQPSTAEEKVLL\n>c\nMKILAAGIAGLLLQPSAEEKV\n" });
        var lines = res.Files["params.out"].Split('\n');
        Assert.Multiple(() =>
        {
            Assert.That(res.ExitCode, Is.Zero);
            Assert.That(lines[3], Is.EqualTo("ARNDCQEGHILKMFPSTWYV"));
            Assert.That(lines, Has.Length.EqualTo(3 + 1 + 20 + 1 + 1));
            Assert.That(lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries), Has.Length.EqualTo(5));
            Assert.That(res.Stderr, Does.Contain("Trained parameter set:"));
        });
    }

    [Test]
    public void Validation_RefusesWhatTheProgramRefuses()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => SequenceAligner.ProbConsRun(Seqs, new ProbConsOptions { ConsistencyReps = 6 }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ProbConsRun(Seqs, new ProbConsOptions { IterativeRefinementReps = 1001 }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ProbConsRun(Seqs, new ProbConsOptions { PreTrainingReps = 21 }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ProbConsRun(Seqs, new ProbConsOptions { GapOpenPenalty = 1 }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ProbConsRun(Seqs, new ProbConsOptions { Cutoff = 1.5f }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ProbConsRun(["MKV", "MK1"]), "Unknown character encountered");
            Assert.Throws<ArgumentException>(() => SequenceAligner.ProbConsRun(["MKV", ""]));
            Assert.Throws<ArgumentException>(() => SequenceAligner.ProbConsRun([]));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.ProbConsRun(null!));
            Assert.Throws<ArgumentException>(() => ProbConsParameters.Parse("0.5 0.5\n"));
        });
    }

    [Test]
    public void SingleSequence_IsReturnedAsItIs()
    {
        var r = SequenceAligner.ProbConsRun(["MKVLAAGIV"]);
        Assert.Multiple(() =>
        {
            Assert.That(r.Alignment.AlignedSequences, Is.EqualTo(new[] { "MKVLAAGIV" }));
            Assert.That(r.GuideTree, Is.EqualTo("S1"));
        });
    }

    [Test]
    public void Program_ClustalOutputHasTheConservationLine()
    {
        var res = ProbConsProgram.Run(["-clustalw", "-ir", "0", "in.fa"], new Dictionary<string, string> { ["in.fa"] = ">a one\nMKVLAAGIVG\n>b\nMKVLAAGIVG\n" });
        Assert.Multiple(() =>
        {
            Assert.That(res.Stdout, Is.EqualTo("PROBCONS version 1.12 multiple sequence alignment\n\na    MKVLAAGIVG\nb    MKVLAAGIVG\n     **********\n"));
            Assert.That(res.Stderr, Does.StartWith("\nPROBCONS version 1.12 - align multiple protein sequences and print to standard output\nWritten by Chuong Do\n"));
        });
    }

    [Test]
    public void Program_NoFileGivesTheDebianPatchMessage_AndNoArgumentsTheUsage()
    {
        var none = ProbConsProgram.Run(["-c", "1"]);
        var usage = ProbConsProgram.Run([]);
        Assert.Multiple(() =>
        {
            Assert.That(none.ExitCode, Is.Zero);
            Assert.That(none.Stderr, Does.EndWith("Zero sequences delivered to function DoAlign()\n"));
            Assert.That(usage.ExitCode, Is.EqualTo(1));
            Assert.That(usage.Stderr, Does.Contain("probcons [OPTION]... [MFAFILE]..."));
            Assert.That(ProbConsProgram.Run(["-e", "x.fa"]).Stderr, Does.EndWith("ERROR: Training emissions (-e) requires training (-t)\n"));
        });
    }

    [Test]
    public void NumberReading_FollowsStrtolBaseZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ProbConsProgram.Parse(["-c", "0x2"]).Options.ConsistencyReps, Is.EqualTo(2), "hexadecimal");
            Assert.That(ProbConsProgram.Parse(["-c", "3abc"]).Options.ConsistencyReps, Is.EqualTo(3), "trailing characters are ignored");
            Assert.That(ProbConsProgram.Parse(["-ir", "010"]).Options.IterativeRefinementReps, Is.EqualTo(8), "octal");
            Assert.Throws<ArgumentException>(() => ProbConsProgram.Parse(["-c", "x"]));
            Assert.Throws<ArgumentException>(() => ProbConsProgram.Parse(["-c", "6"]));
        });
    }
}
