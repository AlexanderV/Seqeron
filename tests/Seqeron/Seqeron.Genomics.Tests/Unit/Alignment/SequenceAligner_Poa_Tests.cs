using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-26, B09 F128): POA 2.0 port (<see cref="SequenceAligner.PoaRun"/>, <see cref="PoaScoreMatrix"/>, <see cref="PoaProgram"/>). The expected values are the output of the
/// installed Debian poa 2.0+20060928 (the committed oracle of 359 command lines is in <see cref="SequenceAligner_PoaOracle_Tests"/>).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_Poa_Tests
{
    private static readonly PoaSequence[] Seqs =
    {
        new("a", "one", "MKVLAAGIVGLLLAQPSTAEEKV"), new("b", null, "MKVLAGIVGLLAQPSTAEEKVLL"), new("c", null, "MKILAAGIAGLLLQPSAEEKV"), new("d", null, "MRVLAAGLVGLLLAQSTAEEK"),
    };

    [Test]
    public void Default_IterativeLocalAlignment_RowsAreThePrograms()
    {
        var r = SequenceAligner.PoaRun(Seqs);
        Assert.Multiple(() =>
        {
            Assert.That(r.Rows, Is.EqualTo(new[] { "MKVLAAGIVGLLLAQPSTAEEKV--", "MKVLA-GIVGLL-AQPSTAEEKVLL", "MKILAAGIAGLL-LQPS-AEEKV--", "MRVLAAGLVGLLLAQ-STAEEK---" }), "poa -read_fasta t1.fa");
            Assert.That(r.Names, Is.EqualTo(new[] { "a", "b", "c", "d" }));
            Assert.That(r.Order, Is.EqualTo(new[] { 0, 1, 2, 3 }));
            Assert.That(r.Pir, Does.StartWith(">a one\nMKVLAAGIVGLLLAQPSTAEEKV..\n>b untitled\nMKVLA.GIVGLL.AQPSTAEEKVLL\n"));
            Assert.That(r.Clustal, Does.StartWith("CLUSTAL W (1.74) multiple sequence alignment\n\n\na" + new string(' ', 35) + "MKVLAAGIVGLLLAQPSTAEEKV--\n"));
            Assert.That(r.Po, Does.StartWith("VERSION=LPO.0.1\nNAME=a\nTITLE=one\nLENGTH=30\nSOURCECOUNT=4\nSOURCENAME=a\nSOURCEINFO=23 0 1 -1 one\n"));
        });
    }

    [Test]
    public void HeaviestBundle_AddsConsensusSequences()
    {
        var r = SequenceAligner.PoaRun(Seqs, null, new PoaOptions { HeaviestBundle = true });
        Assert.Multiple(() =>
        {
            Assert.That(r.Names, Is.EqualTo(new[] { "a", "b", "c", "d", "CONSENS0", "CONSENS1" }));
            Assert.That(r.Rows[4], Is.EqualTo("MKVLAAGIVGLLLAQPSTAEEKVLL"), "poa -hb: the heaviest path");
            Assert.That(r.Rows[5], Is.EqualTo("MKILAAGIAGLL-LQPS-AEEKV--"));
            Assert.That(r.Pir, Does.Contain(">CONSENS0 consensus produced by heaviest_bundle, containing 3 seqs\n"));
            Assert.That(r.Order, Is.EqualTo(new[] { 0, 1, 2, 3, -1, -1 }));
        });
    }

    [Test]
    public void ProgressivePreserveOrder_PrintsThePermutation_AndKeepsInputOrder()
    {
        var r = SequenceAligner.PoaRun(Seqs, null, new PoaOptions { Global = true, Progressive = true, PreserveSequenceOrder = true, Verbose = true });
        Assert.Multiple(() =>
        {
            Assert.That(r.Stdout, Does.EndWith("0 2 3 1 \n"), "poa -do_global -do_progressive -preserve_seqorder -v");
            Assert.That(r.Stdout, Does.StartWith("X-Gap Penalties (Open, Aff1, Aff2; LTrunc, LDecay): 12 6 6 16 0\nX-Gap Penalties (0, 1, 2, ...): 12 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 ... \n"));
            Assert.That(r.Names, Is.EqualTo(new[] { "a", "b", "c", "d" }));
            Assert.That(r.Rows[1], Is.EqualTo("MKVLA-GIVGLLA-QPSTAEEKVLL"));
        });
    }

    [Test]
    public void ScoreMatrix_ParsesTheProgramsFormat_AndPrints()
    {
        var m = PoaScoreMatrix.Parse("# test\nGAP-PENALTIES=10 4 1\nGAP-TRUNCATION-LENGTH=6\n  a c\na 5 -4\nc -4 5\n");
        Assert.Multiple(() =>
        {
            Assert.That(m.Symbols, Is.EqualTo("ac"));
            Assert.That(m.MaxGapLength, Is.EqualTo(6));
            Assert.That(m.GapPenaltyX, Is.EqualTo(new[] { 10, 4, 4, 4, 4, 4, 1, 0 }), "open, extension x (trunc - 1), second extension, special state");
            Assert.That(m.Print("ac"), Is.EqualTo("   a  c\na  5 -4\nc -4  5\n"));
            Assert.Throws<ArgumentException>(() => PoaScoreMatrix.Parse("A C\nA 1 2\nB 3 4\n"), "a row symbol that is not in the symbol line");
            Assert.Throws<ArgumentException>(() => PoaScoreMatrix.Parse("A C\nA 1\nC 3 4\n"), "a missing score");
        });
    }

    [Test]
    public void Letters_OutsideTheMatrix_BecomeItsFirstSymbol_ToUpperKeepsTheLetters()
    {
        var seqs = new[] { new PoaSequence("x", null, "mkv"), new PoaSequence("y", null, "MKV") };
        var plain = SequenceAligner.PoaRun(seqs);
        var upper = SequenceAligner.PoaRun(seqs, null, new PoaOptions { Case = PoaCase.Upper });
        Assert.Multiple(() =>
        {
            Assert.That(plain.Rows, Is.EqualTo(new[] { "AAA---", "---MKV" }), "poa: lower-case letters are not symbols of blosum80.mat, limit_residues maps them to its first symbol A");
            Assert.That(upper.Rows, Is.EqualTo(new[] { "MKV", "MKV" }), "poa -toupper");
        });
    }

    [Test]
    public void Validation_RefusesWhatTheProgramRefuses()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => SequenceAligner.PoaRun([]));
            Assert.Throws<ArgumentException>(() => SequenceAligner.PoaRun([new PoaSequence("x", null, "   ")]));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.PoaRun(null!));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignPoa(["MKV", ""]));
            var unknown = SequenceAligner.PoaRun(Seqs, null, new PoaOptions { PairScores = "a nosuch 5\n" });
            Assert.That(unknown.Rows, Is.EqualTo(new[] { "MKVLAAGIVGLLLAQPSTAEEKV" }), "an unknown name in the pair scores: the program merges nothing and writes the first sequence alone, exit status 0");
        });
    }

    [Test]
    public void MultipleAlignPoa_ReturnsRowsInInputOrder()
    {
        var r = SequenceAligner.MultipleAlignPoa(Seqs.Select(s => s.Residues));
        Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "MKVLAAGIVGLLLAQPSTAEEKV--", "MKVLA-GIVGLL-AQPSTAEEKVLL", "MKILAAGIAGLL-LQPS-AEEKV--", "MRVLAAGLVGLLLAQ-STAEEK---" }));
    }

    [Test]
    public void Program_Usage_AndMsaInputIsRegistered()
    {
        Assert.Multiple(() =>
        {
            var usage = PoaProgram.Run([]);
            Assert.That(usage.ExitCode, Is.EqualTo(255), "exit (-1)");
            Assert.That(usage.Stderr, Does.StartWith("\nUsage: poa [OPTIONS] MATRIXFILE\n"));
            Assert.Throws<NotSupportedException>(() => PoaProgram.Run(["-read_msa", "x.msa", "m.mat"], new Dictionary<string, string> { ["m.mat"] = "A\nA 1\n" }));
        });
    }

    [Test]
    public void FastaReader_FollowsTheProgram()
    {
        var seqs = PoaProgram.ReadFasta("# comment\n>s1 first title\nMKV LAAG\n*ignored\nIVG\n>s2\nMKV\n>s3 \n", out var comment);
        Assert.Multiple(() =>
        {
            Assert.That(comment, Is.EqualTo(" comment"));
            Assert.That(seqs.Select(s => (s.Name, s.Title, s.Residues)), Is.EqualTo(new[] { ("s1", "first title", "MKV LAAGIVG"), ("s2", "untitled", "MKV") }), "a record without residues is dropped");
        });
    }
}
