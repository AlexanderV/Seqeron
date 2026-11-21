using System.Text;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-26c, B09 F232): known answers of SPOA 4.0.8 (<see cref="SpoaProgram"/>, <see cref="SpoaGraph"/>, <see cref="SpoaEngine"/>) on a three-read set; the values are the output of the installed Debian
/// <c>spoa</c> 4.0.8-1 (<c>spoa -r 4 u.fa</c>, <c>spoa -l 1 -r 1 u.fa</c>, <c>spoa -l 2 -r 2 u.fa</c>); the committed oracle (<see cref="SpoaProgram_Oracle_Tests"/>) holds the other command lines.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SpoaProgram_Tests
{
    private const string Reads = ">r1\nACGTACGTAC\n>r2\nACGTTCGTAC\n>r3\nACGACGTAC\n";

    private static SpoaRunResult Run(params string[] args) => SpoaProgram.Run(args, p => p == "u.fa" ? Encoding.Latin1.GetBytes(Reads) : null);

    [Test]
    public void Gfa_OfTheThreeReads_AsTheProgram()
    {
        const string expected = "H\tVN:Z:1.0\nS\t1\tA\tic:Z:true\nL\t1\t+\t2\t+\tOM\tew:f:6\tic:Z:true\nS\t2\tC\tic:Z:true\nL\t2\t+\t3\t+\tOM\tew:f:6\tic:Z:true\n"
            + "S\t3\tG\tic:Z:true\nL\t3\t+\t4\t+\tOM\tew:f:4\tic:Z:true\nL\t3\t+\t5\t+\tOM\tew:f:2\tic:Z:true\nS\t4\tT\tic:Z:true\nL\t4\t+\t5\t+\tOM\tew:f:2\tic:Z:true\nL\t4\t+\t11\t+\tOM\tew:f:2\n"
            + "S\t5\tA\tic:Z:true\nL\t5\t+\t6\t+\tOM\tew:f:4\tic:Z:true\nS\t6\tC\tic:Z:true\nL\t6\t+\t7\t+\tOM\tew:f:6\tic:Z:true\nS\t7\tG\tic:Z:true\nL\t7\t+\t8\t+\tOM\tew:f:6\tic:Z:true\n"
            + "S\t8\tT\tic:Z:true\nL\t8\t+\t9\t+\tOM\tew:f:6\tic:Z:true\nS\t9\tA\tic:Z:true\nL\t9\t+\t10\t+\tOM\tew:f:6\tic:Z:true\nS\t10\tC\tic:Z:true\nS\t11\tT\nL\t11\t+\t6\t+\tOM\tew:f:2\n"
            + "P\tr1\t1+,2+,3+,4+,5+,6+,7+,8+,9+,10+\t*\nP\tr2\t1+,2+,3+,4+,11+,6+,7+,8+,9+,10+\t*\nP\tr3\t1+,2+,3+,5+,6+,7+,8+,9+,10+\t*\nP\tConsensus\t1+,2+,3+,4+,5+,6+,7+,8+,9+,10+\t*\n";
        var r = Run("-r", "4", "u.fa");
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(0));
            Assert.That(r.StdOut, Is.EqualTo(expected));
        });
    }

    [Test]
    public void GlobalMsa_And_SemiGlobalMsaWithConsensus_AsTheProgram()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Run("-l", "1", "-r", "1", "u.fa").StdOut, Is.EqualTo(">r1\nACGTACGTAC\n>r2\nACGTTCGTAC\n>r3\nACG-ACGTAC\n"));
            Assert.That(Run("-l", "2", "-r", "2", "u.fa").StdOut, Is.EqualTo(">r1\nACGTACGTAC\n>r2\nACGTTCGTAC\n>r3\nACG-ACGTAC\n>Consensus\nACGTACGTAC\n"));
        });
    }

    [Test]
    public void Version_Help_AndErrors_AsTheProgram()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Run("--version").StdOut, Is.EqualTo("4.0.8\n"));
            Assert.That(Run("--version").ExitCode, Is.EqualTo(0));
            Assert.That(Run("-h").StdOut, Does.StartWith("usage: spoa [options ...] <sequences>\n"));
            Assert.That(Run("-x").StdErr, Is.EqualTo("spoa: invalid option -- 'x'\n"));
            Assert.That(Run("--bogus").StdErr, Is.EqualTo("spoa: unrecognized option '--bogus'\n"));
            Assert.That(Run("missing.fa").StdErr, Is.EqualTo("[bioparser::Parser::Create] error: unable to open file missing.fa\n"));
            Assert.That(Run("u.txt").StdErr, Does.StartWith("[spoa::CreateParser] error: file u.txt has unsupported format extension"));
            Assert.That(Run("-g", "1", "u.fa").StdErr, Is.EqualTo("[spoa::AlignmentEngine::Create] error: gap opening penalty must be non-positive!\n"));
        });
    }

    [Test]
    public void InvalidFasta_AbortsLikeTheCppRuntime()
    {
        var r = SpoaProgram.Run(["a.fa"], _ => Encoding.Latin1.GetBytes("ACGT\n"));
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(134));
            Assert.That(r.StdErr, Is.EqualTo("terminate called after throwing an instance of 'std::invalid_argument'\n  what():  [bioparser::FastaParser] error: invalid file format\n"));
        });
    }

    [Test]
    public void DotFile_IsWrittenToTheNamedPath()
    {
        var r = Run("-d", "g.dot", "u.fa");
        Assert.Multiple(() =>
        {
            Assert.That(r.StdOut, Is.EqualTo(">Consensus LN:i:10\nACGTACGTAC\n"));
            Assert.That(r.Files["g.dot"], Does.StartWith("digraph 3 {\n  graph [rankdir = LR]\n  0[label = \"0 - A\", style = filled, fillcolor = goldenrod1]\n"));
        });
    }
}
