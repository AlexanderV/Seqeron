using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9e (B09 F206): <c>gaps</c>, <c>delta2blocks</c>, <c>mummer-annotate</c> and the script <c>run-mummer1</c> of MUMmer 3.23 (TestData/Alignment/mummer_gaps_oracle.jsonl.gz, scripts/alignment/generate_mummer_gaps_oracle.py).
/// Fourteen crafted annotate runs read memory in front of the program's buffer (valgrind: "Invalid read of size 1 ... 13 bytes before a block of size 10,000 alloc'd"; three of them also change under MALLOC_PERTURB_): the result of the program
/// is undefined there, they are recorded and refused by the port, not compared.
/// </summary>
[TestFixture]
public class MummerGaps_Oracle_Tests
{
    [Test]
    public void GapsDelta2blocksAnnotateRunMummer1_EqualTheInstalledPrograms_AllRows()
    {
        var (failures, c) = MummerGapsOracleCheck.CheckAll();
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(c["gaps"], Is.GreaterThanOrEqualTo(200));
        Assert.That(c["gapsCrash"], Is.EqualTo(1));
        Assert.That(c["delta2blocks"], Is.GreaterThanOrEqualTo(200));
        Assert.That(c["delta2blocksCrash"], Is.GreaterThan(5));
        Assert.That(c["annotate"], Is.GreaterThanOrEqualTo(200));
        Assert.That(c["annotateUndefined"], Is.GreaterThanOrEqualTo(10));
        Assert.That(c["runmummer1"], Is.GreaterThan(30));
        Assert.That(c["runmummer1Fail"], Is.GreaterThan(30));
    }

    [Test]
    public void Gaps_ReadsPartialTriples_AndRefusesWhatLoopsOrCrashesTheProgram()
    {
        Assert.Throws<FormatException>(() => MummerGaps.Run("g.fa", false, string.Empty));                 // SIGSEGV
        Assert.Throws<FormatException>(() => MummerGaps.Run("g.fa", false, "> q\n   1 2 3\n"));            // scanf never consumes the text: an endless loop
        Assert.That(MummerGaps.Run("g.fa", true, "3 9 20\n"), Is.EqualTo("> g.fa reverse Consistent matches\n       3        9     20    none      -      -\n> g.fa reverse Other matches\n"));
    }

    [Test]
    public void Annotate_AWindowBeforeTheFirstBase_IsUndefinedInTheProgram()
    {
        const string fasta = ">ref\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n";
        const string gaps = "> ref.fa  Consistent matches\n       1        1     30    none      -      -\n       5        5     30    none      1      1\n";
        Assert.Throws<InvalidOperationException>(() => MummerAnnotate.Run(gaps, fasta, new Dictionary<string, string> { ["ref.fa"] = fasta }));
        Assert.Throws<FormatException>(() => MummerAnnotate.Run(gaps, fasta, new Dictionary<string, string>()));
    }

    [Test]
    public void RunMummer1_WithoutTheReverseArgument_FailsLikeTheScript()
    {
        const string fasta = ">r\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n";
        Assert.Throws<FormatException>(() => MummerRunPrograms.RunMummer1(fasta, fasta, "r.fa", "p", reverse: false));
        Assert.Throws<FormatException>(() => MummerRunPrograms.RunMummer1(fasta, fasta, "r.fa", string.Empty));
    }
}
