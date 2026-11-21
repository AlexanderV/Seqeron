using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9e (B09 F204): the <c>mummer</c> command of MUMmer 3.23 on random reference / query files (TestData/Alignment/mummer_cli_oracle.jsonl.gz, scripts/alignment/generate_mummer_cli_oracle.py): every matching mode,
/// <c>-n -l -b -r -s -c -F -L</c>, several query files, multi-record references, empty records, odd characters, CR LF line ends and the runs on which the program fails.
/// </summary>
[TestFixture]
public class MummerCli_Oracle_Tests
{
    [Test]
    public void Mummer_OutputEqualsTheProgram_AllRuns()
    {
        var (failures, c) = MummerCliOracleCheck.CheckAll();
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(c["rows"], Is.GreaterThanOrEqualTo(300));
        Assert.That(c["errors"], Is.GreaterThan(10));
        Assert.That(c["Mum"], Is.GreaterThan(50));
        Assert.That(c["MumReference"], Is.GreaterThan(50));
        Assert.That(c["MaxMatch"], Is.GreaterThan(50));
        Assert.That(c["nonEmpty"], Is.GreaterThan(200));
        Assert.That(c["multiRef"], Is.GreaterThan(50));
        Assert.That(c["strings"], Is.GreaterThan(50));
        Assert.That(c["reverse"], Is.GreaterThan(100));
        Assert.That(c["multiQuery"], Is.GreaterThan(50));
    }

    [Test]
    public void Mummer_OptionConflicts_AreRefused()
    {
        const string fa = ">a\nacgtacgtacgtacgtacgtacgtacgt\n";
        Assert.Throws<ArgumentException>(() => MummerProgram.Run(fa, [fa], new MummerOptions { Both = true, OnlyReverse = true }));
        Assert.Throws<ArgumentException>(() => MummerProgram.Run(fa, [fa], new MummerOptions { ShowReversePositions = true }));
        Assert.Throws<ArgumentException>(() => MummerProgram.Run(fa, [fa], new MummerOptions { MinLength = 0 }));
        Assert.Throws<ArgumentException>(() => MummerProgram.Run(fa, [], null));
        Assert.Throws<FormatException>(() => MummerProgram.Run(string.Empty, [fa], null));
    }
}
