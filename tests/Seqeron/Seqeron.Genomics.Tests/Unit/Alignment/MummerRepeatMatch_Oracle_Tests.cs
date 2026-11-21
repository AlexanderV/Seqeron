using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9e3 (B09 F211): <c>repeat-match</c> and the script <c>exact-tandems</c> of MUMmer 3.23 (TestData/Alignment/repeatmatch_oracle.jsonl.gz, scripts/alignment/generate_repeatmatch_oracle.py): the installed script
/// (broken: <c>awk -f $</c>) and the intended pipeline run by hand. An empty file or a file without a header line leaves the program's sequence buffer uninitialised (valgrind: use of uninitialised values; the output changes under
/// MALLOC_PERTURB_): recorded, refused by the port, not compared.
/// </summary>
[TestFixture]
public class MummerRepeatMatch_Oracle_Tests
{
    [Test]
    public void RepeatMatch_AndExactTandems_EqualTheInstalledPrograms_AllRows()
    {
        var (failures, c) = MummerRepeatMatchOracleCheck.CheckAll();
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(c["rm"], Is.GreaterThanOrEqualTo(200));
        Assert.That(c["rmStops"], Is.GreaterThan(10));
        Assert.That(c["rmUndefined"], Is.EqualTo(2));
        Assert.That(c["rmAbort"], Is.GreaterThanOrEqualTo(3));
        Assert.That(c["et"], Is.GreaterThanOrEqualTo(200));
        Assert.That(c["etx"], Is.GreaterThanOrEqualTo(200));
    }

    [Test]
    public void RepeatMatch_ListsTheMatchesInTheOrderOfTheSuffixTreeOfTheProgram()
    {
        // repeat-match -n 12 t.fa (MUMmer 3.23): the order is that of the program's own tree (here: by the first position, the copies in the order the tree holds them)
        const string fasta = ">tandem demo\nacactgctcactccaaccccggcccctgagtccgaggagaggatcacagtctggatcaca\ngtctggatcacagtctggatcacagtctggatcacagtctgggtgcttcagagtatgtat\naccactgggtggatcaggatcaggatcaaggatacggcggagggcacg\n";
        var r = MummerRepeatMatch.Run(fasta, new MummerRepeatMatchOptions { MinMatchLength = 12 });
        Assert.That(r.Text, Is.EqualTo("Long Exact Matches:\n   Start1     Start2    Length\n       41         53        50\n       41         65        38\n       41         77        26\n       41         89        14\n      131        137        12\n"));
        Assert.That(r.Messages, Is.EqualTo("Genome Length = 168   Used 248 internal nodes\n"));
        Assert.That(r.ExitCode, Is.Zero);
        // the reverse-complement strand: 'r' after the second position; a palindrome stops the program
        Assert.That(MummerRepeatMatch.Run(">pal\nacgtagcta\n", new MummerRepeatMatchOptions { MinMatchLength = 3 }).Text, Is.EqualTo("Long Exact Matches:\n   Start1     Start2    Length\n        4          9r        6\n        1          4r        4\n"));
        var palindrome = MummerRepeatMatch.Run(">p\nacgtacgt\n", new MummerRepeatMatchOptions { MinMatchLength = 3 });
        Assert.That((palindrome.Text, palindrome.ExitCode), Is.EqualTo(("*** Genome is exact palindrome ***\n", 255)));
        Assert.Throws<ArgumentException>(() => MummerRepeatMatch.Run(fasta, new MummerRepeatMatchOptions { MinMatchLength = 0 }));
        Assert.Throws<FormatException>(() => MummerRepeatMatch.Run(string.Empty));
        Assert.Throws<FormatException>(() => MummerRepeatMatch.Run(">" + new string('x', 500) + "\nacgtacgtacgtacgtac\n"));
    }

    [Test]
    public void ExactTandems_TheInstalledScriptIsBroken_TheIntendedPipelineWritesTheTable()
    {
        const string fasta = ">tandem demo\nacactgctcactccaaccccggcccctgagtccgaggagaggatcacagtctggatcaca\ngtctggatcacagtctggatcacagtctggatcacagtctgggtgcttcagagtatgtat\naccactgggtggatcaggatcaggatcaaggatacggcggagggcacg\n";
        var intended = MummerExactTandems.Run(fasta, 6);
        Assert.That(intended.Text, Is.EqualTo("   Start   Extent  UnitLen     Copies\n      41       62       12        5.2\n     131       18        6        3.0\n"));
        var installed = MummerExactTandems.RunInstalled(fasta, 6);
        Assert.That(installed.Text, Is.EqualTo("Finding matches\nTandem repeats\n"));
        Assert.That(installed.ExitCode, Is.EqualTo(2));
        Assert.That(installed.Log, Is.EqualTo("Genome Length = 168   Used 129 internal nodes\n" + MummerExactTandems.InstalledAwkError));
        Assert.That(installed.LeftBehind, Is.EqualTo(new[] { "$$.tmp.matches" }));
        Assert.That(MummerExactTandems.Run(fasta, 0).Text, Is.EqualTo("   Start   Extent  UnitLen     Copies\n"));
    }
}
