using NUnit.Framework;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9e / A3-9f (B09 F208, F209): the perl programs of MUMmer 3.23 - <c>nucmer2xfig</c> (TestData/Alignment/nucmer2xfig_oracle.jsonl.gz, scripts/alignment/generate_nucmer2xfig_oracle.py) and the file naming and
/// standard error of <c>nucmer -p|--prefix</c> (TestData/Alignment/nucmer_prefix_oracle.jsonl.gz, scripts/alignment/generate_nucmer_prefix_oracle.py).
/// </summary>
[TestFixture]
public class MummerPerl_Oracle_Tests
{
    [Test]
    public void Nucmer2xfig_OutputEqualsTheScript_AllRows()
    {
        var (failures, c) = MummerPerlOracleCheck.CheckNucmer2Xfig();
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(c["equal"], Is.GreaterThanOrEqualTo(200));
        Assert.That(c["labels"], Is.GreaterThan(100));
        Assert.That(c["malformed"], Is.GreaterThan(5));
    }

    [Test]
    public void NucmerPrefix_FilesAndStandardError_EqualTheScript_AllRows()
    {
        var (failures, c) = MummerPerlOracleCheck.CheckNucmerPrefix();
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(c["equal"], Is.GreaterThanOrEqualTo(100));
        Assert.That(c["customPrefix"], Is.GreaterThan(40));
        Assert.That(c["nodelta"], Is.GreaterThan(10));
        Assert.That(c["coords"], Is.GreaterThan(20));
        Assert.That(c["subdir"], Is.GreaterThan(5));
    }
}
