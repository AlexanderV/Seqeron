using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9e2 (B09 F210): <c>combineMUMs</c> and the script <c>run-mummer3</c> of MUMmer 3.23 (TestData/Alignment/combinemums_oracle.jsonl.gz, scripts/alignment/generate_combinemums_oracle.py).
/// </summary>
[TestFixture]
public class MummerCombineMums_Oracle_Tests
{
    [Test]
    public void CombineMums_AndRunMummer3_EqualTheInstalledPrograms_AllRows()
    {
        var (failures, c) = MummerCombineMumsOracleCheck.CheckAll();
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(c["cmums"], Is.GreaterThanOrEqualTo(200));
        Assert.That(c["runmummer3"], Is.GreaterThanOrEqualTo(200));
    }
}
