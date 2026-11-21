using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>A3-9e (B09 F207): <c>delta2maf</c> of the Debian MUMmer 3.23 package on nucmer deltas of multi-record, both-strand inputs (TestData/Alignment/delta2maf_oracle.jsonl.gz, scripts/alignment/generate_delta2maf_oracle.py).</summary>
[TestFixture]
public class MummerDelta2Maf_Oracle_Tests
{
    [Test]
    public void Delta2Maf_OutputEqualsTheProgram_AllRows()
    {
        var (failures, c) = MummerDelta2MafOracleCheck.CheckAll();
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(c["equal"], Is.GreaterThanOrEqualTo(200));
        Assert.That(c["errors"], Is.GreaterThan(20));
        Assert.That(c["withIds"], Is.GreaterThan(50));
        Assert.That(c["sorted"], Is.GreaterThan(50));
        Assert.That(c["widths"], Is.GreaterThan(20));
        Assert.That(c["manyAlignments"], Is.GreaterThan(50));
    }

    [Test]
    public void Delta2Maf_PromerDeltaFile_IsRefused_AndTheHeaderWarningsAreTheProgram()
    {
        const string fasta = ">r\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n";
        const string query = ">q\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n";
        Assert.Throws<NotSupportedException>(() => MummerDelta2Maf.Run("r.fa q.fa\nPROMER\n>r q 40 40\n1 9 1 9 0 0 0\n0\n", fasta, query));
        const string delta = "r.fa q.fa\nNUCMER\n>r q 40 40\n1 40 1 40 0 0 0\n0\n";
        var r = MummerDelta2Maf.Run(delta, fasta, query, options: new MummerDelta2MafOptions { ScreenWidth = 5, MatrixType = 9, SortByQuery = true, SortByReference = true });
        Assert.That(r.Messages, Is.EqualTo("WARNING: invalid screen width 100000, using default\nWARNING: invalid matrix type 9, using default\nWARNING: both -r and -q were passed, -q ignored\n"));
    }
}
