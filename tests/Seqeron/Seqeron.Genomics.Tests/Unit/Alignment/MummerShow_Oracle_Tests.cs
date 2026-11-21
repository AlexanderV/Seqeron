using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9b (B09 F200 - F202): the delta programs of MUMmer 3.23 on delta files written by nucmer from rearranged, duplicated, inverted, mutated multi-contig inputs (TestData/Alignment/mummer_show_oracle.jsonl.gz,
/// scripts/alignment/generate_mummer_show_oracle.py; the programs ran with <c>time ()</c> fixed to 1700000000 so that <c>srand (time (NULL))</c> and the btab date are reproducible): <c>delta-filter</c>, <c>show-coords</c>,
/// <c>show-aligns</c>, <c>show-snps</c>, <c>show-diff</c>, <c>show-tiling</c> (with its -t / -u / -p files) and <c>dnadiff</c> (all files; the hash-ordered lines of the report and of .unref / .unqry are compared as sets of lines).
/// </summary>
[TestFixture]
public class MummerShow_Oracle_Tests
{
    [Test]
    public void DeltaPrograms_OutputEqualsTheProgram_AllRuns()
    {
        var (failures, cases, perProgram, f) = MummerShowOracleCheck.CheckAll("mummer_show_oracle.jsonl.gz");
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(cases, Is.GreaterThanOrEqualTo(150));
        foreach (string program in new[] { "delta-filter", "show-coords", "show-aligns", "show-snps", "show-diff", "show-tiling", "dnadiff" })
            Assert.That(perProgram.GetValueOrDefault(program), Is.GreaterThanOrEqualTo(200), program);
        Assert.That(f["tilingWithPath"], Is.GreaterThan(100));
        Assert.That(f["xmlLinks"], Is.GreaterThan(3));
        Assert.That(f["lasRuns"], Is.GreaterThan(3));
        Assert.That(f["snpSelect"], Is.GreaterThan(5));
        Assert.That(f["filterLis"], Is.GreaterThan(100));
    }
}
