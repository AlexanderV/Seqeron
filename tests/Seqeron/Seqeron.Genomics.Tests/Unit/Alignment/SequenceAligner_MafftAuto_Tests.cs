using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c part 4, B09 F287): <see cref="MafftOptions.Auto"/> = <c>mafft --auto</c> of MAFFT 7.505 (the script chooses L-INS-i for fewer than 100 sequences of fewer than 3 000 letters, two iterations of the local pair
/// library for fewer than 200 of fewer than 1 000, FFT-NS-i with two iterations for fewer than 500, FFT-NS-2 for fewer than 20 000) against the installed MAFFT 7.505 (Debian 7.505-1; scripts/alignment/generate_mafft_auto_oracle.py,
/// mafft_auto_oracle.tsv.gz). Every committed run must give the program's rows and its standard error (SHA-256 of the whole text, the full text for every 25th run).
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftAuto_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgressive_Tests.Row>> Oracle = new(() => SequenceAligner_MafftProgressive_Tests.LoadOracle("mafft_auto_oracle.tsv.gz"));

    public static IEnumerable<TestCaseData> Strategies() => Oracle.Value.Select(r => r.Strategy).Distinct().Select(s => new TestCaseData(s).SetArgDisplayNames(s));

    [Test]
    public void Oracle_HasEveryBranchOfTheRule()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(2000));
            foreach (var g in rows.GroupBy(r => r.Strategy))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(rows.Count(r => r.Seqs.Length >= 100 && r.Seqs.Length < 200), Is.GreaterThanOrEqualTo(10), "L-INS-i with 2 iterations");
            Assert.That(rows.Count(r => r.Seqs.Length >= 200 && r.Seqs.Length < 500), Is.GreaterThanOrEqualTo(8), "FFT-NS-i");
            Assert.That(rows.Count(r => r.Seqs.Length >= 500), Is.GreaterThanOrEqualTo(5), "FFT-NS-2");
            Assert.That(rows.Count(r => r.Seqs.Length < 100), Is.GreaterThanOrEqualTo(2000), "L-INS-i");
            Assert.That(rows.Count(r => r.Rc != 0), Is.GreaterThanOrEqualTo(5), "the plain command stops on letters outside the alphabet");
            Assert.That(rows.Count(r => r.Stderr.Length > 0), Is.GreaterThanOrEqualTo(80), "full texts of the standard error");
        });
    }

    [TestCaseSource(nameof(Strategies))]
    public void RowsAndStandardError_EqualMafft(string strategy) => SequenceAligner_MafftProgressive_Tests.CompareStrategy(Oracle.Value, strategy);

    [Test]
    public void Rule_PicksTheStrategyOfTheScript()
    {
        string[] Family(int n) => Enumerable.Range(0, n).Select(i => "MKVLAAGIVGLLLAQ".Remove(i % 12, 1)).ToArray();
        var o = new MafftOptions { AnySymbol = false, Auto = true, Quiet = false };
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.RunMafft(Family(5), o).StandardError, Does.Contain("L-INS-i (Probably most accurate, very slow)\n Iterative refinement method (<16) with LOCAL pairwise alignment information"));
            Assert.That(SequenceAligner.RunMafft(Family(100), o).StandardError, Does.Contain("L-INS-i (Probably most accurate, very slow)\n Iterative refinement method (<2) with LOCAL"));
            Assert.That(SequenceAligner.RunMafft(Family(200), o).StandardError, Does.Contain("FFT-NS-i (Standard)\n Iterative refinement method (max. 2 iterations)"));
            Assert.That(SequenceAligner.RunMafft(Family(500), o).StandardError, Does.Contain("FFT-NS-2 (Fast but rough)"));
        });
    }
}
