using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c part 1, B09 F284): <see cref="SequenceAligner.MultipleAlignMafft"/> / <see cref="SequenceAligner.RunMafft"/> with <c>--maxiterate</c> (FFT-NS-i / NW-NS-i of MAFFT 7.505:
/// <c>disttbfast</c>, then <c>dndpre</c> and <c>dvtditr</c>, the tree-dependent iterative refinement) against the installed MAFFT 7.505 (Debian 7.505-1; scripts/alignment/generate_mafft_iterative_oracle.py,
/// mafft_iterative_oracle.tsv.gz). Every committed run must give the program's rows and its standard error (SHA-256 of the whole text, the full text for every 25th run).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftIterative_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgressive_Tests.Row>> Oracle = new(() => SequenceAligner_MafftProgressive_Tests.LoadOracle("mafft_iterative_oracle.tsv.gz"));

    public static IEnumerable<TestCaseData> Strategies() => Oracle.Value.Select(r => r.Strategy).Distinct().Select(s => new TestCaseData(s).SetArgDisplayNames(s));

    [Test]
    public void Oracle_HasEveryStrategy()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(5000));
            foreach (var g in rows.GroupBy(r => r.Strategy))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(rows.Select(r => r.Strategy).Distinct().Count(), Is.GreaterThanOrEqualTo(20));
            Assert.That(rows.Count(r => r.Rc != 0), Is.GreaterThanOrEqualTo(5), "the plain command stops on letters outside the alphabet");
            Assert.That(rows.Count(r => r.Stderr.Length > 0), Is.GreaterThanOrEqualTo(150), "full texts of the standard error");
            Assert.That(rows.Max(r => r.Seqs.Length), Is.GreaterThanOrEqualTo(40));
            Assert.That(rows.Count(r => r.Seqs.Length == 2), Is.GreaterThanOrEqualTo(200), "two sequences");
        });
    }

    [TestCaseSource(nameof(Strategies))]
    public void RowsAndStandardError_EqualMafft(string strategy) => SequenceAligner_MafftProgressive_Tests.CompareStrategy(Oracle.Value, strategy);

    [Test]
    public void Threads_WithTheIterativeStep_AreRefused_ZeroThreadsAreTheSingleThreadProgram()
    {
        var seqs = new[] { "MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ", "MKVLAAGIVGLLAQ" };
        var plain = new MafftOptions { AnySymbol = false, MaxIterate = 2 };
        Assert.Multiple(() =>
        {
            Assert.Throws<NotSupportedException>(() => SequenceAligner.MultipleAlignMafft(seqs, plain with { Threads = 4 }), "non-deterministic in the program");
            Assert.That(SequenceAligner.MultipleAlignMafft(seqs, plain with { Threads = 4, ThreadsIterate = 0 }).AlignedSequences, Is.EqualTo(SequenceAligner.MultipleAlignMafft(seqs, plain).AlignedSequences));
        });
    }

    [Test]
    public void Strategy_Text_NamesTheIterations()
    {
        var seqs = new[] { "MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ" };
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.RunMafft(seqs, new MafftOptions { AnySymbol = false, MaxIterate = 2 }).StandardError, Does.Contain("Strategy:\n FFT-NS-i (Standard)\n Iterative refinement method (max. 2 iterations)\n"));
            Assert.That(SequenceAligner.RunMafft(seqs, new MafftOptions { AnySymbol = false, MaxIterate = 1000 }).StandardError, Does.Contain("FFT-NS-i (Accurate but slow)\n Iterative refinement method (max. 16 iterations)"), "the script lets 16 iterations run at most");
            Assert.That(SequenceAligner.RunMafft(seqs, new MafftOptions { AnySymbol = false, MaxIterate = 2, Quiet = true }).StandardError, Is.Empty);
        });
    }
}
