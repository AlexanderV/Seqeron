using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c part 3, B09 F286): <see cref="SequenceAligner.MultipleAlignMafft"/> / <see cref="SequenceAligner.RunMafft"/> with a pair library and <c>--maxiterate</c> (G-INS-i, L-INS-i, E-INS-i of MAFFT 7.505:
/// <c>tbfast</c> writes the files hat3 and hat2, <c>dvtditr -l</c> refines every branch of the tree with <c>Falign_localhom</c> / <c>partA__align</c> and the importance of the pair-library blocks) against the installed MAFFT 7.505
/// (Debian 7.505-1; scripts/alignment/generate_mafft_pairlib_iterative_oracle.py, mafft_pairlib_iterative_oracle.tsv.gz). Every committed run must give the program's rows and its standard error (SHA-256 of the whole text,
/// the full text for every 25th run).
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)] // the strategy groups are independent (every call has its own state)
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftPairLibraryIterative_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgressive_Tests.Row>> Oracle = new(() => SequenceAligner_MafftProgressive_Tests.LoadOracle("mafft_pairlib_iterative_oracle.tsv.gz"));

    public static IEnumerable<TestCaseData> Strategies() => Oracle.Value.Select(r => r.Strategy).Distinct().Select(s => new TestCaseData(s).SetArgDisplayNames(s));

    [Test]
    public void Oracle_HasEveryStrategy()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(9000));
            foreach (var g in rows.GroupBy(r => r.Strategy))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            foreach (var tag in new[] { "G-", "L-", "E-" })
                Assert.That(rows.Count(r => r.Strategy.StartsWith(tag, StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(3000), tag);
            Assert.That(rows.Select(r => r.Strategy).Distinct().Count(), Is.GreaterThanOrEqualTo(60));
            Assert.That(rows.Count(r => r.Rc != 0), Is.GreaterThanOrEqualTo(10), "the plain command stops on letters outside the alphabet");
            Assert.That(rows.Count(r => r.Stderr.Length > 0), Is.GreaterThanOrEqualTo(150), "full texts of the standard error");
            Assert.That(rows.Max(r => r.Seqs.Length), Is.GreaterThanOrEqualTo(14));
            Assert.That(rows.Count(r => r.Seqs.Length == 2), Is.GreaterThanOrEqualTo(600), "two sequences");
        });
    }

    [TestCaseSource(nameof(Strategies))]
    public void RowsAndStandardError_EqualMafft(string strategy) => SequenceAligner_MafftProgressive_Tests.CompareStrategy(Oracle.Value, strategy);

    [Test]
    public void Threads_WithTheIterativeStep_AreRefused_ZeroThreadsAreTheSingleThreadProgram()
    {
        var seqs = new[] { "MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ", "MKVLAAGIVGLLAQ" };
        var o = new MafftOptions { AnySymbol = false, PairLibrary = MafftPairLibrary.Global, MaxIterate = 2 };
        Assert.Multiple(() =>
        {
            Assert.Throws<NotSupportedException>(() => SequenceAligner.MultipleAlignMafft(seqs, o with { Threads = 4 }), "non-deterministic in the program");
            Assert.That(SequenceAligner.MultipleAlignMafft(seqs, o with { Threads = 4, ThreadsIterate = 0 }).AlignedSequences, Is.EqualTo(SequenceAligner.MultipleAlignMafft(seqs, o).AlignedSequences));
        });
    }

    [Test]
    public void Strategy_Text_NamesThePairLibraryAndTheIterations()
    {
        var seqs = new[] { "MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ" };
        var o = new MafftOptions { AnySymbol = false, MaxIterate = 5 };
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.RunMafft(seqs, o with { PairLibrary = MafftPairLibrary.Global }).StandardError, Does.Contain("Strategy:\n G-INS-i (Suitable for sequences of similar lengths, very slow)\n Iterative refinement method (<5) with GLOBAL pairwise alignment information\n"));
            Assert.That(SequenceAligner.RunMafft(seqs, o with { PairLibrary = MafftPairLibrary.Local }).StandardError, Does.Contain("L-INS-i (Probably most accurate, very slow)\n Iterative refinement method (<5) with LOCAL pairwise alignment information"));
            Assert.That(SequenceAligner.RunMafft(seqs, o with { PairLibrary = MafftPairLibrary.Genaf, MaxIterate = 1000 }).StandardError, Does.Contain("E-INS-i (Suitable for sequences with long unalignable regions, very slow)\n Iterative refinement method (<16) with LOCAL pairwise alignment with generalized affine gap costs (Altschul 1998)"), "the script lets 16 iterations run at most");
        });
    }
}
