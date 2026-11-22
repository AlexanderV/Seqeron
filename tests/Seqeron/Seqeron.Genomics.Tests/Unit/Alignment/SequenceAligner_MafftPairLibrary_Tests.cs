using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c part 2, B09 F285): <see cref="SequenceAligner.MultipleAlignMafft"/> / <see cref="SequenceAligner.RunMafft"/> with a pair library (<c>--globalpair</c> G-INS-1, <c>--localpair</c> L-INS-1,
/// <c>--genafpair</c> E-INS-1 of MAFFT 7.505: <c>tbfast</c> with <c>pairlocalalign</c>: all pairs aligned by <c>Galign11</c> / <c>Lalign11</c> / <c>genalign11</c>, the guide tree from their distances, the importance of
/// the aligned blocks as constraints of <c>A__align</c>) against the installed MAFFT 7.505 (Debian 7.505-1; scripts/alignment/generate_mafft_pairlib_oracle.py, mafft_pairlib_oracle.tsv.gz). Every committed run
/// must give the program's rows and its standard error (SHA-256 of the whole text, the full text for every 25th run).
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)] // 57 strategy groups of up to 520 runs each: the runs are independent (every call has its own state)
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftPairLibrary_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgressive_Tests.Row>> Oracle = new(() => SequenceAligner_MafftProgressive_Tests.LoadOracle("mafft_pairlib_oracle.tsv.gz"));

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
            Assert.That(rows.Select(r => r.Strategy).Distinct().Count(), Is.GreaterThanOrEqualTo(57), "19 groups for each of G-INS-1, L-INS-1, E-INS-1");
            foreach (var tag in new[] { "G-", "L-", "E-" })
                Assert.That(rows.Count(r => r.Strategy.StartsWith(tag, StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(3000), tag);
            Assert.That(rows.Count(r => r.Rc != 0), Is.GreaterThanOrEqualTo(10), "the plain command stops on letters outside the alphabet");
            Assert.That(rows.Count(r => r.Stderr.Length > 0), Is.GreaterThanOrEqualTo(150), "full texts of the standard error");
            Assert.That(rows.Max(r => r.Seqs.Length), Is.GreaterThanOrEqualTo(21));
            Assert.That(rows.Count(r => r.Seqs.Length == 2), Is.GreaterThanOrEqualTo(600), "two sequences");
        });
    }

    [TestCaseSource(nameof(Strategies))]
    public void RowsAndStandardError_EqualMafft(string strategy) => SequenceAligner_MafftProgressive_Tests.CompareStrategy(Oracle.Value, strategy);

    // The degenerate case found by the oracle: with --fmodel a nucleotide that does not occur gives a matrix whose offset-free scores are positive while the (negative) distances of the pair library shift
    // the matrix of every step of tbfast (makedynamicmtx); the program then keeps the sequences apart. Hand run of the installed mafft 7.505: mafft --globalpair --fmodel --maxiterate 0.
    [Test]
    public void NegativeDistances_ShiftTheMatrixOfEveryStep()
    {
        var seqs = new[] { "CATACACCATAC", "CATACACCATAC", "CATACACCACAC", "CATACACCATAC", "CATACACCATAC", "CATACACCATAAC" };
        var rows = SequenceAligner.MultipleAlignMafft(seqs, new MafftOptions { AnySymbol = false, PairLibrary = MafftPairLibrary.Global, FrequencyModel = true }).AlignedSequences;
        Assert.That(rows, Is.EqualTo(new[]
        {
            "--------------------------------catacaccatac-----------------------",
            "---------------------catacaccatac----------------------------------",
            "-------------------------------------------catacaccacac------------",
            "c----------atacaccatac---------------------------------------------",
            "catacaccatac-------------------------------------------------------",
            "------------------------------------------------------catacaccataac",
        }));
    }

    [Test]
    public void Threads_DoNotChangeTheAlignment()
    {
        // the program is deterministic with several threads (verified on the committed sets: the 'thread' groups); the port is a single-thread run
        foreach (var r in Oracle.Value.Where(x => x.Strategy.EndsWith("-thread", StringComparison.Ordinal) && x.Rc == 0).Take(30))
        {
            var o = SequenceAligner_MafftProgressive_Tests.ParseOptions(r.Options);
            Assert.That(SequenceAligner.MultipleAlignMafft(r.Seqs, o).AlignedSequences, Is.EqualTo(SequenceAligner.MultipleAlignMafft(r.Seqs, o with { Threads = 0 }).AlignedSequences), r.Id);
        }
    }

    [Test]
    public void Contracts()
    {
        var three = new[] { "MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ" };
        var o = new MafftOptions { AnySymbol = false, PairLibrary = MafftPairLibrary.Global };
        Assert.Multiple(() =>
        {
            Assert.Throws<NotSupportedException>(() => SequenceAligner.MultipleAlignMafft(three, o with { MaxIterate = 2 }), "G-INS-i / L-INS-i / E-INS-i: the refinement with the pair library is not ported yet");
            Assert.Throws<NotSupportedException>(() => SequenceAligner.MultipleAlignMafft(three, o with { WeightI = 0.0 }), "--weighti 0: no constraints");
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMafft([.. three, "MK*U+"], o), "letters outside the alphabet stop the plain command");
            Assert.Throws<OperationCanceledException>(() => SequenceAligner.MultipleAlignMafft(three, o with { CancellationToken = new CancellationToken(true) }));
            Assert.That(SequenceAligner.RunMafft(three, o with { Quiet = true }).StandardError, Is.Empty);
        });
    }

    [Test]
    public void Strategy_Text_NamesThePairLibrary()
    {
        var seqs = new[] { "MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ" };
        var o = new MafftOptions { AnySymbol = false };
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.RunMafft(seqs, o with { PairLibrary = MafftPairLibrary.Global }).StandardError, Does.Contain("Strategy:\n G-INS-1 (Not tested.)\n Progressive method incorporating GLOBAL pairwise alignment information\n"));
            Assert.That(SequenceAligner.RunMafft(seqs, o with { PairLibrary = MafftPairLibrary.Local }).StandardError, Does.Contain("Strategy:\n L-INS-1 (Not tested.)\n Progressive method incorporating LOCAL pairwise alignment information\n"));
            Assert.That(SequenceAligner.RunMafft(seqs, o with { PairLibrary = MafftPairLibrary.Genaf }).StandardError, Does.Contain("Strategy:\n E-INS-1 (Not tested.)\n ?\n").And.Contain("Parameters for the E-INS-i option have been changed in version 7.243 (2015 Jun)."));
            Assert.That(SequenceAligner.RunMafft(seqs, o with { PairLibrary = MafftPairLibrary.Genaf, OldGenafParameters = true }).StandardError, Does.Contain("Obsolete parameters used for this calculation."));
        });
    }
}
