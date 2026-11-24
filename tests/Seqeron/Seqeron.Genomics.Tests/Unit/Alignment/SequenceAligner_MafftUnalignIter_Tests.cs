using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s5, B09 WP74, F345): <c>mafft --globalpair --unalignlevel x</c> / <c>--allowshift</c> WITH the iterative refinement (<c>--maxiterate n</c>) against the installed Debian mafft 7.505-1
/// (<c>scripts/alignment/generate_mafft_wp74_oracle.py unalign-iter</c>): return code, standard output, the files and the whole standard error. <c>dvtditr -s</c> computes the largest class of distance
/// (<c>calcmaxdistclass</c>: "ndistclass = 10, maxdistclass = n" on standard error), makes one scoring matrix per class (<c>makescoringmatrices</c>: the matrix shifted by the offset of the distance of the class), and at every branch of the guide tree
/// puts every pair of the two groups into the class of the sum of the distances of its sequences from the branch (<c>classifypairs</c>, <c>USEDISTONTREE</c>); <c>partA__align_variousdist</c> then adds the match scores of the profiles of the classes (the weights
/// of the class let through the pairs of other classes too: "Found a complex grouping." and the subtraction of those pairs by <c>match_calc_del</c>).
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftUnalignIter_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp74_unaligniter_oracle.tsv"));

    [Test]
    public void Oracle_Shape()
    {
        var all = Cases.Value;
        Assert.Multiple(() =>
        {
            Assert.That(all.Count, Is.GreaterThanOrEqualTo(200));
            var ok = all.Where(c => c.Rc == 0).ToList();
            Assert.That(ok.Count(c => c.Args.Contains("--allowshift")), Is.GreaterThanOrEqualTo(70), "--allowshift");
            Assert.That(ok.Count(c => c.Args.Contains("--unalignlevel")), Is.GreaterThanOrEqualTo(120), "--unalignlevel x");
            Assert.That(ok.Count(c => c.Args.Contains("--maxiterate")), Is.EqualTo(ok.Count), "every successful set is iterative");
            Assert.That(ok.Count(c => c.Args.Contains("--skipiterate")), Is.GreaterThanOrEqualTo(5), "--skipiterate");
            Assert.That(ok.Count(c => c.Args.Contains("--shiftpenalty")), Is.GreaterThanOrEqualTo(10), "with a shift penalty of its own");
            Assert.That(ok.Count(c => c.Args.Contains("--treein")), Is.GreaterThanOrEqualTo(5), "a user tree");
            Assert.That(ok.Count(c => c.Args.Contains("--core")), Is.GreaterThanOrEqualTo(8), "--core");
            Assert.That(ok.Count(c => c.Args.Contains("--weighti")), Is.GreaterThanOrEqualTo(10), "--weighti");
            Assert.That(all.Count(c => c.Rc != 0 && c.Err.Contains("supported only with the --globalpair option", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(10), "the other methods");
            Assert.That(all.Count(c => c.Rc != 0 && c.Err.Contains("Iterative refinment is not supported for --large", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(10), "--large with the iterative refinement");
        });
    }

    [Test]
    public void Oracle_EqualsMafft()
    {
        var failures = new List<string>();
        foreach (var c in Cases.Value)
        {
            if (SequenceAligner_MafftProgram_Tests.Verify(c) is { } failure)
                failures.Add(failure);
        }
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void AllowShiftWithIteration_ChangesTheAlignment()
    {
        // installed mafft 7.505-1: the same three DNA sequences with --globalpair --maxiterate 3 --reorder, with and without --allowshift
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = ">s0\nACTTGCAGATATGAAGCGTTCGGCGCTTCG\n>s1\nACAAGGCGATGAAACTTTCA\n>s2\nACAACGCAACGTACATATCA\n" };
        var shifted = MafftProgram.Run(["--quiet", "--globalpair", "--allowshift", "--maxiterate", "3", "--reorder", SequenceAligner_MafftProgram_Tests.InputName], files);
        var plain = MafftProgram.Run(["--quiet", "--globalpair", "--maxiterate", "3", "--reorder", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(shifted.ExitCode, Is.Zero);
        Assert.That(plain.ExitCode, Is.Zero);
        Assert.That(shifted.Stdout, Is.EqualTo(">s0\nacttgcagatatgaagcgttcggcgcttcg\n>s1\nac----------aaggcgatgaaactttca\n>s2\nac----------aacgcaacgtacatatca\n"));
        Assert.That(plain.Stdout, Is.EqualTo(">s0\nacttgcagatatgaagcgttcggcgcttcg\n>s1\n----------acaaggcgatgaaactttca\n>s2\n----------acaacgcaacgtacatatca\n"));
    }

    [Test]
    public void UnalignLevelWithIteration_ReportsTheClassesOfDistance()
    {
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = ">s0\nACTTGCAGATATGAAGCGTTCGGCGCTTCG\n>s1\nACAAGGCGATGAAACTTTCA\n>s2\nACAACGCAACGTACATATCA\n" };
        var r = MafftProgram.Run(["--globalpair", "--unalignlevel", "0.5", "--maxiterate", "2", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.Zero);
        Assert.That(r.Stderr, Does.Contain("ndistclass = 10, maxdistclass = 6\n"));
    }
}
