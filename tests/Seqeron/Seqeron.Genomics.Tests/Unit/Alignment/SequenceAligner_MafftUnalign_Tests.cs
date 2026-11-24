using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s5, B09 WP74, F344): <c>mafft --globalpair --unalignlevel x</c> and <c>--allowshift</c> (= <c>--unalignlevel 0.8 --shiftpenalty 2</c>, <c>-:</c>) without the iterative refinement against the installed Debian mafft 7.505-1
/// (<c>scripts/alignment/generate_mafft_wp74_oracle.py unalign</c>): return code, standard output, the files and the whole standard error. The scoring matrix of an alignment step is shifted by the offset
/// <c>min( 0, dist / 2 - x ) * 600</c> that belongs to the distance of the two groups (<c>makedynamicmtx</c>, <c>dist2offset</c>): <c>tbfast</c> does it for every step of its guide tree (distance from the tips), and the pair stage
/// (<c>pairlocalalign -u</c>, <c>tbfast</c> with the user tree: the distances along the tree) aligns a pair whose distance is large enough a second time with the shifted matrix and keeps the blocks of the second alignment.
/// The script allows both options only with <c>--globalpair</c> and a pair library weight of 1 or more; every other method stops with its message.
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftUnalign_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp74_unalign_oracle.tsv"));

    [Test]
    public void Oracle_Shape()
    {
        var all = Cases.Value;
        Assert.Multiple(() =>
        {
            Assert.That(all.Count, Is.GreaterThanOrEqualTo(200));
            var ok = all.Where(c => c.Rc == 0).ToList();
            Assert.That(ok.Count(c => c.Args.Contains("--allowshift")), Is.GreaterThanOrEqualTo(80), "--allowshift");
            Assert.That(ok.Count(c => c.Args.Contains("--unalignlevel")), Is.GreaterThanOrEqualTo(100), "--unalignlevel x");
            Assert.That(ok.Count(c => c.Args.Contains("--shiftpenalty")), Is.GreaterThanOrEqualTo(15), "with a shift penalty of its own");
            Assert.That(ok.Count(c => c.Args.Contains("--large") || c.Args.Contains("--memsavetree") || c.Args.Contains("--stepadd")), Is.GreaterThanOrEqualTo(15), "nodepair");
            Assert.That(ok.Count(c => c.Args.Contains("--treein")), Is.GreaterThanOrEqualTo(4), "the distances along a user tree");
            Assert.That(ok.Count(c => c.Args.Contains("--weighti")), Is.GreaterThanOrEqualTo(10), "--weighti");
            Assert.That(ok.Count(c => c.Args.Contains("--treeout") || c.Args.Contains("--distout")), Is.GreaterThanOrEqualTo(20), "tree / distance files");
            Assert.That(all.Count(c => c.Rc != 0 && c.Err.Contains("supported only with the --globalpair option", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(10), "the other methods");
            Assert.That(all.Count(c => c.Rc != 0 && c.Err.Contains("Do not simultaneously specify --similaritylevel and --unalignlevel", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(3), "--similaritylevel");
            Assert.That(all.Count(c => c.Rc != 0 && c.Err.Contains("must be between 0.0 and 1.0", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(1), "x above 1");
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
    public void AllowShift_ChangesTheAlignment()
    {
        // installed mafft 7.505-1: the same three protein sequences with --globalpair --bl 80 and with --globalpair --allowshift --bl 80
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = ">s0\nAKMMYRIHIMYCYDKKGEW\n>s1\nTSKMYRIHPMYCYAKKGEW\n>s2\nTAKMMYRIMHIMYCYDKDGEW\n" };
        var shifted = MafftProgram.Run(["--quiet", "--globalpair", "--allowshift", "--bl", "80", SequenceAligner_MafftProgram_Tests.InputName], files);
        var plain = MafftProgram.Run(["--quiet", "--globalpair", "--bl", "80", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(shifted.ExitCode, Is.Zero);
        Assert.That(plain.ExitCode, Is.Zero);
        Assert.That(shifted.Stdout, Is.EqualTo(">s0\n-AKMMYRI-HIMYCYDKKGEW\n>s1\nTSK-MYRI-HPMYCYAKKGEW\n>s2\nTAKMMYRIMHIMYCYDKDGEW\n"));
        Assert.That(plain.Stdout, Is.EqualTo(">s0\n-AKMMYRI-HIMYCYDKKGEW\n>s1\n-TSKMYRI-HPMYCYAKKGEW\n>s2\nTAKMMYRIMHIMYCYDKDGEW\n"));
    }

    [Test]
    public void UnalignLevel_LeavesTheDistantPartUnaligned()
    {
        // installed mafft 7.505-1: DNA, --unalignlevel 0.1 --shiftpenalty 5 against --globalpair --shiftpenalty 5
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = ">s0\nGTAGTTGCTGCAATTTGTATCGTATTGCTT\nGGATGGTAAT\n>s1\nAAGGACGGGGAATGTGTGTT\n>s2\nAAGGACGGGGAGTTTTTGAT\n" };
        var level = MafftProgram.Run(["--quiet", "--globalpair", "--unalignlevel", "0.1", "--shiftpenalty", "5", SequenceAligner_MafftProgram_Tests.InputName], files);
        var plain = MafftProgram.Run(["--quiet", "--globalpair", "--shiftpenalty", "5", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(level.Stdout, Is.EqualTo(">s0\ngtagttgctgcaatttgtatcgtattgcttggatggtaat\n>s1\naaggacggggaatgtg--------------------tgtt\n>s2\naaggacggggagtttt--------------------tgat\n"));
        Assert.That(plain.Stdout, Is.EqualTo(">s0\ngtagttgctgcaatttgtatcgtattgcttggatggtaat\n>s1\naaggacg-gggaatgtgtgt-------------------t\n>s2\naaggacg-gggagtttttga-------------------t\n"));
    }

    [Test]
    public void UnalignLevelAbove1_StopsTheScript()
    {
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nMKVLAAGIVG\n>b\nMKVLAGIVG\n" };
        var r = MafftProgram.Run(["--globalpair", "--unalignlevel", "1.5", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.EqualTo(1));
        Assert.That(r.Stderr, Is.EqualTo("\nThe 'unalignlevel' parameter must be between 0.0 and 1.0\n\n"));
    }
}
