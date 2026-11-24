using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit R14 / A2-19c-hybrid, B09 WP78, F359): <c>mafft --hybridpair</c> without a fragment option (strategy Hybrid-INS-1) against the installed Debian mafft 7.505-1
/// (<c>scripts/alignment/generate_mafft_wp78_oracle.py</c>): return code, standard output, the files and the whole standard error of every run.
/// <list type="bullet">
/// <item>The script runs three programs (<c>mafft.tmpl</c> 2615 - 2624): <c>pairlocalalign -Y</c> (the local alignments of the old sequences against the new ones, <c>hat3</c>: with nothing added the file is empty),
/// <c>disttbfast -E 1 -T -y</c> (the 6-mer distances of all sequences, <c>hat2</c>, no alignment) and <c>tbfast</c> that loads <c>hat3</c> and <c>hat2</c> (the UPGMA tree of the distances, the progressive alignment with the blocks as constraints).</item>
/// <item>Two runs of the program end with a segmentation fault (valgrind on the installed binaries) and are refused: every <c>tbfast -S</c> (<c>--scoreout</c>), and <c>pairlocalalign -Y</c> with threads when nothing is added.</item>
/// </list>
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftHybridPair_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp78_oracle.tsv"));

    public static IEnumerable<TestCaseData> Groups() => Cases.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        var groups = Cases.Value.GroupBy(c => c.Group).ToDictionary(g => g.Key, g => g.ToList());
        Assert.Multiple(() =>
        {
            foreach (var g in new[] { "hybrid", "hybrid-opt", "hybrid-iter", "hybrid-tree", "hybrid-add" })
            {
                Assert.That(groups.ContainsKey(g), Is.True, g);
                if (groups.TryGetValue(g, out var list))
                    Assert.That(list.Count, Is.GreaterThanOrEqualTo(200), g);
            }
            Assert.That(Cases.Value.Count(c => c.Args.Contains("--hybridpair")), Is.GreaterThanOrEqualTo(1000), "--hybridpair in (nearly) every case (the last pair option wins in the others)");
            Assert.That(Cases.Value.Count(c => c.Rc == 0 && c.Err.Contains("Loading 'hat2' ... done.", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(20), "tbfast loads hat2 in the stored standard error (the oracle keeps the text of every 25th run and of the short ones)");
            Assert.That(Cases.Value.Count(c => c.Rc == 0 && c.Err.Contains("Hybrid-INS-1", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(10), "the strategy text");
        });
    }

    [TestCaseSource(nameof(Groups))]
    public void Group_EqualsMafft(string group)
    {
        var failures = new List<string>();
        int n = 0;
        foreach (var c in Cases.Value.Where(x => x.Group == group))
        {
            n++;
            if (SequenceAligner_MafftProgram_Tests.Verify(c) is { } failure)
                failures.Add(failure);
        }

        Assert.That(n, Is.GreaterThan(0));
        Assert.That(failures, Is.Empty, failures.Count + " of " + n + " differ; first:\n" + string.Join("\n", failures.Take(8)));
    }

    private const string FourProteins = ">a\nMKVLAAGIVGLLLAQPSAWAEEKT\n>b\nMKVLAGIVGLLAQPSAWAEEKTQ\n>c\nMKILAAGVVGLLLAQPSWAEEKT\n>d\nMKVLAAGIVGLLLAQSAWEKT\n";

    private static MafftProgramResult Run(string text, params string[] options) =>
        MafftProgram.Run([.. options, SequenceAligner_MafftProgram_Tests.InputName], new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = text });

    [Test]
    public void HybridPair_FourProteins_RunsTheThreeProgramsAndNamesTheStrategy()
    {
        // installed mafft 7.505-1: mafft --hybridpair in.fa: pairlocalalign (alg=Y), disttbfast (6-mer distances), tbfast ("Loading 'hat3'", "Loading 'hat2'"), "Hybrid-INS-1 (Not tested.)"
        var r = Run(FourProteins, "--hybridpair");
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(0));
            Assert.That(r.Stderr, Does.Contain("alg=Y, model=BLOSUM62"));
            Assert.That(r.Stderr, Does.Contain("Loading 'hat3' ... \ndone.\n"));
            Assert.That(r.Stderr, Does.Contain("Loading 'hat2' ... done.\nConstructing a UPGMA tree ... \n"));
            Assert.That(r.Stderr, Does.Contain("Strategy:\n Hybrid-INS-1 (Not tested.)\n ?\n"));
            Assert.That(r.Stdout, Does.StartWith(">a\n"));
        });
    }

    [Test]
    public void HybridPair_WithoutAddedSequences_ReadsAnEmptyPairLibrary_SoTheAlignmentEqualsTheDistanceOnlyRun()
    {
        // pairlocalalign -Y aligns old against new sequences only: with nothing added it writes no block, so the alignment is that of the 6-mer distances and the same progressive stage
        var hybrid = Run(FourProteins, "--hybridpair", "--quiet");
        var plain = Run(FourProteins, "--6merpair", "--quiet");
        Assert.That(hybrid.Stdout, Is.EqualTo(plain.Stdout));
    }

    [TestCase("--scoreout")]
    [TestCase("--localpair", "--scoreout")]
    [TestCase("--globalpair", "--scoreout")]
    public void ScoreoutWithTbfast_IsRefusedAsTheProgramReadsAReleasedTable(params string[] options)
    {
        // installed mafft 7.505-1 (valgrind on tbfast): sumofpairsscore runs after freeconstants(): "Invalid read of size 8", address 0x308 = NULL + 97 * 8, segmentation fault
        string[] all = options.Contains("--localpair") || options.Contains("--globalpair") ? options : ["--hybridpair", .. options];
        var ex = Assert.Throws<NotSupportedException>(() => Run(FourProteins, all))!;
        Assert.That(ex.Message, Does.Contain("freeconstants"));
    }

    [TestCase("--thread", "2")]
    [TestCase("--thread", "-1")]
    public void ThreadsWithoutAddedSequences_AreRefusedAsPairlocalalignReadsAJobThatDoesNotExist(params string[] options)
    {
        // installed mafft 7.505-1 (valgrind on pairlocalalign -Y): "Invalid read of size 1", address 0x0 in a thread, segmentation fault
        var ex = Assert.Throws<NotSupportedException>(() => Run(FourProteins, ["--hybridpair", .. options]))!;
        Assert.That(ex.Message, Does.Contain("pairlocalalign -Y"));
    }

    [Test]
    public void MissingLaraparamsFile_TheShellPrintsTheMessageOfCat_BeforeTheInstallText()
    {
        // installed mafft 7.505-1: mafft --larapair --laraparams nofile in.fa: "cat: nofile: No such file or directory", then the text "Install LaRA", exit status 1
        var r = Run(FourProteins, "--larapair", "--laraparams", "nofile");
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(1));
            Assert.That(r.Stdout, Is.Empty);
            Assert.That(r.Stderr, Does.StartWith("cat: nofile: No such file or directory\n\n== Install LaRA ="));
        });
    }

    [Test]
    public void MissingSeedtableFile_TheShellPrintsTheMessageOfCat_AndTheRunGoesOn()
    {
        // installed mafft 7.505-1: mafft --seedtable nofile in.fa prints the message of cat first and aligns the input (the table is empty)
        var r = Run(FourProteins, "--seedtable", "nofile");
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(0));
            Assert.That(r.Stderr, Does.StartWith("cat: nofile: No such file or directory\n"));
            Assert.That(r.Stdout, Does.StartWith(">a\n"));
        });
    }
}
