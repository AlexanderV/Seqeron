using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s5, B09 WP77, F357): <c>mafft --sp</c> (the programs run with <c>-@</c>: <c>D__align</c> and <c>D__align_variousdist</c> of <c>Dalignmm.c</c>) against the installed Debian mafft 7.505-1
/// (<c>scripts/alignment/generate_mafft_wp77_oracle.py</c>): return code, standard output, the files and the whole standard error of every run (the program prints the progress, the score and the check of the score of every step).
/// <list type="bullet">
/// <item>The dynamic program weights the gap penalty of every cell by the probability that the gap is new: the table of the lengths of the runs of gaps of both groups gives, for the cell, the probability that a gap of the other group
/// continues the gap that the move opens (<c>calcpfac_gap_noidatend</c>, <c>calcpfac_gapex_noidatend</c>, <c>calcpfacnoidatend</c>); the score of the result is checked against the pairwise scores of its rows (the program stops with exit status 1 when they differ).</item>
/// <item><c>D__align_ls</c> takes the group with more columns behind a run of gaps as the first one; two single sequences are aligned by <c>G__align11</c>.</item>
/// <item>With a pair library (<c>constraint</c>) the importance of the blocks enters the match scores and the program adds up the importance on the path (<c>Atracking_localhom</c>); with <c>--unalignlevel</c> every pair of sequences of the two groups
/// has the scoring matrix of its distance class (<c>D__align_variousdist</c>).</item>
/// </list>
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftSp_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp77_oracle.tsv"));

    public static IEnumerable<TestCaseData> Groups() => Cases.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        var groups = Cases.Value.GroupBy(c => c.Group).ToDictionary(g => g.Key, g => g.ToList());
        Assert.Multiple(() =>
        {
            foreach (var g in new[] { "sp", "sp-opt", "sp-iter", "sp-pair", "sp-tree", "sp-unalign", "sp-add" })
            {
                Assert.That(groups.ContainsKey(g), Is.True, g);
                if (groups.TryGetValue(g, out var list))
                    Assert.That(list.Count, Is.GreaterThanOrEqualTo(200), g);
            }
            Assert.That(Cases.Value.All(c => c.Args.Contains("--sp")), Is.True, "every case runs --sp");
            Assert.That(Cases.Value.Count(c => c.Rc == 0 && c.Err.Contains("sequences, len=", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(10), "the progress of D__align in the stored standard error");
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

    [Test]
    public void Sp_FourProteins_EqualsTheInstalledMafft()
    {
        // installed mafft 7.505-1: mafft --sp in.fa (the standard error is the progress of D__align: v1,v2, wm, kenzan, pairscore, diff)
        var r = MafftProgram.Run(["--sp", SequenceAligner_MafftProgram_Tests.InputName], new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = FourProteins });
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(0));
            Assert.That(r.Stdout, Is.EqualTo(">a\nMKVLAAGIVGLLLAQPSAWAEEKT-\n>b\nMKVL-AGIVG-LLAQPSAWAEEKTQ\n>c\nMKILAAGVVGLLLAQPS-WAEEKT-\n>d\nMKVLAAGIVGLLLAQSA-W--EKT-\n"));
            Assert.That(r.Stderr, Does.Contain("2 x 1 sequences, len=24, 21\nwm = 9365.500000\n*impmatch = 0.000000\nkenzan = 3 -> 1.500000\npairscore = 9365.500000\npairscore-nogappairscore = -1375.500000\npairscore-nogappairscore / penalty = 1.500000\ndiff = -0.000000\n"));
            Assert.That(r.Stderr, Does.Contain("3 x 1 sequences, len=24, 23\nwm = 7528.354502\n"));
            Assert.That(r.Stderr, Does.Contain("alg=d, model=BLOSUM62, 1.53, -0.00, -0.00, noshift, amax=0.0"));
        });
    }

    [Test]
    public void SpWithParttree_SplittbfastDoesNotKnowTheOptionAndStops()
    {
        // installed mafft 7.505-1: mafft --sp --parttree in.fa: "illegal option @", "options: Check source file !", exit status 1
        var r = MafftProgram.Run(["--sp", "--parttree", SequenceAligner_MafftProgram_Tests.InputName], new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = FourProteins });
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(1));
            Assert.That(r.Stderr, Does.StartWith("illegal option @\noptions: Check source file !\n"));
        });
    }

    [Test]
    public void SpAdd_TheCheckReadsTheMatrixBeforeItsFirstRow_IsRefusedWithTheProof()
    {
        // installed mafft 7.505-1: this run ends with a segmentation fault of disttbfast (valgrind: "Invalid read of size 8 ... Address is 8 bytes before a block of size 208 alloc'd"): naivepairscore11_dynmtx indexes the matrix with
        // amino_n['='] = -1, the new gap character of --add, in the check of D__align
        const string aligned = ">a0\nSWCQ-EVHSCRENCGKMQNEHQLHDFGY---FKVLPMRPNEV-CSHCRPMPTWEAM-QLN\nFHG-NSRWD-NSF-SGCTQR----IEGV-MGK---MT-TMEPQTWVFWN-GLSKTLDPP-\nKIEICDDSREEHNLFRKCLMTK-\n"
            + ">a1\n-CCQ-E-HSCRE--GKMQNENT-HDFWHM--FKVLP-RPEFVQWHCNIPMPFWEAMQQLN\nFAC-NQY-NQNLF-SFCRRRQ---YEGV--YGKQ-MTWPMPPQTWVFGN-GVSKTLCIPC\nKI-IDDDSREK-NFK--CLHTK-\n"
            + ">a2\n-SVQPEVHSCREN-GKMQTEHELHDFMHAY-FKVYP-RP-VM-CHCCWPPSSWEGM-QLN\nFYG-NYR-D-NSF-SFCRRRQT--YG-V-MGQKQCMQ-PGTPQTWVFWNFGVSKTLCPP-\nKI-I-DDSREEHVNF--CCHTKD\n";
        const string added = ">n0\nSRQEVWNCRNGKMQNETLHYDFLHAYFKVLPRPEVCHCCIPMPWEAMYQNFCPNYRDNSFSFCRRNRYEGVYGKQMTGMPPQTWLFWNPVSKTLCPPKIIDDDSREEHCFKCLHTK\n"
            + ">n1\nSQEVHSCRENGKMQNEHQLHDFLHAYGFKVLPRPEMCHCCIPMPSWAMQLNKFCGNYRDNSFSFRYRQRYEGVMYEGKQMTPMPPQTWVFWVGVSATLCIPKIIDDDSIEEHNGGECLHYTK\n";
        var ex = Assert.Throws<NotSupportedException>(() => MafftProgram.Run(["--sp", "--add", "new.fa", SequenceAligner_MafftProgram_Tests.InputName],
            new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = aligned, ["new.fa"] = added }))!;
        Assert.That(ex.Message, Does.Contain("mtx[-1]"));
    }

    [TestCase("--maxiterate", "2", "--scoreout")]
    [TestCase("--sp", "--maxiterate", "2", "--scoreout")]
    public void MaxiterateWithScoreout_IsRefusedAsTheProgramReadsAFreedBlock(params string[] options)
    {
        // installed mafft 7.505-1: dvtditr -S ends with a segmentation fault (valgrind: invalid read of a freed block of 48 bytes), with and without --sp
        var ex = Assert.Throws<NotSupportedException>(() => MafftProgram.Run([.. options, SequenceAligner_MafftProgram_Tests.InputName], new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = FourProteins }))!;
        Assert.That(ex.Message, Does.Contain("freed"));
    }

    [TestCase("--sp", "--globalpair", "--unalignlevel", "0.3", "--maxiterate", "2")]
    [TestCase("--sp", "--globalpair", "--allowshift", "--maxiterate", "2")]
    public void SpWithTwoSequencesAndUnalignlevel_TakeTheDistancesFromTheOnlyNodeOfTheTree(params string[] options)
    {
        // two sequences: distFromABranch gives the two branch lengths of the tree (len[0][0], len[0][1]), so the pair has a distance class; the iteration prints one line per step
        var r = MafftProgram.Run([.. options, SequenceAligner_MafftProgram_Tests.InputName], new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nMKVLAAGIVGLLLAQPSAWAEEKT\n>b\nMKVLAGIVGLLAQPSAWAEEKTQ\n" });
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(0));
            Assert.That(r.Stderr, Does.Contain("1 x 1 sequences"));
        });
    }
}
