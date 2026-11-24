using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s4c / A2-19c-s5, B09 WP73, F334 ff.): the options of the script <c>mafft</c> that WP73 ported, against the installed Debian mafft 7.505-1
/// (<c>scripts/alignment/generate_mafft_wp73_oracle.py</c>): return code, standard output, the files and the whole standard error of every run.
/// <list type="bullet">
/// <item><c>nan</c> - a guide tree (<c>--treein</c>) with NaN or infinite branch lengths: the weights of the sequences become NaN, every comparison of the dynamic programs is false and the trace-back reads the cell of the shared
/// matrix <c>commonIP</c> that an earlier dynamic program left behind (<c>SequenceAligner.CommonIpState</c>).</item>
/// <item><c>core-anysymbol</c> - <c>--core --anysymbol</c> with a forced <c>--nuc</c> / <c>--amino</c>: <c>setcore</c> reads the alignment again, guesses the alphabet from its letters (upper case for protein) and stops with
/// "Illeagal character" on a letter that is not in it.</item>
/// <item><c>add-pairlib</c> / <c>add-large</c> - <c>--add</c> / <c>--addprofile</c> with a pair library (<c>tbfast -K -I n</c>, the error boxes of <c>tbfast</c>) and with <c>--large</c> / <c>--memsavetree</c> / <c>--stepadd</c> (<c>nodepair</c> + <c>tbfast -U -I n</c>).</item>
/// <item><c>smoothing</c> - <c>--add --smoothing</c> (<c>restorecommongapssmoothly</c>, <c>insertnewgaps_bothorders</c>: the letters next to the inserted columns are moved, the message "n out of m have been smoothed").</item>
/// </list>
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftWp73_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp73_oracle.tsv"));

    public static IEnumerable<TestCaseData> Groups() => Cases.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        Assert.Multiple(() =>
        {
            foreach (var g in Cases.Value.GroupBy(c => c.Group))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            var nan = Cases.Value.Where(c => c.Group == "nan").ToList();
            Assert.That(nan.Count(c => c.Files.Values.Any(v => v.Contains("nan", StringComparison.OrdinalIgnoreCase))), Is.GreaterThanOrEqualTo(150), "NaN branch lengths");
            Assert.That(nan.Count(c => c.Files.Values.Any(v => v.Contains("inf", StringComparison.OrdinalIgnoreCase))), Is.GreaterThanOrEqualTo(150), "infinite branch lengths");
            Assert.That(nan.Count(c => c.Args.Contains("--genafpair") || c.Args.Contains("--oldgenafpair")), Is.GreaterThanOrEqualTo(30), "the case the port got wrong before: the pair library of --genafpair");
            Assert.That(nan.Count(c => c.Rc != 0), Is.GreaterThanOrEqualTo(40), "the stops: dvtditr (mscore / tscore became NaN) and the tree file (glibc scanf reads nan, inf, infinity but not nan(0x12))");
            Assert.That(nan.Count(c => c.Err.Contains("Format error (2)", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(15), "the tree file is read with sscanf");
            Assert.That(nan.Count(c => c.Args.Contains("--weighti")), Is.GreaterThanOrEqualTo(30), "no pair library kept: L__align11_noalign does not touch commonIP");
            var core = Cases.Value.Where(c => c.Group == "core-anysymbol").ToList();
            Assert.That(core.Count(c => c.Rc == 0 && c.Out.Length == 0), Is.GreaterThanOrEqualTo(5), "setcore stops on a letter of replaceu that its guessed alphabet does not know (\"Illeagal character\", an empty output)");
            Assert.That(core.Count(c => c.Out.Length > 0), Is.GreaterThanOrEqualTo(150), "a core is written");
            var pairLibrary = Cases.Value.Where(c => c.Group == "add-pairlib").ToList();
            Assert.That(pairLibrary.Count(c => c.Rc == 0 && c.Out.Length > 0), Is.GreaterThanOrEqualTo(250), "an alignment");
            Assert.That(pairLibrary.Count(c => c.Args.Contains("--addprofile")), Is.GreaterThanOrEqualTo(50), "--addprofile");
            Assert.That(pairLibrary.Count(c => c.Rc != 0), Is.GreaterThanOrEqualTo(10), "the stops of --addprofile (the boxes of tbfast)");
            Assert.That(pairLibrary.Count(c => c.Args.Contains("--keeplength")), Is.GreaterThanOrEqualTo(20), "--keeplength");
            Assert.That(pairLibrary.Count(c => c.OutFiles.Keys.Any(k => k.EndsWith(".map", StringComparison.Ordinal))), Is.GreaterThanOrEqualTo(10), "--mapout");
            var large = Cases.Value.Where(c => c.Group == "add-large").ToList();
            Assert.That(large.Count(c => c.Rc == 0 && c.Out.Length > 0), Is.GreaterThanOrEqualTo(120), "an alignment");
            Assert.That(large.Count(c => c.Err.Contains("Iterative refinment is not supported", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(50), "the script refuses --maxiterate with --large / --memsavetree / --stepadd");
            Assert.That(large.Count(c => c.Args.Contains("--keeplength")), Is.GreaterThanOrEqualTo(20), "--keeplength");
            var smoothing = Cases.Value.Where(c => c.Group == "smoothing").ToList();
            Assert.That(smoothing.Count(c => c.Args.Contains("--keeplength") || c.Args.Contains("--mapout")), Is.GreaterThanOrEqualTo(15), "--smoothing with --keeplength / --mapout");
            Assert.That(smoothing.Count(c => c.Args.Contains("--addprofile")), Is.GreaterThanOrEqualTo(20), "--addprofile with --smoothing");
            Assert.That(smoothing.Count(c => c.Args.Any(a => a is "--localpair" or "--globalpair" or "--genafpair" or "--large")), Is.GreaterThanOrEqualTo(70), "tbfast and nodepair");
        });
    }

    [Test]
    public void CommonIp_IsAllocatedZeroedOnGrowthAndKeepsTheCellsOfEarlierPrograms()
    {
        var c = new SequenceAligner.CommonIpState();
        var a = c.IpFor(5, 7);
        Assert.That(a, Has.Length.EqualTo(15), "AllocateIntMtx( max( orlgth1, commonAlloc1 ) + 10, max( orlgth2, commonAlloc2 ) + 10 )");
        Assert.That(a[0], Has.Length.EqualTo(17));
        a[3][4] = 9;
        Assert.That(c.IpFor(4, 6), Is.SameAs(a), "no program needs more rows or columns: the matrix is shared and not cleared");
        Assert.That(c.IpFor(4, 6)[3][4], Is.EqualTo(9), "the cell an earlier program wrote");
        var grown = c.IpFor(6, 6);
        Assert.That(grown, Is.Not.SameAs(a), "orlgth1 > commonAlloc1: a new matrix");
        Assert.That(grown[3][4], Is.EqualTo(0), "calloc");
        Assert.That(grown, Has.Length.EqualTo(16));
        Assert.That(grown[0], Has.Length.EqualTo(17), "the columns keep max( orlgth2, commonAlloc2 )");
        c.FreeIp();
        Assert.That(c.Alloc1 + c.Alloc2, Is.Zero, "the switch to the memory-saving program");
        Assert.That(c.IpFor(2, 2)[3][4], Is.EqualTo(0));
    }

    [Test]
    public void DpCapacity_FollowsTheStaticArraysOfTheDynamicPrograms()
    {
        var cap = new SequenceAligner.DpCapacity();
        cap.Fit(10, 20);
        Assert.That((cap.Orlgth1, cap.Orlgth2), Is.EqualTo((13, 26)), "ll = max( (int)(1.3 * lgth), orlgth ) + 100; orlgth = ll - 100");
        cap.Fit(10, 25);
        Assert.That((cap.Orlgth1, cap.Orlgth2), Is.EqualTo((13, 26)), "both lengths fit: unchanged");
        cap.Fit(30, 25);
        Assert.That((cap.Orlgth1, cap.Orlgth2), Is.EqualTo((39, 32)), "a program that has to grow recomputes both capacities: max( (int)(1.3 * 25), 26 )");
        cap.Reset();
        Assert.That((cap.Orlgth1, cap.Orlgth2), Is.EqualTo((0, 0)), "A__align( NULL ... )");
    }

    [Test]
    public void NanBranchLength_FollowsTheStaleTraceBackCellOfTheEarlierAlignment()
    {
        // installed mafft 7.505-1: the weights are NaN, every comparison is false, ijp[lgth1][lgth2] keeps the value of the alignment of the first step (the trace-back ends with a gap behind the second sequence)
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">s0\nAGTGAAA\n>s1\nAATGAA\n>s2\nAATGAA\n",
            ["t.txt"] = "1 2 0.1 0.2\n1 3 0.1 nan\n",
        };
        var r = MafftProgram.Run(["--quiet", "--treein", "t.txt", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.Zero);
        Assert.That(r.Stdout, Is.EqualTo(">s0\nagtgaaa\n>s1\naatgaa-\n>s2\naatgaa-\n"));
        // --genafpair with two sequences and NaN: the pair library is made in the same process; the real alignment puts the two sequences next to each other
        var two = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">s0\nMKVLAAGIVGLLLA\n>s1\nMKVLAGIVGLLA\n",
            ["t.txt"] = "1 2 nan 0.1\n",
        };
        var g = MafftProgram.Run(["--quiet", "--genafpair", "--treein", "t.txt", SequenceAligner_MafftProgram_Tests.InputName], two);
        Assert.That(g.Stdout, Is.EqualTo(">s0\n----------MKVLAAGIVGLLLA\n>s1\nMKVLAGIVGLLA------------\n"));
    }

    [Test]
    public void AddWithAPairLibrary_KeepsTheColumnsOfTheOldAlignment()
    {
        // installed mafft 7.505-1: --add new.fa --localpair (tbfast -K -I 1): the old rows keep their columns, the new row follows the gap pattern of the old ones
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">s0\nP-WDVENWCCHCDGWDITVQC\n>s1\nP-CDVENWCPHCDIQPIDVQG\n>s2\nPCWFVENWCTHPHYFDHKVSC\n",
            ["new.fa"] = ">n0\nMCDVENWCTVCDQQTIDYQC\n",
        };
        var r = MafftProgram.Run(["--quiet", "--add", "new.fa", "--localpair", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.Stdout, Is.EqualTo(">s0\nP-WDVENWCCHCDGWDITVQC\n>s1\nP-CDVENWCPHCDIQPIDVQG\n>s2\nPCWFVENWCTHPHYFDHKVSC\n>n0\nM-CDVENWCTVCDQQTIDYQC\n"));
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
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }
}
