using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s4c, B09 WP72, F326 ff.): the guide trees of MAFFT 7.505 that are not the UPGMA of the full matrix - <c>--treein</c> (<c>loadtree</c>), <c>--pileup</c>, <c>--randomchain</c> (<c>createchain</c>),
/// <c>--large</c> / <c>--memsavetree</c> / <c>--memsavetreex</c> / <c>--stepadd</c> / <c>--youngestlinkage</c> (<c>compacttree</c>), <c>--large</c> with a pair library (<c>nodepair</c>) -, the memory-saving alignment
/// (<c>--memsave</c>, <c>MSalignmm</c>, <c>Falign_udpari_long</c>, the switch for groups of more than 30 000 / 10 000 columns) and <c>--auto</c> for 20 000 sequences, against the installed Debian mafft 7.505-1
/// (<c>scripts/alignment/generate_mafft_compacttree_oracle.py</c>): return code, standard output, the files and the standard error (SHA-256) of every run.
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftCompactTree_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Trees = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_compacttree_oracle.tsv"));
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> MemSave = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_memsave_oracle.tsv"));
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Auto20k = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_auto20k_oracle.tsv"));

    private static readonly string[] HeavyGroups = ["longcols", "longparttree", "auto20k"];

    public static IEnumerable<TestCaseData> TreeGroups() => Trees.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    public static IEnumerable<TestCaseData> MemSaveGroups() => MemSave.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        Assert.Multiple(() =>
        {
            foreach (var g in Trees.Value.GroupBy(c => c.Group))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            foreach (var g in MemSave.Value.Where(c => !HeavyGroups.Contains(c.Group)).GroupBy(c => c.Group))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(MemSave.Value.Count(c => c.Group == "longcols"), Is.GreaterThanOrEqualTo(3), "groups of more than 30 000 columns (MSalignmm)");
            Assert.That(MemSave.Value.Count(c => c.Group == "longparttree"), Is.GreaterThanOrEqualTo(4), "groups of more than 10 000 columns in splittbfast");
            Assert.That(Auto20k.Value.Single().Input.Count(ch => ch == '>'), Is.EqualTo(20000), "--auto reaches the memory-saving tree");
            Assert.That(Trees.Value.Count(c => c.OutFiles.Keys.Any(k => k.EndsWith(".tree", StringComparison.Ordinal))), Is.GreaterThanOrEqualTo(300), "guide tree files");
            Assert.That(Trees.Value.Count(c => c.Rc != 0), Is.GreaterThanOrEqualTo(150), "the refusals of the program (trees that are no trees)");
            Assert.That(Trees.Value.Count(c => c.Err.Length > 0), Is.GreaterThanOrEqualTo(150), "full texts of the standard error");
        });
    }

    [TestCaseSource(nameof(TreeGroups))]
    public void GuideTreeFamily_EqualsMafft(string group) => Check(Trees.Value.Where(c => c.Group == group));

    [TestCaseSource(nameof(MemSaveGroups))]
    public void MemorySavingAlignment_EqualsMafft(string group) => Check(MemSave.Value.Where(c => c.Group == group));

    [Test]
    public void AutoForTwentyThousandSequences_EqualsMafft() => Check(Auto20k.Value);

    private static void Check(IEnumerable<SequenceAligner_MafftProgram_Tests.Case> cases)
    {
        var failures = new List<string>();
        int n = 0;
        foreach (var c in cases)
        {
            n++;
            if (SequenceAligner_MafftProgram_Tests.Verify(c) is { } failure)
                failures.Add(failure);
        }
        Assert.That(n, Is.GreaterThan(0));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    // ---------------------------------------------------------------- the C library functions the programs use

    [Test]
    public void CScanner_ReadsNumbersAsGlibcSscanf()
    {
        Assert.Multiple(() =>
        {
            var s = new SequenceAligner.MafftRun.CScanner("  12 -7 +3 x");
            Assert.That(s.Int(out int a) && a == 12, Is.True);
            Assert.That(s.Int(out int b) && b == -7, Is.True);
            Assert.That(s.Int(out int c) && c == 3, Is.True);
            Assert.That(s.Int(out _), Is.False, "a letter is no integer");
            // %lf: the collected characters go to strtod, the rest of a malformed exponent is consumed ("1e" is read as 1)
            var d = new SequenceAligner.MafftRun.CScanner("1e .5 2.5e-1 0x1p-3 inf -inf nan");
            Assert.That(d.Double(out double v1) && v1 == 1.0, Is.True);
            Assert.That(d.Double(out double v2) && v2 == 0.5, Is.True);
            Assert.That(d.Double(out double v3) && v3 == 0.25, Is.True);
            Assert.That(d.Double(out double v4) && v4 == 0.125, Is.True);
            Assert.That(d.Double(out double v5) && double.IsPositiveInfinity(v5), Is.True);
            Assert.That(d.Double(out double v6) && double.IsNegativeInfinity(v6), Is.True);
            Assert.That(d.Double(out double v7) && double.IsNaN(v7), Is.True);
            Assert.That(new SequenceAligner.MafftRun.CScanner(".e1").Double(out _), Is.False);
            Assert.That(new SequenceAligner.MafftRun.CScanner("1 2").Int(out int big) && big == 1, Is.True);
        });
    }

    [Test]
    public void FormatOfNaNAndInfinity_FollowsPrintf()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.MafftGlobals.F(double.PositiveInfinity, 7, 5), Is.EqualTo("    inf"));
            Assert.That(SequenceAligner.MafftGlobals.F(double.NegativeInfinity, 7, 5), Is.EqualTo("   -inf"));
            Assert.That(SequenceAligner.MafftGlobals.F(BitConverter.Int64BitsToDouble(0x7FF8000000000000L), 7, 5), Is.EqualTo("    nan"));
            Assert.That(SequenceAligner.MafftGlobals.F(double.NaN, 7, 5), Is.EqualTo("   -nan"), "the NaN of an operation on x86 has its sign bit set");
            Assert.That(SequenceAligner.MafftGlobals.D5(7), Is.EqualTo("    7"));
            Assert.That(SequenceAligner.MafftGlobals.D5(12345), Is.EqualTo(" 12345"), "% 5d: the blank of the sign is kept");
        });
    }

    [Test]
    public void Contracts()
    {
        string fasta = ">a\nMKVLAAGIVGLLLA\n>b\nMKVLAGIVGLLA\n>c\nMKVIAAGLVGLLLA\n>d\nMKVLAAGIVGLLLG\n";
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = fasta, ["t.txt"] = "1 2 0.1 0.2\n1 3 0.1 0.3\n1 4 0.2 0.1\n" };
        Assert.Multiple(() =>
        {
            // the program writes no distance matrix for these trees: its cp message names a random temporary directory (--distout), or it stops with a segmentation fault that names a process id
            Assert.Throws<NotSupportedException>(() => MafftProgram.Run(["--pileup", "--distout", SequenceAligner_MafftProgram_Tests.InputName], files));
            Assert.Throws<NotSupportedException>(() => MafftProgram.Run(["--large", "--distout", SequenceAligner_MafftProgram_Tests.InputName], files));
            Assert.Throws<NotSupportedException>(() => MafftProgram.Run(["--treein", "t.txt", "--distout", SequenceAligner_MafftProgram_Tests.InputName], files));
            // a tree that names a missing sequence stops the program with its message
            var bad = new Dictionary<string, string>(files) { ["t.txt"] = "1 9 0.1 0.2\n" };
            var r = MafftProgram.Run(["--treein", "t.txt", SequenceAligner_MafftProgram_Tests.InputName], bad);
            Assert.That(r.ExitCode, Is.EqualTo(1));
            Assert.That(r.Stderr, Does.Contain("Check the guide tree."));
            // --pileup adds the sequences in the order of the input: every sequence is a row, the rows keep the input order
            var chain = MafftProgram.Run(["--quiet", "--pileup", SequenceAligner_MafftProgram_Tests.InputName], files);
            Assert.That(chain.ExitCode, Is.EqualTo(0));
            Assert.That(chain.Stdout.Split('>', StringSplitOptions.RemoveEmptyEntries).Select(x => x.Split('\n')[0]), Is.EqualTo(new[] { "a", "b", "c", "d" }));
            // a user tree that is valid gives a run
            var treeRun = MafftProgram.Run(["--quiet", "--treein", "t.txt", SequenceAligner_MafftProgram_Tests.InputName], files);
            Assert.That(treeRun.ExitCode, Is.EqualTo(0));
        });
    }

    [Test]
    public void LibraryOptions_MatchTheCommandLine()
    {
        var seqs = new[] { "MKVLAAGIVGLLLAQWE", "MKVLAGIVGLLAQWE", "MKVIAAGLVGLLLAQWE", "MKVLAAGIVGLLLGQWE", "MKVLSAGIVGLLAQWE" };
        string fasta = string.Concat(seqs.Select((s, i) => ">s" + i.ToString(CultureInfo.InvariantCulture) + "\n" + s + "\n"));
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = fasta };
        foreach (var (method, option) in new[] { (MafftGuideTree.Pileup, "--pileup"), (MafftGuideTree.RandomChain, "--randomchain"), (MafftGuideTree.MemSaveTree, "--large"), (MafftGuideTree.StepAdd, "--stepadd"), (MafftGuideTree.YoungestLinkage, "--youngestlinkage") })
        {
            var viaApi = SequenceAligner.RunMafft(seqs, new MafftOptions { GuideTree = method, Retree = method is MafftGuideTree.Pileup or MafftGuideTree.RandomChain ? 1 : 2, Quiet = true, AnySymbol = false });
            var viaScript = MafftProgram.Run(["--quiet", option, SequenceAligner_MafftProgram_Tests.InputName], files);
            var rows = viaScript.Stdout.Split('>', StringSplitOptions.RemoveEmptyEntries).Select(x => string.Concat(x.Split('\n').Skip(1))).ToArray();
            Assert.That(viaApi.Alignment.AlignedSequences, Is.EqualTo(rows), option);
        }
    }
}
