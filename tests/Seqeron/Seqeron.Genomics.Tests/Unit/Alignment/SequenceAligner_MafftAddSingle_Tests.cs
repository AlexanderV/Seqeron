using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-addsingle, B09 WP76, F352 - F354): <c>mafft --addfragments / --addfull / --addlong / --addtotop</c> (the program <c>addsingle</c> and the programs of its chains: <c>dndpre -M 2</c>,
/// <c>pairlocalalign -Y / -L</c>, <c>disttbfast -T -y</c>) against the installed Debian mafft 7.505-1 (<c>scripts/alignment/generate_mafft_wp75_oracle.py</c>): return code, standard output, the files and the whole standard error of every run.
/// <list type="bullet">
/// <item><c>frag-6mer</c> (<c>--6merpair</c> / <c>--10merpair</c>: <c>addsingle -W n</c>, the k-mer distances of the program, <c>addonetip</c> / <c>addonetip2top</c> (<c>--addtotop</c>), <c>treebase</c> with the new columns kept as '=', the gap lists of all new
/// sequences, <c>dorealignment_tree</c> for the residues that several sequences put into one place, <c>--keeplength</c>, <c>--mapout</c>, <c>--smoothing</c>, <c>--treeout</c>, <c>--distout</c>), <c>frag-6mermulti</c> (<c>-d</c>: <c>dndpre</c> of the old
/// sequences); <c>frag-multi</c> (the default, <c>--multipair</c>: <c>dndpre -M 2</c>, <c>pairlocalalign -Y</c>, the blocks of <c>hat3</c> as the importance of the alignment step, <c>c</c> in the progress), <c>frag-hybrid</c> (<c>disttbfast -E 1 -T -y</c> overwrites
/// <c>hat2</c>), <c>frag-local</c> (<c>pairlocalalign -L</c>, the full table of blocks), <c>frag-codon</c> (<c>--codonpos</c> / <c>--codonscore</c>: the gap penalties of every column of the reference and the score of the codon pairs in <c>G__align11psg</c>),
/// <c>frag-treein</c> (<c>loadtop</c>), <c>frag-opt</c> (further options of the script).</item>
/// <item>A single reference sequence with <c>-t</c> prints the pointer in front of its array of trees: the first line of the tree file is "(null);".</item>
/// </list>
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftAddSingle_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp75_oracle.tsv"));

    public static IEnumerable<TestCaseData> Groups() => Cases.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        Assert.Multiple(() =>
        {
            foreach (var g in Cases.Value.GroupBy(c => c.Group))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(Cases.Value.Select(c => c.Group).Distinct().Order(), Is.EqualTo(new[] { "frag-6mer", "frag-6mermulti", "frag-codon", "frag-hybrid", "frag-local", "frag-multi", "frag-opt", "frag-treein" }));
            Assert.That(Cases.Value.Count(c => c.Rc == 0 && c.Args.Contains("--addtotop")), Is.GreaterThanOrEqualTo(150), "addonetip2top");
            Assert.That(Cases.Value.Count(c => c.Rc == 0 && c.Args.Contains("--keeplength")), Is.GreaterThanOrEqualTo(60), "--keeplength");
            Assert.That(Cases.Value.Count(c => c.Rc == 0 && c.Args.Contains("--codonpos")), Is.GreaterThanOrEqualTo(150), "--codonpos");
            Assert.That(Cases.Value.Count(c => c.Rc == 0 && c.Args.Contains("--codonscore")), Is.GreaterThanOrEqualTo(40), "--codonscore");
            Assert.That(Cases.Value.Count(c => c.OutFiles.Keys.Any(k => k.EndsWith(".tree", StringComparison.Ordinal))), Is.GreaterThanOrEqualTo(50), "--treeout");
            Assert.That(Cases.Value.Count(c => c.Args.Contains("--treein") && c.Rc == 0), Is.GreaterThanOrEqualTo(150), "loadtop");
            Assert.That(Cases.Value.Count(c => c.Rc != 0), Is.GreaterThanOrEqualTo(60), "the stops of the script and of the programs");
        });
    }

    [Test]
    public void SingleReferenceWithTreeOut_TreeFileStartsWithNullPointer()
    {
        // installed mafft 7.505-1: one reference sequence and -t: fixed_musclesupg_double_realloc_nobk_halfmtx_treeout prints tree[-1], a NULL pointer ("(null)" of glibc)
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">a0\nMKVLAAGIVGLLLAQPVSA\n",
            ["new.fa"] = ">f0\nAAGIVGLLLAQ\n>f1\nKVLAAGIAGLL\n",
        };
        var r = MafftProgram.Run(["--quiet", "--addfragments", "new.fa", "--6merpair", "--treeout", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.EqualTo(0));
        Assert.That(r.Files[SequenceAligner_MafftProgram_Tests.InputName + ".tree"], Does.StartWith("(null);\n\n       2: f0\n          nearest sequence: 1\n"));
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
}
