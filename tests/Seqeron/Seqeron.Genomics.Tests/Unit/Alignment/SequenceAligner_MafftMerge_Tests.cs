using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s5, B09 WP74, F340): <c>mafft --merge table</c> against the installed Debian mafft 7.505-1 (<c>scripts/alignment/generate_mafft_wp74_oracle.py</c> merge, merge-pair, merge-iter, merge-err, merge-tree):
/// return code, standard output, the files (<c>GuideTree</c>, <c>input.tree</c>, <c>input.hat2</c>) and the whole standard error of every run.
/// <list type="bullet">
/// <item><c>merge</c> - <c>disttbfast -H 0</c>: <c>readsubalignmentstable</c>, <c>increaseintergroupdistanceshalfmtx</c>, the constrained UPGMA tree
/// (<c>fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained</c>: a join that would break a subalignment is marked inconsistent), the steps inside a subalignment are not aligned (<c>mergeoralign</c> 'n'), the rows of a subalignment keep their columns.</item>
/// <item><c>merge-pair</c> - <c>tbfast -H 0</c> with a pair library (the same tree and steps, the error boxes of <c>tbfast</c>); <c>merge-iter</c> - <c>dvtditr -H 0</c>
/// (<c>fixed_supg_double_treeout_constrained</c>, the branches inside a subalignment are not refined); <c>merge-tree</c> - with <c>--treein</c> (no constrained tree, the check of the groups).</item>
/// <item><c>merge-err</c> - the messages and stops: a table line that is not numbers, a sequence in two groups, a number above the number of sequences, unaligned subalignments ("must be aligned"), a subalignment that is no group of the
/// guide tree ("does not seem to form a monophyletic cluster": disttbfast stops with exit status 0 and the script copies the file <c>GuideTree</c>; tbfast and dvtditr stop with 1), comment-only and empty tables.</item>
/// </list>
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftMerge_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp74_merge_oracle.tsv"));

    public static IEnumerable<TestCaseData> Groups() => Cases.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        Assert.Multiple(() =>
        {
            foreach (var g in Cases.Value.GroupBy(c => c.Group))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(Cases.Value.Select(c => c.Group).Distinct().Order(), Is.EqualTo(new[] { "merge", "merge-err", "merge-iter", "merge-pair", "merge-tree" }));
            // the standard error is stored in full for every 25th case and for the short ones only: the counts below are of the stored texts, of the files and of the arguments
            var merge = Cases.Value.Where(c => c.Group == "merge").ToList();
            Assert.That(merge.Count(c => c.Rc == 0 && c.Out.Length > 0), Is.GreaterThanOrEqualTo(300), "an alignment");
            Assert.That(merge.Count(c => c.Files.Values.Any(v => v.Split('\n').Any(l => l.StartsWith('-') || l.Contains(" -", StringComparison.Ordinal) || l.Contains("\t-", StringComparison.Ordinal)))), Is.GreaterThanOrEqualTo(30), "a negative number: the columns of gaps are kept");
            Assert.That(merge.Count(c => c.OutFiles.Keys.Any(k => k.EndsWith(".tree", StringComparison.Ordinal))), Is.GreaterThanOrEqualTo(20), "--treeout (the file infile.tree of the constrained tree)");
            Assert.That(merge.Count(c => c.OutFiles.ContainsKey("GuideTree")), Is.GreaterThanOrEqualTo(5), "a subalignment that is no group of the tree: disttbfast stops with status 0, the script copies GuideTree");
            Assert.That(merge.Count(c => c.Err.Contains("Checking subalignment", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(8));
            var pair = Cases.Value.Where(c => c.Group == "merge-pair").ToList();
            Assert.That(pair.Count(c => c.Args.Contains("--localpair") || c.Args.Contains("--globalpair") || c.Args.Contains("--genafpair") || c.Args.Contains("--oldgenafpair")), Is.GreaterThanOrEqualTo(200), "tbfast");
            var iter = Cases.Value.Where(c => c.Group == "merge-iter").ToList();
            Assert.That(iter.Count(c => c.Err.Contains("dvtditr", StringComparison.Ordinal) && c.Err.Contains("Checking subalignment", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(5), "dvtditr checks the subalignments");
            Assert.That(iter.Count(c => c.OutFiles.Keys.Any(k => k.EndsWith(".hat2", StringComparison.Ordinal))), Is.GreaterThanOrEqualTo(10), "--distout");
            var err = Cases.Value.Where(c => c.Group == "merge-err").ToList();
            Assert.That(err.Count(c => c.Rc != 0), Is.GreaterThanOrEqualTo(100), "the stops of the programs");
            Assert.That(err.Count(c => c.Err.Contains("must be aligned", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(3), "unaligned subalignment");
            Assert.That(err.Count(c => c.Err.Contains("appears in different groups", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(30), "a sequence in two groups");
            Assert.That(err.Count(c => c.Err.Contains("does not exist in the input sequence file", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(10), "a number above the number of sequences");
            Assert.That(err.Count(c => c.Err.Contains("Format error?", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(3), "a table line that is not numbers");
            Assert.That(err.Count(c => c.Err.Contains("forced to be a monophyletic cluster", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(2), "the constrained tree had to break a join");
            Assert.That(err.Count(c => c.Err.Contains("At least 2 sequences should be input!", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(1), "dndpre stops on the empty file of the stopped disttbfast");
            var tree = Cases.Value.Where(c => c.Group == "merge-tree").ToList();
            Assert.That(tree.Count(c => c.Rc == 0), Is.GreaterThanOrEqualTo(200), "--merge with --treein");
            Assert.That(tree.Count(c => c.OutFiles.ContainsKey("GuideTree")), Is.GreaterThanOrEqualTo(10), "a subalignment that is no group of the given tree");
            Assert.That(tree.Count(c => c.Err.Contains("cp: cannot stat 'infile.tree'", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(1), "the copy of the tree file fails when no program wrote it");
            Assert.That(tree.Count(c => c.Err.Contains("does not seem to form a monophyletic cluster", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(1), "disttbfast / dvtditr");
            Assert.That(tree.Count(c => c.Err.Contains("does not form a monophyletic cluster", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(1), "tbfast");
        });
    }

    [Test]
    public void TwoSubalignments_KeepTheirColumns()
    {
        // installed mafft 7.505-1: "mafft --quiet --merge t.txt in.fa" - sequences 1 + 2 and 3 + 4 are subalignments (the rows of each keep their columns: "MKV-LAAGIVG" / "MKVALAAGIVG" stay one alignment), the two groups are aligned to each other
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nMKV-LAAGIVG\n>b\nMKVALAAGIVG\n>c\nMRVLLGIVGLL\n>d\nMRV-LGIVGLL\n",
            ["t.txt"] = "1 2\n3 4\n",
        };
        var r = MafftProgram.Run(["--quiet", "--merge", "t.txt", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.Zero);
        Assert.That(r.Stdout, Is.EqualTo(">a\nMKV-LAAGIVG--\n>b\nMKVALAAGIVG--\n>c\nMRV--LLGIVGLL\n>d\nMRV---LGIVGLL\n"));
        Assert.That(r.Stderr, Is.Empty);
    }

    [Test]
    public void SubalignmentThatIsNoGroupOfTheTree_StopsDisttbfastWithStatusZeroAndLeavesTheGuideTree()
    {
        // installed mafft 7.505-1: sequences 1 and 4 are far from each other; the constrained tree joins them first, the checking of the groups of the progressive stage finds no group {1, 3} - here: a table for --treein
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nMKVLAAGIVG\n>b\nMKVLAGIVG\n>c\nMKVIAAGLVG\n",
            ["t.txt"] = "1 3\n",
            ["tree.txt"] = "1 2 0.1 0.1\n1 3 0.1 0.1\n",
        };
        var r = MafftProgram.Run(["--merge", "t.txt", "--treein", "tree.txt", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.Zero, "disttbfast stops after its message (goto chudan) with exit status 0 ...");
        Assert.That(r.Stdout, Is.Empty, "... and writes no alignment");
        Assert.That(r.Stderr, Does.Contain("cp: cannot stat 'infile.tree': No such file or directory"));
        Assert.That(r.Stderr, Does.Contain("# Subalignment 1 does not seem to form a monophyletic cluster"));
        Assert.That(r.Files.ContainsKey("GuideTree"), Is.False, "the tree file was not written (--treein without --treeout): cp fails");
    }

    [Test]
    public void TableNumberFarBehindTheBlock_IsRefused()
    {
        // valgrind on the installed disttbfast: "Invalid read of size 4 ... 16 bytes after a block of size 16 alloc'd" (readsubalignmentstable reads tab01[p] before it tests p); the byte behind the padding of the block belongs
        // to the next heap chunk (a size field, then data of other allocations): what the program prints there depends on the heap layout
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nMKVLAAGIVG\n>b\nMKVLAGIVG\n>c\nMKVIAAGLVG\n>d\nMRVLAAGIVG\n",
            ["t.txt"] = "1 2\n1000 1\n",
        };
        Assert.Throws<NotSupportedException>(() => MafftProgram.Run(["--merge", "t.txt", SequenceAligner_MafftProgram_Tests.InputName], files));
        // a number inside the padding of the block (calloc( 4, 4 ) = 16 bytes in a 32-byte chunk: slots 4 .. 5) reads zero and gives the message of the program
        files["t.txt"] = "1 2\n5 1\n";
        var r = MafftProgram.Run(["--merge", "t.txt", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.EqualTo(1));
        Assert.That(r.Stderr, Does.EndWith("Sequence 5 does not exist in the input sequence file.\n"));
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
