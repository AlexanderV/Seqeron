using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s5, B09 WP76, F355 / F356): <c>mafft --adjustdirection</c> / <c>--adjustdirectionaccurately</c> (the programs <c>makedirectionlist -m</c> and <c>setdirection</c>) and <c>mafft --anchors</c> (<c>readexternalanchors</c>, <c>pickpairanch</c>, <c>Falign_givenanchors</c>)
/// against the installed Debian mafft 7.505-1 (<c>scripts/alignment/generate_mafft_wp76_oracle.py</c> dir, dir-add, dir-merge, anchors): return code, standard output, the files and the whole standard error of every run.
/// <list type="bullet">
/// <item><c>makedirectionlist</c> compares every sequence (every new one with <c>-I n</c>) with the sequences before it (in the order of the "contrast": the sequences whose own 6-mers differ most from those of their reverse complement come first):
/// the mean number of shared 6-mers (<c>--adjustdirection</c>, mode a) or the local alignment score (<c>--adjustdirectionaccurately</c>, <c>-d</c>) of the sequence and of its reverse complement with the reference sequences decides between _F_ and _R_;
/// if the first sequence ends up _R_ all labels are turned. <c>setdirection</c> reverse complements the _R_ sequences and puts _R_ in front of their titles (with <c>--anysymbol</c> the sequences of <c>orig</c> are reversed too, their titles stay);
/// with <c>--merge</c> the sequences of a subalignment must agree. A protein input is not changed ("Not necessary!").</item>
/// </list>
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftDirection_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp76_oracle.tsv"));

    public static IEnumerable<TestCaseData> Groups() => Cases.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        Assert.Multiple(() =>
        {
            foreach (var g in Cases.Value.GroupBy(c => c.Group))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            var dir = Cases.Value.Where(c => c.Group.StartsWith("dir", StringComparison.Ordinal)).ToList();
            Assert.That(dir.Count(c => c.Rc == 0 && c.Out.Contains(">_R_", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(250), "reverse complemented sequences with _R_ in their titles");
            Assert.That(dir.Count(c => c.Args.Contains("--adjustdirectionaccurately") && c.Rc == 0), Is.GreaterThanOrEqualTo(100), "the local alignment score");
            Assert.That(dir.Count(c => c.Args.Contains("--adjustdirection") && c.Rc == 0), Is.GreaterThanOrEqualTo(100), "the 6-mers");
            Assert.That(dir.Count(c => c.Args.Contains("--anysymbol") && c.Rc == 0), Is.GreaterThanOrEqualTo(15), "orig");
            Assert.That(Cases.Value.Count(c => c.Group == "dir-merge" && c.Rc == 0 && c.Out.Contains(">_R_", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(40), "--merge with consistent directions and turned subalignments");
            Assert.That(Cases.Value.Count(c => c.Group == "dir-merge" && c.Rc != 0), Is.GreaterThanOrEqualTo(40), "the stops of setdirection and of the script");
            var anchors = Cases.Value.Where(c => c.Group == "anchors").ToList();
            Assert.That(anchors.Count(c => c.Rc == 0), Is.GreaterThanOrEqualTo(180), "alignments with external anchors");
            Assert.That(anchors.Count(c => c.Rc != 0 && c.Err.Contains("The sequence id of the first column must be less than the second.", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(3), "the check of readexternalanchors");
            Assert.That(anchors.Count(c => c.Files["anchors.txt"].Split('\n').Length > 8), Is.GreaterThanOrEqualTo(150), "files with many anchors");
        });
    }

    [Test]
    public void DirMerge_TurnedMemberOfASubalignment_StopsWithTheMessageOfSetdirection()
    {
        // the standard error of those runs is stored as a SHA-256 only: the port, which equals the oracle by that hash in Group_EqualsMafft, is asked for the text
        var stopped = Cases.Value.Where(c => c.Group == "dir-merge" && c.Rc != 0).Select(SequenceAligner_MafftProgram_Tests.RunCase)
            .Count(r => r.Stderr.Contains("Direction of nucleotide sequences seems to be inconsistent", StringComparison.Ordinal) && r.Stderr.Contains("but their directions seem to be different.", StringComparison.Ordinal));
        Assert.That(stopped, Is.GreaterThanOrEqualTo(10));
    }

    [Test]
    public void ReverseComplementedMember_GetsRPrefixAndIsAlignedAsForward()
    {
        // installed mafft 7.505-1: mafft --adjustdirection in.fa  (the second sequence is the reverse complement of the first one)
        const string input = ">a\nATGGCGTACGTTAGCCTAGGATCCGATTACAGGCT\n>b\nAGCCTGTAATCGGATCCTAGGCTAACGTACGCCAT\n>c\nATGGCGTACGTAAGCCTAGGATCCGATTACAGGCT\n";
        var r = MafftProgram.Run(["--quiet", "--adjustdirection", SequenceAligner_MafftProgram_Tests.InputName], new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = input });
        Assert.That(r.ExitCode, Is.EqualTo(0));
        Assert.That(r.Stdout, Does.Contain(">_R_b\n"));
        Assert.That(r.Stdout, Does.Contain(">a\n").And.Contain(">c\n"));
        var rows = r.Stdout.Split('>', StringSplitOptions.RemoveEmptyEntries).Select(x => x.Split('\n')[1]).ToArray();
        Assert.That(rows.Select(x => x.Length).Distinct().Count(), Is.EqualTo(1), "an alignment");
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
