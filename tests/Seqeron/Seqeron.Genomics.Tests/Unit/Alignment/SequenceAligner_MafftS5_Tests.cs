using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s5, B09 WP72, F330 ff.): the options of the script <c>mafft</c> that needed their own programs - <c>--skipiterate</c> (<c>dvtditr -E</c>: <c>generatesubalignmentstable</c>), <c>--core</c> (<c>setcore</c>), <c>--aamatrix</c> (<c>loadaamtx</c>), <c>--add</c> / <c>--addprofile</c> with <c>--maxiterate</c> (<c>dvtditr -K</c>) and <c>--addprofile</c> alone -
/// against the installed Debian mafft 7.505-1 (<c>scripts/alignment/generate_mafft_s5_oracle.py</c>): return code, standard output, the files and the whole standard error of every run.
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftS5_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_s5_oracle.tsv"));

    public static IEnumerable<TestCaseData> Groups() => Cases.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        Assert.Multiple(() =>
        {
            foreach (var g in Cases.Value.GroupBy(c => c.Group))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            var skip = Cases.Value.Where(c => c.Group == "skipiterate").ToList();
            Assert.That(skip.Count(c => c.Err.Contains("Iterative refinment was not done", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(30), "the value is above every branch: nothing is refined");
            Assert.That(skip.Count(c => c.Err.Contains("Checking subalignment", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(30), "subalignments are fixed");
            var core = Cases.Value.Where(c => c.Group == "core").ToList();
            Assert.That(core.Count(c => c.Err.Contains("Make hat2.", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(30), "no file hat2: an empty output");
            Assert.That(core.Count(c => c.Out.Length > 0), Is.GreaterThanOrEqualTo(100), "a core is written");
            var aa = Cases.Value.Where(c => c.Group == "aamatrix").ToList();
            Assert.That(aa.Count(c => c.Err.Contains("no rescale", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(20), "norescale");
            Assert.That(aa.Count(c => c.Err.Contains("not found in the first 20 letters", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(5), "a letter is missing in the first line");
            Assert.That(aa.Count(c => c.Err.Contains("Format error in aa matrix", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(3), "a row has too few numbers");
            var addIter = Cases.Value.Where(c => c.Group == "additer").ToList();
            var addProfile = Cases.Value.Where(c => c.Group == "addprofile").ToList();
            Assert.That(addIter.Count(c => c.Rc == 0 && c.Out.Length > 0), Is.GreaterThanOrEqualTo(150), "additer: an alignment");
            Assert.That(addProfile.Count(c => c.Rc == 0 && c.Out.Length > 0), Is.GreaterThanOrEqualTo(60), "addprofile: an alignment");
            Assert.That(addProfile.Count(c => c.Err.Contains("There is no appropriate position to add", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(50), "addprofile: the added sequences are no group of the guide tree");
            Assert.That(addIter.Count(c => c.Err.Contains("There is no appropriate position to add", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(1), "additer: the same stop");
            Assert.That(Cases.Value.Count(c => c.Group == "addprofile" && c.Err.Contains("must be aligned", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(5), "the added sequences are not an alignment");
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
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }
}
