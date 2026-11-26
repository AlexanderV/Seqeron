using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-rna, B09 WP83, F374): the RNA options of the script <c>mafft</c> (<c>--mccaskill</c>, <c>--rnaalifold</c>, <c>--scarnapair</c>, <c>--dafs</c> / <c>--dafspair</c>, <c>--contrafold</c>) against the second run of the Docker package
/// <c>handoff/A2-19c-rna</c> (MAFFT 7.505 built with its extensions, <c>/opt/mafft-ext</c>, ViennaRNA 2.5.1, linux/amd64; <c>scripts/alignment/convert_rna_returned.py</c>, <c>mafft_rna_oracle.tsv.gz</c>: 150 sets, every one run twice, the repetitions equal).
/// The oracle is run through <see cref="MafftProgram.Run(IReadOnlyList{string}, MafftInstallation, IReadOnlyDictionary{string, string}?, string?, CancellationToken)"/> with <see cref="MafftInstallation.WithExtensions"/>: return code, standard output and the whole standard error byte for byte.
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftRna_Tests
{
    private static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Oracle = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_rna_oracle.tsv"));

    /// <summary>Groups whose runs the port reproduces (the stops of the script and of the first program that fails).</summary>
    public static IEnumerable<TestCaseData> StopGroups() =>
        new[] { "dafs-stop", "contrafold-stop", "scarnapair-stop", "mccaskill-stop" }.Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    /// <summary>Groups that need the base-pairing probabilities of mxscarnamod / McCaskill or the profile scores of RNAalifold.</summary>
    public static IEnumerable<TestCaseData> ComputedGroups() =>
        new[] { "mccaskill", "rnaalifold", "rnaalifold-abort" }.Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    private static MafftProgramResult Run(SequenceAligner_MafftProgram_Tests.Case c, MafftInstallation installation)
        => MafftProgram.Run([.. c.Args, SequenceAligner_MafftProgram_Tests.InputName], installation, new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = c.Input });

    /// <summary>100 further <c>--rnaalifold</c> sets (10 option lines) run in the sandbox with the installed mafft 7.505 and RNAalifold 2.5.1, the versions of the archive (the 25 sets of the archive give equal output here: WP83, F374).</summary>
    private static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> SandboxRnaalifold = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_rnaalifold_oracle.tsv"));

    [Test]
    public void SandboxRnaalifoldOracle_Shape_AndTheRunsAreRefusedNotGuessed()
    {
        var rows = SandboxRnaalifold.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.EqualTo(100));
            Assert.That(rows.Select(r => r.Rc), Is.All.EqualTo(0));
            Assert.That(rows.Select(r => r.Out), Is.All.Not.Empty);
            Assert.That(rows.Count(r => System.Text.RegularExpressions.Regex.IsMatch(r.Err, @"^\d+-\d+, \d+\.\d+$", System.Text.RegularExpressions.RegexOptions.Multiline)), Is.GreaterThan(50), "the base pair lists of RNAalifold -p in the standard error");
            Assert.That(rows.Select(r => string.Join(' ', r.Args)).Distinct().Count(), Is.EqualTo(10));
        });
        foreach (var c in rows.Take(10))
            Assert.Throws<NotSupportedException>(() => Run(c, MafftInstallation.WithExtensions), c.Id);
    }

    [Test]
    public void Oracle_Shape()
    {
        var rows = Oracle.Value;
        var by = rows.GroupBy(r => r.Group).ToDictionary(g => g.Key, g => g.Count());
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.EqualTo(150));
            Assert.That(by, Is.EquivalentTo(new Dictionary<string, int>
            {
                ["mccaskill"] = 73, ["mccaskill-stop"] = 1, ["rnaalifold"] = 25, ["rnaalifold-abort"] = 1,
                ["scarnapair-stop"] = 18, ["dafs-stop"] = 16, ["contrafold-stop"] = 16,
            }));
            Assert.That(rows.Where(r => r.Group is "mccaskill" or "rnaalifold").Select(r => r.Rc), Is.All.EqualTo(0));
            Assert.That(rows.Where(r => r.Group.EndsWith("-stop", StringComparison.Ordinal)).Select(r => r.Rc), Is.All.EqualTo(1));
            Assert.That(rows.Where(r => r.Group.EndsWith("-stop", StringComparison.Ordinal)).Select(r => r.Out), Is.All.Empty);
            Assert.That(rows.Where(r => r.Group is "mccaskill" or "rnaalifold").Select(r => r.Out), Is.All.Not.Empty);
            Assert.That(rows.Select(r => r.Err.Length), Is.All.GreaterThan(0));
        });
    }

    [TestCaseSource(nameof(StopGroups))]
    public void StopsOfTheExtendedInstallation_EqualMafft(string group)
    {
        var failures = new List<string>();
        int n = 0;
        foreach (var c in Oracle.Value.Where(x => x.Group == group))
        {
            n++;
            MafftProgramResult r;
            try
            {
                r = Run(c, MafftInstallation.WithExtensions);
            }
            catch (Exception e)
            {
                failures.Add(c.Id + " (" + string.Join(' ', c.Args) + "): " + e.GetType().Name + " " + e.Message);
                continue;
            }
            if (r.ExitCode != c.Rc || r.Stdout != c.Out || r.Stderr != c.Err)
                failures.Add(c.Id + " (" + string.Join(' ', c.Args) + "): rc " + r.ExitCode + " / " + c.Rc + ", standard error " + (r.Stderr == c.Err ? "equal" : "differs: " + r.Stderr.Replace("\n", "|")));
        }
        Assert.That(n, Is.GreaterThan(0));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    /// <summary>
    /// The runs that need the computed base-pairing probabilities are not answered by the port (it throws instead of printing a different alignment): the day the port computes them, this test is replaced by the byte-for-byte comparison of the group.
    /// </summary>
    [TestCaseSource(nameof(ComputedGroups))]
    public void ComputedRuns_AreRefusedNotGuessed(string group)
    {
        foreach (var c in Oracle.Value.Where(x => x.Group == group))
            Assert.Throws<NotSupportedException>(() => Run(c, MafftInstallation.WithExtensions), c.Id);
    }

    [TestCase("--dafs --localpair")]
    [TestCase("--dafspair")]
    [TestCase("--contrafold --globalpair")]
    [TestCase("--mccaskill --globalpair")]
    [TestCase("--scarnapair")]
    public void DebianInstallation_RefusesEveryRnaOption(string options)
    {
        var args = options.Split(' ').Append("in.fa").ToList();
        var r = MafftProgram.Run(args, new Dictionary<string, string> { ["in.fa"] = ">a\nACGUACGUACGUAC\n>b\nACGUACGGACGUAC\n" });
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(1));
            Assert.That(r.Stderr, Does.Contain("== Install "));
            Assert.That(r.Stderr, Does.Not.Contain("/opt/mafft-ext"));
        });
    }
}
