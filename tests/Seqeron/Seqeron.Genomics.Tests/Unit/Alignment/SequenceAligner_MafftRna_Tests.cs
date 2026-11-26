using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-rna, B09 WP83 / WP86, F374 / F377): the RNA options of the script <c>mafft</c> (<c>--mccaskill</c>, <c>--rnaalifold</c>, <c>--scarnapair</c>, <c>--dafs</c> / <c>--dafspair</c>, <c>--contrafold</c>) against the third run of the Docker package
/// <c>handoff/A2-19c-rna</c> (MAFFT 7.505 built with its extensions, <c>/opt/mafft-ext</c>, ViennaRNA 2.5.1, CONTRAfold 2.02, DAFS 0.0.3, linux/amd64; <c>scripts/alignment/convert_rna_returned.py</c>, <c>mafft_rna_oracle.tsv.gz</c>: 150 sets, every one run twice, the repetitions equal after the process numbers are normalised).
/// The oracle is run through <see cref="MafftProgram.Run(IReadOnlyList{string}, MafftInstallation, IReadOnlyDictionary{string, string}?, string?, CancellationToken)"/> with <see cref="MafftInstallation.WithExtensions"/>: return code, standard output and the whole standard error byte for byte.
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftRna_Tests
{
    private static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Oracle = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_rna_oracle.tsv"));

    /// <summary>Groups whose runs the port reproduces: the stops of the script and of the first program that fails, and the 16 complete runs of <c>--contrafold</c> (contrafoldwrap + CONTRAfold 2.02 + the RNA step of tbfast; F377).</summary>
    public static IEnumerable<TestCaseData> StopGroups() =>
        new[] { "dafs-error", "dafspair-stop", "scarnapair-stop", "mccaskill-stop", "contrafold" }.Select(g => new TestCaseData(g).SetArgDisplayNames(g));

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

    /// <summary>
    /// 96 runs of <c>mafft --contrafold</c> and 42 of <c>mafft --mccaskill</c> with a structure weight (<c>--weightr</c>, <c>--weightm</c>, <c>--rop</c>, three pair libraries, output options) against the installed Debian mafft 7.505 whose <c>contrafold</c> (<c>mxscarnamod</c>) is replaced by a program that answers with the made-up
    /// posteriors (pair lines) of <c>canned</c> (<c>scripts/alignment/generate_mafft_contrafold_oracle.py</c>, F377): the real <c>contrafoldwrap</c> and <c>tbfast</c> (the importance matrix of <c>foldrna</c>, <c>Lalignmm_hmout</c>) are the oracle of the port's wrapper, hat4 and RNA step.
    /// </summary>
    private static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> CannedContrafold = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_contrafold_oracle.tsv"));

    private static MafftInstallation CannedInstallation(SequenceAligner_MafftProgram_Tests.Case c)
    {
        using var doc = System.Text.Json.JsonDocument.Parse(c.Files["canned"]);
        var table = doc.RootElement.EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
        if (c.Group.StartsWith("mccaskill", StringComparison.Ordinal))
            return new MafftInstallation("/usr/lib/mafft/lib/mafft", true, false, false) { MxscarnaPredictor = seq => table[seq.ToUpperInvariant()] };
        return new MafftInstallation("/usr/lib/mafft/lib/mafft", false, false, true) { ContrafoldPredictor = (seq, cutoff) => cutoff == 0.01 ? table[seq] : throw new InvalidOperationException("cutoff") };
    }

    [Test]
    public void CannedContrafold_Shape()
    {
        var by = CannedContrafold.Value.GroupBy(r => r.Group).ToDictionary(g => g.Key, g => g.Count());
        Assert.That(by, Is.EquivalentTo(new Dictionary<string, int> { ["contrafold-canned"] = 96, ["contrafold-canned-iterate"] = 6, ["mccaskill-canned"] = 42, ["mccaskill-canned-iterate"] = 3 }));
        Assert.That(CannedContrafold.Value.Select(r => r.Rc), Is.All.EqualTo(0));
    }

    [Test]
    public void CannedContrafold_EqualsMafft()
    {
        var failures = new List<string>();
        foreach (var c in CannedContrafold.Value.Where(x => x.Group is "contrafold-canned" or "mccaskill-canned"))
        {
            var r = MafftProgram.Run([.. c.Args, SequenceAligner_MafftProgram_Tests.InputName], CannedInstallation(c), new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = c.Input });
            if (r.ExitCode != c.Rc || r.Stdout != c.Out || r.Stderr != c.Err)
                failures.Add(c.Id + " (" + string.Join(' ', c.Args) + "): rc " + r.ExitCode + ", standard output " + (r.Stdout == c.Out ? "equal" : "differs") + ", standard error " + (r.Stderr == c.Err ? "equal" : "differs"));
        }
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void CannedContrafold_WithIterativeRefinement_IsRefusedNotGuessed()
    {
        foreach (var c in CannedContrafold.Value.Where(x => x.Group is "contrafold-canned-iterate" or "mccaskill-canned-iterate"))
            Assert.Throws<NotSupportedException>(() => MafftProgram.Run([.. c.Args, SequenceAligner_MafftProgram_Tests.InputName], CannedInstallation(c), new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = c.Input }), c.Id);
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
                ["scarnapair-stop"] = 18, ["dafs-error"] = 8, ["dafspair-stop"] = 8, ["contrafold"] = 16,
            }));
            Assert.That(rows.Where(r => r.Group is "mccaskill" or "rnaalifold" or "contrafold").Select(r => r.Rc), Is.All.EqualTo(0));
            Assert.That(rows.Where(r => r.Group.EndsWith("-stop", StringComparison.Ordinal)).Select(r => r.Rc), Is.All.EqualTo(1));
            Assert.That(rows.Where(r => r.Group.EndsWith("-stop", StringComparison.Ordinal)).Select(r => r.Out), Is.All.Empty);
            Assert.That(rows.Where(r => r.Group is "mccaskill" or "rnaalifold" or "contrafold").Select(r => r.Out), Is.All.Not.Empty);
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

    /// <summary>Run 2 of the package had the extensions only: the script's "Install DAFS" / "Install CONTRAfold" texts (the 32 sets of run 2 stopped here; checked on one set per option).</summary>
    [TestCase("--dafs --localpair", "\n== Install DAFS===========================================================\nDAFS (Sato et al. Journal 2012 issue:page) is required.\nhttp://www.ncrna.org/ \n==========================================================================\n\n")]
    [TestCase("--dafspair", "\n== Install DAFS===========================================================\nDAFS (Sato et al. Journal 2012 issue:page) is required.\nhttp://www.ncrna.org/ \n==========================================================================\n\n")]
    [TestCase("--contrafold --globalpair", "\n== Install CONTRAfold ===================================================\nThis option requires CONTRAfold\n(Do et al. http://contra.stanford.edu/contrafold/).\nThe executable 'contrafold' have to be installed into \nthe /opt/mafft-ext/libexec/mafft directory. \n=========================================================================\n\n")]
    public void ExtensionsOnlyInstallation_StopsAtTheInstallTextOfRun2(string options, string expectedStderr)
    {
        var args = options.Split(' ').Append("in.fa").ToList();
        var r = MafftProgram.Run(args, MafftInstallation.WithExtensionsOnly, new Dictionary<string, string> { ["in.fa"] = ">s0\nCCCGUCGGUUUAGAGGAUACGCUGG\n>s1\nGUCGUCGGUUUCGAUGAACCGCUGG\n" });
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(1));
            Assert.That(r.Stdout, Is.Empty);
            Assert.That(r.Stderr, Is.EqualTo(expectedStderr));
        });
    }

    /// <summary>
    /// <c>--dafs</c> with a protein input: <c>mccaskillwrap -G</c> stops in front of DAFS ("nuc only"), like <c>--mccaskill</c> (derived from the source, <c>mccaskillwrap.c</c> main; the oracle has no protein <c>--dafs</c> set).
    /// </summary>
    [Test]
    public void DafsWithProtein_StopsWithNucOnly()
    {
        var r = MafftProgram.Run(["--dafs", "--localpair", "in.fa"], MafftInstallation.WithExtensions, new Dictionary<string, string> { ["in.fa"] = ">a\nMKVLAAGIVGLLLAQPSWA\n>b\nMKVLAAGIVGLLLSQPSWA\n" });
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(1));
            Assert.That(r.Stderr, Does.EndWith("inputfile = infile\nnuc only\n"));
            Assert.That(r.Stderr, Does.StartWith("nthread = "));
        });
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
