using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9e4 (B09 F213): the perl script <c>mummerplot</c> of MUMmer 3.23 (TestData/Alignment/mummerplot_oracle.jsonl.gz, scripts/alignment/generate_mummerplot_oracle.py): the data files, the gnuplot script, the filtered delta file, the messages and the
/// arguments of the gnuplot call (a stub gnuplot records them; the real gnuplot 6.0 was run by hand on samples: png / ps come out). Rows whose layout depends on perl's hash order (sequences left out of a layout, equal lengths) are
/// recorded and not compared; the exit status of a die is the stale errno of the last system call (25, 255 ...): only zero / non-zero is compared.
/// </summary>
[TestFixture]
public class MummerPlot_Oracle_Tests
{
    [Test]
    public void Mummerplot_EqualsTheInstalledScript_AllRows()
    {
        var (failures, c) = MummerPlotOracleCheck.CheckAll();
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(c["mp"], Is.GreaterThanOrEqualTo(200));
        Assert.That(c["mpDies"], Is.GreaterThanOrEqualTo(3));
    }

    private const string Delta = "ref.fa qry.fa\nNUCMER\n>ref qry 260 240\n11 250 1 240 3 3 0\n131\n-51\n0\n";

    [Test]
    public void Mummerplot_WritesTheDataFilesAndTheGnuplotScript()
    {
        // mummerplot -p demo --png x.delta (MUMmer 3.23 with gnuplot 6.0): one alignment of 98.76 percent similarity
        var files = new Dictionary<string, string> { ["ref.fa"] = ">ref\n" + new string('a', 260) + "\n", ["qry.fa"] = ">qry\n" + new string('c', 240) + "\n" };
        var r = MummerPlot.Run(Delta, files, new MummerPlotOptions { Prefix = "demo", Png = true, MatchName = "x.delta" });
        Assert.That(r.Files["demo.fplot"], Is.EqualTo("#-- forward hits sorted by %sim\n0 0 0\n0 0 0\n\n\n11 1 98.7551867219917\n250 240 98.7551867219917\n\n\n"));
        Assert.That(r.Files["demo.rplot"], Is.EqualTo("#-- reverse hits sorted by %sim\n0 0 0\n0 0 0\n\n\n"));
        Assert.That(r.Files["demo.gp"], Is.EqualTo("set terminal png tiny size 800,800\nset output \"demo.png\"\nset size 1,1\nset grid\nunset key\nset border 15\nset tics scale 0\nset xlabel \"ref\"\nset ylabel \"qry\"\n"
            + "set format \"%.0f\"\nset mouse format \"%.0f\"\nset mouse mouseformat \"[%.0f, %.0f]\"\nset xrange [1:260]\nset yrange [1:240]\nset style line 1  lt 1 lw 3 pt 6 ps 1\nset style line 2  lt 3 lw 3 pt 6 ps 1\n"
            + "set style line 3  lt 2 lw 3 pt 6 ps 1\nplot \\\n \"demo.fplot\" title \"FWD\" w lp ls 1, \\\n \"demo.rplot\" title \"REV\" w lp ls 2\n"));
        Assert.That(r.Stdout, Is.EqualTo("gnuplot 6.0 patchlevel 0\n"));
        Assert.That(r.Stderr, Is.EqualTo("Reading delta file x.delta\nWriting plot files demo.fplot, demo.rplot\nWriting gnuplot script demo.gp\nRendering plot demo.png\n"));
        Assert.That(r.GnuplotCommand, Is.EqualTo("gnuplot demo.gp"));
        Assert.That(r.ExitCode, Is.Zero);
    }

    [Test]
    public void Mummerplot_RefusesWhatTheScriptDiesOn_AndWarnsLikeTheScript()
    {
        var files = new Dictionary<string, string> { ["ref.fa"] = ">ref\n" + new string('a', 260) + "\n", ["qry.fa"] = ">qry\n" + new string('c', 240) + "\n" };
        var bad = MummerPlot.Run(Delta, files, new MummerPlotOptions { Terminal = "foo" });
        Assert.That(bad.ExitCode, Is.Not.Zero);
        Assert.That(bad.Stderr, Is.EqualTo("ERROR: Invalid terminal type, foo\n"));
        Assert.That(MummerPlot.Run("foo bar\nbaz\n", files).Stderr, Is.EqualTo("ERROR: Could not read out.delta, Unrecognized file type\n"));
        Assert.That(MummerPlot.Run("", files).Stderr, Is.EqualTo("ERROR: Could not read out.delta, File is empty\n"));
        var multi = MummerPlot.Run(Delta, files, new MummerPlotOptions { Png = true, Postscript = true, Large = true, Small = true });
        Assert.That(multi.Stderr, Does.StartWith("WARNING: Multiple terminals not allowed, using 'postscript'\nWARNING: Multiple sizes now allowed, using 'small'\n"));
        // -f / -l / -S need a delta file: a cluster file warns and ignores them
        var cluster = MummerPlot.Run("ref.fa qry.fa\nNUCMER\n>ref qry 260 240\n 1  1\n      11        1     60     -      -\n      72       62     69     1      1\n", files, new MummerPlotOptions { Filter = true, MatchName = "x.cluster" });
        Assert.That(cluster.Stderr, Does.StartWith("WARNING: -f -l -S only work with delta input\nReading cluster file x.cluster\n"));
    }
}
