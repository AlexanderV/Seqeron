using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9e4 (B09 F212): the perl script <c>mapview</c> of MUMmer 3.23 (TestData/Alignment/mapview_oracle.jsonl.gz, scripts/alignment/generate_mapview_oracle.py): the FIG files, standard output / error, the fig2dev command lines
/// of <c>-f pdf / ps</c> (a stub fig2dev records the command and the FIG file it is given) and the side file of <c>-g</c>. Rows whose output depends on perl's hash order (genes with equal 5' ends) are recorded and not compared.
/// </summary>
[TestFixture]
public class MummerMapview_Oracle_Tests
{
    [Test]
    public void Mapview_EqualsTheInstalledScript_AllRows()
    {
        var (failures, c) = MummerMapviewOracleCheck.CheckAll();
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(c["mv"], Is.GreaterThanOrEqualTo(150));
        Assert.That(c["mvConvert"], Is.GreaterThanOrEqualTo(30));
        Assert.That(c["mv"] + c["mvConvert"], Is.GreaterThanOrEqualTo(200));
    }

    [Test]
    public void Mapview_DrawsTheAlignmentsAndConvertsWithFig2dev()
    {
        // mapview -n 1 -p demo x.coords (MUMmer 3.23): one alignment of 98 percent identity over 240 of 260 bases
        const string coords = "ref.fa qry.fa\nNUCMER\n\n    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  | [TAGS]\n"
            + "==========================================================================================================\n      11      250  |        1      240  |      240      240  |    98.76  |      260      240  | ref\tqry\n";
        var r = MummerMapview.Run(coords, options: new MummerMapviewOptions { FileCount = 1, Prefix = "demo" });
        Assert.That(r.Figs.Keys, Is.EqualTo(new[] { "demo_0.fig" }));
        Assert.That(r.Figs["demo_0.fig"], Does.StartWith("#FIG 3.2\nLandscape\nCenter\nInches\nLetter  \n100.00\nMultiple\n-2\n1200 2\n2 1 0 10 11 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 0 3050 17 3050\n"));
        Assert.That(r.Figs["demo_0.fig"], Does.Contain("2 1 0 2 4 0 40 0 -1 0.000 0 0 -1 0 0 2\n\t0 3330 17 3330\n"));
        Assert.That(r.Stdout, Is.Empty);
        Assert.That(r.ExitCode, Is.Zero);
        var pdf = MummerMapview.Run(coords, options: new MummerMapviewOptions { FileCount = 1, Prefix = "demo", Format = "pdf", Magnification = 2, Verbose = true });
        Assert.That(pdf.Figs.Keys, Is.EqualTo(new[] { "tmpdemo_0.fig" }));
        Assert.That(pdf.Commands, Is.EqualTo(new[] { "fig2dev -L pdf -x 100 -m 2 -M tmpdemo_0.fig demo_0.pdf" }));
        Assert.That(pdf.Outputs, Is.EqualTo(new[] { "demo_0.pdf" }));
        Assert.That(pdf.Stdout, Is.EqualTo("demo_0.pdf\t range : 0\t259 \n"));
        // the script goes on after wrong bounds; a UTR file without a CDS file stops it
        Assert.That(MummerMapview.Run(coords, options: new MummerMapviewOptions { X1 = 5, X2 = 3 }).Stdout, Does.StartWith("ERROR3 : wrong range coords : x1 >= x2 !!!\n"));
        var noCds = MummerMapview.Run(coords, "other\tx\tsingle-exon\t5\t50\t.\t+\t.\tg1\n", options: new MummerMapviewOptions());
        Assert.That(noCds.ExitCode, Is.EqualTo(1));
        Assert.That(noCds.Stderr, Does.Contain("ERROR in the input files ! The reference seq ID can't be found in GFF files !"));
    }
}
