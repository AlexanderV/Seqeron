using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// B09 WP67 (F308, R14-13): the MCP tool <c>plotcon</c>. Expected values: what the installed EMBOSS 6.6.0 plotcon writes (<c>-graph data</c>), and the first window derived by hand: three sequences MKTAYIAKQRQISFVKSH / MKTAYLAKQRQISFVRSH / MKSAYIAK-RQISFVKSH,
/// winsize 4: the columns 1-4 are M M M (BLOSUM62 5, three pairs), K K K (5), T T S (T/T 5, T/S 1, T/S 1) and A A A (4): each pair is counted for both sequences, so the sums are 2 x 15, 2 x 15, 2 x 7, 2 x 12 = 98, divided by 3 x 2 x 4 x 4 = 96 gives 1.020833 at the middle of the window (point 3).
/// </summary>
[TestFixture]
public class PlotconToolTests
{
    private const string Alignment = ">a\nMKTAYIAKQRQISFVKSH\n>b\nMKTAYLAKQRQISFVRSH\n>c\nMKSAYIAK-RQISFVKSH\n";

    [Test]
    public void Plotcon_WritesWhatTheInstalledProgramWrites()
    {
        var r = AnalysisTools.Plotcon(Alignment);
        Assert.That(r.Text, Is.EqualTo("##2D Plot\n##Title Similarity Plot of Aligned Sequences\n##Graphs 1\n##Number 1\n##Points 18\n##XminA 0.000000 XmaxA 0.000000 YminA 0.000000 YmaxA 0.000000\n##Xmin 0.000000 Xmax 18.000000 Ymin 0.000000 Ymax 1.187500\n##ScaleXmin 0.000000 ScaleXmax 0.000000 ScaleYmin 0.000000 ScaleYmax 0.000000\n##Maintitle \n##Xtitle Relative Residue Position\n##Ytitle Similarity\n1.000000\t0.000000\n2.000000\t0.000000\n3.000000\t1.020833\n4.000000\t1.145833\n5.000000\t1.000000\n6.000000\t1.104167\n7.000000\t1.166667\n8.000000\t0.833333\n9.000000\t0.979167\n10.000000\t1.041667\n11.000000\t0.979167\n12.000000\t1.125000\n13.000000\t1.187500\n14.000000\t1.125000\n15.000000\t1.062500\n16.000000\t1.062500\n17.000000\t1.187500\n18.000000\t0.000000\n##DataObjects\n##Number 0\n##GraphObjects\n##Number 0\n"));
        var scored = AnalysisTools.Plotcon(Alignment, winsize: 3, scorefile: "EBLOSUM30");
        Assert.That(scored.Text, Does.Contain("##Xmin 0.000000 Xmax 18.000000 Ymin 0.000000 Ymax 2.444444\n"));
        Assert.That(scored.Text, Does.Contain("\n11.000000\t2.444444\n"));
    }

    [Test]
    public void Plotcon_OneSequence_PrintsTheNaNOfTheProgram_AndRefusesBadInput()
    {
        // 0 / 0 in single precision has the sign bit set; ajFmt writes it as -000nan
        var r = AnalysisTools.Plotcon(">a\nMKTAYIAK\n");
        Assert.That(r.Text, Does.Contain("##Xmin 0.000000 Xmax 8.000000 Ymin -000nan Ymax -000nan\n").And.Contain("\n1.000000\t-000nan\n"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Plotcon(""));
        var e = Assert.Throws<ArgumentException>(() => AnalysisTools.Plotcon(Alignment, scorefile: "NOSUCHMATRIX"));
        Assert.That(e!.Message, Does.Contain("Unable to read matrix 'NOSUCHMATRIX'"));
        var w = Assert.Throws<ArgumentException>(() => AnalysisTools.Plotcon(Alignment, winsize: 0, qualifiers: "-sbegin x"));
        Assert.That(w!.Message, Does.Contain("Invalid integer value"));
    }
}
