using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// R14-1 (B09 F288): the EMBOSS formats of the core aligners (global_align / local_align / semi_global_align) take the residues before and after each row from the ROWS, not from an end coordinate whose
/// convention differed per producer. Truth: the installed EMBOSS 6.6.0 `water` (same alignment ACGT / ACGT of a = TTACGTTT, b = ACGT) and the fitting rows of the query inside the reference.
/// water -aformat3 markx10 prints ">>>a, 8 nt vs b, 4 nt", al_start 3 al_stop 6 / 1 4; pair "a  3 ACGT  6" / "b  1 ACGT  4"; sam "b 0 a 3 255 4M".
/// </summary>
[TestFixture]
public class AlignmentResultFormatOffsetsTests
{
    private const string A = "TTACGTTT";
    private const string B = "ACGT";
    private static readonly string[] Names = ["a", "b"];

    [Test]
    public void LocalAlign_Markx10_HasTheLengthsAndCoordinatesOfWater()
    {
        string text = AlignmentTools.LocalAlign(A, B, format: "emboss:markx10", names: Names).Formatted!;
        Assert.That(text, Does.Contain(">>>a, 8 nt vs b, 4 nt"));
        Assert.That(text, Does.Contain(">a ..\n; sq_len: 8\n; sq_type: D\n; al_start: 3\n; al_stop: 6\n; al_display_start: 3\nACGT\n"));
        Assert.That(text, Does.Contain(">b ..\n; sq_len: 4\n; sq_type: D\n; al_start: 1\n; al_stop: 4\n; al_display_start: 1\nACGT\n"));
    }

    [Test]
    public void LocalAlign_Pair_HasTheCoordinatesOfWater()
    {
        string text = AlignmentTools.LocalAlign(A, B, format: "emboss:pair", names: Names).Formatted!;
        Assert.That(text, Does.Contain("a                  3 ACGT      6\n"));
        Assert.That(text, Does.Contain("b                  1 ACGT      4\n"));
    }

    [Test]
    public void LocalAlign_Sam_HasTheStartOfWater()
    {
        string text = AlignmentTools.LocalAlign(A, B, format: "emboss:sam", names: Names).Formatted!;
        Assert.That(text, Does.Contain("b\t0\ta\t3\t255\t4M\t"));
    }

    [Test]
    public void LocalAlign_NoPositiveRegion_HasNoNegativeOffset()
    {
        string text = AlignmentTools.LocalAlign("AAAA", "CCCC", format: "emboss:markx10", names: Names).Formatted!;
        Assert.That(text, Does.Contain(">>>a, 4 nt vs b, 4 nt"));
        Assert.That(text, Does.Not.Contain("-1"));
    }

    [Test]
    public void SemiGlobalAlign_Pair_ListsTheWholeReferenceFromOne()
    {
        // fitting: the rows carry the whole reference (free flanks as gaps in the query), so the reference runs 1 .. 8, not 3 .. 10
        string text = AlignmentTools.SemiGlobalAlign(B, A, format: "emboss:pair", names: ["q", "ref"]).Formatted!;
        Assert.That(text, Does.Contain("ref                1 TTACGTTT      8\n"));
        Assert.That(text, Does.Contain("q                  1 --ACGT--      4\n"));
    }

    [Test]
    public void SemiGlobalAlign_Markx10_HasTheWholeLengths()
    {
        string text = AlignmentTools.SemiGlobalAlign(B, A, format: "emboss:markx10", names: ["q", "ref"]).Formatted!;
        Assert.That(text, Does.Contain(">>>q, 4 nt vs ref, 8 nt"));
        Assert.That(text, Does.Contain("; al_start: 1\n; al_stop: 8\n"));
    }

    [Test]
    public void GlobalAlign_Pair_Needle_Water_Stretcher_AgreeOnTheFullRows()
    {
        foreach (string text in new[]
        {
            AlignmentTools.GlobalAlign(A, B, format: "emboss:pair", names: Names).Formatted!,
            AlignmentTools.NeedleAlign(A, B, format: "pair", names: Names).Formatted!,
            AlignmentTools.StretcherAlign(A, B, format: "pair", names: Names).Formatted!,
        })
        {
            Assert.That(text, Does.Contain("a                  1 "), text);
            Assert.That(text, Does.Contain("      8\n"), text);
            Assert.That(text, Does.Contain("b                  1 ").And.Contain("      4\n"), text);
        }
    }

    [Test]
    public void DualAffineAlign_Local_Markx10_HasTheCoordinatesOfWater()
    {
        string text = AlignmentTools.DualAffineAlign(A, B, mode: "local", format: "emboss:markx10", names: Names).Formatted!;
        Assert.That(text, Does.Contain(">>>a, 8 nt vs b, 4 nt"));
        Assert.That(text, Does.Contain("; al_start: 3\n; al_stop: 6\n"));
    }
}
