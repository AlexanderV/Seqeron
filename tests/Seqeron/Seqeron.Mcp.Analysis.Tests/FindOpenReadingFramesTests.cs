using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_open_reading_frames</c> MCP tool.
/// Expected values from the ATG-to-first-in-frame-stop definition (length divisible
/// by 3, frame 1-3, codonCount = length/3), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindOpenReadingFramesTests
{
    [Test]
    public void FindOpenReadingFrames_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindOpenReadingFrames("ATGTAA", 6));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames("", 6));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames(null!, 6));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames("XYZ", 6));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindOpenReadingFrames("ATGTAA", 0));
    }

    [Test]
    public void FindOpenReadingFrames_Binding_InvokesSuccessfully()
    {
        // Minimal forward ORF ATG-TAA (RC = TTACAT has no ATG).
        var minimal = AnalysisTools.FindOpenReadingFrames("ATGTAA", 6).Items;
        Assert.Multiple(() =>
        {
            Assert.That(minimal, Has.Length.EqualTo(1));
            Assert.That(minimal[0].Sequence, Is.EqualTo("ATGTAA"));
            Assert.That(minimal[0].Position, Is.EqualTo(0));
            Assert.That(minimal[0].Frame, Is.EqualTo(1));
            Assert.That(minimal[0].IsReverseComplement, Is.False);
            Assert.That(minimal[0].Length, Is.EqualTo(6));
            Assert.That(minimal[0].CodonCount, Is.EqualTo(2));
        });

        // ATG-AAA-TAG (RC = CTATTTCAT has no ATG).
        var extra = AnalysisTools.FindOpenReadingFrames("ATGAAATAG", 9).Items;
        Assert.Multiple(() =>
        {
            Assert.That(extra, Has.Length.EqualTo(1));
            Assert.That(extra[0].Sequence, Is.EqualTo("ATGAAATAG"));
            Assert.That(extra[0].Length, Is.EqualTo(9));
            Assert.That(extra[0].CodonCount, Is.EqualTo(3));
        });
    }

    // GENOMIC-ORF-001 (B09 F21-F26): options and coordinates. Expected values from the independent Python oracle
    // (scripts/analysis/generate_orf_oracle.py: Biopython finder == NCBI COrf port == EMBOSS getorf).

    [Test]
    public void FindOpenReadingFrames_InputStrandCoordinates_AndTerminationFlag()
    {
        var fwd = AnalysisTools.FindOpenReadingFrames("ATGTAA", 6).Items.Single();
        Assert.Multiple(() =>
        {
            Assert.That(fwd.Start, Is.EqualTo(0));
            Assert.That(fwd.End, Is.EqualTo(5));
            Assert.That(fwd.IsTerminated, Is.True);
        });

        // CC + TTAGGGCAT + GG: reverse complement CC ATGCCCTAA GG, ORF at scanned-strand position 2, input span 2..10.
        var rev = AnalysisTools.FindOpenReadingFrames("CCTTAGGGCATGG", 9).Items.Single();
        Assert.Multiple(() =>
        {
            Assert.That(rev.IsReverseComplement, Is.True);
            Assert.That(rev.Sequence, Is.EqualTo("ATGCCCTAA"));
            Assert.That(rev.Position, Is.EqualTo(2));
            Assert.That(rev.Start, Is.EqualTo(2));
            Assert.That(rev.End, Is.EqualTo(10));
        });
    }

    [Test]
    public void FindOpenReadingFrames_GeneticCode_ChoosesTheStopCodons()
    {
        // Table 2 (vertebrate mitochondrial): TGA = Trp (reads through), AGA = stop. Biopython 1.88 CodonTable.
        Assert.That(AnalysisTools.FindOpenReadingFrames("ATGTGAAAATAA", 1).Items.Single().Sequence, Is.EqualTo("ATGTGA"));
        Assert.That(AnalysisTools.FindOpenReadingFrames("ATGTGAAAATAA", 1, genetic_code: 2).Items.Single().Sequence, Is.EqualTo("ATGTGAAAATAA"));
        Assert.That(AnalysisTools.FindOpenReadingFrames("ATGCCCAGAGGG", 1).Items, Is.Empty);
        Assert.That(AnalysisTools.FindOpenReadingFrames("ATGCCCAGAGGG", 1, genetic_code: 2).Items.Single().Sequence, Is.EqualTo("ATGCCCAGA"));

        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames("ATGTAA", 1, genetic_code: 27), "dual-coding stops");
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames("ATGTAA", 1, genetic_code: 99), "unknown table");
    }

    [Test]
    public void FindOpenReadingFrames_StartCodonModes()
    {
        Assert.That(AnalysisTools.FindOpenReadingFrames("TTGAAATAA", 1).Items, Is.Empty);
        Assert.That(AnalysisTools.FindOpenReadingFrames("TTGAAATAA", 1, startCodons: "geneticcodestarts").Items.Single().Sequence, Is.EqualTo("TTGAAATAA"));
        Assert.That(AnalysisTools.FindOpenReadingFrames("GTGAAATAA", 1, startCodons: "GeneticCodeStarts").Items, Is.Empty, "GTG is not a table-1 start");
        Assert.That(AnalysisTools.FindOpenReadingFrames("GTGAAATAA", 1, genetic_code: 11, startCodons: "GeneticCodeStarts").Items.Single().Sequence, Is.EqualTo("GTGAAATAA"));
        // CCC AAA TAA GGG: stop-to-stop region 0..8 (any sense codon, longest per stop).
        var any = AnalysisTools.FindOpenReadingFrames("CCCAAATAAGGG", 1, startCodons: "AnySenseCodon", nested: false).Items.Single();
        Assert.That((any.Start, any.End), Is.EqualTo((0, 8)));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames("ATGTAA", 1, startCodons: "Everything"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames("ATGTAA", 1, startCodons: "1"), "numeric strings are not enum names");
    }

    [Test]
    public void FindOpenReadingFrames_NestedOpenEndedLengthConventionAndStrands()
    {
        Assert.That(AnalysisTools.FindOpenReadingFrames("ATGGGGATGCCCTAA", 1).Items.Count(o => !o.IsReverseComplement), Is.EqualTo(2));
        var outer = AnalysisTools.FindOpenReadingFrames("ATGGGGATGCCCTAA", 1, nested: false).Items.Single(o => !o.IsReverseComplement);
        Assert.That(outer.Position, Is.EqualTo(0));

        Assert.That(AnalysisTools.FindOpenReadingFrames("ATGAAACCCGG", 1).Items, Is.Empty);
        var open = AnalysisTools.FindOpenReadingFrames("ATGAAACCCGG", 1, includeOpenEnded: true).Items.Single();
        Assert.Multiple(() =>
        {
            Assert.That(open.Sequence, Is.EqualTo("ATGAAACCC"));
            Assert.That(open.IsTerminated, Is.False);
            Assert.That((open.Start, open.End), Is.EqualTo((0, 8)));
        });

        // ATGAAATAA: 9 nt with stop, 6 nt without (NCBI ORFfinder -ml convention).
        Assert.That(AnalysisTools.FindOpenReadingFrames("ATGAAATAA", 9).Items, Has.Length.EqualTo(1));
        Assert.That(AnalysisTools.FindOpenReadingFrames("ATGAAATAA", 9, minLengthExcludesStop: true).Items, Is.Empty);
        Assert.That(AnalysisTools.FindOpenReadingFrames("ATGAAATAA", 6, minLengthExcludesStop: true).Items, Has.Length.EqualTo(1));

        Assert.That(AnalysisTools.FindOpenReadingFrames("TTAGGGGGGCAT", 1).Items, Has.Length.EqualTo(1));
        Assert.That(AnalysisTools.FindOpenReadingFrames("TTAGGGGGGCAT", 1, bothStrands: false).Items, Is.Empty);
    }
}
