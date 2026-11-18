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

    // B09 A1-21: IUPAC input is accepted; NNN is read through, TAR (TAA|TAG) is a stop, an ambiguous codon never starts an ORF
    // (oracle orf_ambiguous_oracle.tsv: Biopython definition + EMBOSS getorf 6.6.0).
    [Test]
    public void FindOpenReadingFrames_IupacInput_AmbiguousCodonRules()
    {
        var through = AnalysisTools.FindOpenReadingFrames("ATGNNNTAA", 1, bothStrands: false).Items.Single();
        Assert.Multiple(() =>
        {
            Assert.That(through.Sequence, Is.EqualTo("ATGNNNTAA"));
            Assert.That((through.Start, through.End), Is.EqualTo((0, 8)));
            Assert.That(AnalysisTools.FindOpenReadingFrames("ATGTARCCCTAA", 1, bothStrands: false).Items.Single().Sequence, Is.EqualTo("ATGTAR"));
            Assert.That(AnalysisTools.FindOpenReadingFrames("NTGCCCTAA", 1, bothStrands: false).Items, Is.Empty);
        });
        var rev = AnalysisTools.FindOpenReadingFrames("TTANNNCAT", 1).Items.Single(i => i.IsReverseComplement);
        Assert.That((rev.Sequence, rev.Start, rev.End), Is.EqualTo(("ATGNNNTAA", 0, 8)));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames("ATGXTAA", 1));
    }

    // B09 A1-18 / A1-19: NCBI ORFfinder 0.4.3 options measured on the binary (orffinder_cli_oracle.tsv, docs/Evidence/GENOMIC-ORF-001-Evidence.md).
    // strand_equal_span_s2_12 (48 nt), -s 2 -ml 30: the plus ORF 7..45 and the minus ORF 7..45 are identical spans; -n true drops the plus one.
    [Test]
    public void FindOpenReadingFrames_OrfFinderOptions_IgnoreNestedStrandAndMinimumLength()
    {
        const string seq = "TAATTAGCTAAAGCTGTCCACCATGCTGTCCTTCCAGGTAGCTAATTA";
        string Spans(FindOpenReadingFramesResult r) => string.Join(";", r.Items.Select(i => $"{(i.IsReverseComplement ? '-' : '+')},{i.Start},{i.End}"));
        var all = AnalysisTools.FindOpenReadingFrames(seq, 30, startCodons: "AnySenseCodon", nested: false, includeOpenEnded: true,
            minLengthExcludesStop: true, orfFinderMinLength: true);
        Assert.That(Spans(all), Is.EqualTo("+,3,44;+,7,45;-,3,44;-,2,40;-,7,45"));
        var dropped = AnalysisTools.FindOpenReadingFrames(seq, 30, startCodons: "AnySenseCodon", nested: false, includeOpenEnded: true,
            minLengthExcludesStop: true, orfFinderMinLength: true, ignoreNested: true);
        Assert.That(Spans(dropped), Is.EqualTo("-,3,44;-,2,40;-,7,45"));
        // oracle flags: -,3,44 complete; -,2,40 open at 3' (2 nt left); -,7,45 begins at offset 2 of the minus strand (5'-partial)
        Assert.That(dropped.Items.Select(i => i.IsStartPartial).ToArray(), Is.EqualTo(new[] { false, false, true }));
        Assert.That(dropped.Items.Select(i => i.IsTerminated).ToArray(), Is.EqualTo(new[] { true, false, true }));

        var minus = AnalysisTools.FindOpenReadingFrames(seq, 30, startCodons: "AnySenseCodon", nested: false, includeOpenEnded: true,
            minLengthExcludesStop: true, orfFinderMinLength: true, strand: "MINUS");
        Assert.That(Spans(minus), Is.EqualTo("-,3,44;-,2,40;-,7,45"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames(seq, 30, bothStrands: false, strand: "plus"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames(seq, 30, strand: "left"));

        // -ml: values below 30 are changed to 30 (also 0 and negative), whole codons: 30 coding nt pass -ml 31.
        string orf = "ATG" + string.Concat(Enumerable.Repeat("GCC", 9)) + "TAA";
        Assert.That(AnalysisTools.FindOpenReadingFrames(orf, 31, minLengthExcludesStop: true, orfFinderMinLength: true, bothStrands: false).Items, Has.Length.EqualTo(1));
        Assert.That(AnalysisTools.FindOpenReadingFrames(orf, 33, minLengthExcludesStop: true, orfFinderMinLength: true, bothStrands: false).Items, Is.Empty);
        Assert.That(AnalysisTools.FindOpenReadingFrames(orf, -5, minLengthExcludesStop: true, orfFinderMinLength: true, bothStrands: false).Items, Has.Length.EqualTo(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindOpenReadingFrames(orf, 0));
    }

    [Test]
    public void FindOpenReadingFrames_OrfFinderOptions_RegionAndDualCodingTables()
    {
        string orf = "ATG" + string.Concat(Enumerable.Repeat("GCC", 12)) + "TAA";
        string seq = new string('C', 20) + orf + new string('C', 20);          // ORF at 20..61

        var exact = AnalysisTools.FindOpenReadingFrames(seq, 30, bothStrands: false, regionStart: 21, regionEnd: 62).Items.Single();
        Assert.Multiple(() =>
        {
            Assert.That((exact.Start, exact.End), Is.EqualTo((20, 61)), "input coordinates");
            Assert.That(exact.Position, Is.EqualTo(0), "relative to the region");
            Assert.That(exact.IsStartPartial, Is.True);
            Assert.That(exact.IsTerminated, Is.True);
        });
        var cut = AnalysisTools.FindOpenReadingFrames(seq, 30, includeOpenEnded: true, bothStrands: false, regionEnd: 50).Items.Single();
        Assert.That((cut.Start, cut.End, cut.IsTerminated), Is.EqualTo((20, 49, false)));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindOpenReadingFrames(seq, 30, regionStart: 50, regionEnd: 40));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindOpenReadingFrames(seq, 30, regionStart: -5, regionEnd: 100));

        // -g 27: rejected by default, accepted with allowDualCodingStopTables (TGA ends the ORF).
        string tga = "ATG" + string.Concat(Enumerable.Repeat("GCC", 12)) + "TGA";
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindOpenReadingFrames(tga, 30, genetic_code: 27));
        var dual = AnalysisTools.FindOpenReadingFrames(tga, 30, genetic_code: 27, bothStrands: false, allowDualCodingStopTables: true).Items.Single();
        Assert.That((dual.Start, dual.End, dual.IsTerminated), Is.EqualTo((0, 41, true)));
    }
}
