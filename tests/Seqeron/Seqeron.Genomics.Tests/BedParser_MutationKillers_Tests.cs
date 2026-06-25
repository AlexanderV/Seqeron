using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.IO;
using static Seqeron.Genomics.IO.BedParser;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Targeted mutation-killing tests for BedParser.cs (checklist 04 row 66, PARSE-BED-001).
/// The canonical suite left the UCSC BED12 field parsing/validation, the 0-based half-open
/// coordinate semantics, the inclusive filters, and the interval algebra (intersect, subtract,
/// merge, coverage, block expansion) under-pinned. These pin the exact UCSC spec behaviour so
/// the boundary/logical/null-coalescing mutants diverge.
/// </summary>
[TestFixture]
public class BedParser_MutationKillers_Tests
{
    // A valid BED12 line per UCSC: starts[0]==0; starts[n-1]+sizes[n-1]==chromEnd-chromStart; non-overlapping.
    private const string Bed12Line = "chr1\t100\t200\tfeat1\t500\t+\t100\t200\t255,0,0\t2\t20,30\t0,70";

    // ── ParseLine: every BED12 column parsed exactly ──────────────────────────────────

    [Test]
    public void Parse_Bed12_ParsesAllTwelveColumnsExactly()
    {
        var rec = Parse(Bed12Line).Single();

        rec.Chrom.Should().Be("chr1");
        rec.ChromStart.Should().Be(100);
        rec.ChromEnd.Should().Be(200);
        rec.Name.Should().Be("feat1");
        rec.Score.Should().Be(500);
        rec.Strand.Should().Be('+');
        rec.ThickStart.Should().Be(100);
        rec.ThickEnd.Should().Be(200);
        rec.ItemRgb.Should().Be("255,0,0");
        rec.BlockCount.Should().Be(2);
        rec.BlockSizes.Should().Equal(20, 30);
        rec.BlockStarts.Should().Equal(0, 70);
        rec.Length.Should().Be(100);            // ChromEnd - ChromStart
        rec.HasBlocks.Should().BeTrue();
    }

    [Test]
    public void Parse_Bed3_HasNoOptionalFields()
    {
        var rec = Parse("chr2\t50\t75").Single();
        rec.Chrom.Should().Be("chr2");
        rec.ChromStart.Should().Be(50);
        rec.ChromEnd.Should().Be(75);
        rec.Name.Should().BeNull();
        rec.Score.Should().BeNull();
        rec.Strand.Should().BeNull();
        rec.HasBlocks.Should().BeFalse();
    }

    [Test]
    public void Parse_ScoreClampedToUcscRange()
    {
        // BED score domain is 0..1000; 5000 must clamp to 1000.
        Parse("chr1\t0\t10\tn\t5000\t+").Single().Score.Should().Be(1000);
    }

    [Test]
    public void Parse_StartGreaterThanEnd_IsRejected()
    {
        // Per UCSC spec chromStart <= chromEnd; a reversed interval yields no record.
        Parse("chr1\t100\t50").Should().BeEmpty();
    }

    [Test]
    public void Parse_Bed12_FirstBlockStartNotZero_IsRejected()
    {
        // Per UCSC spec the first blockStart must be 0.
        Parse("chr1\t100\t200\tn\t0\t+\t100\t200\t0,0,0\t2\t20,30\t5,70").Should().BeEmpty();
    }

    [Test]
    public void Parse_Bed12_LastBlockDoesNotReachEnd_IsRejected()
    {
        // Per UCSC spec final blockStart + final blockSize must equal chromEnd - chromStart.
        Parse("chr1\t100\t200\tn\t0\t+\t100\t200\t0,0,0\t2\t20,30\t0,60").Should().BeEmpty();
    }

    [Test]
    public void Parse_Bed12_OverlappingBlocks_IsRejected()
    {
        // Blocks may not overlap: starts[1]=10 < starts[0]+sizes[0]=20.
        Parse("chr1\t100\t200\tn\t0\t+\t100\t200\t0,0,0\t2\t20,80\t0,10").Should().BeEmpty();
    }

    // ── Filters: inclusive bounds ──────────────────────────────────────────────────────

    [Test]
    public void FilterByLength_BoundsAreInclusive()
    {
        var recs = new[]
        {
            new BedRecord("c", 0, 5),    // len 5
            new BedRecord("c", 0, 10),   // len 10 (== min, kept)
            new BedRecord("c", 0, 20),   // len 20 (== max, kept)
            new BedRecord("c", 0, 21),   // len 21 (> max, dropped)
        };

        FilterByLength(recs, 10, 20).Select(r => r.Length).Should().Equal(10, 20);
    }

    [Test]
    public void FilterByScore_BoundsAreInclusive()
    {
        var recs = new[]
        {
            new BedRecord("c", 0, 1, Score: 99),
            new BedRecord("c", 0, 1, Score: 100), // == min, kept
            new BedRecord("c", 0, 1, Score: 500), // == max, kept
            new BedRecord("c", 0, 1, Score: 501), // > max, dropped
            new BedRecord("c", 0, 1, Score: null) // no score, dropped
        };

        FilterByScore(recs, 100, 500).Select(r => r.Score!.Value).Should().Equal(100, 500);
    }

    [Test]
    public void FilterByRegion_UsesHalfOpenOverlap()
    {
        var recs = new[]
        {
            new BedRecord("chr1", 0, 20, Name: "a"),   // overlaps query [20,40)? End(20) > start(20)? no ⇒ dropped
            new BedRecord("chr1", 20, 25, Name: "b"),  // ChromStart(20) < end(40) ⇒ kept
            new BedRecord("chr1", 40, 50, Name: "c"),  // ChromStart(40) < end(40)? no ⇒ dropped
        };

        FilterByRegion(recs, "chr1", 20, 40).Select(r => r.Name).Should().Equal("b");
    }

    // ── Intersect: half-open intersection with field coalescing (a ?? b) ───────────────

    [Test]
    public void Intersect_ComputesOverlapAndCoalescesFieldsLeftThenRight()
    {
        var a = new[] { new BedRecord("chr1", 10, 30, Name: "A", Score: null, Strand: '+') };
        var b = new[] { new BedRecord("chr1", 20, 40, Name: "B", Score: 5, Strand: '-') };

        var res = Intersect(a, b).Single();

        res.ChromStart.Should().Be(20);     // max(10,20)
        res.ChromEnd.Should().Be(30);       // min(30,40)
        res.Name.Should().Be("A");          // a.Name ?? b.Name
        res.Score.Should().Be(5);           // a.Score(null) ?? b.Score
        res.Strand.Should().Be('+');        // a.Strand ?? b.Strand
    }

    [Test]
    public void Intersect_CoalescesFromBWhenAFieldsMissing()
    {
        var a = new[] { new BedRecord("chr1", 10, 30, Name: null, Score: 7, Strand: null) };
        var b = new[] { new BedRecord("chr1", 20, 40, Name: "B2", Score: 9, Strand: '-') };

        var res = Intersect(a, b).Single();

        res.Name.Should().Be("B2");         // a.Name(null) ?? b.Name
        res.Score.Should().Be(7);           // a.Score ?? b.Score
        res.Strand.Should().Be('-');        // a.Strand(null) ?? b.Strand
    }

    [Test]
    public void Intersect_NonOverlapping_YieldsNothing()
    {
        var a = new[] { new BedRecord("chr1", 0, 10) };
        var b = new[] { new BedRecord("chr1", 10, 20) }; // touches at 10 but half-open ⇒ no overlap
        Intersect(a, b).Should().BeEmpty();
    }

    [Test]
    public void Intersect_SkipsBTouchingAtStartBoundary()
    {
        var a = new[] { new BedRecord("chr1", 10, 30) };
        // b1 ends exactly at a.ChromStart (10) ⇒ `b.ChromEnd <= a.ChromStart` continue (no overlap);
        // b2 overlaps. Kills `<= a.ChromStart` → `< a.ChromStart`, which would emit a spurious [10,10].
        var b = new[] { new BedRecord("chr1", 5, 10), new BedRecord("chr1", 20, 40) };

        var res = Intersect(a, b).ToList();
        res.Should().ContainSingle();
        res[0].ChromStart.Should().Be(20);
        res[0].ChromEnd.Should().Be(30);
    }

    // ── Subtract: removes B intervals, splitting A as needed ───────────────────────────

    [Test]
    public void Subtract_SplitsIntervalAroundInteriorHole()
    {
        var a = new[] { new BedRecord("chr1", 0, 100) };
        var b = new[] { new BedRecord("chr1", 40, 60) };

        var res = Subtract(a, b).ToList();

        res.Should().HaveCount(2);
        res[0].ChromStart.Should().Be(0);
        res[0].ChromEnd.Should().Be(40);    // left piece [0,40)
        res[1].ChromStart.Should().Be(60);
        res[1].ChromEnd.Should().Be(100);   // right piece [60,100)
    }

    [Test]
    public void Subtract_BIntervalsTouchingEdges_TrimWithoutZeroLengthPieces()
    {
        var a = new[] { new BedRecord("chr1", 10, 50) };
        // b1 starts exactly at a.start (10); b2 ends exactly at a.end (50). The result must be the
        // single interior interval [20,40] — no zero-length [10,10] or [50,50] pieces. This kills
        // `b.ChromStart > start` → `>= start` and `b.ChromEnd < end` → `<= end`.
        var b = new[] { new BedRecord("chr1", 10, 20), new BedRecord("chr1", 40, 50) };

        var res = Subtract(a, b).ToList();
        res.Should().ContainSingle();
        res[0].ChromStart.Should().Be(20);
        res[0].ChromEnd.Should().Be(40);
    }

    [Test]
    public void Subtract_NoOverlappingBOnChrom_YieldsAUnchanged()
    {
        var a = new[] { new BedRecord("chr1", 0, 100, Name: "keep") };
        var b = new[] { new BedRecord("chr2", 0, 100) }; // different chromosome

        Subtract(a, b).Single().Name.Should().Be("keep");
    }

    // ── MergeOverlapping: merges intervals that overlap/touch (half-open contact) ───────

    [Test]
    public void MergeOverlapping_MergesTouchingIntervals()
    {
        var recs = new[]
        {
            new BedRecord("chr1", 0, 10),
            new BedRecord("chr1", 10, 20), // ChromStart(10) <= current.ChromEnd(10) ⇒ merge
        };

        var merged = MergeOverlapping(recs).Single();
        merged.ChromStart.Should().Be(0);
        merged.ChromEnd.Should().Be(20);
    }

    // ── Block operations (BED12) ───────────────────────────────────────────────────────

    [Test]
    public void ExpandBlocks_ProducesAbsoluteExonCoordinates()
    {
        var rec = Parse(Bed12Line).Single();
        var blocks = ExpandBlocks(rec).ToList();

        blocks.Should().HaveCount(2);
        blocks[0].ChromStart.Should().Be(100);  // start + blockStarts[0]
        blocks[0].ChromEnd.Should().Be(120);    // + blockSizes[0]
        blocks[0].Name.Should().Be("feat1_block1");
        blocks[1].ChromStart.Should().Be(170);  // start + blockStarts[1]
        blocks[1].ChromEnd.Should().Be(200);    // + blockSizes[1]
        blocks[1].Name.Should().Be("feat1_block2");
    }

    [Test]
    public void GetTotalBlockLength_SumsBlockSizes()
    {
        GetTotalBlockLength(Parse(Bed12Line).Single()).Should().Be(50); // 20 + 30
    }

    [Test]
    public void GetIntrons_ReturnsGapsBetweenBlocks()
    {
        var rec = Parse(Bed12Line).Single();
        var introns = GetIntrons(rec).ToList();

        introns.Should().HaveCount(1);
        introns[0].ChromStart.Should().Be(120); // end of block 1
        introns[0].ChromEnd.Should().Be(170);   // start of block 2
        introns[0].Name.Should().Be("feat1_intron1");
    }

    // ── CalculateCoverage: per-position depth over a queried region ────────────────────

    [Test]
    public void CalculateCoverage_ProducesStepwiseDepthProfile()
    {
        var recs = new[]
        {
            new BedRecord("chr1", 0, 10),
            new BedRecord("chr1", 5, 15),
        };

        var cov = CalculateCoverage(recs, "chr1", 0, 20).ToList();

        cov.Should().Equal((0, 1), (5, 2), (10, 1), (15, 0));
    }

    // ── ExtractSequence: 0-based half-open; minus strand reverse-complements ───────────

    [Test]
    public void ExtractSequence_PlusStrand_UsesZeroBasedHalfOpen()
    {
        var rec = new BedRecord("chr1", 1, 4);
        ExtractSequence(rec, "ACGTAC").Should().Be("CGT"); // slice [1,4)
    }

    [Test]
    public void ExtractSequence_MinusStrand_ReverseComplements()
    {
        var rec = new BedRecord("chr1", 1, 4, Strand: '-');
        ExtractSequence(rec, "ACGTAC").Should().Be("ACG"); // revcomp of "CGT"
    }

    // ── GenomicInterval.Overlaps: half-open contact does NOT overlap ───────────────────

    [Test]
    public void GenomicInterval_Overlaps_IsHalfOpen()
    {
        var x = new GenomicInterval("chr1", 0, 10);
        x.Overlaps(new GenomicInterval("chr1", 5, 15)).Should().BeTrue();   // genuine overlap
        x.Overlaps(new GenomicInterval("chr1", 10, 20)).Should().BeFalse(); // touch at 10 ⇒ no overlap
        x.Overlaps(new GenomicInterval("chr2", 5, 15)).Should().BeFalse();  // different chromosome
    }

    // ── Writing: column count grows with format level ──────────────────────────────────

    [Test]
    public void WriteToStream_Bed6_EmitsExactlySixTabSeparatedColumns()
    {
        var rec = new BedRecord("chr1", 100, 200, Name: "feat", Score: 500, Strand: '+');
        using var sw = new StringWriter();
        WriteToStream(sw, new[] { rec }, BedFormat.BED6);

        sw.ToString().TrimEnd('\n').Should().Be("chr1\t100\t200\tfeat\t500\t+");
    }

    [Test]
    public void WriteToStream_Bed4_EmitsNameButNotScoreOrStrand()
    {
        var rec = new BedRecord("chr1", 100, 200, Name: "feat", Score: 500, Strand: '+');
        using var sw = new StringWriter();
        WriteToStream(sw, new[] { rec }, BedFormat.BED4);
        // format >= BED4 ⇒ Name column present; format >= BED5/BED6 are false ⇒ no score/strand.
        sw.ToString().TrimEnd('\n').Should().Be("chr1\t100\t200\tfeat");
    }

    [Test]
    public void WriteToStream_Bed5_EmitsNameAndScoreButNotStrand()
    {
        var rec = new BedRecord("chr1", 100, 200, Name: "feat", Score: 500, Strand: '+');
        using var sw = new StringWriter();
        WriteToStream(sw, new[] { rec }, BedFormat.BED5);
        sw.ToString().TrimEnd('\n').Should().Be("chr1\t100\t200\tfeat\t500");
    }

    [Test]
    public void WriteToStream_Bed12_EmitsAllTwelveColumns()
    {
        var rec = Parse(Bed12Line).Single();
        using var sw = new StringWriter();
        WriteToStream(sw, new[] { rec }, BedFormat.BED12);

        sw.ToString().TrimEnd('\n').Should().Be(Bed12Line);
    }
}
