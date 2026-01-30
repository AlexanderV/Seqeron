using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for the BED format parser.
/// </summary>
[TestFixture]
public class BedParserTests
{
    #region Test Data

    private const string SimpleBed3 = @"chr1	100	200
chr1	300	400
chr2	500	600";

    private const string SimpleBed6 = @"chr1	100	200	feature1	500	+
chr1	300	400	feature2	800	-
chr2	500	600	feature3	300	.";

    private const string SimpleBed12 = @"chr1	1000	5000	gene1	900	+	1100	4900	255,0,0	3	100,200,300	0,1000,3700";

    private const string BedWithHeaders = @"track name=test
browser position chr1:1-1000
chr1	100	200	feature1
chr1	300	400	feature2";

    #endregion

    #region Basic Parsing Tests

    [Test]
    public void Parse_BED3_ReturnsCorrectRecords()
    {
        var records = BedParser.Parse(SimpleBed3).ToList();

        Assert.That(records, Has.Count.EqualTo(3));
        Assert.That(records[0].Chrom, Is.EqualTo("chr1"));
        Assert.That(records[0].ChromStart, Is.EqualTo(100));
        Assert.That(records[0].ChromEnd, Is.EqualTo(200));
    }

    [Test]
    public void Parse_BED6_ReturnsCorrectRecords()
    {
        var records = BedParser.Parse(SimpleBed6).ToList();

        Assert.That(records, Has.Count.EqualTo(3));
        Assert.That(records[0].Name, Is.EqualTo("feature1"));
        Assert.That(records[0].Score, Is.EqualTo(500));
        Assert.That(records[0].Strand, Is.EqualTo('+'));
    }

    [Test]
    public void Parse_BED12_ReturnsCorrectRecords()
    {
        var records = BedParser.Parse(SimpleBed12).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].BlockCount, Is.EqualTo(3));
        Assert.That(records[0].BlockSizes, Has.Length.EqualTo(3));
        Assert.That(records[0].BlockStarts, Has.Length.EqualTo(3));
    }

    [Test]
    public void Parse_EmptyContent_ReturnsEmpty()
    {
        var records = BedParser.Parse("").ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_NullContent_ReturnsEmpty()
    {
        var records = BedParser.Parse((string)null!).ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_SkipsTrackAndBrowserLines()
    {
        var records = BedParser.Parse(BedWithHeaders).ToList();
        Assert.That(records, Has.Count.EqualTo(2));
    }

    [Test]
    public void Parse_SkipsComments()
    {
        const string bed = @"# This is a comment
chr1	100	200";

        var records = BedParser.Parse(bed).ToList();
        Assert.That(records, Has.Count.EqualTo(1));
    }

    #endregion

    #region Record Property Tests

    [Test]
    public void BedRecord_Length_CalculatedCorrectly()
    {
        var records = BedParser.Parse(SimpleBed3).ToList();

        Assert.That(records[0].Length, Is.EqualTo(100)); // 200 - 100
        Assert.That(records[2].Length, Is.EqualTo(100)); // 600 - 500
    }

    [Test]
    public void BedRecord_HasBlocks_TrueForBED12()
    {
        var records = BedParser.Parse(SimpleBed12).ToList();
        Assert.That(records[0].HasBlocks, Is.True);
    }

    [Test]
    public void BedRecord_HasBlocks_FalseForBED6()
    {
        var records = BedParser.Parse(SimpleBed6).ToList();
        Assert.That(records[0].HasBlocks, Is.False);
    }

    [Test]
    public void Parse_Score_ClampedTo1000()
    {
        const string bed = "chr1\t100\t200\tname\t1500\t+";
        var records = BedParser.Parse(bed).ToList();

        Assert.That(records[0].Score, Is.EqualTo(1000));
    }

    #endregion

    #region Filtering Tests

    [Test]
    public void FilterByChrom_ReturnsMatchingChromosome()
    {
        var records = BedParser.Parse(SimpleBed3).ToList();
        var chr1 = BedParser.FilterByChrom(records, "chr1").ToList();

        Assert.That(chr1, Has.Count.EqualTo(2));
        Assert.That(chr1.All(r => r.Chrom == "chr1"), Is.True);
    }

    [Test]
    public void FilterByChrom_CaseInsensitive()
    {
        var records = BedParser.Parse(SimpleBed3).ToList();
        var chr1 = BedParser.FilterByChrom(records, "CHR1").ToList();

        Assert.That(chr1, Has.Count.EqualTo(2));
    }

    [Test]
    public void FilterByRegion_ReturnsOverlapping()
    {
        var records = BedParser.Parse(SimpleBed3).ToList();
        var inRegion = BedParser.FilterByRegion(records, "chr1", 150, 350).ToList();

        Assert.That(inRegion, Has.Count.EqualTo(2)); // Both chr1 records overlap
    }

    [Test]
    public void FilterByStrand_ReturnsMatchingStrand()
    {
        var records = BedParser.Parse(SimpleBed6).ToList();
        var plusStrand = BedParser.FilterByStrand(records, '+').ToList();

        Assert.That(plusStrand, Has.Count.EqualTo(1));
        Assert.That(plusStrand[0].Name, Is.EqualTo("feature1"));
    }

    [Test]
    public void FilterByLength_FiltersCorrectly()
    {
        const string bed = @"chr1	100	150	short
chr1	100	300	long";

        var records = BedParser.Parse(bed).ToList();
        var filtered = BedParser.FilterByLength(records, minLength: 100).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
        Assert.That(filtered[0].Name, Is.EqualTo("long"));
    }

    [Test]
    public void FilterByScore_FiltersCorrectly()
    {
        var records = BedParser.Parse(SimpleBed6).ToList();
        var highScore = BedParser.FilterByScore(records, minScore: 400).ToList();

        Assert.That(highScore, Has.Count.EqualTo(2));
    }

    #endregion

    #region Interval Operations Tests

    [Test]
    public void ToIntervals_ConvertsCorrectly()
    {
        var records = BedParser.Parse(SimpleBed3).ToList();
        var intervals = BedParser.ToIntervals(records).ToList();

        Assert.That(intervals, Has.Count.EqualTo(3));
        Assert.That(intervals[0].Start, Is.EqualTo(100));
        Assert.That(intervals[0].End, Is.EqualTo(200));
    }

    [Test]
    public void GenomicInterval_Overlaps_DetectsOverlap()
    {
        var i1 = new BedParser.GenomicInterval("chr1", 100, 200);
        var i2 = new BedParser.GenomicInterval("chr1", 150, 250);
        var i3 = new BedParser.GenomicInterval("chr1", 300, 400);

        Assert.That(i1.Overlaps(i2), Is.True);
        Assert.That(i1.Overlaps(i3), Is.False);
    }

    [Test]
    public void GenomicInterval_Overlaps_DifferentChroms_False()
    {
        var i1 = new BedParser.GenomicInterval("chr1", 100, 200);
        var i2 = new BedParser.GenomicInterval("chr2", 100, 200);

        Assert.That(i1.Overlaps(i2), Is.False);
    }

    [Test]
    public void GenomicInterval_Intersect_ReturnsIntersection()
    {
        var i1 = new BedParser.GenomicInterval("chr1", 100, 200);
        var i2 = new BedParser.GenomicInterval("chr1", 150, 250);

        var intersection = i1.Intersect(i2);

        Assert.That(intersection, Is.Not.Null);
        Assert.That(intersection!.Value.Start, Is.EqualTo(150));
        Assert.That(intersection.Value.End, Is.EqualTo(200));
    }

    [Test]
    public void MergeOverlapping_MergesCorrectly()
    {
        const string bed = @"chr1	100	200
chr1	150	250
chr1	400	500";

        var records = BedParser.Parse(bed).ToList();
        var merged = BedParser.MergeOverlapping(records).ToList();

        Assert.That(merged, Has.Count.EqualTo(2));
        Assert.That(merged[0].ChromStart, Is.EqualTo(100));
        Assert.That(merged[0].ChromEnd, Is.EqualTo(250)); // Merged
    }

    [Test]
    public void Intersect_ReturnsIntersections()
    {
        const string bedA = @"chr1	100	300";
        const string bedB = @"chr1	200	400";

        var a = BedParser.Parse(bedA).ToList();
        var b = BedParser.Parse(bedB).ToList();

        var intersections = BedParser.Intersect(a, b).ToList();

        Assert.That(intersections, Has.Count.EqualTo(1));
        Assert.That(intersections[0].ChromStart, Is.EqualTo(200));
        Assert.That(intersections[0].ChromEnd, Is.EqualTo(300));
    }

    [Test]
    public void Subtract_RemovesOverlapping()
    {
        const string bedA = @"chr1	100	300";
        const string bedB = @"chr1	150	200";

        var a = BedParser.Parse(bedA).ToList();
        var b = BedParser.Parse(bedB).ToList();

        var result = BedParser.Subtract(a, b).ToList();

        Assert.That(result, Has.Count.EqualTo(2));
        Assert.That(result[0].ChromEnd, Is.EqualTo(150));
        Assert.That(result[1].ChromStart, Is.EqualTo(200));
    }

    [Test]
    public void ExpandIntervals_ExpandsCorrectly()
    {
        const string bed = "chr1\t100\t200\tname\t0\t+";

        var records = BedParser.Parse(bed).ToList();
        var expanded = BedParser.ExpandIntervals(records, upstream: 50, downstream: 100).ToList();

        Assert.That(expanded[0].ChromStart, Is.EqualTo(50)); // 100 - 50
        Assert.That(expanded[0].ChromEnd, Is.EqualTo(300)); // 200 + 100
    }

    [Test]
    public void ExpandIntervals_MinusStrand_ExpandsCorrectly()
    {
        const string bed = "chr1\t100\t200\tname\t0\t-";

        var records = BedParser.Parse(bed).ToList();
        var expanded = BedParser.ExpandIntervals(records, upstream: 50, downstream: 100).ToList();

        // For minus strand, upstream is at the end
        Assert.That(expanded[0].ChromStart, Is.EqualTo(0)); // 100 - 100, clamped to 0
        Assert.That(expanded[0].ChromEnd, Is.EqualTo(250)); // 200 + 50
    }

    #endregion

    #region Block Operations Tests

    [Test]
    public void ExpandBlocks_ExpandsBED12()
    {
        var records = BedParser.Parse(SimpleBed12).ToList();
        var expanded = BedParser.ExpandBlocks(records[0]).ToList();

        Assert.That(expanded, Has.Count.EqualTo(3));
    }

    [Test]
    public void GetTotalBlockLength_CalculatesCorrectly()
    {
        var records = BedParser.Parse(SimpleBed12).ToList();
        var totalLength = BedParser.GetTotalBlockLength(records[0]);

        Assert.That(totalLength, Is.EqualTo(600)); // 100 + 200 + 300
    }

    [Test]
    public void GetIntrons_ReturnsIntronRegions()
    {
        var records = BedParser.Parse(SimpleBed12).ToList();
        var introns = BedParser.GetIntrons(records[0]).ToList();

        Assert.That(introns, Has.Count.EqualTo(2)); // 3 blocks = 2 introns
    }

    #endregion

    #region Statistics Tests

    [Test]
    public void CalculateStatistics_ReturnsCorrectStats()
    {
        var records = BedParser.Parse(SimpleBed3).ToList();
        var stats = BedParser.CalculateStatistics(records);

        Assert.That(stats.RecordCount, Is.EqualTo(3));
        Assert.That(stats.TotalBases, Is.EqualTo(300)); // 3 x 100
        Assert.That(stats.MinLength, Is.EqualTo(100));
        Assert.That(stats.MaxLength, Is.EqualTo(100));
    }

    [Test]
    public void CalculateStatistics_ChromosomeCounts_Correct()
    {
        var records = BedParser.Parse(SimpleBed3).ToList();
        var stats = BedParser.CalculateStatistics(records);

        Assert.That(stats.ChromosomeCounts["chr1"], Is.EqualTo(2));
        Assert.That(stats.ChromosomeCounts["chr2"], Is.EqualTo(1));
    }

    [Test]
    public void CalculateStatistics_Empty_ReturnsZeros()
    {
        var stats = BedParser.CalculateStatistics(Array.Empty<BedParser.BedRecord>());

        Assert.That(stats.RecordCount, Is.EqualTo(0));
        Assert.That(stats.TotalBases, Is.EqualTo(0));
    }

    #endregion

    #region Writing Tests

    [Test]
    public void WriteToStream_BED6_ValidOutput()
    {
        var records = BedParser.Parse(SimpleBed6).ToList();
        using var writer = new StringWriter();

        BedParser.WriteToStream(writer, records, BedParser.BedFormat.BED6);
        var output = writer.ToString();

        Assert.That(output, Does.Contain("chr1"));
        Assert.That(output, Does.Contain("feature1"));
    }

    [Test]
    public void WriteAndRead_Roundtrip_PreservesData()
    {
        var original = BedParser.Parse(SimpleBed6).ToList();
        using var writer = new StringWriter();

        BedParser.WriteToStream(writer, original, BedParser.BedFormat.BED6);
        var output = writer.ToString();

        var parsed = BedParser.Parse(output).ToList();

        Assert.That(parsed.Count, Is.EqualTo(original.Count));
        Assert.That(parsed[0].Chrom, Is.EqualTo(original[0].Chrom));
        Assert.That(parsed[0].ChromStart, Is.EqualTo(original[0].ChromStart));
    }

    #endregion

    #region Utility Tests

    [Test]
    public void Sort_SortsByPosition()
    {
        const string bed = @"chr2	100	200
chr1	300	400
chr1	100	200";

        var records = BedParser.Parse(bed).ToList();
        var sorted = BedParser.Sort(records).ToList();

        Assert.That(sorted[0].Chrom, Is.EqualTo("chr1"));
        Assert.That(sorted[0].ChromStart, Is.EqualTo(100));
        Assert.That(sorted[1].Chrom, Is.EqualTo("chr1"));
        Assert.That(sorted[1].ChromStart, Is.EqualTo(300));
        Assert.That(sorted[2].Chrom, Is.EqualTo("chr2"));
    }

    [Test]
    public void CalculateCoverage_ReturnsDepthPerPosition()
    {
        const string bed = @"chr1	100	200
chr1	150	250";

        var records = BedParser.Parse(bed).ToList();
        var coverage = BedParser.CalculateCoverage(records, "chr1", 100, 250).ToList();

        Assert.That(coverage.Count, Is.GreaterThan(0));
    }

    [Test]
    public void ExtractSequence_ReturnsCorrectSequence()
    {
        var record = new BedParser.BedRecord("chr1", 4, 8);
        var reference = "ACGTACGTACGT";

        var sequence = BedParser.ExtractSequence(record, reference);

        Assert.That(sequence, Is.EqualTo("ACGT"));
    }

    [Test]
    public void ExtractSequence_MinusStrand_ReturnsReverseComplement()
    {
        var record = new BedParser.BedRecord("chr1", 0, 4, null, null, '-');
        var reference = "ACGT";

        var sequence = BedParser.ExtractSequence(record, reference);

        Assert.That(sequence, Is.EqualTo("ACGT")); // RC of ACGT is ACGT
    }

    #endregion

    #region Edge Cases

    [Test]
    public void Parse_SpaceSeparated_ParsesCorrectly()
    {
        const string bed = "chr1 100 200";
        var records = BedParser.Parse(bed).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
    }

    [Test]
    public void Parse_InvalidCoordinates_Skips()
    {
        const string bed = @"chr1	abc	200
chr1	100	200";

        var records = BedParser.Parse(bed).ToList();
        Assert.That(records, Has.Count.EqualTo(1));
    }

    [Test]
    public void Parse_TooFewColumns_Skips()
    {
        const string bed = @"chr1	100
chr1	100	200";

        var records = BedParser.Parse(bed).ToList();
        Assert.That(records, Has.Count.EqualTo(1));
    }

    #endregion

    #region File I/O Tests

    [Test]
    public void ParseFile_NonexistentFile_ReturnsEmpty()
    {
        var records = BedParser.ParseFile("nonexistent.bed").ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void ParseFile_ValidFile_ParsesRecords()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            File.WriteAllText(tempFile, SimpleBed6);
            var records = BedParser.ParseFile(tempFile).ToList();

            Assert.That(records, Has.Count.EqualTo(3));
        }
        finally
        {
            File.Delete(tempFile);
        }
    }

    #endregion
}
