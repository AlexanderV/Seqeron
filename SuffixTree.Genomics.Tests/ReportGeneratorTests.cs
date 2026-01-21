using NUnit.Framework;
using SuffixTree.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

[TestFixture]
public class ReportGeneratorTests
{
    #region Builder Tests

    [Test]
    public void CreateBuilder_ReturnsBuilder()
    {
        var builder = ReportGenerator.CreateBuilder("Test Report");

        Assert.That(builder, Is.Not.Null);
    }

    [Test]
    public void Build_CreatesReportWithTitle()
    {
        var report = ReportGenerator.CreateBuilder("Test Report").Build();

        Assert.That(report.Title, Is.EqualTo("Test Report"));
    }

    [Test]
    public void WithDescription_SetsDescription()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .WithDescription("Test description")
            .Build();

        Assert.That(report.Description, Is.EqualTo("Test description"));
    }

    [Test]
    public void AddMetadata_AddsToReport()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddMetadata("Key1", "Value1")
            .AddMetadata("Key2", 42)
            .Build();

        Assert.That(report.Metadata.Count, Is.EqualTo(2));
        Assert.That(report.Metadata["Key1"], Is.EqualTo("Value1"));
    }

    [Test]
    public void AddSection_AddsSectionToReport()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddSection("Section 1", "Content 1")
            .AddSection("Section 2", "Content 2")
            .Build();

        Assert.That(report.Sections.Count, Is.EqualTo(2));
    }

    [Test]
    public void AddSummary_CreatesSummarySection()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddSummary("Summary", "Summary content")
            .Build();

        Assert.That(report.Sections[0].Type, Is.EqualTo(ReportGenerator.ReportSectionType.Summary));
    }

    [Test]
    public void AddTable_CreatesTableSection()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddTable("Results",
                new[] { "Column1", "Column2" },
                new List<IReadOnlyList<string>> { new[] { "A", "B" } })
            .Build();

        Assert.That(report.Sections[0].Type, Is.EqualTo(ReportGenerator.ReportSectionType.Table));
    }

    [Test]
    public void AddChart_CreatesChartSection()
    {
        var chart = new ReportGenerator.ChartConfig(
            "bar", "Test Chart",
            new[] { "A", "B", "C" },
            new[] { 1.0, 2.0, 3.0 });

        var report = ReportGenerator.CreateBuilder("Test")
            .AddChart("Chart", chart)
            .Build();

        Assert.That(report.Sections[0].Type, Is.EqualTo(ReportGenerator.ReportSectionType.Chart));
    }

    [Test]
    public void AddSequence_CreatesSequenceSection()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddSequence("DNA Sequence", "ACGTACGTACGT")
            .Build();

        Assert.That(report.Sections[0].Type, Is.EqualTo(ReportGenerator.ReportSectionType.Sequence));
    }

    [Test]
    public void AddAlignment_CreatesAlignmentSection()
    {
        var sequences = new[]
        {
            ("Seq1", "ACGT"),
            ("Seq2", "ACGT")
        };

        var report = ReportGenerator.CreateBuilder("Test")
            .AddAlignment("Alignment", sequences)
            .Build();

        Assert.That(report.Sections[0].Type, Is.EqualTo(ReportGenerator.ReportSectionType.Alignment));
    }

    [Test]
    public void AddStatistics_CreatesStatisticsSection()
    {
        var stats = new Dictionary<string, double>
        {
            ["Mean"] = 0.5,
            ["StdDev"] = 0.1
        };

        var report = ReportGenerator.CreateBuilder("Test")
            .AddStatistics("Stats", stats)
            .Build();

        Assert.That(report.Sections[0].Type, Is.EqualTo(ReportGenerator.ReportSectionType.Statistics));
    }

    #endregion

    #region HTML Generation Tests

    [Test]
    public void Generate_Html_ReturnsValidHtml()
    {
        var report = ReportGenerator.CreateBuilder("Test Report")
            .WithDescription("Test description")
            .AddSection("Section 1", "Content")
            .Build();

        var html = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Html);

        Assert.That(html, Does.StartWith("<!DOCTYPE html>"));
        Assert.That(html, Does.Contain("<title>Test Report</title>"));
        Assert.That(html, Does.Contain("Test description"));
        Assert.That(html, Does.Contain("</html>"));
    }

    [Test]
    public void Generate_Html_IncludesTable()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddTable("Data",
                new[] { "Name", "Value" },
                new List<IReadOnlyList<string>> { new[] { "A", "1" }, new[] { "B", "2" } })
            .Build();

        var html = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Html);

        Assert.That(html, Does.Contain("<table>"));
        Assert.That(html, Does.Contain("<th>Name</th>"));
        Assert.That(html, Does.Contain("<td>A</td>"));
    }

    [Test]
    public void Generate_Html_EncodesSpecialCharacters()
    {
        var report = ReportGenerator.CreateBuilder("Test <>&")
            .AddSection("Section", "Content with <script> and & chars")
            .Build();

        var html = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Html);

        Assert.That(html, Does.Contain("&lt;"));
        Assert.That(html, Does.Contain("&gt;"));
        Assert.That(html, Does.Contain("&amp;"));
    }

    [Test]
    public void Generate_Html_IncludesChart()
    {
        var chart = new ReportGenerator.ChartConfig(
            "bar", "Test",
            new[] { "A" },
            new[] { 10.0 });

        var report = ReportGenerator.CreateBuilder("Test")
            .AddChart("Chart", chart)
            .Build();

        var html = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Html);

        Assert.That(html, Does.Contain("<svg"));
        Assert.That(html, Does.Contain("<rect"));
    }

    #endregion

    #region JSON Generation Tests

    [Test]
    public void Generate_Json_ReturnsValidJson()
    {
        var report = ReportGenerator.CreateBuilder("Test Report")
            .AddMetadata("Key", "Value")
            .AddSection("Section", "Content")
            .Build();

        var json = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Json);

        Assert.That(json, Does.Contain("\"title\""));
        Assert.That(json, Does.Contain("Test Report"));
        Assert.That(json, Does.Contain("\"sections\""));
    }

    [Test]
    public void Generate_Json_IncludesMetadata()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddMetadata("Count", 42)
            .Build();

        var json = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Json);

        Assert.That(json, Does.Contain("\"metadata\""));
        Assert.That(json, Does.Contain("Count"));
    }

    #endregion

    #region Markdown Generation Tests

    [Test]
    public void Generate_Markdown_ReturnsValidMarkdown()
    {
        var report = ReportGenerator.CreateBuilder("Test Report")
            .WithDescription("Description")
            .AddSection("Section 1", "Content 1")
            .Build();

        var md = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Markdown);

        Assert.That(md, Does.StartWith("# Test Report"));
        Assert.That(md, Does.Contain("## Section 1"));
    }

    [Test]
    public void Generate_Markdown_IncludesTable()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddTable("Data",
                new[] { "A", "B" },
                new List<IReadOnlyList<string>> { new[] { "1", "2" } })
            .Build();

        var md = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Markdown);

        Assert.That(md, Does.Contain("| A | B |"));
        Assert.That(md, Does.Contain("| --- | --- |"));
        Assert.That(md, Does.Contain("| 1 | 2 |"));
    }

    [Test]
    public void Generate_Markdown_SequenceInCodeBlock()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddSequence("DNA", "ACGT")
            .Build();

        var md = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Markdown);

        Assert.That(md, Does.Contain("```"));
    }

    #endregion

    #region PlainText Generation Tests

    [Test]
    public void Generate_PlainText_ReturnsFormattedText()
    {
        var report = ReportGenerator.CreateBuilder("Test Report")
            .AddSection("Section 1", "Content")
            .Build();

        var text = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.PlainText);

        Assert.That(text, Does.Contain("TEST REPORT"));
        Assert.That(text, Does.Contain("SECTION 1"));
    }

    [Test]
    public void Generate_PlainText_FormatsTable()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddTable("Data",
                new[] { "Name", "Value" },
                new List<IReadOnlyList<string>> { new[] { "Test", "123" } })
            .Build();

        var text = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.PlainText);

        Assert.That(text, Does.Contain("Name"));
        Assert.That(text, Does.Contain("Value"));
        Assert.That(text, Does.Contain("Test"));
    }

    #endregion

    #region CSV Generation Tests

    [Test]
    public void Generate_Csv_ReturnsValidCsv()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddTable("Data",
                new[] { "Col1", "Col2" },
                new List<IReadOnlyList<string>> { new[] { "A", "B" } })
            .Build();

        var csv = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Csv);

        Assert.That(csv, Does.Contain("Col1,Col2"));
        Assert.That(csv, Does.Contain("A,B"));
    }

    [Test]
    public void Generate_Csv_HandlesSpecialCharacters()
    {
        var report = ReportGenerator.CreateBuilder("Test")
            .AddTable("Data",
                new[] { "Name" },
                new List<IReadOnlyList<string>> { new[] { "Value, with comma" } })
            .Build();

        var csv = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Csv);

        Assert.That(csv, Does.Contain("\"Value, with comma\""));
    }

    [Test]
    public void Generate_Csv_IncludesStatistics()
    {
        var stats = new Dictionary<string, double>
        {
            ["Metric1"] = 1.5,
            ["Metric2"] = 2.5
        };

        var report = ReportGenerator.CreateBuilder("Test")
            .AddStatistics("Stats", stats)
            .Build();

        var csv = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Csv);

        Assert.That(csv, Does.Contain("Metric,Value"));
        Assert.That(csv, Does.Contain("Metric1"));
    }

    #endregion

    #region Template Tests

    [Test]
    public void CreateSequenceAnalysisReport_CreatesValidReport()
    {
        var composition = new ReportGenerator.CompositionResult(0.25, 0.25, 0.25, 0.25, 0.5);

        var report = ReportGenerator.CreateSequenceAnalysisReport(
            "TestSeq",
            "ACGTACGT",
            composition);

        Assert.That(report.Title, Does.Contain("TestSeq"));
        Assert.That(report.Sections.Count, Is.GreaterThan(0));
    }

    [Test]
    public void CreateSequenceAnalysisReport_IncludesMotifs()
    {
        var composition = new ReportGenerator.CompositionResult(0.25, 0.25, 0.25, 0.25, 0.5);
        var motifs = new[] { ("TATA", 10), ("GC-box", 25) };

        var report = ReportGenerator.CreateSequenceAnalysisReport(
            "TestSeq",
            "ACGT",
            composition,
            motifs);

        var tableSection = report.Sections.FirstOrDefault(s => s.Type == ReportGenerator.ReportSectionType.Table);
        Assert.That(tableSection.Title, Does.Contain("Motif"));
    }

    [Test]
    public void CreateComparisonReport_CreatesValidReport()
    {
        var sequences = new[]
        {
            ("Seq1", "ACGTACGT"),
            ("Seq2", "ACGTACGT")
        };

        var report = ReportGenerator.CreateComparisonReport("Comparison", sequences, 1.0);

        Assert.That(report.Title, Is.EqualTo("Comparison"));
        Assert.That(report.Metadata.ContainsKey("Number of Sequences"), Is.True);
    }

    #endregion

    #region Edge Cases

    [Test]
    public void Generate_EmptyReport_ProducesOutput()
    {
        var report = ReportGenerator.CreateBuilder("Empty").Build();

        var html = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Html);
        var json = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Json);
        var md = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Markdown);

        Assert.That(html, Is.Not.Empty);
        Assert.That(json, Is.Not.Empty);
        Assert.That(md, Is.Not.Empty);
    }

    [Test]
    public void Generate_LongSequence_FormatsCorrectly()
    {
        var longSeq = new string('A', 200);
        var report = ReportGenerator.CreateBuilder("Test")
            .AddSequence("Long Sequence", longSeq)
            .Build();

        var html = ReportGenerator.Generate(report, ReportGenerator.OutputFormat.Html);

        Assert.That(html, Does.Contain("AAAA"));
    }

    [Test]
    public void Build_SetsGeneratedAtTimestamp()
    {
        var before = System.DateTime.UtcNow;
        var report = ReportGenerator.CreateBuilder("Test").Build();
        var after = System.DateTime.UtcNow;

        Assert.That(report.GeneratedAt, Is.InRange(before, after));
    }

    #endregion
}
