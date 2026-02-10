using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace Seqeron.Genomics.Reports;

/// <summary>
/// Generates HTML and JSON reports for genomic analysis results.
/// Supports various analysis types and customizable output formats.
/// </summary>
public static class ReportGenerator
{
    #region Report Types

    /// <summary>Output format for reports</summary>
    public enum OutputFormat
    {
        Html,
        Json,
        Markdown,
        PlainText,
        Csv
    }

    /// <summary>Report section definition</summary>
    public readonly record struct ReportSection(
        string Title,
        string Content,
        ReportSectionType Type,
        IReadOnlyDictionary<string, object>? Data = null);

    /// <summary>Types of report sections</summary>
    public enum ReportSectionType
    {
        Summary,
        Table,
        Chart,
        Sequence,
        Alignment,
        Statistics,
        Custom
    }

    /// <summary>Complete report structure</summary>
    public readonly record struct Report(
        string Title,
        string Description,
        DateTime GeneratedAt,
        IReadOnlyList<ReportSection> Sections,
        IReadOnlyDictionary<string, object> Metadata);

    /// <summary>Table data for reports</summary>
    public readonly record struct TableData(
        IReadOnlyList<string> Headers,
        IReadOnlyList<IReadOnlyList<string>> Rows);

    /// <summary>Chart configuration</summary>
    public readonly record struct ChartConfig(
        string Type,
        string Title,
        IReadOnlyList<string> Labels,
        IReadOnlyList<double> Values,
        IReadOnlyDictionary<string, string>? Options = null);

    #endregion

    #region Report Building

    /// <summary>
    /// Creates a new report builder.
    /// </summary>
    public static ReportBuilder CreateBuilder(string title)
    {
        return new ReportBuilder(title);
    }

    /// <summary>Fluent builder for creating reports</summary>
    public class ReportBuilder
    {
        private readonly string _title;
        private string _description = "";
        private readonly List<ReportSection> _sections = new();
        private readonly Dictionary<string, object> _metadata = new();

        internal ReportBuilder(string title)
        {
            _title = title;
        }

        public ReportBuilder WithDescription(string description)
        {
            _description = description;
            return this;
        }

        public ReportBuilder AddMetadata(string key, object value)
        {
            _metadata[key] = value;
            return this;
        }

        public ReportBuilder AddSection(string title, string content, ReportSectionType type = ReportSectionType.Custom)
        {
            _sections.Add(new ReportSection(title, content, type));
            return this;
        }

        public ReportBuilder AddSummary(string title, string content)
        {
            _sections.Add(new ReportSection(title, content, ReportSectionType.Summary));
            return this;
        }

        public ReportBuilder AddTable(string title, TableData table)
        {
            var data = new Dictionary<string, object> { ["table"] = table };
            _sections.Add(new ReportSection(title, "", ReportSectionType.Table, data));
            return this;
        }

        public ReportBuilder AddTable(string title, IReadOnlyList<string> headers,
            IReadOnlyList<IReadOnlyList<string>> rows)
        {
            return AddTable(title, new TableData(headers, rows));
        }

        public ReportBuilder AddChart(string title, ChartConfig chart)
        {
            var data = new Dictionary<string, object> { ["chart"] = chart };
            _sections.Add(new ReportSection(title, "", ReportSectionType.Chart, data));
            return this;
        }

        public ReportBuilder AddSequence(string title, string sequence, int lineLength = 60)
        {
            var formattedSeq = FormatSequence(sequence, lineLength);
            _sections.Add(new ReportSection(title, formattedSeq, ReportSectionType.Sequence));
            return this;
        }

        public ReportBuilder AddAlignment(string title, IEnumerable<(string Name, string Sequence)> alignedSequences)
        {
            var alignment = FormatAlignment(alignedSequences);
            _sections.Add(new ReportSection(title, alignment, ReportSectionType.Alignment));
            return this;
        }

        public ReportBuilder AddStatistics(string title, IReadOnlyDictionary<string, double> stats)
        {
            var data = new Dictionary<string, object> { ["statistics"] = stats };
            var content = FormatStatistics(stats);
            _sections.Add(new ReportSection(title, content, ReportSectionType.Statistics, data));
            return this;
        }

        public Report Build()
        {
            return new Report(_title, _description, DateTime.UtcNow, _sections, _metadata);
        }
    }

    #endregion

    #region Generation Methods

    /// <summary>
    /// Generates a report in the specified format.
    /// </summary>
    public static string Generate(Report report, OutputFormat format)
    {
        return format switch
        {
            OutputFormat.Html => GenerateHtml(report),
            OutputFormat.Json => GenerateJson(report),
            OutputFormat.Markdown => GenerateMarkdown(report),
            OutputFormat.PlainText => GeneratePlainText(report),
            OutputFormat.Csv => GenerateCsv(report),
            _ => throw new ArgumentException($"Unsupported format: {format}")
        };
    }

    /// <summary>
    /// Saves a report to a file.
    /// </summary>
    public static void SaveToFile(Report report, string filePath, OutputFormat? format = null)
    {
        var outputFormat = format ?? InferFormat(filePath);
        var content = Generate(report, outputFormat);
        File.WriteAllText(filePath, content, Encoding.UTF8);
    }

    private static OutputFormat InferFormat(string filePath)
    {
        var ext = Path.GetExtension(filePath).ToLowerInvariant();
        return ext switch
        {
            ".html" or ".htm" => OutputFormat.Html,
            ".json" => OutputFormat.Json,
            ".md" or ".markdown" => OutputFormat.Markdown,
            ".txt" => OutputFormat.PlainText,
            ".csv" => OutputFormat.Csv,
            _ => OutputFormat.Html
        };
    }

    #endregion

    #region HTML Generation

    private static string GenerateHtml(Report report)
    {
        var sb = new StringBuilder();

        sb.AppendLine("<!DOCTYPE html>");
        sb.AppendLine("<html lang=\"en\">");
        sb.AppendLine("<head>");
        sb.AppendLine("  <meta charset=\"UTF-8\">");
        sb.AppendLine("  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">");
        sb.AppendLine($"  <title>{HtmlEncode(report.Title)}</title>");
        sb.AppendLine("  <style>");
        sb.AppendLine(GetCssStyles());
        sb.AppendLine("  </style>");
        sb.AppendLine("</head>");
        sb.AppendLine("<body>");
        sb.AppendLine("  <div class=\"container\">");

        // Header
        sb.AppendLine("    <header>");
        sb.AppendLine($"      <h1>{HtmlEncode(report.Title)}</h1>");
        if (!string.IsNullOrEmpty(report.Description))
        {
            sb.AppendLine($"      <p class=\"description\">{HtmlEncode(report.Description)}</p>");
        }
        sb.AppendLine($"      <p class=\"timestamp\">Generated: {report.GeneratedAt:yyyy-MM-dd HH:mm:ss} UTC</p>");
        sb.AppendLine("    </header>");

        // Metadata
        if (report.Metadata.Count > 0)
        {
            sb.AppendLine("    <section class=\"metadata\">");
            sb.AppendLine("      <h2>Analysis Information</h2>");
            sb.AppendLine("      <dl>");
            foreach (var (key, value) in report.Metadata)
            {
                sb.AppendLine($"        <dt>{HtmlEncode(key)}</dt>");
                sb.AppendLine($"        <dd>{HtmlEncode(value?.ToString() ?? "")}</dd>");
            }
            sb.AppendLine("      </dl>");
            sb.AppendLine("    </section>");
        }

        // Sections
        foreach (var section in report.Sections)
        {
            sb.AppendLine(GenerateHtmlSection(section));
        }

        sb.AppendLine("  </div>");
        sb.AppendLine("</body>");
        sb.AppendLine("</html>");

        return sb.ToString();
    }

    private static string GenerateHtmlSection(ReportSection section)
    {
        var sb = new StringBuilder();
        var cssClass = section.Type.ToString().ToLowerInvariant();

        sb.AppendLine($"    <section class=\"{cssClass}\">");
        sb.AppendLine($"      <h2>{HtmlEncode(section.Title)}</h2>");

        switch (section.Type)
        {
            case ReportSectionType.Table:
                sb.AppendLine(GenerateHtmlTable(section));
                break;
            case ReportSectionType.Chart:
                sb.AppendLine(GenerateHtmlChart(section));
                break;
            case ReportSectionType.Sequence:
                sb.AppendLine($"      <pre class=\"sequence\">{HtmlEncode(section.Content)}</pre>");
                break;
            case ReportSectionType.Alignment:
                sb.AppendLine($"      <pre class=\"alignment\">{HtmlEncode(section.Content)}</pre>");
                break;
            case ReportSectionType.Statistics:
                sb.AppendLine(GenerateHtmlStatistics(section));
                break;
            default:
                sb.AppendLine($"      <div class=\"content\">{HtmlEncode(section.Content)}</div>");
                break;
        }

        sb.AppendLine("    </section>");
        return sb.ToString();
    }

    private static string GenerateHtmlTable(ReportSection section)
    {
        if (section.Data?.TryGetValue("table", out var tableObj) != true || tableObj is not TableData table)
            return "";

        var sb = new StringBuilder();
        sb.AppendLine("      <table>");
        sb.AppendLine("        <thead>");
        sb.AppendLine("          <tr>");
        foreach (var header in table.Headers)
        {
            sb.AppendLine($"            <th>{HtmlEncode(header)}</th>");
        }
        sb.AppendLine("          </tr>");
        sb.AppendLine("        </thead>");
        sb.AppendLine("        <tbody>");
        foreach (var row in table.Rows)
        {
            sb.AppendLine("          <tr>");
            foreach (var cell in row)
            {
                sb.AppendLine($"            <td>{HtmlEncode(cell)}</td>");
            }
            sb.AppendLine("          </tr>");
        }
        sb.AppendLine("        </tbody>");
        sb.AppendLine("      </table>");
        return sb.ToString();
    }

    private static string GenerateHtmlChart(ReportSection section)
    {
        if (section.Data?.TryGetValue("chart", out var chartObj) != true || chartObj is not ChartConfig chart)
            return "";

        // Generate a simple SVG bar chart
        var sb = new StringBuilder();
        var maxValue = chart.Values.Count > 0 ? chart.Values.Max() : 1;
        var barWidth = 40;
        var gap = 10;
        var height = 200;
        var width = (barWidth + gap) * chart.Values.Count + 50;

        sb.AppendLine($"      <svg class=\"chart\" width=\"{width}\" height=\"{height + 50}\" viewBox=\"0 0 {width} {height + 50}\">");

        for (int i = 0; i < chart.Values.Count; i++)
        {
            var barHeight = maxValue > 0 ? (chart.Values[i] / maxValue) * height : 0;
            var x = 30 + i * (barWidth + gap);
            var y = height - barHeight;

            sb.AppendLine($"        <rect x=\"{x}\" y=\"{y}\" width=\"{barWidth}\" height=\"{barHeight}\" fill=\"#4CAF50\"/>");

            if (i < chart.Labels.Count)
            {
                sb.AppendLine($"        <text x=\"{x + barWidth / 2}\" y=\"{height + 20}\" text-anchor=\"middle\" font-size=\"10\">{HtmlEncode(chart.Labels[i])}</text>");
            }

            sb.AppendLine($"        <text x=\"{x + barWidth / 2}\" y=\"{y - 5}\" text-anchor=\"middle\" font-size=\"10\">{chart.Values[i]:F1}</text>");
        }

        sb.AppendLine("      </svg>");
        return sb.ToString();
    }

    private static string GenerateHtmlStatistics(ReportSection section)
    {
        if (section.Data?.TryGetValue("statistics", out var statsObj) != true ||
            statsObj is not IReadOnlyDictionary<string, double> stats)
            return $"      <div class=\"content\">{HtmlEncode(section.Content)}</div>";

        var sb = new StringBuilder();
        sb.AppendLine("      <dl class=\"statistics\">");
        foreach (var (key, value) in stats)
        {
            sb.AppendLine($"        <dt>{HtmlEncode(key)}</dt>");
            sb.AppendLine($"        <dd>{value:F4}</dd>");
        }
        sb.AppendLine("      </dl>");
        return sb.ToString();
    }

    private static string GetCssStyles()
    {
        return @"
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; line-height: 1.6; color: #333; max-width: 1200px; margin: 0 auto; padding: 20px; }
    .container { background: #fff; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); padding: 30px; }
    header { border-bottom: 2px solid #4CAF50; padding-bottom: 20px; margin-bottom: 30px; }
    h1 { color: #2E7D32; margin: 0 0 10px 0; }
    h2 { color: #388E3C; border-bottom: 1px solid #ddd; padding-bottom: 10px; }
    .description { color: #666; font-size: 1.1em; }
    .timestamp { color: #999; font-size: 0.9em; }
    section { margin-bottom: 30px; }
    table { width: 100%; border-collapse: collapse; margin: 15px 0; }
    th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
    th { background: #4CAF50; color: white; }
    tr:nth-child(even) { background: #f9f9f9; }
    tr:hover { background: #f5f5f5; }
    pre { background: #f4f4f4; padding: 15px; border-radius: 4px; overflow-x: auto; font-family: 'Courier New', monospace; }
    .sequence { font-size: 14px; letter-spacing: 1px; }
    .alignment { font-size: 12px; }
    dl { display: grid; grid-template-columns: 1fr 2fr; gap: 10px; }
    dt { font-weight: bold; color: #555; }
    dd { margin: 0; }
    .chart { display: block; margin: 20px auto; }
    .statistics dt { background: #f0f0f0; padding: 8px; border-radius: 4px; }
    .statistics dd { padding: 8px; }";
    }

    private static string HtmlEncode(string text)
    {
        if (string.IsNullOrEmpty(text)) return "";
        return text
            .Replace("&", "&amp;")
            .Replace("<", "&lt;")
            .Replace(">", "&gt;")
            .Replace("\"", "&quot;")
            .Replace("'", "&#39;");
    }

    #endregion

    #region JSON Generation

    private static string GenerateJson(Report report)
    {
        using var stream = new MemoryStream();
        using (var writer = new Utf8JsonWriter(stream, new JsonWriterOptions { Indented = true }))
        {
            writer.WriteStartObject();
            writer.WriteString("title", report.Title);
            writer.WriteString("description", report.Description);
            writer.WriteString("generatedAt", report.GeneratedAt.ToString("o", CultureInfo.InvariantCulture));

            // Metadata
            if (report.Metadata != null)
            {
                writer.WriteStartObject("metadata");
                foreach (var kvp in report.Metadata)
                {
                    WriteJsonValue(writer, kvp.Key, kvp.Value);
                }
                writer.WriteEndObject();
            }

            // Sections
            writer.WriteStartArray("sections");
            foreach (var s in report.Sections)
            {
                writer.WriteStartObject();
                writer.WriteString("title", s.Title);
                writer.WriteString("content", s.Content);
                writer.WriteString("type", s.Type.ToString());

                if (s.Data != null)
                {
                    writer.WriteStartObject("data");
                    foreach (var kvp in s.Data)
                    {
                        WriteJsonValue(writer, kvp.Key, kvp.Value);
                    }
                    writer.WriteEndObject();
                }

                writer.WriteEndObject();
            }
            writer.WriteEndArray();

            writer.WriteEndObject();
        }

        return Encoding.UTF8.GetString(stream.ToArray());
    }

    private static void WriteJsonValue(Utf8JsonWriter writer, string propertyName, object? value)
    {
        if (value is null) { writer.WriteNull(propertyName); return; }
        switch (value)
        {
            case string s: writer.WriteString(propertyName, s); break;
            case int i: writer.WriteNumber(propertyName, i); break;
            case long l: writer.WriteNumber(propertyName, l); break;
            case double d: writer.WriteNumber(propertyName, d); break;
            case float f: writer.WriteNumber(propertyName, f); break;
            case decimal dec: writer.WriteNumber(propertyName, dec); break;
            case bool b: writer.WriteBoolean(propertyName, b); break;
            case DateTime dt: writer.WriteString(propertyName, dt.ToString("o", CultureInfo.InvariantCulture)); break;
            default: writer.WriteString(propertyName, value.ToString()); break;
        }
    }

    #endregion

    #region Markdown Generation

    private static string GenerateMarkdown(Report report)
    {
        var sb = new StringBuilder();

        sb.AppendLine($"# {report.Title}");
        sb.AppendLine();

        if (!string.IsNullOrEmpty(report.Description))
        {
            sb.AppendLine($"*{report.Description}*");
            sb.AppendLine();
        }

        sb.AppendLine($"**Generated:** {report.GeneratedAt:yyyy-MM-dd HH:mm:ss} UTC");
        sb.AppendLine();

        // Metadata
        if (report.Metadata.Count > 0)
        {
            sb.AppendLine("## Analysis Information");
            sb.AppendLine();
            foreach (var (key, value) in report.Metadata)
            {
                sb.AppendLine($"- **{key}:** {value}");
            }
            sb.AppendLine();
        }

        // Sections
        foreach (var section in report.Sections)
        {
            sb.AppendLine($"## {section.Title}");
            sb.AppendLine();

            switch (section.Type)
            {
                case ReportSectionType.Table:
                    sb.AppendLine(GenerateMarkdownTable(section));
                    break;
                case ReportSectionType.Sequence:
                case ReportSectionType.Alignment:
                    sb.AppendLine("```");
                    sb.AppendLine(section.Content);
                    sb.AppendLine("```");
                    break;
                case ReportSectionType.Statistics:
                    sb.AppendLine(GenerateMarkdownStatistics(section));
                    break;
                default:
                    sb.AppendLine(section.Content);
                    break;
            }
            sb.AppendLine();
        }

        return sb.ToString();
    }

    private static string GenerateMarkdownTable(ReportSection section)
    {
        if (section.Data?.TryGetValue("table", out var tableObj) != true || tableObj is not TableData table)
            return "";

        var sb = new StringBuilder();

        // Headers
        sb.AppendLine("| " + string.Join(" | ", table.Headers) + " |");
        sb.AppendLine("| " + string.Join(" | ", table.Headers.Select(_ => "---")) + " |");

        // Rows
        foreach (var row in table.Rows)
        {
            sb.AppendLine("| " + string.Join(" | ", row) + " |");
        }

        return sb.ToString();
    }

    private static string GenerateMarkdownStatistics(ReportSection section)
    {
        if (section.Data?.TryGetValue("statistics", out var statsObj) != true ||
            statsObj is not IReadOnlyDictionary<string, double> stats)
            return section.Content;

        var sb = new StringBuilder();
        foreach (var (key, value) in stats)
        {
            sb.AppendLine($"- **{key}:** {value:F4}");
        }
        return sb.ToString();
    }

    #endregion

    #region Plain Text Generation

    private static string GeneratePlainText(Report report)
    {
        var sb = new StringBuilder();
        var line = new string('=', 60);
        var subline = new string('-', 40);

        sb.AppendLine(line);
        sb.AppendLine(report.Title.ToUpperInvariant());
        sb.AppendLine(line);
        sb.AppendLine();

        if (!string.IsNullOrEmpty(report.Description))
        {
            sb.AppendLine(report.Description);
            sb.AppendLine();
        }

        sb.AppendLine($"Generated: {report.GeneratedAt:yyyy-MM-dd HH:mm:ss} UTC");
        sb.AppendLine();

        // Metadata
        if (report.Metadata.Count > 0)
        {
            sb.AppendLine("ANALYSIS INFORMATION");
            sb.AppendLine(subline);
            foreach (var (key, value) in report.Metadata)
            {
                sb.AppendLine($"  {key}: {value}");
            }
            sb.AppendLine();
        }

        // Sections
        foreach (var section in report.Sections)
        {
            sb.AppendLine(section.Title.ToUpperInvariant());
            sb.AppendLine(subline);

            switch (section.Type)
            {
                case ReportSectionType.Table:
                    sb.AppendLine(GeneratePlainTextTable(section));
                    break;
                case ReportSectionType.Statistics:
                    sb.AppendLine(GeneratePlainTextStatistics(section));
                    break;
                default:
                    sb.AppendLine(section.Content);
                    break;
            }
            sb.AppendLine();
        }

        return sb.ToString();
    }

    private static string GeneratePlainTextTable(ReportSection section)
    {
        if (section.Data?.TryGetValue("table", out var tableObj) != true || tableObj is not TableData table)
            return "";

        var sb = new StringBuilder();
        var colWidths = new int[table.Headers.Count];

        // Calculate column widths
        for (int i = 0; i < table.Headers.Count; i++)
        {
            colWidths[i] = table.Headers[i].Length;
            foreach (var row in table.Rows)
            {
                if (i < row.Count)
                    colWidths[i] = Math.Max(colWidths[i], row[i].Length);
            }
        }

        // Headers
        for (int i = 0; i < table.Headers.Count; i++)
        {
            sb.Append(table.Headers[i].PadRight(colWidths[i] + 2));
        }
        sb.AppendLine();
        sb.AppendLine(new string('-', colWidths.Sum() + colWidths.Length * 2));

        // Rows
        foreach (var row in table.Rows)
        {
            for (int i = 0; i < row.Count && i < colWidths.Length; i++)
            {
                sb.Append(row[i].PadRight(colWidths[i] + 2));
            }
            sb.AppendLine();
        }

        return sb.ToString();
    }

    private static string GeneratePlainTextStatistics(ReportSection section)
    {
        if (section.Data?.TryGetValue("statistics", out var statsObj) != true ||
            statsObj is not IReadOnlyDictionary<string, double> stats)
            return section.Content;

        var sb = new StringBuilder();
        var maxKeyLen = stats.Keys.Max(k => k.Length);
        foreach (var (key, value) in stats)
        {
            sb.AppendLine($"  {key.PadRight(maxKeyLen)} : {value:F4}");
        }
        return sb.ToString();
    }

    #endregion

    #region CSV Generation

    private static string GenerateCsv(Report report)
    {
        var sb = new StringBuilder();

        // Report header as comments
        sb.AppendLine($"# {report.Title}");
        sb.AppendLine($"# Generated: {report.GeneratedAt:yyyy-MM-dd HH:mm:ss} UTC");
        sb.AppendLine();

        // Find and output tables
        foreach (var section in report.Sections.Where(s => s.Type == ReportSectionType.Table))
        {
            if (section.Data?.TryGetValue("table", out var tableObj) != true || tableObj is not TableData table)
                continue;

            sb.AppendLine($"# {section.Title}");
            sb.AppendLine(string.Join(",", table.Headers.Select(CsvEncode)));
            foreach (var row in table.Rows)
            {
                sb.AppendLine(string.Join(",", row.Select(CsvEncode)));
            }
            sb.AppendLine();
        }

        // Output statistics as table
        foreach (var section in report.Sections.Where(s => s.Type == ReportSectionType.Statistics))
        {
            if (section.Data?.TryGetValue("statistics", out var statsObj) != true ||
                statsObj is not IReadOnlyDictionary<string, double> stats)
                continue;

            sb.AppendLine($"# {section.Title}");
            sb.AppendLine("Metric,Value");
            foreach (var (key, value) in stats)
            {
                sb.AppendLine($"{CsvEncode(key)},{value:F6}");
            }
            sb.AppendLine();
        }

        return sb.ToString();
    }

    private static string CsvEncode(string value)
    {
        if (string.IsNullOrEmpty(value)) return "";
        if (value.Contains(',') || value.Contains('"') || value.Contains('\n'))
        {
            return "\"" + value.Replace("\"", "\"\"") + "\"";
        }
        return value;
    }

    #endregion

    #region Formatting Helpers

    private static string FormatSequence(string sequence, int lineLength)
    {
        if (string.IsNullOrEmpty(sequence))
            return "";

        var sb = new StringBuilder();
        for (int i = 0; i < sequence.Length; i += lineLength)
        {
            var lineNum = (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8);
            var line = sequence.Substring(i, Math.Min(lineLength, sequence.Length - i));

            // Add spaces every 10 characters for readability
            var formatted = new StringBuilder();
            for (int j = 0; j < line.Length; j++)
            {
                if (j > 0 && j % 10 == 0)
                    formatted.Append(' ');
                formatted.Append(line[j]);
            }

            sb.AppendLine($"{lineNum}  {formatted}");
        }
        return sb.ToString();
    }

    private static string FormatAlignment(IEnumerable<(string Name, string Sequence)> alignedSequences)
    {
        var seqs = alignedSequences.ToList();
        if (seqs.Count == 0)
            return "";

        var maxNameLen = seqs.Max(s => s.Name.Length);
        var sb = new StringBuilder();
        var blockSize = 60;

        var maxLen = seqs.Max(s => s.Sequence.Length);

        for (int pos = 0; pos < maxLen; pos += blockSize)
        {
            foreach (var (name, sequence) in seqs)
            {
                var paddedName = name.PadRight(maxNameLen);
                var block = pos < sequence.Length
                    ? sequence.Substring(pos, Math.Min(blockSize, sequence.Length - pos))
                    : "";
                sb.AppendLine($"{paddedName}  {block}");
            }

            // Consensus line
            var consensusLine = new char[Math.Min(blockSize, maxLen - pos)];
            for (int i = 0; i < consensusLine.Length; i++)
            {
                var colPos = pos + i;
                var chars = seqs.Where(s => colPos < s.Sequence.Length).Select(s => s.Sequence[colPos]).ToList();
                consensusLine[i] = chars.Count > 0 && chars.All(c => c == chars[0]) ? '*' : ' ';
            }
            sb.AppendLine($"{new string(' ', maxNameLen)}  {new string(consensusLine)}");
            sb.AppendLine();
        }

        return sb.ToString();
    }

    private static string FormatStatistics(IReadOnlyDictionary<string, double> stats)
    {
        var sb = new StringBuilder();
        foreach (var (key, value) in stats)
        {
            sb.AppendLine($"{key}: {value:F4}");
        }
        return sb.ToString();
    }

    #endregion

    #region Preset Report Templates

    /// <summary>
    /// Simple composition result for report generation.
    /// </summary>
    public readonly record struct CompositionResult(
        double A,
        double C,
        double G,
        double T,
        double GcContent);

    /// <summary>
    /// Creates a sequence analysis report.
    /// </summary>
    public static Report CreateSequenceAnalysisReport(
        string sequenceName,
        string sequence,
        CompositionResult composition,
        IEnumerable<(string Name, int Position)>? motifs = null)
    {
        var builder = CreateBuilder($"Sequence Analysis: {sequenceName}")
            .WithDescription("Comprehensive analysis of sequence composition and features")
            .AddMetadata("Sequence Length", sequence.Length)
            .AddMetadata("Analysis Type", "Composition Analysis");

        // Composition statistics
        builder.AddStatistics("Nucleotide Composition", new Dictionary<string, double>
        {
            ["A Content"] = composition.A,
            ["C Content"] = composition.C,
            ["G Content"] = composition.G,
            ["T Content"] = composition.T,
            ["GC Content"] = composition.GcContent
        });

        // Sequence display
        builder.AddSequence("Sequence", sequence);

        // Motifs table if provided
        if (motifs?.Any() == true)
        {
            var rows = motifs.Select(m => new List<string> { m.Name, m.Position.ToString(CultureInfo.InvariantCulture) } as IReadOnlyList<string>).ToList();
            builder.AddTable("Detected Motifs", new[] { "Motif", "Position" }, rows);
        }

        return builder.Build();
    }

    /// <summary>
    /// Creates a comparison report for multiple sequences.
    /// </summary>
    public static Report CreateComparisonReport(
        string title,
        IEnumerable<(string Name, string Sequence)> sequences,
        double? identityScore = null)
    {
        var seqList = sequences.ToList();
        var builder = CreateBuilder(title)
            .WithDescription("Comparison of multiple sequences")
            .AddMetadata("Number of Sequences", seqList.Count);

        if (identityScore.HasValue)
        {
            builder.AddMetadata("Overall Identity", $"{identityScore.Value:P2}");
        }

        // Summary table
        var rows = seqList.Select(s =>
            new List<string> { s.Name, s.Sequence.Length.ToString(CultureInfo.InvariantCulture) } as IReadOnlyList<string>
        ).ToList();
        builder.AddTable("Sequences", new[] { "Name", "Length" }, rows);

        // Alignment view
        builder.AddAlignment("Sequence Alignment", seqList);

        return builder.Build();
    }

    #endregion
}
