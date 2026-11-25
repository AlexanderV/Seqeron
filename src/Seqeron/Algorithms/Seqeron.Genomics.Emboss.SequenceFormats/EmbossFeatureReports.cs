// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using System.Text;
using System.Globalization;
using Seqeron.Genomics.Infrastructure;

/// <summary>The reports of the programs that report the feature tables of the feature library (<c>ajReportWrite</c> with a table and a sequence).</summary>
public static class EmbossFeatureReports
{
    /// <summary>The options of a report as <c>acdSetReport</c> sets them.</summary>
    public static EmbossReportOptions Options(string program, string defaultFormat, EmbossReportQualifiers? rq, IReadOnlyList<string>? commandLine, DateTime? runDate, string reportFile, EmbossFeatureReadContext context, bool footer)
    {
        rq ??= new EmbossReportQualifiers();
        string formatText = rq.Format ?? defaultFormat;
        if (!EmbossReportFormats.TryParse(formatText, out var format))
            throw new EmbossReportException($"Unknown report format '{formatText}'");
        return new EmbossReportOptions
        {
            Program = program,
            CommandLine = commandLine ?? [],
            RunDate = runDate ?? DateTime.Now,
            ReportFile = reportFile,
            Format = format,
            FormatText = formatText,
            Tags = EmbossReportTagLists.TagList(program),
            ShowScore = rq.ShowScore ?? true,
            ShowStrand = rq.ShowStrand ?? true,
            ShowAccession = rq.ShowAccession,
            ShowDescription = rq.ShowDescription,
            ShowUsa = rq.ShowUsa,
            MaxAll = rq.MaxAll,
            MaxSeq = rq.MaxSeq,
            Precision = 3,
            WriteFooter = footer,
            WriteEmptyTables = false,
            Context = context,
        };
    }

    /// <summary>The ACD validation of the report (<c>acdSetReport</c>): a format that needs more tags than the program has (diffseq) is an error and the program dies.</summary>
    public static bool ValidateTags(EmbossReportOptions options, EmbossFeatureReadContext context)
    {
        var flags = EmbossReportFormats.Flags(options.Format);
        if (flags.MinTags <= options.Tags.Count)
            return false;
        context.Error($"Report format '{options.FormatText ?? EmbossReportFormats.Name(options.Format)}' needs {flags.MinTags} tags, has only {options.Tags.Count}");
        return true;
    }

    /// <summary>The report table of a feature table: the full table for the feature formats, the flat list of its top-level features for the others.</summary>
    public static EmbossReport Build(EmbossSeq seq, EmbossSeqFeatureTable table, IReadOnlyList<EmbossReportTag> tags)
    {
        var features = table.Features.Select(f => new EmbossFeature(f.Type, f.Source, f.Start, f.End, f.Score, f.Strand, f.Frame, [.. f.Tags.Select(t => new EmbossFeatureTag(t.Name, t.Value))],
            f.Flags.HasFlag(EmbossLocationBits.BetweenSeq), f.Protein) { FromTable = true }).ToList();
        var sequence = new EmbossReportSequence(seq.Name, MakeUsa(seq), seq.Acc, seq.Desc, seq.Seq.ToString(), EmbossSeqRange.GetBegin(seq), EmbossSeqRange.GetEnd(seq), seq.Offset, seq.Reversed, seq.Type == "N");
        bool protein = table.Type == "P" || (table.Type.Length == 0 && seq.Type == "P");            // ajReportSetType: a table with no type takes the type of the sequence
        return new EmbossReport(sequence, features, protein, table.Type, unchecked((int)table.Start), unchecked((int)table.End), unchecked((int)table.Offset), table.Length, (int)table.Groups, table.FormatName ?? "<null>", tags)
        {
            FeatureTable = table,
            TypeUnknown = table.Type.Length == 0,
        };
    }

    /// <summary>Writes the report file; a format the table is not written in is an error of the library ("Report format ... not supported ...") and the file stays empty.</summary>
    public static string Write(EmbossReportOptions options, IReadOnlyList<EmbossReport> tables, EmbossFeatureReadContext context, Func<int, (long Sequences, long Length)>? totals, out bool failed)
    {
        failed = false;
        try
        {
            return EmbossReportWriter.WriteFile(options, tables, totals);
        }
        catch (EmbossReportException e)
        {
            // ajReportWrite refuses every table of the run: one error for each, and no hit is reported
            for (int i = 0; i < Math.Max(1, tables.Count); i++)
                context.Error(e.Message);
            failed = true;
            return "";
        }
    }

    // seqMakeUsa: the USA of a sequence (format::file:entry and the range)
    public static string MakeUsa(EmbossSeq t)
    {
        var sb = new System.Text.StringBuilder();
        if (t.Db.Length > 0)
            sb.Append(t.Db).Append("-id:").Append(t.Entryname);
        else if (t.Entryname.Length > 0)
            sb.Append(t.Formatstr).Append("::").Append(t.Filename).Append(':').Append(t.Entryname);
        else
            sb.Append(t.Formatstr).Append("::").Append(t.Filename);
        var range = new System.Text.StringBuilder("[");
        if (t.Reversed)
        {
            if (t.End != 0)
                range.Append(-t.End);
            range.Append(':');
            if (t.Begin != 0)
                range.Append(-t.Begin);
            range.Append(":r");
        }
        else
        {
            if (t.Begin != 0)
                range.Append(t.Begin);
            range.Append(':');
            if (t.End != 0)
                range.Append(t.End);
        }

        range.Append(']');
        if (range.Length > 3)
            sb.Append(range);
        return sb.ToString();
    }
}
