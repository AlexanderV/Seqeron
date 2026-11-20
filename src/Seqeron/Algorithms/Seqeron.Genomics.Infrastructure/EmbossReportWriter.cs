using System.Text;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>The EMBOSS 6.6 feature-report formats (<c>-rformat</c>, <c>ajreport.c</c>) the report layer writes. The remaining names of the program's list (seqtable, table, embl, genbank, ... ) are not written yet (B09.md A5-8).</summary>
public enum EmbossReportFormat
{
    /// <summary>GFF version 3 (<c>-rformat gff</c>).</summary>
    Gff,
    /// <summary>Tab-separated columns with a header line (<c>-rformat excel</c>).</summary>
    Excel,
}

/// <summary>One feature of an EMBOSS report in the form the formats print it (audit A5-8, B09 F176).</summary>
/// <param name="ExcelStart">First coordinate of the <c>excel</c> report (the program's own, not reordered).</param>
/// <param name="ExcelEnd">Second coordinate of the <c>excel</c> report.</param>
/// <param name="GffStart">First coordinate of the <c>gff</c> report.</param>
/// <param name="GffEnd">Second coordinate of the <c>gff</c> report.</param>
/// <param name="Score">The score as printed.</param>
/// <param name="ExcelStrand">The strand column of the <c>excel</c> report (<c>+</c> / <c>-</c>).</param>
/// <param name="GffStrand">The strand column of the <c>gff</c> report (<c>+</c> / <c>-</c> / <c>.</c>).</param>
/// <param name="GffPhase">The phase column of the <c>gff</c> report (<c>.</c>, or the frame minus one).</param>
/// <param name="ExcelColumns">The program's own columns after the strand of the <c>excel</c> report.</param>
/// <param name="GffAttributes">The attributes after <c>ID=name.N;</c> of the <c>gff</c> report, one per <c>;</c>.</param>
public sealed record EmbossReportFeature(
    long ExcelStart, long ExcelEnd, long GffStart, long GffEnd, string Score, string ExcelStrand, string GffStrand, string GffPhase, IReadOnlyList<string> ExcelColumns, IReadOnlyList<string> GffAttributes);

/// <summary>What a program reports about one sequence: the header values and the features (in the order the program writes them).</summary>
/// <param name="Program">The program name (the gff source column).</param>
/// <param name="SequenceName">The sequence name.</param>
/// <param name="RegionEnd">The end of <c>##sequence-region</c> (the sequence length for most programs; fuzztran lets the features decide, see <see cref="EmbossReportWriter.UnsignedRegionEnd"/>).</param>
/// <param name="GffType">The SO feature type (<c>nucleotide_motif</c>, <c>polypeptide_motif</c>, <c>sequence_feature</c>, <c>repeat_region</c>).</param>
/// <param name="Protein"><c>#!Type Protein</c> instead of <c>DNA</c>.</param>
/// <param name="ExcelHeader">The excel header columns after <c>Strand</c>.</param>
/// <param name="Features">The features.</param>
/// <param name="WriteWhenEmpty">The program writes the header of the report even when it found nothing (etandem, equicktandem); fuzznuc, fuzzpro and fuzztran write nothing.</param>
public sealed record EmbossReport(string Program, string SequenceName, long RegionEnd, string GffType, bool Protein, IReadOnlyList<string> ExcelHeader, IReadOnlyList<EmbossReportFeature> Features, bool WriteWhenEmpty = false);

/// <summary>
/// The one writer of the EMBOSS feature-report texts (<c>-rformat gff</c> and <c>excel</c> so far), shared by fuzznuc, fuzzpro, fuzztran, etandem and equicktandem (audit A5-8, B09 F176): the layout of
/// <c>ajreport.c</c> / <c>ajfeat.c</c> as measured on the installed 6.6.0. A report without a feature is empty (the program prints only a warning).
/// </summary>
public static class EmbossReportWriter
{
    /// <summary>The end of <c>##sequence-region</c> of a report whose feature table has no length (fuzztran): the largest coordinate read as an unsigned 32-bit integer, so a report with negative reverse-frame coordinates names a number near 2^32.</summary>
    public static long UnsignedRegionEnd(IEnumerable<EmbossReportFeature> features)
    {
        long end = 0;
        foreach (var f in features)
            end = Math.Max(end, Math.Max(unchecked((uint)f.GffStart), unchecked((uint)f.GffEnd)));
        return end;
    }

    // GFF3 attribute escaping (ajfeat.c): the separators of the attribute list inside a value are percent-encoded; the first '=' of an attribute (key = value) is the separator itself
    private static string Escape(string attribute)
    {
        int eq = attribute.IndexOf('=');
        string key = eq < 0 ? attribute : attribute[..(eq + 1)];
        string value = eq < 0 ? "" : attribute[(eq + 1)..];
        var sb = new StringBuilder(key);
        foreach (char c in value)
        {
            switch (c)
            {
                case '%': sb.Append("%25"); break;
                case ',': sb.Append("%2C"); break;
                case ';': sb.Append("%3B"); break;
                case '=': sb.Append("%3D"); break;
                case '&': sb.Append("%26"); break;
                case '\t': sb.Append("%09"); break;
                case '\n': sb.Append("%0A"); break;
                case '\r': sb.Append("%0D"); break;
                default: sb.Append(c); break;
            }
        }

        return sb.ToString();
    }

    /// <summary>Writes the report in <paramref name="format"/>; <paramref name="date"/> is the <c>#!Date</c> of the gff header (the date of the run).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="report"/> is null.</exception>
    public static string Write(EmbossReport report, EmbossReportFormat format, DateTime date)
    {
        ArgumentNullException.ThrowIfNull(report);
        if (report.Features.Count == 0 && !report.WriteWhenEmpty)
            return string.Empty;
        var sb = new StringBuilder();
        if (format == EmbossReportFormat.Excel)
        {
            sb.Append("SeqName\tStart\tEnd\tScore\tStrand");
            foreach (string c in report.ExcelHeader)
                sb.Append('\t').Append(c);
            sb.Append('\n');
            foreach (var f in report.Features)
            {
                sb.Append(report.SequenceName).Append('\t').Append(f.ExcelStart).Append('\t').Append(f.ExcelEnd).Append('\t').Append(f.Score).Append('\t').Append(f.ExcelStrand);
                foreach (string c in f.ExcelColumns)
                    sb.Append('\t').Append(c);
                sb.Append('\n');
            }

            return sb.ToString();
        }

        sb.Append("##gff-version 3\n").Append("##sequence-region ").Append(report.SequenceName).Append(" 1 ").Append(report.RegionEnd).Append('\n');
        sb.Append("#!Date ").Append(date.ToString("yyyy-MM-dd", System.Globalization.CultureInfo.InvariantCulture)).Append('\n');
        sb.Append("#!Type ").Append(report.Protein ? "Protein" : "DNA").Append('\n').Append("#!Source-version EMBOSS 6.6.0.0\n");
        int n = 0;
        foreach (var f in report.Features)
        {
            n++;
            sb.Append(report.SequenceName).Append('\t').Append(report.Program).Append('\t').Append(report.GffType).Append('\t').Append(f.GffStart).Append('\t').Append(f.GffEnd).Append('\t').Append(f.Score == "0" ? "." : f.Score)       // a zero score is undefined in GFF
                .Append('\t').Append(f.GffStrand).Append('\t').Append(f.GffPhase).Append("\tID=").Append(report.SequenceName).Append('.').Append(n);
            foreach (string a in f.GffAttributes)
                sb.Append(';').Append(Escape(a));
            sb.Append('\n');
        }

        return sb.ToString();
    }
}
