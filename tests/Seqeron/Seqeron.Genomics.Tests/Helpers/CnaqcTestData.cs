using System.Globalization;

namespace Seqeron.Genomics.Tests.Helpers;

/// <summary>
/// Loads the seeded synthetic CNAqc datasets (TestData/CNAqc/cnaqc_D*.tsv: VAF, Major, minor, NV, DP) used to lock
/// the KDE / peakPick ports and <c>AnalyzePurityPeaks</c> against CNAqc 1.1.5 run in R (ONCO-PURITY-001, FIN-B24 F33/F34).
/// </summary>
internal static class CnaqcTestData
{
    internal readonly record struct Row(double Vaf, int Major, int Minor, int Nv, int Dp);

    internal static IReadOnlyList<Row> Load(string dataset)
    {
        string resource = $"Seqeron.Genomics.Tests.TestData.CNAqc.cnaqc_{dataset}.tsv";
        using Stream stream = typeof(CnaqcTestData).Assembly.GetManifestResourceStream(resource)
            ?? throw new InvalidOperationException($"Embedded resource '{resource}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        reader.ReadLine(); // header
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            string[] f = line.Split('\t');
            rows.Add(new Row(
                double.Parse(f[0], CultureInfo.InvariantCulture),
                int.Parse(f[1], CultureInfo.InvariantCulture),
                int.Parse(f[2], CultureInfo.InvariantCulture),
                int.Parse(f[3], CultureInfo.InvariantCulture),
                int.Parse(f[4], CultureInfo.InvariantCulture)));
        }

        return rows;
    }

    /// <summary>Raw text of an embedded TestData/CNAqc file (e.g. the R reference outputs of FIN-B24 F62).</summary>
    internal static string Text(string file)
    {
        string resource = $"Seqeron.Genomics.Tests.TestData.CNAqc.{file}";
        using Stream stream = typeof(CnaqcTestData).Assembly.GetManifestResourceStream(resource)
            ?? throw new InvalidOperationException($"Embedded resource '{resource}' not found.");
        using var reader = new StreamReader(stream);
        return reader.ReadToEnd();
    }

    /// <summary>
    /// Subclonal segment set S1–S3 (cnaqc_S*_segments.tsv: chr, from, to, k1, k2, CCF; cnaqc_S*_mutations.tsv: 1-based
    /// segment, VAF, NV, DP) as <see cref="Seqeron.Genomics.Oncology.OncologyAnalyzer.SubclonalPeakSegment"/>s.
    /// </summary>
    internal static List<Seqeron.Genomics.Oncology.OncologyAnalyzer.SubclonalPeakSegment> Segments(string dataset)
    {
        static List<string[]> Rows(string text) => text.Split('\n').Skip(1).Where(l => l.Length > 0).Select(l => l.TrimEnd('\r').Split('\t')).ToList();
        var mutations = Rows(Text($"cnaqc_{dataset}_mutations.tsv"));
        var segments = Rows(Text($"cnaqc_{dataset}_segments.tsv"));
        return segments.Select((f, i) =>
        {
            int[] k1 = f[3].Split(':').Select(v => int.Parse(v, CultureInfo.InvariantCulture)).ToArray();
            int[] k2 = f[4].Split(':').Select(v => int.Parse(v, CultureInfo.InvariantCulture)).ToArray();
            return new Seqeron.Genomics.Oncology.OncologyAnalyzer.SubclonalPeakSegment(
                f[0], (long)double.Parse(f[1], CultureInfo.InvariantCulture), (long)double.Parse(f[2], CultureInfo.InvariantCulture),
                k1[0], k1[1], k2[0], k2[1], double.Parse(f[5], CultureInfo.InvariantCulture),
                mutations.Where(m => int.Parse(m[0], CultureInfo.InvariantCulture) == i + 1)
                    .Select(m => double.Parse(m[1], CultureInfo.InvariantCulture)).ToList());
        }).ToList();
    }

    internal static double[] Vafs(string dataset, int major, int minor) =>
        Load(dataset).Where(r => r.Major == major && r.Minor == minor).Select(r => r.Vaf).ToArray();
}
