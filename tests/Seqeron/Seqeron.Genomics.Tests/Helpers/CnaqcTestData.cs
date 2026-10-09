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

    internal static double[] Vafs(string dataset, int major, int minor) =>
        Load(dataset).Where(r => r.Major == major && r.Minor == minor).Select(r => r.Vaf).ToArray();
}
