// ALIGN-MULTI-001 (B09 audit A3-9a, F150) — GenomicAnalyzer.FindNucmerClusters = the `mummer | mgaps` pipeline of MUMmer 3.23 nucmer.
// Oracle: scripts/alignment/generate_nucmer_clusters_oracle.py — 700 related random pairs (substitutions, indels, repeats, inverted segments, N), mummer -mum | -mumreference | -maxmatch -b -l L -n piped into
// mgaps -l c -s g -d D -f F, the clusters of both strands. Evidence: docs/Evidence/ALIGN-MULTI-001-Evidence.md.

using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_NucmerClusters_Tests
{
    private sealed record Row(string Reference, string Query, string Mode, int MinMatch, NucmerClusterOptions Options, string Clusters);

    private static readonly Lazy<List<Row>> OracleRows = new(ReadOracle);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.nucmer_clusters_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            var options = new NucmerClusterOptions
            {
                MinClusterLength = long.Parse(f[4], CultureInfo.InvariantCulture),
                MaxSeparation = long.Parse(f[5], CultureInfo.InvariantCulture),
                DiagonalDifference = int.Parse(f[6], CultureInfo.InvariantCulture),
                DiagonalFraction = double.Parse(f[7], CultureInfo.InvariantCulture),
            };
            rows.Add(new Row(f[0], f[1], f[2], int.Parse(f[3], CultureInfo.InvariantCulture), options, f.Length > 8 ? f[8] : string.Empty));
        }
        return rows;
    }

    private static CommonRegionMode ModeOf(string m) => m switch
    {
        "mum" => CommonRegionMode.MaximalUnique,
        "mumreference" => CommonRegionMode.MaximalUniqueReference,
        _ => CommonRegionMode.MaximalExactMatches,
    };

    private static string Format(IReadOnlyList<StrandedNucmerCluster> clusters) =>
        string.Join("|", clusters.Select(c => (c.IsReverseComplement ? "R/" : "F/") + string.Join(";", c.Cluster.Matches.Select(m => $"{m.Start1}:{m.Start2}:{m.Length}:{m.Overlap}"))));

    [Test]
    public void FindNucmerClusters_EqualsMummerPipedIntoMgaps_AllRows()
    {
        var rows = OracleRows.Value;
        Assert.That(rows, Has.Count.EqualTo(700));
        var failures = new List<string>();
        int withClusters = 0, reverse = 0, multiMatch = 0, maxmatch = 0, mumref = 0, mum = 0, overlapped = 0;
        foreach (var r in rows)
        {
            var actual = Format(GenomicAnalyzer.FindNucmerClusters(r.Reference, r.Query, r.MinMatch, ModeOf(r.Mode), r.Options));
            if (actual != r.Clusters)
                failures.Add($"{r.Mode} l={r.MinMatch}\n expected {r.Clusters}\n actual   {actual}");
            if (r.Clusters.Length > 0) withClusters++;
            if (r.Clusters.Contains("R/")) reverse++;
            if (r.Clusters.Contains(';')) multiMatch++;
            if (r.Mode == "maxmatch") maxmatch++;
            if (r.Mode == "mumreference") mumref++;
            if (r.Mode == "mum") mum++;
            if (r.Clusters.Split('|', ';', ',').Any(m => m.Length > 0 && !m.EndsWith(":0", StringComparison.Ordinal))) overlapped++;
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
        Assert.That(withClusters, Is.GreaterThan(600));
        Assert.That(reverse, Is.GreaterThan(200));
        Assert.That(multiMatch, Is.GreaterThan(300));
        Assert.That(maxmatch, Is.GreaterThan(100));
        Assert.That(mumref, Is.GreaterThan(100));
        Assert.That(mum, Is.GreaterThan(300));
        Assert.That(overlapped, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void ForwardOnly_AndDefaults()
    {
        var r = OracleRows.Value.First(x => x.Clusters.Contains("F/") && x.Clusters.Contains("R/"));
        var forward = GenomicAnalyzer.FindNucmerClusters(r.Reference, r.Query, r.MinMatch, ModeOf(r.Mode), r.Options, StrandSelection.Forward);
        Assert.That(forward.All(c => !c.IsReverseComplement), Is.True);
        Assert.That(Format(forward), Is.EqualTo(string.Join("|", r.Clusters.Split('|').Where(c => c.StartsWith("F/", StringComparison.Ordinal)))));
        Assert.That(GenomicAnalyzer.FindNucmerClusters("ACGT", "ACGT"), Is.Empty);
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindNucmerClusters(null!, "ACGT"));
    }
}
