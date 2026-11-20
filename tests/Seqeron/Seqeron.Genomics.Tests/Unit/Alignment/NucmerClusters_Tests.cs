using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A3-9a, B09 F150): <see cref="NucmerClusters.Cluster"/> = the <c>mgaps</c> program of MUMmer 3.23 (the clustering between <c>mummer</c> and <c>postnuc</c> inside <c>nucmer</c>).
/// Oracle (TestData/Alignment/mgaps_oracle.tsv, scripts/alignment/generate_mgaps_oracle.py): 758 match lists - the forward-strand matches of <c>mummer -maxmatch -n -l 6..12</c> on related random pairs
/// (repeats, indels, overlapping matches on different diagonals) and shuffled synthetic lists - with random <c>-l -s -d -f [-e]</c>, and the clusters the installed <c>mgaps</c> prints for them.
/// </summary>
[TestFixture]
public class NucmerClusters_Tests
{
    private sealed record Row(string Kind, NucmerClusterOptions Options, (long, long, long)[] Matches, string Clusters);

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.mgaps_oracle.tsv";
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
                MinClusterLength = long.Parse(f[1], CultureInfo.InvariantCulture),
                MaxSeparation = long.Parse(f[2], CultureInfo.InvariantCulture),
                DiagonalDifference = int.Parse(f[3], CultureInfo.InvariantCulture),
                DiagonalFraction = double.Parse(f[4], CultureInfo.InvariantCulture),
                UseExtent = f[5] == "1",
            };
            var matches = f[6].Split(',').Select(m => m.Split(':')).Select(p => (long.Parse(p[0], CultureInfo.InvariantCulture), long.Parse(p[1], CultureInfo.InvariantCulture), long.Parse(p[2], CultureInfo.InvariantCulture))).ToArray();
            rows.Add(new Row(f[0], options, matches, f[7]));
        }
        return rows;
    }

    private static string Format(IReadOnlyList<NucmerCluster> clusters) =>
        string.Join("|", clusters.Select(c => string.Join(";", c.Matches.Select(m => $"{m.Start1}:{m.Start2}:{m.Length}:{m.Overlap}"))));

    [Test]
    public void Cluster_EqualsMgaps_AllRows()
    {
        var rows = Oracle.Value;
        Assert.That(rows, Has.Count.EqualTo(758));
        var failures = new List<string>();
        int withClusters = 0, multi = 0, overlapped = 0, extent = 0, mummer = 0, many = 0, unions = 0;
        foreach (var r in rows)
        {
            var actual = Format(NucmerClusters.Cluster(r.Matches, r.Options));
            if (actual != r.Clusters)
                failures.Add($"{r.Kind} l={r.Options.MinClusterLength} s={r.Options.MaxSeparation} d={r.Options.DiagonalDifference} f={r.Options.DiagonalFraction} e={r.Options.UseExtent}\n{string.Join(",", r.Matches)}\n expected {r.Clusters}\n actual   {actual}");
            if (r.Clusters.Length > 0) withClusters++;
            if (r.Clusters.Contains('|')) multi++;
            if (r.Clusters.Split('|', ';', ',').Any(m => !m.EndsWith(":0", StringComparison.Ordinal))) overlapped++;
            if (r.Options.UseExtent) extent++;
            if (r.Kind == "mummer") mummer++;
            if (r.Matches.Length > 20) many++;
            if (r.Clusters.Contains(';')) unions++;
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
        Assert.That(withClusters, Is.GreaterThan(500));
        Assert.That(multi, Is.GreaterThan(150));
        Assert.That(overlapped, Is.GreaterThan(30));
        Assert.That(extent, Is.GreaterThan(80));
        Assert.That(mummer, Is.GreaterThan(300));
        Assert.That(many, Is.GreaterThan(150));
        Assert.That(unions, Is.GreaterThan(300));
    }

    [Test]
    public void WorkedValues_RunOnTheBinary()
    {
        // mummer -mum -b -l 20 -n | mgaps -l 65 -s 90 -d 5 -f 0.12 on two shared segments of 250 and 240 bases separated by a 20 base insertion in the query:
        // "51 1 250 none - -" and "321 271 240 none 20 20" are two matches of ONE cluster (separation 20 on both sequences)
        var c = NucmerClusters.Cluster([(51, 1, 250), (321, 271, 240)]);
        Assert.That(c, Has.Count.EqualTo(1));
        Assert.That(c[0].Matches, Is.EqualTo(new[] { new NucmerClusterMatch(51, 1, 250, 0), new NucmerClusterMatch(321, 271, 240, 0) }));
        // two matches farther apart than -s are separate clusters; each must reach -c
        Assert.That(NucmerClusters.Cluster([(1, 1, 70), (500, 500, 70)]), Has.Count.EqualTo(2));
        Assert.That(NucmerClusters.Cluster([(1, 1, 30), (500, 500, 30)]), Is.Empty);
        // a coordinate shift does not change the clusters
        Assert.That(NucmerClusters.Cluster([(0, 0, 70), (80, 82, 70)])[0].Matches.Select(m => (m.Length, m.Overlap)),
            Is.EqualTo(NucmerClusters.Cluster([(1, 1, 70), (81, 83, 70)])[0].Matches.Select(m => (m.Length, m.Overlap))));
    }

    [Test]
    public void Validation()
    {
        Assert.Throws<ArgumentNullException>(() => NucmerClusters.Cluster(null!));
        Assert.Throws<ArgumentException>(() => NucmerClusters.Cluster([(1, 1, 0)]));
        Assert.Throws<ArgumentException>(() => NucmerClusters.Cluster([(1, 1, 10), (1, 1, 20)]));
        Assert.That(NucmerClusters.Cluster([]), Is.Empty);
    }
}
