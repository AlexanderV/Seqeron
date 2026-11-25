using System.Globalization;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-STATS-001 (audit A1-9, B09 F43): statistics of multiple alignments — per-sequence comparison with a reference / the
/// consensus (EMBOSS infoalign 6.6.0 <c>infoalign_Compare</c>), pairwise distance / identity matrices (EMBOSS distmat 6.6.0
/// uncorrected + Jukes-Cantor; Biopython 1.88 <c>DistanceCalculator('identity')</c>), ClustalW 2.1 conservation line.
/// Oracle: TestData/Alignment/msa_stats_oracle.tsv (scripts/alignment/generate_msa_stats_oracle.py): 150 random gapped DNA/IUPAC and
/// protein alignments (3-8 rows, leading/trailing/internal gaps).
/// </summary>
[TestFixture]
public class MultipleAlignmentStatistics_Tests
{
    private sealed record Row(string Kind, string[] Rows, string Param, string Expected);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.msa_stats_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1].Split(','), f[2], f[3]));
        }
        return rows;
    });

    private static double Parse(string s)
    {
        if (s.Contains("nan", StringComparison.OrdinalIgnoreCase)) return double.NaN;
        if (s.Equals("inf", StringComparison.OrdinalIgnoreCase)) return double.PositiveInfinity;
        return double.Parse(s, CultureInfo.InvariantCulture);
    }

    [Test]
    public void Oracle_ClustalConservationLine_DnaAndProtein()
    {
        var failures = new List<string>();
        int strong = 0, weak = 0, star = 0;
        foreach (var r in Oracle.Value.Where(r => r.Kind == "clustal"))
        {
            var type = r.Param == "protein" ? ConsensusResidueType.Protein : ConsensusResidueType.Nucleotide;
            string got = MultipleAlignmentStatistics.ClustalConservationLine(r.Rows, type);
            string want = r.Expected.Replace('|', ' ').TrimEnd();
            strong += got.Count(c => c == ':');
            weak += got.Count(c => c == '.');
            star += got.Count(c => c == '*');
            if (got.TrimEnd() != want) failures.Add($"{r.Param}: {got.TrimEnd()} expected {want} in {string.Join(',', r.Rows)}");
        }
        Assert.Multiple(() =>
        {
            Assert.That(strong, Is.GreaterThan(30));
            Assert.That(weak, Is.GreaterThan(30));
            Assert.That(star, Is.GreaterThan(500));
            Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
        });
    }

    /// <summary>Biopython 1.88 <c>DistanceCalculator('identity')</c> compares characters (skip_letters = ()), so it equals the
    /// gap-weight-1 identity (matches / alignment length) plus the columns where both rows have '-'.</summary>
    [Test]
    public void Oracle_PairwiseIdentity_EqualsBiopythonDistanceCalculatorIdentity()
    {
        var failures = new List<string>();
        int matrices = 0, bothGapColumns = 0;
        foreach (var r in Oracle.Value.Where(r => r.Kind == "bioid"))
        {
            var id = MultipleAlignmentStatistics.CalculatePairwiseIdentityMatrix(r.Rows, gapWeight: 1);
            var expected = r.Expected.Split(',').Select(Parse).ToArray();
            int len = r.Rows[0].Length, idx = 0;
            matrices++;
            for (int i = 0; i < r.Rows.Length; i++)
            {
                for (int j = i + 1; j < r.Rows.Length; j++, idx++)
                {
                    int bothGap = Enumerable.Range(0, len).Count(k => r.Rows[i][k] == '-' && r.Rows[j][k] == '-');
                    bothGapColumns += bothGap;
                    double got = id[i, j] + 100.0 * bothGap / len;
                    if (Math.Abs(got - expected[idx]) > 1e-4 || id[i, j] != id[j, i])
                        failures.Add($"[{i},{j}] expected {expected[idx]} got {got}");
                }
                Assert.That(id[i, i], Is.EqualTo(100.0));
            }
        }
        Assert.That(matrices, Is.GreaterThanOrEqualTo(70));
        Assert.That(bothGapColumns, Is.GreaterThan(0), "gap/gap columns exercised");
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
    }

    /// <summary>ClustalW 2.1 -convert of this protein alignment prints the conservation line "*:* ** ** *..".</summary>
    [Test]
    public void ClustalConservationLine_ProteinExample_IsTheClustalwLine()
    {
        var rows = new[] { "MKVSLAKAG-HST", "MKVALAWAGQHSA", "MRVGLADAGNHNV", "MKVCLAEAG-HTA" };
        Assert.That(MultipleAlignmentStatistics.ClustalConservationLine(rows, ConsensusResidueType.Protein), Is.EqualTo("*:* ** ** *.."));
    }
}
