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
    public void Oracle_Infoalign_PerSequenceStatistics_AgainstReferenceAndConsensus()
    {
        var failures = new List<string>();
        int checkedRows = 0, consensusRuns = 0, similarCounted = 0;
        foreach (var r in Oracle.Value.Where(r => r.Kind == "infoalign"))
        {
            int refNo = int.Parse(r.Param, CultureInfo.InvariantCulture);
            var stats = MultipleAlignmentStatistics.CalculateSequenceStatistics(r.Rows, refNo == 0 ? null : refNo - 1);
            var expected = r.Expected.Split(';');
            if (refNo == 0) consensusRuns++;
            for (int i = 0; i < expected.Length; i++)
            {
                var e = expected[i].Split(',');
                var a = stats[i];
                string got = $"{a.SequenceLength},{a.AlignmentLength},{a.GapRuns},{a.GapLength},{a.Identities},{a.Similarities},{a.Differences}";
                string want = string.Join(',', e.Take(7));
                checkedRows++;
                if (a.Similarities > 0) similarCounted++;
                if (got != want || Math.Abs(a.PercentChange - Parse(e[7])) > 1e-4)
                    failures.Add($"ref {refNo} row {i}: {got} {a.PercentChange:F6} expected {expected[i]} in {string.Join(',', r.Rows)}");
            }
        }
        Assert.Multiple(() =>
        {
            Assert.That(checkedRows, Is.GreaterThan(3000));
            Assert.That(consensusRuns, Is.GreaterThanOrEqualTo(150), "consensus reference exercised (cons via MotifFinder.GenerateEmbossConsensus)");
            Assert.That(similarCounted, Is.GreaterThan(100), "similar (positive-score) residues exercised");
            Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
        });
    }

    [Test]
    public void Oracle_Distmat_Uncorrected_And_JukesCantor_WithGapWeight()
    {
        var failures = new List<string>();
        int matrices = 0, nanEntries = 0;
        foreach (var r in Oracle.Value.Where(r => r.Kind == "distmat"))
        {
            var p = r.Param.Split(':');
            var method = p[0] == "0" ? MsaDistanceMethod.Uncorrected : MsaDistanceMethod.JukesCantor;
            double gw = double.Parse(p[1], CultureInfo.InvariantCulture);
            var d = MultipleAlignmentStatistics.CalculateDistanceMatrix(r.Rows, method, ConsensusResidueType.Auto, gw);
            var expected = r.Expected.Split(',').Select(Parse).ToArray();
            int n = r.Rows.Length, idx = 0;
            matrices++;
            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++, idx++)
                {
                    double got = d[i, j] * 100;
                    if (double.IsNaN(expected[idx]))
                    {
                        nanEntries++;
                        if (!double.IsNaN(got)) failures.Add($"{r.Param} [{i},{j}] expected nan got {got}");
                    }
                    else if (double.IsInfinity(expected[idx]))
                    {
                        if (!double.IsPositiveInfinity(got)) failures.Add($"{r.Param} [{i},{j}] expected inf got {got}");
                    }
                    // distmat computes in single precision (float D, float log): two decimals plus a 1e-4 relative allowance near saturation
                    else if (Math.Abs(got - expected[idx]) > 0.0051 + 1e-4 * Math.Abs(expected[idx]) || d[i, j] != d[j, i])
                    {
                        failures.Add($"{r.Param} [{i},{j}] expected {expected[idx]} got {got:F4} in {string.Join(',', r.Rows)}");
                    }
                }
                if (d[i, i] != 0) failures.Add("non-zero diagonal");
            }
        }
        Assert.Multiple(() =>
        {
            Assert.That(matrices, Is.GreaterThanOrEqualTo(600));
            Assert.That(nanEntries, Is.GreaterThan(0), "Jukes-Cantor saturation (log of a non-positive number) exercised");
            Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
        });
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

    [Test]
    public void Calculate_ColumnsGapsAndSummary()
    {
        var rows = new[] { "ACGT-ACGTAC", "ACGTTACG-AC", "AGGT-ACGTAC", "A-GTTTCGTAA" };
        var s = MultipleAlignmentStatistics.Calculate(rows, ConsensusResidueType.Nucleotide);
        Assert.Multiple(() =>
        {
            Assert.That(s.SequenceCount, Is.EqualTo(4));
            Assert.That(s.AlignmentLength, Is.EqualTo(11));
            Assert.That(s.Columns[0].Conservation, Is.EqualTo('*'));
            Assert.That(s.Columns[0].MostCommonResidue, Is.EqualTo('A'));
            Assert.That(s.Columns[0].MostCommonCount, Is.EqualTo(4));
            Assert.That(s.Columns[1].Counts['C'], Is.EqualTo(2));
            Assert.That(s.Columns[1].GapCount, Is.EqualTo(1));
            Assert.That(s.Columns[1].GapFraction, Is.EqualTo(0.25));
            Assert.That(s.Columns[4].GapCount, Is.EqualTo(2));
            Assert.That(s.Columns[4].MostCommonResidue, Is.EqualTo('T'));
            Assert.That(s.GapPercent, Is.EqualTo(100.0 * 4 / 44).Within(1e-12));
            Assert.That(s.GapFreeColumns, Is.EqualTo(8));
            // distmat 6.6.0 uncorrected (x100, printed with two decimals; the program computes in single precision): s1/s2 0.00, s1/s3 10.00, s1/s4 22.22, s2/s3 11.11, s2/s4 22.22, s3/s4 22.22
            var d = MultipleAlignmentStatistics.CalculateDistanceMatrix(rows);
            Assert.That(d[0, 1] * 100, Is.EqualTo(0).Within(1e-9));
            Assert.That(d[0, 2] * 100, Is.EqualTo(10.0).Within(0.005));
            Assert.That(d[1, 2] * 100, Is.EqualTo(11.11).Within(0.005));
            Assert.That(d[2, 3] * 100, Is.EqualTo(22.22).Within(0.005));
            Assert.That(s.PairwiseIdentity[0, 2], Is.EqualTo(90.0).Within(1e-9));
            // Biopython identity (gap columns count as differences): s1/s3 9 matches / 11 columns
            Assert.That(MultipleAlignmentStatistics.CalculatePairwiseIdentityMatrix(rows, 1)[0, 2], Is.EqualTo(100.0 * 9 / 11).Within(1e-9));
        });
        // infoalign -refseq 1 (EMBOSS 6.6.0): s2 10 11 1 1 9 0 1 18.181818
        var info = MultipleAlignmentStatistics.CalculateSequenceStatistics(rows, 0)[1];
        Assert.Multiple(() =>
        {
            Assert.That((info.SequenceLength, info.AlignmentLength, info.GapRuns, info.GapLength, info.Identities, info.Similarities, info.Differences),
                Is.EqualTo((10, 11, 1, 1, 9, 0, 1)));
            Assert.That(info.PercentChange, Is.EqualTo(18.181818).Within(1e-5));
        });
    }

    [Test]
    public void MultipleAlignmentResult_Overload_And_Validation()
    {
        var result = new MultipleAlignmentResult(new[] { "AC-T", "ACGT" }, "ACGT", 0);
        Assert.Multiple(() =>
        {
            Assert.That(MultipleAlignmentStatistics.Calculate(result).SequenceCount, Is.EqualTo(2));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentStatistics.Calculate(new[] { "AC", "ACG" }));
            Assert.Throws<ArgumentNullException>(() => MultipleAlignmentStatistics.Calculate((MultipleAlignmentResult)null!));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentStatistics.CalculateDistanceMatrix(new[] { "ACGT" }));
            Assert.Throws<ArgumentOutOfRangeException>(() => MultipleAlignmentStatistics.CalculateSequenceStatistics(new[] { "AC", "AG" }, 5));
            Assert.That(MultipleAlignmentStatistics.Calculate(Array.Empty<string>()).AlignmentLength, Is.EqualTo(0));
            Assert.That(MultipleAlignmentStatistics.Calculate(new[] { "ACGT" }).MeanPairwiseIdentity, Is.EqualTo(100.0));
        });
    }
}
