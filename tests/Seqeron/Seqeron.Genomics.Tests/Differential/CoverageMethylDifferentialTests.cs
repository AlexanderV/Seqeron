// 08_DIFFERENTIAL_TESTING rows 142, 185. Independent oracles: a manual best-match read placement +
// per-position depth, and a manual per-CpG bisulfite methylation level (methylated / total).

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class CoverageMethylDifferentialTests
{
    private const double Tol = 1e-12;

    // ---- Row 142: ASSEMBLY-COVER-001 — CalculateCoverage vs manual placement + depth ----

    private static int[] CoverageOracle(string reference, string[] reads, int minOverlap)
    {
        var cov = new int[reference.Length];
        foreach (var read in reads)
        {
            int bestPos = -1, bestScore = minOverlap - 1;
            for (int pos = 0; pos <= reference.Length - read.Length; pos++)
            {
                int m = 0;
                for (int i = 0; i < read.Length; i++)
                    if (char.ToUpperInvariant(reference[pos + i]) == char.ToUpperInvariant(read[i])) m++;
                if (m > bestScore) { bestScore = m; bestPos = pos; }
            }
            if (bestPos >= 0)
                for (int i = bestPos; i < bestPos + read.Length && i < reference.Length; i++) cov[i]++;
        }
        return cov;
    }

    [Test]
    [Category("ASSEMBLY-COVER-001")]
    public void CalculateCoverage_MatchesManualPlacement()
    {
        const string reference = "ACGTACGTAC";
        var reads = new[] { "ACGT", "GTAC", "ACGTAC" };
        const int minOverlap = 3;
        Assert.That(SequenceAssembler.CalculateCoverage(reference, reads, minOverlap),
            Is.EqualTo(CoverageOracle(reference, reads, minOverlap)));
    }

    // ---- Row 185: EPIGEN-METHYL-001 — per-CpG methylation level = methylated / total ----

    [Test]
    [Category("EPIGEN-METHYL-001")]
    public void MethylationFromBisulfite_MatchesMethylatedOverTotal()
    {
        const string reference = "ACGAACGT"; // CpG cytosines at positions 1 and 5
        var reads = new (string, int)[] { ("CG", 1), ("TG", 1), ("CG", 5) };

        // Independent oracle.
        var cpg = new List<int>();
        for (int i = 0; i < reference.Length - 1; i++)
            if (reference[i] == 'C' && reference[i + 1] == 'G') cpg.Add(i);
        var data = cpg.ToDictionary(p => p, _ => (meth: 0, total: 0));
        foreach (var (read, start) in reads)
            for (int i = 0; i < read.Length && start + i < reference.Length - 1; i++)
            {
                int refPos = start + i;
                if (!data.ContainsKey(refPos)) continue;
                if (read[i] == 'C') data[refPos] = (data[refPos].meth + 1, data[refPos].total + 1);
                else if (read[i] == 'T') data[refPos] = (data[refPos].meth, data[refPos].total + 1);
            }
        var expected = cpg.Where(p => data[p].total > 0)
            .Select(p => (p, (double)data[p].meth / data[p].total, data[p].total)).ToList();

        var actual = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(reference, reads)
            .Select(s => (s.Position, s.MethylationLevel, s.Coverage)).ToList();

        Assert.That(actual.Count, Is.EqualTo(expected.Count));
        for (int k = 0; k < actual.Count; k++)
        {
            Assert.That(actual[k].Position, Is.EqualTo(expected[k].p), $"pos[{k}]");
            Assert.That(actual[k].MethylationLevel, Is.EqualTo(expected[k].Item2).Within(Tol), $"level[{k}]");
            Assert.That(actual[k].Coverage, Is.EqualTo(expected[k].Item3), $"cov[{k}]");
        }
    }
}
