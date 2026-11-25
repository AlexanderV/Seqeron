// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
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

    [Test]
    public void Oracle_Infoalign_PerSequenceStatistics_AgainstReferenceAndConsensus()
    {
        var failures = new List<string>();
        int checkedRows = 0, consensusRuns = 0, similarCounted = 0;
        foreach (var r in Oracle.Value.Where(r => r.Kind == "infoalign"))
        {
            int refNo = int.Parse(r.Param, CultureInfo.InvariantCulture);
            var stats = EmbossAlignmentStatistics.CalculateSequenceStatistics(r.Rows, refNo == 0 ? null : refNo - 1);
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
            Assert.That(consensusRuns, Is.GreaterThanOrEqualTo(150), "consensus reference exercised (cons via EmbossConsensus.GenerateEmbossConsensus)");
            Assert.That(similarCounted, Is.GreaterThan(100), "similar (positive-score) residues exercised");
            Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
        });
    }

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
}
