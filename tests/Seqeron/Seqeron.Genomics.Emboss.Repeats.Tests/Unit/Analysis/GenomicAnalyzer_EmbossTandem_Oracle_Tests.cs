// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Repeats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
// GENOMIC-TANDEM-001 (B09 A3-8) — EMBOSS etandem / equicktandem 6.6.0 ports (EmbossRepeatPrograms.FindTandemRepeatsEtandem / FindTandemRepeatsEquicktandem).
// Oracle: scripts/analysis/generate_emboss_tandem_oracle.py — every row is the original-format output (-origfile) of the installed binary on one (sequence, parameters) case.
// Evidence: docs/Evidence/GENOMIC-TANDEM-001-Evidence.md.

using Seqeron.Genomics.Emboss.Repeats;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_EmbossTandem_Oracle_Tests
{
    private sealed record Row(string Program, string Params, string Sequence, string Hits);

    // One parsed copy per fixture (audit A3-15): the oracle was re-read and re-parsed on every call, i.e. in every test.
    private static readonly Lazy<List<Row>> OracleRows = new(ReadOracleRows);

    private static List<Row> ReadOracle() => OracleRows.Value;

    private static List<Row> ReadOracleRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.emboss_tandem_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], f[2], f[3]));
        }
        return rows;
    }

    // The programs' original output line: "%6d %10d %10d %2d %3d %5.1f consensus" / "%6d %10d %10d %2d %3d" -> compared field by field.
    private static string FormatEtandem(EmbossTandemRepeat h)
        => string.Join(",", h.Score, h.Start, h.End, h.RepeatSize, h.RepeatCount, CFormat.Fixed(h.Identity, 1), h.Consensus);

    private static string FormatEquick(EmbossTandemRepeat h)
        => string.Join(",", h.Score, h.Start, h.End, h.RepeatSize, h.RepeatCount);

    private static string Join(IEnumerable<string> rows)
    {
        string s = string.Join("|", rows);
        return s.Length == 0 ? "-" : s;
    }

    [Test]
    public void FindTandemRepeatsEtandem_MatchesEtandem_AllRows()
    {
        var rows = ReadOracle().Where(r => r.Program == "etandem").ToList();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(1000));
        int nonEmpty = 0, mismatchRows = 0, uniformRows = 0, multiHit = 0, lower = 0, withN = 0;
        foreach (var r in rows)
        {
            var p = r.Params.Split(',').Select(int.Parse).ToArray();
            var hits = EmbossRepeatPrograms.FindTandemRepeatsEtandem(r.Sequence, p[0], p[1], p[2], p[3] == 1, p[4] == 1);
            Assert.That(Join(hits.Select(FormatEtandem)), Is.EqualTo(r.Hits), $"etandem {r.Params} {r.Sequence}");
            if (hits.Count > 0) nonEmpty++;
            if (hits.Count > 1) multiHit++;
            if (p[3] == 1) mismatchRows++;
            if (p[4] == 1) uniformRows++;
            if (r.Sequence.Any(char.IsLower)) lower++;
            if (r.Sequence.Contains('N') || r.Sequence.Contains('n')) withN++;
        }

        Assert.That(nonEmpty, Is.GreaterThan(700));
        Assert.That(multiHit, Is.GreaterThan(300));
        Assert.That(mismatchRows, Is.GreaterThan(300));
        Assert.That(uniformRows, Is.GreaterThan(250));
        Assert.That(lower, Is.GreaterThan(200));
        Assert.That(withN, Is.GreaterThan(200));
    }

    [Test]
    public void FindTandemRepeatsEquicktandem_MatchesEquicktandem_AllRows()
    {
        var rows = ReadOracle().Where(r => r.Program == "equicktandem").ToList();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(1000));
        int nonEmpty = 0, multiHit = 0;
        foreach (var r in rows)
        {
            var p = r.Params.Split(',').Select(int.Parse).ToArray();
            var hits = EmbossRepeatPrograms.FindTandemRepeatsEquicktandem(r.Sequence, p[0], p[1]);
            Assert.That(Join(hits.Select(FormatEquick)), Is.EqualTo(r.Hits), $"equicktandem {r.Params} {r.Sequence}");
            if (hits.Count > 0) nonEmpty++;
            if (hits.Count > 1) multiHit++;
        }

        Assert.That(nonEmpty, Is.GreaterThan(500));
        Assert.That(multiHit, Is.GreaterThan(200));
    }

    [Test]
    public void Etandem_WorkedExample_ReportsTheTwoArraysHighestScoreFirst()
    {
        // run by hand: etandem -minrepeat 2 -maxrepeat 8 -threshold 10 on this sequence printed
        //   "    33         61         96  3  12 100.0 cag" and "    26        137        172  6   6  94.4 acgtta"
        const string seq = "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACACGTTAACGTTAACGTTAACGTTAACGTTAacgtNN";
        var hits = EmbossRepeatPrograms.FindTandemRepeatsEtandem(seq, 2, 8, 10);
        var rows = hits.Select(h => (h.Score, h.Start, h.End, h.RepeatSize, h.RepeatCount, CFormat.Fixed(h.Identity, 1), h.Consensus)).ToArray();
        Assert.That(rows, Is.EqualTo(new[] { (33, 61, 96, 3, 12, "100.0", "cag"), (26, 137, 172, 6, 6, "94.4", "acgtta") }));
    }

    [Test]
    public void Equicktandem_ExactRepeat_ScoreIsLengthMinusGap()
    {
        // a perfect array of 8 copies of ACG (equicktandem -maxrepeat 6 -threshold 10 printed "21 1 24 3 8"): gap 3 gives 21 matches; gap 6 is masked
        var hits = EmbossRepeatPrograms.FindTandemRepeatsEquicktandem("ACGACGACGACGACGACGACGACG", 6, 10);
        Assert.That(hits, Has.Count.EqualTo(1));
        Assert.That((hits[0].Score, hits[0].Position, hits[0].RepeatSize, hits[0].RepeatCount), Is.EqualTo((21, 0, 3, 8)));
    }

    [Test]
    public void Arguments_AreValidated()
    {
        Assert.Throws<ArgumentNullException>(() => EmbossRepeatPrograms.FindTandemRepeatsEtandem(null!, 2, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossRepeatPrograms.FindTandemRepeatsEtandem("ACGT", 1, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossRepeatPrograms.FindTandemRepeatsEtandem("ACGT", 3, 2));
        Assert.Throws<ArgumentNullException>(() => EmbossRepeatPrograms.FindTandemRepeatsEquicktandem(null!, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossRepeatPrograms.FindTandemRepeatsEquicktandem("ACGT", 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossRepeatPrograms.FindTandemRepeatsEquicktandem("ACGT", 3, 0));
        Assert.That(EmbossRepeatPrograms.FindTandemRepeatsEtandem("", 2, 5), Is.Empty);
        Assert.That(EmbossRepeatPrograms.FindTandemRepeatsEquicktandem("", 5), Is.Empty);
    }
}
