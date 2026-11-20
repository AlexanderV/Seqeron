// GENOMIC-TANDEM-001 (B09 A2-34) — GenomicAnalyzer.FindTandemRepeats against MUMmer 3.23 `exact-tandems`
// (repeat-match -t -n L | sort | tandem-repeat.awk), and the canonical-rotation option.
// Oracle: scripts/analysis/generate_exact_tandems_oracle.py (3630 (sequence, L) rows; the least rotation of each unit
// is computed there by brute force). Evidence: docs/Evidence/GENOMIC-TANDEM-001-Evidence.md.
//
// exact-tandems lists, per start position, the tandem repeat with the SMALLEST unit length whose maximal match
// (extent - unit) is at least L and that has >= 2 complete copies. FindTandemRepeats lists every maximal run with a
// primitive unit; the two agree after keeping, per start, the run with the smallest unit among those with extent - unit >= L.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_ExactTandems_Oracle_Tests
{
    private sealed record Row(string Sequence, int MinMatch, List<(int Start, int Unit, int Extent, string Canonical)> Hits);

    // One parsed copy per fixture (audit A3-15): the oracle was re-read and re-parsed on every call, i.e. in every test.
    private static readonly Lazy<List<Row>> OracleRows = new(ReadOracleRows);

    private static List<Row> ReadOracle() => OracleRows.Value;

    private static List<Row> ReadOracleRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.exact_tandems_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            var hits = new List<(int, int, int, string)>();
            if (f[2] != "-")
                foreach (string h in f[2].Split(';'))
                {
                    var p = h.Split(':');
                    hits.Add((int.Parse(p[0]), int.Parse(p[1]), int.Parse(p[2]), p[3]));
                }
            rows.Add(new Row(f[0], int.Parse(f[1]), hits));
        }
        return rows;
    }

    // Extent of the maximal period-p interval that starts at a (the run S[a..e) with S[x] = S[x+p]).
    private static int Extent(string s, int a, int p)
    {
        int e = a + p;
        while (e < s.Length && s[e] == s[e - p]) e++;
        return e - a;
    }

    [Test]
    public void FindTandemRepeats_MatchesExactTandems_AllRows()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.EqualTo(3630));
        int nonEmpty = 0;
        foreach (var row in rows)
        {
            var ours = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(row.Sequence), 1, 2)
                .Select(t => (Start: t.Position, Unit: t.Unit.Length, Extent: Extent(row.Sequence, t.Position, t.Unit.Length)))
                .Where(t => t.Extent - t.Unit >= row.MinMatch)
                .GroupBy(t => t.Start)
                .Select(g => g.OrderBy(t => t.Unit).First())
                .OrderBy(t => t.Start)
                .ToList();
            var expected = row.Hits.Select(h => (h.Start, h.Unit, h.Extent)).ToList();
            Assert.That(ours, Is.EqualTo(expected), $"{row.Sequence} L={row.MinMatch}");
            if (expected.Count > 0) nonEmpty++;
        }
        Assert.That(nonEmpty, Is.GreaterThan(3000));
    }

    [Test]
    public void CanonicalRotation_IsLeastRotationOfTheOracleUnit_AllRows()
    {
        int checkedUnits = 0;
        foreach (var row in ReadOracle().Where(r => r.MinMatch == 1))
        {
            var runs = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(row.Sequence), 1, 2, int.MaxValue, canonicalRotation: true)
                .ToDictionary(t => (t.Position, t.Unit.Length));
            foreach (var h in row.Hits)
            {
                var t = runs[(h.Start, h.Unit)];
                Assert.That(t.Unit, Is.EqualTo(h.Canonical), $"{row.Sequence} @{h.Start}");
                Assert.That(t.OriginalUnit, Is.EqualTo(row.Sequence.Substring(h.Start, h.Unit)));
                Assert.That(t.FullSequence, Is.EqualTo(row.Sequence.Substring(h.Start, h.Unit * t.Repetitions)));
                checkedUnits++;
            }
        }
        Assert.That(checkedUnits, Is.GreaterThan(5000));
    }

    [Test]
    public void CanonicalRotation_OnlyRewritesTheUnit()
    {
        const string s = "GATATGATGATGCCAGTAGTAGT";
        var plain = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(s), 2, 2, int.MaxValue).ToList();
        var canon = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(s), 2, 2, int.MaxValue, canonicalRotation: true).ToList();
        Assert.That(canon.Select(t => (t.Position, t.Repetitions)), Is.EqualTo(plain.Select(t => (t.Position, t.Repetitions))));
        Assert.That(plain.Select(t => t.Unit), Does.Contain("TGA").Or.Contain("ATG"));
        Assert.That(canon.Select(t => t.Unit), Is.All.Matches<string>(u => IsLeastRotation(u)));
        Assert.That(GenomicAnalyzer.FindTandemRepeats(new DnaSequence(s), 2, 2, int.MaxValue, canonicalRotation: false).Select(t => t.Unit),
            Is.EqualTo(plain.Select(t => t.Unit)));
    }

    [Test]
    public void CanonicalRotation_StringOverload_WithN_MatchesDnaOverloadOnAcgt()
    {
        var withN = GenomicAnalyzer.FindTandemRepeats("ACNACNACNTGATGATGAT", 2, 2, int.MaxValue, canonicalRotation: true).ToList();
        Assert.That(withN.Select(t => (t.Unit, t.Position, t.Repetitions)), Is.EqualTo(new[] { ("ATG", 9, 3) }));
        Assert.That(withN[0].OriginalUnit, Is.EqualTo("TGA"));
        Assert.That(withN[0].Rotation, Is.EqualTo(2));
    }

    private static bool IsLeastRotation(string u)
        => Enumerable.Range(0, u.Length).All(i => string.CompareOrdinal(u, u[i..] + u[..i]) <= 0);
}
