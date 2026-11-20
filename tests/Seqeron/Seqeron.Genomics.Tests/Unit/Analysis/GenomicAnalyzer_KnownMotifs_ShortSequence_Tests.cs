// GENOMIC-MOTIFS-001 (audit A2-18, B09 F91) - the "overhanging hit" of EMBOSS fuzznuc 6.6.0 on a sequence shorter than the pattern, decided against the binary.
// Oracle: TestData/Analysis/fuzznuc_short_sequence_oracle.tsv (scripts/analysis/generate_fuzznuc_short_sequence_oracle.py, 300 runs of fuzznuc -pmismatch k on a sequence of 1 .. m-1 bases).
// Decision: fuzznuc's reports there are an artefact of its matcher (coordinates beyond the sequence end, reverse-strand lines with end < start, IUPAC patterns reporting nothing),
// not a defined search semantics: FindKnownMotifHits reports no hit for a sequence shorter than the motif, whatever the mismatch bound.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_KnownMotifs_ShortSequence_Tests
{
    private sealed record Row(string Sequence, string Pattern, int K, bool Complement, string[] Hits);

    // One parsed copy per fixture (audit A3-15): the oracle was re-read and re-parsed on every call, i.e. in every test.
    private static readonly Lazy<List<Row>> OracleRows = new(ReadOracleRows);

    private static List<Row> ReadOracle() => OracleRows.Value;

    private static List<Row> ReadOracleRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.fuzznuc_short_sequence_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], int.Parse(f[2]), f[3] == "Y", f[4] == "-" ? [] : f[4].Split(';')));
        }

        return rows;
    }

    [Test]
    public void Oracle_ShowsTheArtefacts_OverhangingEndsAndReversedCoordinates()
    {
        var rows = ReadOracle();
        int overhanging = 0, reversed = 0, withHits = 0, iupacWithHits = 0;
        foreach (var r in rows)
        {
            if (r.Hits.Length > 0)
            {
                withHits++;
                if (r.Pattern.Any(c => "RYN".Contains(c)))
                    iupacWithHits++;
            }

            foreach (var h in r.Hits)
            {
                var f = h.Split(',');
                int start = int.Parse(f[0]), end = int.Parse(f[1]);
                if (end > r.Sequence.Length)
                    overhanging++;
                if (end < start)
                    reversed++;
            }
        }

        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.EqualTo(300));
            Assert.That(rows.All(r => r.Sequence.Length < r.Pattern.Length), Is.True, "every sequence is shorter than its pattern");
            Assert.That(withHits, Is.GreaterThan(30), "fuzznuc does report 'hits' there (49 of 300 runs)");
            Assert.That(overhanging, Is.GreaterThan(30), "hit ends beyond the end of the sequence (GFF end > sequence length): 41 hits");
            Assert.That(reversed, Is.GreaterThan(5), "reverse-strand lines with end < start (e.g. 2..-1): 10 hits");
        });
    }

    [Test]
    public void FindKnownMotifHits_SequenceShorterThanMotif_ReportsNoHit_AtAnyMismatchBound()
    {
        foreach (var r in ReadOracle())
        {
            var hits = GenomicAnalyzer.FindKnownMotifHits(new DnaSequence(r.Sequence), [r.Pattern], bothStrands: r.Complement, allowIupac: true, maxMismatches: r.K);
            Assert.That(hits.TryGetValue(r.Pattern, out var found) ? found : [], Is.Empty, $"{r.Sequence} / {r.Pattern} k={r.K} complement={r.Complement}: fuzznuc printed {string.Join(';', r.Hits)}");
        }
    }
}
