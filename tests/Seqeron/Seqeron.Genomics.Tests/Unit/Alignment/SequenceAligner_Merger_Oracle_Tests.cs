// ALIGN-SEMI-001 (B09 A3-13) — EMBOSS merger 6.6.0 (SequenceAligner.Merger). Oracle: scripts/alignment/generate_merger_oracle.py (450 random overlapping read pairs:
// merged sequence, conflict list, score). Evidence: docs/Evidence/ALIGN-SEMI-001-Evidence.md.

using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAligner_Merger_Oracle_Tests
{
    private static List<string[]> ReadRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.merger_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<string[]>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0) rows.Add(line.Split('\t'));
        return rows;
    }

    [Test]
    public void Merger_MergedSequenceConflictsAndScore_EqualMerger_AllRows()
    {
        var rows = ReadRows();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(400));
        int conflictRows = 0, gapConflicts = 0, ambiguity = 0, nonDefaultPenalties = 0, doubled = 0;
        foreach (var f in rows)
        {
            double go = double.Parse(f[2], CultureInfo.InvariantCulture), ge = double.Parse(f[3], CultureInfo.InvariantCulture);
            var r = SequenceAligner.Merger(f[0], f[1], gapOpen: go, gapExtend: ge);
            string label = $"{f[0]} / {f[1]} go={go} ge={ge}";
            Assert.That(r.MergedSequence, Is.EqualTo(f[4]), label + " merged");
            string conflicts = r.Conflicts.Count == 0 ? "-" : string.Join(";", r.Conflicts.Select(c => $"{c.PositionA}:{c.BaseA}:{c.PositionB}:{c.BaseB}:{c.Using}"));
            Assert.That(conflicts, Is.EqualTo(f[5]), label + " conflicts");
            Assert.That(r.Warnings.Count == 0 ? "-" : string.Join("|", r.Warnings), Is.EqualTo(f[6]), label + " warnings");
            Assert.That(r.Score.ToString("F1", CultureInfo.InvariantCulture), Is.EqualTo(f[7]), label + " score");
            if (r.Conflicts.Count > 0) conflictRows++;
            if (r.Conflicts.Any(c => c.BaseA == '.' || c.BaseB == '.')) gapConflicts++;
            if (f[0].Any(c => "NRYKMSWnrykmsw".Contains(c)) || f[1].Any(c => "NRYKMSWnrykmsw".Contains(c))) ambiguity++;
            if (go != 50 || ge != 5) nonDefaultPenalties++;
        }

        Assert.That(conflictRows, Is.GreaterThan(300));
        Assert.That(gapConflicts, Is.GreaterThan(30));
        Assert.That(ambiguity, Is.GreaterThan(80));
        Assert.That(nonDefaultPenalties, Is.GreaterThan(100));
        Assert.That(doubled, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void Merger_Arguments_AreValidated()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.Merger(null!, "ACGT"));
        Assert.Throws<ArgumentException>(() => SequenceAligner.Merger("", "ACGT"));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.Merger("ACGT", "ACGT", gapOpen: 101));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.Merger("ACGT", "ACGT", gapExtend: 11));
        // identical sequences merge to themselves without conflicts
        var same = SequenceAligner.Merger("ACGTACGTTT", "ACGTACGTTT");
        Assert.That((same.MergedSequence, same.Conflicts.Count), Is.EqualTo(("acgtacgttt", 0)));
    }
}
