using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-STATS-001 (audit A1-9, B09 F42): BLAST-style integer percentages and counts of <see cref="AlignmentStatistics"/>.
/// Definition: NCBI C++ toolkit <c>CAlignFormatUtil::GetPercentMatch</c> (src/objtools/align_format/align_format_util.cpp, opened):
/// <c>(int)(0.5 + 100·n/d)</c>, 100 only if n == d, otherwise at most 99. Oracle: NCBI BLAST+ 2.12.0 blastn / blastp
/// (TestData/Alignment/blast_stats_oracle.tsv, 320 best HSPs, scripts/alignment/generate_blast_stats_oracle.py): tabular
/// <c>length nident mismatch gapopen gaps positive pident</c> and the pairwise-report line
/// <c> Identities = a/L (p%), Positives = b/L (q%), Gaps = g/L (r%)</c>.
/// </summary>
[TestFixture]
public class SequenceAligner_BlastStatistics_Tests
{
    private sealed record Row(string Program, string Q, string S, int Length, int Nident, int Mismatch, int GapOpen, int Gaps,
        int Positive, string Pident, int IdPct, int? PosPct, int GapPct);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.blast_stats_oracle.tsv";
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
            int I(string x) => int.Parse(x, CultureInfo.InvariantCulture);
            rows.Add(new Row(f[0], f[1], f[2], I(f[3]), I(f[4]), I(f[5]), I(f[6]), I(f[7]), I(f[8]), f[9], I(f[10]),
                f[11] == "NA" ? null : I(f[11]), I(f[12])));
        }
        return rows;
    }

    private static AlignmentStatistics Stats(Row r)
    {
        var ar = new AlignmentResult(r.Q, r.S, 0, AlignmentType.Local, 0, 0, 0, 0);
        return r.Program == "blastp"
            ? SequenceAligner.CalculateStatistics(ar, SubstitutionMatrix.Blosum62)
            : SequenceAligner.CalculateStatistics(ar, new ScoringMatrix(1, -2, 0, 0));
    }

    [Test]
    public void Oracle_Counts_And_IntegerPercentages_EqualBlastPlus_AllRows()
    {
        var rows = ReadOracle();
        var failures = new List<string>();
        foreach (var r in rows)
        {
            var s = Stats(r);
            string tag = $"{r.Program} {r.Nident}/{r.Length}";
            void Check(string what, object expected, object actual)
            {
                if (!Equals(expected, actual)) failures.Add($"{tag}: {what} expected {expected} got {actual}");
            }

            Check("length", r.Length, s.AlignmentLength);
            Check("nident", r.Nident, s.Matches);
            Check("mismatch", r.Mismatch, s.Mismatches);
            Check("gapopen", r.GapOpen, s.GapOpens);
            Check("gaps", r.Gaps, s.Gaps);
            Check("positive", r.Positive, s.BlastPositives);
            Check("pident", r.Pident, s.BlastTabularPidentText);
            Check("Identities %", r.IdPct, s.BlastPercentIdentity);
            Check("Gaps %", r.GapPct, s.BlastPercentGaps);
            if (r.PosPct is { } pp) Check("Positives %", pp, s.BlastPercentPositives);

            string expectedLine = r.Program == "blastp"
                ? $" Identities = {r.Nident}/{r.Length} ({r.IdPct}%), Positives = {r.Positive}/{r.Length} ({r.PosPct}%), Gaps = {r.Gaps}/{r.Length} ({r.GapPct}%)"
                : $" Identities = {r.Nident}/{r.Length} ({r.IdPct}%), Gaps = {r.Gaps}/{r.Length} ({r.GapPct}%)";
            Check("report line", expectedLine, s.ToBlastString(r.Program == "blastp"));
        }
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(300));
            Assert.That(rows.Count(r => r.Program == "blastp"), Is.GreaterThanOrEqualTo(150));
            Assert.That(rows.Count(r => r.GapOpen > 0), Is.GreaterThan(100), "gapped HSPs");
            // the "never 100% unless identical" cap is exercised: pident >= 99.5 but not 100 is shown as 99
            Assert.That(rows.Count(r => double.Parse(r.Pident, CultureInfo.InvariantCulture) is >= 99.5 and < 100 && r.IdPct == 99), Is.GreaterThan(10));
            Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
        });
    }

    /// <summary>Definition (align_format_util.cpp): round half up, 100 only for n == d, never above 99 otherwise.</summary>
    [TestCase(0, 10, 0)]
    [TestCase(1, 2, 50)]
    [TestCase(1, 8, 13)]      // 12.5 -> 13 (half up)
    [TestCase(12, 30, 40)]    // BLAST text example "Identities = 12/30 (40%)"
    [TestCase(199, 200, 99)]  // 99.5 would round to 100: capped
    [TestCase(995, 1000, 99)]
    [TestCase(200, 200, 100)]
    [TestCase(0, 0, 0)]       // empty alignment
    public void BlastPercent_FollowsGetPercentMatch(int n, int d, int expected) =>
        Assert.That(AlignmentStatistics.BlastPercent(n, d), Is.EqualTo(expected));

    [Test]
    public void ToBlastString_NucleotideAndProtein_Format()
    {
        var s = new AlignmentStatistics(Matches: 12, Mismatches: 15, Gaps: 3, AlignmentLength: 30, Identity: 40, Similarity: 0,
            GapPercent: 10, Positives: 14, GapOpens: 2, BlastPositives: 14);
        Assert.Multiple(() =>
        {
            Assert.That(s.ToBlastString(), Is.EqualTo(" Identities = 12/30 (40%), Gaps = 3/30 (10%)"));
            Assert.That(s.ToBlastString(true), Is.EqualTo(" Identities = 12/30 (40%), Positives = 14/30 (47%), Gaps = 3/30 (10%)"));
            Assert.That(s.BlastTabularPidentText, Is.EqualTo("40.000"));
            Assert.That(AlignmentStatistics.Empty.BlastTabularPidentText, Is.EqualTo("0.000"));
        });
    }

    /// <summary>BLAST counts an identical pair as an identity even when its own score is not positive (X/X in BLOSUM62), whereas the
    /// EMBOSS similarity count does not: BlastPositives = identities + non-identical positives.</summary>
    [Test]
    public void BlastPositives_CountsIdentityPairsWithNonPositiveSelfScore()
    {
        var ar = new AlignmentResult("AXL", "AXI", 0, AlignmentType.Global, 0, 0, 2, 2);
        var s = SequenceAligner.CalculateStatistics(ar, SubstitutionMatrix.Blosum62);
        Assert.Multiple(() =>
        {
            Assert.That(s.Matches, Is.EqualTo(2));
            Assert.That(s.Positives, Is.EqualTo(2), "A/A and L/I (EMBOSS similarity; X/X has score -1)");
            Assert.That(s.BlastPositives, Is.EqualTo(3), "A/A, X/X (identity) and L/I");
        });
    }
}
