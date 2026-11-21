using System.Globalization;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-LOCAL-001 (audit A6-3, B09 F178): <see cref="SequenceAligner.MatcherProgram"/> on every <c>-gapopen</c> / <c>-gapextend</c> / <c>-alternatives</c> the installed EMBOSS 6.6.0 <c>matcher</c> accepts, in particular <c>-gapopen</c>
/// BELOW <c>-gapextend</c> (<c>matcher a b -gapopen 0 -gapextend 4</c>: score 105 and 4 gaps for the auditor's pair; the first gap residue is cheaper than a further one), 0, 1, 100, and the values the ACD resets (negative: warning
/// "integer value out of range -1 less than (reset to) 0") or refuses ("Invalid integer value"). Oracle: TestData/Alignment/matcher_gap_oracle.tsv (scripts/alignment/generate_matcher_gap_oracle.py): alignment lists (score, coordinates, rows) and the program's warnings.
/// </summary>
[TestFixture]
public class SequenceAligner_MatcherGap_Tests
{
    private sealed record Row(int Line, string Kind, string GapOpen, string GapExtend, string Alternatives, string Seq1, string Seq2, string Status, string Alignments);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.matcher_gap_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        reader.ReadLine();
        int line = 1;
        string? text;
        while ((text = reader.ReadLine()) is not null)
        {
            line++;
            if (text.Length == 0 || text[0] == '#')
                continue;
            var f = text.Split('\t');
            rows.Add(new Row(line, f[0], f[1], f[2], f[3], f[4], f[5], f[6], f.Length > 7 ? f[7] : ""));
        }

        return rows;
    });

    private static readonly Regex Plain = new(@"^[+-]?\d+(\.\d+)?$", RegexOptions.CultureInvariant);

    private static string Key(EmbossPairAlignment a)
    {
        int Residues(string row) => row.Count(c => c != '-');
        return string.Create(CultureInfo.InvariantCulture, $"{(int)Math.Round(a.Score)}|{a.First.Offset}|{a.Second.Offset}|{a.First.Offset + Residues(a.First.Row) - 1}|{a.Second.Offset + Residues(a.Second.Row) - 1}|{a.First.Row}|{a.Second.Row}");
    }

    private static EmbossProgramResult Run(Row r)
    {
        var forced = r.Kind == "dna" ? new EmbossSequenceQualifiers { Nucleotide = true } : new EmbossSequenceQualifiers { Protein = true };
        return SequenceAligner.MatcherProgram(new EmbossSequence("a", r.Seq1, "", forced), new EmbossSequence("b", r.Seq2, "", forced),
            r.Kind == "dna" ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62, double.Parse(r.GapOpen, CultureInfo.InvariantCulture), double.Parse(r.GapExtend, CultureInfo.InvariantCulture),
            int.Parse(r.Alternatives, CultureInfo.InvariantCulture));
    }

    /// <summary>The rows the numeric API can express: the texts are plain decimal numbers ("abc" / "1e2" / "" are not numbers a caller can pass).</summary>
    private static IEnumerable<Row> Rows() => Oracle.Value.Where(r => Plain.IsMatch(r.GapOpen) && Plain.IsMatch(r.GapExtend) && int.TryParse(r.Alternatives, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out _));

    [Test]
    public void MatcherProgram_EqualsEmbossMatcher_OnEveryGapPair_IncludingOpenBelowExtend()
    {
        var failures = new List<string>();
        int ok = 0, died = 0, openBelowExtend = 0, warned = 0;
        foreach (var r in Rows())
        {
            EmbossProgramResult? got = null;
            Exception? error = null;
            try
            {
                got = Run(r);
            }
            catch (Exception e) when (e is ArgumentException or InvalidOperationException)
            {
                error = e;
            }

            if (r.Status == "DIED")
            {
                died++;
                if (error is null)
                    failures.Add($"line {r.Line} {r.Kind} -gapopen {r.GapOpen} -gapextend {r.GapExtend} -alternatives {r.Alternatives}: the program died, the port returned {got!.Alignments.Count} alignments");
                continue;
            }

            if (error is not null)
            {
                failures.Add($"line {r.Line} {r.Kind} -gapopen {r.GapOpen} -gapextend {r.GapExtend} -alternatives {r.Alternatives}: the port refused: {error.Message}");
                continue;
            }

            var parts = r.Alignments.Split('~');
            string expectedList = parts[0];
            string expectedWarnings = parts.Length > 1 ? parts[1] : "";
            string gotList = string.Join(';', got!.Alignments.Select(Key));
            string gotWarnings = string.Join('|', got.StandardError.Split('\n', StringSplitOptions.RemoveEmptyEntries));
            if (gotList != expectedList || gotWarnings != expectedWarnings)
                failures.Add($"line {r.Line} {r.Kind} -gapopen {r.GapOpen} -gapextend {r.GapExtend} -alternatives {r.Alternatives} {r.Seq1} {r.Seq2}\n   expected {expectedList} ~ {expectedWarnings}\n   actual   {gotList} ~ {gotWarnings}");
            else
                ok++;
            if (double.Parse(r.GapOpen, CultureInfo.InvariantCulture) < double.Parse(r.GapExtend, CultureInfo.InvariantCulture))
                openBelowExtend++;
            if (expectedWarnings.Length > 0)
                warned++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count}:\n" + string.Join("\n", failures.Take(4)));
        Assert.That(ok, Is.GreaterThan(600));
        Assert.That(openBelowExtend, Is.GreaterThan(200), "rows with -gapopen below -gapextend");
        Assert.That(died, Is.GreaterThan(20));
        Assert.That(warned, Is.GreaterThan(5), "rows with a reset warning");
    }

    [Test]
    public void MatcherProgram_OpenBelowExtend_TheAuditorsPair_IsAnAlignmentNotAnException()
    {
        // A6-3: matcher a.fa b.fa -gapopen 0 -gapextend 4 was an ArgumentException ("GapOpen must be a finite score <= 0 ... got 4") in the port; a gap of length k costs 0 + (k - 1) * 4
        var first = new EmbossSequence("a", "ACGTACGTAAACGTTTGCACCCGGTACGTA");
        var second = new EmbossSequence("b", "ACGTACGTACGTTTGCACCCGGGTACGTA");
        foreach (int open in new[] { 0, 1, 2, 3, 4 })
            Assert.That(SequenceAligner.MatcherProgram(first, second, null, open, 4).Alignments, Has.Count.EqualTo(1), $"-gapopen {open}");
        Assert.That(SequenceAligner.MatcherProgram(first, second, null, 0, 4).Alignments[0].Score, Is.GreaterThanOrEqualTo(SequenceAligner.MatcherProgram(first, second, null, 4, 4).Alignments[0].Score));
    }

    [Test]
    public void MatcherProgram_IntegerOptions_AreValidatedLikeTheAcd()
    {
        var a = new EmbossSequence("a", "ACGTACGTAAACGTTTGCACCCGGTACGTA");
        var b = new EmbossSequence("b", "ACGTACGTACGTTTGCACCCGGGTACGTA");
        // a non-integer is "Invalid integer value" (the program dies); NaN, infinity and a value beyond 64 bits are not integers either (1e308 was an OverflowException from the decimal conversion)
        foreach (double bad in new[] { 1.5, 0.5, double.NaN, double.PositiveInfinity, double.NegativeInfinity, 1e308, -1e308, 9.3e18 })
        {
            Assert.That(Assert.Throws<ArgumentException>(() => SequenceAligner.MatcherProgram(a, b, null, bad, 4))!.ParamName, Is.EqualTo("gapOpen"), $"{bad}");
            Assert.That(Assert.Throws<ArgumentException>(() => SequenceAligner.MatcherProgram(a, b, null, 16, bad))!.ParamName, Is.EqualTo("gapExtend"), $"{bad}");
        }

        // below the minimum: the program warns and goes on with the minimum
        var reset = SequenceAligner.MatcherProgram(a, b, null, -1, -3, 0);
        var zero = SequenceAligner.MatcherProgram(a, b, null, 0, 0, 1);
        Assert.That(reset.Alignments.Select(Key), Is.EqualTo(zero.Alignments.Select(Key)));
        Assert.That(reset.StandardError, Is.EqualTo("Warning: integer value out of range 0 less than (reset to) 1\nWarning: integer value out of range -1 less than (reset to) 0\nWarning: integer value out of range -3 less than (reset to) 0\n"));
        Assert.That(zero.StandardError, Is.Empty);
        // 100 is legal (no maximum in matcher.acd); a value that would overflow the 32-bit scores of the program is refused by the range check
        Assert.That(SequenceAligner.MatcherProgram(a, b, null, 100, 10).Alignments, Is.Not.Empty);
        Assert.Throws<ArgumentException>(() => SequenceAligner.MatcherProgram(a, b, null, int.MaxValue, 4));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MatcherProgram(a, b, null, 16, 1e12));
    }
}
