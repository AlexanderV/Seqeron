using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit A2-29b, B09 F125): the ksw2 extension functions <c>ksw_extz</c>, <c>ksw_extz2_sse</c>, <c>ksw_extd</c>, <c>ksw_extd2_sse</c> (<see cref="Ksw2Aligner"/>) against an independent
/// Gotoh program over the band (TestData/Alignment/ksw2_exact_oracle.tsv, scripts/alignment/generate_ksw2_oracle.py): max, mqe, mte and score of 330 cases per variant, and every CIGAR re-scored along its own path.
/// The Z-drop / approximate cases need the real C program (docs/Validation/review-2026-09/handoff/A2-29b-ksw2).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class Ksw2Aligner_Tests
{
    private sealed record Row(int Id, string Variant, Ksw2Option Flag, int W, int EndBonus, int Q, int E, int Q2, int E2, int M, sbyte[] Mat, byte[] Query, byte[] Target, int Max, int Mqe, int Mte, int Score);

    private static int Parse(string s) => s == "NEG" ? Ksw2Aligner.NegInf : int.Parse(s, CultureInfo.InvariantCulture);

    private static byte[] Digits(string s) => s == "-" ? [] : [.. s.Select(c => (byte)(c - '0'))];

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.ksw2_exact_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(int.Parse(f[0], CultureInfo.InvariantCulture), f[1], (Ksw2Option)int.Parse(f[2], CultureInfo.InvariantCulture), int.Parse(f[3], CultureInfo.InvariantCulture),
                int.Parse(f[4], CultureInfo.InvariantCulture), int.Parse(f[5], CultureInfo.InvariantCulture), int.Parse(f[6], CultureInfo.InvariantCulture), int.Parse(f[7], CultureInfo.InvariantCulture),
                int.Parse(f[8], CultureInfo.InvariantCulture), int.Parse(f[9], CultureInfo.InvariantCulture), [.. f[10].Split(',').Select(x => sbyte.Parse(x, CultureInfo.InvariantCulture))],
                Digits(f[11]), Digits(f[12]), Parse(f[13]), Parse(f[14]), Parse(f[15]), Parse(f[16])));
        }

        return list;
    });

    private static Ksw2Result Run(Row r) => r.Variant switch
    {
        "extz" => Ksw2Aligner.ExtendZ(r.Query, r.Target, r.M, r.Mat, r.Q, r.E, r.W, -1, r.Flag),
        "extz2" => Ksw2Aligner.ExtendZ2(r.Query, r.Target, r.M, r.Mat, r.Q, r.E, r.W, -1, r.EndBonus, r.Flag),
        "extd" => Ksw2Aligner.ExtendD(r.Query, r.Target, r.M, r.Mat, r.Q, r.E, r.Q2, r.E2, r.W, -1, r.Flag),
        _ => Ksw2Aligner.ExtendD2(r.Query, r.Target, r.M, r.Mat, r.Q, r.E, r.Q2, r.E2, r.W, -1, r.EndBonus, r.Flag),
    };

    // the residue-pair score the function uses: the SIMD functions without GENERIC use match / mismatch / wildcard (the wildcard is -e, or -e2 for the dual function, when mat[m*m-1] is 0)
    private static int PairScore(Row r, int t, int q)
    {
        bool simd = r.Variant.EndsWith('2');
        if (simd && (r.Flag & Ksw2Option.GenericScore) == 0)
        {
            if (t == r.M - 1 || q == r.M - 1)
            {
                int nn = r.Mat[(r.M * r.M) - 1];
                return nn == 0 ? -(r.Variant == "extd2" ? r.E2 : r.E) : nn;
            }

            return t == q ? r.Mat[0] : r.Mat[1];
        }

        return r.Mat[(t * r.M) + q];
    }

    private static int GapCost(Row r, int k) => r.Variant.StartsWith("extd", StringComparison.Ordinal) ? Math.Min(r.Q + (k * r.E), r.Q2 + (k * r.E2)) : r.Q + (k * r.E);

    // walks the CIGAR from (0, 0): returns (score, target residues, query residues consumed)
    private static (int Score, int T, int Q) Rescore(Row r, IReadOnlyList<(char Op, int Len)> cigar)
    {
        int t = 0;
        int q = 0;
        int score = 0;
        foreach (var (op, len) in cigar)
        {
            if (op == 'M')
            {
                for (int i = 0; i < len; i++)
                    score += PairScore(r, r.Target[t++], r.Query[q++]);
            }
            else
            {
                score -= GapCost(r, len);
                if (op == 'D')
                    t += len;
                else
                    q += len;
            }
        }

        return (score, t, q);
    }

    private static IReadOnlyList<(char Operation, int Length)> Forward(Ksw2Result res, Row r) =>
        (r.Flag & Ksw2Option.ReverseCigar) != 0 ? [.. res.Cigar.Reverse()] : res.Cigar;

    private static void CheckVariant(string variant)
    {
        var failures = new List<string>();
        var rows = Oracle.Value.Where(r => r.Variant == variant).ToList();
        foreach (var r in rows)
        {
            var res = Run(r);
            string head = $"#{r.Id} {variant} flag={(int)r.Flag} w={r.W} q={r.Q} e={r.E} q2={r.Q2} e2={r.E2} eb={r.EndBonus} len={r.Query.Length}x{r.Target.Length}";
            if (res.Max != r.Max || res.Mqe != r.Mqe || res.Mte != r.Mte || res.Score != r.Score)
                failures.Add($"{head}: max/mqe/mte/score got {res.Max}/{res.Mqe}/{res.Mte}/{res.Score} expected {r.Max}/{r.Mqe}/{r.Mte}/{r.Score}");
            if ((r.Flag & Ksw2Option.ScoreOnly) != 0)
            {
                if (res.Cigar.Count != 0)
                    failures.Add($"{head}: score-only run returned a CIGAR");
                continue;
            }

            // the path: to the corner, to the query end (end bonus) or to the best cell
            bool simd = variant.EndsWith('2');
            bool ext = (r.Flag & Ksw2Option.ExtensionOnly) != 0;
            int endT;
            int endQ;
            int expected;
            bool reach = false;
            if (!res.ZDropped && !ext)
            {
                (endT, endQ, expected) = (r.Target.Length - 1, r.Query.Length - 1, res.Score);
            }
            else if (simd && !res.ZDropped && ext && res.Mqe + r.EndBonus > res.Max)
            {
                (endT, endQ, expected, reach) = (res.MqeT, r.Query.Length - 1, res.Mqe, true);
            }
            else if (res.MaxT >= 0 && res.MaxQ >= 0)
            {
                (endT, endQ, expected) = (res.MaxT, res.MaxQ, res.Max);
            }
            else
            {
                if (res.Cigar.Count != 0)
                    failures.Add($"{head}: CIGAR without a path end");
                continue;
            }

            if (res.ReachEnd != reach)
                failures.Add($"{head}: reach_end {res.ReachEnd} expected {reach}");
            var (score, t, q) = Rescore(r, Forward(res, r));
            if (t != endT + 1 || q != endQ + 1)
                failures.Add($"{head}: CIGAR {res.CigarString} consumes {t}/{q}, expected {endT + 1}/{endQ + 1}");
            else if (score != expected)
                failures.Add($"{head}: CIGAR {res.CigarString} scores {score}, expected {expected}");
        }

        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(300));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}\n" + string.Join("\n", failures.Take(25)));
    }

    [Test]
    public void ExtZ_MatchesTheIndependentDpOnAllCases() => CheckVariant("extz");

    [Test]
    public void ExtZ2Sse_MatchesTheIndependentDpOnAllCases() => CheckVariant("extz2");

    [Test]
    public void ExtD_MatchesTheIndependentDpOnAllCases() => CheckVariant("extd");

    [Test]
    public void ExtD2Sse_MatchesTheIndependentDpOnAllCases() => CheckVariant("extd2");

    [Test]
    public void IdenticalSequences_GiveOneMatchRun()
    {
        byte[] s = Ksw2Aligner.EncodeNucleotides("ACGTACGTTGCA");
        var mat = Ksw2Aligner.SimpleMatrix(5, 2, -4, 1);
        var r = Ksw2Aligner.ExtendZ2(s, s, 5, mat, 4, 2, -1, -1, 0);
        Assert.That(r.CigarString, Is.EqualTo("12M"));
        Assert.That(r.Score, Is.EqualTo(24));
        Assert.That(r.Max, Is.EqualTo(24));
    }

    [Test]
    public void AmbiguousBase_ScoresTheWildcardPenalty()
    {
        // one N against A: -ambiguous (minimap2 -N), the rest matches
        var mat = Ksw2Aligner.SimpleMatrix(5, 2, -4, 1);
        var r = Ksw2Aligner.ExtendD2(Ksw2Aligner.EncodeNucleotides("ACGNACGT"), Ksw2Aligner.EncodeNucleotides("ACGTACGT"), 5, mat, 4, 2, 24, 1, -1, -1, 0);
        Assert.That(r.Score, Is.EqualTo((7 * 2) - 1));
        Assert.That(r.CigarString, Is.EqualTo("8M"));
    }

    [Test]
    public void ZDrop_StopsTheExtensionAfterTheDivergence()
    {
        // 30 matching bases followed by unrelated ones: the extension stops (zdropped) with the best cell at the end of the match
        string common = "ACGTTGCAACGTAGCTAGCTAGGATCCATCGA";
        byte[] q = Ksw2Aligner.EncodeNucleotides(common + new string('A', 60));
        byte[] t = Ksw2Aligner.EncodeNucleotides(common + new string('C', 60));
        var mat = Ksw2Aligner.SimpleMatrix(5, 2, -4, 1);
        var r = Ksw2Aligner.ExtendD2(q, t, 5, mat, 4, 2, 24, 1, -1, 10, 0, Ksw2Option.ExtensionOnly);
        Assert.That(r.ZDropped, Is.True);
        Assert.That(r.Max, Is.EqualTo(common.Length * 2));
        Assert.That(r.MaxQ, Is.EqualTo(common.Length - 1));
        Assert.That(r.CigarString, Is.EqualTo(common.Length + "M"));
    }

    [Test]
    public void EndBonus_PrefersReachingTheQueryEnd()
    {
        // the query ends with a mismatch: without the bonus the extension stops before it, with a bonus larger than the loss it is taken to the end
        var mat = Ksw2Aligner.SimpleMatrix(5, 2, -4, 1);
        byte[] q = Ksw2Aligner.EncodeNucleotides("ACGTACGTAC" + "A");
        byte[] t = Ksw2Aligner.EncodeNucleotides("ACGTACGTAC" + "C" + "GGGG");
        var without = Ksw2Aligner.ExtendZ2(q, t, 5, mat, 4, 2, -1, 100, 0, Ksw2Option.ExtensionOnly);
        var with = Ksw2Aligner.ExtendZ2(q, t, 5, mat, 4, 2, -1, 100, 30, Ksw2Option.ExtensionOnly);
        Assert.That(without.ReachEnd, Is.False);
        Assert.That(without.CigarString, Is.EqualTo("10M"));
        Assert.That(with.ReachEnd, Is.True);
        Assert.That(with.CigarString, Is.EqualTo("11M"));
    }
}
