using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-STATS-001 (audit R14-7, B09 F300): alignment significance statistics - Karlin-Altschul lambda / K / H, the NCBI gapped parameter tables, bit score, E-value (Spouge finite-size correction for the protein
/// matrices, length-adjusted Karlin-Altschul for blastn), P-value, shuffle Z-score. Oracle: BLAST+ 2.12.0 (<c>blastp -comp_based_stats 0 -seg no -subject</c>, <c>blastn -task blastn -dust no -subject</c>,
/// <c>-outfmt 5</c>): TestData/Alignment/blast_significance_oracle.tsv.gz (scripts/alignment/generate_significance_oracle.py), 2 283 HSPs of seeded related random pairs over 18 protein matrix / gap-cost
/// combinations and 15 blastn reward / penalty / gap-cost combinations; the ungapped parameters equal the NCBI table rows (BLOSUM62 0.3176 / 0.134 / 0.4012) computed from the score distribution.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-STATS-001")]
public class AlignmentSignificance_Tests
{
    private sealed record Row(string Kind, string MatrixOrReward, int Penalty, int Existence, int Extension, int QueryLength, int SubjectLength, int Score, double Bits, double EValue);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        using var stream = TestDataResource.Open("Seqeron.Genomics.Tests.TestData.Alignment.blast_significance_oracle.tsv")
            ?? throw new InvalidOperationException("blast_significance_oracle.tsv not embedded");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var rows = new List<Row>();
        string? line;
        var inv = CultureInfo.InvariantCulture;
        while ((line = reader.ReadLine()) is not null)
        {
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], int.Parse(f[2], inv), int.Parse(f[3], inv), int.Parse(f[4], inv), int.Parse(f[5], inv), int.Parse(f[6], inv), int.Parse(f[7], inv), double.Parse(f[8], inv), double.Parse(f[9], inv)));
        }

        return rows;
    });

    private static void AssertClose(double actual, double expected, double relative, string what)
        => Assert.That(Math.Abs(actual - expected), Is.LessThanOrEqualTo(relative * Math.Abs(expected) + 1e-300), $"{what}: {actual:R} vs {expected:R}");

    [Test]
    public void ProteinHsps_BitScoreAndSpougeEValue_EqualBlastPlus()
    {
        var rows = Oracle.Value.Where(r => r.Kind == "protein").ToList();
        Assert.That(rows.Count, Is.GreaterThan(1000));
        foreach (var r in rows)
        {
            var sig = AlignmentSignificance.SpougeSignificance(r.Score, r.MatrixOrReward, r.Existence, r.Extension, r.QueryLength, r.SubjectLength);
            string id = $"{r.MatrixOrReward} {r.Existence}/{r.Extension} score {r.Score} {r.QueryLength}x{r.SubjectLength}";
            Assert.That(sig.BitScore, Is.EqualTo(r.Bits).Within(0.0006), id);
            AssertClose(sig.EValue, r.EValue, 2e-5, "E-value " + id);
            Assert.That(sig.Method, Is.EqualTo("SpougeFiniteSizeCorrection"));
        }
    }

    [Test]
    public void NucleotideHsps_BitScoreAndLengthAdjustedEValue_EqualBlastPlus()
    {
        var rows = Oracle.Value.Where(r => r.Kind == "nucleotide").ToList();
        Assert.That(rows.Count, Is.GreaterThan(800));
        foreach (var r in rows)
        {
            int reward = int.Parse(r.MatrixOrReward, CultureInfo.InvariantCulture);
            var sig = AlignmentSignificance.ForNucleotide(r.Score, new ScoringMatrix(reward, r.Penalty, -r.Existence, -r.Extension), r.QueryLength, r.SubjectLength);
            string id = $"{reward}/{r.Penalty} {r.Existence}/{r.Extension} score {r.Score} {r.QueryLength}x{r.SubjectLength}";
            Assert.That(sig.BitScore, Is.EqualTo(r.Bits).Within(0.0006), id);
            AssertClose(sig.EValue, r.EValue, 2e-5, "E-value " + id);
            Assert.That(sig.Method, Is.EqualTo("LengthAdjustedKarlinAltschul"));
        }
    }

    [TestCase("BLOSUM62", 0.3176, 0.134, 0.4012)]
    [TestCase("BLOSUM45", 0.2291, 0.0924, 0.2514)]
    [TestCase("BLOSUM50", 0.2318, 0.112, 0.3362)]
    [TestCase("BLOSUM90", 0.3346, 0.19, 0.7547)]
    [TestCase("PAM30", 0.3400, 0.283, 1.7536)]
    [TestCase("PAM70", 0.3345, 0.229, 1.0291)]
    [TestCase("PAM250", 0.2252, 0.0870, 0.2223)]
    public void UngappedParameters_EqualTheNcbiRowOfTheMatrix(string name, double lambda, double k, double h)
    {
        // blast_stat.c: the first row of each table (the "ideal" block: Robinson & Robinson background); K of the table is printed to 3 digits
        var p = AlignmentSignificance.Ungapped(SubstitutionMatrix.Get(name));
        Assert.That(p.Lambda, Is.EqualTo(lambda).Within(6e-4), name);
        Assert.That(p.K, Is.EqualTo(k).Within(0.002), name);
        Assert.That(p.H, Is.EqualTo(h).Within(0.0006), name);
    }

    [Test]
    public void Ungapped_BlastnSchemes_AreTheNcbiUngappedRows()
    {
        // blastn_values_*: the first rows with existence 0 / extension 0 of 1/-5, 1/-4 and 1/-3 are the ungapped values (1/-2 and 2/-3 list the linear-gap values there, which differ from the ungapped ones: 1/-2 ungapped is the root of e^l / 4 + 3 e^(-2 l) / 4 = 1)
        foreach (var (reward, penalty, lambda, k, h) in new[] { (1, -5, 1.39, 0.747, 1.38), (1, -4, 1.383, 0.738, 1.36), (1, -3, 1.374, 0.711, 1.31) })
        {
            var p = AlignmentSignificance.UngappedNucleotide(reward, penalty);
            Assert.That(p.Lambda, Is.EqualTo(lambda).Within(0.006), $"{reward}/{penalty}");
            Assert.That(p.K, Is.EqualTo(k).Within(0.01), $"{reward}/{penalty}");
            Assert.That(p.H, Is.EqualTo(h).Within(0.01), $"{reward}/{penalty}");
        }

        Assert.That(AlignmentSignificance.UngappedNucleotide(1, -2).Lambda, Is.EqualTo(1.3327057628).Within(1e-9), "root of exp(l) / 4 + 3 exp(-2 l) / 4 = 1");
        // the closed forms: match 1 / mismatch -1 on a uniform alphabet is the +-1 walk: lambda = ln 3, K = 1/3
        var pm = AlignmentSignificance.UngappedNucleotide(1, -1);
        Assert.That(pm.Lambda, Is.EqualTo(Math.Log(3)).Within(1e-9));
        Assert.That(pm.K, Is.EqualTo(1.0 / 3).Within(1e-9));
    }

    [Test]
    public void Lambda_SolvesTheDefiningEquation_ForAnyDistribution()
    {
        // Karlin & Altschul: sum p_s exp(lambda s) = 1 (the worked example of NCBI's comment: scores -2 / 0 / 3 with probabilities 0.7 / 0.1 / 0.2 give lambda 0.330, K 0.154)
        var p = AlignmentSignificance.FromScoreDistribution([(-2, 0.7), (0, 0.1), (3, 0.2)]);
        Assert.That(p.Lambda, Is.EqualTo(0.330).Within(0.0005));
        Assert.That(p.K, Is.EqualTo(0.154).Within(0.0005));
        double residual = (0.7 * Math.Exp(-2 * p.Lambda)) + 0.1 + (0.2 * Math.Exp(3 * p.Lambda)) - 1;
        Assert.That(residual, Is.EqualTo(0).Within(1e-6));
        Assert.Throws<ArgumentException>(() => AlignmentSignificance.FromScoreDistribution([(1, 0.5), (-1, 0.5)]), "expected score 0");
        Assert.Throws<ArgumentException>(() => AlignmentSignificance.FromScoreDistribution([(-1, 0.5), (-2, 0.5)]), "no positive score");
    }

    [Test]
    public void TabulatedGappedParameters_BLOSUM62_11_1_AreTheDocumentedOnes()
    {
        Assert.That(AlignmentSignificance.TryGappedProtein("BLOSUM62", 11, 1, out var p), Is.True);
        Assert.That((p.Lambda, p.K, p.H, p.Alpha, p.Beta), Is.EqualTo((0.267, 0.041, 0.14, 1.9, -30.0)));
        Assert.That(AlignmentSignificance.TryGappedProtein("EBLOSUM62", 11, 1, out _), Is.True, "an EMBOSS name is accepted");
        Assert.That(AlignmentSignificance.TryGappedProtein("BLOSUM62", 7, 1, out _), Is.False, "not tabulated");
        Assert.That(AlignmentSignificance.TryGappedProtein("GONNET", 11, 1, out _), Is.False);
        Assert.That(AlignmentSignificance.TryGappedNucleotide(2, -3, 5, 2, out var n, out bool down), Is.True);
        Assert.That((n.Lambda, n.K, n.H, down), Is.EqualTo((0.625, 0.41, 0.78, true)));
        Assert.That(AlignmentSignificance.TryGappedNucleotide(4, -6, 10, 4, out var scaled, out _), Is.True, "reward / penalty are divided by their gcd (4/-6 = 2/-3, costs 10/4 = 5/2): lambda is halved");
        Assert.That(scaled.Lambda, Is.EqualTo(0.625 / 2).Within(1e-12));
        Assert.That(AlignmentSignificance.TryGappedNucleotide(7, -3, 5, 2, out _, out _), Is.False);
    }

    [Test]
    public void LengthAdjustment_And_PValue_FollowTheNcbiDefinitions()
    {
        var p = new KarlinAltschulParameters(0.625, 0.41, 0.78, 0.8, -2);
        Assert.That(AlignmentSignificance.LengthAdjustment(p, 430, 880), Is.EqualTo(13), "blastn -task blastn, the pair of the oracle");
        Assert.That(KarlinAltschulParameters.PValue(1e-10), Is.EqualTo(1e-10 - 5e-21).Within(1e-30), "-expm1(-E) = E - E^2 / 2 for a small E (no cancellation)");
        Assert.That(KarlinAltschulParameters.PValue(1), Is.EqualTo(1 - Math.Exp(-1)).Within(1e-15));
        Assert.That(KarlinAltschulParameters.PValue(1.24e-33), Is.EqualTo(1.24e-33).Within(1e-45), "a 1e-33 E-value keeps its digits (Math.Exp(x) - 1 would give 0)");
        Assert.That(p.BitScore(100), Is.EqualTo(((0.625 * 100) - Math.Log(0.41)) / Math.Log(2)).Within(1e-12));
        Assert.That(AlignmentSignificance.LengthAdjustment(p, 2, 2), Is.EqualTo(0), "K m n below max(m, n): no room for an adjustment (NCBI returns 0)");
    }

    [Test]
    public void ForMatrix_UsesTheNcbiMatrixOnly_TheEmbossBlosum80HasOtherUnits()
    {
        var s = AlignmentSignificance.ForMatrix(256, new MatrixScoring(SubstitutionMatrix.Get("BLOSUM62"), -11, -1), 150, 160);
        Assert.That(s.EValue, Is.EqualTo(1.24008e-33).Within(1e-38));
        Assert.That(s.BitScore, Is.EqualTo(103.219).Within(0.001));
        // the EMBOSS EBLOSUM80 has A/A = 7, the NCBI BLOSUM80 5: the NCBI table (lambda 0.343) is not its table
        Assert.Throws<ArgumentException>(() => AlignmentSignificance.ForMatrix(100, new MatrixScoring(SubstitutionMatrix.Get("BLOSUM80"), -10, -1), 150, 160));
        Assert.Throws<ArgumentException>(() => AlignmentSignificance.ForMatrix(100, new MatrixScoring(SubstitutionMatrix.Get("BLOSUM62"), -7, -1), 150, 160), "gap costs not in the table");
    }

    [Test]
    public void ShuffleZScore_IsReproducible_AndLargeForARealMatch()
    {
        string core = "ACGTTGCAAGGCTTAACCGGTTAACGTACGGTACCATGCA";
        string a = "TTTT" + core + "GGGG", b = "CCCCAAAA" + core + "TTTTCCCC";
        var z1 = AlignmentSignificance.ShuffleZScore(a, b, shuffles: 200, seed: 7);
        var z2 = AlignmentSignificance.ShuffleZScore(a, b, shuffles: 200, seed: 7);
        Assert.That(z1, Is.EqualTo(z2));
        Assert.That(z1.Score, Is.GreaterThanOrEqualTo(core.Length));
        Assert.That(z1.ZScore, Is.GreaterThan(5), "the shared 40-mer is far above the shuffled scores");
        Assert.That(z1.StandardDeviation, Is.GreaterThan(0));
        var unrelated = AlignmentSignificance.ShuffleZScore("ACGTACGTAC", "ACGTACGTAC", shuffles: 50, seed: 3);
        Assert.That(unrelated.ZScore, Is.GreaterThan(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentSignificance.ShuffleZScore(a, b, shuffles: 1));
    }
}
