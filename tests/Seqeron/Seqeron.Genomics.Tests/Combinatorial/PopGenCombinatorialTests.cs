namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the PopGen area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("PopGen")]
public class PopGenCombinatorialTests
{
    private static string DiverseDna(int n, uint seed)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    private static char NextBase(char c) => "ACGT"[("ACGT".IndexOf(c) + 1) % 4];

    /// <summary>A deterministic haploid sample: each sequence carries a few substitutions of a base.</summary>
    private static List<char[]> MakePopulation(int nSeqs, int seqLen, uint seed = 0xD1Eu)
    {
        string baseSeq = DiverseDna(seqLen, seed);
        var pop = new List<char[]>();
        for (int i = 0; i < nSeqs; i++)
        {
            var chars = baseSeq.ToCharArray();
            for (int k = 0; k <= i; k++) { int p = (k * 7 + i * 3) % seqLen; chars[p] = NextBase(chars[p]); }
            pop.Add(chars);
        }
        return pop;
    }

    // ── Independent re-derivations of the population-genetics estimators (ground truth) ──
    private static int SegregatingSites(List<char[]> pop)
    {
        int len = pop[0].Length, s = 0;
        for (int pos = 0; pos < len; pos++)
        {
            char first = pop[0][pos];
            if (pop.Any(x => x[pos] != first)) s++;
        }
        return s;
    }

    private static double Pi(List<char[]> pop)
    {
        int n = pop.Count, len = pop[0].Length;
        double totalDiff = 0; int comparisons = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                int d = 0;
                for (int k = 0; k < len; k++) if (pop[i][k] != pop[j][k]) d++;
                totalDiff += d; comparisons++;
            }
        return totalDiff / (comparisons * len);
    }

    private static double Theta(int s, int n, int len)
    {
        if (n < 2 || len <= 0) return 0;
        double a1 = 0; for (int i = 1; i < n; i++) a1 += 1.0 / i;
        return (double)s / (a1 * len);
    }

    private static double TajimaD(double kHat, int s, int n)
    {
        if (s == 0 || n < 3) return 0;
        double a1 = 0, a2 = 0; for (int i = 1; i < n; i++) { a1 += 1.0 / i; a2 += 1.0 / (i * i); }
        double watterson = (double)s / a1;
        double b1 = (n + 1.0) / (3 * (n - 1));
        double b2 = 2.0 * (n * n + n + 3) / (9.0 * n * (n - 1));
        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1);
        double e1 = c1 / a1, e2 = c2 / (a1 * a1 + a2);
        double variance = e1 * s + e2 * s * (s - 1);
        return variance <= 0 ? 0 : (kHat - watterson) / Math.Sqrt(variance);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: POP-DIV-001 — Nucleotide diversity estimators (PopGen)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 44.
    // Spec: tests/TestSpecs/POP-DIV-001.md (canonical CalculateDiversityStatistics et al.).
    // Dimensions: nSeqs(3) × seqLen(3) × method(3: π/θ/D). Grid 3×3×3 = 27.
    //
    // Model: π (Nei & Li 1979) is the mean per-site pairwise difference; Watterson's θ
    // (1975) is S/(a₁·L) with a₁=Σ1/i; Tajima's D (1989) is (k̂ − S/a₁)/√V̂. The combined
    // CalculateDiversityStatistics must report each estimator equal to its standalone routine
    // and to an independent re-derivation of the published formula.
    //
    // The combinatorial point: sample size, length and estimator choice interact (a₁ depends
    // on n, normalisation on L, D on all of S, k̂, n); every cell is checked against the
    // closed-form definition so no estimator silently diverges.
    // ═══════════════════════════════════════════════════════════════════════

    public enum DiversityEstimator { Pi, Theta, TajimaD }

    [Test, Combinatorial]
    public void PopDiv_Estimators_MatchClosedForm_AcrossSampleAndLength(
        [Values(3, 4, 5)] int nSeqs,
        [Values(20, 40, 80)] int seqLen,
        [Values(DiversityEstimator.Pi, DiversityEstimator.Theta, DiversityEstimator.TajimaD)] DiversityEstimator estimator)
    {
        var pop = MakePopulation(nSeqs, seqLen);
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(pop);
        int s = SegregatingSites(pop);
        stats.SegregratingSites.Should().Be(s, "segregating-site count is independently reproduced");

        switch (estimator)
        {
            case DiversityEstimator.Pi:
                stats.NucleotideDiversity.Should().BeGreaterThanOrEqualTo(0);
                stats.NucleotideDiversity.Should().BeApproximately(Pi(pop), 1e-12);
                stats.NucleotideDiversity.Should().Be(PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(pop));
                break;
            case DiversityEstimator.Theta:
                stats.WattersonTheta.Should().BeGreaterThanOrEqualTo(0);
                stats.WattersonTheta.Should().BeApproximately(Theta(s, nSeqs, seqLen), 1e-12);
                stats.WattersonTheta.Should().Be(
                    PopulationGeneticsAnalyzer.CalculateWattersonTheta(s, nSeqs, seqLen));
                break;
            default:
                double kHat = Pi(pop) * seqLen;
                stats.TajimasD.Should().BeApproximately(TajimaD(kHat, s, nSeqs), 1e-9);
                stats.TajimasD.Should().Be(
                    PopulationGeneticsAnalyzer.CalculateTajimasD(stats.NucleotideDiversity * seqLen, s, nSeqs));
                break;
        }
    }

    /// <summary>
    /// Interaction witness: a monomorphic sample has zero segregating sites, so π, θ and
    /// Tajima's D are all zero regardless of sample size or length.
    /// </summary>
    [Test]
    public void PopDiv_MonomorphicSample_AllZero()
    {
        var pop = Enumerable.Range(0, 5).Select(_ => new string('A', 30).ToCharArray()).ToList();
        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(pop);

        stats.SegregratingSites.Should().Be(0);
        stats.NucleotideDiversity.Should().Be(0);
        stats.WattersonTheta.Should().Be(0);
        stats.TajimasD.Should().Be(0);
    }

    /// <summary>
    /// Interaction witness: Tajima's D sign tracks the allele-frequency spectrum — a single
    /// intermediate-frequency site (excess π) gives D &gt; 0, whereas a singleton (excess rare
    /// variants) gives D &lt; 0. Both share S = 1 and n = 4.
    /// </summary>
    [Test]
    public void PopDiv_TajimaD_SignTracksFrequencySpectrum()
    {
        // One site, balanced 2:2 → intermediate frequency → D > 0.
        var balanced = new List<char[]>
        {
            "ACGT".ToCharArray(), "ACGT".ToCharArray(), "TCGT".ToCharArray(), "TCGT".ToCharArray(),
        };
        PopulationGeneticsAnalyzer.CalculateDiversityStatistics(balanced).TajimasD.Should().BePositive();

        // One site, singleton 1:3 → excess rare variants → D < 0.
        var singleton = new List<char[]>
        {
            "TCGT".ToCharArray(), "ACGT".ToCharArray(), "ACGT".ToCharArray(), "ACGT".ToCharArray(),
        };
        PopulationGeneticsAnalyzer.CalculateDiversityStatistics(singleton).TajimasD.Should().BeNegative();
    }

    /// <summary>
    /// Worked example: π is the mean per-site pairwise difference — for {AAAA, AAAA, AATA}
    /// that is 2 differences over 3 pairs × 4 sites = 1/6.
    /// </summary>
    [Test]
    public void PopDiv_NucleotideDiversity_WorkedExample()
    {
        var pop = new List<char[]> { "AAAA".ToCharArray(), "AAAA".ToCharArray(), "AATA".ToCharArray() };
        PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(pop).Should().BeApproximately(2.0 / 12.0, 1e-12);
    }
}
