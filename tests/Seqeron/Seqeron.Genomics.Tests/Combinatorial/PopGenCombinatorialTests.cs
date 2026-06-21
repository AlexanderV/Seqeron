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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: POP-FST-001 — Population differentiation Fst (PopGen)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 46.
    // Spec: tests/TestSpecs/POP-FST-001.md (canonical CalculateFst / CalculateFStatistics).
    // Dimensions: nPops(3) × nSamples(3) × method(2: Wright/WC). Grid 3×3×2 = 18.
    //
    // Model (Wright 1965; Nei Gst): Fst measures the among-population share of genetic
    // variance, in [0,1]. The implementation offers two estimators — Wright's variance form
    // σ²_S/p̄(1−p̄) (CalculateFst) and the heterozygosity form 1−H_S/H_T (CalculateFStatistics).
    // The "WC" axis maps to the heterozygosity-based estimator (the other implemented form).
    //
    // The combinatorial point: number of populations, sample size and estimator interact, yet
    // both estimators must lie in [0,1], give 0 for identical populations, 1 for a fixed
    // difference, and match their closed-form definitions; the Wright pairwise matrix over
    // nPops must be symmetric with a zero diagonal.
    // ═══════════════════════════════════════════════════════════════════════

    public enum FstEstimator { Wright, Heterozygosity }

    private static readonly double[] BaseFreqs = { 0.5, 0.2, 0.8, 0.4, 0.6 };

    private static (double Freq, int N)[] PopProfile(int popIndex, int nSamples) =>
        BaseFreqs.Select(b => (Math.Clamp(b + popIndex * 0.12, 0.0, 1.0), nSamples)).ToArray();

    private static double IndepWrightFst((double Freq, int N)[] a, (double Freq, int N)[] b)
    {
        double num = 0, den = 0;
        for (int i = 0; i < a.Length; i++)
        {
            double p1 = a[i].Freq, p2 = b[i].Freq; int n1 = a[i].N, n2 = b[i].N;
            double pBar = (n1 * p1 + n2 * p2) / (n1 + n2);
            num += ((p1 - pBar) * (p1 - pBar) * n1 + (p2 - pBar) * (p2 - pBar) * n2) / (n1 + n2);
            den += pBar * (1 - pBar);
        }
        return den > 0 ? num / den : 0;
    }

    private static IEnumerable<(int, int, int, int, double, double)> HetData((double Freq, int N)[] a, (double Freq, int N)[] b)
    {
        for (int i = 0; i < a.Length; i++)
        {
            int h1 = (int)Math.Round(2 * a[i].Freq * (1 - a[i].Freq) * a[i].N);
            int h2 = (int)Math.Round(2 * b[i].Freq * (1 - b[i].Freq) * b[i].N);
            yield return (h1, a[i].N, h2, b[i].N, a[i].Freq, b[i].Freq);
        }
    }

    private static double IndepHetFst((double Freq, int N)[] a, (double Freq, int N)[] b)
    {
        double obs = 0, exp = 0, tot = 0; int totalN = 0;
        foreach (var (h1, n1, h2, n2, p1, p2) in HetData(a, b))
        {
            double pBar = (n1 * p1 + n2 * p2) / (n1 + n2);
            obs += h1 + h2;
            exp += 2 * p1 * (1 - p1) * n1 + 2 * p2 * (1 - p2) * n2;
            tot += 2 * pBar * (1 - pBar) * (n1 + n2);
            totalN += n1 + n2;
        }
        double hs = totalN > 0 ? exp / totalN : 0, ht = totalN > 0 ? tot / totalN : 0;
        return ht > 0 ? 1 - hs / ht : 0;
    }

    private static double Fst(FstEstimator method, (double Freq, int N)[] a, (double Freq, int N)[] b) =>
        method == FstEstimator.Wright
            ? PopulationGeneticsAnalyzer.CalculateFst(a, b)
            : PopulationGeneticsAnalyzer.CalculateFStatistics("a", "b", HetData(a, b)).Fst;

    [Test, Combinatorial]
    public void PopFst_InRangeMatchesFormula_AndMatrixWellFormed(
        [Values(2, 3, 4)] int nPops,
        [Values(10, 30, 50)] int nSamples,
        [Values(FstEstimator.Wright, FstEstimator.Heterozygosity)] FstEstimator method)
    {
        var pops = Enumerable.Range(0, nPops).Select(p => PopProfile(p, nSamples)).ToArray();

        // The most-diverged pair (0, nPops-1) differs, so Fst is positive but bounded.
        double fst = Fst(method, pops[0], pops[nPops - 1]);
        fst.Should().BeInRange(0.0, 1.0 + 1e-9, "Fst is a variance proportion");
        double expected = method == FstEstimator.Wright
            ? IndepWrightFst(pops[0], pops[nPops - 1])
            : IndepHetFst(pops[0], pops[nPops - 1]);
        fst.Should().BeApproximately(expected, 1e-9, "Fst equals its closed-form definition");

        // A population against itself shows no differentiation.
        Fst(method, pops[0], pops[0]).Should().BeApproximately(0.0, 1e-12, "no structure within one population");

        if (method == FstEstimator.Wright)
        {
            var matrix = PopulationGeneticsAnalyzer.CalculatePairwiseFst(
                pops.Select((p, i) => ($"P{i}", (IReadOnlyList<(double, int)>)p)));
            for (int i = 0; i < nPops; i++)
            {
                matrix[i, i].Should().Be(0, "self Fst is zero");
                for (int j = i + 1; j < nPops; j++)
                    matrix[i, j].Should().Be(matrix[j, i], "the Fst matrix is symmetric");
            }
        }
    }

    /// <summary>
    /// Interaction witness: a fixed difference (one population fixed for the allele, the other
    /// for the alternative) is maximal differentiation — Fst = 1 under both estimators.
    /// </summary>
    [Test]
    public void PopFst_FixedDifference_IsOne()
    {
        var fixedA = new[] { (1.0, 40), (1.0, 40) };
        var fixedB = new[] { (0.0, 40), (0.0, 40) };

        PopulationGeneticsAnalyzer.CalculateFst(fixedA, fixedB).Should().BeApproximately(1.0, 1e-12);
        PopulationGeneticsAnalyzer.CalculateFStatistics("a", "b", HetData(fixedA, fixedB)).Fst
            .Should().BeApproximately(1.0, 1e-12);
    }

    /// <summary>
    /// Interaction witness: Fst increases with divergence — populations farther apart in
    /// allele frequency are more differentiated (Wright estimator).
    /// </summary>
    [Test]
    public void PopFst_IncreasesWithDivergence()
    {
        var p0 = PopProfile(0, 40);
        double near = PopulationGeneticsAnalyzer.CalculateFst(p0, PopProfile(1, 40));
        double far = PopulationGeneticsAnalyzer.CalculateFst(p0, PopProfile(3, 40));
        far.Should().BeGreaterThan(near, "greater allele-frequency divergence raises Fst");
    }
}
