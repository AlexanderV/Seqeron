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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: POP-LD-001 — Linkage disequilibrium (PopGen)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 47.
    // Spec: tests/TestSpecs/POP-LD-001.md (canonical CalculateLD).
    // Dimensions: nLoci(3) × nSamples(3) × metric(2: D'/r²). Grid 3×3×2 = 18.
    //
    // Model (Lewontin 1964 D'; Hill & Robertson 1968 r²): from diploid genotypes (0/1/2)
    // CalculateLD reports r² = Cor(X₁,X₂)² and D' = |Cov/2| / D_max, both in [0,1]. They are
    // distinct measures — D' can be 1 (complete LD) while r² < 1 when allele frequencies differ.
    //
    // The combinatorial point: number of loci, sample size and metric choice interact, yet
    // each reported statistic stays in [0,1] and equals its closed-form definition for every
    // locus pair, sample size and metric.
    // ═══════════════════════════════════════════════════════════════════════

    public enum LdMetric { DPrime, RSquared }

    private static double IndepRSquared(int[] x, int[] y)
    {
        int n = x.Length; double mx = x.Average(), my = y.Average();
        double cov = 0, vx = 0, vy = 0;
        for (int i = 0; i < n; i++) { cov += (x[i] - mx) * (y[i] - my); vx += (x[i] - mx) * (x[i] - mx); vy += (y[i] - my) * (y[i] - my); }
        cov /= n; vx /= n; vy /= n;
        double r2 = vx > 0 && vy > 0 ? cov * cov / (vx * vy) : 0;
        return Math.Clamp(r2, 0, 1);
    }

    private static double IndepDPrime(int[] x, int[] y)
    {
        int n = x.Length;
        double p1 = x.Sum() / (2.0 * n), p2 = y.Sum() / (2.0 * n), q1 = 1 - p1, q2 = 1 - p2;
        double mx = x.Average(), my = y.Average();
        double cov = 0; for (int i = 0; i < n; i++) cov += (x[i] - mx) * (y[i] - my);
        cov /= n;
        double d = cov / 2;
        double dMax = d >= 0 ? Math.Min(p1 * q2, q1 * p2) : Math.Min(p1 * p2, q1 * q2);
        double dp = dMax > 1e-10 ? Math.Abs(d) / dMax : 0;
        return Math.Min(dp, 1.0);
    }

    [Test, Combinatorial]
    public void PopLd_StatisticsInRangeAndMatchFormula(
        [Values(3, 4, 5)] int nLoci,
        [Values(12, 30, 60)] int nSamples,
        [Values(LdMetric.DPrime, LdMetric.RSquared)] LdMetric metric)
    {
        // Deterministic 0/1/2 genotypes; locus l = (individual + l) mod 3.
        var loci = Enumerable.Range(0, nLoci)
            .Select(l => Enumerable.Range(0, nSamples).Select(i => (i + l) % 3).ToArray()).ToArray();

        for (int a = 0; a < nLoci; a++)
            for (int b = a + 1; b < nLoci; b++)
            {
                var genos = loci[a].Zip(loci[b], (g1, g2) => (g1, g2));
                var ld = PopulationGeneticsAnalyzer.CalculateLD($"L{a}", $"L{b}", genos, b - a);

                if (metric == LdMetric.RSquared)
                {
                    ld.RSquared.Should().BeInRange(0.0, 1.0);
                    ld.RSquared.Should().BeApproximately(IndepRSquared(loci[a], loci[b]), 1e-9);
                }
                else
                {
                    ld.DPrime.Should().BeInRange(0.0, 1.0);
                    ld.DPrime.Should().BeApproximately(IndepDPrime(loci[a], loci[b]), 1e-9);
                }
                ld.Distance.Should().Be(b - a);
            }
    }

    /// <summary>
    /// Interaction witness: perfectly correlated loci have r² = D' = 1, while independent loci
    /// (zero covariance) have r² = D' = 0.
    /// </summary>
    [Test]
    public void PopLd_PerfectAndZeroLinkage()
    {
        var same = new[] { 0, 1, 2, 1, 0 };
        var perfect = PopulationGeneticsAnalyzer.CalculateLD("a", "b", same.Zip(same, (g1, g2) => (g1, g2)), 1);
        perfect.RSquared.Should().BeApproximately(1.0, 1e-9);
        perfect.DPrime.Should().BeApproximately(1.0, 1e-9);

        int[] x = { 0, 2, 0, 2 }, y = { 0, 0, 2, 2 }; // orthogonal → zero covariance
        var none = PopulationGeneticsAnalyzer.CalculateLD("a", "b", x.Zip(y, (g1, g2) => (g1, g2)), 1);
        none.RSquared.Should().BeApproximately(0.0, 1e-9);
        none.DPrime.Should().BeApproximately(0.0, 1e-9);
    }

    /// <summary>
    /// Interaction witness: D' and r² are distinct metrics — complete LD with unequal allele
    /// frequencies gives D' = 1 but r² = 0.5.
    /// </summary>
    [Test]
    public void PopLd_DPrimeAndRSquared_Differ()
    {
        int[] x = { 0, 1, 2, 1 }, y = { 0, 2, 2, 0 };
        var ld = PopulationGeneticsAnalyzer.CalculateLD("a", "b", x.Zip(y, (g1, g2) => (g1, g2)), 1);

        ld.DPrime.Should().BeApproximately(1.0, 1e-9, "complete LD");
        ld.RSquared.Should().BeApproximately(0.5, 1e-9, "but r² < 1 under unequal frequencies");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: POP-ANCESTRY-001 — Supervised ancestry estimation (PopGen)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 207.
    // Spec: tests/TestSpecs/POP-ANCESTRY-001.md (canonical EstimateAncestry). ADVANCED §10.
    // Dimensions: nPops(3) × nMarkers(3) × method(2). Grid 3×3×2 = 18 (full, exhaustive).
    //
    // Model (Alexander 2009 FRAPPE EM; supervised/projection): with reference allele frequencies fixed,
    // the EM returns ancestry fractions q that sum to 1 and lie in [0,1].
    //
    // Axis mapping (documented): nPops → number of reference populations; nMarkers → number of SNPs;
    // method → {Pure (genotypes match pop0's fixed allele-1 frequencies), Admixed (all-heterozygous)}.
    // Engineered references: pop0 carries allele-1 (freq 1.0), all others allele-2 (freq 0.0). The
    // combinatorial point: a pure individual is assigned ~100% to pop0, an admixed individual spreads
    // its ancestry (pop0 dominant but < 100%), and proportions always form a simplex.
    // ═══════════════════════════════════════════════════════════════════════

    public enum AncestryKind { Pure, Admixed }

    [Test, Combinatorial]
    public void PopAncestry_SimplexAndAssignment_AcrossPopsMarkersMethod(
        [Values(2, 3, 4)] int nPops,
        [Values(5, 10, 20)] int nMarkers,
        [Values(AncestryKind.Pure, AncestryKind.Admixed)] AncestryKind method)
    {
        var refs = Enumerable.Range(0, nPops)
            .Select(p => ($"pop{p}", (IReadOnlyList<double>)Enumerable.Repeat(p == 0 ? 1.0 : 0.0, nMarkers).ToList()))
            .ToList();
        int coreGeno = method == AncestryKind.Pure ? 2 : 1; // all allele-1 homozygous vs all heterozygous
        var individual = ("ind", (IReadOnlyList<int>)Enumerable.Repeat(coreGeno, nMarkers).ToList());

        var result = PopulationGeneticsAnalyzer.EstimateAncestry(new[] { individual }, refs).Single();
        var props = result.Proportions;

        props.Values.Sum().Should().BeApproximately(1.0, 1e-6, "ancestry proportions form a simplex");
        props.Values.Should().OnlyContain(v => v >= -1e-9 && v <= 1.0 + 1e-9, "each proportion is in [0,1]");

        if (method == AncestryKind.Pure)
            props["pop0"].Should().BeGreaterThan(0.99, "a pure pop0 individual is assigned to pop0");
        else
        {
            props["pop0"].Should().Be(props.Values.Max(), "the allele-1 carrier favours pop0");
            props["pop0"].Should().BeLessThan(0.99, "an admixed individual is not pure");
        }
    }

    /// <summary>
    /// Interaction witness — an all-heterozygous individual between an allele-1 and an allele-2 source
    /// population is assigned 50/50.
    /// </summary>
    [Test]
    public void PopAncestry_AllHeterozygous_IsFiftyFifty()
    {
        var refs = new[]
        {
            ("A", (IReadOnlyList<double>)Enumerable.Repeat(1.0, 10).ToList()),
            ("B", (IReadOnlyList<double>)Enumerable.Repeat(0.0, 10).ToList()),
        };
        var ind = ("h", (IReadOnlyList<int>)Enumerable.Repeat(1, 10).ToList());

        var props = PopulationGeneticsAnalyzer.EstimateAncestry(new[] { ind }, refs).Single().Proportions;
        props["A"].Should().BeApproximately(0.5, 1e-9);
        props["B"].Should().BeApproximately(0.5, 1e-9);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: POP-ROH-001 — Runs of homozygosity (PopGen)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 208.
    // Spec: tests/TestSpecs/POP-ROH-001.md (canonical FindROH / CalculateInbreedingFromROH). ADVANCED §10.
    // Dimensions: minLen(3) × density(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Marras 2015; PLINK --homozyg): a run of consecutive homozygous SNPs is reported when it
    // has ≥ minSnps SNPs and spans ≥ minLength bp with no gap > maxGap; F_ROH = ΣL_roh / L_auto.
    //
    // Axis mapping (documented): minLen → the minimum run length (bp); density → the SNP spacing (bp).
    // Engineered construct: 10 consecutive homozygous SNPs at the chosen spacing. The combinatorial
    // point: the run is reported exactly when its span (9·spacing) ≥ minLength, with the correct
    // SNP count, and the ROH inbreeding coefficient equals ΣL/genomeLength.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void PopRoh_RunReportedWhenSpanMeetsMinLength_AcrossMinLenAndDensity(
        [Values(10_000, 50_000, 100_000)] int minLen,
        [Values(2_000, 10_000, 20_000)] int spacing)
    {
        var genotypes = Enumerable.Range(0, 10).Select(i => (i * spacing, 0)).ToList(); // 10 homozygous SNPs
        int span = 9 * spacing;

        var roh = PopulationGeneticsAnalyzer.FindROH(genotypes, minSnps: 5, minLength: minLen,
            maxHeterozygotes: 1, maxGap: 1_000_000).ToList();

        bool expectReported = span >= minLen;
        roh.Any().Should().Be(expectReported, "a run is reported iff its span ≥ minLength");
        if (expectReported)
        {
            var run = roh.Should().ContainSingle().Subject;
            run.Start.Should().Be(0);
            run.End.Should().Be(span);
            run.SnpCount.Should().Be(10);

            double f = PopulationGeneticsAnalyzer.CalculateInbreedingFromROH(
                new[] { (run.Start, run.End) }, genomeLength: 1_000_000);
            f.Should().BeApproximately((double)span / 1_000_000, 1e-9, "F_ROH = ΣL_roh / L_auto");
        }
    }

    /// <summary>
    /// Interaction witness — a heterozygous stretch beyond the tolerance breaks a run, and F_ROH is the
    /// summed ROH length over the genome.
    /// </summary>
    [Test]
    public void PopRoh_HeterozygoteBreaksRun_AndInbreeding()
    {
        // 6 homozygous, then two heterozygotes (exceeds maxHeterozygotes 1), then 6 homozygous.
        var genos = new List<(int, int)>();
        for (int i = 0; i < 6; i++) genos.Add((i * 10_000, 0));
        genos.Add((6 * 10_000, 1));
        genos.Add((7 * 10_000, 1));
        for (int i = 8; i < 14; i++) genos.Add((i * 10_000, 0));

        var roh = PopulationGeneticsAnalyzer.FindROH(genos, minSnps: 5, minLength: 10_000,
            maxHeterozygotes: 1, maxGap: 1_000_000).ToList();
        roh.Should().HaveCount(2, "the heterozygote stretch splits the ROH into two");

        PopulationGeneticsAnalyzer.CalculateInbreedingFromROH(new[] { (0, 100), (200, 350) }, 1000)
            .Should().BeApproximately(0.25, 1e-9, "(100 + 150) / 1000");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: POP-SELECT-001 — Selection statistics (EHH / iHS / scan) (PopGen)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 209.
    // Spec: tests/TestSpecs/POP-SELECT-001.md (canonical CalculateEhh / CalculateIHS / ScanForSelection).
    // ADVANCED §10.
    // Dimensions: statistic(3) × windowSize(3) × nPops(2). Grid 3×3×2 = 18 (full, exhaustive).
    //
    // Model (Sabeti 2002 EHH; Voight 2006 iHS): EHH = Σ C(n_h,2)/C(n_c,2) over haplotype classes;
    // iHS = ln(iHH_A/iHH_D); a selection scan reports the proportion of |iHS| > 2 per window.
    //
    // Axis mapping (documented): statistic → the canonical method (EHH/iHS/Scan); windowSize → the
    // haplotype count / scan window; nPops → the haplotype length / sample factor. The combinatorial
    // point: EHH equals its closed form, iHS is deterministic with a valid derived-allele frequency,
    // and the scan's per-window extreme proportion equals an independent recount.
    // ═══════════════════════════════════════════════════════════════════════

    public enum SelectStat { Ehh, Ihs, Scan }

    private static double Choose2(int n) => n * (n - 1) / 2.0;

    [Test, Combinatorial]
    public void PopSelect_StatisticsAreCorrect_AcrossStatWindowPops(
        [Values(SelectStat.Ehh, SelectStat.Ihs, SelectStat.Scan)] SelectStat statistic,
        [Values(4, 6, 8)] int windowSize,
        [Values(2, 3)] int popsFactor)
    {
        switch (statistic)
        {
            case SelectStat.Ehh:
            {
                int hapLen = popsFactor * 2;
                int half = windowSize / 2;
                var haps = Enumerable.Repeat(new string('A', hapLen), half)
                    .Concat(Enumerable.Repeat(new string('C', hapLen), windowSize - half)).ToList();
                double expected = (Choose2(half) + Choose2(windowSize - half)) / Choose2(windowSize);
                PopulationGeneticsAnalyzer.CalculateEhh(haps).Should().BeApproximately(expected, 1e-9,
                    "EHH = Σ C(n_h,2)/C(n_c,2)");
                break;
            }
            case SelectStat.Ihs:
            {
                int markers = windowSize;
                int nHaps = popsFactor + 2; // ≥ 3 so the core is polymorphic
                var positions = Enumerable.Range(0, markers).Select(i => i * 1000).ToList();
                var haps = Enumerable.Range(0, nHaps)
                    .Select(h => new string(Enumerable.Range(0, markers).Select(m => (char)('0' + (h + m) % 2)).ToArray()))
                    .ToList();
                int core = markers / 2;

                var r1 = PopulationGeneticsAnalyzer.CalculateIHS(haps, positions, core);
                var r2 = PopulationGeneticsAnalyzer.CalculateIHS(haps, positions, core);
                r1.Should().Be(r2, "iHS is deterministic");
                r1.DerivedAlleleFrequency.Should().BeInRange(0.0, 1.0, "a derived-allele frequency is in [0,1]");
                break;
            }
            default:
            {
                // Alternating extreme (|score|=3 > 2) and non-extreme (1) scores.
                var scores = Enumerable.Range(0, windowSize * popsFactor)
                    .Select(i => i % 2 == 0 ? 3.0 : 1.0).ToList();
                var windows = PopulationGeneticsAnalyzer.ScanForSelection(scores, windowSize).ToList();

                foreach (var w in windows)
                {
                    int start = w.WindowIndex * windowSize;
                    int bruteExtreme = Enumerable.Range(start, w.SnpCount).Count(i => Math.Abs(scores[i]) > 2.0);
                    w.ExtremeCount.Should().Be(bruteExtreme, "extreme count = #|score|>2 in the window");
                    w.ProportionExtreme.Should().BeApproximately((double)bruteExtreme / w.SnpCount, 1e-12);
                }
                break;
            }
        }
    }

    /// <summary>
    /// Interaction witness — EHH boundary cases: identical haplotypes give EHH 1, all-distinct give 0.
    /// </summary>
    [Test]
    public void PopSelect_EhhBoundaryCases()
    {
        PopulationGeneticsAnalyzer.CalculateEhh(Enumerable.Repeat("ACGT", 5).ToList())
            .Should().Be(1.0, "all-identical haplotypes are in one class ⇒ EHH 1");
        PopulationGeneticsAnalyzer.CalculateEhh(new[] { "AAAA", "CCCC", "GGGG", "TTTT" })
            .Should().Be(0.0, "all-distinct haplotypes ⇒ no shared pairs ⇒ EHH 0");
    }
}
