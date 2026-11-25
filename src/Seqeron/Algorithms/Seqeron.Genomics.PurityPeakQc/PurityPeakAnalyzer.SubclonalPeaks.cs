// SPDX-License-Identifier: GPL-3.0-or-later
// Port of CNAqc 1.1.5 (caravagnalab/CNAqc, "License: GPL (>= 3) + file LICENSE"), BMix 1.0.0 (caravagnalab/BMix, "License: GPL (>= 3)") and the R stats kmeans (Hartigan-Wong, MacQueen).
// Part of Seqeron.Genomics.PurityPeakQc, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.PurityPeakQc;

/// <summary>CNAqc subclonal peak expectations (port of CNAqc 1.1.5).</summary>
public static partial class PurityPeakAnalyzer
{
    /// <summary>
    /// A subclonal copy-number segment for <see cref="AnalyzeSubclonalPurityPeaks"/>: two subclones with karyotypes
    /// Major:minor (cell fraction <see cref="Ccf"/>) and SecondMajor:SecondMinor (1 − Ccf), and the VAFs of the mutations
    /// mapped to it (CNAqc <c>x$cna_subclonal</c> row: <c>karyotype</c>, <c>karyotype_2</c>, <c>CCF</c>, <c>mutations</c>).
    /// Both karyotypes must be simple (1:0, 1:1, 2:0, 2:1, 2:2) — CNAqc drops other subclonal segments at <c>init</c>.
    /// </summary>
    /// <param name="Chromosome">Chromosome (identification only).</param>
    /// <param name="Start">Segment start (identification only).</param>
    /// <param name="End">Segment end (identification only).</param>
    /// <param name="MajorCopyNumber">First subclone Major.</param>
    /// <param name="MinorCopyNumber">First subclone minor.</param>
    /// <param name="SecondMajorCopyNumber">Second subclone Major.</param>
    /// <param name="SecondMinorCopyNumber">Second subclone minor.</param>
    /// <param name="Ccf">Cell fraction of the first subclone, in (0, 1).</param>
    /// <param name="Vafs">VAFs of the segment's mutations (count = CNAqc <c>n</c>).</param>
    public sealed record SubclonalPeakSegment(
        string Chromosome,
        long Start,
        long End,
        int MajorCopyNumber,
        int MinorCopyNumber,
        int SecondMajorCopyNumber,
        int SecondMinorCopyNumber,
        double Ccf,
        IReadOnlyList<double> Vafs);

    /// <summary>Evolutionary model of a subclonal segment (CNAqc <c>model</c>).</summary>
    public enum SubclonalEvolutionModel
    {
        /// <summary>Both subclones descend from the starting state independently (<c>model_id</c> "start -&gt; g1 | g2").</summary>
        Branching,

        /// <summary>One subclone descends from the other (<c>model_id</c> "start -&gt; g1 -&gt; g2").</summary>
        Linear,
    }

    /// <summary>Parameters of <see cref="AnalyzeSubclonalPurityPeaks"/>; defaults are what CNAqc 1.1.5 <c>analyze_peaks</c> passes.</summary>
    public sealed record SubclonalPeakOptions
    {
        /// <summary>Default options.</summary>
        public static SubclonalPeakOptions Default { get; } = new();

        /// <summary>Segments with n &gt; this (strict) are analysed (<c>n_min = min_absolute_karyotype_mutations</c>, 100).</summary>
        public int MinMutations { get; init; } = 100;

        /// <summary>Peak-matching tolerance |x − peak| ≤ ε (<c>epsilon = purity_error</c>, 0.05), in (0, 1).</summary>
        public double Epsilon { get; init; } = 0.05;

        /// <summary>KDE bandwidth multiplier (<c>kernel_adjust</c>, 1), &gt; 0.</summary>
        public double KernelAdjust { get; init; } = 1.0;

        /// <summary>Starting state of the evolution models (<c>starting_state_subclonal_evolution</c>, "1:1"), Major ≥ minor, Major ≥ 1.</summary>
        public (int Major, int Minor) StartingState { get; init; } = (1, 1);

        /// <summary>
        /// R <c>set.seed</c> value for the random 8-letter mutation identifiers CNAqc draws while building the models
        /// (<c>sample(LETTERS, 8, replace = TRUE)</c>). Expected peaks, matches and decisions do not depend on it; it only
        /// decides which of several mutations with an identical expected VAF is reported (copies, role, identifier).
        /// </summary>
        public int Seed { get; init; }

        /// <summary>true reproduces R ≤ 4.3 <c>density</c> values (<c>old.coords = TRUE</c>); default R ≥ 4.4.</summary>
        public bool LegacyDensityCoordinates { get; init; }

        /// <summary>Peak-detection bootstrap replicates per segment (<c>n_bootstrap</c>, default 1), ≥ 1; drawn from the
        /// same R stream after the model identifiers (see <see cref="PurityPeakOptions.BootstrapCount"/>).</summary>
        public int BootstrapCount { get; init; } = 1;
    }

    /// <summary>One expected VAF peak of an evolutionary model (CNAqc <c>peaks_analysis$subclonal$expected_peaks</c> row).</summary>
    /// <param name="ModelId">CNAqc <c>model_id</c>, e.g. "A1B1 -&gt; A1A2B1 | A1B1" (allele genotypes).</param>
    /// <param name="Model">Linear or branching.</param>
    /// <param name="MutationId">The representative mutation's 8-letter identifier (R RNG, see <see cref="SubclonalPeakOptions.Seed"/>).</param>
    /// <param name="FirstCloneCopies">Copies of the mutation in the first-listed clone of the model (<c>n.clone_1</c>).</param>
    /// <param name="SecondCloneCopies">Copies in the second-listed clone (<c>n.clone_2</c>).</param>
    /// <param name="FirstCloneGenotype">First clone's sorted alleles (e.g. "A1A2B1"), null if the mutation is absent there.</param>
    /// <param name="SecondCloneGenotype">Second clone's sorted alleles, null if absent there.</param>
    /// <param name="Shared">CNAqc <c>role = "shared"</c> (present in both clones), else private.</param>
    /// <param name="ExpectedPeak">π·(n1·c1 + n2·c2) / (2(1−π) + π·(c1·ploidy1 + c2·ploidy2)).</param>
    /// <param name="Matched">Some data peak (discarded ones included) lies within ε (|x − peak| ≤ ε).</param>
    public sealed record SubclonalExpectedPeak(
        string ModelId,
        SubclonalEvolutionModel Model,
        string MutationId,
        int FirstCloneCopies,
        int SecondCloneCopies,
        string? FirstCloneGenotype,
        string? SecondCloneGenotype,
        bool Shared,
        double ExpectedPeak,
        bool Matched);

    /// <summary>A model's share of matched expected peaks (CNAqc decision table row).</summary>
    public readonly record struct SubclonalModelScore(string ModelId, SubclonalEvolutionModel Model, double MatchedProportion);

    /// <summary>Per-segment result of <see cref="AnalyzeSubclonalPurityPeaks"/>.</summary>
    /// <param name="Segment">The input segment.</param>
    /// <param name="Density">Gaussian KDE of the segment's VAFs.</param>
    /// <param name="Peaks">KDE data peaks.</param>
    /// <param name="ExpectedPeaks">Expected peaks of all models (branching, then linear first→second, then second→first;
    /// each model's peaks ascending). Empty when no model can reach the karyotypes from the starting state (CNAqc aborts
    /// that segment's expectations).</param>
    /// <param name="Rankings">All models by descending matched proportion (ties by model_id, C collation).</param>
    /// <param name="BestModels">The models with the top proportion — CNAqc <c>peaks_analysis$subclonal$summary</c>.</param>
    public sealed record SubclonalSegmentPeakResult(
        SubclonalPeakSegment Segment,
        KernelDensityEstimate Density,
        IReadOnlyList<PurityDataPeak> Peaks,
        IReadOnlyList<SubclonalExpectedPeak> ExpectedPeaks,
        IReadOnlyList<SubclonalModelScore> Rankings,
        IReadOnlyList<SubclonalModelScore> BestModels);

    /// <summary>
    /// CNAqc peak QC of subclonal simple copy-number segments (CNAqc 1.1.5 <c>analyze_peaks</c> →
    /// <c>analyze_peaks_subclonal</c> with <c>cluster_subclonal_CCF = FALSE</c>, <c>R/equations.R</c>
    /// <c>expectations_subclonal</c>). Each segment with n &gt; MinMutations is analysed on its own:
    /// <list type="number">
    /// <item>Evolution models from the starting state (default 1:1, alleles A1 B1): states evolve by single-allele
    /// amplification / deletion and genome doubling (breadth-first, ploidy ≤ 2·target ploidy) until the target karyotype
    /// appears; distinct allele sets are kept and every allele gains a new mutation at each reached state. Branching:
    /// both clones from the start; linear: clone 1 then clone 2 from it (only when clone 1 has no LOH or both have LOH),
    /// and clone 2 then clone 1 (symmetric condition, cell fraction 1 − CCF).</item>
    /// <item>Expected peaks of a clone pair: for every mutation with copies n1, n2 in the two clones,
    /// π·(n1·CCF + n2·(1 − CCF)) / (2(1−π) + π·(CCF·ploidy1 + (1 − CCF)·ploidy2)), distinct values, ascending;
    /// solutions with identical peak vectors (15 significant digits) are dropped.</item>
    /// <item>Data peaks: KDE peaks of the segment's VAFs (no BMix); an expected peak is matched iff |x − peak| ≤ ε for
    /// some data peak; each model scores prop = matched / #peaks; the top-prop models are the decision.</item>
    /// </list>
    /// A starting state with LOH cannot reach karyotypes without LOH: such a segment gets no models (CNAqc aborts it).
    /// CCF clustering of segments (<c>cluster_subclonal_CCF = TRUE</c>, mclust) is not ported.
    /// </summary>
    /// <param name="segments">Subclonal segments in CNAqc order.</param>
    /// <param name="purity">The purity π ∈ (0, 1] to QC.</param>
    /// <param name="options">Parameters (null = defaults).</param>
    /// <returns>One result per analysed segment, in input order.</returns>
    /// <exception cref="ArgumentNullException">A required argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid purity, option, CCF or VAF.</exception>
    /// <exception cref="ArgumentException">A segment karyotype is not simple, or the starting state is invalid.</exception>
    /// <exception cref="InvalidOperationException">A model search cannot terminate (CNAqc's <c>evolve</c> loops forever,
    /// e.g. from a 2:2 start to 1:0) or a segment yields no KDE peak (R error).</exception>
    public static IReadOnlyList<SubclonalSegmentPeakResult> AnalyzeSubclonalPurityPeaks(
        IEnumerable<SubclonalPeakSegment> segments,
        double purity,
        SubclonalPeakOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(segments);
        options ??= SubclonalPeakOptions.Default;
        if (!(purity > 0.0 && purity <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(purity), purity, "Purity must be in (0, 1].");
        if (!(options.Epsilon > 0.0 && options.Epsilon < 1.0))
            throw new ArgumentOutOfRangeException(nameof(options), options.Epsilon, "Epsilon must be in (0, 1).");
        if (!(options.KernelAdjust > 0.0) || double.IsPositiveInfinity(options.KernelAdjust))
            throw new ArgumentOutOfRangeException(nameof(options), options.KernelAdjust, "KernelAdjust must be finite and > 0.");
        if (options.MinMutations < 0)
            throw new ArgumentOutOfRangeException(nameof(options), options.MinMutations, "MinMutations must be ≥ 0.");
        if (options.BootstrapCount < 1)
            throw new ArgumentOutOfRangeException(nameof(options), options.BootstrapCount, "BootstrapCount must be ≥ 1.");
        (int startMajor, int startMinor) = options.StartingState;
        if (startMajor < 1 || startMinor < 0 || startMinor > startMajor)
            throw new ArgumentException("StartingState must satisfy Major ≥ 1 and 0 ≤ minor ≤ Major.", nameof(options));

        var all = segments.ToList();
        foreach (SubclonalPeakSegment s in all)
        {
            if (s is null) throw new ArgumentException("Segments must not contain null.", nameof(segments));
            if (s.Vafs is null) throw new ArgumentException("Segment VAFs must not be null.", nameof(segments));
            if (!PurityPeakOptions.SimpleClonalKaryotypes.Contains((s.MajorCopyNumber, s.MinorCopyNumber))
                || !PurityPeakOptions.SimpleClonalKaryotypes.Contains((s.SecondMajorCopyNumber, s.SecondMinorCopyNumber)))
                throw new ArgumentException("Subclonal segment karyotypes must be simple (1:0, 1:1, 2:0, 2:1, 2:2).", nameof(segments));
            if (!(s.Ccf > 0.0 && s.Ccf < 1.0))
                throw new ArgumentOutOfRangeException(nameof(segments), s.Ccf, "Segment CCF must be in (0, 1).");
            foreach (double v in s.Vafs)
            {
                if (double.IsNaN(v) || v < 0.0 || v > 1.0)
                    throw new ArgumentOutOfRangeException(nameof(segments), v, "VAF must be in [0, 1].");
            }
        }

        // subclonal_calls = filter(n > n_min); expectations for every segment first (one R RNG stream), then data peaks.
        var calls = all.Where(s => s.Vafs.Count > options.MinMutations).ToList();
        var rng = new RMersenneTwister(options.Seed);
        string start = $"{startMajor}:{startMinor}";
        var expectations = new List<List<CnaqcModelPeak>?>(calls.Count);
        foreach (SubclonalPeakSegment s in calls)
        {
            expectations.Add(CnaqcExpectationsSubclonal(start, s.Ccf,
                $"{s.MajorCopyNumber}:{s.MinorCopyNumber}", $"{s.SecondMajorCopyNumber}:{s.SecondMinorCopyNumber}", purity, rng));
        }

        var results = new List<SubclonalSegmentPeakResult>(calls.Count);
        for (int i = 0; i < calls.Count; i++)
        {
            SubclonalPeakSegment s = calls[i];
            (KernelDensityEstimate density, List<PurityDataPeak> peaks) = SimplePeakDetector(
                s.Vafs, options.KernelAdjust, options.LegacyDensityCoordinates, $"subclonal segment {s.Chromosome}:{s.Start}",
                options.BootstrapCount, rng);
            var expected = new List<SubclonalExpectedPeak>();
            foreach (CnaqcModelPeak e in expectations[i] ?? new List<CnaqcModelPeak>())
            {
                bool hit = peaks.Any(p => Math.Abs(p.X - e.Peak) <= options.Epsilon);
                expected.Add(new SubclonalExpectedPeak(e.ModelId, e.Model, e.MutationId, e.N1, e.N2, e.Genotype1, e.Genotype2,
                    e.Genotype1 is not null && e.Genotype2 is not null, e.Peak, hit));
            }

            // group_by(model_id) (C collation) → prop = sum(matched) / n() → arrange(desc(prop)) (stable).
            var rankings = expected
                .GroupBy(e => e.ModelId, StringComparer.Ordinal)
                .OrderBy(g => g.Key, StringComparer.Ordinal)
                .Select(g => new SubclonalModelScore(g.Key, g.First().Model, (double)g.Count(e => e.Matched) / g.Count()))
                .OrderByDescending(r => r.MatchedProportion)
                .ToList();
            var best = rankings.Where(r => rankings.Count > 0 && r.MatchedProportion == rankings[0].MatchedProportion).ToList();
            results.Add(new SubclonalSegmentPeakResult(s, density, peaks, expected, rankings, best));
        }

        return results;
    }

    // One row of CNAqc expectations_subclonal (get_peaks + model annotation).
    private sealed record CnaqcModelPeak(
        string ModelId, SubclonalEvolutionModel Model, string MutationId, int N1, int N2, string? Genotype1, string? Genotype2, double Peak);

    // A copy state: alleles (e.g. "A1", "B2") in dplyr arrange(allele) order, each with the mutation IDs it carries.
    private sealed record CnaqcAllele(string Name, string[] Mutations);

    private static int CnaqcAlleleIndex(string allele) =>
        int.Parse(allele.Replace("A", string.Empty, StringComparison.Ordinal).Replace("B", string.Empty, StringComparison.Ordinal),
            System.Globalization.CultureInfo.InvariantCulture);

    private static List<CnaqcAllele> CnaqcArrange(IEnumerable<CnaqcAllele> state) =>
        state.OrderBy(a => a.Name, StringComparer.Ordinal).ToList(); // stable, C collation

    // expectations_subclonal(starting, CCF_1, karyotype_1, karyotype_2, purity); null when CNAqc aborts (no model).
    private static List<CnaqcModelPeak>? CnaqcExpectationsSubclonal(
        string starting, double ccf, string karyotype1, string karyotype2, double purity, RMersenneTwister rng)
    {
        static bool AnyZeroPart(string k) => k.Split(':').Any(p => p == "0");
        static bool NoZeroPart(string k) => k.Split(':').All(p => p != "0");
        if (AnyZeroPart(starting) && (NoZeroPart(karyotype1) || NoZeroPart(karyotype2))) return null; // cli_abort

        var rows = new List<CnaqcModelPeak>();
        rows.AddRange(CnaqcBranching(starting, karyotype1, karyotype2, ccf, purity, rng));
        bool zero1 = karyotype1.Contains('0', StringComparison.Ordinal), zero2 = karyotype2.Contains('0', StringComparison.Ordinal);
        if ((zero1 && zero2) || !zero1) rows.AddRange(CnaqcLinear(starting, karyotype1, karyotype2, ccf, purity, rng));
        if ((zero2 && zero1) || !zero2) rows.AddRange(CnaqcLinear(starting, karyotype2, karyotype1, 1 - ccf, purity, rng));
        return rows;
    }

    private static IEnumerable<CnaqcModelPeak> CnaqcBranching(
        string starting, string left, string right, double ccf, double purity, RMersenneTwister rng)
    {
        List<CnaqcAllele> start = CnaqcInitialState(starting, rng);
        List<List<CnaqcAllele>> branchLeft = CnaqcEvolve(start, left, rng);
        List<List<CnaqcAllele>> branchRight = CnaqcEvolve(start, right, rng);
        var solutions = new List<List<CnaqcPeakRow>>();
        foreach (List<CnaqcAllele> l in branchLeft)
        {
            foreach (List<CnaqcAllele> r in branchRight) solutions.Add(CnaqcGetPeaks(l, r, ccf, purity));
        }

        return CnaqcAnnotate(start, CnaqcDistinctSolutions(solutions), SubclonalEvolutionModel.Branching, " | ");
    }

    private static IEnumerable<CnaqcModelPeak> CnaqcLinear(
        string starting, string first, string second, double ccf, double purity, RMersenneTwister rng)
    {
        List<CnaqcAllele> start = CnaqcInitialState(starting, rng);
        List<List<CnaqcAllele>> firstChildren = CnaqcEvolve(start, first, rng);
        var secondChildren = firstChildren.Select(c => CnaqcEvolve(c, second, rng)).ToList();
        var solutions = new List<List<CnaqcPeakRow>>();
        for (int i = 0; i < firstChildren.Count; i++)
        {
            foreach (List<CnaqcAllele> y in secondChildren[i]) solutions.Add(CnaqcGetPeaks(firstChildren[i], y, ccf, purity));
        }

        return CnaqcAnnotate(start, CnaqcDistinctSolutions(solutions), SubclonalEvolutionModel.Linear, " -> ");
    }

    private sealed record CnaqcPeakRow(string MutationId, int N1, int N2, string? Genotype1, string? Genotype2, double Peak);

    // solutions[!duplicated(peak %>% paste(collapse = ';'))]: R as.character(double) keeps 15 significant digits.
    private static List<List<CnaqcPeakRow>> CnaqcDistinctSolutions(List<List<CnaqcPeakRow>> solutions)
    {
        var seen = new HashSet<string>(StringComparer.Ordinal);
        var kept = new List<List<CnaqcPeakRow>>();
        foreach (List<CnaqcPeakRow> s in solutions)
        {
            string id = string.Join(";", s.Select(r => r.Peak.ToString("E14", System.Globalization.CultureInfo.InvariantCulture)));
            if (seen.Add(id)) kept.Add(s);
        }

        return kept;
    }

    private static IEnumerable<CnaqcModelPeak> CnaqcAnnotate(
        List<CnaqcAllele> start, List<List<CnaqcPeakRow>> solutions, SubclonalEvolutionModel model, string separator)
    {
        string initial = CnaqcGenotype(start);
        foreach (List<CnaqcPeakRow> s in solutions)
        {
            // paste0 drops nothing: an empty unique(genotype) vector is recycled to "".
            string g1 = s.Select(r => r.Genotype1).FirstOrDefault(g => g is not null) ?? string.Empty;
            string g2 = s.Select(r => r.Genotype2).FirstOrDefault(g => g is not null) ?? string.Empty;
            string modelId = $"{initial} -> {g1}{separator}{g2}";
            foreach (CnaqcPeakRow r in s)
                yield return new CnaqcModelPeak(modelId, model, r.MutationId, r.N1, r.N2, r.Genotype1, r.Genotype2, r.Peak);
        }
    }

    private static string CnaqcGenotype(List<CnaqcAllele> state) =>
        string.Concat(state.Select(a => a.Name).OrderBy(n => n, StringComparer.Ordinal));

    // as_karyotype: allele-letter counts sorted decreasingly, "a:b" (a single letter → "a:0").
    private static string CnaqcKaryotype(List<CnaqcAllele> state)
    {
        var counts = state.GroupBy(a => a.Name[0]).Select(g => g.Count()).OrderByDescending(c => c).ToList();
        string text = string.Join(":", counts);
        return text.Contains(':', StringComparison.Ordinal) ? text : text + ":0";
    }

    // get_peaks(clone_1, clone_2, CCF_1, purity)
    private static List<CnaqcPeakRow> CnaqcGetPeaks(List<CnaqcAllele> clone1, List<CnaqcAllele> clone2, double ccf, double purity)
    {
        static List<(string Id, int N)> Table(List<CnaqcAllele> state) =>
            state.SelectMany(a => a.Mutations).GroupBy(m => m, StringComparer.Ordinal)
                .Select(g => (g.Key, g.Count())).OrderBy(t => t.Key, StringComparer.Ordinal).ToList();

        List<(string Id, int N)> m1 = Table(clone1), m2 = Table(clone2);
        string g1 = CnaqcGenotype(clone1), g2 = CnaqcGenotype(clone2);
        double denominator = MixtureCopiesPerCell(purity, (ccf * clone1.Count) + ((1 - ccf) * clone2.Count));
        var in2 = m2.ToDictionary(t => t.Id, t => t.N, StringComparer.Ordinal);
        var in1 = new HashSet<string>(m1.Select(t => t.Id), StringComparer.Ordinal);

        var joined = new List<(string Id, int N1, int N2, string? G1, string? G2, double Raw)>();
        foreach ((string id, int n1) in m1) // full_join: x rows (with their y match), then unmatched y rows
        {
            bool shared = in2.TryGetValue(id, out int n2);
            double x1 = n1 * ccf, x2 = shared ? n2 * (1 - ccf) : 0.0;
            joined.Add((id, n1, shared ? n2 : 0, g1, shared ? g2 : null, (x1 + x2) * purity));
        }

        foreach ((string id, int n2) in m2)
        {
            if (!in1.Contains(id)) joined.Add((id, 0, n2, null, g2, (0.0 + (n2 * (1 - ccf))) * purity));
        }

        var seen = new HashSet<double>();
        return joined.Where(j => seen.Add(j.Raw))                                   // distinct(peak)
            .Select(j => new CnaqcPeakRow(j.Id, j.N1, j.N2, j.G1, j.G2, j.Raw / denominator))
            .OrderBy(r => r.Peak)                                                   // arrange(peak), stable
            .ToList();
    }

    // initial_state(target): A1, B1 with one mutation each; evolve to the target (first reached state) unless 1:1.
    private static List<CnaqcAllele> CnaqcInitialState(string target, RMersenneTwister rng)
    {
        var state = CnaqcMutation(new List<CnaqcAllele> { new("A1", Array.Empty<string>()), new("B1", Array.Empty<string>()) }, rng);
        return target == "1:1" ? state : CnaqcEvolve(state, target, rng)[0];
    }

    // mutation(copy_state): every allele gains a fresh identifier not yet used in the state.
    private static List<CnaqcAllele> CnaqcMutation(List<CnaqcAllele> state, RMersenneTwister rng)
    {
        var rows = new List<CnaqcAllele>(state);
        for (int i = 0; i < rows.Count; i++)
        {
            var used = new HashSet<string>(rows.SelectMany(a => a.Mutations), StringComparer.Ordinal);
            string id;
            do
            {
                var chars = new char[8];
                for (int k = 0; k < 8; k++) chars[k] = (char)('A' + rng.UnifIndex(26)); // sample(LETTERS, 8, replace = TRUE)
                id = new string(chars);
            }
            while (used.Contains(id));

            rows[i] = rows[i] with { Mutations = rows[i].Mutations.Append(id).ToArray() };
        }

        return rows;
    }

    // evolve(copy_state, target)
    private static List<List<CnaqcAllele>> CnaqcEvolve(List<CnaqcAllele> state, string target, RMersenneTwister rng)
    {
        if (CnaqcKaryotype(state) == target) return new List<List<CnaqcAllele>> { CnaqcMutation(state, rng) };

        int cap = 2 * target.Split(':').Sum(p => int.Parse(p, System.Globalization.CultureInfo.InvariantCulture));
        var current = new List<List<CnaqcAllele>> { state };
        List<List<CnaqcAllele>> reached;
        while (true)
        {
            var amp = current.SelectMany(CnaqcAmplify).Where(s => s.Count <= cap).ToList();
            var del = current.SelectMany(CnaqcDelete).Where(s => s.Count <= cap).ToList();
            var wgs = current.Select(CnaqcGenomeDouble).Where(s => s.Count <= cap).ToList();
            current = amp.Concat(del).Concat(wgs).ToList();
            if (current.Count == 0)
            {
                throw new InvalidOperationException(
                    $"No evolution path reaches {target} within ploidy {cap} (CNAqc's evolve() does not terminate here).");
            }

            reached = current.Where(s => CnaqcKaryotype(s) == target).ToList();
            if (reached.Count > 0) break;
        }

        var identities = new HashSet<string>(StringComparer.Ordinal);
        return reached.Where(s => identities.Add(CnaqcGenotype(s))).Select(s => CnaqcMutation(s, rng)).ToList();
    }

    // amplify: one state per allele, with a renamed copy of that allele (next free index of its letter).
    private static IEnumerable<List<CnaqcAllele>> CnaqcAmplify(List<CnaqcAllele> state)
    {
        var names = state.Select(a => a.Name).ToList();
        foreach (CnaqcAllele allele in state)
        {
            char letter = allele.Name[0];
            int n = CnaqcAlleleIndex(allele.Name);
            string created = $"{letter}{n + 1}";
            while (names.Contains(created))
            {
                n++;
                created = $"{letter}{n + 1}";
            }

            var copies = state.Where(a => a.Name == allele.Name).Select(a => a with { Name = created });
            yield return CnaqcArrange(state.Concat(copies));
        }
    }

    // delete: one state per allele, without it.
    private static IEnumerable<List<CnaqcAllele>> CnaqcDelete(List<CnaqcAllele> state)
    {
        foreach (CnaqcAllele allele in state)
            yield return CnaqcArrange(state.Where(a => a.Name != allele.Name));
    }

    // genome_double: every allele duplicated under the next free index of its letter.
    private static List<CnaqcAllele> CnaqcGenomeDouble(List<CnaqcAllele> state)
    {
        var copy = new List<CnaqcAllele>(state);
        for (int i = 0; i < copy.Count; i++)
        {
            char letter = state[i].Name[0];
            int n = int.Parse(state[i].Name.AsSpan(1), System.Globalization.CultureInfo.InvariantCulture);
            string created = $"{letter}{n + 1}";
            while (state.Any(a => a.Name == created) || copy.Any(a => a.Name == created))
            {
                n++;
                created = $"{letter}{n + 1}";
            }

            copy[i] = copy[i] with { Name = created };
        }

        return CnaqcArrange(state.Concat(copy));
    }

    // ---- BMix bmixfit (Binomial mixtures) as called by CNAqc mixture_peak_detector (FIN-B24 F63) ----
}
