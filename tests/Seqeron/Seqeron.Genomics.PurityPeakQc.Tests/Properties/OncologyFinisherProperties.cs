// SPDX-License-Identifier: GPL-3.0-or-later
// Tests of the Port of CNAqc 1.1.5, BMix 1.0.0 and R kmeans; part of Seqeron.Genomics.PurityPeakQc.Tests, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the production project).
using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// FIN-B24 heavy tier — property-based tests for the behaviour added by the B24 finisher (F44–F64). Oracles are
/// independent restatements of the reference definitions (never routed through production helpers):
/// <list type="bullet">
/// <item>LICHeE cluster summaries (F42/F43): with zero SDs and robust clusters the per-cluster margins collapse to ε,
/// i.e. the result equals <see cref="OncologyAnalyzer.ReconstructPhylogeny"/>; a bijective id relabelling maps the tree.</item>
/// <item>maftools <c>math.score</c> filters (F56): raising the VAF cutoff can only drop VAFs, so once a sample is skipped
/// (null) every higher cutoff skips it; cutoff 0 with n_min 1 is the unfiltered MATH.</item>
/// <item>Shannon clonal diversity (F57): proportion-only (permutation/scale invariant), 0 ≤ H ≤ ln k.</item>
/// <item>CNAqc peak QC (F34/F62/F64): expected peaks are <c>expected_vaf_fun</c> m·π/(2(1 − π) + π·n) for the QC'd π, so
/// a purity shift moves them exactly along that curve; the verdicts depend on the multiset of VAFs (mutation order is
/// irrelevant); the subclonal R seed only names mutations; BMix auto mode = on iff read counts are supplied.</item>
/// </list>
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Oncology")]
public class OncologyFinisherProperties
{

    private static List<PurityPeakAnalyzer.PurityPeakMutation> SimulateMutations(
        Random rng, double purity, (int Major, int Minor)[] karyotypes, int perKaryotype, bool withReads)
    {
        var list = new List<PurityPeakAnalyzer.PurityPeakMutation>();
        foreach (var (major, minor) in karyotypes)
        {
            for (int i = 0; i < perKaryotype; i++)
            {
                int m = rng.Next(2) == 0 ? 1 : Math.Max(1, major);
                double mean = ExpectedVaf(m, purity, major + minor);
                int depth = 120;
                int alt = 0;
                for (int r = 0; r < depth; r++)
                {
                    if (rng.NextDouble() < mean) alt++;
                }

                alt = Math.Max(1, alt);
                var mutation = new PurityPeakAnalyzer.PurityPeakMutation((double)alt / depth, major, minor);
                list.Add(withReads ? mutation with { AlternateReads = alt, Depth = depth } : mutation);
            }
        }

        return list;
    }

    /// <summary>
    /// F34: for any QC'd purity q the expected peaks are CNAqc <c>expected_vaf_fun</c> m·q/(2(1 − q) + q·(Major + minor))
    /// for m ∈ {1, Major} — a purity shift moves every expected peak along that curve — and the analysis depends only on
    /// the multiset of mutations (score, verdict and matched flags are permutation-invariant; KDE peaks to the 1e−12
    /// rounding of the density sums).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 40)]
    public Property AnalyzePurityPeaks_ExpectedPeaksFollowPurity_PermutationInvariant()
    {
        return Prop.ForAll(PeakScenarioGen(SimplePool).ToArbitrary(), s =>
        {
            var rng = new Random(s.Seed);
            var mutations = SimulateMutations(rng, s.Purity, s.Karyotypes, 120, withReads: false);
            var a = PurityPeakAnalyzer.AnalyzePurityPeaks(mutations, s.QcPurity);
            var b = PurityPeakAnalyzer.AnalyzePurityPeaks(mutations.OrderBy(_ => rng.Next()).ToList(), s.QcPurity);

            bool expected = a.Matches.All(m =>
                (m.Multiplicity == 1 || m.Multiplicity == m.MajorCopyNumber)
                && Math.Abs(m.ExpectedPeak - ExpectedVaf(m.Multiplicity, s.QcPurity, m.MajorCopyNumber + m.MinorCopyNumber)) <= 1e-15);
            bool sameScore = double.IsNaN(a.Score) ? double.IsNaN(b.Score) : Math.Abs(a.Score - b.Score) <= 1e-12;
            bool invariant = sameScore && a.Pass == b.Pass
                             && a.Matches.Select(m => (m.MajorCopyNumber, m.MinorCopyNumber, m.Multiplicity, m.Matched))
                                 .SequenceEqual(b.Matches.Select(m => (m.MajorCopyNumber, m.MinorCopyNumber, m.Multiplicity, m.Matched)));
            return (a.Matches.Count > 0 && expected && invariant)
                .Label($"π={s.Purity}, q={s.QcPurity}, karyotypes=[{string.Join(",", s.Karyotypes)}], score {a.Score} vs {b.Score}, pass {a.Pass} vs {b.Pass}");
        });
    }

    /// <summary>
    /// F62 (<c>analyze_peaks_general</c>): for every analysed complex karyotype the expected peaks are m = 1..Major with
    /// m·q/(2(1 − q) + q·(Major + minor)); matched + mismatched = #expected; the matched proportion and PASS
    /// (prop ≥ 0.5) are invariant to the mutation order.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 40)]
    public Property AnalyzeComplexKaryotypePeaks_ExpectedPeaks_Bookkeeping_PermutationInvariant()
    {
        return Prop.ForAll(PeakScenarioGen(ComplexPool).ToArbitrary(), s =>
        {
            var rng = new Random(s.Seed);
            var mutations = SimulateMutations(rng, s.Purity, s.Karyotypes, 130, withReads: false);
            var a = PurityPeakAnalyzer.AnalyzeComplexKaryotypePeaks(mutations, s.QcPurity);
            var b = PurityPeakAnalyzer.AnalyzeComplexKaryotypePeaks(mutations.OrderBy(_ => rng.Next()).ToList(), s.QcPurity);

            bool ok = a.Ran && a.Karyotypes.Count == s.Karyotypes.Length;
            foreach (var k in a.Karyotypes)
            {
                var exp = k.ExpectedPeaks;
                ok &= exp.Select(e => e.Multiplicity).SequenceEqual(Enumerable.Range(1, k.MajorCopyNumber))
                      && exp.All(e => Math.Abs(e.ExpectedPeak - ExpectedVaf(e.Multiplicity, s.QcPurity, k.MajorCopyNumber + k.MinorCopyNumber)) <= 1e-15)
                      && k.MatchedPeaks + k.MismatchedPeaks == exp.Count
                      && k.MatchedPeaks == exp.Count(e => e.Matched)
                      && k.Pass == (k.MatchedProportion >= 0.5);
            }

            ok &= a.Karyotypes.Select(k => (k.MajorCopyNumber, k.MinorCopyNumber, k.MatchedProportion, k.Pass))
                .SequenceEqual(b.Karyotypes.Select(k => (k.MajorCopyNumber, k.MinorCopyNumber, k.MatchedProportion, k.Pass)));
            return ok.Label($"π={s.Purity}, q={s.QcPurity}, karyotypes=[{string.Join(",", s.Karyotypes)}]");
        });
    }

    /// <summary>
    /// F62 (<c>analyze_peaks_subclonal</c>): the R seed only draws the 8-letter mutation identifiers, so expected-peak
    /// values, matches, rankings and the best models do not depend on it; permuting a segment's VAFs leaves them unchanged.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 30)]
    public Property AnalyzeSubclonalPurityPeaks_SeedAndVafOrder_DoNotChangeDecisions()
    {
        var pairs = new[] { ((2, 1), (1, 1)), ((2, 0), (1, 1)), ((1, 0), (1, 1)), ((2, 1), (2, 0)), ((1, 1), (2, 2)) };
        var arb = (from purityPct in Gen.Choose(40, 95)
                   from ccfPct in Gen.Choose(20, 80)
                   from pair in Gen.Elements(pairs)
                   from seed in Gen.Choose(0, int.MaxValue)
                   from rSeed in Gen.Choose(0, 100_000)
                   select (purityPct / 100.0, ccfPct / 100.0, pair, seed, rSeed)).ToArbitrary();
        return Prop.ForAll(arb, t =>
        {
            (double purity, double ccf, var pair, int seed, int rSeed) = t;
            var rng = new Random(seed);
            ((int ma, int mi), (int ma2, int mi2)) = pair;
            double tumourCopies = ccf * (ma + mi) + (1 - ccf) * (ma2 + mi2);
            var vafs = Enumerable.Range(0, 150).Select(_ =>
            {
                double mean = (rng.Next(2) == 0 ? ccf : 1.0) * purity / (2 * (1 - purity) + purity * tumourCopies);
                return Math.Clamp(mean + (rng.NextDouble() - 0.5) * 0.06, 0.01, 0.99);
            }).ToList();
            var segment = new PurityPeakAnalyzer.SubclonalPeakSegment("1", 1, 1_000_000, ma, mi, ma2, mi2, ccf, vafs);
            var shuffled = segment with { Vafs = vafs.OrderBy(_ => rng.Next()).ToList() };

            string Run(PurityPeakAnalyzer.SubclonalPeakSegment seg, int rs)
            {
                try
                {
                    var r = PurityPeakAnalyzer.AnalyzeSubclonalPurityPeaks(new[] { seg }, purity,
                        PurityPeakAnalyzer.SubclonalPeakOptions.Default with { Seed = rs });
                    return string.Join(";", r.Select(x =>
                        string.Join(",", x.ExpectedPeaks.Select(e => (e.ModelId, e.ExpectedPeak.ToString("R"), e.Matched)).OrderBy(e => e.ToString(), StringComparer.Ordinal))
                        + "|" + string.Join(",", x.Rankings) + "|" + string.Join(",", x.BestModels)));
                }
                catch (InvalidOperationException ex)
                {
                    return "IOE:" + ex.GetType().Name;
                }
            }

            string baseRun = Run(segment, rSeed);
            return (baseRun == Run(segment, rSeed + 1) && baseRun == Run(shuffled, rSeed))
                .Label($"π={purity}, ccf={ccf}, pair={pair}");
        });
    }

    /// <summary>
    /// F64 (CNAqc <c>analyze_peaks</c> always runs BMix when NV/DP exist): with <see cref="PurityPeakAnalyzer.PurityPeakOptions.FitMixturePeaks"/>
    /// null, VAF-only input equals FitMixturePeaks = false and read-count input equals FitMixturePeaks = true (same seed).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 8)]
    public Property AnalyzePurityPeaks_MixtureAutoDefault_FollowsReadCountAvailability()
    {
        return Prop.ForAll(PeakScenarioGen(new[] { (1, 1), (2, 1) }).ToArbitrary(), s =>
        {
            var rng = new Random(s.Seed);
            var withReads = SimulateMutations(rng, s.Purity, s.Karyotypes.Take(1).ToArray(), 110, withReads: true);
            var vafOnly = withReads.Select(m => new PurityPeakAnalyzer.PurityPeakMutation(m.Vaf, m.MajorCopyNumber, m.MinorCopyNumber)).ToList();
            var dflt = PurityPeakAnalyzer.PurityPeakOptions.Default with { Seed = s.Seed % 1000 };

            static string Key(PurityPeakAnalyzer.PurityPeakAnalysis r) =>
                $"{r.Score:R}|{r.Pass}|" + string.Join(";", r.Karyotypes.SelectMany(k => k.Peaks).Select(p => $"{p.X:R},{p.Source},{p.Discarded}"));

            bool vafSide = Key(PurityPeakAnalyzer.AnalyzePurityPeaks(vafOnly, s.QcPurity, dflt))
                           == Key(PurityPeakAnalyzer.AnalyzePurityPeaks(vafOnly, s.QcPurity, dflt with { FitMixturePeaks = false }));
            bool readSide = Key(PurityPeakAnalyzer.AnalyzePurityPeaks(withReads, s.QcPurity, dflt))
                            == Key(PurityPeakAnalyzer.AnalyzePurityPeaks(withReads, s.QcPurity, dflt with { FitMixturePeaks = true }));
            return (vafSide && readSide).Label($"vafSide={vafSide}, readSide={readSide}");
        });
    }


    // -------------------------------------------------------------------------
    // CNAqc peak QC (F34, F62, F64)
    // -------------------------------------------------------------------------

    private static double ExpectedVaf(int m, double purity, int copies) => m * purity / (2 * (1 - purity) + purity * copies);

    private static Gen<(double Purity, double QcPurity, (int Major, int Minor)[] Karyotypes, int Seed)> PeakScenarioGen(
        (int Major, int Minor)[] pool) =>
        from purityPct in Gen.Choose(35, 95)
        from shift in Gen.Choose(-15, 15)
        from count in Gen.Choose(1, 2)
        from picks in Gen.Shuffle(pool)
        from seed in Gen.Choose(0, int.MaxValue)
        select (purityPct / 100.0, Math.Clamp(purityPct + shift, 20, 100) / 100.0, picks.Take(count).ToArray(), seed);


    private static readonly (int, int)[] SimplePool = { (1, 1), (2, 0), (2, 1), (2, 2), (1, 0) };

    private static readonly (int, int)[] ComplexPool = { (3, 1), (3, 0), (4, 2), (3, 2), (4, 1) };
}
