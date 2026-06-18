using FsCheck;
using FsCheck.Fluent;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the Oncology algorithm group.
///
/// This file is the single home for the entire Oncology block of checklist 01
/// (rows #86–#120). Each test unit lives in its own <c>#region</c> so siblings can
/// be appended without disturbing the others. Oracles are derived INDEPENDENTLY
/// from the cited theory/doc (never routed through production constants), so a
/// self-consistent-but-wrong production constant is still caught.
///
/// Test Units: ONCO-IMMUNE-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Oncology")]
public class OncologyProperties
{
    #region ONCO-IMMUNE-001 — Immune Infiltration Estimation

    // -------------------------------------------------------------------------
    // Theory oracles (literal published constants — NOT routed through production)
    // -------------------------------------------------------------------------

    /// <summary>Yoshihara et al. (2013) tumor-purity intercept coefficient (literal).</summary>
    private const double YoshiharaA = 0.6049872018;

    /// <summary>Yoshihara et al. (2013) tumor-purity ESTIMATE-score coefficient (literal).</summary>
    private const double YoshiharaB = 0.0001467884;

    /// <summary>
    /// Independent tumor-purity oracle transcribed LITERALLY from Yoshihara et al. (2013):
    /// <c>P = clamp(cos(a + b · S_ESTIMATE), 0, 1)</c>. Constants are local literals so that a
    /// wrong production constant (shared with the implementation) is caught.
    /// </summary>
    private static double ExpectedTumorPurity(double estimateScore) =>
        Math.Clamp(Math.Cos(YoshiharaA + YoshiharaB * estimateScore), 0.0, 1.0);

    // Small gene-symbol pool: a handful of real ESTIMATE immune-signature genes plus
    // some non-signature filler, so random profiles have variable signature overlap.
    private static readonly string[] ImmunePoolGenes =
    {
        "LCP2", "PTPRC", "CD2", "CD3D", "GZMB", "PRF1", "HLA-B", "CD27",
    };

    private static readonly string[] StromalPoolGenes =
    {
        "DCN", "THBS2", "COL1A2", "FAP", "BGN", "LUM", "VCAM1", "CD200",
    };

    private static readonly string[] FillerGenes =
    {
        "GAPDH", "ACTB", "TUBB", "RPLP0", "FOO1", "BAR2", "BAZ3", "QUX4",
    };

    /// <summary>All genes the random-profile generator may draw from.</summary>
    private static readonly string[] ProfilePoolGenes =
        ImmunePoolGenes.Concat(StromalPoolGenes).Concat(FillerGenes).ToArray();

    /// <summary>
    /// Generates a random expression profile: a non-empty subset of the gene pool, each mapped
    /// to a finite, strictly-positive expression value. Distinct genes only (dictionary keys).
    /// </summary>
    private static Arbitrary<IReadOnlyDictionary<string, double>> ExpressionProfileArbitrary() =>
        (from count in Gen.Choose(1, ProfilePoolGenes.Length)
         from genes in Gen.Elements(ProfilePoolGenes).ArrayOf(count).Select(a => a.Distinct().ToArray())
         from values in Gen.Choose(1, 100_000).Select(v => v / 100.0).ArrayOf(genes.Length)
         select (IReadOnlyDictionary<string, double>)genes
             .Select((g, i) => (g, values[i]))
             .ToDictionary(t => t.g, t => t.Item2))
        .ToArbitrary();

    /// <summary>
    /// Builds a deconvolution scenario with TWO cell types whose marker genes are DISJOINT
    /// (orthogonal signature columns), a mixing fraction w ∈ (0,1), and the mixture
    /// <c>m = w·col(A) + (1−w)·col(B)</c>. Because the columns are orthogonal, NNLS recovers the
    /// mixing weights exactly (up to normalization), giving a true exact-recovery anchor.
    /// </summary>
    private static Arbitrary<(IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>> sig,
                              IReadOnlyDictionary<string, double> mixture,
                              double wA)> OrthogonalMixtureArbitrary() =>
        (from wPermille in Gen.Choose(50, 950) // wA ∈ [0.05, 0.95], away from degenerate ends
         from aGenes in Gen.Choose(2, 4)
         from bGenes in Gen.Choose(2, 4)
         from aVals in Gen.Choose(10, 200).Select(v => v / 10.0).ArrayOf(aGenes)
         from bVals in Gen.Choose(10, 200).Select(v => v / 10.0).ArrayOf(bGenes)
         select BuildOrthogonalMixture(wPermille / 1000.0, aVals, bVals))
        .ToArbitrary();

    private static (IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>> sig,
                    IReadOnlyDictionary<string, double> mixture,
                    double wA)
        BuildOrthogonalMixture(double wA, double[] aVals, double[] bVals)
    {
        var colA = new Dictionary<string, double>();
        for (int i = 0; i < aVals.Length; i++)
            colA[$"A_marker_{i}"] = aVals[i];

        var colB = new Dictionary<string, double>();
        for (int i = 0; i < bVals.Length; i++)
            colB[$"B_marker_{i}"] = bVals[i];

        var sig = new Dictionary<string, IReadOnlyDictionary<string, double>>
        {
            ["CellType_A"] = colA,
            ["CellType_B"] = colB,
        };

        // m = wA·col(A) + (1−wA)·col(B) over the disjoint marker genes.
        var mixture = new Dictionary<string, double>();
        foreach (var (g, v) in colA)
            mixture[g] = wA * v;
        foreach (var (g, v) in colB)
            mixture[g] = (1.0 - wA) * v;

        return (sig, mixture, wA);
    }

    // -------------------------------------------------------------------------
    // ESTIMATE side
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-ESTIMATE-01: <c>EstimateScore == ImmuneScore + StromalScore</c> for arbitrary profiles.
    /// The aggregate ESTIMATE score is, by construction, the additive sum of the two ssGSEA scores
    /// (Yoshihara et al., 2013). Recomputed independently from the two returned component scores.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EstimateInfiltration_EstimateScore_IsSumOfComponents()
    {
        return Prop.ForAll(ExpressionProfileArbitrary(), profile =>
        {
            var r = ImmuneAnalyzer.EstimateInfiltration(profile);
            return (Math.Abs(r.EstimateScore - (r.ImmuneScore + r.StromalScore)) < 1e-9)
                .Label($"EstimateScore {r.EstimateScore} != Immune {r.ImmuneScore} + Stromal {r.StromalScore}");
        });
    }

    /// <summary>
    /// TumorPurity formula + INV-ESTIMATE-02 (R: score ∈ [0,1]): tumor purity equals
    /// <c>clamp(cos(0.6049872018 + 0.0001467884·EstimateScore), 0, 1)</c>, recomputed INDEPENDENTLY
    /// from the returned EstimateScore using literal Yoshihara (2013) constants, and lies in [0,1].
    /// NOTE: the ssGSEA Immune/Stromal/Estimate scores are NOT bounded to [0,1] and are deliberately
    /// not asserted as such — only TumorPurity is the bounded infiltration-derived score.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EstimateInfiltration_TumorPurity_MatchesYoshiharaCosineAndIsInUnitRange()
    {
        return Prop.ForAll(ExpressionProfileArbitrary(), profile =>
        {
            var r = ImmuneAnalyzer.EstimateInfiltration(profile);
            double expected = ExpectedTumorPurity(r.EstimateScore);
            bool matchesFormula = Math.Abs(r.TumorPurity - expected) < 1e-9;
            bool inRange = r.TumorPurity is >= 0.0 and <= 1.0;
            return (matchesFormula && inRange)
                .Label($"purity {r.TumorPurity} vs oracle {expected} (estimate {r.EstimateScore}), inRange={inRange}");
        });
    }

    /// <summary>
    /// INV-ESTIMATE-03 + overlap: <c>OverlappingImmuneGenes</c> / <c>OverlappingStromalGenes</c> equal
    /// the count of the supplied signature genes present in the profile, recomputed independently;
    /// and a custom immune/stromal gene list with ZERO profile overlap yields that enrichment score 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EstimateInfiltration_OverlapCounts_MatchIndependentIntersection()
    {
        // Custom immune signature drawn from the immune pool, custom stromal from the stromal pool,
        // so overlap with a profile is exactly the independently-computed intersection size.
        return Prop.ForAll(ExpressionProfileArbitrary(), profile =>
        {
            IReadOnlyList<string> immune = ImmunePoolGenes;
            IReadOnlyList<string> stromal = StromalPoolGenes;

            var r = ImmuneAnalyzer.EstimateInfiltration(profile, immune, stromal);

            int expectedImmune = immune.Count(g => profile.ContainsKey(g));
            int expectedStromal = stromal.Count(g => profile.ContainsKey(g));

            return (r.OverlappingImmuneGenes == expectedImmune && r.OverlappingStromalGenes == expectedStromal)
                .Label($"overlap immune {r.OverlappingImmuneGenes}/{expectedImmune}, " +
                       $"stromal {r.OverlappingStromalGenes}/{expectedStromal}");
        });
    }

    /// <summary>
    /// INV-ESTIMATE-03 (zero-overlap): custom signature lists sharing NO genes with the profile
    /// produce zero immune AND stromal enrichment scores (the ssGSEA helper returns 0 on an empty
    /// hit set), hence EstimateScore 0 and tumor purity at the formula's score-0 value ≈ 0.8225.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EstimateInfiltration_ZeroSignatureOverlap_YieldsZeroScores()
    {
        return Prop.ForAll(ExpressionProfileArbitrary(), profile =>
        {
            // Signatures use genes guaranteed absent from the profile pool.
            IReadOnlyList<string> immune = new[] { "ABSENT_IMM_1", "ABSENT_IMM_2", "ABSENT_IMM_3" };
            IReadOnlyList<string> stromal = new[] { "ABSENT_STR_1", "ABSENT_STR_2", "ABSENT_STR_3" };

            var r = ImmuneAnalyzer.EstimateInfiltration(profile, immune, stromal);

            return (r.ImmuneScore == 0.0 && r.StromalScore == 0.0 && r.EstimateScore == 0.0
                    && r.OverlappingImmuneGenes == 0 && r.OverlappingStromalGenes == 0
                    && Math.Abs(r.TumorPurity - ExpectedTumorPurity(0.0)) < 1e-9)
                .Label($"non-zero score on zero overlap: immune={r.ImmuneScore}, stromal={r.StromalScore}, " +
                       $"purity={r.TumorPurity}");
        });
    }

    /// <summary>
    /// Edge anchor: an EMPTY expression profile returns zero immune/stromal/ESTIMATE scores and
    /// tumor purity ≈ 0.8225 — asserted via the literal Yoshihara cosine oracle at score 0.
    /// Source: Immune_Infiltration_Estimation.md §3.3 / §6.1.
    /// </summary>
    [Test]
    [Category("Property")]
    public void EstimateInfiltration_EmptyProfile_ZeroScoresAndCosineZeroPurity()
    {
        var r = ImmuneAnalyzer.EstimateInfiltration(new Dictionary<string, double>());
        Assert.Multiple(() =>
        {
            Assert.That(r.ImmuneScore, Is.EqualTo(0.0));
            Assert.That(r.StromalScore, Is.EqualTo(0.0));
            Assert.That(r.EstimateScore, Is.EqualTo(0.0));
            Assert.That(r.OverlappingImmuneGenes, Is.EqualTo(0));
            Assert.That(r.OverlappingStromalGenes, Is.EqualTo(0));
            Assert.That(r.TumorPurity, Is.EqualTo(ExpectedTumorPurity(0.0)).Within(1e-12));
            // Absolute literature value: cos(0.6049872018) ≈ 0.8225.
            Assert.That(r.TumorPurity, Is.EqualTo(0.8225).Within(1e-3));
        });
    }

    /// <summary>
    /// D (determinism): identical input yields an identical <see cref="ImmuneAnalyzer.InfiltrationResult"/>
    /// across repeated calls (all fields bit-identical).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EstimateInfiltration_IsDeterministic()
    {
        return Prop.ForAll(ExpressionProfileArbitrary(), profile =>
        {
            var a = ImmuneAnalyzer.EstimateInfiltration(profile);
            var b = ImmuneAnalyzer.EstimateInfiltration(profile);
            return a.Equals(b).Label($"non-deterministic InfiltrationResult: {a} vs {b}");
        });
    }

    /// <summary>
    /// Edge: a null expression profile throws <see cref="ArgumentNullException"/>.
    /// </summary>
    [Test]
    [Category("Property")]
    public void EstimateInfiltration_NullProfile_Throws() =>
        Assert.Throws<ArgumentNullException>(() => ImmuneAnalyzer.EstimateInfiltration(null!));

    // -------------------------------------------------------------------------
    // NNLS side
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-NNLS-01 (fractions ≥ 0) and P (Σ fractions ≤ 1.0): for arbitrary orthogonal mixtures every
    /// cell fraction is non-negative and the fractions sum to at most 1 (here exactly 1, since there is
    /// overlap with positive fitted mass — normalization divides by the positive total).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DeconvoluteImmuneCells_FractionsNonNegative_AndSumLeqOne()
    {
        return Prop.ForAll(OrthogonalMixtureArbitrary(), s =>
        {
            var r = ImmuneAnalyzer.DeconvoluteImmuneCells(s.mixture, s.sig);
            bool nonNeg = r.CellFractions.Values.All(v => v >= 0.0);
            double sum = r.CellFractions.Values.Sum();
            bool leqOne = sum <= 1.0 + 1e-9;
            bool sumToOne = Math.Abs(sum - 1.0) < 1e-9; // positive fitted mass ⇒ normalized to 1
            return (nonNeg && leqOne && sumToOne)
                .Label($"nonNeg={nonNeg}, sum={sum}");
        });
    }

    /// <summary>
    /// NNLS exact-recovery anchor (the real correctness test): for an orthogonal two-cell-type signature
    /// and mixture <c>m = wA·col(A) + (1−wA)·col(B)</c>, the returned fractions recover
    /// {A: wA, B: 1−wA} after normalization, with Correlation ≈ 1 and Rmse ≈ 0. Orthogonal columns make
    /// the NNLS solution unique and exact, so this validates the solver, not just its bounds.
    /// Source: linear mixture model m = S·f (Immune_Infiltration_Estimation.md §2.B).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DeconvoluteImmuneCells_OrthogonalMixture_RecoversFractionsExactly()
    {
        return Prop.ForAll(OrthogonalMixtureArbitrary(), s =>
        {
            var r = ImmuneAnalyzer.DeconvoluteImmuneCells(s.mixture, s.sig);
            double fa = r.CellFractions["CellType_A"];
            double fb = r.CellFractions["CellType_B"];
            bool recoverA = Math.Abs(fa - s.wA) < 1e-6;
            bool recoverB = Math.Abs(fb - (1.0 - s.wA)) < 1e-6;
            bool corr = Math.Abs(r.Correlation - 1.0) < 1e-6;
            bool rmse = Math.Abs(r.Rmse) < 1e-6;
            return (recoverA && recoverB && corr && rmse)
                .Label($"wA={s.wA}: fa={fa}, fb={fb}, corr={r.Correlation}, rmse={r.Rmse}");
        });
    }

    /// <summary>
    /// Single-cell-type anchor: a mixture equal to a single signature column is fully attributed to that
    /// cell type (fraction == 1) once normalized, with perfect reconstruction.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DeconvoluteImmuneCells_SingleCellType_FractionIsOne()
    {
        var sig = new Dictionary<string, IReadOnlyDictionary<string, double>>
        {
            ["Only_A"] = new Dictionary<string, double> { ["A1"] = 5.0, ["A2"] = 8.0, ["A3"] = 3.0 },
            ["Only_B"] = new Dictionary<string, double> { ["B1"] = 4.0, ["B2"] = 9.0, ["B3"] = 2.0 },
        };
        var mixture = new Dictionary<string, double> { ["A1"] = 5.0, ["A2"] = 8.0, ["A3"] = 3.0 };

        var r = ImmuneAnalyzer.DeconvoluteImmuneCells(mixture, sig);

        Assert.Multiple(() =>
        {
            Assert.That(r.CellFractions["Only_A"], Is.EqualTo(1.0).Within(1e-6));
            Assert.That(r.CellFractions["Only_B"], Is.EqualTo(0.0).Within(1e-6));
            Assert.That(r.Correlation, Is.EqualTo(1.0).Within(1e-6));
            Assert.That(r.Rmse, Is.EqualTo(0.0).Within(1e-6));
        });
    }

    /// <summary>
    /// No-overlap branch: a profile sharing NO genes with the signature matrix yields all-zero
    /// fractions, Correlation == 0, Rmse == 0, OverlappingGenes == 0, and Σ fractions == 0.
    /// Source: Immune_Infiltration_Estimation.md §6.1.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DeconvoluteImmuneCells_NoOverlap_AllZero()
    {
        var sig = new Dictionary<string, IReadOnlyDictionary<string, double>>
        {
            ["CellType_A"] = new Dictionary<string, double> { ["SIG_A1"] = 5.0, ["SIG_A2"] = 8.0 },
            ["CellType_B"] = new Dictionary<string, double> { ["SIG_B1"] = 4.0, ["SIG_B2"] = 9.0 },
        };
        IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>> signature = sig;

        // Profile genes are disjoint from every signature gene.
        var profileGen = (from count in Gen.Choose(1, 5)
                          from vals in Gen.Choose(1, 1000).Select(v => v / 10.0).ArrayOf(count)
                          select (IReadOnlyDictionary<string, double>)Enumerable.Range(0, count)
                              .ToDictionary(i => $"NONSIG_{i}", i => vals[i]))
            .ToArbitrary();

        return Prop.ForAll(profileGen, profile =>
        {
            var r = ImmuneAnalyzer.DeconvoluteImmuneCells(profile, signature);
            bool allZero = r.CellFractions.Values.All(v => v == 0.0);
            return (allZero && r.Correlation == 0.0 && r.Rmse == 0.0 && r.OverlappingGenes == 0
                    && Math.Abs(r.CellFractions.Values.Sum()) < 1e-12)
                .Label($"no-overlap branch broken: corr={r.Correlation}, rmse={r.Rmse}, overlap={r.OverlappingGenes}");
        });
    }

    /// <summary>
    /// INV-NNLS-03: <c>Rmse ≥ 0</c> always, and <c>Correlation ∈ [-1, 1]</c> whenever defined,
    /// across arbitrary orthogonal mixtures.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DeconvoluteImmuneCells_RmseNonNegative_CorrelationInUnitInterval()
    {
        return Prop.ForAll(OrthogonalMixtureArbitrary(), s =>
        {
            var r = ImmuneAnalyzer.DeconvoluteImmuneCells(s.mixture, s.sig);
            // Pearson r is mathematically in [-1, 1]; perfect reconstruction can produce 1 + ε
            // from floating-point rounding, so the bound carries a 1e-9 fp tolerance only.
            return (r.Rmse >= 0.0 && r.Correlation is >= -1.0 - 1e-9 and <= 1.0 + 1e-9)
                .Label($"rmse={r.Rmse}, corr={r.Correlation:R}");
        });
    }

    /// <summary>
    /// D (determinism): identical input yields an identical <see cref="ImmuneAnalyzer.DeconvolutionResult"/>
    /// — scalar fields bit-identical and CellFractions equal entrywise.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DeconvoluteImmuneCells_IsDeterministic()
    {
        return Prop.ForAll(OrthogonalMixtureArbitrary(), s =>
        {
            var a = ImmuneAnalyzer.DeconvoluteImmuneCells(s.mixture, s.sig);
            var b = ImmuneAnalyzer.DeconvoluteImmuneCells(s.mixture, s.sig);

            bool scalarsEqual = a.Correlation.Equals(b.Correlation)
                                && a.Rmse.Equals(b.Rmse)
                                && a.OverlappingGenes == b.OverlappingGenes
                                && a.CellFractions.Count == b.CellFractions.Count;
            bool fractionsEqual = a.CellFractions.All(kvp =>
                b.CellFractions.TryGetValue(kvp.Key, out double v) && v.Equals(kvp.Value));

            return (scalarsEqual && fractionsEqual)
                .Label("non-deterministic DeconvolutionResult");
        });
    }

    /// <summary>
    /// Edge: a null expression profile throws <see cref="ArgumentNullException"/>.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DeconvoluteImmuneCells_NullProfile_Throws() =>
        Assert.Throws<ArgumentNullException>(() => ImmuneAnalyzer.DeconvoluteImmuneCells(null!));

    #endregion
}
