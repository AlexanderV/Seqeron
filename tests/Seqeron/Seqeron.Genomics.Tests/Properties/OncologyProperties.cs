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
/// Test Units: ONCO-IMMUNE-001, ONCO-SOMATIC-001, ONCO-VAF-001
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

    #region ONCO-SOMATIC-001 — Somatic Mutation Calling

    // -------------------------------------------------------------------------
    // Independent oracles (transcribed from Somatic_Mutation_Calling.md §2.2/§2.4,
    // NOT routed through production constants). The decision rule, VAF, and score
    // are recomputed here from literals so a self-consistent-but-wrong production
    // constant is still caught.
    // -------------------------------------------------------------------------

    /// <summary>Documented tumor detection threshold τ_t = 0.05 (literal, Yan et al. 2021).</summary>
    private const double TauT = 0.05;

    /// <summary>Documented normal absence ceiling τ_n = 0.01 (literal, Saunders et al. 2012).</summary>
    private const double TauN = 0.01;

    /// <summary>
    /// INV-05 VAF oracle: f = totalReads == 0 ? 0 : altReads / totalReads. Recomputed independently.
    /// </summary>
    private static double ExpectedVaf(int altReads, int totalReads) =>
        totalReads == 0 ? 0.0 : (double)altReads / totalReads;

    /// <summary>
    /// INV-01/INV-02 status oracle, transcribed literally from the §4.2 decision table:
    /// f_t &lt; τ_t ⇒ NotDetected; f_t ≥ τ_t ∧ f_n ≤ τ_n ⇒ Somatic; f_t ≥ τ_t ∧ f_n &gt; τ_n ⇒ Germline.
    /// </summary>
    private static OncologyAnalyzer.SomaticStatus ExpectedStatus(double fT, double fN, double tauT, double tauN)
    {
        if (fT < tauT)
        {
            return OncologyAnalyzer.SomaticStatus.NotDetected;
        }

        return fN <= tauN
            ? OncologyAnalyzer.SomaticStatus.Somatic
            : OncologyAnalyzer.SomaticStatus.Germline;
    }

    /// <summary>INV-03 score oracle: the allele-frequency separation max(0, f_t − f_n).</summary>
    private static double ExpectedSeparationScore(double fT, double fN) => Math.Max(0.0, fT - fN);

    // -------------------------------------------------------------------------
    // Generators
    // -------------------------------------------------------------------------

    /// <summary>
    /// Generates a read-count pair (alt, total) honouring the contract 0 ≤ alt ≤ total. Includes the
    /// uncovered site (total = 0) so the VAF-0 branch (INV-05) is exercised.
    /// </summary>
    private static Gen<(int alt, int total)> ReadCountGen() =>
        from total in Gen.Choose(0, 500)
        from alt in Gen.Choose(0, total)
        select (alt, total);

    /// <summary>
    /// Generates a contract-valid <see cref="OncologyAnalyzer.VariantObservation"/> (0 ≤ alt ≤ total in
    /// both tumor and normal), with tumor totals drawn ≥ 1 so f_t can land on either side of τ_t.
    /// </summary>
    private static Arbitrary<OncologyAnalyzer.VariantObservation> VariantArbitrary() =>
        (from tTotal in Gen.Choose(1, 500)
         from tAlt in Gen.Choose(0, tTotal)
         from nTotal in Gen.Choose(0, 500)
         from nAlt in Gen.Choose(0, Math.Max(0, nTotal))
         from pos in Gen.Choose(1, 1_000_000)
         select new OncologyAnalyzer.VariantObservation(
             "chr1", pos, "A", "T", tAlt, tTotal, nAlt, nTotal))
        .ToArbitrary();

    /// <summary>
    /// Generates a variant together with random thresholds τ_t, τ_n ∈ [0, 1], so the status oracle is
    /// exercised against arbitrary cutoffs rather than only the defaults.
    /// </summary>
    private static Arbitrary<(OncologyAnalyzer.VariantObservation variant, double tauT, double tauN)>
        VariantWithThresholdsArbitrary() =>
        (from variant in VariantArbitrary().Generator
         from tt in Gen.Choose(0, 1000).Select(x => x / 1000.0)
         from tn in Gen.Choose(0, 1000).Select(x => x / 1000.0)
         select (variant, tt, tn))
        .ToArbitrary();

    /// <summary>Generates a list (0..12) of contract-valid variants, preserving order.</summary>
    private static Arbitrary<OncologyAnalyzer.VariantObservation[]> VariantListArbitrary() =>
        (from n in Gen.Choose(0, 12)
         from arr in VariantArbitrary().Generator.ArrayOf(n)
         select arr)
        .ToArbitrary();

    // -------------------------------------------------------------------------
    // VAF (R: VAF ∈ [0,1], INV-05)
    // -------------------------------------------------------------------------

    /// <summary>
    /// R + INV-05: <c>TumorVaf</c> and <c>NormalVaf</c> equal the independently recomputed
    /// altReads/totalReads (0 when totalReads = 0), and both lie in [0, 1] (since 0 ≤ alt ≤ total).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Somatic_Vaf_MatchesOracleAndInUnitRange()
    {
        return Prop.ForAll(VariantArbitrary(), v =>
        {
            var call = OncologyAnalyzer.Classify(v);
            double expectedT = ExpectedVaf(v.TumorAltReads, v.TumorTotalReads);
            double expectedN = ExpectedVaf(v.NormalAltReads, v.NormalTotalReads);
            bool matches = Math.Abs(call.TumorVaf - expectedT) < 1e-12
                           && Math.Abs(call.NormalVaf - expectedN) < 1e-12;
            bool inRange = call.TumorVaf is >= 0.0 and <= 1.0 && call.NormalVaf is >= 0.0 and <= 1.0;
            return (matches && inRange)
                .Label($"tVaf {call.TumorVaf}/{expectedT}, nVaf {call.NormalVaf}/{expectedN}, inRange={inRange}");
        });
    }

    // -------------------------------------------------------------------------
    // Status decision rule (INV-01, INV-02)
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-01/INV-02: <c>Status</c> matches the independent §4.2 decision oracle EXACTLY, driven over
    /// random read counts AND random thresholds τ_t, τ_n ∈ [0, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Somatic_Status_MatchesDecisionOracle()
    {
        return Prop.ForAll(VariantWithThresholdsArbitrary(), t =>
        {
            var call = OncologyAnalyzer.Classify(t.variant, t.tauT, t.tauN);
            double fT = ExpectedVaf(t.variant.TumorAltReads, t.variant.TumorTotalReads);
            double fN = ExpectedVaf(t.variant.NormalAltReads, t.variant.NormalTotalReads);
            var expected = ExpectedStatus(fT, fN, t.tauT, t.tauN);
            return (call.Status == expected)
                .Label($"status {call.Status} != oracle {expected} (fT={fT}, fN={fN}, τt={t.tauT}, τn={t.tauN})");
        });
    }

    // -------------------------------------------------------------------------
    // P: somatic calls absent in matched normal
    // -------------------------------------------------------------------------

    /// <summary>
    /// P: every call classified <see cref="OncologyAnalyzer.SomaticStatus.Somatic"/> has its matched-normal
    /// VAF ≤ τ_n (absent in normal) AND tumor VAF ≥ τ_t (present in tumor), over random thresholds.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Somatic_SomaticCalls_AbsentInMatchedNormal()
    {
        return Prop.ForAll(VariantWithThresholdsArbitrary(), t =>
        {
            var call = OncologyAnalyzer.Classify(t.variant, t.tauT, t.tauN);
            if (call.Status != OncologyAnalyzer.SomaticStatus.Somatic)
            {
                return true.ToProperty();
            }

            return (call.NormalVaf <= t.tauN && call.TumorVaf >= t.tauT)
                .Label($"somatic but nVaf={call.NormalVaf} (>τn {t.tauN}) or tVaf={call.TumorVaf} (<τt {t.tauT})");
        });
    }

    // -------------------------------------------------------------------------
    // INV-03: SomaticScore
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-03: <c>SomaticCall.SomaticScore == (Status==Somatic ? max(0, f_t − f_n) : 0)</c> and ∈ [0, 1].
    /// The SomaticCall score is gated on the Somatic status (it is 0 for Germline/NotDetected), recomputed
    /// independently from the VAF separation. Random thresholds.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Somatic_SomaticCallScore_IsGatedSeparationAndInUnitRange()
    {
        return Prop.ForAll(VariantWithThresholdsArbitrary(), t =>
        {
            var call = OncologyAnalyzer.Classify(t.variant, t.tauT, t.tauN);
            double fT = ExpectedVaf(t.variant.TumorAltReads, t.variant.TumorTotalReads);
            double fN = ExpectedVaf(t.variant.NormalAltReads, t.variant.NormalTotalReads);
            double expected = call.Status == OncologyAnalyzer.SomaticStatus.Somatic
                ? ExpectedSeparationScore(fT, fN)
                : 0.0;
            bool matches = Math.Abs(call.SomaticScore - expected) < 1e-9;
            bool inRange = call.SomaticScore is >= 0.0 and <= 1.0;
            return (matches && inRange)
                .Label($"score {call.SomaticScore} != {expected} (status {call.Status}), inRange={inRange}");
        });
    }

    /// <summary>
    /// INV-03 (standalone): <c>CalculateSomaticScore(variant) == max(0, f_t − f_n)</c> ALWAYS — the standalone
    /// scorer is status-INDEPENDENT (no NotDetected/Germline gating), unlike the SomaticCall score — and ∈ [0, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Somatic_StandaloneScore_IsUngatedSeparationAndInUnitRange()
    {
        return Prop.ForAll(VariantArbitrary(), v =>
        {
            double score = OncologyAnalyzer.CalculateSomaticScore(v);
            double fT = ExpectedVaf(v.TumorAltReads, v.TumorTotalReads);
            double fN = ExpectedVaf(v.NormalAltReads, v.NormalTotalReads);
            double expected = ExpectedSeparationScore(fT, fN);
            return (Math.Abs(score - expected) < 1e-9 && score is >= 0.0 and <= 1.0)
                .Label($"standalone score {score} != {expected} (fT={fT}, fN={fN})");
        });
    }

    // -------------------------------------------------------------------------
    // M (REMAPPED): the checklist's "higher tumor depth → more confident" is FALSE for this rule —
    // the score is max(0, f_t − f_n), which is DEPTH-INDEPENDENT (scaling alt & total together leaves
    // the VAF, hence the score, unchanged). The honest monotonicity is in the VAF SEPARATION:
    //   • holding the normal and totals fixed, raising tumorAltReads (⇒ higher f_t) never DECREASES the
    //     SomaticScore and never turns a detected (Somatic/Germline) call into NotDetected;
    //   • holding the tumor and totals fixed, raising normalAltReads (⇒ higher f_n) never INCREASES the
    //     standalone separation score.
    // -------------------------------------------------------------------------

    /// <summary>Generator for the tumor-alt monotonicity check: shared totals, two ordered tumor-alt counts.</summary>
    private static Arbitrary<(int tTotal, int tAltLo, int tAltHi, int nTotal, int nAlt)> TumorAltMonotoneArbitrary() =>
        (from tTotal in Gen.Choose(1, 500)
         from a1 in Gen.Choose(0, tTotal)
         from a2 in Gen.Choose(0, tTotal)
         from nTotal in Gen.Choose(0, 500)
         from nAlt in Gen.Choose(0, Math.Max(0, nTotal))
         select (tTotal, Math.Min(a1, a2), Math.Max(a1, a2), nTotal, nAlt))
        .ToArbitrary();

    /// <summary>
    /// M (remapped, tumor side): with normal and totals fixed, increasing tumorAltReads (⇒ higher f_t) never
    /// DECREASES the standalone separation score and never demotes a detected call to NotDetected.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Somatic_Score_MonotoneNonDecreasingInTumorAlt()
    {
        return Prop.ForAll(TumorAltMonotoneArbitrary(), t =>
        {
            var lo = new OncologyAnalyzer.VariantObservation("chr1", 1, "A", "T", t.tAltLo, t.tTotal, t.nAlt, t.nTotal);
            var hi = new OncologyAnalyzer.VariantObservation("chr1", 1, "A", "T", t.tAltHi, t.tTotal, t.nAlt, t.nTotal);

            double scoreLo = OncologyAnalyzer.CalculateSomaticScore(lo);
            double scoreHi = OncologyAnalyzer.CalculateSomaticScore(hi);

            var callLo = OncologyAnalyzer.Classify(lo);
            var callHi = OncologyAnalyzer.Classify(hi);

            bool scoreMonotone = scoreHi >= scoreLo - 1e-12;
            // If the lower-alt call was detected, the higher-alt call (≥ f_t) cannot become NotDetected.
            bool detectionMonotone =
                callLo.Status == OncologyAnalyzer.SomaticStatus.NotDetected
                || callHi.Status != OncologyAnalyzer.SomaticStatus.NotDetected;

            return (scoreMonotone && detectionMonotone)
                .Label($"scoreLo={scoreLo}, scoreHi={scoreHi}, statusLo={callLo.Status}, statusHi={callHi.Status}");
        });
    }

    /// <summary>Generator for the normal-alt monotonicity check: shared tumor & totals, two ordered normal-alt counts.</summary>
    private static Arbitrary<(int tTotal, int tAlt, int nTotal, int nAltLo, int nAltHi)> NormalAltMonotoneArbitrary() =>
        (from tTotal in Gen.Choose(1, 500)
         from tAlt in Gen.Choose(0, tTotal)
         from nTotal in Gen.Choose(1, 500)
         from a1 in Gen.Choose(0, nTotal)
         from a2 in Gen.Choose(0, nTotal)
         select (tTotal, tAlt, nTotal, Math.Min(a1, a2), Math.Max(a1, a2)))
        .ToArbitrary();

    /// <summary>
    /// M (remapped, normal side): with tumor and totals fixed, increasing normalAltReads (⇒ higher f_n) never
    /// INCREASES the standalone separation score (it is non-increasing in f_n).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Somatic_Score_MonotoneNonIncreasingInNormalAlt()
    {
        return Prop.ForAll(NormalAltMonotoneArbitrary(), t =>
        {
            var lo = new OncologyAnalyzer.VariantObservation("chr1", 1, "A", "T", t.tAlt, t.tTotal, t.nAltLo, t.nTotal);
            var hi = new OncologyAnalyzer.VariantObservation("chr1", 1, "A", "T", t.tAlt, t.tTotal, t.nAltHi, t.nTotal);

            double scoreLo = OncologyAnalyzer.CalculateSomaticScore(lo);
            double scoreHi = OncologyAnalyzer.CalculateSomaticScore(hi);

            return (scoreHi <= scoreLo + 1e-12)
                .Label($"normalAlt {t.nAltLo}->{t.nAltHi}: scoreLo={scoreLo}, scoreHi={scoreHi}");
        });
    }

    // -------------------------------------------------------------------------
    // INV-04: FilterGermlineVariants
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-04: <c>FilterGermlineVariants</c> returns EXACTLY the Somatic subset of
    /// <c>CallSomaticMutations</c>, in input order (same length and same calls).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Somatic_FilterGermline_IsSomaticSubsetInOrder()
    {
        return Prop.ForAll(VariantListArbitrary(), variants =>
        {
            var filtered = OncologyAnalyzer.FilterGermlineVariants(variants);
            var expected = OncologyAnalyzer.CallSomaticMutations(variants)
                .Where(c => c.Status == OncologyAnalyzer.SomaticStatus.Somatic)
                .ToList();

            bool sameLength = filtered.Count == expected.Count;
            bool allSomatic = filtered.All(c => c.Status == OncologyAnalyzer.SomaticStatus.Somatic);
            bool sameOrder = sameLength && filtered.Zip(expected, (a, b) => a.Equals(b)).All(x => x);

            return (sameLength && allSomatic && sameOrder)
                .Label($"filtered {filtered.Count} vs expected {expected.Count}, allSomatic={allSomatic}");
        });
    }

    // -------------------------------------------------------------------------
    // Order/count preservation + D (determinism)
    // -------------------------------------------------------------------------

    /// <summary>
    /// D + order/count preservation: <c>CallSomaticMutations</c> output has the same length and order as the
    /// input (its embedded <c>Variant</c> matches positionally), and identical inputs ⇒ identical results
    /// (all fields, via record equality).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Somatic_CallSomaticMutations_PreservesOrderAndIsDeterministic()
    {
        return Prop.ForAll(VariantListArbitrary(), variants =>
        {
            var a = OncologyAnalyzer.CallSomaticMutations(variants);
            var b = OncologyAnalyzer.CallSomaticMutations(variants);

            bool lengthOk = a.Count == variants.Length;
            bool orderOk = a.Select((c, i) => c.Variant.Equals(variants[i])).All(x => x);
            bool deterministic = a.Count == b.Count && a.Zip(b, (x, y) => x.Equals(y)).All(z => z);

            return (lengthOk && orderOk && deterministic)
                .Label($"length={a.Count}/{variants.Length}, orderOk={orderOk}, deterministic={deterministic}");
        });
    }

    // -------------------------------------------------------------------------
    // Edge cases / validation
    // -------------------------------------------------------------------------

    /// <summary>Edge: null variants ⇒ <see cref="ArgumentNullException"/>.</summary>
    [Test]
    [Category("Property")]
    public void Somatic_NullVariants_Throws() =>
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CallSomaticMutations(null!));

    /// <summary>Edge: a tumor threshold outside [0, 1] ⇒ <see cref="ArgumentOutOfRangeException"/>.</summary>
    [TestCase(-0.01)]
    [TestCase(1.01)]
    [Category("Property")]
    public void Somatic_TumorThresholdOutOfRange_Throws(double tauT) =>
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            OncologyAnalyzer.CallSomaticMutations(Array.Empty<OncologyAnalyzer.VariantObservation>(), tauT, 0.01));

    /// <summary>Edge: a normal threshold outside [0, 1] ⇒ <see cref="ArgumentOutOfRangeException"/>.</summary>
    [TestCase(-0.01)]
    [TestCase(1.01)]
    [Category("Property")]
    public void Somatic_NormalThresholdOutOfRange_Throws(double tauN) =>
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            OncologyAnalyzer.CallSomaticMutations(Array.Empty<OncologyAnalyzer.VariantObservation>(), 0.05, tauN));

    /// <summary>
    /// Edge: a variant with negative reads or alt &gt; total ⇒ <see cref="ArgumentOutOfRangeException"/>,
    /// thrown by the VAF computation inside <c>Classify</c>/<c>CallSomaticMutations</c> (not the ctor).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Somatic_NegativeReads_Throws()
    {
        var negative = new OncologyAnalyzer.VariantObservation("chr1", 1, "A", "T", -1, 100, 0, 100);
        var altExceedsTotal = new OncologyAnalyzer.VariantObservation("chr1", 1, "A", "T", 120, 100, 0, 100);

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.Classify(negative));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.Classify(altExceedsTotal));
        });
    }

    /// <summary>Edge: empty input ⇒ empty output.</summary>
    [Test]
    [Category("Property")]
    public void Somatic_EmptyInput_EmptyOutput() =>
        Assert.That(
            OncologyAnalyzer.CallSomaticMutations(Array.Empty<OncologyAnalyzer.VariantObservation>()),
            Is.Empty);

    /// <summary>
    /// Edge: tumor-only (normalTotal = 0 ⇒ f_n = 0) with f_t ≥ τ_t ⇒ Somatic (ℓ_n = 1 analogue, §6.1).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Somatic_TumorOnly_AboveThreshold_IsSomatic()
    {
        var v = new OncologyAnalyzer.VariantObservation("chr1", 100, "A", "T", 30, 100, 0, 0);
        var call = OncologyAnalyzer.Classify(v);
        Assert.Multiple(() =>
        {
            Assert.That(call.NormalVaf, Is.EqualTo(0.0));
            Assert.That(call.Status, Is.EqualTo(OncologyAnalyzer.SomaticStatus.Somatic));
            Assert.That(call.SomaticScore, Is.EqualTo(0.30).Within(1e-9));
        });
    }

    // -------------------------------------------------------------------------
    // Worked-example anchors (Somatic_Mutation_Calling.md §7.1)
    // -------------------------------------------------------------------------

    /// <summary>
    /// §7.1 anchor: (25/100 tumor, 0/100 normal) ⇒ Somatic, score 0.25; (48/100 tumor, 50/100 normal)
    /// ⇒ Germline, score 0.0; (4/100 tumor) ⇒ NotDetected.
    /// </summary>
    [TestCase(25, 100, 0, 100, OncologyAnalyzer.SomaticStatus.Somatic, 0.25)]
    [TestCase(48, 100, 50, 100, OncologyAnalyzer.SomaticStatus.Germline, 0.0)]
    [TestCase(4, 100, 0, 100, OncologyAnalyzer.SomaticStatus.NotDetected, 0.0)]
    [Category("Property")]
    public void Somatic_WorkedExample_Anchors(
        int tAlt, int tTotal, int nAlt, int nTotal,
        OncologyAnalyzer.SomaticStatus expectedStatus, double expectedScore)
    {
        var v = new OncologyAnalyzer.VariantObservation("chr1", 100, "A", "T", tAlt, tTotal, nAlt, nTotal);
        var call = OncologyAnalyzer.Classify(v);
        Assert.Multiple(() =>
        {
            Assert.That(call.Status, Is.EqualTo(expectedStatus));
            Assert.That(call.SomaticScore, Is.EqualTo(expectedScore).Within(1e-9));
        });
    }

    #endregion

    #region ONCO-VAF-001 — Variant Allele Frequency Analysis

    // -------------------------------------------------------------------------
    // Independent oracles (transcribed from Variant_Allele_Frequency.md §2.2,
    // NOT routed through production constants). The empirical VAF, the Wilson
    // score interval (literal z = 1.96), and the purity/ploidy correction are
    // recomputed here from literals so a self-consistent-but-wrong production
    // constant is still caught.
    // -------------------------------------------------------------------------

    /// <summary>
    /// Literal 95% standard-normal quantile z = 1.96 (Wilson 1927, via the binomial-proportion
    /// confidence-interval specification). Declared as a LOCAL literal so the oracle never routes
    /// through the production <c>ZScore95</c> constant.
    /// </summary>
    private const double VafZScore95 = 1.96;

    /// <summary>The only supported (and default) confidence level for the Wilson interval.</summary>
    private const double VafConfidence95 = 0.95;

    /// <summary>
    /// INV-01 empirical-VAF oracle: <c>VAF = totalReads == 0 ? 0 : altReads / totalReads</c>,
    /// recomputed independently (§2.2 / §4.1).
    /// </summary>
    private static double ExpectedVafEmpirical(int altReads, int totalReads) =>
        totalReads == 0 ? 0.0 : (double)altReads / totalReads;

    /// <summary>
    /// Independent Wilson score interval oracle (§2.2) with literal z = 1.96:
    /// <c>center = (p̂ + z²/2n)/(1 + z²/n)</c>,
    /// <c>margin = (z/(1 + z²/n))·√(p̂(1−p̂)/n + z²/4n²)</c>, then bounds clamped to [0, 1].
    /// Returns (lower, upper).
    /// </summary>
    private static (double lower, double upper) ExpectedWilsonInterval(int altReads, int totalReads)
    {
        double n = totalReads;
        double pHat = (double)altReads / totalReads;
        double z2 = VafZScore95 * VafZScore95;
        double denom = 1.0 + z2 / n;
        double center = (pHat + z2 / (2.0 * n)) / denom;
        double margin = (VafZScore95 / denom) * Math.Sqrt(pHat * (1.0 - pHat) / n + z2 / (4.0 * n * n));
        return (Math.Max(0.0, center - margin), Math.Min(1.0, center + margin));
    }

    /// <summary>
    /// INV-04 purity/ploidy-correction oracle (§2.2): <c>adjusted = vaf·(2(1−π) + π·n_tot)/π</c>,
    /// with the normal contribution fixed at the literal diploid 2. Recomputed independently.
    /// </summary>
    private static double ExpectedAdjustedVaf(double vaf, double purity, double ploidy) =>
        vaf * (2.0 * (1.0 - purity) + purity * ploidy) / purity;

    // -------------------------------------------------------------------------
    // Generators
    // -------------------------------------------------------------------------

    /// <summary>
    /// (alt, total) honouring the contract 0 ≤ alt ≤ total, with total = 0 included so the VAF-0 /
    /// no-coverage branch (INV-01 / §6.1) is exercised.
    /// </summary>
    private static Arbitrary<(int alt, int total)> VafCountsArbitrary() =>
        (from total in Gen.Choose(0, 5000)
         from alt in Gen.Choose(0, total)
         select (alt, total))
        .ToArbitrary();

    /// <summary>(alt, total) with total ≥ 1, so the Wilson interval is always defined (CI requires coverage).</summary>
    private static Arbitrary<(int alt, int total)> VafCoveredCountsArbitrary() =>
        (from total in Gen.Choose(1, 5000)
         from alt in Gen.Choose(0, total)
         select (alt, total))
        .ToArbitrary();

    /// <summary>
    /// Monotonicity generator: a shared total ≥ 2 and two ORDERED, DISTINCT alt counts (lo &lt; hi),
    /// both in [0, total], so increasing alt strictly increases the empirical VAF.
    /// </summary>
    private static Arbitrary<(int total, int altLo, int altHi)> VafMonotoneArbitrary() =>
        (from total in Gen.Choose(2, 5000)
         from a1 in Gen.Choose(0, total)
         from a2 in Gen.Choose(0, total)
         where a1 != a2
         select (total, Math.Min(a1, a2), Math.Max(a1, a2)))
        .ToArbitrary();

    /// <summary>Tumor purity π ∈ (0, 1] (in 1/1000 steps, excluding 0).</summary>
    private static Arbitrary<double> PurityArbitrary() =>
        Gen.Choose(1, 1000).Select(x => x / 1000.0).ToArbitrary();

    /// <summary>Tumor total copy number n_tot &gt; 0 (in 1/10 steps).</summary>
    private static Arbitrary<double> PloidyArbitrary() =>
        Gen.Choose(1, 80).Select(x => x / 10.0).ToArbitrary();

    /// <summary>(vaf, purity, ploidy) with vaf ∈ [0, 1], purity ∈ (0, 1], ploidy &gt; 0.</summary>
    private static Arbitrary<(double vaf, double purity, double ploidy)> PurityCorrectionArbitrary() =>
        (from vafMilli in Gen.Choose(0, 1000)
         from purityMilli in Gen.Choose(1, 1000)
         from ploidyTenth in Gen.Choose(1, 80)
         select (vafMilli / 1000.0, purityMilli / 1000.0, ploidyTenth / 10.0))
        .ToArbitrary();

    // -------------------------------------------------------------------------
    // CalculateVAF (R + INV-01, M, R: VAF==0 at no alt reads, D)
    // -------------------------------------------------------------------------

    /// <summary>
    /// R + INV-01: <c>CalculateVAF</c> equals the independently recomputed
    /// <c>totalReads == 0 ? 0 : altReads/totalReads</c> and always lies in [0, 1] (since 0 ≤ alt ≤ total).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Vaf_Empirical_MatchesOracleAndInUnitRange()
    {
        return Prop.ForAll(VafCountsArbitrary(), c =>
        {
            double vaf = OncologyAnalyzer.CalculateVAF(c.alt, c.total);
            double expected = ExpectedVafEmpirical(c.alt, c.total);
            bool matches = Math.Abs(vaf - expected) < 1e-12;
            bool inRange = vaf is >= 0.0 and <= 1.0;
            return (matches && inRange)
                .Label($"vaf {vaf} vs oracle {expected} (alt={c.alt}, total={c.total}), inRange={inRange}");
        });
    }

    /// <summary>
    /// R (§6.1): <c>CalculateVAF</c> returns exactly 0 when altReads == 0 and total &gt; 0 (no alt support).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Vaf_ZeroAltReads_IsZero()
    {
        return Prop.ForAll(Gen.Choose(1, 5000).ToArbitrary(), total =>
        {
            double vaf = OncologyAnalyzer.CalculateVAF(0, total);
            return (vaf == 0.0).Label($"VAF for 0/{total} = {vaf}, expected 0");
        });
    }

    /// <summary>
    /// M (§2.4 corollary): at a FIXED totalReads, increasing altReads STRICTLY increases the VAF —
    /// more alt reads ⇒ higher VAF. Driven over ordered distinct alt counts on a shared total.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Vaf_StrictlyIncreasesWithAltReads()
    {
        return Prop.ForAll(VafMonotoneArbitrary(), t =>
        {
            double lo = OncologyAnalyzer.CalculateVAF(t.altLo, t.total);
            double hi = OncologyAnalyzer.CalculateVAF(t.altHi, t.total);
            return (hi > lo)
                .Label($"VAF not strictly increasing: {t.altLo}/{t.total}={lo} vs {t.altHi}/{t.total}={hi}");
        });
    }

    /// <summary>
    /// D (determinism): identical (alt, total) inputs ⇒ bit-identical <c>CalculateVAF</c> output.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Vaf_Empirical_IsDeterministic()
    {
        return Prop.ForAll(VafCountsArbitrary(), c =>
        {
            double a = OncologyAnalyzer.CalculateVAF(c.alt, c.total);
            double b = OncologyAnalyzer.CalculateVAF(c.alt, c.total);
            return a.Equals(b).Label($"non-deterministic CalculateVAF: {a} vs {b} for {c.alt}/{c.total}");
        });
    }

    /// <summary>§7.1 anchor: <c>CalculateVAF(25, 100) == 0.25</c>.</summary>
    [Test]
    [Category("Property")]
    public void Vaf_WorkedExample_25Of100_IsQuarter() =>
        Assert.That(OncologyAnalyzer.CalculateVAF(25, 100), Is.EqualTo(0.25).Within(1e-12));

    // -------------------------------------------------------------------------
    // Wilson score confidence interval (INV-02, INV-03)
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-02 + INV-03: the returned Wilson bounds match the independent z = 1.96 oracle (§2.2) and obey
    /// <c>0 ≤ Lower ≤ Upper ≤ 1</c>; the <c>Vaf</c> field equals the empirical p̂ and <c>Confidence == 0.95</c>.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Vaf_WilsonInterval_MatchesOracleAndIsBounded()
    {
        return Prop.ForAll(VafCoveredCountsArbitrary(), c =>
        {
            var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(c.alt, c.total);
            var (expLower, expUpper) = ExpectedWilsonInterval(c.alt, c.total);
            double expVaf = ExpectedVafEmpirical(c.alt, c.total);

            bool boundsMatch = Math.Abs(ci.Lower - expLower) < 1e-9 && Math.Abs(ci.Upper - expUpper) < 1e-9;
            bool vafIsPHat = Math.Abs(ci.Vaf - expVaf) < 1e-12;
            bool ordered = ci.Lower <= ci.Upper;
            bool bounded = ci.Lower >= 0.0 && ci.Upper <= 1.0;
            bool conf = ci.Confidence == VafConfidence95;

            return (boundsMatch && vafIsPHat && ordered && bounded && conf)
                .Label($"alt={c.alt}/{c.total}: Lower {ci.Lower}/{expLower}, Upper {ci.Upper}/{expUpper}, " +
                       $"Vaf {ci.Vaf}/{expVaf}, conf {ci.Confidence}");
        });
    }

    /// <summary>
    /// D (determinism): identical inputs ⇒ identical <see cref="OncologyAnalyzer.VafConfidenceInterval"/>
    /// (all fields, via record equality).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Vaf_WilsonInterval_IsDeterministic()
    {
        return Prop.ForAll(VafCoveredCountsArbitrary(), c =>
        {
            var a = OncologyAnalyzer.CalculateVAFConfidenceInterval(c.alt, c.total);
            var b = OncologyAnalyzer.CalculateVAFConfidenceInterval(c.alt, c.total);
            return a.Equals(b).Label($"non-deterministic CI: {a} vs {b} for {c.alt}/{c.total}");
        });
    }

    /// <summary>
    /// §7.1 worked example anchor: 25/100 at 95% (z = 1.96) ⇒ Vaf 0.25, Lower ≈ 0.1754509,
    /// Upper ≈ 0.3430465 (the exact §7.1 numbers, within 1e-6), Confidence 0.95.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Vaf_WilsonInterval_WorkedExample_25Of100()
    {
        var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(25, 100);
        Assert.Multiple(() =>
        {
            Assert.That(ci.Vaf, Is.EqualTo(0.25).Within(1e-12));
            Assert.That(ci.Lower, Is.EqualTo(0.1754509).Within(1e-6));
            Assert.That(ci.Upper, Is.EqualTo(0.3430465).Within(1e-6));
            Assert.That(ci.Confidence, Is.EqualTo(0.95));
        });
    }

    /// <summary>
    /// §6.1 edge (p̂ = 0): alt 0 with coverage ⇒ Wilson Lower == 0 and Upper &gt; 0 (no overshoot,
    /// non-zero width).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Vaf_WilsonInterval_PHatZero_LowerIsZeroUpperPositive()
    {
        var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(0, 100);
        Assert.Multiple(() =>
        {
            Assert.That(ci.Vaf, Is.EqualTo(0.0));
            Assert.That(ci.Lower, Is.EqualTo(0.0));
            Assert.That(ci.Upper, Is.GreaterThan(0.0));
            Assert.That(ci.Upper, Is.LessThanOrEqualTo(1.0));
        });
    }

    /// <summary>
    /// §6.1 edge (p̂ = 1): alt == total ⇒ Wilson Upper == 1 (no overshoot — the upper bound reaches but
    /// does not exceed 1; equals 1 in exact arithmetic, within fp tolerance) and Lower &lt; 1 (non-zero
    /// width). Per INV-03 (§2.4) the upper bound is ≤ 1; the §6.1 "upper = 1" is the exact-arithmetic value
    /// and the implementation's Min-clamp does not lift the sub-ulp fp drift up to a literal 1.0.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Vaf_WilsonInterval_PHatOne_UpperIsOneLowerBelowOne()
    {
        var ci = OncologyAnalyzer.CalculateVAFConfidenceInterval(100, 100);
        Assert.Multiple(() =>
        {
            Assert.That(ci.Vaf, Is.EqualTo(1.0));
            Assert.That(ci.Upper, Is.EqualTo(1.0).Within(1e-12)); // reaches 1 (no overshoot, INV-03)
            Assert.That(ci.Upper, Is.LessThanOrEqualTo(1.0));
            Assert.That(ci.Lower, Is.LessThan(1.0));
            Assert.That(ci.Lower, Is.GreaterThanOrEqualTo(0.0));
        });
    }

    // -------------------------------------------------------------------------
    // AdjustVAFForPurity (INV-04)
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-04 formula: <c>AdjustVAFForPurity</c> equals the independent oracle
    /// <c>vaf·(2(1−π) + π·ploidy)/π</c> (§2.2), over random vaf ∈ [0, 1], π ∈ (0, 1], ploidy &gt; 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Vaf_PurityCorrection_MatchesOracle()
    {
        return Prop.ForAll(PurityCorrectionArbitrary(), t =>
        {
            double adjusted = OncologyAnalyzer.AdjustVAFForPurity(t.vaf, t.purity, t.ploidy);
            double expected = ExpectedAdjustedVaf(t.vaf, t.purity, t.ploidy);
            return (Math.Abs(adjusted - expected) < 1e-9)
                .Label($"adjusted {adjusted} vs oracle {expected} (vaf={t.vaf}, π={t.purity}, ploidy={t.ploidy})");
        });
    }

    /// <summary>
    /// INV-04 round-trip: for a diploid heterozygous clonal SNV the expected VAF is v = π/2, so
    /// <c>AdjustVAFForPurity(π/2, π, 2) == 1</c> for random π ∈ (0, 1] (within 1e-9).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Vaf_PurityCorrection_DiploidHetRoundTripIsOne()
    {
        return Prop.ForAll(PurityArbitrary(), purity =>
        {
            double adjusted = OncologyAnalyzer.AdjustVAFForPurity(purity / 2.0, purity, 2.0);
            return (Math.Abs(adjusted - 1.0) < 1e-9)
                .Label($"round-trip ≠ 1: AdjustVAFForPurity({purity / 2.0}, {purity}, 2) = {adjusted}");
        });
    }

    /// <summary>
    /// D (determinism): identical inputs ⇒ bit-identical <c>AdjustVAFForPurity</c> output.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Vaf_PurityCorrection_IsDeterministic()
    {
        return Prop.ForAll(PurityCorrectionArbitrary(), t =>
        {
            double a = OncologyAnalyzer.AdjustVAFForPurity(t.vaf, t.purity, t.ploidy);
            double b = OncologyAnalyzer.AdjustVAFForPurity(t.vaf, t.purity, t.ploidy);
            return a.Equals(b).Label($"non-deterministic AdjustVAFForPurity: {a} vs {b}");
        });
    }

    /// <summary>§7.1 anchor: <c>AdjustVAFForPurity(0.40, 0.80, 2) == 1.0</c> (clonal het, diploid).</summary>
    [Test]
    [Category("Property")]
    public void Vaf_PurityCorrection_WorkedExample_IsOne() =>
        Assert.That(OncologyAnalyzer.AdjustVAFForPurity(0.40, 0.80, 2.0), Is.EqualTo(1.0).Within(1e-9));

    // -------------------------------------------------------------------------
    // Edge cases / validation (§3.3 / §6.1)
    // -------------------------------------------------------------------------

    /// <summary>Edge: <c>CalculateVAF</c> with altReads &lt; 0 or altReads &gt; totalReads ⇒ ArgumentOutOfRangeException.</summary>
    [Test]
    [Category("Property")]
    public void Vaf_InvalidCounts_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateVAF(-1, 100));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateVAF(120, 100));
        });
    }

    /// <summary>Edge (§6.1): <c>CalculateVAF</c> with totalReads == 0 ⇒ returns 0 (no coverage).</summary>
    [Test]
    [Category("Property")]
    public void Vaf_ZeroCoverage_ReturnsZero() =>
        Assert.That(OncologyAnalyzer.CalculateVAF(0, 0), Is.EqualTo(0.0));

    /// <summary>Edge (§6.1): the confidence interval with totalReads == 0 ⇒ ArgumentOutOfRangeException (undefined).</summary>
    [Test]
    [Category("Property")]
    public void Vaf_ConfidenceInterval_ZeroCoverage_Throws() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateVAFConfidenceInterval(0, 0));

    /// <summary>
    /// Edge (§3.3): a confidence level outside (0, 1), or any level other than the supported 0.95,
    /// ⇒ ArgumentOutOfRangeException.
    /// </summary>
    [TestCase(-0.01)]
    [TestCase(0.0)]
    [TestCase(0.90)]
    [TestCase(0.99)]
    [TestCase(1.0)]
    [TestCase(1.5)]
    [Category("Property")]
    public void Vaf_ConfidenceInterval_UnsupportedConfidence_Throws(double confidence) =>
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            OncologyAnalyzer.CalculateVAFConfidenceInterval(25, 100, confidence));

    /// <summary>Edge (§3.3): <c>AdjustVAFForPurity</c> with vaf outside [0, 1] ⇒ ArgumentOutOfRangeException.</summary>
    [TestCase(-0.01)]
    [TestCase(1.01)]
    [Category("Property")]
    public void Vaf_PurityCorrection_VafOutOfRange_Throws(double vaf) =>
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AdjustVAFForPurity(vaf, 0.8, 2.0));

    /// <summary>Edge (§3.3): <c>AdjustVAFForPurity</c> with purity outside (0, 1] ⇒ ArgumentOutOfRangeException.</summary>
    [TestCase(0.0)]
    [TestCase(-0.5)]
    [TestCase(1.01)]
    [Category("Property")]
    public void Vaf_PurityCorrection_PurityOutOfRange_Throws(double purity) =>
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AdjustVAFForPurity(0.4, purity, 2.0));

    /// <summary>Edge (§3.3): <c>AdjustVAFForPurity</c> with ploidy ≤ 0 ⇒ ArgumentOutOfRangeException.</summary>
    [TestCase(0.0)]
    [TestCase(-1.0)]
    [Category("Property")]
    public void Vaf_PurityCorrection_PloidyNonPositive_Throws(double ploidy) =>
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AdjustVAFForPurity(0.4, 0.8, ploidy));

    #endregion
}
