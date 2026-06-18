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
/// Test Units: ONCO-IMMUNE-001, ONCO-SOMATIC-001, ONCO-VAF-001, ONCO-DRIVER-001,
/// ONCO-ARTIFACT-001, ONCO-ANNOT-001, ONCO-TMB-001
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

    #region ONCO-DRIVER-001 — Driver Mutation Detection (20/20 Rule)

    // -------------------------------------------------------------------------
    // Independent oracles (transcribed from Driver_Mutation_Detection.md §2.2/§2.4,
    // NOT routed through production constants). The 20/20-rule threshold, the
    // recurrent-position min count, and the truncating consequence set are LOCAL
    // literals here, so a self-consistent-but-wrong production constant is caught.
    //
    // DOC AMBIGUITY (resolved against the SOURCE): the doc's INV-02 phrases the
    // dual-pass tie-break as "Oncogene needs f_OG ≥ f_TSG", which would label an
    // exact f_OG == f_TSG tie Oncogene; but §4.1/§5.3/§6 (and the implementation,
    // OncologyAnalyzer.ClassifyGene) treat an EXACT tie with BOTH fractions > 0.20
    // as Ambiguous. We follow the SOURCE: the property generators are constructed so
    // they NEVER emit an exact-tie dual-pass case (we assert dominant-fraction roles
    // only), and the exact-tie corner is pinned in a dedicated [Test] anchor below
    // (Driver_DualPassExactTie_IsAmbiguous) with this note.
    // -------------------------------------------------------------------------

    /// <summary>Literal 20/20-rule driver fraction threshold (strict <c>&gt;</c>), Vogelstein 2013 / Tokheim 2020.</summary>
    private const double DriverThreshold = 0.20;

    /// <summary>Literal recurrent-position minimum count: a position is recurrent at ≥ 2 missense (Miller 2017).</summary>
    private const int RecurrentMinCount = 2;

    /// <summary>
    /// INV-05 fraction oracle, recomputed INDEPENDENTLY of the production grouping:
    /// N = total mutations; recurrentMissenseCount = sum over protein positions with ≥ 2 missense of their
    /// FULL missense count; truncatingCount = #(consequence ∈ {Nonsense, Frameshift, SpliceSite});
    /// f_OG = recurrentMissenseCount / N, f_TSG = truncatingCount / N (both 0 on empty).
    /// </summary>
    private static (double fOg, double fTsg) ExpectedFractions(IReadOnlyList<OncologyAnalyzer.GeneMutation> mutations)
    {
        int n = mutations.Count;
        if (n == 0)
        {
            return (0.0, 0.0);
        }

        int truncatingCount = mutations.Count(IsTruncatingOracle);

        int recurrentMissenseCount = mutations
            .Where(m => m.Consequence == OncologyAnalyzer.MutationConsequence.Missense)
            .GroupBy(m => m.ProteinPosition)
            .Where(g => g.Count() >= RecurrentMinCount)
            .Sum(g => g.Count());

        return ((double)recurrentMissenseCount / n, (double)truncatingCount / n);
    }

    /// <summary>Independent truncating predicate: nonsense, frameshift, or splice-site (Schroeder 2014; Miller 2017).</summary>
    private static bool IsTruncatingOracle(OncologyAnalyzer.GeneMutation m) =>
        m.Consequence is OncologyAnalyzer.MutationConsequence.Nonsense
            or OncologyAnalyzer.MutationConsequence.Frameshift
            or OncologyAnalyzer.MutationConsequence.SpliceSite;

    /// <summary>
    /// INV-02/INV-03 role oracle (full 20/20 rule, transcribed from the SOURCE — NOT the doc's literal
    /// INV-02 phrasing): only f_OG &gt; 0.20 ⇒ Oncogene; only f_TSG &gt; 0.20 ⇒ TumorSuppressor; both &gt; 0.20
    /// with a STRICTLY dominant fraction ⇒ that fraction's role; both &gt; 0.20 with an EXACT tie ⇒ Ambiguous
    /// (the documented dual-pass tie behavior, §4.1/§5.3/§6); neither ⇒ Ambiguous. The dedicated-role property
    /// generators never emit an exact tie, but this oracle is also used by the multi-gene membership property
    /// (which CAN), so it must encode the exact-tie⇒Ambiguous resolution to stay an honest independent oracle.
    /// </summary>
    private static OncologyAnalyzer.DriverGeneRole ExpectedRole(double fOg, double fTsg)
    {
        bool isOg = fOg > DriverThreshold;
        bool isTsg = fTsg > DriverThreshold;

        if (isOg && isTsg)
        {
            if (fOg > fTsg)
            {
                return OncologyAnalyzer.DriverGeneRole.Oncogene;
            }

            return fTsg > fOg
                ? OncologyAnalyzer.DriverGeneRole.TumorSuppressor
                : OncologyAnalyzer.DriverGeneRole.Ambiguous; // exact tie
        }

        if (isOg)
        {
            return OncologyAnalyzer.DriverGeneRole.Oncogene;
        }

        return isTsg
            ? OncologyAnalyzer.DriverGeneRole.TumorSuppressor
            : OncologyAnalyzer.DriverGeneRole.Ambiguous;
    }

    // -------------------------------------------------------------------------
    // Generators
    // -------------------------------------------------------------------------

    /// <summary>Small gene-symbol pool, so IdentifyDriverMutations groups several mutations per gene.</summary>
    private static readonly string[] DriverGenePool = { "KRAS", "TP53", "EGFR", "BRAF", "PTEN" };

    /// <summary>All consequence categories the generator may draw, including non-driver <c>Other</c>.</summary>
    private static readonly OncologyAnalyzer.MutationConsequence[] DriverConsequencePool =
    {
        OncologyAnalyzer.MutationConsequence.Missense,
        OncologyAnalyzer.MutationConsequence.Nonsense,
        OncologyAnalyzer.MutationConsequence.Frameshift,
        OncologyAnalyzer.MutationConsequence.SpliceSite,
        OncologyAnalyzer.MutationConsequence.Other,
    };

    /// <summary>A single random mutation: a gene from the pool, a position in a SMALL range (to force recurrence), any consequence.</summary>
    private static Gen<OncologyAnalyzer.GeneMutation> DriverMutationGen(string gene) =>
        from posObj in Gen.Choose(1, 4) // tiny range ⇒ collisions ⇒ recurrent positions
        from consObj in Gen.Elements(DriverConsequencePool)
        select new OncologyAnalyzer.GeneMutation(gene, posObj, consObj);

    /// <summary>A list (1..12) of mutations for ONE gene — the contract for ClassifyGene/ScoreDriverPotential.</summary>
    private static Arbitrary<OncologyAnalyzer.GeneMutation[]> SingleGeneMutationsArbitrary() =>
        (from gene in Gen.Elements(DriverGenePool)
         from n in Gen.Choose(1, 12)
         from arr in DriverMutationGen(gene).ArrayOf(n)
         select arr)
        .ToArbitrary();

    /// <summary>A list (0..18) of mutations across SEVERAL genes — the contract for IdentifyDriverMutations.</summary>
    private static Arbitrary<OncologyAnalyzer.GeneMutation[]> MultiGeneMutationsArbitrary() =>
        (from n in Gen.Choose(0, 18)
         from arr in (from gene in Gen.Elements(DriverGenePool)
                      from m in DriverMutationGen(gene)
                      select m).ArrayOf(n)
         select arr)
        .ToArbitrary();

    /// <summary>
    /// Builds a PURE-ONCOGENE gene: <paramref name="recurrent"/> missense ALL at one position (≥ 2 ⇒ recurrent)
    /// plus <paramref name="filler"/> non-recurrent, non-truncating <c>Other</c> mutations at distinct positions.
    /// f_OG = recurrent/(recurrent+filler); f_TSG = 0. With recurrent &gt; filler·(0.25) it exceeds 0.20.
    /// </summary>
    private static OncologyAnalyzer.GeneMutation[] BuildPureOncogene(string gene, int recurrent, int filler)
    {
        var list = new List<OncologyAnalyzer.GeneMutation>();
        for (int i = 0; i < recurrent; i++)
        {
            list.Add(new OncologyAnalyzer.GeneMutation(gene, 100, OncologyAnalyzer.MutationConsequence.Missense));
        }

        for (int i = 0; i < filler; i++)
        {
            list.Add(new OncologyAnalyzer.GeneMutation(gene, 200 + i, OncologyAnalyzer.MutationConsequence.Other));
        }

        return list.ToArray();
    }

    /// <summary>
    /// Builds a PURE-TSG gene: <paramref name="truncating"/> nonsense mutations (at distinct positions, so they
    /// contribute to f_TSG only) plus <paramref name="filler"/> non-recurrent <c>Other</c> mutations.
    /// f_TSG = truncating/(truncating+filler); f_OG = 0.
    /// </summary>
    private static OncologyAnalyzer.GeneMutation[] BuildPureTsg(string gene, int truncating, int filler)
    {
        var list = new List<OncologyAnalyzer.GeneMutation>();
        for (int i = 0; i < truncating; i++)
        {
            list.Add(new OncologyAnalyzer.GeneMutation(gene, 300 + i, OncologyAnalyzer.MutationConsequence.Nonsense));
        }

        for (int i = 0; i < filler; i++)
        {
            list.Add(new OncologyAnalyzer.GeneMutation(gene, 500 + i, OncologyAnalyzer.MutationConsequence.Other));
        }

        return list.ToArray();
    }

    /// <summary>
    /// Builds an AMBIGUOUS gene: only <c>Other</c> and single (non-recurrent) missense / single truncating,
    /// so BOTH fractions are 0 (≤ 0.20). Distinct positions ⇒ no recurrence; one truncating over a large N.
    /// </summary>
    private static OncologyAnalyzer.GeneMutation[] BuildAmbiguous(string gene, int filler)
    {
        var list = new List<OncologyAnalyzer.GeneMutation>();
        // One isolated missense (not recurrent) + one truncating, drowned by Other filler.
        list.Add(new OncologyAnalyzer.GeneMutation(gene, 700, OncologyAnalyzer.MutationConsequence.Missense));
        for (int i = 0; i < filler; i++)
        {
            list.Add(new OncologyAnalyzer.GeneMutation(gene, 800 + i, OncologyAnalyzer.MutationConsequence.Other));
        }

        return list.ToArray();
    }

    // -------------------------------------------------------------------------
    // INV-05: fraction oracle
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-05: <c>RecurrentMissenseFraction</c> and <c>TruncatingFraction</c> equal the INDEPENDENTLY recomputed
    /// counts/N (recurrent = missense positions with ≥ 2 missense, full count; truncating = nonsense+frameshift+
    /// splice), within 1e-9; <c>MutationCount</c> == N; both fractions ∈ [0, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Driver_Fractions_MatchIndependentOracle()
    {
        return Prop.ForAll(SingleGeneMutationsArbitrary(), mutations =>
        {
            var c = OncologyAnalyzer.ClassifyGene(mutations);
            var (fOg, fTsg) = ExpectedFractions(mutations);

            bool ogOk = Math.Abs(c.RecurrentMissenseFraction - fOg) < 1e-9;
            bool tsgOk = Math.Abs(c.TruncatingFraction - fTsg) < 1e-9;
            bool countOk = c.MutationCount == mutations.Length;
            bool inRange = c.RecurrentMissenseFraction is >= 0.0 and <= 1.0
                           && c.TruncatingFraction is >= 0.0 and <= 1.0;

            return (ogOk && tsgOk && countOk && inRange)
                .Label($"f_OG {c.RecurrentMissenseFraction}/{fOg}, f_TSG {c.TruncatingFraction}/{fTsg}, " +
                       $"count {c.MutationCount}/{mutations.Length}");
        });
    }

    // -------------------------------------------------------------------------
    // INV-02/INV-03: role (unambiguous cases only — generators exclude exact ties)
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-02/INV-03: <c>Role</c> matches the independent 20/20-rule oracle for the UNAMBIGUOUS cases. The
    /// fraction generator is built so an exact f_OG == f_TSG dual pass is impossible: a pure-oncogene gene has
    /// f_TSG = 0, a pure-TSG gene has f_OG = 0, an ambiguous gene has both 0. Thus the role is fully determined
    /// by the single non-zero criterion (or Ambiguous when neither exceeds 0.20).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Driver_Role_MatchesOracle_NoExactTie()
    {
        // Each scenario yields f_OG=0 XOR f_TSG=0 (or both 0) ⇒ never an exact dual-pass tie over 0.20.
        var scenarioGen =
            from gene in Gen.Elements(DriverGenePool)
            from kind in Gen.Choose(0, 2)
            from major in Gen.Choose(1, 10)
            from filler in Gen.Choose(0, 10)
            select kind switch
            {
                0 => BuildPureOncogene(gene, major + 1, filler), // ≥ 2 missense ⇒ recurrent
                1 => BuildPureTsg(gene, major, filler),
                _ => BuildAmbiguous(gene, filler),
            };

        return Prop.ForAll(scenarioGen.ToArbitrary(), mutations =>
        {
            var c = OncologyAnalyzer.ClassifyGene(mutations);
            var (fOg, fTsg) = ExpectedFractions(mutations);
            var expected = ExpectedRole(fOg, fTsg);
            return (c.Role == expected)
                .Label($"role {c.Role} != oracle {expected} (f_OG={fOg}, f_TSG={fTsg})");
        });
    }

    // -------------------------------------------------------------------------
    // R / INV-04: ScoreDriverPotential
    // -------------------------------------------------------------------------

    /// <summary>
    /// R + INV-04: <c>ScoreDriverPotential</c> equals <c>max(f_OG, f_TSG)</c> (independent oracle, within 1e-9)
    /// and lies in [0, 1] — driver score ≥ 0 always.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Driver_Score_IsMaxFractionAndInUnitRange()
    {
        return Prop.ForAll(SingleGeneMutationsArbitrary(), mutations =>
        {
            double score = OncologyAnalyzer.ScoreDriverPotential(mutations);
            var (fOg, fTsg) = ExpectedFractions(mutations);
            double expected = Math.Max(fOg, fTsg);
            return (Math.Abs(score - expected) < 1e-9 && score is >= 0.0 and <= 1.0)
                .Label($"score {score} != max({fOg},{fTsg})={expected}");
        });
    }

    // -------------------------------------------------------------------------
    // M: more recurrent / hotspot → higher driver likelihood
    // -------------------------------------------------------------------------

    /// <summary>
    /// M (monotone construction): a gene whose mutations are ALL missense at ONE position (f_OG = 1.0 ⇒ score
    /// 1.0) scores ≥ ANY mixed gene of the same size that adds non-recurrent <c>Other</c> filler (which only
    /// lowers the recurrent-missense fraction). Concretely, score(all-recurrent) == 1.0 ≥ score(recurrent+filler),
    /// strictly &gt; whenever filler &gt; 0. Constructed, not a vague claim.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Driver_Score_MoreRecurrent_IsHigher()
    {
        var gen =
            from gene in Gen.Elements(DriverGenePool)
            from recurrent in Gen.Choose(2, 10) // ≥ 2 ⇒ recurrent
            from filler in Gen.Choose(0, 10)
            select (gene, recurrent, filler);

        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            // Pure recurrent: f_OG = 1.0 ⇒ score 1.0.
            var pure = BuildPureOncogene(t.gene, t.recurrent, 0);
            // Diluted: same recurrent block + non-recurrent Other filler ⇒ f_OG < 1 when filler > 0.
            var diluted = BuildPureOncogene(t.gene, t.recurrent, t.filler);

            double pureScore = OncologyAnalyzer.ScoreDriverPotential(pure);
            double dilutedScore = OncologyAnalyzer.ScoreDriverPotential(diluted);

            bool pureIsOne = Math.Abs(pureScore - 1.0) < 1e-9;
            bool monotone = pureScore >= dilutedScore - 1e-12;
            bool strictWhenDiluted = t.filler == 0 || pureScore > dilutedScore + 1e-12;

            return (pureIsOne && monotone && strictWhenDiluted)
                .Label($"pure={pureScore} (expect 1.0), diluted={dilutedScore}, filler={t.filler}");
        });
    }

    // -------------------------------------------------------------------------
    // INV-01 / P: IdentifyDriverMutations subset + MatchCancerHotspots membership
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-01 + P: <c>IdentifyDriverMutations</c> returns a SUBSET of the input in INPUT ORDER, and every
    /// returned mutation is justified — its gene is a driver gene (Oncogene/TumorSuppressor by the independent
    /// oracle) OR its (gene, position) is a known hotspot. Conversely no NON-driver, NON-hotspot mutation is
    /// returned. Driven with a random hotspot set built from the input's own (gene, position) keys.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Driver_Identify_IsJustifiedSubsetInOrder()
    {
        var gen =
            from mutations in MultiGeneMutationsArbitrary().Generator
            // Hotspot set: a random subset of the input's (gene, position) keys.
            from mask in Gen.Choose(0, 1).ArrayOf(mutations.Length)
            select (mutations, mask);

        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var hotspots = t.mutations
                .Where((_, i) => t.mask.Length > i && t.mask[i] == 1)
                .Select(m => (m.Gene, m.ProteinPosition))
                .ToHashSet();
            IReadOnlySet<(string, int)> hotspotSet = hotspots;

            var drivers = OncologyAnalyzer.IdentifyDriverMutations(t.mutations, hotspotSet);

            // Independent oracle: a gene is a driver gene iff its oracle role != Ambiguous.
            var driverGenes = t.mutations
                .GroupBy(m => m.Gene, StringComparer.Ordinal)
                .Where(g =>
                {
                    var (fOg, fTsg) = ExpectedFractions(g.ToList());
                    return ExpectedRole(fOg, fTsg) != OncologyAnalyzer.DriverGeneRole.Ambiguous;
                })
                .Select(g => g.Key)
                .ToHashSet(StringComparer.Ordinal);

            bool IsExpectedDriver(OncologyAnalyzer.GeneMutation m) =>
                driverGenes.Contains(m.Gene) || hotspots.Contains((m.Gene, m.ProteinPosition));

            // Expected = exactly the input mutations satisfying the predicate, in input order.
            var expected = t.mutations.Where(IsExpectedDriver).ToList();

            bool sameLength = drivers.Count == expected.Count;
            bool sameOrder = sameLength && drivers.Zip(expected, (a, b) => a.Equals(b)).All(x => x);
            bool allJustified = drivers.All(IsExpectedDriver);

            return (sameLength && sameOrder && allJustified)
                .Label($"drivers {drivers.Count} vs expected {expected.Count}, order={sameOrder}");
        });
    }

    /// <summary>
    /// P (construction): an Oncogene-classified gene ⇒ ALL its mutations are returned; a hotspot-only mutation
    /// in an Ambiguous gene ⇒ returned IFF its (gene, position) ∈ knownHotspots.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Driver_Identify_OncogeneAllReturned_HotspotGatesAmbiguous()
    {
        // KRAS: pure oncogene (3 missense at pos 100) ⇒ all 3 returned regardless of hotspots.
        var kras = BuildPureOncogene("KRAS", 3, 0);
        // PTEN: ambiguous (1 isolated missense + Other filler) ⇒ only hotspot mutation returned.
        var pten = BuildAmbiguous("PTEN", 3);
        var input = kras.Concat(pten).ToArray();

        // Hotspot only on the PTEN isolated missense at position 700.
        IReadOnlySet<(string, int)> hotspots = new HashSet<(string, int)> { ("PTEN", 700) };

        var drivers = OncologyAnalyzer.IdentifyDriverMutations(input, hotspots);

        Assert.Multiple(() =>
        {
            // All KRAS oncogene mutations present.
            Assert.That(drivers.Count(m => m.Gene == "KRAS"), Is.EqualTo(3));
            // Exactly the single hotspot PTEN mutation present (the Other filler is neither driver nor hotspot).
            Assert.That(drivers.Count(m => m.Gene == "PTEN"), Is.EqualTo(1));
            Assert.That(drivers.Single(m => m.Gene == "PTEN").ProteinPosition, Is.EqualTo(700));
            // Subset of input.
            Assert.That(drivers.All(input.Contains), Is.True);
        });
    }

    /// <summary>
    /// P: <c>MatchCancerHotspots</c> is true ⟺ (gene, position) ∈ the supplied set, over arbitrary mutations
    /// and a random hotspot set drawn from the gene/position space.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Driver_MatchHotspots_IsSetMembership()
    {
        var gen =
            from gene in Gen.Elements(DriverGenePool)
            from pos in Gen.Choose(1, 6)
            from cons in Gen.Elements(DriverConsequencePool)
            from inSetFlag in Gen.Choose(0, 1)
            select (mutation: new OncologyAnalyzer.GeneMutation(gene, pos, cons), inSet: inSetFlag == 1);

        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var set = new HashSet<(string, int)>();
            if (t.inSet)
            {
                set.Add((t.mutation.Gene, t.mutation.ProteinPosition));
            }

            IReadOnlySet<(string, int)> hotspots = set;
            bool matched = OncologyAnalyzer.MatchCancerHotspots(t.mutation, hotspots);
            bool expected = set.Contains((t.mutation.Gene, t.mutation.ProteinPosition));
            return (matched == expected)
                .Label($"matched={matched}, expected={expected}, inSet={t.inSet}");
        });
    }

    // -------------------------------------------------------------------------
    // D: determinism
    // -------------------------------------------------------------------------

    /// <summary>
    /// D (determinism): identical inputs ⇒ identical classification (record equality), identical score, and an
    /// identical driver list (same length, entrywise equal) across repeated calls.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Driver_IsDeterministic()
    {
        return Prop.ForAll(MultiGeneMutationsArbitrary(), mutations =>
        {
            var c1 = OncologyAnalyzer.ClassifyGene(mutations);
            var c2 = OncologyAnalyzer.ClassifyGene(mutations);
            double s1 = OncologyAnalyzer.ScoreDriverPotential(mutations);
            double s2 = OncologyAnalyzer.ScoreDriverPotential(mutations);
            var d1 = OncologyAnalyzer.IdentifyDriverMutations(mutations);
            var d2 = OncologyAnalyzer.IdentifyDriverMutations(mutations);

            bool classOk = c1.Equals(c2);
            bool scoreOk = s1.Equals(s2);
            bool driversOk = d1.Count == d2.Count && d1.Zip(d2, (a, b) => a.Equals(b)).All(x => x);

            return (classOk && scoreOk && driversOk)
                .Label($"classOk={classOk}, scoreOk={scoreOk}, driversOk={driversOk}");
        });
    }

    // -------------------------------------------------------------------------
    // Edge cases / validation (READ-from-source behavior, §3.3 / §6.1)
    // -------------------------------------------------------------------------

    /// <summary>Edge (§3.3): null mutations ⇒ <see cref="ArgumentNullException"/> on every entry point.</summary>
    [Test]
    [Category("Property")]
    public void Driver_NullMutations_Throws() =>
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.ClassifyGene(null!));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.ScoreDriverPotential(null!));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.IdentifyDriverMutations(null!));
        });

    /// <summary>
    /// Edge (SOURCE behavior): <c>IdentifyDriverMutations</c> treats a null <c>knownHotspots</c> as EMPTY
    /// (does NOT throw) — the parameter is nullable with a null-coalescing default in OncologyAnalyzer. This is
    /// the actual implementation behavior; the doc §3.3's blanket "null knownHotspots throws" applies to
    /// <c>MatchCancerHotspots</c> (whose set parameter is NON-nullable and DOES throw — asserted separately).
    /// With no hotspots, an Oncogene gene's mutations are still returned (gene-driven), proving null⇒empty
    /// rather than null⇒throw.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Driver_Identify_NullHotspots_TreatedAsEmpty()
    {
        var kras = BuildPureOncogene("KRAS", 3, 0);
        IReadOnlyList<OncologyAnalyzer.GeneMutation> drivers =
            OncologyAnalyzer.IdentifyDriverMutations(kras, null);
        // No throw; the oncogene's mutations are returned via gene classification, not via hotspots.
        Assert.That(drivers.Count, Is.EqualTo(3));
    }

    /// <summary>Edge (§3.3): <c>MatchCancerHotspots</c> with a null set ⇒ <see cref="ArgumentNullException"/>.</summary>
    [Test]
    [Category("Property")]
    public void Driver_MatchHotspots_NullSet_Throws() =>
        Assert.Throws<ArgumentNullException>(() =>
            OncologyAnalyzer.MatchCancerHotspots(
                new OncologyAnalyzer.GeneMutation("KRAS", 12, OncologyAnalyzer.MutationConsequence.Missense),
                null!));

    /// <summary>Edge (§6.1): an empty mutation list ⇒ Ambiguous, both fractions 0, count 0, score 0.</summary>
    [Test]
    [Category("Property")]
    public void Driver_EmptyList_AmbiguousZeroFractions()
    {
        var c = OncologyAnalyzer.ClassifyGene(Array.Empty<OncologyAnalyzer.GeneMutation>());
        Assert.Multiple(() =>
        {
            Assert.That(c.Role, Is.EqualTo(OncologyAnalyzer.DriverGeneRole.Ambiguous));
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.0));
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(0.0));
            Assert.That(c.MutationCount, Is.EqualTo(0));
            Assert.That(OncologyAnalyzer.ScoreDriverPotential(Array.Empty<OncologyAnalyzer.GeneMutation>()),
                Is.EqualTo(0.0));
        });
    }

    /// <summary>
    /// Edge (§6.1, strict threshold): a truncating fraction EXACTLY 0.20 is NOT a TumorSuppressor — 1 nonsense
    /// among 5 mutations (f_TSG = 0.20) classifies Ambiguous, since the rule is strict <c>&gt;</c> 0.20.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Driver_TruncatingFractionExactly020_NotTumorSuppressor()
    {
        var mutations = new[]
        {
            new OncologyAnalyzer.GeneMutation("TP53", 1, OncologyAnalyzer.MutationConsequence.Nonsense),
            new OncologyAnalyzer.GeneMutation("TP53", 2, OncologyAnalyzer.MutationConsequence.Other),
            new OncologyAnalyzer.GeneMutation("TP53", 3, OncologyAnalyzer.MutationConsequence.Other),
            new OncologyAnalyzer.GeneMutation("TP53", 4, OncologyAnalyzer.MutationConsequence.Other),
            new OncologyAnalyzer.GeneMutation("TP53", 5, OncologyAnalyzer.MutationConsequence.Other),
        };
        var c = OncologyAnalyzer.ClassifyGene(mutations);
        Assert.Multiple(() =>
        {
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.20).Within(1e-12));
            Assert.That(c.Role, Is.Not.EqualTo(OncologyAnalyzer.DriverGeneRole.TumorSuppressor));
            Assert.That(c.Role, Is.EqualTo(OncologyAnalyzer.DriverGeneRole.Ambiguous));
        });
    }

    /// <summary>
    /// Edge (§6.1): a SINGLE missense at a position is NOT recurrent — three distinct-position single missenses
    /// give f_OG = 0 (recurrence requires ≥ 2 at one position).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Driver_SingleMissensePerPosition_NotRecurrent()
    {
        var mutations = new[]
        {
            new OncologyAnalyzer.GeneMutation("EGFR", 10, OncologyAnalyzer.MutationConsequence.Missense),
            new OncologyAnalyzer.GeneMutation("EGFR", 20, OncologyAnalyzer.MutationConsequence.Missense),
            new OncologyAnalyzer.GeneMutation("EGFR", 30, OncologyAnalyzer.MutationConsequence.Missense),
        };
        var c = OncologyAnalyzer.ClassifyGene(mutations);
        Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(0.0));
    }

    /// <summary>
    /// Dual-pass EXACT-TIE corner (DOC AMBIGUITY): a gene with f_OG == f_TSG, BOTH strictly &gt; 0.20, is
    /// Ambiguous per the SOURCE (OncologyAnalyzer.ClassifyGene resolves an exact tie to Ambiguous; §4.1/§5.3/§6),
    /// NOT Oncogene as a literal reading of doc INV-02's "f_OG ≥ f_TSG" would suggest. Construction: 2 recurrent
    /// missense (one position) + 2 nonsense over N=4 ⇒ f_OG = f_TSG = 0.50.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Driver_DualPassExactTie_IsAmbiguous()
    {
        var mutations = new[]
        {
            new OncologyAnalyzer.GeneMutation("BRAF", 600, OncologyAnalyzer.MutationConsequence.Missense),
            new OncologyAnalyzer.GeneMutation("BRAF", 600, OncologyAnalyzer.MutationConsequence.Missense),
            new OncologyAnalyzer.GeneMutation("BRAF", 1, OncologyAnalyzer.MutationConsequence.Nonsense),
            new OncologyAnalyzer.GeneMutation("BRAF", 2, OncologyAnalyzer.MutationConsequence.Nonsense),
        };
        var c = OncologyAnalyzer.ClassifyGene(mutations);
        Assert.Multiple(() =>
        {
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(0.50).Within(1e-12));
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.50).Within(1e-12));
            Assert.That(c.Role, Is.EqualTo(OncologyAnalyzer.DriverGeneRole.Ambiguous));
        });
    }

    /// <summary>
    /// IDH1 anchor (§7.1): 10 missense at codon 132 ⇒ Role Oncogene, RecurrentMissenseFraction 1.0,
    /// TruncatingFraction 0.0, ScoreDriverPotential 1.0.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Driver_Idh1WorkedExample_Anchor()
    {
        var idh1 = Enumerable.Range(0, 10)
            .Select(_ => new OncologyAnalyzer.GeneMutation("IDH1", 132, OncologyAnalyzer.MutationConsequence.Missense))
            .ToArray();
        var c = OncologyAnalyzer.ClassifyGene(idh1);
        Assert.Multiple(() =>
        {
            Assert.That(c.Role, Is.EqualTo(OncologyAnalyzer.DriverGeneRole.Oncogene));
            Assert.That(c.RecurrentMissenseFraction, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(c.TruncatingFraction, Is.EqualTo(0.0));
            Assert.That(OncologyAnalyzer.ScoreDriverPotential(idh1), Is.EqualTo(1.0).Within(1e-12));
        });
    }

    #endregion

    #region ONCO-ARTIFACT-001 — Sequencing Artifact Detection (OxoG / FFPE / Strand Bias)

    // -------------------------------------------------------------------------
    // Independent oracles, transcribed from Sequencing_Artifact_Detection.md
    // (§2.2 substitution classes, §2.4 INV-01..INV-05, §4.2 decision table),
    // NOT routed through production constants. The substitution map, GIV ratio
    // and flag rule are recomputed here from literals so a self-consistent-but-
    // wrong production constant is still caught.
    //
    // Checklist-row remap notes (two row items are bogus and are NOT tested as
    // written):
    //  * "R: strand-bias ∈ [0,1]" is WRONG — FisherStrand FS is Phred-scaled and
    //    UNBOUNDED above (e.g. 108.384). The genuine contract (INV-03) is FS ≥ 0,
    //    which is what we assert; the [0,1] claim is dropped.
    //  * "M: stricter thresholds → ≤ survivors" does NOT map — FilterArtifacts has
    //    no tunable threshold (OxoG uses the fixed 1.5). It is remapped to GIV/damage
    //    monotonicity: for an OxoG substitution, once IsArtifact is true at some GIV
    //    (GIV > 1.5), it stays true at any higher GIV — raising the imbalance can
    //    never turn a flagged OxoG variant back into a survivor.
    // -------------------------------------------------------------------------

    /// <summary>Documented damaged-GIV threshold = 1.5 (literal, Nature Methods 2017 / Chen et al. 2017).</summary>
    private const double ArtifactDamagedGivThreshold = 1.5;

    /// <summary>
    /// INV-04 classification oracle: the substitution map on UPPER-cased (ref, alt).
    /// {C&gt;T, G&gt;A} ⇒ FfpeDeamination; {G&gt;T, C&gt;A} ⇒ OxoG; everything else ⇒ None. The two
    /// artifact classes are disjoint. Recomputed independently from §2.2 / §4.2.
    /// </summary>
    private static OncologyAnalyzer.ArtifactType ExpectedArtifactType(char reference, char alternate)
    {
        char r = char.ToUpperInvariant(reference);
        char a = char.ToUpperInvariant(alternate);
        return (r, a) switch
        {
            ('C', 'T') => OncologyAnalyzer.ArtifactType.FfpeDeamination,
            ('G', 'A') => OncologyAnalyzer.ArtifactType.FfpeDeamination,
            ('G', 'T') => OncologyAnalyzer.ArtifactType.OxoG,
            ('C', 'A') => OncologyAnalyzer.ArtifactType.OxoG,
            _ => OncologyAnalyzer.ArtifactType.None
        };
    }

    /// <summary>
    /// INV-02 GIV oracle: GIV = r2 == 0 ? (r1 == 0 ? 1.0 : +∞) : r1 / r2. Recomputed from §6.1 / §2.2.
    /// </summary>
    private static double ExpectedGiv(int r1, int r2) =>
        r2 == 0 ? (r1 == 0 ? 1.0 : double.PositiveInfinity) : (double)r1 / r2;

    /// <summary>
    /// INV-04 flag oracle (§4.2 decision table): FFPE ⇒ always flagged; OxoG ⇒ flagged iff GIV &gt; 1.5;
    /// None ⇒ never flagged. Recomputed independently from the documented rule.
    /// </summary>
    private static bool ExpectedIsArtifact(OncologyAnalyzer.ArtifactType type, double giv) => type switch
    {
        OncologyAnalyzer.ArtifactType.FfpeDeamination => true,
        OncologyAnalyzer.ArtifactType.OxoG => giv > ArtifactDamagedGivThreshold,
        _ => false
    };

    // -------------------------------------------------------------------------
    // Generators
    // -------------------------------------------------------------------------

    private static readonly char[] ArtifactBases = { 'A', 'C', 'G', 'T' };

    /// <summary>
    /// Generates a full <see cref="OncologyAnalyzer.ArtifactObservation"/> with ref/alt drawn from
    /// {A,C,G,T} (so all four artifact pairs AND non-artifact pairs are reached) and non-negative strand
    /// and read-mate counts. Read-mate counts span the GIV branches (r2 = 0 reachable).
    /// </summary>
    private static Arbitrary<OncologyAnalyzer.ArtifactObservation> ArtifactObservationArbitrary() =>
        (from refIdx in Gen.Choose(0, 3)
         from altIdx in Gen.Choose(0, 3)
         from refFwd in Gen.Choose(0, 400)
         from refRev in Gen.Choose(0, 400)
         from altFwd in Gen.Choose(0, 400)
         from altRev in Gen.Choose(0, 400)
         from r1 in Gen.Choose(0, 400)
         from r2 in Gen.Choose(0, 400)
         select new OncologyAnalyzer.ArtifactObservation(
             ArtifactBases[refIdx], ArtifactBases[altIdx],
             refFwd, refRev, altFwd, altRev, r1, r2))
        .ToArbitrary();

    /// <summary>Generates a list (0..12) of contract-valid observations, preserving order.</summary>
    private static Arbitrary<OncologyAnalyzer.ArtifactObservation[]> ArtifactListArbitrary() =>
        (from n in Gen.Choose(0, 12)
         from arr in ArtifactObservationArbitrary().Generator.ArrayOf(n)
         select arr)
        .ToArbitrary();

    /// <summary>Generates an (r1, r2) pair across all GIV branches, including r2 = 0.</summary>
    private static Arbitrary<(int r1, int r2)> GivPairArbitrary() =>
        (from r1 in Gen.Choose(0, 500)
         from r2 in Gen.Choose(0, 500)
         select (r1, r2))
        .ToArbitrary();

    /// <summary>Generates a balanced strand table (refFwd == refRev, altFwd == altRev) ⇒ FS = 0.</summary>
    private static Arbitrary<(int refSym, int altSym)> BalancedTableArbitrary() =>
        (from refSym in Gen.Choose(0, 300)
         from altSym in Gen.Choose(0, 300)
         select (refSym, altSym))
        .ToArbitrary();

    /// <summary>Generates an arbitrary non-negative strand table for the FS ≥ 0 universal property.</summary>
    private static Arbitrary<(int refFwd, int refRev, int altFwd, int altRev)> StrandTableArbitrary() =>
        (from refFwd in Gen.Choose(0, 300)
         from refRev in Gen.Choose(0, 300)
         from altFwd in Gen.Choose(0, 300)
         from altRev in Gen.Choose(0, 300)
         select (refFwd, refRev, altFwd, altRev))
        .ToArbitrary();

    // -------------------------------------------------------------------------
    // Classification (INV-04) — substitution map + flag rule
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-04: <c>ClassifyArtifact</c>.Type matches the independent substitution map EXACTLY (the two
    /// artifact classes disjoint, None for everything else), and <c>IsArtifact</c> matches the §4.2 flag
    /// rule (FFPE always; OxoG iff GIV &gt; 1.5; None never), driven over random ref/alt/counts.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_Classification_MatchesSubstitutionAndFlagOracle()
    {
        return Prop.ForAll(ArtifactObservationArbitrary(), o =>
        {
            var call = OncologyAnalyzer.ClassifyArtifact(o);
            var expectedType = ExpectedArtifactType(o.ReferenceAllele, o.AlternateAllele);
            double expectedGiv = ExpectedGiv(o.AltReadsR1, o.AltReadsR2);
            bool expectedFlag = ExpectedIsArtifact(expectedType, expectedGiv);
            return (call.Type == expectedType && call.IsArtifact == expectedFlag)
                .Label($"{o.ReferenceAllele}>{o.AlternateAllele}: type {call.Type}/{expectedType}, " +
                       $"flag {call.IsArtifact}/{expectedFlag} (giv={expectedGiv})");
        });
    }

    /// <summary>
    /// INV-04 (case-insensitivity, §3.3): lower-cased bases classify identically to their upper-cased
    /// counterparts — same Type and IsArtifact.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_Classification_IsCaseInsensitive()
    {
        return Prop.ForAll(ArtifactObservationArbitrary(), o =>
        {
            var upper = OncologyAnalyzer.ClassifyArtifact(o);
            var lower = OncologyAnalyzer.ClassifyArtifact(o with
            {
                ReferenceAllele = char.ToLowerInvariant(o.ReferenceAllele),
                AlternateAllele = char.ToLowerInvariant(o.AlternateAllele)
            });
            return (upper.Type == lower.Type && upper.IsArtifact == lower.IsArtifact)
                .Label($"upper {upper.Type}/{upper.IsArtifact} vs lower {lower.Type}/{lower.IsArtifact}");
        });
    }

    // -------------------------------------------------------------------------
    // GIV (INV-02) — ratio oracle + non-negativity
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-02: <c>CalculateGivScore</c> equals the independently recomputed ratio (1.0 when both 0,
    /// +∞ when only r2 = 0, else r1/r2), is always ≥ 0, and equal positive counts ⇒ exactly 1.0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_Giv_MatchesRatioOracle()
    {
        return Prop.ForAll(GivPairArbitrary(), p =>
        {
            double giv = OncologyAnalyzer.CalculateGivScore(p.r1, p.r2);
            double expected = ExpectedGiv(p.r1, p.r2);
            bool matches = double.IsPositiveInfinity(expected)
                ? double.IsPositiveInfinity(giv)
                : Math.Abs(giv - expected) < 1e-12;
            bool nonNegative = giv >= 0.0;
            bool equalCountsOne = p.r1 != p.r2 || p.r1 == 0 || Math.Abs(giv - 1.0) < 1e-12;
            return (matches && nonNegative && equalCountsOne)
                .Label($"giv {giv} vs {expected} (r1={p.r1}, r2={p.r2})");
        });
    }

    // -------------------------------------------------------------------------
    // FisherStrand FS (INV-03, INV-05) — FS ≥ 0, balanced ⇒ 0, monotone segregation
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-03: <c>CalculateStrandBias</c> (FS) is ALWAYS ≥ 0 over arbitrary non-negative tables (Phred of
    /// a p ≤ 1). NOTE: this is the genuine contract; the bogus checklist claim "strand-bias ∈ [0,1]" is
    /// dropped — FS is unbounded above.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_StrandBias_IsNonNegative()
    {
        return Prop.ForAll(StrandTableArbitrary(), t =>
        {
            double fs = OncologyAnalyzer.CalculateStrandBias(t.refFwd, t.refRev, t.altFwd, t.altRev);
            return (fs >= 0.0)
                .Label($"FS={fs} for [{t.refFwd},{t.refRev},{t.altFwd},{t.altRev}]");
        });
    }

    /// <summary>
    /// INV-03: a symmetric/balanced table (refFwd == refRev AND altFwd == altRev, incl. all-zero) has no
    /// strand bias (p = 1) ⇒ FS = 0 (within 1e-9).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_StrandBias_BalancedTableIsZero()
    {
        return Prop.ForAll(BalancedTableArbitrary(), t =>
        {
            double fs = OncologyAnalyzer.CalculateStrandBias(t.refSym, t.refSym, t.altSym, t.altSym);
            return (Math.Abs(fs) < 1e-9)
                .Label($"FS={fs} for balanced [{t.refSym},{t.refSym},{t.altSym},{t.altSym}]");
        });
    }

    /// <summary>
    /// INV-05: greater strand segregation ⇒ FS non-decreasing. Family of fully-segregated tables
    /// [n, 0, 0, n] for growing n; the Fisher p decreases (more extreme), so FS must be non-decreasing.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_StrandBias_MonotoneInSegregation()
    {
        return Prop.ForAll(Gen.Choose(1, 60).ToArbitrary(), n =>
        {
            double fsLo = OncologyAnalyzer.CalculateStrandBias(n, 0, 0, n);
            double fsHi = OncologyAnalyzer.CalculateStrandBias(n + 1, 0, 0, n + 1);
            return (fsHi >= fsLo - 1e-9)
                .Label($"FS({n})={fsLo} > FS({n + 1})={fsHi}");
        });
    }

    // -------------------------------------------------------------------------
    // M (remapped) — GIV/damage monotonicity for OxoG
    // -------------------------------------------------------------------------

    /// <summary>
    /// M (remapped from the bogus "stricter thresholds → ≤ survivors"): for a fixed OxoG substitution,
    /// raising the GIV (more read-1 imbalance) can never un-flag the variant — if it is flagged at a lower
    /// GIV (GIV &gt; 1.5) it stays flagged at any higher GIV. We hold r2 fixed and let r1 grow.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_OxoG_FlagMonotoneInGiv()
    {
        var gen = (from r2 in Gen.Choose(1, 200)
                   from r1Lo in Gen.Choose(0, 600)
                   from delta in Gen.Choose(0, 600)
                   select (r2, r1Lo, r1Hi: r1Lo + delta)).ToArbitrary();
        return Prop.ForAll(gen, t =>
        {
            // G>T is the canonical OxoG substitution.
            var lo = new OncologyAnalyzer.ArtifactObservation('G', 'T', 10, 10, 5, 5, t.r1Lo, t.r2);
            var hi = lo with { AltReadsR1 = t.r1Hi };
            var callLo = OncologyAnalyzer.ClassifyArtifact(lo);
            var callHi = OncologyAnalyzer.ClassifyArtifact(hi);
            // Higher GIV (r1Hi ≥ r1Lo) ⇒ if flagged at lo, still flagged at hi.
            bool monotone = !callLo.IsArtifact || callHi.IsArtifact;
            return monotone
                .Label($"r2={t.r2}: lo(r1={t.r1Lo}) flag={callLo.IsArtifact}, hi(r1={t.r1Hi}) flag={callHi.IsArtifact}");
        });
    }

    // -------------------------------------------------------------------------
    // P (INV-01) — FilterArtifacts is a non-flagged subset in input order
    // -------------------------------------------------------------------------

    /// <summary>
    /// P / INV-01: <c>FilterArtifacts</c> returns EXACTLY the non-flagged inputs in input order — i.e. it
    /// equals <c>input.Where(x =&gt; !ClassifyArtifact(x).IsArtifact)</c>. The result is a subset of the input,
    /// preserves order, and every survivor is non-flagged.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_Filter_IsNonFlaggedSubsetInOrder()
    {
        return Prop.ForAll(ArtifactListArbitrary(), variants =>
        {
            var survivors = OncologyAnalyzer.FilterArtifacts(variants);
            var expected = variants.Where(x => !OncologyAnalyzer.ClassifyArtifact(x).IsArtifact).ToArray();
            bool matchesOracle = survivors.SequenceEqual(expected);
            bool allNonFlagged = survivors.All(x => !OncologyAnalyzer.ClassifyArtifact(x).IsArtifact);
            return (matchesOracle && allNonFlagged)
                .Label($"survivors={survivors.Count}, expected={expected.Length}, allNonFlagged={allNonFlagged}");
        });
    }

    // -------------------------------------------------------------------------
    // DetectOxoGArtifacts — exactly the OxoG calls with GIV > 1.5
    // -------------------------------------------------------------------------

    /// <summary>
    /// <c>DetectOxoGArtifacts</c> returns exactly the calls whose Type is OxoG AND GIV &gt; 1.5 (the OxoG
    /// flag rule), in input order, recomputed independently from the substitution map and GIV oracle.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_DetectOxoG_ReturnsExactlyDamagedOxoG()
    {
        return Prop.ForAll(ArtifactListArbitrary(), variants =>
        {
            var detected = OncologyAnalyzer.DetectOxoGArtifacts(variants);
            var expected = variants
                .Where(o => ExpectedArtifactType(o.ReferenceAllele, o.AlternateAllele)
                                == OncologyAnalyzer.ArtifactType.OxoG
                            && ExpectedGiv(o.AltReadsR1, o.AltReadsR2) > ArtifactDamagedGivThreshold)
                .Select(o => OncologyAnalyzer.ClassifyArtifact(o))
                .ToArray();
            return detected.SequenceEqual(expected)
                .Label($"detected={detected.Count}, expected={expected.Length}");
        });
    }

    // -------------------------------------------------------------------------
    // D (determinism)
    // -------------------------------------------------------------------------

    /// <summary>
    /// D: identical inputs ⇒ identical results — <c>ClassifyArtifact</c> returns an equal <c>ArtifactCall</c>
    /// (all fields; +∞ == +∞ via record equality) across two independent calls.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_Classification_IsDeterministic()
    {
        return Prop.ForAll(ArtifactObservationArbitrary(), o =>
        {
            var a = OncologyAnalyzer.ClassifyArtifact(o);
            var b = OncologyAnalyzer.ClassifyArtifact(o);
            return a.Equals(b).Label($"{a} != {b}");
        });
    }

    /// <summary>D: <c>FilterArtifacts</c> is deterministic — two runs over the same list yield equal results.</summary>
    [FsCheck.NUnit.Property]
    public Property Artifact_Filter_IsDeterministic()
    {
        return Prop.ForAll(ArtifactListArbitrary(), variants =>
        {
            var a = OncologyAnalyzer.FilterArtifacts(variants);
            var b = OncologyAnalyzer.FilterArtifacts(variants);
            return a.SequenceEqual(b).Label($"a={a.Count}, b={b.Count}");
        });
    }

    // -------------------------------------------------------------------------
    // Edge cases and worked-example anchors (§6.1, §7.1)
    // -------------------------------------------------------------------------

    /// <summary>§6.1: GIV with r2 = 0, r1 = 0 ⇒ 1.0 (no imbalance evidence).</summary>
    [Test]
    [Category("Property")]
    public void Artifact_Giv_BothZero_IsOne()
    {
        Assert.That(OncologyAnalyzer.CalculateGivScore(0, 0), Is.EqualTo(1.0));
    }

    /// <summary>§6.1: GIV with r2 = 0, r1 &gt; 0 ⇒ +∞ (maximal one-sided imbalance).</summary>
    [Test]
    [Category("Property")]
    public void Artifact_Giv_OnlyR2Zero_IsPositiveInfinity()
    {
        Assert.That(OncologyAnalyzer.CalculateGivScore(7, 0), Is.EqualTo(double.PositiveInfinity));
    }

    /// <summary>§6.1: equal positive counts ⇒ GIV = 1.0 (balanced library).</summary>
    [Test]
    [Category("Property")]
    public void Artifact_Giv_EqualPositiveCounts_IsOne()
    {
        Assert.That(OncologyAnalyzer.CalculateGivScore(50, 50), Is.EqualTo(1.0).Within(1e-12));
    }

    /// <summary>§6.1: a balanced strand table ⇒ FS = 0 (p = 1).</summary>
    [Test]
    [Category("Property")]
    public void Artifact_StrandBias_BalancedAnchor_IsZero()
    {
        Assert.That(OncologyAnalyzer.CalculateStrandBias(20, 20, 10, 10), Is.EqualTo(0.0).Within(1e-9));
    }

    /// <summary>§7.1 numerical walk-through: table [20,0,0,20] ⇒ FS ≈ 108.384 (within 1e-2).</summary>
    [Test]
    [Category("Property")]
    public void Artifact_StrandBias_WorkedExampleAnchor()
    {
        double fs = OncologyAnalyzer.CalculateStrandBias(20, 0, 0, 20);
        Assert.That(fs, Is.EqualTo(108.384).Within(1e-2));
    }

    /// <summary>§6.1: a non-artifact substitution (A&gt;G) ⇒ Type None, never flagged.</summary>
    [Test]
    [Category("Property")]
    public void Artifact_NonArtifactSubstitution_IsNoneNotFlagged()
    {
        var o = new OncologyAnalyzer.ArtifactObservation('A', 'G', 10, 10, 10, 10, 100, 1);
        var call = OncologyAnalyzer.ClassifyArtifact(o);
        Assert.Multiple(() =>
        {
            Assert.That(call.Type, Is.EqualTo(OncologyAnalyzer.ArtifactType.None));
            Assert.That(call.IsArtifact, Is.False);
        });
    }

    /// <summary>
    /// §7.1 worked example: G&gt;T with r1 = 200, r2 = 100 ⇒ GIV = 2.0 (&gt; 1.5), Type OxoG, IsArtifact true.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Artifact_OxoGWorkedExample_Anchor()
    {
        var v = new OncologyAnalyzer.ArtifactObservation('G', 'T', 20, 18, 12, 10, 200, 100);
        var call = OncologyAnalyzer.ClassifyArtifact(v);
        Assert.Multiple(() =>
        {
            Assert.That(call.GivScore, Is.EqualTo(2.0).Within(1e-12));
            Assert.That(call.Type, Is.EqualTo(OncologyAnalyzer.ArtifactType.OxoG));
            Assert.That(call.IsArtifact, Is.True);
        });
    }

    /// <summary>INV-01: an empty variant list ⇒ empty filter result.</summary>
    [Test]
    [Category("Property")]
    public void Artifact_Filter_EmptyList_IsEmpty()
    {
        Assert.That(OncologyAnalyzer.FilterArtifacts(Array.Empty<OncologyAnalyzer.ArtifactObservation>()),
            Is.Empty);
    }

    /// <summary>§6.1 API contract: a null variant list ⇒ <c>ArgumentNullException</c>.</summary>
    [Test]
    [Category("Property")]
    public void Artifact_Filter_NullList_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.FilterArtifacts(null!));
    }

    /// <summary>§6.1 API contract: <c>DetectOxoGArtifacts</c> with a null list ⇒ <c>ArgumentNullException</c>.</summary>
    [Test]
    [Category("Property")]
    public void Artifact_DetectOxoG_NullList_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.DetectOxoGArtifacts(null!));
    }

    /// <summary>§3.3 validation: a negative GIV count ⇒ <c>ArgumentOutOfRangeException</c>.</summary>
    [Test]
    [Category("Property")]
    public void Artifact_Giv_NegativeCount_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateGivScore(-1, 5));
    }

    /// <summary>§3.3 validation: a negative strand-table cell ⇒ <c>ArgumentOutOfRangeException</c>.</summary>
    [Test]
    [Category("Property")]
    public void Artifact_StrandBias_NegativeCount_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateStrandBias(-1, 0, 0, 0));
    }

    #endregion

    #region ONCO-ANNOT-001 — Cancer Variant Annotation (AMP/ASCO/CAP 2017 Tiers)

    // -------------------------------------------------------------------------
    // Independent theory oracle — the AMP/ASCO/CAP 2017 Level→Tier cascade
    // (Li MM et al. 2017, J Mol Diagn 19(1):4–23, Figure 2 / Tables 6–7).
    //
    // NOTE on the checklist row: the loose row text reads "R: consequence ∈ SO
    // vocabulary", but Cancer_Variant_Annotation.md classifies into the four-member
    // AMP/ASCO/CAP VariantTier vocabulary, not Sequence-Ontology consequences. We
    // therefore map "R" to the documented contract: tier ∈ {TierI,TierII,TierIII,
    // TierIV} (INV-01, total/exhaustive). The remaining rows (P: annotation preserves
    // the variant; D: deterministic) are honoured verbatim.
    //
    // The 0.01 benign cutoff is a LOCAL literal here (not routed through the
    // production constant) so a wrong production constant is still caught.
    // -------------------------------------------------------------------------

    /// <summary>
    /// AMP/ASCO/CAP 2017 benign population-frequency cutoff (1%), transcribed LITERALLY from the
    /// guideline ("the work group recommends using 1% (0.01) as a primary cutoff"). Kept as a local
    /// literal so a wrong production constant is not silently mirrored into the oracle.
    /// </summary>
    private const double AnnotBenignMafCutoff = 0.01;

    /// <summary>
    /// Independent Level→Tier cascade transcribed from Cancer_Variant_Annotation.md §4.1 / §4.2
    /// (Li et al. 2017, Figure 2): Level A or B ⇒ Tier I; Level C or D ⇒ Tier II; otherwise (no level)
    /// MAF ≥ 0.01 OR no cancer association ⇒ Tier IV; otherwise Tier III. Recomputed without calling the
    /// production classifier, so a wrong production rule is caught.
    /// </summary>
    private static OncologyAnalyzer.VariantTier ExpectedTier(
        OncologyAnalyzer.ClinicalEvidenceLevel level, double maf, bool hasCancerAssociation)
    {
        if (level is OncologyAnalyzer.ClinicalEvidenceLevel.A or OncologyAnalyzer.ClinicalEvidenceLevel.B)
        {
            return OncologyAnalyzer.VariantTier.TierI_StrongClinicalSignificance;
        }

        if (level is OncologyAnalyzer.ClinicalEvidenceLevel.C or OncologyAnalyzer.ClinicalEvidenceLevel.D)
        {
            return OncologyAnalyzer.VariantTier.TierII_PotentialClinicalSignificance;
        }

        if (maf >= AnnotBenignMafCutoff || !hasCancerAssociation)
        {
            return OncologyAnalyzer.VariantTier.TierIV_BenignOrLikelyBenign;
        }

        return OncologyAnalyzer.VariantTier.TierIII_UnknownClinicalSignificance;
    }

    /// <summary>The four members of the controlled <see cref="OncologyAnalyzer.VariantTier"/> vocabulary.</summary>
    private static readonly OncologyAnalyzer.VariantTier[] AnnotTierVocabulary =
    {
        OncologyAnalyzer.VariantTier.TierI_StrongClinicalSignificance,
        OncologyAnalyzer.VariantTier.TierII_PotentialClinicalSignificance,
        OncologyAnalyzer.VariantTier.TierIII_UnknownClinicalSignificance,
        OncologyAnalyzer.VariantTier.TierIV_BenignOrLikelyBenign,
    };

    // -------------------------------------------------------------------------
    // Generators
    // -------------------------------------------------------------------------

    /// <summary>All five evidence levels, including None, so every cascade branch is reachable.</summary>
    private static readonly OncologyAnalyzer.ClinicalEvidenceLevel[] AnnotEvidenceLevels =
    {
        OncologyAnalyzer.ClinicalEvidenceLevel.None,
        OncologyAnalyzer.ClinicalEvidenceLevel.A,
        OncologyAnalyzer.ClinicalEvidenceLevel.B,
        OncologyAnalyzer.ClinicalEvidenceLevel.C,
        OncologyAnalyzer.ClinicalEvidenceLevel.D,
    };

    /// <summary>Small gene-symbol pool for variant generation and COSMIC-key collisions.</summary>
    private static readonly string[] AnnotGenes = { "BRAF", "KRAS", "EGFR", "TP53", "PIK3CA" };

    /// <summary>Small protein-change pool (HGVS p.-notation).</summary>
    private static readonly string[] AnnotProteinChanges = { "p.V600E", "p.G12D", "p.L858R", "p.R175H" };

    /// <summary>
    /// Generates a contract-valid MAF in [0, 1] with deliberate emphasis on the 0.01 cutoff boundary:
    /// the exact cutoff (0.01), just below (0.0099), 0.0, 1.0, and a dense band of values straddling
    /// the threshold, plus a uniform sweep over the full range.
    /// </summary>
    private static Gen<double> AnnotMafGen() =>
        Gen.Frequency(
            (1, Gen.Constant(0.0)),
            (1, Gen.Constant(1.0)),
            (1, Gen.Constant(AnnotBenignMafCutoff)),                // exactly 0.01 ⇒ Tier IV
            (1, Gen.Constant(0.0099)),                              // just below ⇒ III/IV by assoc
            (2, Gen.Choose(0, 200).Select(p => p / 10000.0)),       // dense band [0, 0.02]
            (2, Gen.Choose(0, 10000).Select(p => p / 10000.0)));    // uniform [0, 1]

    /// <summary>Generates a fully-populated, contract-valid variant spanning all branches of the cascade.</summary>
    private static Arbitrary<OncologyAnalyzer.CancerVariantAnnotationInput> AnnotVariantArbitrary() =>
        (from gene in Gen.Elements(AnnotGenes)
         from change in Gen.Elements(AnnotProteinChanges)
         from level in Gen.Elements(AnnotEvidenceLevels)
         from maf in AnnotMafGen()
         from assoc in Gen.Elements(true, false)
         select new OncologyAnalyzer.CancerVariantAnnotationInput(gene, change, level, maf, assoc))
        .ToArbitrary();

    /// <summary>Generates a batch (0..12) of contract-valid variants, order-preserving.</summary>
    private static Arbitrary<OncologyAnalyzer.CancerVariantAnnotationInput[]> AnnotVariantListArbitrary() =>
        (from n in Gen.Choose(0, 12)
         from arr in AnnotVariantArbitrary().Generator.ArrayOf(n)
         select arr)
        .ToArbitrary();

    /// <summary>
    /// Generates a variant restricted to evidence level A or B (so the oracle expects Tier I regardless
    /// of MAF / association), paired with a random MAF and association flag to prove INV-02 independence.
    /// </summary>
    private static Arbitrary<OncologyAnalyzer.CancerVariantAnnotationInput> AnnotTierIVariantArbitrary() =>
        (from gene in Gen.Elements(AnnotGenes)
         from change in Gen.Elements(AnnotProteinChanges)
         from level in Gen.Elements(
             OncologyAnalyzer.ClinicalEvidenceLevel.A, OncologyAnalyzer.ClinicalEvidenceLevel.B)
         from maf in AnnotMafGen()
         from assoc in Gen.Elements(true, false)
         select new OncologyAnalyzer.CancerVariantAnnotationInput(gene, change, level, maf, assoc))
        .ToArbitrary();

    /// <summary>Generates a variant restricted to evidence level C or D (oracle: Tier II regardless of MAF/assoc).</summary>
    private static Arbitrary<OncologyAnalyzer.CancerVariantAnnotationInput> AnnotTierIIVariantArbitrary() =>
        (from gene in Gen.Elements(AnnotGenes)
         from change in Gen.Elements(AnnotProteinChanges)
         from level in Gen.Elements(
             OncologyAnalyzer.ClinicalEvidenceLevel.C, OncologyAnalyzer.ClinicalEvidenceLevel.D)
         from maf in AnnotMafGen()
         from assoc in Gen.Elements(true, false)
         select new OncologyAnalyzer.CancerVariantAnnotationInput(gene, change, level, maf, assoc))
        .ToArbitrary();

    // -------------------------------------------------------------------------
    // Tier classification (INV-01, INV-02, INV-03 / R: tier ∈ vocabulary)
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-01/INV-03 (R: tier ∈ the four-member VariantTier vocabulary): <c>ClassifyVariantTier</c> matches
    /// the independent Level→Tier cascade EXACTLY over generators spanning every evidence level (incl. None),
    /// MAF ∈ [0,1] with boundary emphasis (0.01 / 0.0099), and both association flags. The cascade is total,
    /// so this also proves the result is always one of the four tiers.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annot_ClassifyTier_MatchesCascadeOracle()
    {
        return Prop.ForAll(AnnotVariantArbitrary(), v =>
        {
            var tier = OncologyAnalyzer.ClassifyVariantTier(v);
            var expected = ExpectedTier(v.EvidenceLevel, v.PopulationMaf, v.HasCancerAssociation);
            return (tier == expected)
                .Label($"level={v.EvidenceLevel}, maf={v.PopulationMaf}, assoc={v.HasCancerAssociation}: " +
                       $"got {tier}, expected {expected}");
        });
    }

    /// <summary>
    /// INV-01 (exhaustive/total): for any contract-valid variant, <c>ClassifyVariantTier</c> returns a value
    /// that is a member of the four-element <see cref="OncologyAnalyzer.VariantTier"/> vocabulary.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annot_ClassifyTier_AlwaysInControlledVocabulary()
    {
        return Prop.ForAll(AnnotVariantArbitrary(), v =>
        {
            var tier = OncologyAnalyzer.ClassifyVariantTier(v);
            return AnnotTierVocabulary.Contains(tier).Label($"tier {tier} not in VariantTier vocabulary");
        });
    }

    /// <summary>
    /// INV-02 (Level A/B ⇒ Tier I, independent of MAF / association): driven with RANDOM MAF and association,
    /// a level-A/B variant is always Tier I — proving the tier depends only on the evidence level here.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annot_LevelAB_AlwaysTierI_RegardlessOfMafAndAssociation()
    {
        return Prop.ForAll(AnnotTierIVariantArbitrary(), v =>
        {
            var tier = OncologyAnalyzer.ClassifyVariantTier(v);
            return (tier == OncologyAnalyzer.VariantTier.TierI_StrongClinicalSignificance)
                .Label($"level={v.EvidenceLevel}, maf={v.PopulationMaf}, assoc={v.HasCancerAssociation}: got {tier}");
        });
    }

    /// <summary>
    /// INV-02 (Level C/D ⇒ Tier II, independent of MAF / association): driven with RANDOM MAF and association,
    /// a level-C/D variant is always Tier II.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annot_LevelCD_AlwaysTierII_RegardlessOfMafAndAssociation()
    {
        return Prop.ForAll(AnnotTierIIVariantArbitrary(), v =>
        {
            var tier = OncologyAnalyzer.ClassifyVariantTier(v);
            return (tier == OncologyAnalyzer.VariantTier.TierII_PotentialClinicalSignificance)
                .Label($"level={v.EvidenceLevel}, maf={v.PopulationMaf}, assoc={v.HasCancerAssociation}: got {tier}");
        });
    }

    // -------------------------------------------------------------------------
    // Batch annotation (P: annotation preserves the variant + input order)
    // -------------------------------------------------------------------------

    /// <summary>
    /// P (annotation preserves the variant / order, §3.2): <c>AnnotateCancerVariants</c> returns exactly one
    /// annotation per input, in input order, each preserving the input variant VERBATIM (Gene, ProteinChange,
    /// EvidenceLevel, PopulationMaf, HasCancerAssociation) with <c>Tier == ClassifyVariantTier(thatInput)</c>.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annot_Annotate_PreservesVariantAndOrderWithMatchingTier()
    {
        return Prop.ForAll(AnnotVariantListArbitrary(), variants =>
        {
            var annotations = OncologyAnalyzer.AnnotateCancerVariants(variants);

            bool countOk = annotations.Count == variants.Length;
            bool allOk = countOk && variants
                .Select((v, i) => annotations[i].Variant.Equals(v)
                    && annotations[i].Tier == OncologyAnalyzer.ClassifyVariantTier(v))
                .All(ok => ok);

            return allOk.Label($"count={annotations.Count}/{variants.Length}, perItem preserved+tier-matched={allOk}");
        });
    }

    // -------------------------------------------------------------------------
    // COSMIC lookup (exact ordinal (gene, proteinChange) key)
    // -------------------------------------------------------------------------

    /// <summary>
    /// §3.3 / §5.2: <c>GetCOSMICAnnotation</c> returns the catalog value on an exact (Gene, ProteinChange)
    /// hit and <c>null</c> on a miss. Built by seeding the catalog with the variant's own key, so a hit is
    /// guaranteed and the returned value equals the seeded id.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annot_Cosmic_ReturnsSeededValueOnExactHit()
    {
        return Prop.ForAll(AnnotVariantArbitrary(), v =>
        {
            var catalog = new Dictionary<(string Gene, string ProteinChange), string>
            {
                [(v.Gene, v.ProteinChange)] = "COSV-SEED",
            };
            return (OncologyAnalyzer.GetCOSMICAnnotation(v, catalog) == "COSV-SEED")
                .Label($"hit on ({v.Gene},{v.ProteinChange})");
        });
    }

    /// <summary>
    /// §6.1 (catalog miss ⇒ null): against an EMPTY catalog, <c>GetCOSMICAnnotation</c> always returns null.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annot_Cosmic_ReturnsNullOnMiss()
    {
        var empty = new Dictionary<(string Gene, string ProteinChange), string>();
        return Prop.ForAll(AnnotVariantArbitrary(), v =>
            (OncologyAnalyzer.GetCOSMICAnnotation(v, empty) == null)
                .Label($"expected null for ({v.Gene},{v.ProteinChange})"));
    }

    /// <summary>
    /// §3.3 ordinal (case-sensitive) key equality: a catalog keyed by an UPPER-cased gene does not match a
    /// variant whose gene differs only in case ⇒ <c>null</c> (no case-folding).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Annot_Cosmic_IsOrdinalCaseSensitive()
    {
        // Genes in the pool are upper-case; a lower-cased gene is a distinct ordinal key.
        return Prop.ForAll(AnnotVariantArbitrary(), v =>
        {
            string lowered = v.Gene.ToLowerInvariant();
            return Prop.When(!string.Equals(lowered, v.Gene, StringComparison.Ordinal), () =>
            {
                var catalog = new Dictionary<(string Gene, string ProteinChange), string>
                {
                    [(v.Gene, v.ProteinChange)] = "COSV-SEED",
                };
                var lowerVariant = v with { Gene = lowered };
                return (OncologyAnalyzer.GetCOSMICAnnotation(lowerVariant, catalog) == null)
                    .Label($"case-differing gene '{lowered}' vs key '{v.Gene}' should miss");
            });
        });
    }

    // -------------------------------------------------------------------------
    // D — determinism
    // -------------------------------------------------------------------------

    /// <summary>D (determinism): identical input ⇒ identical tier, annotation list, and COSMIC result.</summary>
    [FsCheck.NUnit.Property]
    public Property Annot_Deterministic()
    {
        return Prop.ForAll(AnnotVariantListArbitrary(), variants =>
        {
            var tiers1 = variants.Select(OncologyAnalyzer.ClassifyVariantTier).ToArray();
            var tiers2 = variants.Select(OncologyAnalyzer.ClassifyVariantTier).ToArray();

            var ann1 = OncologyAnalyzer.AnnotateCancerVariants(variants);
            var ann2 = OncologyAnalyzer.AnnotateCancerVariants(variants);

            var catalog = new Dictionary<(string Gene, string ProteinChange), string>
            {
                [("BRAF", "p.V600E")] = "COSV56056643",
            };
            var cos1 = variants.Select(v => OncologyAnalyzer.GetCOSMICAnnotation(v, catalog)).ToArray();
            var cos2 = variants.Select(v => OncologyAnalyzer.GetCOSMICAnnotation(v, catalog)).ToArray();

            return (tiers1.SequenceEqual(tiers2)
                    && ann1.SequenceEqual(ann2)
                    && cos1.SequenceEqual(cos2))
                .Label($"n={variants.Length}: tiers/annotations/cosmic stable across repeats");
        });
    }

    // -------------------------------------------------------------------------
    // Validation and worked-example anchors (§3.3, §6.1, §7.1)
    // -------------------------------------------------------------------------

    /// <summary>§3.3: PopulationMaf NaN ⇒ <c>ArgumentOutOfRangeException</c>.</summary>
    [Test]
    [Category("Property")]
    public void Annot_Maf_NaN_Throws()
    {
        var v = new OncologyAnalyzer.CancerVariantAnnotationInput(
            "BRAF", "p.V600E", OncologyAnalyzer.ClinicalEvidenceLevel.None, double.NaN, true);
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyVariantTier(v));
    }

    /// <summary>§3.3: PopulationMaf &lt; 0 ⇒ <c>ArgumentOutOfRangeException</c>.</summary>
    [Test]
    [Category("Property")]
    public void Annot_Maf_Negative_Throws()
    {
        var v = new OncologyAnalyzer.CancerVariantAnnotationInput(
            "BRAF", "p.V600E", OncologyAnalyzer.ClinicalEvidenceLevel.None, -0.01, true);
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyVariantTier(v));
    }

    /// <summary>§3.3: PopulationMaf &gt; 1 ⇒ <c>ArgumentOutOfRangeException</c>.</summary>
    [Test]
    [Category("Property")]
    public void Annot_Maf_AboveOne_Throws()
    {
        var v = new OncologyAnalyzer.CancerVariantAnnotationInput(
            "BRAF", "p.V600E", OncologyAnalyzer.ClinicalEvidenceLevel.None, 1.5, true);
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyVariantTier(v));
    }

    /// <summary>§3.3: a null variants batch ⇒ <c>ArgumentNullException</c>.</summary>
    [Test]
    [Category("Property")]
    public void Annot_Annotate_NullBatch_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.AnnotateCancerVariants(null!));
    }

    /// <summary>§3.3: an empty batch ⇒ an empty annotation list.</summary>
    [Test]
    [Category("Property")]
    public void Annot_Annotate_EmptyBatch_IsEmpty()
    {
        Assert.That(
            OncologyAnalyzer.AnnotateCancerVariants(
                Array.Empty<OncologyAnalyzer.CancerVariantAnnotationInput>()),
            Is.Empty);
    }

    /// <summary>§3.3: a null COSMIC catalog ⇒ <c>ArgumentNullException</c>.</summary>
    [Test]
    [Category("Property")]
    public void Annot_Cosmic_NullCatalog_Throws()
    {
        var v = new OncologyAnalyzer.CancerVariantAnnotationInput(
            "BRAF", "p.V600E", OncologyAnalyzer.ClinicalEvidenceLevel.A, 0.0, true);
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.GetCOSMICAnnotation(v, null!));
    }

    /// <summary>§6.1 anchor: Level A but high MAF ⇒ Tier I (categorized by evidence level).</summary>
    [Test]
    [Category("Property")]
    public void Annot_LevelA_HighMaf_IsTierI()
    {
        var v = new OncologyAnalyzer.CancerVariantAnnotationInput(
            "BRAF", "p.V600E", OncologyAnalyzer.ClinicalEvidenceLevel.A, 0.95, false);
        Assert.That(OncologyAnalyzer.ClassifyVariantTier(v),
            Is.EqualTo(OncologyAnalyzer.VariantTier.TierI_StrongClinicalSignificance));
    }

    /// <summary>§6.1 anchor: no level, MAF exactly 0.01 ⇒ Tier IV (cutoff is inclusive).</summary>
    [Test]
    [Category("Property")]
    public void Annot_NoLevel_MafExactlyCutoff_IsTierIV()
    {
        var v = new OncologyAnalyzer.CancerVariantAnnotationInput(
            "TP53", "p.R175H", OncologyAnalyzer.ClinicalEvidenceLevel.None, 0.01, true);
        Assert.That(OncologyAnalyzer.ClassifyVariantTier(v),
            Is.EqualTo(OncologyAnalyzer.VariantTier.TierIV_BenignOrLikelyBenign));
    }

    /// <summary>§6.1 anchor: no level, MAF 0.0099, cancer association ⇒ Tier III (below cutoff, associated).</summary>
    [Test]
    [Category("Property")]
    public void Annot_NoLevel_JustBelowCutoff_WithAssociation_IsTierIII()
    {
        var v = new OncologyAnalyzer.CancerVariantAnnotationInput(
            "TP53", "p.R175H", OncologyAnalyzer.ClinicalEvidenceLevel.None, 0.0099, true);
        Assert.That(OncologyAnalyzer.ClassifyVariantTier(v),
            Is.EqualTo(OncologyAnalyzer.VariantTier.TierIII_UnknownClinicalSignificance));
    }

    /// <summary>§6.1 anchor: no level, low MAF, no cancer association ⇒ Tier IV.</summary>
    [Test]
    [Category("Property")]
    public void Annot_NoLevel_LowMaf_NoAssociation_IsTierIV()
    {
        var v = new OncologyAnalyzer.CancerVariantAnnotationInput(
            "TP53", "p.R175H", OncologyAnalyzer.ClinicalEvidenceLevel.None, 0.0, false);
        Assert.That(OncologyAnalyzer.ClassifyVariantTier(v),
            Is.EqualTo(OncologyAnalyzer.VariantTier.TierIV_BenignOrLikelyBenign));
    }

    /// <summary>
    /// §7.1 worked example: BRAF p.V600E Level A ⇒ Tier I, and the COSMIC lookup returns the catalog id.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Annot_BrafV600E_WorkedExample_Anchor()
    {
        var v = new OncologyAnalyzer.CancerVariantAnnotationInput(
            "BRAF", "p.V600E", OncologyAnalyzer.ClinicalEvidenceLevel.A, 0.0, true);
        var catalog = new Dictionary<(string Gene, string ProteinChange), string>
        {
            [("BRAF", "p.V600E")] = "COSV56056643",
        };
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.ClassifyVariantTier(v),
                Is.EqualTo(OncologyAnalyzer.VariantTier.TierI_StrongClinicalSignificance));
            Assert.That(OncologyAnalyzer.GetCOSMICAnnotation(v, catalog), Is.EqualTo("COSV56056643"));
        });
    }

    #endregion

    #region ONCO-TMB-001 — Tumor Mutational Burden

    // -------------------------------------------------------------------------
    // Theory oracle (literal published constant — NOT routed through production)
    // -------------------------------------------------------------------------

    /// <summary>
    /// FDA TMB-High cutoff transcribed LITERALLY from Marcus et al. (2021) — pembrolizumab approved for
    /// solid tumors with "TMB ≥ 10 mut/Mb", inclusive boundary. Kept as a LOCAL literal (not
    /// <c>OncologyAnalyzer.TmbHighThreshold</c>) so a wrong production constant is still caught.
    /// </summary>
    private const double FdaTmbHighCutoff = 10.0;

    /// <summary>Generates a mutation count in [0, 100_000] (≥ 0; spans 0 and large counts).</summary>
    private static Arbitrary<int> TmbCountArbitrary() =>
        Gen.Choose(0, 100_000).ToArbitrary();

    /// <summary>
    /// Generates a finite, strictly-positive target region size in Mb. Drawn from integer milli-Mb in
    /// [1, 5_000] (i.e. 0.001 .. 5.0 Mb) so values stay in a sane range that keeps the n/r quotient exact
    /// to within 1e-12 and avoids floating-point noise.
    /// </summary>
    private static Arbitrary<double> TmbRegionMbArbitrary() =>
        Gen.Choose(1, 5_000).Select(milli => milli / 1_000.0).ToArbitrary();

    /// <summary>
    /// Generates a TMB value ≥ 0 spanning the classification boundary: exact 10.0, just below (9.999),
    /// just above (10.001), 0, plus random finite values straddling the cutoff.
    /// </summary>
    private static Arbitrary<double> TmbValueArbitrary() =>
        Gen.OneOf(
            Gen.Constant(0.0),
            Gen.Constant(9.999),
            Gen.Constant(10.0),
            Gen.Constant(10.001),
            Gen.Choose(0, 20_000).Select(milli => milli / 1_000.0)) // 0.000 .. 20.000
        .ToArbitrary();

    /// <summary>
    /// Generates a list of <see cref="OncologyAnalyzer.SomaticCall"/> with MIXED statuses
    /// (Somatic / Germline / NotDetected), so the overload's Somatic-only counting can be verified
    /// against an independent count. Length 0..40. Numeric fields are arbitrary-but-valid placeholders;
    /// only <c>Status</c> drives the TMB count.
    /// </summary>
    private static Arbitrary<OncologyAnalyzer.SomaticCall[]> SomaticCallListArbitrary() =>
        (from n in Gen.Choose(0, 40)
         from statuses in Gen.Elements(
                 OncologyAnalyzer.SomaticStatus.Somatic,
                 OncologyAnalyzer.SomaticStatus.Germline,
                 OncologyAnalyzer.SomaticStatus.NotDetected)
             .ArrayOf(n)
         select statuses
             .Select((s, i) => new OncologyAnalyzer.SomaticCall(
                 new OncologyAnalyzer.VariantObservation(
                     "chr1", i + 1, "A", "T", 30, 60, 0, 60),
                 TumorVaf: 0.5,
                 NormalVaf: 0.0,
                 Status: s,
                 SomaticScore: 0.9))
             .ToArray())
        .ToArbitrary();

    // -------------------------------------------------------------------------
    // P + INV-01: TMB = count / Mb (exact ratio, mutationCount==0 ⇒ 0)
    // -------------------------------------------------------------------------

    /// <summary>
    /// P + INV-01 (§2.2/§2.4 — TMB = mutationCount / targetRegionMb): for every n ≥ 0 and finite r > 0,
    /// <c>CalculateTMB(n, r)</c> equals the independently computed quotient <c>n / (double)r</c> to within
    /// 1e-12. The oracle is plain division, NOT routed through any production helper.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tmb_Calculate_EqualsCountOverMb()
    {
        return Prop.ForAll(TmbCountArbitrary(), TmbRegionMbArbitrary(), (n, r) =>
        {
            double actual = OncologyAnalyzer.CalculateTMB(n, r);
            double expected = n / r;
            return (Math.Abs(actual - expected) <= 1e-12)
                .Label($"CalculateTMB({n},{r})={actual}, expected n/r={expected}");
        });
    }

    /// <summary>
    /// P + INV-01 boundary (§6.1 — mutationCount = 0 ⇒ TMB = 0): zero somatic mutations give exactly 0
    /// mut/Mb for any finite r > 0 (quotient of zero).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tmb_Calculate_ZeroCount_IsZero()
    {
        return Prop.ForAll(TmbRegionMbArbitrary(), r =>
        {
            double actual = OncologyAnalyzer.CalculateTMB(0, r);
            return (actual == 0.0).Label($"CalculateTMB(0,{r})={actual}, expected 0");
        });
    }

    // -------------------------------------------------------------------------
    // R + INV-02: TMB ≥ 0
    // -------------------------------------------------------------------------

    /// <summary>
    /// R + INV-02 (§2.4 — TMB ≥ 0 for n ≥ 0, r > 0): the quotient of a non-negative count by a positive
    /// region is always non-negative.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tmb_Calculate_IsNonNegative()
    {
        return Prop.ForAll(TmbCountArbitrary(), TmbRegionMbArbitrary(), (n, r) =>
        {
            double actual = OncologyAnalyzer.CalculateTMB(n, r);
            return (actual >= 0.0).Label($"CalculateTMB({n},{r})={actual} is negative");
        });
    }

    // -------------------------------------------------------------------------
    // M + INV-03: monotonicity (strict in count at fixed r; strict-decreasing in r at fixed count>0)
    // -------------------------------------------------------------------------

    /// <summary>
    /// M + INV-03 (§2.4 — non-decreasing in count): at a FIXED region size, increasing the mutation count
    /// strictly increases TMB. Driven over (n, r) with the comparison count n+1; division by a positive
    /// constant is strictly increasing, so the inequality is strict.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tmb_Calculate_StrictlyIncreasingInCount()
    {
        return Prop.ForAll(TmbCountArbitrary(), TmbRegionMbArbitrary(), (n, r) =>
        {
            double lo = OncologyAnalyzer.CalculateTMB(n, r);
            double hi = OncologyAnalyzer.CalculateTMB(n + 1, r);
            return (hi > lo).Label($"TMB({n},{r})={lo} should be < TMB({n + 1},{r})={hi}");
        });
    }

    /// <summary>
    /// M + INV-03 (§2.4 — non-increasing in r): at a FIXED count &gt; 0, increasing the region size
    /// strictly decreases TMB. The count is forced ≥ 1 (a 0 count gives 0 for any r, no strict change),
    /// and the larger region is r + 1 Mb.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tmb_Calculate_StrictlyDecreasingInRegion()
    {
        var positiveCount = Gen.Choose(1, 100_000).ToArbitrary();
        return Prop.ForAll(positiveCount, TmbRegionMbArbitrary(), (n, r) =>
        {
            double smallRegion = OncologyAnalyzer.CalculateTMB(n, r);
            double largeRegion = OncologyAnalyzer.CalculateTMB(n, r + 1.0);
            return (largeRegion < smallRegion)
                .Label($"TMB({n},{r})={smallRegion} should be > TMB({n},{r + 1.0})={largeRegion}");
        });
    }

    // -------------------------------------------------------------------------
    // INV-04: classification at the inclusive ≥ 10 cutoff
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-04 (§2.4/§4.2 — ClassifyTMB(tmb) = High ⇔ tmb ≥ 10, inclusive): for random tmb ≥ 0 (including
    /// the exact boundary 10.0 ⇒ High, just below ⇒ Low, and 0 ⇒ Low), the production result equals the
    /// independent rule recomputed with a LOCAL literal 10.0 cutoff.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tmb_Classify_MatchesInclusiveCutoffOracle()
    {
        return Prop.ForAll(TmbValueArbitrary(), tmb =>
        {
            var actual = OncologyAnalyzer.ClassifyTMB(tmb);
            var expected = tmb >= FdaTmbHighCutoff
                ? OncologyAnalyzer.TmbStatus.High
                : OncologyAnalyzer.TmbStatus.Low;
            return (actual == expected)
                .Label($"ClassifyTMB({tmb})={actual}, expected {expected} (cutoff {FdaTmbHighCutoff})");
        });
    }

    // -------------------------------------------------------------------------
    // SomaticCall overload: counts ONLY Somatic-status calls
    // -------------------------------------------------------------------------

    /// <summary>
    /// SomaticCall overload (§3.1/§5.1 — only <c>Somatic</c> calls counted): for a mixed list of
    /// Somatic / Germline / NotDetected calls, <c>CalculateTMB(calls, r)</c> equals the INDEPENDENT count
    /// of Somatic-status entries divided by r (within 1e-12) — proving Germline and NotDetected are
    /// excluded from the burden.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tmb_CalculateFromCalls_CountsOnlySomatic()
    {
        return Prop.ForAll(SomaticCallListArbitrary(), TmbRegionMbArbitrary(), (calls, r) =>
        {
            double actual = OncologyAnalyzer.CalculateTMB(calls, r);
            int somatic = calls.Count(c => c.Status == OncologyAnalyzer.SomaticStatus.Somatic);
            double expected = somatic / r;
            return (Math.Abs(actual - expected) <= 1e-12)
                .Label($"CalculateTMB(calls[{calls.Length}],{r})={actual}, somatic={somatic}, expected {expected}");
        });
    }

    // -------------------------------------------------------------------------
    // D: determinism (identical inputs ⇒ identical TMB / status)
    // -------------------------------------------------------------------------

    /// <summary>
    /// D (deterministic): repeated calls on identical inputs yield identical TMB (both overloads) and
    /// identical classification — no hidden state or nondeterminism.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tmb_IsDeterministic()
    {
        return Prop.ForAll(SomaticCallListArbitrary(), TmbCountArbitrary(), TmbRegionMbArbitrary(),
            (calls, n, r) =>
            {
                bool intStable = OncologyAnalyzer.CalculateTMB(n, r) == OncologyAnalyzer.CalculateTMB(n, r);
                bool callsStable = OncologyAnalyzer.CalculateTMB(calls, r) == OncologyAnalyzer.CalculateTMB(calls, r);
                double tmb = OncologyAnalyzer.CalculateTMB(n, r);
                bool classifyStable = OncologyAnalyzer.ClassifyTMB(tmb) == OncologyAnalyzer.ClassifyTMB(tmb);
                return (intStable && callsStable && classifyStable)
                    .Label($"non-deterministic: int={intStable}, calls={callsStable}, classify={classifyStable}");
            });
    }

    // -------------------------------------------------------------------------
    // Validation / edge cases (§3.3 / §6.1) and anchors (§7.1)
    // -------------------------------------------------------------------------

    /// <summary>§3.3: negative mutationCount ⇒ ArgumentOutOfRangeException.</summary>
    [Test]
    [Category("Property")]
    public void Tmb_Calculate_NegativeCount_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateTMB(-1, 1.1));
    }

    /// <summary>§3.3/§6.1: non-finite or ≤ 0 targetRegionMb ⇒ ArgumentOutOfRangeException (TMB undefined).</summary>
    [TestCase(0.0)]
    [TestCase(-1.0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    [TestCase(double.NegativeInfinity)]
    [Category("Property")]
    public void Tmb_Calculate_InvalidRegion_Throws(double region)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateTMB(11, region));
    }

    /// <summary>§3.3: null calls collection ⇒ ArgumentNullException.</summary>
    [Test]
    [Category("Property")]
    public void Tmb_CalculateFromCalls_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.CalculateTMB((IEnumerable<OncologyAnalyzer.SomaticCall>)null!, 1.1));
    }

    /// <summary>§3.3: negative or non-finite tmb ⇒ ArgumentOutOfRangeException from ClassifyTMB.</summary>
    [TestCase(-0.0001)]
    [TestCase(-1.0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    [TestCase(double.NegativeInfinity)]
    [Category("Property")]
    public void Tmb_Classify_InvalidValue_Throws(double tmb)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyTMB(tmb));
    }

    /// <summary>§7.1 worked example: 11 mutations / 1.1 Mb = 10.0 mut/Mb ⇒ TMB-High (inclusive cutoff).</summary>
    [Test]
    [Category("Property")]
    public void Tmb_WorkedExample_11Over1Point1_Is10_High()
    {
        double tmb = OncologyAnalyzer.CalculateTMB(11, 1.1);
        Assert.Multiple(() =>
        {
            Assert.That(tmb, Is.EqualTo(10.0).Within(1e-12));
            Assert.That(OncologyAnalyzer.ClassifyTMB(tmb), Is.EqualTo(OncologyAnalyzer.TmbStatus.High));
        });
    }

    /// <summary>§6.1 boundary anchor: tmb = 10.0 ⇒ High (inclusive); tmb = 9.999 ⇒ Low.</summary>
    [TestCase(10.0, OncologyAnalyzer.TmbStatus.High)]
    [TestCase(9.999, OncologyAnalyzer.TmbStatus.Low)]
    [TestCase(0.0, OncologyAnalyzer.TmbStatus.Low)]
    [Category("Property")]
    public void Tmb_Classify_BoundaryAnchors(double tmb, OncologyAnalyzer.TmbStatus expected)
    {
        Assert.That(OncologyAnalyzer.ClassifyTMB(tmb), Is.EqualTo(expected));
    }

    #endregion

    #region ONCO-MSI-001 — Microsatellite Instability (MSI) Detection

    // -------------------------------------------------------------------------
    // Independent oracles (transcribed from Microsatellite_Instability_Detection.md
    // §2.2/§2.4, NOT routed through production constants). The score, the continuous
    // 0.20 cutoff and the Bethesda count rule are recomputed here from LOCAL literals
    // so a self-consistent-but-wrong production constant is still caught.
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-03 cutoff transcribed as a LOCAL literal (MSIsensor2: "msi high: msi score &gt;= 20%").
    /// Deliberately NOT <see cref="OncologyAnalyzer.MsiHighScoreThreshold"/> so the test is an
    /// independent oracle rather than a tautology against the production constant.
    /// </summary>
    private const double MsiHighCutoffOracle = 0.20;

    /// <summary>
    /// INV-01/INV-02 score oracle: score = u / n (exact rational, computed in double). Recomputed
    /// independently of <see cref="OncologyAnalyzer.CalculateMSIScore"/>.
    /// </summary>
    private static double ExpectedMsiScore(int unstableLoci, int totalLoci) =>
        (double)unstableLoci / totalLoci;

    /// <summary>
    /// INV-03 continuous-status oracle, transcribed literally: MSI-H iff score ≥ 0.20 (inclusive),
    /// else MSS. The continuous classifier never yields MSI_Low.
    /// </summary>
    private static OncologyAnalyzer.MsiStatus ExpectedContinuousStatus(double score) =>
        score >= MsiHighCutoffOracle
            ? OncologyAnalyzer.MsiStatus.MSI_High
            : OncologyAnalyzer.MsiStatus.MSS;

    /// <summary>
    /// INV-04 Bethesda oracle, transcribed literally from Boland et al. (1998): depends ONLY on the
    /// absolute unstable-marker count — ≥2 ⇒ MSI-H, exactly 1 ⇒ MSI-L, 0 ⇒ MSS — irrespective of the
    /// panel size (confirmed against <c>ClassifyBethesdaPanel</c>, which branches on the count alone).
    /// </summary>
    private static OncologyAnalyzer.MsiStatus ExpectedBethesdaStatus(int unstableMarkers) =>
        unstableMarkers >= 2
            ? OncologyAnalyzer.MsiStatus.MSI_High
            : unstableMarkers == 1
                ? OncologyAnalyzer.MsiStatus.MSI_Low
                : OncologyAnalyzer.MsiStatus.MSS;

    // -------------------------------------------------------------------------
    // Generators
    // -------------------------------------------------------------------------

    /// <summary>
    /// Generates a contract-valid loci pair (u, n) with n ≥ 1 and 0 ≤ u ≤ n, covering the degenerate
    /// extremes u = 0 (score 0) and u = n (score 1).
    /// </summary>
    private static Arbitrary<(int unstable, int total)> MsiLociArbitrary() =>
        (from total in Gen.Choose(1, 500)
         from unstable in Gen.Choose(0, total)
         select (unstable, total))
        .ToArbitrary();

    /// <summary>
    /// Generates, at a FIXED total, two ordered unstable counts u_lo ≤ u_hi (both ≤ total) for the
    /// monotonicity check.
    /// </summary>
    private static Arbitrary<(int total, int uLo, int uHi)> MsiMonotoneArbitrary() =>
        (from total in Gen.Choose(1, 500)
         from a in Gen.Choose(0, total)
         from b in Gen.Choose(0, total)
         select (total, Math.Min(a, b), Math.Max(a, b)))
        .ToArbitrary();

    /// <summary>
    /// Generates an MSI score as a fraction in [0,1], biased to include the 0.20 boundary and points just
    /// below/above it (so the inclusive cutoff is exercised), plus the extremes 0 and 1.
    /// </summary>
    private static Arbitrary<double> MsiScoreArbitrary() =>
        Gen.OneOf(
            Gen.Choose(0, 1000).Select(x => x / 1000.0),               // uniform over [0,1]
            Gen.Elements(0.0, 0.199999, 0.2, 0.200001, 1.0))          // boundary anchors
        .ToArbitrary();

    /// <summary>
    /// Generates a contract-valid Bethesda count pair (u, total) with total ≥ 1 and 0 ≤ u ≤ total,
    /// over panel sizes both equal to and different from the classical 5 (so the count-only rule is probed).
    /// </summary>
    private static Arbitrary<(int unstable, int total)> BethesdaArbitrary() =>
        (from total in Gen.Choose(1, 20)
         from unstable in Gen.Choose(0, total)
         select (unstable, total))
        .ToArbitrary();

    /// <summary>Generates a non-empty per-locus boolean flag list (length 1..50) for <c>DetectMSI</c>.</summary>
    private static Arbitrary<bool[]> MsiFlagsArbitrary() =>
        (from n in Gen.Choose(1, 50)
         from arr in Gen.Elements(true, false).ArrayOf(n)
         select arr)
        .ToArbitrary();

    // -------------------------------------------------------------------------
    // R + INV-01/INV-02: score = u/n exactly, and ∈ [0,1]
    // -------------------------------------------------------------------------

    /// <summary>
    /// R + INV-01/INV-02: <c>CalculateMSIScore(u,n)</c> equals the independently recomputed u/n (within
    /// 1e-12) and lies in [0,1] for every contract-valid 0 ≤ u ≤ n, n ≥ 1.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Msi_Score_EqualsFractionAndInUnitRange()
    {
        return Prop.ForAll(MsiLociArbitrary(), p =>
        {
            double score = OncologyAnalyzer.CalculateMSIScore(p.unstable, p.total);
            double expected = ExpectedMsiScore(p.unstable, p.total);
            bool matches = Math.Abs(score - expected) < 1e-12;
            bool inRange = score is >= 0.0 and <= 1.0;
            return (matches && inRange)
                .Label($"score {score} != u/n {expected} (u={p.unstable}, n={p.total}), inRange={inRange}");
        });
    }

    // -------------------------------------------------------------------------
    // M: more unstable loci → strictly higher score
    // -------------------------------------------------------------------------

    /// <summary>
    /// M: at a FIXED total, a strictly greater unstable count yields a strictly greater score
    /// (u_lo &lt; u_hi ⇒ score(u_lo) &lt; score(u_hi)); equal counts give equal scores.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Msi_Score_StrictlyMonotoneInUnstableCount()
    {
        return Prop.ForAll(MsiMonotoneArbitrary(), t =>
        {
            double lo = OncologyAnalyzer.CalculateMSIScore(t.uLo, t.total);
            double hi = OncologyAnalyzer.CalculateMSIScore(t.uHi, t.total);
            bool ok = t.uLo < t.uHi ? hi > lo : hi == lo;
            return ok.Label($"n={t.total}: score({t.uLo})={lo}, score({t.uHi})={hi}");
        });
    }

    /// <summary>
    /// M (end-to-end): adding one more <c>true</c> flag to a list of the same length (one stable flag
    /// flipped to unstable) strictly raises <c>DetectMSI(...).Score</c>.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Msi_DetectMSI_MoreTrueFlagsRaisesScore()
    {
        var gen =
            from n in Gen.Choose(1, 50)
            from k in Gen.Choose(0, n - 1)                  // current unstable count, leaving ≥1 stable to flip
            select (n, k);

        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            // base list: k unstable, (n-k) stable; flipped list: one extra unstable.
            var baseFlags = new bool[t.n];
            for (int i = 0; i < t.k; i++)
            {
                baseFlags[i] = true;
            }

            var moreFlags = (bool[])baseFlags.Clone();
            moreFlags[t.k] = true; // flip one stable → unstable (same total length)

            double baseScore = OncologyAnalyzer.DetectMSI(baseFlags).Score;
            double moreScore = OncologyAnalyzer.DetectMSI(moreFlags).Score;
            return (moreScore > baseScore)
                .Label($"n={t.n}, k={t.k}: base={baseScore}, more={moreScore}");
        });
    }

    // -------------------------------------------------------------------------
    // P + INV-03: continuous status = MSI-H iff score ≥ 0.20 (inclusive); never MSI_Low
    // -------------------------------------------------------------------------

    /// <summary>
    /// P + INV-03: <c>ClassifyMSIStatus(score)</c> equals the independent literal-cutoff oracle
    /// (score ≥ 0.20 ⇒ MSI_High else MSS) over random scores in [0,1] including the inclusive 0.20
    /// boundary, and is NEVER MSI_Low.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Msi_ContinuousStatus_MatchesInclusiveCutoffOracle()
    {
        return Prop.ForAll(MsiScoreArbitrary(), score =>
        {
            var actual = OncologyAnalyzer.ClassifyMSIStatus(score);
            var expected = ExpectedContinuousStatus(score);
            bool notLow = actual != OncologyAnalyzer.MsiStatus.MSI_Low;
            return (actual == expected && notLow)
                .Label($"score={score}: actual={actual}, oracle={expected}");
        });
    }

    // -------------------------------------------------------------------------
    // INV-04: Bethesda categorical = f(count only): ≥2→MSI-H, 1→MSI-L, 0→MSS
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-04: <c>ClassifyBethesdaPanel(u,total)</c> equals the independent count-only oracle
    /// (≥2 ⇒ MSI-H, 1 ⇒ MSI-L, 0 ⇒ MSS) over random valid (u,total). Confirms the result depends on the
    /// absolute unstable count alone, irrespective of the panel size.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Msi_Bethesda_MatchesCountOnlyOracle()
    {
        return Prop.ForAll(BethesdaArbitrary(), p =>
        {
            var actual = OncologyAnalyzer.ClassifyBethesdaPanel(p.unstable, p.total);
            var expected = ExpectedBethesdaStatus(p.unstable);
            return (actual == expected)
                .Label($"u={p.unstable}, total={p.total}: actual={actual}, oracle={expected}");
        });
    }

    /// <summary>
    /// INV-04 (count independence): for the SAME unstable count, two different (valid) panel sizes give the
    /// SAME Bethesda status — proving the classification ignores the total.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Msi_Bethesda_IsIndependentOfPanelSize()
    {
        var gen =
            from u in Gen.Choose(0, 20)
            from t1 in Gen.Choose(u <= 0 ? 1 : u, 25)
            from t2 in Gen.Choose(u <= 0 ? 1 : u, 25)
            select (u, t1, t2);

        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var a = OncologyAnalyzer.ClassifyBethesdaPanel(t.u, t.t1);
            var b = OncologyAnalyzer.ClassifyBethesdaPanel(t.u, t.t2);
            return (a == b).Label($"u={t.u}: total {t.t1}→{a}, total {t.t2}→{b}");
        });
    }

    // -------------------------------------------------------------------------
    // DetectMSI end-to-end: counts, score, status all consistent
    // -------------------------------------------------------------------------

    /// <summary>
    /// DetectMSI end-to-end: <c>UnstableLoci</c> = #true, <c>TotalLoci</c> = count, <c>Score</c> = #true/count
    /// (within 1e-12), and <c>Status</c> = <c>ClassifyMSIStatus(Score)</c> — each recomputed independently.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Msi_DetectMSI_FieldsAreConsistent()
    {
        return Prop.ForAll(MsiFlagsArbitrary(), flags =>
        {
            int trueCount = flags.Count(f => f);
            int total = flags.Length;
            var r = OncologyAnalyzer.DetectMSI(flags);

            double expectedScore = (double)trueCount / total;
            var expectedStatus = ExpectedContinuousStatus(expectedScore);

            bool ok = r.UnstableLoci == trueCount
                      && r.TotalLoci == total
                      && Math.Abs(r.Score - expectedScore) < 1e-12
                      && r.Status == expectedStatus;
            return ok.Label(
                $"u={r.UnstableLoci}/{trueCount}, n={r.TotalLoci}/{total}, score={r.Score}/{expectedScore}, status={r.Status}/{expectedStatus}");
        });
    }

    // -------------------------------------------------------------------------
    // D: determinism — identical inputs ⇒ identical results for all four methods
    // -------------------------------------------------------------------------

    /// <summary>
    /// D: <c>CalculateMSIScore</c>, <c>ClassifyMSIStatus</c>, <c>ClassifyBethesdaPanel</c> and <c>DetectMSI</c>
    /// each return identical results on repeated identical inputs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Msi_AllMethods_AreDeterministic()
    {
        var gen =
            from loci in MsiLociArbitrary().Generator
            from score in MsiScoreArbitrary().Generator
            from beth in BethesdaArbitrary().Generator
            from flags in MsiFlagsArbitrary().Generator
            select (loci, score, beth, flags);

        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            bool scoreDet = OncologyAnalyzer.CalculateMSIScore(t.loci.unstable, t.loci.total)
                            == OncologyAnalyzer.CalculateMSIScore(t.loci.unstable, t.loci.total);
            bool contDet = OncologyAnalyzer.ClassifyMSIStatus(t.score)
                           == OncologyAnalyzer.ClassifyMSIStatus(t.score);
            bool bethDet = OncologyAnalyzer.ClassifyBethesdaPanel(t.beth.unstable, t.beth.total)
                           == OncologyAnalyzer.ClassifyBethesdaPanel(t.beth.unstable, t.beth.total);
            bool detectDet = OncologyAnalyzer.DetectMSI(t.flags) == OncologyAnalyzer.DetectMSI(t.flags);
            return (scoreDet && contDet && bethDet && detectDet)
                .Label($"score={scoreDet}, cont={contDet}, beth={bethDet}, detect={detectDet}");
        });
    }

    // -------------------------------------------------------------------------
    // Validation / edge cases (§3.3, §6.1) and worked-example anchors (§7.1)
    // -------------------------------------------------------------------------

    /// <summary>§3.3: <c>CalculateMSIScore</c> rejects totalLoci ≤ 0, unstableLoci &lt; 0, or unstableLoci &gt; totalLoci.</summary>
    [TestCase(0, 0)]    // totalLoci ≤ 0
    [TestCase(1, 0)]    // totalLoci ≤ 0
    [TestCase(0, -1)]   // totalLoci ≤ 0
    [TestCase(-1, 10)]  // unstableLoci < 0
    [TestCase(11, 10)]  // unstableLoci > totalLoci
    [Category("Property")]
    public void Msi_CalculateMSIScore_RejectsInvalidLoci(int unstableLoci, int totalLoci)
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.CalculateMSIScore(unstableLoci, totalLoci));
    }

    /// <summary>§3.3: <c>ClassifyMSIStatus</c> rejects non-finite scores or scores outside [0,1].</summary>
    [TestCase(-0.0001)]
    [TestCase(1.0001)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    [TestCase(double.NegativeInfinity)]
    [Category("Property")]
    public void Msi_ClassifyMSIStatus_RejectsOutOfRangeScore(double score)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyMSIStatus(score));
    }

    /// <summary>§3.3: <c>ClassifyBethesdaPanel</c> rejects totalMarkers ≤ 0, unstable &lt; 0, or unstable &gt; total.</summary>
    [TestCase(0, 0)]
    [TestCase(1, 0)]
    [TestCase(0, -1)]
    [TestCase(-1, 5)]
    [TestCase(6, 5)]
    [Category("Property")]
    public void Msi_ClassifyBethesdaPanel_RejectsInvalidCounts(int unstableMarkers, int totalMarkers)
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.ClassifyBethesdaPanel(unstableMarkers, totalMarkers));
    }

    /// <summary>§3.3: <c>DetectMSI</c> rejects a null sequence.</summary>
    [Test]
    [Category("Property")]
    public void Msi_DetectMSI_NullThrows()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.DetectMSI(null!));
    }

    /// <summary>§3.3/§6.1: <c>DetectMSI</c> rejects an empty sequence (no valid loci).</summary>
    [Test]
    [Category("Property")]
    public void Msi_DetectMSI_EmptyThrows()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.DetectMSI(Array.Empty<bool>()));
    }

    /// <summary>§6.1 boundary anchors: score 0.20 ⇒ MSI-H (inclusive); 0.0 ⇒ MSS; 0.199999 ⇒ MSS; 1.0 ⇒ MSI-H.</summary>
    [TestCase(0.20, OncologyAnalyzer.MsiStatus.MSI_High)]
    [TestCase(0.0, OncologyAnalyzer.MsiStatus.MSS)]
    [TestCase(0.199999, OncologyAnalyzer.MsiStatus.MSS)]
    [TestCase(1.0, OncologyAnalyzer.MsiStatus.MSI_High)]
    [Category("Property")]
    public void Msi_ClassifyMSIStatus_BoundaryAnchors(double score, OncologyAnalyzer.MsiStatus expected)
    {
        Assert.That(OncologyAnalyzer.ClassifyMSIStatus(score), Is.EqualTo(expected));
    }

    /// <summary>§7.1 worked example: 6 unstable of 20 loci ⇒ score 0.30 ⇒ MSI-High.</summary>
    [Test]
    [Category("Property")]
    public void Msi_DetectMSI_WorkedExample_6Of20_IsMsiHigh()
    {
        var flags = new bool[20];
        for (int i = 0; i < 6; i++)
        {
            flags[i] = true;
        }

        var r = OncologyAnalyzer.DetectMSI(flags);
        Assert.Multiple(() =>
        {
            Assert.That(r.UnstableLoci, Is.EqualTo(6));
            Assert.That(r.TotalLoci, Is.EqualTo(20));
            Assert.That(r.Score, Is.EqualTo(0.30).Within(1e-12));
            Assert.That(r.Status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSI_High));
        });
    }

    /// <summary>§7.1 / §6.1 Bethesda anchors: (2,5) ⇒ MSI-H, (1,5) ⇒ MSI-L, (0,5) ⇒ MSS.</summary>
    [TestCase(2, 5, OncologyAnalyzer.MsiStatus.MSI_High)]
    [TestCase(1, 5, OncologyAnalyzer.MsiStatus.MSI_Low)]
    [TestCase(0, 5, OncologyAnalyzer.MsiStatus.MSS)]
    [Category("Property")]
    public void Msi_ClassifyBethesdaPanel_Anchors(int unstable, int total, OncologyAnalyzer.MsiStatus expected)
    {
        Assert.That(OncologyAnalyzer.ClassifyBethesdaPanel(unstable, total), Is.EqualTo(expected));
    }

    #endregion

    #region ONCO-HRD-001 — Homologous Recombination Deficiency (HRD) Score

    // -------------------------------------------------------------------------
    // Independent oracles (transcribed from HRD_Score.md §2.2/§2.4, NOT routed
    // through production constants). The unweighted sum and the ≥ 42 cutoff are
    // recomputed here from LOCAL literals so a self-consistent-but-wrong
    // production constant (e.g. a drifted threshold) is still caught.
    //   INV-01 HRD = LOH + TAI + LST   (unweighted sum, Telli 2016)
    //   INV-02 sum is order-independent (integer addition; "unweighted")
    //   INV-03 status = HrdHigh iff score ≥ 42 (boundary inclusive)
    //   INV-04 score ≥ 0 and each component ≥ 0 (event counts)
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-03 local literal cutoff (Telli et al. 2016 "HRD score ≥42"), transcribed independently of
    /// the production constant <see cref="OncologyAnalyzer.HrdHighScoreThreshold"/> so threshold drift is caught.
    /// </summary>
    private const int HrdHighCutoff = 42;

    /// <summary>INV-01 score oracle: the unweighted sum of the three genomic-scar component counts.</summary>
    private static int ExpectedHrdScore(int loh, int tai, int lst) => loh + tai + lst;

    /// <summary>INV-03 status oracle: HrdHigh iff score ≥ 42 (inclusive), recomputed from the local literal.</summary>
    private static OncologyAnalyzer.HrdStatus ExpectedHrdStatus(int score) =>
        score >= HrdHighCutoff ? OncologyAnalyzer.HrdStatus.HrdHigh : OncologyAnalyzer.HrdStatus.HrdNegative;

    // -------------------------------------------------------------------------
    // Generators
    // -------------------------------------------------------------------------

    /// <summary>
    /// Generates a triple of non-negative genomic-scar component counts (INV-04). The 0..60 range per
    /// component spans both sides of the 42 cutoff for the summed score, including all-zero (near-diploid).
    /// </summary>
    private static Arbitrary<(int loh, int tai, int lst)> HrdComponentTripleArbitrary() =>
        (from loh in Gen.Choose(0, 60)
         from tai in Gen.Choose(0, 60)
         from lst in Gen.Choose(0, 60)
         select (loh, tai, lst))
        .ToArbitrary();

    /// <summary>
    /// Generates a non-negative HRD score (INV-04), drawn 0..120 so it straddles the 42 cutoff densely,
    /// with the exact boundaries 41 and 42 injected so the inclusive comparison is always exercised.
    /// </summary>
    private static Arbitrary<int> HrdScoreArbitrary() =>
        Gen.OneOf(Gen.Choose(0, 120), Gen.Elements(41, 42)).ToArbitrary();

    // -------------------------------------------------------------------------
    // P + INV-01 (additive) / INV-04 (R: score ≥ 0)
    // -------------------------------------------------------------------------

    /// <summary>
    /// P + INV-01 + INV-04 (R): <c>CalculateHRDScore(loh, tai, lst)</c> equals the independently recomputed
    /// unweighted sum <c>loh + tai + lst</c>, and the result is ≥ 0 for non-negative components.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hrd_Score_IsUnweightedSumAndNonNegative()
    {
        return Prop.ForAll(HrdComponentTripleArbitrary(), t =>
        {
            int score = OncologyAnalyzer.CalculateHRDScore(t.loh, t.tai, t.lst);
            int expected = ExpectedHrdScore(t.loh, t.tai, t.lst);
            return (score == expected && score >= 0)
                .Label($"score {score} != sum {expected} (loh={t.loh}, tai={t.tai}, lst={t.lst})");
        });
    }

    // -------------------------------------------------------------------------
    // INV-02 (commutative / order-independent)
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-02: <c>CalculateHRDScore</c> is invariant under every permutation of its three arguments
    /// (unweighted integer addition is commutative). All 6 permutations of (loh, tai, lst) must agree.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hrd_Score_IsOrderIndependent()
    {
        return Prop.ForAll(HrdComponentTripleArbitrary(), t =>
        {
            int s = OncologyAnalyzer.CalculateHRDScore(t.loh, t.tai, t.lst);
            int[] perms =
            {
                OncologyAnalyzer.CalculateHRDScore(t.loh, t.lst, t.tai),
                OncologyAnalyzer.CalculateHRDScore(t.tai, t.loh, t.lst),
                OncologyAnalyzer.CalculateHRDScore(t.tai, t.lst, t.loh),
                OncologyAnalyzer.CalculateHRDScore(t.lst, t.loh, t.tai),
                OncologyAnalyzer.CalculateHRDScore(t.lst, t.tai, t.loh),
            };
            return perms.All(p => p == s)
                .Label($"permutations disagree with {s} (loh={t.loh}, tai={t.tai}, lst={t.lst})");
        });
    }

    // -------------------------------------------------------------------------
    // M (more genomic scars → higher score)
    // -------------------------------------------------------------------------

    /// <summary>Generator for the monotonicity check: a base triple plus a strictly-positive bump (1..30).</summary>
    private static Arbitrary<(int loh, int tai, int lst, int bump)> HrdBumpArbitrary() =>
        (from loh in Gen.Choose(0, 60)
         from tai in Gen.Choose(0, 60)
         from lst in Gen.Choose(0, 60)
         from bump in Gen.Choose(1, 30)
         select (loh, tai, lst, bump))
        .ToArbitrary();

    /// <summary>
    /// M: increasing ANY single component (others fixed) by a strictly-positive amount strictly increases the
    /// score — more genomic scars ⇒ higher HRD. Checked independently for each of the three components.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hrd_Score_StrictlyIncreasesWhenAnyComponentIncreases()
    {
        return Prop.ForAll(HrdBumpArbitrary(), t =>
        {
            int baseScore = OncologyAnalyzer.CalculateHRDScore(t.loh, t.tai, t.lst);
            int upLoh = OncologyAnalyzer.CalculateHRDScore(t.loh + t.bump, t.tai, t.lst);
            int upTai = OncologyAnalyzer.CalculateHRDScore(t.loh, t.tai + t.bump, t.lst);
            int upLst = OncologyAnalyzer.CalculateHRDScore(t.loh, t.tai, t.lst + t.bump);
            return (upLoh > baseScore && upTai > baseScore && upLst > baseScore)
                .Label($"not strictly increasing: base={baseScore}, +loh={upLoh}, +tai={upTai}, +lst={upLst} (bump={t.bump})");
        });
    }

    // -------------------------------------------------------------------------
    // INV-03 (classification at the inclusive 42 cutoff)
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-03: <c>ClassifyHRDStatus(score)</c> equals <c>score ≥ 42 ? HrdHigh : HrdNegative</c>, recomputed from
    /// the LOCAL literal 42 (boundary inclusive). Driven over random scores ≥ 0 incl. the boundaries 41 and 42.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hrd_Status_MatchesInclusive42Cutoff()
    {
        return Prop.ForAll(HrdScoreArbitrary(), score =>
        {
            var status = OncologyAnalyzer.ClassifyHRDStatus(score);
            var expected = ExpectedHrdStatus(score);
            return (status == expected)
                .Label($"status {status} != oracle {expected} (score={score}, cutoff={HrdHighCutoff})");
        });
    }

    // -------------------------------------------------------------------------
    // DetectHRD end-to-end
    // -------------------------------------------------------------------------

    /// <summary>
    /// End-to-end (INV-01 + INV-03): <c>DetectHRD</c> returns the input components unchanged, a Score equal to the
    /// independent sum, and a Status equal to the independent inclusive-42 classification of that score.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hrd_DetectHRD_ComposesSumAndClassification()
    {
        return Prop.ForAll(HrdComponentTripleArbitrary(), t =>
        {
            var result = OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(t.loh, t.tai, t.lst));
            int expectedScore = ExpectedHrdScore(t.loh, t.tai, t.lst);
            var expectedStatus = ExpectedHrdStatus(expectedScore);
            bool componentsEqual = result.Components.Loh == t.loh
                                   && result.Components.Tai == t.tai
                                   && result.Components.Lst == t.lst;
            return (componentsEqual && result.Score == expectedScore && result.Status == expectedStatus)
                .Label($"components/score/status mismatch: {result} vs sum {expectedScore}/{expectedStatus} " +
                       $"(loh={t.loh}, tai={t.tai}, lst={t.lst})");
        });
    }

    // -------------------------------------------------------------------------
    // D (determinism)
    // -------------------------------------------------------------------------

    /// <summary>
    /// D (determinism): identical inputs yield identical results across repeated calls for all three methods
    /// (<c>CalculateHRDScore</c>, <c>ClassifyHRDStatus</c>, <c>DetectHRD</c>) — pure, side-effect-free.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hrd_AllMethods_AreDeterministic()
    {
        return Prop.ForAll(HrdComponentTripleArbitrary(), t =>
        {
            int s1 = OncologyAnalyzer.CalculateHRDScore(t.loh, t.tai, t.lst);
            int s2 = OncologyAnalyzer.CalculateHRDScore(t.loh, t.tai, t.lst);
            var c1 = OncologyAnalyzer.ClassifyHRDStatus(s1);
            var c2 = OncologyAnalyzer.ClassifyHRDStatus(s1);
            var d1 = OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(t.loh, t.tai, t.lst));
            var d2 = OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(t.loh, t.tai, t.lst));
            return (s1 == s2 && c1 == c2 && d1.Equals(d2))
                .Label($"non-deterministic: score {s1}/{s2}, status {c1}/{c2}, result {d1} vs {d2}");
        });
    }

    // -------------------------------------------------------------------------
    // Validation / edge cases (§3.3, §6.1) — explicit anchors and throw contracts
    // -------------------------------------------------------------------------

    /// <summary>§7.1 worked example: (20, 15, 12) ⇒ score 47 ⇒ HrdHigh (47 ≥ 42).</summary>
    [Test]
    [Category("Property")]
    public void Hrd_DetectHRD_WorkedExample_20_15_12_Is47AndHrdHigh()
    {
        var result = OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(20, 15, 12));
        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(47));
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdHigh));
            Assert.That(result.Components.Loh, Is.EqualTo(20));
            Assert.That(result.Components.Tai, Is.EqualTo(15));
            Assert.That(result.Components.Lst, Is.EqualTo(12));
        });
    }

    /// <summary>
    /// §6.1 classification anchors: score 42 ⇒ HrdHigh (inclusive boundary), 41 ⇒ HrdNegative,
    /// 0 ⇒ HrdNegative (all-zero near-diploid).
    /// </summary>
    [TestCase(42, OncologyAnalyzer.HrdStatus.HrdHigh)]
    [TestCase(41, OncologyAnalyzer.HrdStatus.HrdNegative)]
    [TestCase(0, OncologyAnalyzer.HrdStatus.HrdNegative)]
    [Category("Property")]
    public void Hrd_ClassifyHRDStatus_BoundaryAnchors(int score, OncologyAnalyzer.HrdStatus expected)
    {
        Assert.That(OncologyAnalyzer.ClassifyHRDStatus(score), Is.EqualTo(expected));
    }

    /// <summary>§6.1 all-components-zero anchor: (0, 0, 0) ⇒ score 0 ⇒ HrdNegative.</summary>
    [Test]
    [Category("Property")]
    public void Hrd_DetectHRD_AllZeroComponents_IsZeroAndHrdNegative()
    {
        var result = OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(0, 0, 0));
        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(0));
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdNegative));
        });
    }

    /// <summary>§3.3: any negative component to <c>CalculateHRDScore</c> throws ArgumentOutOfRangeException.</summary>
    [TestCase(-1, 0, 0)]
    [TestCase(0, -1, 0)]
    [TestCase(0, 0, -1)]
    [Category("Property")]
    public void Hrd_CalculateHRDScore_NegativeComponent_Throws(int loh, int tai, int lst)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateHRDScore(loh, tai, lst));
    }

    /// <summary>§3.3: any negative component to <c>DetectHRD</c> throws ArgumentOutOfRangeException.</summary>
    [TestCase(-1, 0, 0)]
    [TestCase(0, -1, 0)]
    [TestCase(0, 0, -1)]
    [Category("Property")]
    public void Hrd_DetectHRD_NegativeComponent_Throws(int loh, int tai, int lst)
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(loh, tai, lst)));
    }

    /// <summary>§3.3: a negative score to <c>ClassifyHRDStatus</c> throws ArgumentOutOfRangeException.</summary>
    [TestCase(-1)]
    [TestCase(-42)]
    [Category("Property")]
    public void Hrd_ClassifyHRDStatus_NegativeScore_Throws(int score)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyHRDStatus(score));
    }

    #endregion

    #region ONCO-LOH-001 — Loss of Heterozygosity (HRD-LOH)

    // -------------------------------------------------------------------------
    // Independent scarHRD/Abkevich oracle, transcribed from
    // Loss_Of_Heterozygosity.md §2.2/§2.4 (INV-01..INV-06) and §4.1, NOT routed
    // through any production helper. The LOH predicate, the whole-chromosome
    // exclusion, the adjacent-same-state merge and the strict 15 Mb count are all
    // recomputed here from literals, so a self-consistent-but-wrong implementation
    // is still caught.
    //
    // Checklist-row refinement (the row wording is loose / partly bogus):
    //   * "P: LOH ⟺ BAF → 0/1" — the implementation takes NO BAF input. The real
    //     contract (INV-03 / §2.2) is: a segment is LOH iff cn_minor == 0 AND
    //     cn_major != 0. That is what we test.
    //   * "M: lower BAF-dev threshold → ≥ LOH" — there is NO tunable threshold; the
    //     region-size limit is the FIXED constant 15,000,000 bp. This monotonicity
    //     item does NOT apply and is dropped; in its place we pin the fixed strict
    //     15 Mb boundary (INV-04, §6.1).
    // -------------------------------------------------------------------------

    /// <summary>Documented strict minimum HRD-LOH region length, 15 Mb (literal, Abkevich 2012 / scarHRD <c>sizelimitLOH = 15e6</c>).</summary>
    private const long LohMinRegionLengthBp = 15_000_000L;

    /// <summary>The handful of chromosomes used by the generators (kept small so collisions and merges are exercised).</summary>
    private static readonly string[] LohChromosomes = { "1", "2", "3", "X" };

    /// <summary>
    /// INV-03 LOH predicate, transcribed literally from §2.2: a segment is LOH iff the minor-allele copy
    /// number is 0 and the major-allele copy number is non-zero (homozygous deletion 0|0 is NOT LOH;
    /// heterozygous retention minor≠0 is NOT LOH).
    /// </summary>
    private static bool OracleIsLoh(OncologyAnalyzer.AlleleSpecificSegment s)
        => s.MinorCopyNumber == 0 && s.MajorCopyNumber != 0;

    /// <summary>
    /// Independent scarHRD HRD-LOH oracle (INV-01,03,04,05,06). Reproduces §4.1 from scratch: group by
    /// chromosome; EXCLUDE any chromosome whose every segment has cn_minor == 0 (whole-chromosome LOH,
    /// scarHRD <c>chrDel</c>); within each remaining chromosome sort by start and merge adjacent segments of
    /// the SAME LOH state when the gap is ≤ 1 bp; count merged LOH-state runs whose length (end − start) is
    /// STRICTLY &gt; 15,000,000 bp. Returns both the count and the surviving (chr,start,end) runs.
    /// </summary>
    private static (int score, List<(string chr, long start, long end)> regions) OracleHrdLoh(
        IEnumerable<OncologyAnalyzer.AlleleSpecificSegment> segments)
    {
        // Group preserving nothing about order — the oracle is set-based per chromosome (INV-06).
        var byChr = new Dictionary<string, List<OncologyAnalyzer.AlleleSpecificSegment>>(StringComparer.Ordinal);
        foreach (var s in segments)
        {
            if (!byChr.TryGetValue(s.Chromosome, out var list))
            {
                list = new List<OncologyAnalyzer.AlleleSpecificSegment>();
                byChr[s.Chromosome] = list;
            }

            list.Add(s);
        }

        var regions = new List<(string chr, long start, long end)>();
        foreach (var kv in byChr)
        {
            var group = kv.Value;

            // Whole-chromosome LOH: every segment has minor == 0 → excluded entirely (INV-05).
            if (group.All(s => s.MinorCopyNumber == 0))
            {
                continue;
            }

            // Sort by start (ties by end) then merge adjacent same-LOH-state runs (gap ≤ 1 bp).
            var sorted = group.OrderBy(s => s.Start).ThenBy(s => s.End).ToList();
            var runs = new List<(bool loh, long start, long end)>();
            foreach (var s in sorted)
            {
                bool loh = OracleIsLoh(s);
                if (runs.Count == 0)
                {
                    runs.Add((loh, s.Start, s.End));
                    continue;
                }

                var last = runs[^1];
                bool sameState = last.loh == loh;
                bool adjacent = s.Start - last.end <= 1L;
                if (sameState && adjacent)
                {
                    runs[^1] = (last.loh, last.start, Math.Max(last.end, s.End));
                }
                else
                {
                    runs.Add((loh, s.Start, s.End));
                }
            }

            foreach (var run in runs)
            {
                // INV-04: strict > 15 Mb; INV-03: only LOH-state runs counted.
                if (run.loh && (run.end - run.start) > LohMinRegionLengthBp)
                {
                    regions.Add((kv.Key, run.start, run.end));
                }
            }
        }

        return (regions.Count, regions);
    }

    /// <summary>
    /// Independent LOH-fraction oracle (INV-02, §5.2): (Σ LOH-segment lengths on <paramref name="chromosome"/>)
    /// ÷ (Σ all covered segment lengths on that chromosome), with NO 15 Mb filter and NO whole-chromosome
    /// exclusion. A chromosome absent from the input (zero covered length) yields 0.0.
    /// </summary>
    private static double OracleLohFraction(
        IEnumerable<OncologyAnalyzer.AlleleSpecificSegment> segments, string chromosome)
    {
        long total = 0;
        long loh = 0;
        foreach (var s in segments)
        {
            if (!string.Equals(s.Chromosome, chromosome, StringComparison.Ordinal))
            {
                continue;
            }

            total += s.End - s.Start;
            if (OracleIsLoh(s))
            {
                loh += s.End - s.Start;
            }
        }

        return total == 0 ? 0.0 : (double)loh / total;
    }

    // -------------------------------------------------------------------------
    // Generators — non-overlapping ascending segments over a few chromosomes,
    // with varied major/minor CN (LOH / het / homozygous-deletion) and lengths
    // straddling the 15 Mb boundary so the strict filter is exercised on both
    // sides.
    // -------------------------------------------------------------------------

    /// <summary>
    /// Generates a single chromosome's segment list: contiguous, ascending, non-overlapping segments whose
    /// individual lengths straddle the 15 Mb boundary, with copy numbers spanning LOH (minor 0, major≥1),
    /// heterozygous (minor≥1) and homozygous deletion (0|0). Gaps between segments are 0 or 1 bp so the
    /// merge step is exercised.
    /// </summary>
    private static Gen<OncologyAnalyzer.AlleleSpecificSegment[]> SegmentsForChromosomeGen(string chr) =>
        from n in Gen.Choose(0, 5)
        from lengthsMb in Gen.Choose(1, 40).ArrayOf(n)   // 1..40 Mb per segment (straddles 15 Mb)
        from gaps in Gen.Choose(0, 1).ArrayOf(n)          // 0 or 1 bp gaps → adjacency/merge
        from majors in Gen.Choose(0, 4).ArrayOf(n)        // major CN incl. 0 (homozygous deletion)
        from minors in Gen.Choose(0, 3).ArrayOf(n)        // minor CN incl. 0 (LOH candidate)
        select BuildChromosomeSegments(chr, lengthsMb, gaps, majors, minors);

    private static OncologyAnalyzer.AlleleSpecificSegment[] BuildChromosomeSegments(
        string chr, int[] lengthsMb, int[] gaps, int[] majors, int[] minors)
    {
        var result = new OncologyAnalyzer.AlleleSpecificSegment[lengthsMb.Length];
        long cursor = 0;
        for (int i = 0; i < lengthsMb.Length; i++)
        {
            long start = cursor + (i == 0 ? 0 : gaps[i]);
            long end = start + (long)lengthsMb[i] * 1_000_000L;
            result[i] = new OncologyAnalyzer.AlleleSpecificSegment(chr, start, end, majors[i], minors[i]);
            cursor = end;
        }

        return result;
    }

    /// <summary>
    /// Generates a full multi-chromosome segment set (contract-valid: End &gt; Start, CN ≥ 0) over the
    /// generator's small chromosome pool. Some chromosomes may be whole-chromosome LOH, some empty.
    /// </summary>
    private static Arbitrary<OncologyAnalyzer.AlleleSpecificSegment[]> SegmentSetArbitrary() =>
        (from c0 in SegmentsForChromosomeGen(LohChromosomes[0])
         from c1 in SegmentsForChromosomeGen(LohChromosomes[1])
         from c2 in SegmentsForChromosomeGen(LohChromosomes[2])
         from c3 in SegmentsForChromosomeGen(LohChromosomes[3])
         select c0.Concat(c1).Concat(c2).Concat(c3).ToArray())
        .ToArbitrary();

    /// <summary>FsCheck-friendly shuffle of a segment array (drives INV-06 order independence).</summary>
    private static Gen<OncologyAnalyzer.AlleleSpecificSegment[]> ShuffleGen(
        OncologyAnalyzer.AlleleSpecificSegment[] source) =>
        Gen.Choose(int.MinValue, int.MaxValue).ArrayOf(source.Length)
            .Select(keys => source
                .Select((seg, idx) => (seg, key: keys.Length == 0 ? 0 : keys[idx]))
                .OrderBy(t => t.key)
                .Select(t => t.seg)
                .ToArray());

    // -------------------------------------------------------------------------
    // HRD-LOH score vs. independent oracle (INV-01, INV-03, INV-04, INV-05)
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-01/03/04/05: <c>CalculateHrdLohScore</c> and <c>DetectLOH(...).Score</c> both equal the independent
    /// scarHRD oracle (LOH = minor 0 &amp; major≠0; whole-chromosome LOH excluded; adjacent same-state merge;
    /// strict &gt; 15 Mb), and the score is ≥ 0 (INV-01).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Loh_HrdScore_MatchesScarHrdOracle()
    {
        return Prop.ForAll(SegmentSetArbitrary(), segs =>
        {
            int expected = OracleHrdLoh(segs).score;
            int score = OncologyAnalyzer.CalculateHrdLohScore(segs);
            int detectScore = OncologyAnalyzer.DetectLOH(segs).Score;
            return (score == expected && detectScore == expected && score >= 0)
                .Label($"score={score}, detect={detectScore}, oracle={expected}");
        });
    }

    /// <summary>
    /// INV-06: the HRD-LOH score is invariant under input-segment-order shuffles (per-chromosome aggregation
    /// is set-based). The shuffled score equals the original score equals the oracle.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Loh_HrdScore_OrderIndependent()
    {
        return Prop.ForAll(SegmentSetArbitrary(), segs =>
            Prop.ForAll(ShuffleGen(segs).ToArbitrary(), shuffled =>
            {
                int baseline = OncologyAnalyzer.CalculateHrdLohScore(segs);
                int permuted = OncologyAnalyzer.CalculateHrdLohScore(shuffled);
                int oracle = OracleHrdLoh(segs).score;
                return (baseline == permuted && permuted == oracle)
                    .Label($"baseline={baseline}, permuted={permuted}, oracle={oracle}");
            }));
    }

    // -------------------------------------------------------------------------
    // Per-region invariants (R, INV-03/04/05)
    // -------------------------------------------------------------------------

    /// <summary>
    /// R + INV-03/04/05: every region returned by <c>DetectLOH</c> has Start &lt; End, length (End − Start)
    /// strictly &gt; 15 Mb, a non-empty chromosome that is NOT whole-chromosome-LOH, and matches a region the
    /// independent oracle also produced (LOH state on a kept chromosome). <c>Length</c> equals End − Start.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Loh_DetectedRegions_SatisfyContract()
    {
        return Prop.ForAll(SegmentSetArbitrary(), segs =>
        {
            var result = OncologyAnalyzer.DetectLOH(segs);
            var oracleRegions = OracleHrdLoh(segs).regions
                .Select(r => (r.chr, r.start, r.end))
                .ToHashSet();

            foreach (var region in result.Regions)
            {
                bool positions = region.Start < region.End;
                bool lengthField = region.Length == region.End - region.Start;
                bool strict15 = region.End - region.Start > LohMinRegionLengthBp;
                bool inOracle = oracleRegions.Contains((region.Chromosome, region.Start, region.End));
                if (!(positions && lengthField && strict15 && inOracle))
                {
                    return false.Label(
                        $"bad region chr={region.Chromosome} [{region.Start},{region.End}] len={region.Length} " +
                        $"pos={positions} lenField={lengthField} strict15={strict15} inOracle={inOracle}");
                }
            }

            // Regions.Count must equal Score (and the oracle count).
            return (result.Regions.Count == result.Score && result.Score == oracleRegions.Count)
                .Label($"count={result.Regions.Count}, score={result.Score}, oracle={oracleRegions.Count}");
        });
    }

    // -------------------------------------------------------------------------
    // CalculateLOHFraction (INV-02, §5.2)
    // -------------------------------------------------------------------------

    /// <summary>
    /// INV-02 + §5.2: <c>CalculateLOHFraction</c> lies in [0, 1] and equals the independent fraction oracle
    /// (Σ LOH lengths / Σ covered lengths on the chromosome; no size filter, no whole-chromosome exclusion)
    /// for every chromosome in the pool. Absent chromosomes (zero coverage) give 0.0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Loh_Fraction_MatchesOracleAndInUnitRange()
    {
        return Prop.ForAll(SegmentSetArbitrary(), segs =>
        {
            foreach (var chr in LohChromosomes.Append("absent-chromosome"))
            {
                double frac = OncologyAnalyzer.CalculateLOHFraction(segs, chr);
                double expected = OracleLohFraction(segs, chr);
                bool inRange = frac is >= 0.0 and <= 1.0;
                bool matches = Math.Abs(frac - expected) < 1e-12;
                if (!(inRange && matches))
                {
                    return false.Label($"chr={chr}: frac={frac}, oracle={expected}, inRange={inRange}");
                }
            }

            return true.ToProperty();
        });
    }

    // -------------------------------------------------------------------------
    // D (determinism)
    // -------------------------------------------------------------------------

    /// <summary>
    /// D: identical input ⇒ identical result. Two independent calls over the same segments yield the same
    /// score, the same ordered region list, and the same per-chromosome fraction.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Loh_Deterministic()
    {
        return Prop.ForAll(SegmentSetArbitrary(), segs =>
        {
            var a = OncologyAnalyzer.DetectLOH(segs);
            var b = OncologyAnalyzer.DetectLOH(segs);
            bool sameScore = a.Score == b.Score;
            bool sameRegions = a.Regions.SequenceEqual(b.Regions);
            bool sameFraction = LohChromosomes.All(chr =>
                OncologyAnalyzer.CalculateLOHFraction(segs, chr)
                    == OncologyAnalyzer.CalculateLOHFraction(segs, chr));
            return (sameScore && sameRegions && sameFraction)
                .Label($"sameScore={sameScore}, sameRegions={sameRegions}, sameFraction={sameFraction}");
        });
    }

    // -------------------------------------------------------------------------
    // Validation / edge anchors (§3.3, §6.1)
    // -------------------------------------------------------------------------

    /// <summary>§3.3: null segments ⇒ ArgumentNullException on all three entry points.</summary>
    [Test]
    [Category("Property")]
    public void Loh_NullSegments_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.DetectLOH(null!));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CalculateHrdLohScore(null!));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CalculateLOHFraction(null!, "1"));
        });
    }

    /// <summary>§3.3: a null chromosome to <c>CalculateLOHFraction</c> ⇒ ArgumentNullException.</summary>
    [Test]
    [Category("Property")]
    public void Loh_NullChromosome_Throws()
    {
        var segs = new[] { new OncologyAnalyzer.AlleleSpecificSegment("1", 0, 1_000_000, 1, 0) };
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CalculateLOHFraction(segs, null!));
    }

    /// <summary>§3.3 / §6.1: a segment with End ≤ Start ⇒ ArgumentException.</summary>
    [TestCase(100L, 100L)]
    [TestCase(100L, 50L)]
    [Category("Property")]
    public void Loh_NonPositiveLength_Throws(long start, long end)
    {
        var segs = new[] { new OncologyAnalyzer.AlleleSpecificSegment("1", start, end, 1, 0) };
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.DetectLOH(segs));
    }

    /// <summary>§3.3 / §6.1: a negative copy number ⇒ ArgumentException.</summary>
    [TestCase(-1, 0)]
    [TestCase(0, -1)]
    [Category("Property")]
    public void Loh_NegativeCopyNumber_Throws(int major, int minor)
    {
        var segs = new[] { new OncologyAnalyzer.AlleleSpecificSegment("1", 0, 20_000_000, major, minor) };
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.DetectLOH(segs));
    }

    /// <summary>§6.1: empty input ⇒ score 0 and fraction 0.</summary>
    [Test]
    [Category("Property")]
    public void Loh_EmptyInput_ScoreAndFractionZero()
    {
        var empty = Array.Empty<OncologyAnalyzer.AlleleSpecificSegment>();
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CalculateHrdLohScore(empty), Is.EqualTo(0));
            Assert.That(OncologyAnalyzer.DetectLOH(empty).Score, Is.EqualTo(0));
            Assert.That(OncologyAnalyzer.CalculateLOHFraction(empty, "1"), Is.EqualTo(0.0));
        });
    }

    /// <summary>
    /// INV-04 / §6.1: an LOH region of length EXACTLY 15,000,000 bp is NOT counted (strict &gt;); 15,000,001 IS
    /// counted. A trailing heterozygous segment keeps the chromosome out of the whole-chromosome-LOH
    /// exclusion (INV-05), so the size filter alone is under test.
    /// </summary>
    [TestCase(15_000_000L, 0)]
    [TestCase(15_000_001L, 1)]
    [Category("Property")]
    public void Loh_FifteenMbStrictBoundary(long lohEnd, int expectedScore)
    {
        var segs = new[]
        {
            new OncologyAnalyzer.AlleleSpecificSegment("1", 0, lohEnd, 1, 0),
            new OncologyAnalyzer.AlleleSpecificSegment("1", lohEnd + 1, lohEnd + 5_000_000, 1, 1),
        };
        Assert.That(OncologyAnalyzer.CalculateHrdLohScore(segs), Is.EqualTo(expectedScore));
    }

    /// <summary>INV-03 / §6.1: homozygous deletion (minor 0, major 0) is NOT LOH ⇒ not counted, even &gt; 15 Mb.</summary>
    [Test]
    [Category("Property")]
    public void Loh_HomozygousDeletion_NotCounted()
    {
        var segs = new[]
        {
            new OncologyAnalyzer.AlleleSpecificSegment("1", 0, 20_000_000, 0, 0),
            new OncologyAnalyzer.AlleleSpecificSegment("1", 20_000_000, 60_000_000, 1, 1),
        };
        Assert.That(OncologyAnalyzer.CalculateHrdLohScore(segs), Is.EqualTo(0));
    }

    /// <summary>INV-03 / §6.1: heterozygous retention (minor ≠ 0) is NOT LOH ⇒ not counted, even &gt; 15 Mb.</summary>
    [Test]
    [Category("Property")]
    public void Loh_HeterozygousRetained_NotCounted()
    {
        var segs = new[]
        {
            new OncologyAnalyzer.AlleleSpecificSegment("1", 0, 30_000_000, 2, 1),
            new OncologyAnalyzer.AlleleSpecificSegment("1", 30_000_000, 60_000_000, 1, 1),
        };
        Assert.That(OncologyAnalyzer.CalculateHrdLohScore(segs), Is.EqualTo(0));
    }

    /// <summary>
    /// INV-05 / §6.1: a whole-chromosome-LOH chromosome (every segment minor 0) is excluded — its &gt; 15 Mb
    /// LOH run is NOT counted, while an equivalent run on a non-whole-chromosome-LOH chromosome IS counted.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Loh_WholeChromosomeLoh_NotCounted()
    {
        // chr1: every segment LOH → whole-chromosome LOH → excluded.
        // chr2: a 20 Mb LOH run plus a het segment → NOT whole-chromosome → counted once.
        var segs = new[]
        {
            new OncologyAnalyzer.AlleleSpecificSegment("1", 0, 20_000_000, 1, 0),
            new OncologyAnalyzer.AlleleSpecificSegment("1", 20_000_000, 60_000_000, 2, 0),
            new OncologyAnalyzer.AlleleSpecificSegment("2", 0, 20_000_000, 1, 0),
            new OncologyAnalyzer.AlleleSpecificSegment("2", 20_000_000, 60_000_000, 1, 1),
        };
        Assert.That(OncologyAnalyzer.CalculateHrdLohScore(segs), Is.EqualTo(1));
    }

    /// <summary>
    /// §7.1 worked example: 20 Mb LOH on "1" (counted) + het on "1" + 10 Mb LOH on "2" (≤ 15 Mb, not counted)
    /// ⇒ score 1; CalculateLOHFraction("1") = 20M / 60M = 0.3333…
    /// </summary>
    [Test]
    [Category("Property")]
    public void Loh_WorkedExample_Section71()
    {
        var segments = new[]
        {
            new OncologyAnalyzer.AlleleSpecificSegment("1", 0, 20_000_000, 1, 0),
            new OncologyAnalyzer.AlleleSpecificSegment("1", 20_000_000, 60_000_000, 1, 1),
            new OncologyAnalyzer.AlleleSpecificSegment("2", 0, 10_000_000, 2, 0),
        };
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CalculateHrdLohScore(segments), Is.EqualTo(1));
            Assert.That(
                OncologyAnalyzer.CalculateLOHFraction(segments, "1"),
                Is.EqualTo(20_000_000.0 / 60_000_000.0).Within(1e-12));
        });
    }

    #endregion

    #region ONCO-SIG-001 — SBS-96 Trinucleotide Context Catalog

    // =========================================================================
    // ONCO-SIG-001 — SBS-96 Trinucleotide Context Catalog
    // Doc: docs/algorithms/Oncology/SBS96_Trinucleotide_Context_Catalog.md
    // Source: Alexandrov et al. (2013) Nature 500:415-421; COSMIC SBS96;
    //         Bergstrom et al. (2019) BMC Genomics 20:685.
    // Checklist row #96: R exactly 96 SBS channels; P each SNV → one
    //   pyrimidine-centred trinucleotide; P channel counts sum = #SNVs;
    //   D deterministic.
    // Invariants (doc §2.4): INV-01 every channel's ref base is a pyrimidine
    //   (C/T); INV-02 channel space is exactly 96 distinct labels; INV-03
    //   Σ counts = #classifiable variants; INV-04 folding is reverse-complement
    //   (a purine-ref variant and its pyrimidine-strand form share a channel).
    //
    // Every expected value below is derived from independent oracles re-built
    // here from the doc (the pyrimidine fold rule §2.2, the 6×4×4 channel set
    // §2.3, and the group-by partition §2.4), never from the implementation.
    // =========================================================================

    /// <summary>The four DNA bases, used to build oracle generators and the 96-channel set independently.</summary>
    private static readonly char[] SbsBases = { 'A', 'C', 'G', 'T' };

    /// <summary>
    /// The six pyrimidine substitutions of the SBS-96 classification (doc §2.2): C&gt;A, C&gt;G, C&gt;T,
    /// T&gt;A, T&gt;C, T&gt;G. Declared independently of the implementation.
    /// </summary>
    private static readonly (char Ref, char Alt)[] SbsPyrimidineSubstitutions =
    {
        ('C', 'A'), ('C', 'G'), ('C', 'T'), ('T', 'A'), ('T', 'C'), ('T', 'G')
    };

    /// <summary>Watson-Crick complement A↔T, C↔G (doc §4.2). Independent oracle helper.</summary>
    private static char SbsComplement(char b) => b switch
    {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => throw new ArgumentOutOfRangeException(nameof(b), b, "non-ACGT base in oracle"),
    };

    /// <summary>
    /// INDEPENDENT pyrimidine-strand fold oracle (doc §2.2 / INV-01, INV-04). Upper-cases the bases, then:
    /// if ref ∈ {C,T} the mutation is already on the pyrimidine strand → <c>"{5'}[{ref}&gt;{alt}]{3'}"</c>;
    /// if ref ∈ {A,G} it is reverse-complemented → <c>"{comp(3')}[{comp(ref)}&gt;{comp(alt)}]{comp(5')}"</c>.
    /// This is recomputed from the doc, NOT copied from <see cref="OncologyAnalyzer.ClassifySbsContext"/>.
    /// </summary>
    private static string FoldChannelOracle(char five, char reference, char alt, char three)
    {
        char f = char.ToUpperInvariant(five);
        char r = char.ToUpperInvariant(reference);
        char a = char.ToUpperInvariant(alt);
        char t = char.ToUpperInvariant(three);

        if (r is 'A' or 'G')
        {
            return $"{SbsComplement(t)}[{SbsComplement(r)}>{SbsComplement(a)}]{SbsComplement(f)}";
        }

        return $"{f}[{r}>{a}]{t}";
    }

    /// <summary>
    /// INDEPENDENT 96-channel set oracle (doc §2.3 / INV-02): all <c>5'[s&gt;e]3'</c> over the six pyrimidine
    /// substitutions × 4 5'-bases × 4 3'-bases. Built here, not read from the implementation.
    /// </summary>
    private static HashSet<string> BuildExpected96ChannelSet()
    {
        var set = new HashSet<string>(StringComparer.Ordinal);
        foreach (var (r, a) in SbsPyrimidineSubstitutions)
        {
            foreach (char five in SbsBases)
            {
                foreach (char three in SbsBases)
                {
                    set.Add($"{five}[{r}>{a}]{three}");
                }
            }
        }

        return set;
    }

    /// <summary>Generator for one valid SNV (5', ref, alt, 3') over A/C/G/T with ref ≠ alt.</summary>
    private static Gen<(char Five, char Ref, char Alt, char Three)> SnvGen() =>
        from five in Gen.Elements(SbsBases)
        from reference in Gen.Elements(SbsBases)
        from altOffset in Gen.Choose(1, 3) // ref ≠ alt: pick one of the 3 other bases deterministically
        from three in Gen.Elements(SbsBases)
        let alt = SbsBases[(Array.IndexOf(SbsBases, reference) + altOffset) % 4]
        select (five, reference, alt, three);

    private static Arbitrary<(char Five, char Ref, char Alt, char Three)> SnvArbitrary() =>
        SnvGen().ToArbitrary();

    /// <summary>Generator for a multiset (list, length 0–40) of valid SNVs — drives the catalog properties.</summary>
    private static Arbitrary<(char Five, char Ref, char Alt, char Three)[]> SnvListArbitrary() =>
        (from n in Gen.Choose(0, 40)
         from arr in SnvGen().ArrayOf(n)
         select arr)
        .ToArbitrary();

    // -------------------------------------------------------------------------
    // R + INV-02 — exactly 96 distinct channels equal to the recomputed set
    // -------------------------------------------------------------------------

    /// <summary>
    /// R + INV-02: <see cref="OncologyAnalyzer.EnumerateSbs96Channels"/> has exactly 96 entries, all DISTINCT,
    /// and as a set equals the independently recomputed 6×4×4 channel set. Also pins
    /// <c>Sbs96ChannelCount == 96</c>. This is a deterministic fact, asserted as a single [Test].
    /// </summary>
    [Test]
    [Category("Property")]
    public void Sbs96_Enumerate_IsExactly96DistinctChannels_EqualToRecomputedSet()
    {
        var channels = OncologyAnalyzer.EnumerateSbs96Channels();
        var expected = BuildExpected96ChannelSet();
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.Sbs96ChannelCount, Is.EqualTo(96), "Sbs96ChannelCount");
            Assert.That(channels, Has.Count.EqualTo(96), "enumerated count");
            Assert.That(channels.Distinct(StringComparer.Ordinal).Count(), Is.EqualTo(96), "distinct count");
            Assert.That(expected, Has.Count.EqualTo(96), "oracle set size");
            Assert.That(
                new HashSet<string>(channels, StringComparer.Ordinal).SetEquals(expected),
                Is.True,
                "enumerated channels must equal the recomputed 6×4×4 set");
        });
    }

    /// <summary>
    /// R + INV-02: every <see cref="OncologyAnalyzer.ClassifySbsContext"/> result for an arbitrary valid SNV is
    /// a member of the independently recomputed 96-channel set.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Sbs96_Classify_AlwaysYieldsAChannelInThe96Set()
    {
        var expected = BuildExpected96ChannelSet();
        return Prop.ForAll(SnvArbitrary(), snv =>
        {
            string channel = OncologyAnalyzer.ClassifySbsContext(snv.Five, snv.Ref, snv.Alt, snv.Three);
            return expected.Contains(channel)
                .Label($"{snv.Five}[{snv.Ref}>{snv.Alt}]{snv.Three} → {channel} not in 96-set");
        });
    }

    // -------------------------------------------------------------------------
    // P + INV-01 / INV-04 — each SNV folds to one pyrimidine-centred channel
    // -------------------------------------------------------------------------

    /// <summary>
    /// P + INV-01/INV-04: <see cref="OncologyAnalyzer.ClassifySbsContext"/> equals the INDEPENDENT pyrimidine
    /// fold oracle for every valid (5', ref, alt, 3'); the reference base in the result (the char before '&gt;')
    /// is ALWAYS a pyrimidine C/T (INV-01).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Sbs96_Classify_MatchesFoldOracle_AndRefIsPyrimidine()
    {
        return Prop.ForAll(SnvArbitrary(), snv =>
        {
            string actual = OncologyAnalyzer.ClassifySbsContext(snv.Five, snv.Ref, snv.Alt, snv.Three);
            string expected = FoldChannelOracle(snv.Five, snv.Ref, snv.Alt, snv.Three);
            char refBase = actual[2]; // "5'[R>A]3'": index 2 is R
            bool pyrimidine = refBase is 'C' or 'T';
            return (actual == expected && pyrimidine)
                .Label($"{snv.Five}[{snv.Ref}>{snv.Alt}]{snv.Three}: actual={actual}, oracle={expected}, refPyr={pyrimidine}");
        });
    }

    /// <summary>
    /// P (case-insensitivity, doc §3.1 / §6.1): lower-cased inputs classify identically to upper-cased inputs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Sbs96_Classify_IsCaseInsensitive()
    {
        return Prop.ForAll(SnvArbitrary(), snv =>
        {
            string upper = OncologyAnalyzer.ClassifySbsContext(snv.Five, snv.Ref, snv.Alt, snv.Three);
            string lower = OncologyAnalyzer.ClassifySbsContext(
                char.ToLowerInvariant(snv.Five),
                char.ToLowerInvariant(snv.Ref),
                char.ToLowerInvariant(snv.Alt),
                char.ToLowerInvariant(snv.Three));
            return (upper == lower).Label($"upper={upper} != lower={lower}");
        });
    }

    /// <summary>
    /// INV-04 (strand equivalence): a variant and its reverse-complement strand form map to the SAME channel —
    /// <c>Classify(5',ref,alt,3') == Classify(comp(3'),comp(ref),comp(alt),comp(5'))</c>, over arbitrary valid SNVs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Sbs96_Classify_FoldsRevCompToSameChannel_INV04()
    {
        return Prop.ForAll(SnvArbitrary(), snv =>
        {
            char fU = char.ToUpperInvariant(snv.Five);
            char rU = char.ToUpperInvariant(snv.Ref);
            char aU = char.ToUpperInvariant(snv.Alt);
            char tU = char.ToUpperInvariant(snv.Three);

            string forward = OncologyAnalyzer.ClassifySbsContext(fU, rU, aU, tU);
            string revComp = OncologyAnalyzer.ClassifySbsContext(
                SbsComplement(tU), SbsComplement(rU), SbsComplement(aU), SbsComplement(fU));
            return (forward == revComp)
                .Label($"{fU}[{rU}>{aU}]{tU} → {forward} ; revcomp → {revComp}");
        });
    }

    // -------------------------------------------------------------------------
    // P (Σ counts = #SNVs) + INV-03 — catalog partition / group-by oracle
    // -------------------------------------------------------------------------

    /// <summary>
    /// P + INV-03: <see cref="OncologyAnalyzer.Build96ContextCatalog"/> has exactly 96 keys equal to the
    /// enumerated channel set (all present, zero-count included), Σ values == #variants, and each channel's
    /// count equals the number of input variants that fold to it under the INDEPENDENT group-by oracle.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Sbs96_Catalog_MatchesGroupByOracle_AndSumsToVariantCount()
    {
        var channelSet = BuildExpected96ChannelSet();
        return Prop.ForAll(SnvListArbitrary(), variants =>
        {
            var catalog = OncologyAnalyzer.Build96ContextCatalog(
                variants.Select(v => (v.Five, v.Ref, v.Alt, v.Three)));

            // Independent group-by oracle: fold each variant, tally per channel.
            var oracle = channelSet.ToDictionary(c => c, _ => 0, StringComparer.Ordinal);
            foreach (var v in variants)
            {
                oracle[FoldChannelOracle(v.Five, v.Ref, v.Alt, v.Three)]++;
            }

            bool keysMatch = catalog.Count == 96
                             && new HashSet<string>(catalog.Keys, StringComparer.Ordinal).SetEquals(channelSet);
            bool sumMatches = catalog.Values.Sum() == variants.Length;
            bool entriesMatch = oracle.All(kv => catalog[kv.Key] == kv.Value);

            return (keysMatch && sumMatches && entriesMatch)
                .Label($"keys={keysMatch}, sum={catalog.Values.Sum()}/{variants.Length}, entries={entriesMatch}");
        });
    }

    /// <summary>
    /// P (edge, doc §6.1): an empty variant collection yields all 96 channels present, each with count 0.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Sbs96_Catalog_Empty_HasAll96ChannelsAtZero()
    {
        var catalog = OncologyAnalyzer.Build96ContextCatalog(
            Array.Empty<(char, char, char, char)>());
        var channelSet = BuildExpected96ChannelSet();
        Assert.Multiple(() =>
        {
            Assert.That(catalog, Has.Count.EqualTo(96));
            Assert.That(
                new HashSet<string>(catalog.Keys, StringComparer.Ordinal).SetEquals(channelSet),
                Is.True);
            Assert.That(catalog.Values, Is.All.EqualTo(0));
            Assert.That(catalog.Values.Sum(), Is.EqualTo(0));
        });
    }

    // -------------------------------------------------------------------------
    // D — determinism
    // -------------------------------------------------------------------------

    /// <summary>
    /// D: identical inputs ⇒ identical classification (same SNV classified twice yields the same channel).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Sbs96_Classify_IsDeterministic()
    {
        return Prop.ForAll(SnvArbitrary(), snv =>
        {
            string a = OncologyAnalyzer.ClassifySbsContext(snv.Five, snv.Ref, snv.Alt, snv.Three);
            string b = OncologyAnalyzer.ClassifySbsContext(snv.Five, snv.Ref, snv.Alt, snv.Three);
            return (a == b).Label($"{a} != {b}");
        });
    }

    /// <summary>
    /// D: identical variant collections ⇒ identical catalogs (entrywise over all 96 channels).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Sbs96_Catalog_IsDeterministic()
    {
        return Prop.ForAll(SnvListArbitrary(), variants =>
        {
            var c1 = OncologyAnalyzer.Build96ContextCatalog(variants.Select(v => (v.Five, v.Ref, v.Alt, v.Three)));
            var c2 = OncologyAnalyzer.Build96ContextCatalog(variants.Select(v => (v.Five, v.Ref, v.Alt, v.Three)));
            bool equal = c1.Count == c2.Count && c1.All(kv => c2.TryGetValue(kv.Key, out int v) && v == kv.Value);
            return equal.Label("catalogs differ entrywise");
        });
    }

    // -------------------------------------------------------------------------
    // Validation / edge cases (doc §3.3, §6.1)
    // -------------------------------------------------------------------------

    /// <summary>Validation: a non-ACGT base in ANY of the four positions raises <see cref="ArgumentException"/>.</summary>
    [TestCase('X', 'C', 'A', 'A')]
    [TestCase('A', 'X', 'A', 'A')]
    [TestCase('A', 'C', 'X', 'A')]
    [TestCase('A', 'C', 'A', 'X')]
    [TestCase('N', 'C', 'A', 'A')]
    [TestCase('1', 'C', 'A', 'A')]
    [Category("Property")]
    public void Sbs96_Classify_NonAcgtBase_Throws(char five, char reference, char alt, char three)
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ClassifySbsContext(five, reference, alt, three));
    }

    /// <summary>Validation (doc §3.3): reference == alternate is not a substitution ⇒ <see cref="ArgumentException"/>.</summary>
    [TestCase('A', 'C', 'C', 'A')]
    [TestCase('T', 'G', 'G', 'A')]
    [TestCase('A', 'A', 'A', 'A')]
    [Category("Property")]
    public void Sbs96_Classify_RefEqualsAlt_Throws(char five, char reference, char alt, char three)
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ClassifySbsContext(five, reference, alt, three));
    }

    /// <summary>Validation (doc §3.3): <see cref="OncologyAnalyzer.Build96ContextCatalog"/>(null) ⇒ <see cref="ArgumentNullException"/>.</summary>
    [Test]
    [Category("Property")]
    public void Sbs96_Catalog_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.Build96ContextCatalog(null!));
    }

    /// <summary>Anchor (doc §7.1): a G&gt;T at 5'-TGA-3' folds to <c>T[C&gt;A]A</c>; lower-case input is identical.</summary>
    [TestCase('T', 'G', 'T', 'A', "T[C>A]A")]
    [TestCase('t', 'g', 't', 'a', "T[C>A]A")]
    [TestCase('A', 'C', 'A', 'A', "A[C>A]A")]
    [Category("Property")]
    public void Sbs96_Classify_WorkedExample_Section71(char five, char reference, char alt, char three, string expected)
    {
        Assert.That(
            OncologyAnalyzer.ClassifySbsContext(five, reference, alt, three),
            Is.EqualTo(expected));
    }

    /// <summary>
    /// Anchor (doc §7.1 catalog walk-through): A[C&gt;A]A and the folded T[C&gt;A]A each count 1, all other
    /// 94 channels count 0; total = 2.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Sbs96_Catalog_WorkedExample_Section71()
    {
        var catalog = OncologyAnalyzer.Build96ContextCatalog(new[]
        {
            ('A', 'C', 'A', 'A'), // A[C>A]A (pyrimidine, unchanged)
            ('T', 'G', 'T', 'A'), // folds to T[C>A]A
        });
        Assert.Multiple(() =>
        {
            Assert.That(catalog, Has.Count.EqualTo(96));
            Assert.That(catalog["A[C>A]A"], Is.EqualTo(1));
            Assert.That(catalog["T[C>A]A"], Is.EqualTo(1));
            Assert.That(catalog.Values.Sum(), Is.EqualTo(2));
            Assert.That(
                catalog.Where(kv => kv.Key is not ("A[C>A]A" or "T[C>A]A")).All(kv => kv.Value == 0),
                Is.True,
                "all other 94 channels must be 0");
        });
    }

    #endregion

    #region ONCO-SIG-002 — Mutational Signature Fitting (NNLS Refitting + Cosine Similarity)

    // -------------------------------------------------------------------------
    // Independent oracles (derived from the doc, NOT from the implementation):
    //   - cosine(a,b) = dot(a,b) / (‖a‖·‖b‖); 0 when either norm is 0  (doc §2.2, INV-01)
    //   - reconstruction[k] = Σ_j signatures[j][k]·exposures[j]        (doc §2.2, R = S·x)
    //   - residual error E(x) = ‖S·x − d‖₂²                            (doc §2.2 NNLS objective)
    // These are recomputed here from scratch so the tests assert the theory,
    // not the code's own output.
    // -------------------------------------------------------------------------

    private const double SigClosedFormTolerance = 1e-9;   // closed-form algebra
    private const double SigSolverRecoveryTolerance = 1e-6; // NNLS recovery / optimality

    /// <summary>Independent cosine oracle: dot(a,b)/(‖a‖·‖b‖); 0.0 when either Euclidean norm is 0.</summary>
    private static double OracleCosine(IReadOnlyList<double> a, IReadOnlyList<double> b)
    {
        double dot = 0.0, na = 0.0, nb = 0.0;
        for (int i = 0; i < a.Count; i++)
        {
            dot += a[i] * b[i];
            na += a[i] * a[i];
            nb += b[i] * b[i];
        }

        if (na == 0.0 || nb == 0.0)
        {
            return 0.0;
        }

        return dot / (Math.Sqrt(na) * Math.Sqrt(nb));
    }

    /// <summary>Independent linear-combination oracle for the reconstruction S·x.</summary>
    private static double[] OracleReconstruct(IReadOnlyList<IReadOnlyList<double>> signatures, IReadOnlyList<double> x)
    {
        int channels = signatures[0].Count;
        var r = new double[channels];
        for (int j = 0; j < signatures.Count; j++)
        {
            for (int k = 0; k < channels; k++)
            {
                r[k] += signatures[j][k] * x[j];
            }
        }

        return r;
    }

    /// <summary>Independent residual oracle E(x) = ‖S·x − d‖₂² (the NNLS objective).</summary>
    private static double OracleResidual(
        IReadOnlyList<IReadOnlyList<double>> signatures, IReadOnlyList<double> x, IReadOnlyList<double> d)
    {
        double[] r = OracleReconstruct(signatures, x);
        double sum = 0.0;
        for (int k = 0; k < r.Length; k++)
        {
            double diff = r[k] - d[k];
            sum += diff * diff;
        }

        return sum;
    }

    // ---- generators -----------------------------------------------------------

    /// <summary>Non-negative finite double in a sane magnitude (0 .. 100, milli-resolution).</summary>
    private static Gen<double> NonNegDoubleGen() => Gen.Choose(0, 100_000).Select(v => v / 1000.0);

    /// <summary>A non-negative double vector of the given length.</summary>
    private static Gen<double[]> NonNegVectorGen(int length) => NonNegDoubleGen().ArrayOf(length);

    /// <summary>A pair of non-negative vectors sharing a common (small) length, length ≥ 1.</summary>
    private static Arbitrary<(double[] a, double[] b)> NonNegVectorPairArbitrary() =>
        (from n in Gen.Choose(1, 8)
         from a in NonNegVectorGen(n)
         from b in NonNegVectorGen(n)
         select (a, b)).ToArbitrary();

    /// <summary>A single non-zero non-negative vector (forces at least one positive entry).</summary>
    private static Arbitrary<double[]> NonZeroVectorArbitrary() =>
        (from n in Gen.Choose(1, 8)
         from baseVec in NonNegVectorGen(n)
         from idx in Gen.Choose(0, n - 1)
         from bump in Gen.Choose(1, 100_000).Select(v => v / 1000.0)
         select Bump(baseVec, idx, bump)).ToArbitrary();

    private static double[] Bump(double[] v, int idx, double bump)
    {
        var c = (double[])v.Clone();
        c[idx] += bump; // guarantees a strictly positive entry ⇒ non-zero norm
        return c;
    }

    /// <summary>
    /// A general fitting problem: n channels, k signatures (non-negative finite entries) and a
    /// non-negative catalog of length n — for invariants that must hold on ARBITRARY problems.
    /// </summary>
    private static Arbitrary<(double[][] signatures, double[] catalog)> FitProblemArbitrary() =>
        (from n in Gen.Choose(1, 6)
         from k in Gen.Choose(1, 4)
         from sigs in NonNegVectorGen(n).ArrayOf(k)
         from catalog in NonNegVectorGen(n)
         select (sigs, catalog)).ToArbitrary();

    /// <summary>
    /// A CONSTRUCTIBLE problem with a known answer: well-separated (standard-basis-scaled) signatures
    /// S and a known non-negative exposure x, with d = S·x. Because each signature occupies a distinct
    /// channel block scaled positively, S has full column rank and the unique NNLS optimum equals x.
    /// </summary>
    private static Arbitrary<(double[][] signatures, double[] exposures, double[] catalog)> KnownAnswerArbitrary() =>
        (from k in Gen.Choose(1, 4)
         from scales in Gen.Choose(1, 50).Select(v => v / 10.0).ArrayOf(k) // strictly positive scales
         from xRaw in Gen.Choose(0, 50_000).Select(v => v / 1000.0).ArrayOf(k)
         select BuildKnownAnswer(scales, xRaw)).ToArbitrary();

    private static (double[][] signatures, double[] exposures, double[] catalog) BuildKnownAnswer(
        double[] scales, double[] x)
    {
        int k = scales.Length;
        // Signature j is scales[j]·e_j over k channels ⇒ S diagonal, full rank, columns orthogonal.
        var sigs = new double[k][];
        for (int j = 0; j < k; j++)
        {
            sigs[j] = new double[k];
            sigs[j][j] = scales[j];
        }

        double[] d = OracleReconstruct(sigs, x);
        return (sigs, x, d);
    }

    private static IReadOnlyList<IReadOnlyList<double>> AsSignatures(double[][] sigs) =>
        sigs.Select(s => (IReadOnlyList<double>)s).ToArray();

    // ===== CosineSimilarity ===================================================

    /// <summary>
    /// INV-01: <c>CosineSimilarity(a,b) ∈ [0,1]</c> for non-negative vectors and equals the independent
    /// dot/norm oracle within 1e-9; a zero-norm vector yields exactly 0.0. (doc §2.2, §2.4 INV-01, §6.1)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CosineSimilarity_MatchesDotNormOracle_AndIsInUnitRange()
    {
        return Prop.ForAll(NonNegVectorPairArbitrary(), pair =>
        {
            double actual = OncologyAnalyzer.CosineSimilarity(pair.a, pair.b);
            double oracle = OracleCosine(pair.a, pair.b);
            bool matches = Math.Abs(actual - oracle) < SigClosedFormTolerance;
            // For non-negative inputs the cosine is in [0,1] (allow fp slack at the bounds).
            bool inRange = actual >= -SigClosedFormTolerance && actual <= 1.0 + SigClosedFormTolerance;
            return (matches && inRange).Label($"cos {actual} vs oracle {oracle}, inRange={inRange}");
        });
    }

    /// <summary>
    /// INV-02: <c>CosineSimilarity(a,a) = 1</c> (within 1e-9) for any non-zero vector a. (doc §2.4 INV-02)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CosineSimilarity_OfVectorWithItself_IsOne()
    {
        return Prop.ForAll(NonZeroVectorArbitrary(), a =>
        {
            double self = OncologyAnalyzer.CosineSimilarity(a, a);
            return (Math.Abs(self - 1.0) < SigClosedFormTolerance).Label($"cos(a,a) = {self} ≠ 1");
        });
    }

    /// <summary>
    /// INV-03 (scale invariance): <c>CosineSimilarity(a, k·b) = CosineSimilarity(a,b)</c> for any k &gt; 0.
    /// The cosine of an angle is invariant under positive scaling. (doc §2.4 INV-03)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CosineSimilarity_IsInvariantUnderPositiveScaling()
    {
        var arb = (from pair in NonNegVectorPairArbitrary().Generator
                   from kMilli in Gen.Choose(1, 100_000)
                   select (pair.a, pair.b, k: kMilli / 1000.0)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            double baseCos = OncologyAnalyzer.CosineSimilarity(t.a, t.b);
            double[] scaled = t.b.Select(x => x * t.k).ToArray();
            double scaledCos = OncologyAnalyzer.CosineSimilarity(t.a, scaled);
            return (Math.Abs(baseCos - scaledCos) < SigClosedFormTolerance)
                .Label($"cos(a,b)={baseCos} ≠ cos(a,{t.k}·b)={scaledCos}");
        });
    }

    /// <summary>
    /// Anchor: orthogonal non-negative vectors (disjoint support) ⇒ cosine 0; a zero vector ⇒ 0.0.
    /// (doc §6.1 orthogonal/zero-norm cases)
    /// </summary>
    [Test]
    [Category("Property")]
    public void CosineSimilarity_OrthogonalAndZeroVectors_AreZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CosineSimilarity(new double[] { 1, 0 }, new double[] { 0, 1 }),
                Is.EqualTo(0.0).Within(SigClosedFormTolerance), "orthogonal ⇒ 0");
            Assert.That(OncologyAnalyzer.CosineSimilarity(new double[] { 0, 0 }, new double[] { 3, 4 }),
                Is.EqualTo(0.0), "zero-norm ⇒ 0.0");
            Assert.That(OncologyAnalyzer.CosineSimilarity(new double[] { 3, 4 }, new double[] { 0, 0 }),
                Is.EqualTo(0.0), "zero-norm ⇒ 0.0");
        });
    }

    // ===== ReconstructCatalog =================================================

    /// <summary>
    /// Reconstruction oracle: <c>Reconstruction[k] = Σ_j signatures[j][k]·exposures[j]</c> within 1e-9,
    /// recomputed independently. (doc §2.2 R = S·x)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReconstructCatalog_EqualsLinearCombinationOracle()
    {
        var arb = (from n in Gen.Choose(1, 6)
                   from k in Gen.Choose(1, 4)
                   from sigs in NonNegVectorGen(n).ArrayOf(k)
                   from exposures in NonNegVectorGen(k)
                   select (sigs, exposures)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            IReadOnlyList<double> actual = OncologyAnalyzer.ReconstructCatalog(AsSignatures(t.sigs), t.exposures);
            double[] oracle = OracleReconstruct(AsSignatures(t.sigs), t.exposures);
            bool ok = actual.Count == oracle.Length;
            for (int k = 0; ok && k < oracle.Length; k++)
            {
                ok &= Math.Abs(actual[k] - oracle[k]) < SigClosedFormTolerance;
            }

            return ok.Label("ReconstructCatalog ≠ Σ_j signatures[j][k]·exposures[j]");
        });
    }

    // ===== FitSignatures ======================================================

    /// <summary>
    /// R + INV-04: every fitted exposure is non-negative (NNLS constraint x ≥ 0; allow −1e-9 fp slack).
    /// (doc §2.4 INV-04, checklist R: exposures ≥ 0)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FitSignatures_Exposures_AreNonNegative()
    {
        return Prop.ForAll(FitProblemArbitrary(), t =>
        {
            var fit = OncologyAnalyzer.FitSignatures(t.catalog, AsSignatures(t.signatures));
            return fit.Exposures.All(e => e >= -SigClosedFormTolerance)
                .Label($"negative exposure in [{string.Join(",", fit.Exposures)}]");
        });
    }

    /// <summary>
    /// P + INV-06: <c>NormalizedExposures = Exposures / Σ Exposures</c> (sum to 1 within 1e-9) when the
    /// total is positive; otherwise all zero. The doc-derived "Σ proportions = 1" form, NOT the loose
    /// "Σ exposures = #mutations". (doc §2.4 INV-06; checklist P: normalised 1.0)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FitSignatures_NormalizedExposures_AreExposuresOverTheirSum()
    {
        return Prop.ForAll(FitProblemArbitrary(), t =>
        {
            var fit = OncologyAnalyzer.FitSignatures(t.catalog, AsSignatures(t.signatures));
            double sum = fit.Exposures.Sum();

            bool ok;
            if (sum > 0.0)
            {
                ok = Math.Abs(fit.NormalizedExposures.Sum() - 1.0) < SigClosedFormTolerance;
                for (int j = 0; ok && j < fit.Exposures.Count; j++)
                {
                    ok &= Math.Abs(fit.NormalizedExposures[j] - fit.Exposures[j] / sum) < SigClosedFormTolerance;
                }
            }
            else
            {
                ok = fit.NormalizedExposures.All(p => p == 0.0);
            }

            return ok.Label($"normalized {string.Join(",", fit.NormalizedExposures)} (Σexp={sum})");
        });
    }

    /// <summary>
    /// Reconstruction + quality: the result's <c>Reconstruction</c> equals
    /// <c>ReconstructCatalog(signatures, Exposures)</c> and <c>ReconstructionCosineSimilarity</c> equals
    /// <c>CosineSimilarity(catalog, Reconstruction)</c>, both within 1e-9. (doc §2.2/§3.2)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FitSignatures_ReconstructionAndCosine_AreSelfConsistent()
    {
        return Prop.ForAll(FitProblemArbitrary(), t =>
        {
            var sigs = AsSignatures(t.signatures);
            var fit = OncologyAnalyzer.FitSignatures(t.catalog, sigs);

            double[] oracleRecon = OracleReconstruct(sigs, fit.Exposures);
            bool reconOk = fit.Reconstruction.Count == oracleRecon.Length;
            for (int k = 0; reconOk && k < oracleRecon.Length; k++)
            {
                reconOk &= Math.Abs(fit.Reconstruction[k] - oracleRecon[k]) < SigClosedFormTolerance;
            }

            double oracleCos = OracleCosine(t.catalog, fit.Reconstruction);
            bool cosOk = Math.Abs(fit.ReconstructionCosineSimilarity - oracleCos) < SigClosedFormTolerance;

            return (reconOk && cosOk)
                .Label($"recon/cos mismatch: cos={fit.ReconstructionCosineSimilarity} vs {oracleCos}");
        });
    }

    /// <summary>
    /// Known-answer recovery (standard basis): with signature j = e_j (S = identity over n channels),
    /// <c>FitSignatures(catalog, basis).Exposures == catalog</c> exactly and the reconstruction cosine
    /// is 1 for a non-zero catalog. The identity decomposition is unique. (doc §7.1; INV-07)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FitSignatures_StandardBasis_RecoversCatalogExactly()
    {
        return Prop.ForAll(NonZeroVectorArbitrary(), catalog =>
        {
            int n = catalog.Length;
            var basis = new double[n][];
            for (int j = 0; j < n; j++)
            {
                basis[j] = new double[n];
                basis[j][j] = 1.0; // e_j
            }

            var fit = OncologyAnalyzer.FitSignatures(catalog, AsSignatures(basis));

            bool exposuresOk = fit.Exposures.Count == n;
            for (int j = 0; exposuresOk && j < n; j++)
            {
                exposuresOk &= Math.Abs(fit.Exposures[j] - catalog[j]) < SigClosedFormTolerance;
            }

            bool cosOk = Math.Abs(fit.ReconstructionCosineSimilarity - 1.0) < SigClosedFormTolerance;
            return (exposuresOk && cosOk)
                .Label($"basis recovery failed: exposures={string.Join(",", fit.Exposures)} cos={fit.ReconstructionCosineSimilarity}");
        });
    }

    /// <summary>
    /// Known-answer recovery (well-conditioned S, d = S·x): with full-column-rank (diagonal-scaled)
    /// signatures and a known non-negative x, fitting d = S·x recovers x within ~1e-6 and cosine ≈ 1.
    /// The unconstrained LS optimum is already non-negative, so the NNLS optimum equals it (INV-07).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FitSignatures_WellConditionedKnownX_RecoversExposures()
    {
        return Prop.ForAll(KnownAnswerArbitrary(), t =>
        {
            var sigs = AsSignatures(t.signatures);
            var fit = OncologyAnalyzer.FitSignatures(t.catalog, sigs);

            bool xOk = fit.Exposures.Count == t.exposures.Length;
            for (int j = 0; xOk && j < t.exposures.Length; j++)
            {
                xOk &= Math.Abs(fit.Exposures[j] - t.exposures[j]) < SigSolverRecoveryTolerance;
            }

            // Cosine ≈ 1 only when d ≠ 0 (a zero x gives a zero catalog ⇒ cosine 0 by convention).
            double dNorm = Math.Sqrt(t.catalog.Sum(v => v * v));
            bool cosOk = dNorm == 0.0
                ? fit.ReconstructionCosineSimilarity == 0.0
                : Math.Abs(fit.ReconstructionCosineSimilarity - 1.0) < SigSolverRecoveryTolerance;

            return (xOk && cosOk)
                .Label($"recovered {string.Join(",", fit.Exposures)} vs x {string.Join(",", t.exposures)} cos={fit.ReconstructionCosineSimilarity}");
        });
    }

    /// <summary>
    /// M ("better fit → lower reconstruction error") = NNLS optimality (INV-05/INV-07). The fitted
    /// solution's residual ‖S·x_nnls − d‖² is (a) ≤ ‖d‖² because x = 0 is feasible, and (b) no worse
    /// than the residual of ANY randomly generated non-negative exposure vector y (within 1e-6):
    /// a strictly better fit can never have a higher reconstruction error. (doc §2.4 INV-05/INV-07)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FitSignatures_Residual_IsNoWorseThanAnyFeasibleExposure()
    {
        var arb = (from problem in FitProblemArbitrary().Generator
                   from y in NonNegVectorGen(problem.signatures.Length)
                   select (problem.signatures, problem.catalog, y)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            var sigs = AsSignatures(t.signatures);
            var fit = OncologyAnalyzer.FitSignatures(t.catalog, sigs);

            double residualNnls = OracleResidual(sigs, fit.Exposures, t.catalog);
            double dNormSquared = t.catalog.Sum(v => v * v);
            double residualY = OracleResidual(sigs, t.y, t.catalog);

            // (a) x = 0 feasible ⇒ optimum residual ≤ ‖d‖²; (b) optimum ≤ any feasible y's residual.
            bool boundedByD = residualNnls <= dNormSquared + SigSolverRecoveryTolerance;
            bool optimal = residualNnls <= residualY + SigSolverRecoveryTolerance;

            return (boundedByD && optimal)
                .Label($"nnls residual {residualNnls} > ‖d‖²={dNormSquared} or > feasible-y residual {residualY}");
        });
    }

    /// <summary>
    /// D (determinism): identical inputs produce an identical <see cref="SignatureFitResult"/> — every
    /// field equal entrywise. (checklist D: deterministic)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FitSignatures_IsDeterministic()
    {
        return Prop.ForAll(FitProblemArbitrary(), t =>
        {
            var sigs1 = AsSignatures(t.signatures);
            var sigs2 = AsSignatures(t.signatures);
            var a = OncologyAnalyzer.FitSignatures(t.catalog, sigs1);
            var b = OncologyAnalyzer.FitSignatures(t.catalog, sigs2);

            bool ok = a.Exposures.SequenceEqual(b.Exposures)
                      && a.NormalizedExposures.SequenceEqual(b.NormalizedExposures)
                      && a.Reconstruction.SequenceEqual(b.Reconstruction)
                      && a.ReconstructionCosineSimilarity.Equals(b.ReconstructionCosineSimilarity);

            return ok.Label("FitSignatures is not deterministic for identical inputs");
        });
    }

    // ===== Validation / edge cases ============================================

    /// <summary>Validation: null vectors to <c>CosineSimilarity</c> throw <see cref="ArgumentNullException"/>.</summary>
    [Test]
    [Category("Property")]
    public void CosineSimilarity_NullArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CosineSimilarity(null!, new double[] { 1 }));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CosineSimilarity(new double[] { 1 }, null!));
        });
    }

    /// <summary>Validation: empty or length-mismatched cosine vectors throw <see cref="ArgumentException"/>.</summary>
    [Test]
    [Category("Property")]
    public void CosineSimilarity_EmptyOrMismatchedVectors_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() =>
                OncologyAnalyzer.CosineSimilarity(Array.Empty<double>(), Array.Empty<double>()));
            Assert.Throws<ArgumentException>(() =>
                OncologyAnalyzer.CosineSimilarity(new double[] { 1, 2 }, new double[] { 1 }));
        });
    }

    /// <summary>Validation: null catalog/signatures (or a null signature vector) throw <see cref="ArgumentNullException"/>.</summary>
    [Test]
    [Category("Property")]
    public void FitSignatures_NullArguments_Throw()
    {
        var sigs = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.FitSignatures(null!, sigs));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.FitSignatures(new double[] { 1, 1 }, null!));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.ReconstructCatalog(null!, new double[] { 1 }));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.ReconstructCatalog(sigs, null!));
        });
    }

    /// <summary>
    /// Validation: empty, ragged, or dimension-mismatched inputs throw <see cref="ArgumentException"/>.
    /// (doc §3.3)
    /// </summary>
    [Test]
    [Category("Property")]
    public void FitSignatures_EmptyRaggedOrMismatched_Throw()
    {
        var ragged = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1, 2 } };
        var good = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        Assert.Multiple(() =>
        {
            // No signatures.
            Assert.Throws<ArgumentException>(() =>
                OncologyAnalyzer.FitSignatures(new double[] { 1 }, Array.Empty<IReadOnlyList<double>>()));
            // Ragged signature matrix.
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.FitSignatures(new double[] { 1, 0 }, ragged));
            // Catalog length ≠ channel count.
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.FitSignatures(new double[] { 1, 2, 3 }, good));
            // ReconstructCatalog: exposure count ≠ signature count.
            Assert.Throws<ArgumentException>(() =>
                OncologyAnalyzer.ReconstructCatalog(good, new double[] { 1, 2, 3 }));
        });
    }

    /// <summary>
    /// Edge (doc §6.1): a zero catalog d = 0 yields all-zero exposures, all-zero normalized exposures,
    /// and a zero reconstruction (the only feasible minimiser of ‖S·x‖², x ≥ 0).
    /// </summary>
    [Test]
    [Category("Property")]
    public void FitSignatures_ZeroCatalog_YieldsAllZeros()
    {
        var sigs = new IReadOnlyList<double>[] { new double[] { 1, 2 }, new double[] { 3, 1 } };
        var fit = OncologyAnalyzer.FitSignatures(new double[] { 0, 0 }, sigs);
        Assert.Multiple(() =>
        {
            Assert.That(fit.Exposures, Is.All.EqualTo(0.0));
            Assert.That(fit.NormalizedExposures, Is.All.EqualTo(0.0));
            Assert.That(fit.Reconstruction, Is.All.EqualTo(0.0));
        });
    }

    /// <summary>
    /// Anchor (doc §7.1 API example): catalog [3,5] with standard-basis signatures sig1=[1,0], sig2=[0,1]
    /// ⇒ exposures [3,5], reconstruction [3,5], cosine 1.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FitSignatures_Anchor_BasisCatalog35()
    {
        var sigs = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        var fit = OncologyAnalyzer.FitSignatures(new double[] { 3, 5 }, sigs);
        Assert.Multiple(() =>
        {
            Assert.That(fit.Exposures[0], Is.EqualTo(3.0).Within(SigClosedFormTolerance));
            Assert.That(fit.Exposures[1], Is.EqualTo(5.0).Within(SigClosedFormTolerance));
            Assert.That(fit.Reconstruction[0], Is.EqualTo(3.0).Within(SigClosedFormTolerance));
            Assert.That(fit.Reconstruction[1], Is.EqualTo(5.0).Within(SigClosedFormTolerance));
            Assert.That(fit.ReconstructionCosineSimilarity, Is.EqualTo(1.0).Within(SigClosedFormTolerance));
        });
    }

    /// <summary>
    /// Anchor (doc §7.1 walk-through): catalog [0,1] with sig1=[1,0], sig2=[1,1]. The unconstrained LS
    /// gives x₁ &lt; 0, so sig1 clamps to 0 and sig2 refits to 0.5 ⇒ exposures [0, 0.5],
    /// reconstruction [0.5, 0.5].
    /// </summary>
    [Test]
    [Category("Property")]
    public void FitSignatures_Anchor_ClampedNnls01()
    {
        var sigs = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 1, 1 } };
        var fit = OncologyAnalyzer.FitSignatures(new double[] { 0, 1 }, sigs);
        Assert.Multiple(() =>
        {
            Assert.That(fit.Exposures[0], Is.EqualTo(0.0).Within(SigSolverRecoveryTolerance));
            Assert.That(fit.Exposures[1], Is.EqualTo(0.5).Within(SigSolverRecoveryTolerance));
            Assert.That(fit.Reconstruction[0], Is.EqualTo(0.5).Within(SigSolverRecoveryTolerance));
            Assert.That(fit.Reconstruction[1], Is.EqualTo(0.5).Within(SigSolverRecoveryTolerance));
        });
    }

    #endregion

    #region ONCO-SIG-003 — Signature Exposure Bootstrap Confidence Intervals

    // -------------------------------------------------------------------------
    // Theory (Senkin 2021 MSA; Huang 2018; Efron 1979 percentile; Hyndman & Fan 1996 type-7):
    //   • Parametric bootstrap: each replicate resamples the integer catalog as a multinomial draw of
    //     N = Σ catalog mutations with pₖ = catalogₖ/N, then refits by NNLS.            (TestSpec §1.2.1–2)
    //   • Per-signature interval = [½(1−c), 1−½(1−c)] percentiles of the replicate exposures. (§1.2.3)
    //   • Point estimate = NNLS exposure of the OBSERVED (un-resampled) catalog.            (§1.2.5, INV-5)
    //   • Corner cases: N = 0 ⇒ all-zero; a single non-zero channel ⇒ deterministic collapse. (§1.3)
    //
    // GUARANTEED invariants are asserted on ARBITRARY problems: per-signature ordering (INV-2),
    // non-negativity (INV-1), determinism under a fixed seed (INV-4), the point-estimate contract
    // (INV-5), and percentile-bound monotonicity in the confidence level. The bracketing
    // "Lower ≤ point ≤ Upper" (checklist R) and "bootstrap mean ≈ point estimate" (checklist P)
    // are only mathematically EXACT on the deterministic-collapse construction (single non-zero
    // channel) and the N = 0 case — a percentile bracket of a skewed stochastic bootstrap need not
    // contain its own mean — so they are proven there as exact equalities, with one concrete
    // non-degenerate anchor (default 1000 replicates) verifying the statistical bracketing.
    // -------------------------------------------------------------------------

    private const double SigBootstrapTolerance = 1e-9;

    /// <summary>Reproducible RNG seed reused across the randomized bootstrap properties (INV-4).</summary>
    private const int SigBootstrapSeed = 42;

    /// <summary>Modest replicate count keeping randomized property runs fast yet non-degenerate.</summary>
    private const int SigBootstrapReplicates = 60;

    /// <summary>
    /// An arbitrary bootstrap problem: n channels (1..6), k signatures (non-negative finite vectors), and
    /// a non-negative INTEGER catalog of length n. Integer counts (not proportions) are required because
    /// the multinomial resample needs the total N = Σ catalog (TestSpec §1.4.2).
    /// </summary>
    private static Arbitrary<(double[][] signatures, int[] catalog)> BootstrapProblemArbitrary() =>
        (from n in Gen.Choose(1, 6)
         from k in Gen.Choose(1, 4)
         from sigs in NonNegVectorGen(n).ArrayOf(k)
         from catalog in Gen.Choose(0, 40).ArrayOf(n)
         select (sigs, catalog)).ToArbitrary();

    /// <summary>
    /// A deterministic-collapse problem: standard-basis signatures over n channels (so NNLS recovers the
    /// catalog exactly, exposure[j] = catalog[j]) and a catalog with exactly ONE positive channel. A
    /// multinomial draw of N over a single non-zero probability is deterministic ⇒ every resample equals
    /// the observed catalog ⇒ every replicate exposure equals the point estimate. (TestSpec §1.3)
    /// </summary>
    private static Arbitrary<(double[][] signatures, int[] catalog)> DeterministicCollapseArbitrary() =>
        (from n in Gen.Choose(1, 5)
         from hot in Gen.Choose(0, n - 1)
         from count in Gen.Choose(1, 500)
         select BuildCollapseProblem(n, hot, count)).ToArbitrary();

    private static (double[][] signatures, int[] catalog) BuildCollapseProblem(int n, int hot, int count)
    {
        // Signature j = e_j over n channels ⇒ S = identity ⇒ NNLS exposure[j] = catalog[j].
        var sigs = new double[n][];
        for (int j = 0; j < n; j++)
        {
            sigs[j] = new double[n];
            sigs[j][j] = 1.0;
        }

        var catalog = new int[n];
        catalog[hot] = count; // single non-zero channel ⇒ deterministic multinomial resample
        return (sigs, catalog);
    }

    /// <summary>
    /// INV-5 (contract): exactly one interval per signature, in signature order, and each interval's
    /// <c>PointEstimate</c> equals the NNLS exposure of the OBSERVED (un-resampled) catalog — recomputed
    /// independently via <see cref="OncologyAnalyzer.FitSignatures(IReadOnlyList{double}, IReadOnlyList{IReadOnlyList{double}})"/>
    /// (the documented point estimate, Senkin 2021; Huang 2018). The reported confidence is the one requested.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BootstrapExposures_PointEstimate_IsObservedNnlsFit_OneIntervalPerSignatureInOrder()
    {
        return Prop.ForAll(BootstrapProblemArbitrary(), t =>
        {
            var sigs = AsSignatures(t.signatures);
            double[] observed = t.catalog.Select(c => (double)c).ToArray();
            IReadOnlyList<double> fit = OncologyAnalyzer.FitSignatures(observed, sigs).Exposures;

            var intervals = OncologyAnalyzer.BootstrapExposures(
                t.catalog, sigs, SigBootstrapReplicates, OncologyAnalyzer.DefaultBootstrapConfidence, SigBootstrapSeed);

            bool ok = intervals.Count == t.signatures.Length;
            for (int j = 0; ok && j < intervals.Count; j++)
            {
                ok &= Math.Abs(intervals[j].PointEstimate - fit[j]) < SigBootstrapTolerance;
                ok &= Math.Abs(intervals[j].Confidence - OncologyAnalyzer.DefaultBootstrapConfidence) < SigBootstrapTolerance;
            }

            return ok.Label(
                $"point estimate ≠ observed NNLS fit, or interval count {intervals.Count} ≠ #signatures {t.signatures.Length}");
        });
    }

    /// <summary>
    /// INV-1 + INV-2: on an arbitrary catalog every interval is ordered (<c>Lower ≤ Upper</c>) and all of
    /// <c>Lower</c>, <c>Upper</c>, <c>Mean</c>, <c>PointEstimate</c> are non-negative — resampled counts are
    /// ≥ 0 and the NNLS constraint forces x ≥ 0, while the lower percentile probability ½(1−c) is strictly
    /// below the upper 1−½(1−c) (Efron 1979). (checklist R; TestSpec INV-1/INV-2)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BootstrapExposures_IntervalsOrdered_AndAllBoundsNonNegative()
    {
        return Prop.ForAll(BootstrapProblemArbitrary(), t =>
        {
            var intervals = OncologyAnalyzer.BootstrapExposures(
                t.catalog, AsSignatures(t.signatures), SigBootstrapReplicates, 0.95, SigBootstrapSeed);

            bool ok = true;
            foreach (var ci in intervals)
            {
                ok &= ci.Lower <= ci.Upper + SigBootstrapTolerance;
                ok &= ci.Lower >= -SigBootstrapTolerance;
                ok &= ci.Upper >= -SigBootstrapTolerance;
                ok &= ci.Mean >= -SigBootstrapTolerance;
                ok &= ci.PointEstimate >= -SigBootstrapTolerance;
            }

            return ok.Label("an interval was not ordered (Lower ≤ Upper) or a bound/mean/point was negative");
        });
    }

    /// <summary>
    /// Percentile-bound monotonicity in the confidence level: re-using the SAME replicate distribution
    /// (identical seed), a wider confidence c widens each interval — the lower bound moves down (or stays)
    /// and the upper bound moves up (or stays). Direct consequence of the type-7 quantile being monotone
    /// non-decreasing in its probability, with ½(1−c) decreasing and 1−½(1−c) increasing as c grows.
    /// (TestSpec S1; Hyndman &amp; Fan 1996)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BootstrapExposures_WiderConfidence_WidensOrKeepsInterval()
    {
        return Prop.ForAll(BootstrapProblemArbitrary(), t =>
        {
            var sigs = AsSignatures(t.signatures);
            var narrow = OncologyAnalyzer.BootstrapExposures(t.catalog, sigs, SigBootstrapReplicates, 0.50, SigBootstrapSeed);
            var wide = OncologyAnalyzer.BootstrapExposures(t.catalog, sigs, SigBootstrapReplicates, 0.95, SigBootstrapSeed);

            bool ok = narrow.Count == wide.Count;
            for (int j = 0; ok && j < narrow.Count; j++)
            {
                ok &= wide[j].Lower <= narrow[j].Lower + SigBootstrapTolerance; // lower bound non-increasing in c
                ok &= wide[j].Upper >= narrow[j].Upper - SigBootstrapTolerance; // upper bound non-decreasing in c
            }

            return ok.Label("a wider confidence level produced a strictly narrower interval");
        });
    }

    /// <summary>
    /// INV-4 (determinism): identical <c>(catalog, signatures, replicates, confidence, seed)</c> arguments
    /// produce element-wise identical intervals (every field equal). The fixed RNG seed makes the
    /// multinomial resampling reproducible. (TestSpec INV-4)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BootstrapExposures_SameArguments_AreDeterministic()
    {
        return Prop.ForAll(BootstrapProblemArbitrary(), t =>
        {
            var sigs = AsSignatures(t.signatures);
            var a = OncologyAnalyzer.BootstrapExposures(t.catalog, sigs, SigBootstrapReplicates, 0.95, SigBootstrapSeed);
            var b = OncologyAnalyzer.BootstrapExposures(t.catalog, sigs, SigBootstrapReplicates, 0.95, SigBootstrapSeed);

            bool ok = a.Count == b.Count;
            for (int j = 0; ok && j < a.Count; j++)
            {
                ok &= a[j].Equals(b[j]); // record struct ⇒ structural, field-wise equality
            }

            return ok.Label("BootstrapExposures is not deterministic for identical (…, seed) arguments");
        });
    }

    /// <summary>
    /// Checklist P + R, made EXACT (TestSpec §1.3 single-non-zero-channel collapse): with standard-basis
    /// signatures and a catalog whose single non-zero channel carries the whole mutation load, the
    /// multinomial resample is deterministic, so every replicate exposure equals the point estimate. Hence
    /// for every signature j the bootstrap <c>Mean</c>, <c>Lower</c> and <c>Upper</c> all equal the
    /// <c>PointEstimate</c> = <c>catalog[j]</c>, and trivially <c>Lower ≤ point ≤ Upper</c>.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BootstrapExposures_DeterministicCollapse_MeanAndBoundsEqualPointEstimate()
    {
        return Prop.ForAll(DeterministicCollapseArbitrary(), t =>
        {
            var intervals = OncologyAnalyzer.BootstrapExposures(
                t.catalog, AsSignatures(t.signatures), SigBootstrapReplicates, 0.95, SigBootstrapSeed);

            bool ok = intervals.Count == t.signatures.Length;
            for (int j = 0; ok && j < intervals.Count; j++)
            {
                double expected = t.catalog[j]; // basis recovery, constant across all resamples
                var ci = intervals[j];
                ok &= Math.Abs(ci.PointEstimate - expected) < SigBootstrapTolerance;
                ok &= Math.Abs(ci.Mean - expected) < SigBootstrapTolerance;
                ok &= Math.Abs(ci.Lower - expected) < SigBootstrapTolerance;
                ok &= Math.Abs(ci.Upper - expected) < SigBootstrapTolerance;
                ok &= ci.Lower <= ci.PointEstimate + SigBootstrapTolerance
                      && ci.PointEstimate <= ci.Upper + SigBootstrapTolerance;
            }

            return ok.Label("deterministic-collapse interval did not degenerate to the exact point estimate");
        });
    }

    /// <summary>
    /// Edge (TestSpec §1.3, N = 0): a zero-mutation catalog gives an all-zero resample every replicate, so
    /// every signature's interval is degenerate at 0 — point, mean, lower and upper all 0.
    /// </summary>
    [Test]
    [Category("Property")]
    public void BootstrapExposures_ZeroCatalog_AllIntervalFieldsZero()
    {
        var sigs = new IReadOnlyList<double>[] { new double[] { 1, 0, 0 }, new double[] { 0, 1, 0 } };
        var intervals = OncologyAnalyzer.BootstrapExposures(new[] { 0, 0, 0 }, sigs, replicates: 100, confidence: 0.95, seed: 42);

        Assert.That(intervals, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            foreach (var ci in intervals)
            {
                Assert.That(ci.PointEstimate, Is.EqualTo(0.0).Within(SigBootstrapTolerance));
                Assert.That(ci.Mean, Is.EqualTo(0.0).Within(SigBootstrapTolerance));
                Assert.That(ci.Lower, Is.EqualTo(0.0).Within(SigBootstrapTolerance));
                Assert.That(ci.Upper, Is.EqualTo(0.0).Within(SigBootstrapTolerance));
            }
        });
    }

    /// <summary>
    /// Statistical anchor (checklist R + P; TestSpec C1): on a concrete non-degenerate catalog with the
    /// default 1000 replicates, the percentile interval brackets the point estimate
    /// (<c>Lower ≤ point ≤ Upper</c>) for every signature, and the bootstrap mean is close to the point
    /// estimate. With standard-basis signatures the replicate exposure of channel j is the resampled count
    /// ~ Binomial(N, catalogⱼ/N), whose mean N·pⱼ equals the observed count (Senkin 2021).
    /// </summary>
    [Test]
    [Category("Property")]
    public void BootstrapExposures_NonDegenerateDefaultReplicates_BracketsPointEstimate()
    {
        var sigs = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        var catalog = new[] { 40, 10 }; // N = 50, signature 0 dominant
        var intervals = OncologyAnalyzer.BootstrapExposures(catalog, sigs); // defaults: 1000 reps, 0.95, seed 42

        Assert.That(intervals, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            foreach (var ci in intervals)
            {
                Assert.That(ci.Lower, Is.LessThanOrEqualTo(ci.PointEstimate + SigBootstrapTolerance),
                    "Lower bound must not exceed the point estimate.");
                Assert.That(ci.Upper, Is.GreaterThanOrEqualTo(ci.PointEstimate - SigBootstrapTolerance),
                    "Upper bound must not fall below the point estimate.");
            }

            Assert.That(intervals[0].Mean, Is.EqualTo(40.0).Within(1.5),
                "Bootstrap mean of the dominant signature ≈ its point estimate (N·p = 50·0.8).");
            Assert.That(intervals[1].Mean, Is.EqualTo(10.0).Within(1.5),
                "Bootstrap mean of the minor signature ≈ its point estimate (N·p = 50·0.2).");
        });
    }

    /// <summary>
    /// Validation: <c>replicates &lt; 1</c> or a confidence outside the open interval (0, 1) throws
    /// <see cref="ArgumentOutOfRangeException"/>. (TestSpec failure modes)
    /// </summary>
    [Test]
    [Category("Property")]
    public void BootstrapExposures_InvalidReplicatesOrConfidence_Throw()
    {
        var sigs = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.BootstrapExposures(new[] { 1, 1 }, sigs, replicates: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.BootstrapExposures(new[] { 1, 1 }, sigs, confidence: 0.0));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.BootstrapExposures(new[] { 1, 1 }, sigs, confidence: 1.0));
        });
    }

    /// <summary>
    /// Validation: a null catalog throws <see cref="ArgumentNullException"/>; a negative count or a catalog
    /// whose length differs from the signature channel count throws <see cref="ArgumentException"/>.
    /// (TestSpec failure modes)
    /// </summary>
    [Test]
    [Category("Property")]
    public void BootstrapExposures_NullNegativeOrMismatchedCatalog_Throw()
    {
        var sigs = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.BootstrapExposures(null!, sigs));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.BootstrapExposures(new[] { -1, 2 }, sigs));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.BootstrapExposures(new[] { 1, 2, 3 }, sigs));
        });
    }

    #endregion
}
