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
}
