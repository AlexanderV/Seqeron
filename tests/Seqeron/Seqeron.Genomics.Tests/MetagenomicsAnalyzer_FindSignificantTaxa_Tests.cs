// META-TAXA-001 — Significant Taxa Detection (Mann–Whitney U / Wilcoxon rank-sum)
// Evidence: docs/Evidence/META-TAXA-001-Evidence.md
// TestSpec: tests/TestSpecs/META-TAXA-001.md
// Source: Mann HB, Whitney DR (1947). Ann. Math. Statist. 18(1):50–60;
//         SciPy scipy.stats.mannwhitneyu documentation; Xia & Sun (2017) Genes & Diseases 4(3).

using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class MetagenomicsAnalyzer_FindSignificantTaxa_Tests
{
    // SciPy reference example (mannwhitneyu docs): U1=17, U2=3, σ=sqrt(200/12).
    private static readonly double[] ScipyX = { 19, 22, 16, 29, 24 };
    private static readonly double[] ScipyY = { 20, 11, 17, 12 };
    private const double ErfPTolerance = 1e-6; // A&S 7.1.26 erf approximation, |ε| ≤ 1.5e-7

    #region MannWhitneyU

    // M1 — SciPy mannwhitneyu(x=[19,22,16,29,24], y=[20,11,17,12]) → U1=17.0, U2=3.0.
    [Test]
    public void MannWhitneyU_ScipyExample_ReturnsU17And3()
    {
        var r = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY);

        Assert.Multiple(() =>
        {
            Assert.That(r.U1, Is.EqualTo(17.0).Within(1e-10),
                "SciPy reference: U1 for x=[19,22,16,29,24] is 17.0.");
            Assert.That(r.U2, Is.EqualTo(3.0).Within(1e-10),
                "SciPy reference: U2 = n1*n2 - U1 = 20 - 17 = 3.0.");
        });
    }

    // M2 — σ_U = sqrt(n1*n2*(n1+n2+1)/12) = sqrt(200/12); z (no cc) = (17-10)/σ = 1.7146428199482247.
    [Test]
    public void MannWhitneyU_ScipyExample_NoContinuity_MatchesSigmaAndZ()
    {
        var r = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY, useContinuityCorrection: false);

        // σ = sqrt(5*4*10/12) = 4.08248290463863 ; z = (max(17,3)-10)/σ
        const double expectedZ = 7.0 / 4.08248290463863; // = 1.7146428199482247
        Assert.That(r.Z, Is.EqualTo(expectedZ).Within(1e-10),
            "z = (|U - m_U|)/σ_U with m_U=10, σ_U=sqrt(200/12) per Mann & Whitney (1947).");
    }

    // M3 — SciPy asymptotic two-tailed p without continuity = 0.0864107329737.
    [Test]
    public void MannWhitneyU_ScipyExample_NoContinuity_MatchesPValue()
    {
        var r = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY, useContinuityCorrection: false);

        Assert.That(r.PValue, Is.EqualTo(0.0864107329737).Within(ErfPTolerance),
            "SciPy mannwhitneyu asymptotic p (no continuity) = 0.0864107329737; matched within erf approx tolerance.");
    }

    // M4 — SciPy default (continuity on): z=(7-0.5)/σ=1.5921683328090657, p=0.11134688653314041.
    [Test]
    public void MannWhitneyU_ScipyExample_WithContinuity_MatchesZAndPValue()
    {
        var r = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY); // default cc = true

        const double expectedZ = 6.5 / 4.08248290463863; // = 1.5921683328090657
        Assert.Multiple(() =>
        {
            Assert.That(r.Z, Is.EqualTo(expectedZ).Within(1e-10),
                "Continuity correction reduces |U - m_U| by 0.5 (SciPy default): z = (7-0.5)/σ.");
            Assert.That(r.PValue, Is.EqualTo(0.11134688653314041).Within(ErfPTolerance),
                "SciPy mannwhitneyu asymptotic p (with continuity) = 0.11134688653314041.");
        });
    }

    // M5 — Mann & Whitney tortoise/hare worked example: U_T=11, U_H=25 (ranks 12,6,5,4,3,2 vs rest).
    [Test]
    public void MannWhitneyU_TortoiseHare_ReturnsU11And25()
    {
        // Pooled finish order (rank 12 = first): T H H H H H T T T T T H
        // Tortoise positions (1-based from the rank-12 end) give rank sum R_T = 12+6+5+4+3+2 = 32.
        // Build values so the pooled sort reproduces those ranks. Larger value = better rank.
        // Assign distinct ascending values to ranks 1..12; tortoise occupies ranks {12,6,5,4,3,2}.
        var values = Enumerable.Range(1, 12).Select(i => (double)i).ToArray(); // value == rank
        int[] tortoiseRanks = { 12, 6, 5, 4, 3, 2 };
        int[] hareRanks = { 11, 10, 9, 8, 7, 1 };
        var tort = tortoiseRanks.Select(rk => values[rk - 1]).ToArray();
        var hare = hareRanks.Select(rk => values[rk - 1]).ToArray();

        var r = MetagenomicsAnalyzer.MannWhitneyU(tort, hare);

        Assert.Multiple(() =>
        {
            Assert.That(r.U1, Is.EqualTo(11.0).Within(1e-10),
                "U_T = R_T - n(n+1)/2 = 32 - 21 = 11 (Mann & Whitney tortoise/hare example).");
            Assert.That(r.U2, Is.EqualTo(25.0).Within(1e-10),
                "U_H = n1*n2 - U_T = 36 - 11 = 25.");
        });
    }

    // M6 — INV-01: U1 + U2 = n1*n2 for arbitrary inputs (property-based).
    [Test]
    public void MannWhitneyU_AnyInput_USumEqualsProductOfSizes()
    {
        var rng = new Random(20260613); // fixed seed: deterministic
        Assert.Multiple(() =>
        {
            for (int trial = 0; trial < 50; trial++)
            {
                int n1 = rng.Next(1, 12), n2 = rng.Next(1, 12);
                var g1 = Enumerable.Range(0, n1).Select(_ => Math.Round(rng.NextDouble() * 10, 1)).ToArray();
                var g2 = Enumerable.Range(0, n2).Select(_ => Math.Round(rng.NextDouble() * 10, 1)).ToArray();

                var r = MetagenomicsAnalyzer.MannWhitneyU(g1, g2);

                Assert.That(r.U1 + r.U2, Is.EqualTo((double)n1 * n2).Within(1e-9),
                    $"INV-01: U1+U2 must equal n1*n2 = {n1 * n2} (trial {trial}).");
            }
        });
    }

    // M7 — Ties: midrank + tie-corrected σ. Group1=[1,1,2], Group2=[1,2,2]; derive σ by hand.
    [Test]
    public void MannWhitneyU_WithTies_UsesTieCorrectedSigma()
    {
        // Pooled values: 1,1,1 (rank 1-3 → midrank 2), 2,2,2 (rank 4-6 → midrank 5). n1=n2=3, n=6.
        // R1 (group1 = [1,1,2]): two 1s (rank 2 each) + one 2 (rank 5) = 9.
        // U1 = 9 - 3*4/2 = 9 - 6 = 3 ; U2 = 9 - 3 = 6.
        // Tie groups: t=3 (value 1), t=3 (value 2) → Σ(t³-t) = 24 + 24 = 48.
        // variance = n1*n2/12 * ((n+1) - Σ/(n(n-1))) = 9/12 * (7 - 48/30) = 0.75*(7-1.6)=0.75*5.4=4.05.
        var g1 = new double[] { 1, 1, 2 };
        var g2 = new double[] { 1, 2, 2 };
        var r = MetagenomicsAnalyzer.MannWhitneyU(g1, g2, useContinuityCorrection: false);

        const double expectedSigma = 2.0124611797498108; // sqrt(4.05)
        const double expectedZ = 1.5 / expectedSigma;      // (max(3,6)-4.5)=1.5
        Assert.Multiple(() =>
        {
            Assert.That(r.U1, Is.EqualTo(3.0).Within(1e-10), "U1 = R1 - n1(n1+1)/2 = 9 - 6 = 3 with midranks.");
            Assert.That(r.U2, Is.EqualTo(6.0).Within(1e-10), "U2 = n1*n2 - U1 = 9 - 3 = 6.");
            Assert.That(r.Z, Is.EqualTo(expectedZ).Within(1e-10),
                "Tie-corrected variance 4.05 → σ=sqrt(4.05); z=(|6-4.5|)/σ (Mann & Whitney tie correction).");
        });
    }

    // M8 / INV-05 — all-tied groups → σ=0 → p=1, z=0.
    [Test]
    public void MannWhitneyU_IdenticalGroups_PValueIsOne()
    {
        var g1 = new double[] { 5, 5, 5 };
        var g2 = new double[] { 5, 5, 5 };
        var r = MetagenomicsAnalyzer.MannWhitneyU(g1, g2);

        Assert.Multiple(() =>
        {
            Assert.That(r.Z, Is.EqualTo(0.0).Within(1e-12), "All observations tied: σ→0, z defined as 0.");
            Assert.That(r.PValue, Is.EqualTo(1.0).Within(1e-12), "INV-05: no evidence against H0 → p = 1.");
        });
    }

    // S1 / INV-06 — symmetry: swapping groups leaves p unchanged.
    [Test]
    public void MannWhitneyU_SwapGroups_PValueUnchanged()
    {
        var forward = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY);
        var swapped = MetagenomicsAnalyzer.MannWhitneyU(ScipyY, ScipyX);

        Assert.That(swapped.PValue, Is.EqualTo(forward.PValue).Within(1e-12),
            "INV-06: the two-tailed p-value is symmetric under group swap.");
    }

    // S2 / INV-03 — p-value stays within [0,1].
    [Test]
    public void MannWhitneyU_AnyInput_PValueInUnitInterval()
    {
        var r = MetagenomicsAnalyzer.MannWhitneyU(
            new double[] { 100, 200, 300 }, new double[] { 1, 2, 3 });

        Assert.That(r.PValue, Is.InRange(0.0, 1.0),
            "INV-03: p-value must lie in [0,1] even for fully separated groups.");
    }

    // C1 — continuity correction increases the (two-tailed) p-value (SciPy parameter semantics).
    [Test]
    public void MannWhitneyU_ContinuityCorrection_IncreasesPValue()
    {
        var withCc = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY, useContinuityCorrection: true);
        var noCc = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY, useContinuityCorrection: false);

        Assert.That(withCc.PValue, Is.GreaterThan(noCc.PValue),
            "Continuity correction reduces |U - m_U|, lowering z and thus raising the p-value (0.1113 > 0.0864).");
    }

    [Test]
    public void MannWhitneyU_NullGroup_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => MetagenomicsAnalyzer.MannWhitneyU(null!, ScipyY),
                "Null group1 must throw ArgumentNullException.");
            Assert.Throws<ArgumentNullException>(() => MetagenomicsAnalyzer.MannWhitneyU(ScipyX, null!),
                "Null group2 must throw ArgumentNullException.");
        });
    }

    [Test]
    public void MannWhitneyU_EmptyGroup_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => MetagenomicsAnalyzer.MannWhitneyU(Array.Empty<double>(), ScipyY),
            "Empty group must throw ArgumentException (test undefined for n=0).");
    }

    #endregion

    #region FindSignificantTaxa

    private static IReadOnlyDictionary<string, double> Profile(params (string Taxon, double Abundance)[] xs)
        => xs.ToDictionary(x => x.Taxon, x => x.Abundance);

    // M9 — separated taxon flagged significant, overlapping taxon not (per-taxon rank-sum, Xia & Sun 2017).
    [Test]
    public void FindSignificantTaxa_SeparatedAndOverlapping_FlagsCorrectly()
    {
        // 6 samples per group. TaxonA: group1 low, group2 high (fully separated → smallest p).
        // TaxonB: heavily overlapping → not significant at 0.05.
        var profiles = new List<IReadOnlyDictionary<string, double>>
        {
            Profile(("TaxonA", 1), ("TaxonB", 5)),
            Profile(("TaxonA", 2), ("TaxonB", 6)),
            Profile(("TaxonA", 3), ("TaxonB", 5)),
            Profile(("TaxonA", 10), ("TaxonB", 6)),
            Profile(("TaxonA", 11), ("TaxonB", 5)),
            Profile(("TaxonA", 12), ("TaxonB", 6)),
        };
        var groups = new[] { 1, 1, 1, 2, 2, 2 };

        var result = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, pThreshold: 0.05);

        var a = result.Single(t => t.Taxon == "TaxonA");
        var b = result.Single(t => t.Taxon == "TaxonB");
        Assert.Multiple(() =>
        {
            // Sourced p-values (scipy 1.13.1, mannwhitneyu method='asymptotic', use_continuity=True,
            // alternative='two-sided'), recomputed this validation session:
            //   TaxonA: group1=[1,2,3] vs group2=[10,11,12] (fully separated, U1=0,U2=9) → p=0.08085559837005224
            //   TaxonB: group1=[5,6,5] vs group2=[6,5,6] (overlapping, U=3)              → p=0.6192567541768621
            // Erf approximation tolerance applies (A&S 7.1.26, |ε| ≤ 1.5e-7 ⇒ p within ~1e-6).
            Assert.That(a.PValue, Is.EqualTo(0.08085559837005224).Within(ErfPTolerance),
                "scipy mannwhitneyu([1,2,3],[10,11,12]) two-sided asymptotic p (cc) = 0.08085559837005224.");
            Assert.That(b.PValue, Is.EqualTo(0.6192567541768621).Within(ErfPTolerance),
                "scipy mannwhitneyu([5,6,5],[6,5,6]) two-sided asymptotic p (cc) = 0.6192567541768621.");
            // At p=0.05 neither taxon clears the threshold for n1=n2=3 (smallest attainable cc p is 0.0809).
            Assert.That(a.Significant, Is.False,
                "INV-04: TaxonA p=0.0809 ≥ 0.05 → not significant (3-vs-3 cannot reach 0.05).");
            Assert.That(b.Significant, Is.False,
                "INV-04: TaxonB p=0.619 ≥ 0.05 → not significant.");
            Assert.That(a.PValue, Is.LessThan(b.PValue),
                "Fully separated TaxonA must have a smaller p-value than overlapping TaxonB.");
            Assert.That(result.Select(t => t.PValue), Is.Ordered,
                "Results must be ordered by ascending p-value.");
        });
    }

    // M9b — with enough replicates (4 vs 4) a fully separated taxon clears p<0.05 and is flagged.
    [Test]
    public void FindSignificantTaxa_FullySeparated_FlagsSignificant()
    {
        // group1=[1,2,3,4] vs group2=[10,11,12,13]; scipy mannwhitneyu two-sided asymptotic (cc)
        // p = 0.030382821976577504 < 0.05 (recomputed this session, scipy 1.13.1;
        // hand-check: U1=0,U2=16, μ=8, σ=sqrt(4·4·9/12)=sqrt(12), z=(8-0.5)/σ=2.16506, p=2·SF(z)).
        var profiles = new List<IReadOnlyDictionary<string, double>>
        {
            Profile(("TaxonA", 1)), Profile(("TaxonA", 2)), Profile(("TaxonA", 3)), Profile(("TaxonA", 4)),
            Profile(("TaxonA", 10)), Profile(("TaxonA", 11)), Profile(("TaxonA", 12)), Profile(("TaxonA", 13)),
        };
        var groups = new[] { 1, 1, 1, 1, 2, 2, 2, 2 };

        var a = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, pThreshold: 0.05)
            .Single(t => t.Taxon == "TaxonA");

        Assert.Multiple(() =>
        {
            Assert.That(a.PValue, Is.EqualTo(0.030382821976577504).Within(ErfPTolerance),
                "scipy mannwhitneyu([1,2,3,4],[10,11,12,13]) two-sided asymptotic p (cc) = 0.030382821976577504.");
            Assert.That(a.Significant, Is.True,
                "INV-04: p=0.0285 < 0.05 → significant.");
        });
    }

    // S3 — taxon absent in some profiles is treated as abundance 0 (ASM-03).
    [Test]
    public void FindSignificantTaxa_AbsentTaxon_TreatedAsZeroAbundance()
    {
        // TaxonZ present (high) only in group2; absent (→0) in group1.
        var profiles = new List<IReadOnlyDictionary<string, double>>
        {
            Profile(("TaxonZ", 0)),       // explicit 0
            Profile(("Other", 1)),        // TaxonZ absent → treated as 0
            Profile(("TaxonZ", 50)),
            Profile(("TaxonZ", 60)),
        };
        var groups = new[] { 1, 1, 2, 2 };

        var result = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups);

        var z = result.Single(t => t.Taxon == "TaxonZ");
        // group1 = [0, 0] (one explicit 0, one filled because TaxonZ is absent from that profile),
        // group2 = [50, 60]. The absent-→-0 fill must make this a fully separated 2-vs-2 comparison.
        // Sourced value (scipy 1.13.1, mannwhitneyu([0,0],[50,60]) method='asymptotic',
        // use_continuity=True, alternative='two-sided'), recomputed this session: p = 0.22067136191984682.
        // This exact value can only be produced if the absent TaxonZ was filled with 0 (ASM-03);
        // a NaN/skip or a different fill would not reproduce it, so it genuinely locks the behaviour.
        Assert.That(z.PValue, Is.EqualTo(0.22067136191984682).Within(ErfPTolerance),
            "Absent TaxonZ filled with 0 → mannwhitneyu([0,0],[50,60]) two-sided asymptotic p (cc) = 0.22067136191984682.");
    }

    [Test]
    public void FindSignificantTaxa_NullInputs_Throw()
    {
        var profiles = new List<IReadOnlyDictionary<string, double>> { Profile(("T", 1)) };
        var groups = new[] { 1 };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => MetagenomicsAnalyzer.FindSignificantTaxa(null!, groups),
                "Null profiles must throw.");
            Assert.Throws<ArgumentNullException>(() => MetagenomicsAnalyzer.FindSignificantTaxa(profiles, null!),
                "Null groups must throw.");
        });
    }

    [Test]
    public void FindSignificantTaxa_EmptyProfiles_ReturnsEmpty()
    {
        var result = MetagenomicsAnalyzer.FindSignificantTaxa(
            new List<IReadOnlyDictionary<string, double>>(), Array.Empty<int>());

        Assert.That(result, Is.Empty, "Empty profiles yield an empty result list.");
    }

    [Test]
    public void FindSignificantTaxa_MismatchedLengths_Throws()
    {
        var profiles = new List<IReadOnlyDictionary<string, double>> { Profile(("T", 1)) };
        Assert.Throws<ArgumentException>(
            () => MetagenomicsAnalyzer.FindSignificantTaxa(profiles, new[] { 1, 2 }),
            "profiles and groups length mismatch must throw ArgumentException.");
    }

    [Test]
    public void FindSignificantTaxa_SingleGroupOnly_Throws()
    {
        var profiles = new List<IReadOnlyDictionary<string, double>>
        {
            Profile(("T", 1)), Profile(("T", 2)),
        };
        Assert.Throws<ArgumentException>(
            () => MetagenomicsAnalyzer.FindSignificantTaxa(profiles, new[] { 1, 1 }),
            "Only group 1 present (no group 2) must throw ArgumentException.");
    }

    [Test]
    public void FindSignificantTaxa_InvalidGroupLabel_Throws()
    {
        var profiles = new List<IReadOnlyDictionary<string, double>>
        {
            Profile(("T", 1)), Profile(("T", 2)),
        };
        Assert.Throws<ArgumentException>(
            () => MetagenomicsAnalyzer.FindSignificantTaxa(profiles, new[] { 1, 3 }),
            "Group label other than 1/2 must throw ArgumentException.");
    }

    #endregion
}
