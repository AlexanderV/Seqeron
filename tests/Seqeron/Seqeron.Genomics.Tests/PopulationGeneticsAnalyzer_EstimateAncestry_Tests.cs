// POP-ANCESTRY-001 — Ancestry Estimation (supervised/projection ADMIXTURE, EM with fixed F)
// Evidence: docs/Evidence/POP-ANCESTRY-001-Evidence.md
// TestSpec: tests/TestSpecs/POP-ANCESTRY-001.md
// Source: Alexander DH, Novembre J, Lange K (2009). Genome Research 19(9):1655-1664.
//         Equations 2 (log-likelihood), 4 (FRAPPE EM ancestry update), 5 (convergence).

using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.Population.PopulationGeneticsAnalyzer;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class PopulationGeneticsAnalyzer_EstimateAncestry_Tests
{
    private const double Tol = 1e-10;

    private static List<(string, IReadOnlyList<int>)> Inds(params (string, int[])[] xs) =>
        xs.Select(x => (x.Item1, (IReadOnlyList<int>)x.Item2)).ToList();

    private static List<(string, IReadOnlyList<double>)> Refs(params (string, double[])[] xs) =>
        xs.Select(x => (x.Item1, (IReadOnlyList<double>)x.Item2)).ToList();

    #region EstimateAncestry

    // M1 — Eq. 4 one iteration from uniform start, symmetric 2-pop panel.
    // f_A=(0.8,0.2), f_B=(0.2,0.8), g=(2,0): each SNP mix=0.5, A gets 2*(0.5*0.8/0.5)=1.6, B 0.4;
    // q_A=(1.6+1.6)/(2*2)=0.8, q_B=0.2.
    [Test]
    public void EstimateAncestry_OneIterationSymmetricPanel_MatchesEq4()
    {
        var inds = Inds(("IND1", new[] { 2, 0 }));
        var refs = Refs(("A", new[] { 0.8, 0.2 }), ("B", new[] { 0.2, 0.8 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 1).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.Proportions["A"], Is.EqualTo(0.8).Within(Tol),
                "Eq. 4 one iteration on f_A=(0.8,0.2)/f_B=(0.2,0.8), g=(2,0) gives q_A = 0.8.");
            Assert.That(r.Proportions["B"], Is.EqualTo(0.2).Within(Tol),
                "Complement of q_A; Eq. 4 yields q_B = 0.2.");
        });
    }

    // M2 — Eq. 4, single SNP g=2, f=(0.9,0.1): q_A = (1/2)*[2*(0.5*0.9/0.5)] = 0.9.
    [Test]
    public void EstimateAncestry_SingleSnpHomozygousAllele1_MatchesEq4()
    {
        var inds = Inds(("IND1", new[] { 2 }));
        var refs = Refs(("A", new[] { 0.9 }), ("B", new[] { 0.1 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 1).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.Proportions["A"], Is.EqualTo(0.9).Within(Tol),
                "g=2, f=(0.9,0.1): one Eq. 4 iteration gives q_A = 0.9.");
            Assert.That(r.Proportions["B"], Is.EqualTo(0.1).Within(Tol),
                "q_B = 0.1 by complement.");
        });
    }

    // M3 — Eq. 4, single SNP g=0, f=(0.9,0.1): mirror of M2 via (1-f). q_A=0.1, q_B=0.9.
    [Test]
    public void EstimateAncestry_SingleSnpHomozygousAllele2_MatchesEq4()
    {
        var inds = Inds(("IND1", new[] { 0 }));
        var refs = Refs(("A", new[] { 0.9 }), ("B", new[] { 0.1 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 1).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.Proportions["A"], Is.EqualTo(0.1).Within(Tol),
                "g=0 uses (1-f); A has 1-f=0.1, so q_A = 0.1 after one Eq. 4 iteration.");
            Assert.That(r.Proportions["B"], Is.EqualTo(0.9).Within(Tol),
                "B has 1-f=0.9, so q_B = 0.9.");
        });
    }

    // M4 — INV-01: proportions sum to 1 for an arbitrary informative input.
    [Test]
    public void EstimateAncestry_ArbitraryInput_ProportionsSumToOne()
    {
        var inds = Inds(("IND1", new[] { 2, 1, 0, 1 }));
        var refs = Refs(
            ("A", new[] { 0.7, 0.3, 0.6, 0.2 }),
            ("B", new[] { 0.2, 0.8, 0.1, 0.9 }),
            ("C", new[] { 0.5, 0.4, 0.5, 0.5 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 100).Single();

        Assert.That(r.Proportions.Values.Sum(), Is.EqualTo(1.0).Within(Tol),
            "Eq. 4 normalizes by 2J each step, so Sigma_k q_ik = 1 (constraint, Alexander 2009).");
    }

    // M5 — convergence: diagnostic individual driven to its source population (q_A -> 1).
    [Test]
    public void EstimateAncestry_DiagnosticIndividual_ConvergesToSource()
    {
        var inds = Inds(("IND1", new[] { 2, 0 }));
        var refs = Refs(("A", new[] { 0.8, 0.2 }), ("B", new[] { 0.2, 0.8 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 1000).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.Proportions["A"], Is.EqualTo(1.0).Within(1e-3),
                "EM monotonically ascends to the MLE; this genotype is fully population-A ancestry.");
            Assert.That(r.Proportions["B"], Is.EqualTo(0.0).Within(1e-3),
                "Population-B ancestry vanishes at the maximum likelihood estimate.");
        });
    }

    // M6 — INV-04: identical reference panels are uninformative; uniform q stays uniform.
    [Test]
    public void EstimateAncestry_IdenticalPanels_StaysUniform()
    {
        var inds = Inds(("IND1", new[] { 2, 1, 0 }));
        var refs = Refs(("A", new[] { 0.5, 0.5, 0.5 }), ("B", new[] { 0.5, 0.5, 0.5 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 100).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.Proportions["A"], Is.EqualTo(0.5).Within(Tol),
                "Identical panels make a/b proportional; uniform q is a fixed point of Eq. 4.");
            Assert.That(r.Proportions["B"], Is.EqualTo(0.5).Within(Tol),
                "Same reasoning for population B.");
        });
    }

    // M7 — symmetric heterozygote: f=(0.9,0.1), g=1 keeps a uniform individual uniform.
    [Test]
    public void EstimateAncestry_SymmetricHeterozygote_StaysUniform()
    {
        var inds = Inds(("IND1", new[] { 1 }));
        var refs = Refs(("A", new[] { 0.9 }), ("B", new[] { 0.1 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 1).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.Proportions["A"], Is.EqualTo(0.5).Within(Tol),
                "g=1: a-term gives 0.9 to A, b-term gives 0.1; (0.9+0.1)/2 = 0.5 (Eq. 4).");
            Assert.That(r.Proportions["B"], Is.EqualTo(0.5).Within(Tol),
                "By symmetry q_B = 0.5.");
        });
    }

    // S1 — INV-03: Eq. 2 log-likelihood is non-decreasing across EM iterations (property test).
    [Test]
    public void EstimateAncestry_LogLikelihood_IsNonDecreasingAcrossIterations()
    {
        var inds = Inds(("IND1", new[] { 2, 0, 1, 2, 0 }));
        var refs = Refs(
            ("A", new[] { 0.8, 0.2, 0.5, 0.7, 0.1 }),
            ("B", new[] { 0.1, 0.9, 0.4, 0.2, 0.8 }),
            ("C", new[] { 0.4, 0.4, 0.6, 0.5, 0.5 }));

        double[] LogLik(int iters)
        {
            var p = EstimateAncestry(inds, refs, maxIterations: iters).Single().Proportions;
            double[] q = { p["A"], p["B"], p["C"] };
            double[] fA = { 0.8, 0.2, 0.5, 0.7, 0.1 };
            double[] fB = { 0.1, 0.9, 0.4, 0.2, 0.8 };
            double[] fC = { 0.4, 0.4, 0.6, 0.5, 0.5 };
            int[] g = { 2, 0, 1, 2, 0 };
            double l = 0;
            for (int j = 0; j < g.Length; j++)
            {
                double mix1 = q[0] * fA[j] + q[1] * fB[j] + q[2] * fC[j];
                double mix0 = q[0] * (1 - fA[j]) + q[1] * (1 - fB[j]) + q[2] * (1 - fC[j]);
                l += g[j] * System.Math.Log(mix1) + (2 - g[j]) * System.Math.Log(mix0);
            }
            return new[] { l };
        }

        double l1 = LogLik(1)[0];
        double l2 = LogLik(2)[0];
        double l3 = LogLik(3)[0];

        Assert.Multiple(() =>
        {
            Assert.That(l2, Is.GreaterThanOrEqualTo(l1 - Tol),
                "EM is a monotone ascent algorithm: L after 2 iterations >= L after 1 (Eq. 2/5).");
            Assert.That(l3, Is.GreaterThanOrEqualTo(l2 - Tol),
                "EM ascent continues: L after 3 iterations >= L after 2.");
        });
    }

    // S2 — empty individuals -> empty result.
    [Test]
    public void EstimateAncestry_NoIndividuals_ReturnsEmpty()
    {
        var refs = Refs(("A", new[] { 0.8 }), ("B", new[] { 0.2 }));
        Assert.That(EstimateAncestry(Inds(), refs).ToList(), Is.Empty,
            "No individuals to estimate yields an empty result.");
    }

    // S3 — empty reference panels -> empty result.
    [Test]
    public void EstimateAncestry_NoReferencePanels_ReturnsEmpty()
    {
        var inds = Inds(("IND1", new[] { 2, 0 }));
        Assert.That(EstimateAncestry(inds, Refs()).ToList(), Is.Empty,
            "Without reference panels there is nothing to estimate against.");
    }

    // S4 — individual whose genotype length != panel SNP count is skipped.
    [Test]
    public void EstimateAncestry_MismatchedGenotypeLength_SkipsIndividual()
    {
        var inds = Inds(("SHORT", new[] { 2 }), ("OK", new[] { 2, 0 }));
        var refs = Refs(("A", new[] { 0.8, 0.2 }), ("B", new[] { 0.2, 0.8 }));

        var results = EstimateAncestry(inds, refs, maxIterations: 1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1),
                "The genotype of wrong length cannot be index-aligned to the panels; it is skipped.");
            Assert.That(results[0].IndividualId, Is.EqualTo("OK"),
                "Only the correctly sized individual is returned.");
        });
    }

    // S6 — missing genotype (value outside {0,1,2}) is skipped; result equals the
    // estimate from the informative SNP alone. Genotype (2, -1): SNP1 informative,
    // SNP0 missing => identical to the single-SNP g=2 case (q_A=0.9 after 1 iteration).
    [Test]
    public void EstimateAncestry_MissingGenotype_SkipsThatSnp()
    {
        var indsMissing = Inds(("IND1", new[] { -1, 2 }));
        var refs = Refs(("A", new[] { 0.3, 0.9 }), ("B", new[] { 0.7, 0.1 }));

        var rMissing = EstimateAncestry(indsMissing, refs, maxIterations: 1).Single();

        Assert.Multiple(() =>
        {
            Assert.That(rMissing.Proportions["A"], Is.EqualTo(0.9).Within(Tol),
                "Missing SNP0 contributes no Eq. 2 term; only informative SNP1 (g=2,f=0.9) drives q_A=0.9.");
            Assert.That(rMissing.Proportions["B"], Is.EqualTo(0.1).Within(Tol),
                "q_B = 0.1 from the single informative SNP.");
        });
    }

    // S7 — all genotypes missing -> uniform prior returned (no informative term).
    [Test]
    public void EstimateAncestry_AllGenotypesMissing_ReturnsUniformPrior()
    {
        var inds = Inds(("IND1", new[] { -1, 3 }));
        var refs = Refs(("A", new[] { 0.8, 0.2 }), ("B", new[] { 0.2, 0.8 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 10).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.Proportions["A"], Is.EqualTo(0.5).Within(Tol),
                "With no informative SNP, EM cannot update; the uniform prior 1/K is returned.");
            Assert.That(r.Proportions["B"], Is.EqualTo(0.5).Within(Tol),
                "Uniform prior for K=2.");
        });
    }

    // S5 — INV-04 / result keys: three panels -> exactly those keys, summing to 1.
    [Test]
    public void EstimateAncestry_ThreePanels_KeysMatchAndSumToOne()
    {
        var inds = Inds(("IND1", new[] { 2, 1, 0 }));
        var refs = Refs(
            ("A", new[] { 0.9, 0.5, 0.1 }),
            ("B", new[] { 0.1, 0.5, 0.9 }),
            ("C", new[] { 0.5, 0.5, 0.5 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 50).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.Proportions.Keys, Is.EquivalentTo(new[] { "A", "B", "C" }),
                "Result is keyed by exactly the supplied reference-population ids.");
            Assert.That(r.Proportions.Values.Sum(), Is.EqualTo(1.0).Within(Tol),
                "Constraint Sigma_k q_ik = 1 holds for K = 3.");
        });
    }

    // INV-2 — every proportion lies in [0,1] (constraint q_ik >= 0 with Sigma_k q_ik = 1,
    // Alexander 2009). Exercised on an informative 3-population case.
    [Test]
    public void EstimateAncestry_ArbitraryInput_ProportionsInUnitInterval()
    {
        var inds = Inds(("IND1", new[] { 2, 1, 0, 2, 1 }));
        var refs = Refs(
            ("A", new[] { 0.8, 0.3, 0.6, 0.9, 0.2 }),
            ("B", new[] { 0.1, 0.9, 0.2, 0.1, 0.7 }),
            ("C", new[] { 0.5, 0.4, 0.5, 0.5, 0.5 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 100).Single();

        Assert.Multiple(() =>
        {
            foreach (var v in r.Proportions.Values)
            {
                Assert.That(v, Is.GreaterThanOrEqualTo(0.0 - Tol),
                    "Constraint q_ik >= 0 (Alexander 2009).");
                Assert.That(v, Is.LessThanOrEqualTo(1.0 + Tol),
                    "q_ik <= 1 because q_ik >= 0 and Sigma_k q_ik = 1 (Alexander 2009).");
            }
        });
    }

    // maxIterations=0 boundary — no EM update runs, so the uniform prior 1/K is returned
    // (same code path as the all-missing case; documented contract maxIterations >= 0).
    [Test]
    public void EstimateAncestry_ZeroIterations_ReturnsUniformPrior()
    {
        var inds = Inds(("IND1", new[] { 2, 0 }));
        var refs = Refs(("A", new[] { 0.8, 0.2 }), ("B", new[] { 0.2, 0.8 }));

        var r = EstimateAncestry(inds, refs, maxIterations: 0).Single();

        Assert.Multiple(() =>
        {
            Assert.That(r.Proportions["A"], Is.EqualTo(0.5).Within(Tol),
                "With zero iterations the EM never updates; the uniform prior 1/K is returned.");
            Assert.That(r.Proportions["B"], Is.EqualTo(0.5).Within(Tol),
                "Uniform prior for K=2.");
        });
    }

    // C1 — label invariance (Eq. 2): permuting panels permutes the proportions consistently.
    [Test]
    public void EstimateAncestry_PermutedPanels_PermutesProportions()
    {
        var inds = Inds(("IND1", new[] { 2, 0, 1 }));
        var refsAB = Refs(("A", new[] { 0.8, 0.2, 0.6 }), ("B", new[] { 0.2, 0.8, 0.4 }));
        var refsBA = Refs(("B", new[] { 0.2, 0.8, 0.4 }), ("A", new[] { 0.8, 0.2, 0.6 }));

        var rAB = EstimateAncestry(inds, refsAB, maxIterations: 50).Single();
        var rBA = EstimateAncestry(inds, refsBA, maxIterations: 50).Single();

        Assert.Multiple(() =>
        {
            Assert.That(rBA.Proportions["A"], Is.EqualTo(rAB.Proportions["A"]).Within(Tol),
                "Eq. 2 is invariant under panel relabeling; q_A is unchanged by panel order.");
            Assert.That(rBA.Proportions["B"], Is.EqualTo(rAB.Proportions["B"]).Within(Tol),
                "q_B is likewise unchanged by panel order.");
        });
    }

    #endregion
}
