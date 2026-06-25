// MIRNA-TARGET-001 (McCaskill partition-function folder) — Turner-2004 McCaskill partition
//   function, base-pair / unpaired probabilities, ensemble free energy, and region accessibility (SA).
// Evidence: docs/Evidence/MIRNA-TARGET-001-Evidence.md  (§ SA — structural accessibility)
// TestSpec: tests/TestSpecs/MIRNA-TARGET-001.md          (§ MCC-001..MCC-0xx)
// Source: McCaskill JS (1990). Biopolymers 29(6-7):1105-1119 (PMID 1695107).
//         Lorenz et al. (2011). ViennaRNA Package 2.0, Algorithms Mol Biol 6:26.
//         RNAplfold man page + Bernhart et al. (2006) Bioinformatics 22:614.
//         Turner-2004 nearest-neighbour parameters (NNDB, rna.urmc.rochester.edu/NNDB/turner04).

using System;
using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class RnaSecondaryStructure_UnpairedProbabilities_Tests
{
    // RT = R·T/1000 with R = 1.987 cal/(mol·K), T = 310.15 K → 0.61626805 kcal/mol.
    private const double Rt = 1.987 * 310.15 / 1000.0;

    #region Analytic tiny case (the solid pin)

    // MCC-001 — GAAAC can adopt only the open chain OR the single hairpin closed by G(0)-C(4)
    // over the 3-nt loop "AAA". The Turner hairpin energy is CalculateHairpinLoopEnergy("AAA",
    // 'G','C',false) = 5.4 kcal/mol (G-C closing ⇒ no terminal-AU penalty). Hence, analytically,
    //   w = exp(-5.4/RT),  Z = 1 + w,  P(0,4) = w/Z,  p_unpaired = 1 - P(0,4).
    [Test]
    public void CalculateUnpairedProbabilities_GAAAC_MatchesAnalyticPartitionFunction()
    {
        // Hand-derived from the Turner hairpin energy (verified: CalculateHairpinLoopEnergy = 5.4).
        const double dG = 5.4;
        double w = Math.Exp(-dG / Rt);
        double expectedZ = 1.0 + w;                 // 1.0001565052764922
        double expectedP = w / expectedZ;           // 0.00015648078642340854
        double expectedPu = 1.0 - expectedP;        // 0.9998435192135765

        var r = CalculateUnpairedProbabilities("GAAAC");

        Assert.Multiple(() =>
        {
            Assert.That(r.PartitionFunction, Is.EqualTo(expectedZ).Within(1e-12),
                "Z = 1 + exp(-ΔG_hairpin/RT) over {open chain, single G-C hairpin}");
            Assert.That(r.BasePairProbabilities[(0, 4)], Is.EqualTo(expectedP).Within(1e-12),
                "P(0,4) = exp(-ΔG/RT)/Z");
            Assert.That(r.UnpairedProbabilities[0], Is.EqualTo(expectedPu).Within(1e-12),
                "p_unpaired(0) = 1 - P(0,4)");
            Assert.That(r.UnpairedProbabilities[4], Is.EqualTo(expectedPu).Within(1e-12),
                "p_unpaired(4) = 1 - P(0,4) (the only pair (0,4) involves position 4)");
            Assert.That(r.UnpairedProbabilities[2], Is.EqualTo(1.0).Within(1e-12),
                "interior loop base 2 can never pair ⇒ p_unpaired = 1");
            Assert.That(r.EnsembleFreeEnergy, Is.EqualTo(-Rt * Math.Log(expectedZ)).Within(1e-12),
                "ensemble free energy = -RT·ln Z");
        });
    }

    // MCC-002 — exact analytic value, copied verbatim (would fail any wrong/MFE-as-PF impl).
    [Test]
    public void CalculateUnpairedProbabilities_GAAAC_ExactNumericValues()
    {
        var r = CalculateUnpairedProbabilities("GAAAC");

        Assert.Multiple(() =>
        {
            Assert.That(r.PartitionFunction, Is.EqualTo(1.0001565052764922).Within(1e-15),
                "Z exactly 1.0001565052764922");
            Assert.That(r.BasePairProbabilities[(0, 4)], Is.EqualTo(0.00015648078642340854).Within(1e-15),
                "P(0,4) exactly 0.00015648078642340854");
            Assert.That(r.UnpairedProbabilities[0], Is.EqualTo(0.9998435192135765).Within(1e-15),
                "p_unpaired(0) exactly 0.9998435192135765");
        });
    }

    // MCC-002b — SECOND independent analytic pin (4-nt loop, exercises the terminal-mismatch path
    // that GAAAC's 3-nt loop does NOT). CAAAAG can adopt only the open chain OR the single hairpin
    // closed by C(0)-G(5) over the 4-nt loop "AAAA". The Turner hairpin energy for a size-4 loop is
    //   ΔG = initiation(4) + terminalMismatch[c5=C, mm5=A, mm3=A, c3=G]
    //      = 5.6 + (-1.5)  = 4.1 kcal/mol      (NNDB tm key "CAAG" = -1.5; C-G closing ⇒ no AU end).
    // The closing helix spans the whole 6-nt sequence ⇒ no dangling ends, no terminal-AU penalty, so
    //   w = exp(-4.1/RT),  Z = 1 + w,  P(0,5) = w/Z,  p_unpaired = 1 - P(0,5).
    // Hand-derived from the published Turner-2004 tables (NOT read back from the code):
    //   Z = 1.0012902114608,  P(0,5) = 0.0012885489601637966,  p_unpaired(0) = 0.9987114510398362.
    [Test]
    public void CalculateUnpairedProbabilities_CAAAAG_MatchesAnalyticPartitionFunction()
    {
        const double dG = 5.6 - 1.5;               // initiation(4) + terminal mismatch CAAG
        double w = Math.Exp(-dG / Rt);
        double expectedZ = 1.0 + w;
        double expectedP = w / expectedZ;
        double expectedPu = 1.0 - expectedP;

        var r = CalculateUnpairedProbabilities("CAAAAG");

        Assert.Multiple(() =>
        {
            Assert.That(r.PartitionFunction, Is.EqualTo(expectedZ).Within(1e-12),
                "Z = 1 + exp(-(5.6-1.5)/RT) over {open chain, single C-G hairpin with 4-nt loop}");
            Assert.That(r.PartitionFunction, Is.EqualTo(1.0012902114608).Within(1e-12),
                "Z exactly 1.0012902114608 (analytic)");
            Assert.That(r.BasePairProbabilities[(0, 5)], Is.EqualTo(expectedP).Within(1e-12),
                "P(0,5) = exp(-ΔG/RT)/Z");
            Assert.That(r.BasePairProbabilities[(0, 5)], Is.EqualTo(0.0012885489601637966).Within(1e-12),
                "P(0,5) exactly 0.0012885489601637966 (analytic)");
            Assert.That(r.UnpairedProbabilities[0], Is.EqualTo(expectedPu).Within(1e-12),
                "p_unpaired(0) = 1 - P(0,5)");
            Assert.That(r.UnpairedProbabilities[0], Is.EqualTo(0.9987114510398362).Within(1e-12),
                "p_unpaired(0) exactly 0.9987114510398362 (analytic)");
            Assert.That(r.UnpairedProbabilities[2], Is.EqualTo(1.0).Within(1e-12),
                "interior loop base 2 can never pair ⇒ p_unpaired = 1");
            Assert.That(r.EnsembleFreeEnergy, Is.EqualTo(-Rt * Math.Log(expectedZ)).Within(1e-12),
                "ensemble free energy = -RT·ln Z");
        });
    }

    #endregion

    #region Invariants

    // MCC-003 — every P(i,j) ∈ [0,1]; Σ_j P(i,j) ≤ 1; p_unpaired(i) = 1 - Σ_j P(i,j) ∈ [0,1];
    // Z > 0. Checked across hairpin, multiloop, GC-rich, and mixed sequences.
    [TestCase("GGGAAACCC")]
    [TestCase("GGGAAACCCAAAGGGAAACCC")]
    [TestCase("GCGCGCAAAAGCGCGC")]
    [TestCase("AUGCAUGCAUGCAUGCAUGC")]
    [TestCase("GGGGAAAACCCCUUUUGGGGAAAACCCC")]
    public void CalculateUnpairedProbabilities_Invariants_Hold(string seq)
    {
        var r = CalculateUnpairedProbabilities(seq);
        var sumPair = new double[seq.Length];
        foreach (var kv in r.BasePairProbabilities)
        {
            sumPair[kv.Key.I] += kv.Value;
            sumPair[kv.Key.J] += kv.Value;
        }

        Assert.Multiple(() =>
        {
            Assert.That(r.PartitionFunction, Is.GreaterThan(0.0), "Z > 0");
            foreach (var kv in r.BasePairProbabilities)
                Assert.That(kv.Value, Is.InRange(0.0, 1.0 + 1e-12),
                    $"P{kv.Key} must lie in [0,1]");
            for (int i = 0; i < seq.Length; i++)
            {
                Assert.That(sumPair[i], Is.LessThanOrEqualTo(1.0 + 1e-9),
                    $"Σ_j P({i},j) ≤ 1");
                Assert.That(r.UnpairedProbabilities[i], Is.InRange(-1e-9, 1.0 + 1e-9),
                    $"p_unpaired({i}) ∈ [0,1]");
                // McCaskill consistency: 1 - Σ_j P(i,j) == p_unpaired(i).
                Assert.That(r.UnpairedProbabilities[i], Is.EqualTo(1.0 - sumPair[i]).Within(1e-6),
                    $"p_unpaired({i}) = 1 - Σ_j P({i},j)");
            }
        });
    }

    // MCC-004 — ensemble free energy −RT·ln Z is a LOWER bound on the MFE (the MFE structure is
    // a single term in Z, so Z ≥ exp(-MFE/RT) ⇒ -RT·ln Z ≤ MFE). Lorenz et al. 2011 / statistical
    // mechanics. (Pseudoknot-free MFE folder uses the same Turner-2004 model.)
    [TestCase("GGGAAACCC")]
    [TestCase("GGGAAACCCAAAGGGAAACCC")]
    [TestCase("GCGCGCAAAAGCGCGC")]
    [TestCase("GGGGAAAACCCCUUUUGGGGAAAACCCC")]
    public void CalculateUnpairedProbabilities_EnsembleFreeEnergy_AtMostMfe(string seq)
    {
        var r = CalculateUnpairedProbabilities(seq);
        double mfe = CalculateMinimumFreeEnergy(seq);

        Assert.That(r.EnsembleFreeEnergy, Is.LessThanOrEqualTo(mfe + 1e-9),
            "ensemble free energy -RT·ln Z ≤ MFE (same Turner-2004 model)");
    }

    #endregion

    #region Edge cases

    [Test]
    public void CalculateUnpairedProbabilities_NullOrEmpty_ReturnsUnitPartition()
    {
        var rEmpty = CalculateUnpairedProbabilities("");
        var rNull = CalculateUnpairedProbabilities(null!);

        Assert.Multiple(() =>
        {
            Assert.That(rEmpty.PartitionFunction, Is.EqualTo(1.0), "empty sequence ⇒ Z = 1 (only the empty structure)");
            Assert.That(rEmpty.BasePairProbabilities, Is.Empty, "no pairs in an empty sequence");
            Assert.That(rEmpty.EnsembleFreeEnergy, Is.EqualTo(0.0), "−RT·ln 1 = 0");
            Assert.That(rNull.PartitionFunction, Is.EqualTo(1.0), "null treated as empty ⇒ Z = 1");
        });
    }

    [Test]
    public void CalculateUnpairedProbabilities_TooShortToPair_AllUnpaired()
    {
        // Length < minLoopSize + 2 ⇒ no pair can form ⇒ only the open chain.
        var r = CalculateUnpairedProbabilities("GCGC"); // span 3 < 5

        Assert.Multiple(() =>
        {
            Assert.That(r.PartitionFunction, Is.EqualTo(1.0), "no admissible pair ⇒ Z = 1");
            Assert.That(r.BasePairProbabilities, Is.Empty, "no base pairs");
            Assert.That(r.UnpairedProbabilities, Is.All.EqualTo(1.0), "every base unpaired");
        });
    }

    [Test]
    public void CalculateUnpairedProbabilities_NonPositiveTemperature_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => CalculateUnpairedProbabilities("GGGAAACCC", temperature: 0.0),
            "temperature must be positive (Kelvin)");
    }

    // MCC-007 — single base: no pair can possibly form (span 0 < minLoopSize+2) ⇒ Z = 1, the lone
    // base is unpaired with probability 1, no base pairs, ΔG_ensemble = -RT·ln 1 = 0.
    [TestCase("G")]
    [TestCase("A")]
    public void CalculateUnpairedProbabilities_SingleBase_AllUnpaired(string seq)
    {
        var r = CalculateUnpairedProbabilities(seq);

        Assert.Multiple(() =>
        {
            Assert.That(r.PartitionFunction, Is.EqualTo(1.0), "single base ⇒ Z = 1");
            Assert.That(r.UnpairedProbabilities, Has.Count.EqualTo(1), "one entry per position");
            Assert.That(r.UnpairedProbabilities[0], Is.EqualTo(1.0), "the only base is unpaired");
            Assert.That(r.BasePairProbabilities, Is.Empty, "no base pairs");
            Assert.That(r.EnsembleFreeEnergy, Is.EqualTo(0.0), "−RT·ln 1 = 0");
        });
    }

    // MCC-008 — non-ACGU characters cannot pair (PairType returns 0 for them), so they behave like
    // permanently-unpaired positions. Here only the flanking G…C can pair; the result must stay a
    // finite, deterministic probability distribution (no NaN/throw on the ambiguity codes), and the
    // three central N's must be unpaired with probability 1.
    [Test]
    public void CalculateUnpairedProbabilities_NonAcgu_TreatedAsUnpairable()
    {
        var r1 = CalculateUnpairedProbabilities("GGGNNNCCC");
        var r2 = CalculateUnpairedProbabilities("GGGNNNCCC");

        Assert.Multiple(() =>
        {
            Assert.That(r1.PartitionFunction, Is.EqualTo(r2.PartitionFunction),
                "deterministic across calls");
            Assert.That(double.IsFinite(r1.PartitionFunction) && r1.PartitionFunction >= 1.0, Is.True,
                "Z finite and ≥ 1 even with non-ACGU bases");
            for (int i = 3; i <= 5; i++)
                Assert.That(r1.UnpairedProbabilities[i], Is.EqualTo(1.0).Within(1e-12),
                    $"non-pairable base N at {i} ⇒ p_unpaired = 1");
            foreach (var kv in r1.BasePairProbabilities)
            {
                Assert.That(kv.Key.I, Is.Not.InRange(3, 5), "no pair may involve an N position");
                Assert.That(kv.Key.J, Is.Not.InRange(3, 5), "no pair may involve an N position");
            }
        });
    }

    #endregion

    #region Region accessibility (SA)

    // MCC-005 — for GAAAC the only structure leaving the whole 5-nt window unpaired is the open
    // chain (weight 1), so P(window unpaired) = Z_open/Z = 1/Z exactly.
    [Test]
    public void CalculateRegionUnpairedProbability_GAAAC_WholeWindow_EqualsInverseZ()
    {
        var pf = CalculateUnpairedProbabilities("GAAAC");
        double expected = 1.0 / pf.PartitionFunction; // only the open chain has the window unpaired

        double pu = CalculateRegionUnpairedProbability("GAAAC", windowEnd: 4, windowLength: 5);

        Assert.That(pu, Is.EqualTo(expected).Within(1e-12),
            "P(whole window unpaired) = Z_open/Z = 1/Z (open chain is the only window-unpaired structure)");
    }

    // MCC-006 — accessibility is a probability in [0,1] and 1.0 when no pair can form.
    [Test]
    public void CalculateRegionUnpairedProbability_Bounds()
    {
        double full = CalculateRegionUnpairedProbability("GGGGAAAACCCCUUUUGGGG", windowEnd: 13, windowLength: 14);
        double tooShort = CalculateRegionUnpairedProbability("GCGC", windowEnd: 3, windowLength: 4);

        Assert.Multiple(() =>
        {
            Assert.That(full, Is.InRange(0.0, 1.0), "accessibility ∈ [0,1]");
            Assert.That(tooShort, Is.EqualTo(1.0), "no pair can form ⇒ region unpaired with probability 1");
        });
    }

    // MCC-009 — RNAplfold/Bernhart (2006) accessibility is MONOTONE NON-INCREASING in window
    // length: extending the queried region (forbidding more positions from pairing) can only shrink
    // Z_open, hence Z_open/Z cannot rise. Checked over every length 1..14 ending at a fixed anchor.
    [Test]
    public void CalculateRegionUnpairedProbability_LongerWindow_NeverMoreAccessible()
    {
        const string seq = "GGGGAAAACCCCUUUUGGGG";
        const int anchor = 13;
        double previous = double.PositiveInfinity;
        for (int len = 1; len <= 14; len++)
        {
            double p = CalculateRegionUnpairedProbability(seq, anchor, len);
            Assert.That(p, Is.LessThanOrEqualTo(previous + 1e-12),
                $"P(window len {len} unpaired) must not exceed P(window len {len - 1} unpaired)");
            Assert.That(p, Is.InRange(0.0, 1.0), "accessibility ∈ [0,1]");
            previous = p;
        }
    }

    // MCC-010 — cross-consistency of the two public methods: a length-1 region ending at i is
    // exactly "position i is unpaired", so CalculateRegionUnpairedProbability(s,i,1) must equal the
    // per-base p_unpaired(i) from CalculateUnpairedProbabilities. Both are Z_forbid(i)/Z.
    [TestCase("GGGGAAAACCCCUUUUGGGG")]
    [TestCase("GCGCGCAAAAGCGCGC")]
    public void CalculateRegionUnpairedProbability_Length1_EqualsPerBaseUnpaired(string seq)
    {
        var full = CalculateUnpairedProbabilities(seq);
        for (int i = 0; i < seq.Length; i++)
        {
            double region1 = CalculateRegionUnpairedProbability(seq, i, 1);
            Assert.That(region1, Is.EqualTo(full.UnpairedProbabilities[i]).Within(1e-9),
                $"length-1 region @ {i} = per-base p_unpaired({i})");
        }
    }

    [Test]
    public void CalculateRegionUnpairedProbability_WindowOutOfRange_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => CalculateRegionUnpairedProbability("GGGGAAAACCCC", windowEnd: 2, windowLength: 14),
                "window does not fit (start < 0)");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => CalculateRegionUnpairedProbability("GGGGAAAACCCC", windowEnd: 0, windowLength: 0),
                "window length must be positive");
        });
    }

    #endregion
}
