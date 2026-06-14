// RNA-PARTITION-001 — RNA Partition Function (McCaskill) & Boltzmann Structure Probability
// Evidence: docs/Evidence/RNA-PARTITION-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-PARTITION-001.md
// Source: McCaskill JS (1990). Biopolymers 29(6-7):1105-1119. DOI 10.1002/bip.360290621 (PMID 1695107).
//         Will S, MIT 18.417 McCaskill slides; Freiburg RNA McCaskill teaching tool; ViennaRNA pf_fold.

using System;
using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class RnaSecondaryStructure_PartitionFunction_Tests
{
    // RT = R·T/1000 with R = 1.987 cal/(mol·K), T = 310.15 K → 0.61626805 kcal/mol (ViennaRNA k≈1.987e-3).
    private const double Rt = 1.987 * 310.15 / 1000.0;

    #region CalculatePartitionFunction

    // M1 — Z = 1 when no canonical base pair can form (AAAA has no A-U/G-C/G-U pair).
    // Evidence: McCaskill base case Q=1 (only the empty structure). INV-1.
    [Test]
    public void CalculatePartitionFunction_NoPossiblePair_ReturnsOne()
    {
        var result = CalculatePartitionFunction("AAAA", basePairEnergy: 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(result.PartitionFunction, Is.EqualTo(1.0).Within(1e-10),
                "With no canonical pair, only the empty structure exists, so Z must equal exactly 1.");
            Assert.That(result.BasePairProbabilities, Is.Empty,
                "No pair can form, so the base-pair-probability map must be empty.");
        });
    }

    // M2 — Z = 1 when the only candidate pair has span ≤ minimum hairpin loop (3).
    // GC: pair (0,1) has span 1 ≤ 3 → forbidden. Evidence: min-loop constraint (MIT slides).
    [Test]
    public void CalculatePartitionFunction_PairSpanBelowMinLoop_ReturnsOne()
    {
        var result = CalculatePartitionFunction("GC", basePairEnergy: 0.0);

        Assert.That(result.PartitionFunction, Is.EqualTo(1.0).Within(1e-10),
            "The single G-C pair has span 1 (≤ 3), so it is sterically forbidden and Z must remain 1.");
    }

    // M3 — With E_bp = 0, every Boltzmann weight is exp(0)=1, so Z equals the number of
    // admissible pseudoknot-free structures. GGGGCCCC → 16 (recurrence + exhaustive enumeration).
    // Evidence dataset (independent of implementation). A wrong/ambiguous DP would not give 16.
    [Test]
    public void CalculatePartitionFunction_GGGGCCCC_ZeroEnergy_CountsSixteenStructures()
    {
        var result = CalculatePartitionFunction("GGGGCCCC", basePairEnergy: 0.0);

        Assert.That(result.PartitionFunction, Is.EqualTo(16.0).Within(1e-10),
            "At E_bp=0, Z counts admissible structures of GGGGCCCC, which is exactly 16 (independently enumerated).");
    }

    // M4 — GGGAAACCC → 20 admissible structures at E_bp = 0. Evidence dataset.
    [Test]
    public void CalculatePartitionFunction_GGGAAACCC_ZeroEnergy_CountsTwentyStructures()
    {
        var result = CalculatePartitionFunction("GGGAAACCC", basePairEnergy: 0.0);

        Assert.That(result.PartitionFunction, Is.EqualTo(20.0).Within(1e-10),
            "At E_bp=0, Z counts admissible structures of GGGAAACCC, which is exactly 20 (independently enumerated).");
    }

    // M5 — Single-pair sequence: GAAAAC has exactly one admissible pair (0,5). Z = 2
    // ({empty},{(0,5)}); the pair's probability is 1/2. Evidence: external P formula (MIT mccaskill2).
    [Test]
    public void CalculatePartitionFunction_SingleAdmissiblePair_ProbabilityIsOneHalf()
    {
        var result = CalculatePartitionFunction("GAAAAC", basePairEnergy: 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(result.PartitionFunction, Is.EqualTo(2.0).Within(1e-10),
                "Two structures exist (empty and the single pair), so Z must be exactly 2.");
            Assert.That(result.BasePairProbabilities[(0, 5)], Is.EqualTo(0.5).Within(1e-10),
                "The pair (0,5) appears in 1 of 2 equally weighted structures, so P[0,5] = 0.5.");
            Assert.That(result.BasePairProbabilities, Has.Count.EqualTo(1),
                "Only one pair is admissible, so exactly one probability entry must be present.");
        });
    }

    // M6 — Exact base-pair-probability spectrum for GGGAAACCC at E_bp=0 (Z=20).
    // P[0,8]=6/20, P[0,7]=3/20, P[0,6]=1/20. Evidence: external P formula; counts independently derived.
    [Test]
    public void CalculatePartitionFunction_GGGAAACCC_ZeroEnergy_HasExactPairProbabilities()
    {
        var result = CalculatePartitionFunction("GGGAAACCC", basePairEnergy: 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(result.BasePairProbabilities[(0, 8)], Is.EqualTo(0.30).Within(1e-10),
                "Pair (0,8) occurs in 6 of 20 equally weighted structures → P = 0.30.");
            Assert.That(result.BasePairProbabilities[(0, 7)], Is.EqualTo(0.15).Within(1e-10),
                "Pair (0,7) occurs in 3 of 20 structures → P = 0.15.");
            Assert.That(result.BasePairProbabilities[(0, 6)], Is.EqualTo(0.05).Within(1e-10),
                "Pair (0,6) occurs in 1 of 20 structures → P = 0.05.");
        });
    }

    // M7 — Every base-pair probability must lie in [0,1] for a folding sequence (INV-2).
    [Test]
    public void CalculatePartitionFunction_FoldingSequence_AllProbabilitiesInUnitInterval()
    {
        var result = CalculatePartitionFunction("GGGAAACCC", basePairEnergy: -2.0);

        Assert.That(result.BasePairProbabilities.Values,
            Is.All.InRange(0.0, 1.0),
            "Each base-pair probability is a Boltzmann probability and must lie within [0,1].");
    }

    // M10 — null sequence throws ArgumentNullException (contract).
    [Test]
    public void CalculatePartitionFunction_NullSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => CalculatePartitionFunction(null!),
            "A null sequence is invalid input and must raise ArgumentNullException.");
    }

    // M11 — empty sequence → Z = 1 (only the empty structure), no pairs. INV-1.
    [Test]
    public void CalculatePartitionFunction_EmptySequence_ReturnsOneNoPairs()
    {
        var result = CalculatePartitionFunction("");

        Assert.Multiple(() =>
        {
            Assert.That(result.PartitionFunction, Is.EqualTo(1.0).Within(1e-10),
                "An empty sequence has exactly one (empty) structure with weight 1, so Z = 1.");
            Assert.That(result.BasePairProbabilities, Is.Empty,
                "An empty sequence has no base pairs.");
        });
    }

    // S1 — Z strictly increases as E_bp decreases (more favourable pairing). INV-4.
    [Test]
    public void CalculatePartitionFunction_LowerEnergy_StrictlyIncreasesPartitionFunction()
    {
        double z0 = CalculatePartitionFunction("GGGAAACCC", basePairEnergy: 0.0).PartitionFunction;
        double zNeg1 = CalculatePartitionFunction("GGGAAACCC", basePairEnergy: -1.0).PartitionFunction;
        double zNeg2 = CalculatePartitionFunction("GGGAAACCC", basePairEnergy: -2.0).PartitionFunction;

        Assert.Multiple(() =>
        {
            Assert.That(zNeg1, Is.GreaterThan(z0),
                "More favourable pairing (E_bp=-1 vs 0) increases every pair's Boltzmann weight, so Z must grow.");
            Assert.That(zNeg2, Is.GreaterThan(zNeg1),
                "Lowering E_bp further (−2 vs −1) must increase Z again (strict monotonicity).");
        });
    }

    // S2 — non-positive temperature throws ArgumentOutOfRangeException (contract).
    [Test]
    public void CalculatePartitionFunction_NonPositiveTemperature_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => CalculatePartitionFunction("GGGAAACCC", temperature: 0.0),
            "Absolute temperature must be positive (Kelvin); zero is physically invalid.");
    }

    // C1 — Property-based invariant test for the O(n³) DP: for many seeded random RNAs,
    // Z ≥ 1 and every base-pair probability is in [0,1] (INV-1, INV-2). Also serves as the
    // performance baseline check (length 300 case runs in-suite). Deterministic seed.
    [Test]
    public void CalculatePartitionFunction_RandomSequences_PreservesInvariants()
    {
        var rng = new Random(20260614);
        for (int trial = 0; trial < 40; trial++)
        {
            int len = rng.Next(0, 31);
            string seq = GenerateRandomRna(len, rng);

            var result = CalculatePartitionFunction(seq, basePairEnergy: -1.0);

            Assert.That(result.PartitionFunction, Is.GreaterThanOrEqualTo(1.0 - 1e-12),
                $"INV-1: Z must be ≥ 1 for any sequence (trial {trial}, seq '{seq}').");
            Assert.That(result.BasePairProbabilities.Values, Is.All.InRange(0.0, 1.0),
                $"INV-2: every base-pair probability must be in [0,1] (trial {trial}, seq '{seq}').");
        }
    }

    // Performance baseline: a length-300 RNA partition function completes quickly (O(n³)).
    [Test]
    public void CalculatePartitionFunction_Length300_CompletesUnderTimeBudget()
    {
        string seq = GenerateRandomRna(300, new Random(7));
        var sw = System.Diagnostics.Stopwatch.StartNew();

        var result = CalculatePartitionFunction(seq);

        sw.Stop();
        Assert.Multiple(() =>
        {
            Assert.That(result.PartitionFunction, Is.GreaterThanOrEqualTo(1.0 - 1e-12),
                "Z must remain ≥ 1 even at n=300.");
            Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(10.0),
                "The O(n³) DP on n=300 must complete well within the 10 s baseline budget.");
        });
    }

    #endregion

    #region CalculateStructureProbability

    // M8 — When structure energy equals ensemble energy, the structure is the whole
    // ensemble, so p = 1. Evidence: p = exp(−βE)/Z; here numerator = denominator. INV-5.
    [Test]
    public void CalculateStructureProbability_EnergyEqualsEnsemble_ReturnsOne()
    {
        double p = CalculateStructureProbability(structureEnergy: -5.0, ensembleEnergy: -5.0);

        Assert.That(p, Is.EqualTo(1.0).Within(1e-10),
            "When E_struct = E_ensemble the structure carries all the weight, so p must equal 1.");
    }

    // M9 — Boltzmann ratio exp(−(E_s − E_ens)/RT). For E_s=-5, E_ens=-6 → exp(-1/RT).
    // RT = 1.987·310.15/1000. Evidence: McCaskill p=exp(−βE)/Z; ViennaRNA β=1/kT.
    [Test]
    public void CalculateStructureProbability_KnownEnergies_ReturnsBoltzmannRatio()
    {
        double expected = Math.Exp(-1.0 / Rt); // 0.197370910785...

        double p = CalculateStructureProbability(structureEnergy: -5.0, ensembleEnergy: -6.0);

        Assert.That(p, Is.EqualTo(expected).Within(1e-10),
            "p = exp(−(E_s−E_ens)/RT) = exp(−1/RT) ≈ 0.19737 for E_s=−5, E_ens=−6 kcal/mol.");
    }

    #endregion

    #region GenerateRandomRna

    // M12 — seeded GenerateRandomRna is deterministic; correct length and RNA alphabet. INV-6.
    [Test]
    public void GenerateRandomRna_SameSeed_IsDeterministicAndValidAlphabet()
    {
        string a = GenerateRandomRna(50, new Random(42));
        string b = GenerateRandomRna(50, new Random(42));

        Assert.Multiple(() =>
        {
            Assert.That(a, Is.EqualTo(b),
                "Two generators seeded identically must produce identical sequences (determinism).");
            Assert.That(a, Has.Length.EqualTo(50),
                "Generated sequence length must equal the requested length.");
            Assert.That(a.All(c => "ACGU".Contains(c)), Is.True,
                "Generated RNA must contain only A, C, G, U bases.");
        });
    }

    // S3 — length 0 yields empty string; length 100 yields 100 RNA bases. INV-6.
    [Test]
    public void GenerateRandomRna_BoundaryLengths_ProduceCorrectOutput()
    {
        string empty = GenerateRandomRna(0, new Random(1));
        string hundred = GenerateRandomRna(100, new Random(1));

        Assert.Multiple(() =>
        {
            Assert.That(empty, Is.Empty, "Length 0 must produce an empty sequence.");
            Assert.That(hundred, Has.Length.EqualTo(100), "Length 100 must produce 100 bases.");
            Assert.That(hundred.All(c => "ACGU".Contains(c)), Is.True,
                "All generated bases must be valid RNA nucleotides.");
        });
    }

    // Delegate smoke test: parameterless-seed overload returns a valid sequence of the right length.
    [Test]
    public void GenerateRandomRna_DefaultOverload_ReturnsValidSequence()
    {
        string s = GenerateRandomRna(20);

        Assert.Multiple(() =>
        {
            Assert.That(s, Has.Length.EqualTo(20), "Default overload must honour the requested length.");
            Assert.That(s.All(c => "ACGU".Contains(c)), Is.True, "Default overload must produce valid RNA bases.");
        });
    }

    #endregion
}
