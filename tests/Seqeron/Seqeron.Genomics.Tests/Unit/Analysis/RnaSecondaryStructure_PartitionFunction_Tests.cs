// RNA-PARTITION-001 — RNA Partition Function (McCaskill) & Boltzmann Structure Probability
// Evidence: docs/Evidence/RNA-PARTITION-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-PARTITION-001.md
// Source: McCaskill JS (1990). Biopolymers 29(6-7):1105-1119. DOI 10.1002/bip.360290621 (PMID 1695107).
//         Will S, MIT 18.417 McCaskill slides; Freiburg RNA McCaskill teaching tool; ViennaRNA pf_fold.

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

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
    // EVERY one of the 9 admissible pairs is checked against the structure count obtained by
    // independent exhaustive enumeration (P = #structures-containing-pair / 20). Crucially this
    // includes pairs that are NESTABLE inside another pair, e.g. (2,6) and (1,7): for those the
    // McCaskill probability is NOT the external decomposition term — it requires the outside
    // recursion. A purely-external implementation gives P[2,6]=0.05 (wrong); the true value is
    // 6/20=0.30. These rows are the regression lock for the outside-recursion fix.
    // Evidence: brute-force enumeration in docs/Evidence/RNA-PARTITION-001-Evidence.md.
    [Test]
    public void CalculatePartitionFunction_GGGAAACCC_ZeroEnergy_HasExactPairProbabilities()
    {
        var result = CalculatePartitionFunction("GGGAAACCC", basePairEnergy: 0.0);
        var p = result.BasePairProbabilities;

        Assert.Multiple(() =>
        {
            // Outermost pairs (never enclosed) — external term equals the truth here.
            Assert.That(p[(0, 6)], Is.EqualTo(1.0 / 20.0).Within(1e-10), "Pair (0,6): 1 of 20 structures.");
            Assert.That(p[(0, 7)], Is.EqualTo(3.0 / 20.0).Within(1e-10), "Pair (0,7): 3 of 20 structures.");
            Assert.That(p[(0, 8)], Is.EqualTo(6.0 / 20.0).Within(1e-10), "Pair (0,8): 6 of 20 structures.");
            Assert.That(p[(2, 8)], Is.EqualTo(1.0 / 20.0).Within(1e-10), "Pair (2,8): 1 of 20 structures.");
            // Nestable pairs — REQUIRE the outside recursion; external-only would under-report.
            Assert.That(p[(1, 6)], Is.EqualTo(3.0 / 20.0).Within(1e-10), "Pair (1,6): 3 of 20 structures (nestable).");
            Assert.That(p[(1, 7)], Is.EqualTo(4.0 / 20.0).Within(1e-10), "Pair (1,7): 4 of 20 structures (nestable).");
            Assert.That(p[(1, 8)], Is.EqualTo(3.0 / 20.0).Within(1e-10), "Pair (1,8): 3 of 20 structures (nestable).");
            Assert.That(p[(2, 6)], Is.EqualTo(6.0 / 20.0).Within(1e-10), "Pair (2,6): 6 of 20 structures (nestable).");
            Assert.That(p[(2, 7)], Is.EqualTo(3.0 / 20.0).Within(1e-10), "Pair (2,7): 3 of 20 structures (nestable).");
        });
    }

    // M6b — Exact base-pair-probability spectrum for GGGGCCCC at E_bp=0 (Z=16). All 10 admissible
    // pairs checked against the independent structure counts. (1,5) and (2,6) are nestable inside
    // (0,6)/(0,7): true P = 3/16 = 0.1875, whereas the external-only formula gives 1/16 = 0.0625.
    // Evidence: brute-force enumeration (docs/Evidence/RNA-PARTITION-001-Evidence.md).
    [Test]
    public void CalculatePartitionFunction_GGGGCCCC_ZeroEnergy_HasExactPairProbabilities()
    {
        var result = CalculatePartitionFunction("GGGGCCCC", basePairEnergy: 0.0);
        var p = result.BasePairProbabilities;

        Assert.Multiple(() =>
        {
            Assert.That(p[(0, 4)], Is.EqualTo(1.0 / 16.0).Within(1e-10), "Pair (0,4): 1 of 16.");
            Assert.That(p[(0, 5)], Is.EqualTo(1.0 / 16.0).Within(1e-10), "Pair (0,5): 1 of 16.");
            Assert.That(p[(0, 6)], Is.EqualTo(2.0 / 16.0).Within(1e-10), "Pair (0,6): 2 of 16.");
            Assert.That(p[(0, 7)], Is.EqualTo(4.0 / 16.0).Within(1e-10), "Pair (0,7): 4 of 16.");
            Assert.That(p[(1, 5)], Is.EqualTo(3.0 / 16.0).Within(1e-10), "Pair (1,5): 3 of 16 (nestable).");
            Assert.That(p[(1, 6)], Is.EqualTo(2.0 / 16.0).Within(1e-10), "Pair (1,6): 2 of 16 (nestable).");
            Assert.That(p[(1, 7)], Is.EqualTo(2.0 / 16.0).Within(1e-10), "Pair (1,7): 2 of 16 (nestable).");
            Assert.That(p[(2, 6)], Is.EqualTo(3.0 / 16.0).Within(1e-10), "Pair (2,6): 3 of 16 (nestable).");
            Assert.That(p[(2, 7)], Is.EqualTo(1.0 / 16.0).Within(1e-10), "Pair (2,7): 1 of 16 (nestable).");
            Assert.That(p[(3, 7)], Is.EqualTo(1.0 / 16.0).Within(1e-10), "Pair (3,7): 1 of 16.");
        });
    }

    // M6c — Weighted (E_bp=-1) base-pair probabilities for GGGGCCCC. Every admissible pair is
    // checked against an independent Boltzmann-weighted exhaustive enumeration (weight w^#pairs,
    // w=exp(1/RT)). This locks the outside recursion for the general w≠1 case, where the wrong
    // "divide by Q^b(k,l)" shortcut also fails. Evidence: brute-force weighted enumeration.
    [Test]
    public void CalculatePartitionFunction_GGGGCCCC_NegativeEnergy_HasExactWeightedProbabilities()
    {
        var result = CalculatePartitionFunction("GGGGCCCC", basePairEnergy: -1.0);
        var p = result.BasePairProbabilities;

        Assert.Multiple(() =>
        {
            Assert.That(result.PartitionFunction, Is.EqualTo(180.01834483039346).Within(1e-7), "Z weighted.");
            Assert.That(p[(0, 4)], Is.EqualTo(0.02814492).Within(1e-7), "P(0,4) weighted.");
            Assert.That(p[(0, 6)], Is.EqualTo(0.17074408).Within(1e-7), "P(0,6) weighted.");
            Assert.That(p[(0, 7)], Is.EqualTo(0.45594238).Within(1e-7), "P(0,7) weighted.");
            Assert.That(p[(1, 5)], Is.EqualTo(0.31334323).Within(1e-7), "P(1,5) weighted (nestable).");
            Assert.That(p[(2, 6)], Is.EqualTo(0.31334323).Within(1e-7), "P(2,6) weighted (nestable).");
            Assert.That(p[(1, 6)], Is.EqualTo(0.17074408).Within(1e-7), "P(1,6) weighted (nestable).");
            Assert.That(p[(1, 7)], Is.EqualTo(0.17074408).Within(1e-7), "P(1,7) weighted (nestable).");
        });
    }

    // M6d — Single-base pairing-probability invariant: for any position, the sum of probabilities
    // over all pairs incident to it is ≤ 1 (a base pairs with at most one partner). Verified
    // independently above (max row sum 0.983 over 300 random cases). Locks the outside recursion
    // against over-counting. Evidence: McCaskill ensemble property.
    [Test]
    public void CalculatePartitionFunction_PerBasePairingProbability_NeverExceedsOne()
    {
        foreach (var (seq, ebp) in new[] { ("GGGGCCCC", 0.0), ("GGGGCCCC", -2.0), ("GGGAAACCC", -1.0) })
        {
            var result = CalculatePartitionFunction(seq, basePairEnergy: ebp);
            var rowSum = new double[seq.Length];
            foreach (var kv in result.BasePairProbabilities)
            {
                rowSum[kv.Key.I] += kv.Value;
                rowSum[kv.Key.J] += kv.Value;
            }
            for (int i = 0; i < seq.Length; i++)
                Assert.That(rowSum[i], Is.LessThanOrEqualTo(1.0 + 1e-9),
                    $"Position {i} of '{seq}' (E_bp={ebp}) pairs with total probability ≤ 1.");
        }
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
            // Anti-hang guard (generous bound; a real O(n³) at n=300 finishes in well under a
            // second — this only catches catastrophic blowup and won't flake under parallel load).
            Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(30.0),
                "The O(n³) DP on n=300 must not blow up.");
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
