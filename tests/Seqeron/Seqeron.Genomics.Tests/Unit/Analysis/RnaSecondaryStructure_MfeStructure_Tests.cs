// RNA-STRUCT-001 — Secondary Structure Prediction (MFE-optimal structure via DP traceback)
// Evidence: docs/Evidence/RNA-STRUCT-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-STRUCT-001.md
// Source: Zuker M, Stiegler P (1981). Nucleic Acids Res. 9(1):133-148 (W/V matrices, MFE traceback).
//         MIT 6.047 Lecture 08 (Washietl), Fig. 13 — F/C/M/M¹ recurrence decomposition.
//         Nussinov R, Jacobson AB (1980). PNAS 77(11):6309-6313 (DP traceback principle).
//
// These tests cover CalculateMfeStructure / PredictStructureMfe — the DP-traceback structure that
// is the optimal counterpart of the greedy PredictStructure. The defining contract (verified here)
// is that the reconstructed structure's free energy EQUALS the scalar MFE the DP already returns,
// because the traceback follows the SAME V/W/WM recurrences. The scalar MFE itself is validated
// against exact NNDB Turner 2004 values by RnaSecondaryStructure_MinimumFreeEnergy_Tests.

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_MfeStructure_Tests
{
    #region CalculateMfeStructure — energy consistency (the core invariant)

    // INV — Zuker-Stiegler traceback consistency: the reconstructed structure's energy MUST equal
    // the scalar MFE for the SAME input and parameters. This is the property the whole fix exists
    // to guarantee, and it would FAIL for the old greedy PredictStructure (whose energy is a sum of
    // independent stem-loops, not the DP optimum). Tested across hairpin, multi-stem and multiloop.
    [TestCase("GGGAAACCC")]
    [TestCase("GGGGAAAACCCC")]
    [TestCase("CACAAAAAAAUGUG")]
    [TestCase("GCAGCAAAAGCGC")]
    [TestCase("GGGAAACCCAAAGGGAAACCC")]
    [TestCase("GGGGAAAACCCCUUUUGGGGAAAACCCC")]
    [TestCase("GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA")]
    public void CalculateMfeStructure_ReconstructedEnergy_EqualsScalarMfe(string rna)
    {
        double scalar = CalculateMinimumFreeEnergy(rna);

        var structure = CalculateMfeStructure(rna);

        Assert.That(structure.FreeEnergy, Is.EqualTo(scalar).Within(1e-9),
            $"Traceback structure energy must equal the scalar DP MFE for '{rna}' " +
            "(same recurrences ⇒ same energy; a divergence would mean the traceback is not " +
            "consistent with the DP).");
    }

    // M1 — Simple hairpin GGGAAACCC: the only stabilizing fold is the 3-bp stem closing a 3-nt loop.
    // Optimal dot-bracket is (((...))) and the energy matches the NNDB-validated scalar MFE (-1.12).
    [Test]
    public void CalculateMfeStructure_SimpleHairpin_ReturnsNestedStemAndLoop()
    {
        var structure = CalculateMfeStructure("GGGAAACCC");

        Assert.Multiple(() =>
        {
            Assert.That(structure.DotBracket, Is.EqualTo("(((...)))"),
                "MFE structure of GGGAAACCC is the 3-bp hairpin with a 3-nt loop.");
            Assert.That(structure.FreeEnergy, Is.EqualTo(-1.12).Within(1e-9),
                "Energy equals the NNDB-validated scalar MFE (2 GC/GC stacks -3.42·... per RNA-MFE-001).");
            Assert.That(structure.BasePairs, Is.EquivalentTo(new[] { (0, 8), (1, 7), (2, 6) }),
                "Pairs are the three antiparallel WC pairs of the stem (5'<3', sorted).");
        });
    }

    // M2 — Four-pair GC hairpin GGGGAAAACCCC → ((((....)))) (4-bp stem, 4-nt loop), energy -5.28.
    [Test]
    public void CalculateMfeStructure_FourPairHairpin_ReturnsFourBpStem()
    {
        var structure = CalculateMfeStructure("GGGGAAAACCCC");

        Assert.Multiple(() =>
        {
            Assert.That(structure.DotBracket, Is.EqualTo("((((....))))"),
                "MFE structure is the 4-bp hairpin closing a 4-nt loop.");
            Assert.That(structure.FreeEnergy, Is.EqualTo(-5.28).Within(1e-9),
                "Energy equals the NNDB-validated scalar MFE for this hairpin (RNA-MFE-001).");
        });
    }

    // M3 — GGGAAACCC|AAA|GGGAAACCC. The outer GGG (0-2) is complementary to the final CCC (18-20),
    // so the DP nests an inner hairpin inside an outer stem — a single nested fold that is more
    // stable than two side-by-side hairpins (an extra GC stack stabilizes the outer helix). The DP
    // finds this global optimum; energy equals the NNDB-validated scalar MFE.
    [Test]
    public void CalculateMfeStructure_NestableSequence_FindsNestedFold()
    {
        var structure = CalculateMfeStructure("GGGAAACCCAAAGGGAAACCC");

        Assert.Multiple(() =>
        {
            Assert.That(structure.DotBracket, Is.EqualTo("(((...(((...)))...)))"),
                "Outer stem (0-2 with 18-20) enclosing an inner hairpin (6-8 with 12-14).");
            Assert.That(ValidateDotBracket(structure.DotBracket), Is.True,
                "Notation must be balanced/well-nested.");
            Assert.That(structure.FreeEnergy, Is.EqualTo(CalculateMinimumFreeEnergy("GGGAAACCCAAAGGGAAACCC")).Within(1e-9),
                "Energy equals the scalar DP MFE.");
        });
    }

    // M4 — Optimality over greedy: GGGGAAAACCCCUUUUGGGGAAAACCCC. The greedy PredictStructure picks a
    // single hairpin and misses the deeper fold; the DP traceback finds the NESTED two-helix optimum
    // ((((....((((....))))....)))) with strictly lower (more negative) energy. This is the exact
    // defect RNA-STRUCT-001 fixes: optimal structure ≠ greedy structure.
    [Test]
    public void CalculateMfeStructure_OutperformsGreedy_FindsNestedOptimum()
    {
        const string rna = "GGGGAAAACCCCUUUUGGGGAAAACCCC";

        var optimal = CalculateMfeStructure(rna);
        var greedy = PredictStructure(rna, minStemLength: 3, minLoopSize: 3, maxLoopSize: 10);

        Assert.Multiple(() =>
        {
            Assert.That(optimal.DotBracket, Is.EqualTo("((((....((((....))))....))))"),
                "DP traceback recovers the nested two-helix MFE structure.");
            Assert.That(optimal.FreeEnergy, Is.EqualTo(-13.76).Within(1e-9),
                "Optimal energy equals the scalar DP MFE.");
            Assert.That(optimal.FreeEnergy, Is.LessThan(greedy.MinimumFreeEnergy),
                "The DP optimum is strictly more stable than the greedy stem-loop sum.");
        });
    }

    // M5 — Multiloop / branched structure (tRNA-like) exercises the WM (multiloop) traceback path.
    // The reconstructed structure must be valid, full-length, and energy-consistent with the scalar.
    [Test]
    public void CalculateMfeStructure_TrnaLike_BranchedStructureIsConsistent()
    {
        const string trna =
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA";

        var structure = CalculateMfeStructure(trna);

        Assert.Multiple(() =>
        {
            Assert.That(structure.DotBracket.Length, Is.EqualTo(trna.Length),
                "Dot-bracket length equals sequence length.");
            Assert.That(ValidateDotBracket(structure.DotBracket), Is.True,
                "Branched (multiloop) notation must be balanced.");
            Assert.That(structure.FreeEnergy, Is.EqualTo(CalculateMinimumFreeEnergy(trna)).Within(1e-9),
                "Branched-structure energy equals the scalar DP MFE.");
            Assert.That(structure.BasePairs.Count, Is.GreaterThan(0),
                "A 73-nt tRNA folds into multiple helices, not the open chain.");
        });
    }

    #endregion

    #region CalculateMfeStructure — notation/pair consistency and invariants

    // M6 — The base pairs returned and the pairs implied by the dot-bracket must be identical sets,
    // so the two representations never disagree.
    [Test]
    public void CalculateMfeStructure_BasePairs_MatchDotBracketPairs()
    {
        var structure = CalculateMfeStructure("GGGAAACCCAAAGGGAAACCC");

        var fromPairs = new HashSet<(int, int)>(structure.BasePairs);
        var fromNotation = new HashSet<(int, int)>(ParseDotBracket(structure.DotBracket));

        Assert.That(fromPairs.SetEquals(fromNotation), Is.True,
            "BasePairs and the pairs parsed from DotBracket must be the same set.");
    }

    // M7 — Every returned pair must be a valid WC or wobble pair (the DP only pairs CanPair bases).
    [Test]
    public void CalculateMfeStructure_AllPairs_AreWatsonCrickOrWobble()
    {
        var structure = CalculateMfeStructure(
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA");

        Assert.That(structure.BasePairs.All(p =>
                CanPair(structure.Sequence[p.Position1], structure.Sequence[p.Position2])),
            Is.True, "Every optimal base pair must be a canonical WC or G-U wobble pair.");
    }

    // M8 — Optimal pseudoknot-free DP can never produce crossing pairs: i<k<j<l never occurs.
    [Test]
    public void CalculateMfeStructure_OptimalStructure_HasNoCrossingPairs()
    {
        var structure = CalculateMfeStructure("GGGGAAAACCCCUUUUGGGGAAAACCCC");

        var pairs = structure.BasePairs.OrderBy(p => p.Position1).ToList();
        bool anyCrossing = false;
        for (int a = 0; a < pairs.Count; a++)
            for (int b = a + 1; b < pairs.Count; b++)
            {
                var (i, j) = pairs[a];
                var (k, l) = pairs[b];
                if (i < k && k < j && j < l) anyCrossing = true; // crossing ⇒ pseudoknot
            }

        Assert.That(anyCrossing, Is.False,
            "A pseudoknot-free DP optimum must have only nested or disjoint pairs (no crossings).");
    }

    #endregion

    #region CalculateMfeStructure — edge cases

    // E1 — Null sequence ⇒ empty structure, no pairs, ΔG = 0 (matches CalculateMinimumFreeEnergy).
    [Test]
    public void CalculateMfeStructure_Null_ReturnsEmptyStructure()
    {
        var structure = CalculateMfeStructure(null!);

        Assert.Multiple(() =>
        {
            Assert.That(structure.DotBracket, Is.EqualTo(""), "Null ⇒ empty dot-bracket.");
            Assert.That(structure.BasePairs, Is.Empty, "Null ⇒ no base pairs.");
            Assert.That(structure.FreeEnergy, Is.EqualTo(0.0).Within(1e-12), "Null ⇒ ΔG = 0.");
        });
    }

    // E2 — Empty sequence ⇒ empty structure, ΔG = 0.
    [Test]
    public void CalculateMfeStructure_Empty_ReturnsEmptyStructure()
    {
        var structure = CalculateMfeStructure("");

        Assert.Multiple(() =>
        {
            Assert.That(structure.DotBracket, Is.EqualTo(""), "Empty ⇒ empty dot-bracket.");
            Assert.That(structure.BasePairs, Is.Empty, "Empty ⇒ no base pairs.");
            Assert.That(structure.FreeEnergy, Is.EqualTo(0.0).Within(1e-12), "Empty ⇒ ΔG = 0.");
        });
    }

    // E3 — Too short to hold a hairpin (len < minLoopSize+2): all-dots, no pairs, ΔG = 0.
    [Test]
    public void CalculateMfeStructure_TooShort_ReturnsAllDots()
    {
        var structure = CalculateMfeStructure("GC");

        Assert.Multiple(() =>
        {
            Assert.That(structure.DotBracket, Is.EqualTo(".."),
                "Too short ⇒ all positions unpaired.");
            Assert.That(structure.BasePairs, Is.Empty, "Too short ⇒ no base pairs.");
            Assert.That(structure.FreeEnergy, Is.EqualTo(0.0).Within(1e-12), "Too short ⇒ ΔG = 0.");
        });
    }

    // E4 — A homopolymer has no complementary bases: the open chain (all dots, ΔG = 0) is optimal.
    [Test]
    public void CalculateMfeStructure_PolyA_ReturnsOpenChain()
    {
        var structure = CalculateMfeStructure("AAAAAAAAAAAA");

        Assert.Multiple(() =>
        {
            Assert.That(structure.DotBracket, Is.EqualTo("............"),
                "Poly-A cannot pair ⇒ open chain.");
            Assert.That(structure.BasePairs, Is.Empty, "Poly-A ⇒ no base pairs.");
            Assert.That(structure.FreeEnergy, Is.EqualTo(0.0).Within(1e-12), "Poly-A ⇒ ΔG = 0.");
        });
    }

    // E5 — Case insensitivity and DNA: lowercase RNA and the T-as-U convention fold identically to
    // the uppercase RNA form (T and U are interchangeable for folding).
    [Test]
    public void CalculateMfeStructure_LowercaseAndDna_FoldLikeUppercaseRna()
    {
        var reference = CalculateMfeStructure("GGGAAACCC");
        var lower = CalculateMfeStructure("gggaaaccc");
        var dna = CalculateMfeStructure("GGGTTTCCC"); // TTT read as UUU loop

        Assert.Multiple(() =>
        {
            Assert.That(lower.DotBracket, Is.EqualTo(reference.DotBracket),
                "Lowercase input folds like uppercase.");
            Assert.That(lower.FreeEnergy, Is.EqualTo(reference.FreeEnergy).Within(1e-9),
                "Lowercase energy equals uppercase energy.");
            Assert.That(dna.DotBracket, Is.EqualTo(reference.DotBracket),
                "DNA (T) folds like the RNA (U) form.");
        });
    }

    #endregion

    #region PredictStructureMfe — optimal-path overload

    // M9 — PredictStructureMfe returns the DP optimum in the SecondaryStructure shape, with the same
    // dot-bracket and energy as CalculateMfeStructure, and no pseudoknots (DP is pseudoknot-free).
    [Test]
    public void PredictStructureMfe_ReturnsOptimalStructureInSecondaryStructureShape()
    {
        const string rna = "GGGGAAAACCCCUUUUGGGGAAAACCCC";
        var direct = CalculateMfeStructure(rna);

        var prediction = PredictStructureMfe(rna);

        Assert.Multiple(() =>
        {
            Assert.That(prediction.DotBracket, Is.EqualTo(direct.DotBracket),
                "PredictStructureMfe dot-bracket equals CalculateMfeStructure.");
            Assert.That(prediction.MinimumFreeEnergy, Is.EqualTo(direct.FreeEnergy).Within(1e-9),
                "PredictStructureMfe energy equals the DP MFE.");
            Assert.That(prediction.BasePairs.Count, Is.EqualTo(direct.BasePairs.Count),
                "Same number of base pairs as the traceback result.");
            Assert.That(prediction.Pseudoknots, Is.Empty,
                "The pseudoknot-free DP optimum contains no crossing pairs.");
        });
    }

    // M10 — PredictStructureMfe base pairs carry the correct bases and pair types.
    [Test]
    public void PredictStructureMfe_BasePairs_CarryCorrectBasesAndTypes()
    {
        var prediction = PredictStructureMfe("GGGAAACCC");

        Assert.Multiple(() =>
        {
            Assert.That(prediction.BasePairs.Count, Is.EqualTo(3), "3-bp stem.");
            foreach (var bp in prediction.BasePairs)
            {
                Assert.That(bp.Type, Is.EqualTo(BasePairType.WatsonCrick),
                    "All three G-C pairs are Watson-Crick.");
                Assert.That(prediction.Sequence[bp.Position1], Is.EqualTo(bp.Base1));
                Assert.That(prediction.Sequence[bp.Position2], Is.EqualTo(bp.Base2));
            }
        });
    }

    #endregion

    #region Property-based test (O(n³) algorithm — required by Definition of Done)

    // P1 — For random RNA sequences, the traceback structure's energy ALWAYS equals the scalar MFE,
    // the dot-bracket is ALWAYS balanced and full-length, and pairs are ALWAYS canonical. Fixed seed
    // ⇒ deterministic. This is the invariant property of a DP traceback (Zuker-Stiegler 1981).
    [Test]
    public void CalculateMfeStructure_RandomSequences_AlwaysConsistentValidCanonical()
    {
        var rng = new Random(20260622); // fixed, documented seed → deterministic

        Assert.Multiple(() =>
        {
            for (int t = 0; t < 50; t++)
            {
                int len = 12 + rng.Next(0, 40); // 12..51 nt
                string rna = GenerateRandomRna(len, rng, gcContent: 0.5);

                double scalar = CalculateMinimumFreeEnergy(rna);
                var structure = CalculateMfeStructure(rna);

                Assert.That(structure.FreeEnergy, Is.EqualTo(scalar).Within(1e-9),
                    $"[{rna}] traceback energy must equal scalar MFE.");
                Assert.That(structure.DotBracket.Length, Is.EqualTo(rna.Length),
                    $"[{rna}] dot-bracket length must equal sequence length.");
                Assert.That(ValidateDotBracket(structure.DotBracket), Is.True,
                    $"[{rna}] dot-bracket must be balanced/well-nested.");
                Assert.That(structure.BasePairs.All(p =>
                        CanPair(structure.Sequence[p.Position1], structure.Sequence[p.Position2])),
                    Is.True, $"[{rna}] every pair must be a canonical WC/wobble pair.");
            }
        });
    }

    #endregion
}
