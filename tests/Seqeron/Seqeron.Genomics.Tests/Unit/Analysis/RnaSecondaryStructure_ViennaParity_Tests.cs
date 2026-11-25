// RNA-STRUCT-001 — Secondary Structure Prediction: Turner 2004 interior loops + d2 loop model vs ViennaRNA
// Evidence: docs/Evidence/RNA-STRUCT-001-Evidence.md (section "Validation 2026-11: ViennaRNA parity")
// TestSpec: tests/TestSpecs/RNA-STRUCT-001.md (VP-001 … VP-008)
// Sources: NNDB Turner 2004 (turner_2004_int11_dg.txt, int21_dg.txt, int22_dg.txt, internal_loop_dg.txt,
//          terminal_mismatch_dg.txt, bulge loops page); ViennaRNA 2.7.2 sdist (eval/eval_internal.c
//          vrna_E_internal, eval/eval_exterior.c vrna_E_exterior_stem, params/default.c).
//
// Every expected number below was produced by the reference implementation ViennaRNA 2.7.2 (Python `RNA`,
// md.dangles = 2) loaded with the NNDB-precision Turner 2004 parameter file described in the Evidence
// (ViennaRNA's default file rounds the NNDB stacks to 0.1 kcal/mol and the terminal AU penalty 0.45 → 0.5,
// and differs from NNDB in 27 mismatch entries); the default-file value is quoted where it differs.

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_ViennaParity_Tests
{
    // VP-001: interior loops — int21 (both orientations), int22, 1×n, 2×3, generic, 'N' — equal
    // ViennaRNA eval_int_loop. The comment holds the probe sequence (outer pair = 1..n, inner hairpin AAAAA).
    [TestCase(1, 2, 'C', 'G', 'C', 'G', 'A', 'A', 'A', 'G', 0.80)]   // CACAAAAAGGAG   int21
    [TestCase(1, 2, 'G', 'C', 'A', 'U', 'U', 'U', 'U', 'C', 2.20)]   // GUAAAAAAUCUC   int21
    [TestCase(2, 1, 'A', 'U', 'G', 'C', 'G', 'C', 'A', 'C', 3.00)]   // AGAGAAAAACCU   int21 (2×1 orientation)
    [TestCase(1, 2, 'U', 'G', 'G', 'U', 'C', 'A', 'C', 'A', 3.70)]   // UCGAAAAAUAAG   int21, GU closures
    [TestCase(2, 2, 'C', 'G', 'G', 'C', 'G', 'G', 'A', 'A', -0.30)]  // CGAGAAAAACAGG  int22
    [TestCase(2, 2, 'A', 'U', 'U', 'A', 'U', 'U', 'U', 'U', 0.60)]   // AUUUAAAAAAUUU  int22
    [TestCase(2, 2, 'G', 'C', 'C', 'G', 'G', 'A', 'G', 'A', -1.90)]  // GGGCAAAAAGAAC  int22 tandem GA
    [TestCase(2, 2, 'G', 'U', 'A', 'U', 'A', 'G', 'N', 'G', 2.70)]   // GANAAAAAAUGGU  int22 with unknown base
    [TestCase(1, 1, 'G', 'C', 'G', 'C', 'N', 'A', 'N', 'A', 0.90)]   // GNGAAAAACAC    int11 with unknown base
    [TestCase(1, 4, 'C', 'G', 'A', 'U', 'A', 'G', 'A', 'G', 4.50)]   // CAAAAAAAUGAAGG 1×4: init(5)+asym 1.8+AU 0.7
    [TestCase(2, 3, 'A', 'U', 'G', 'C', 'G', 'A', 'A', 'C', 2.10)]   // AGAGAAAAACCAAU 2×3, 5'RG/3'YA −1.2 (default file 2.20)
    [TestCase(2, 3, 'G', 'C', 'C', 'G', 'G', 'A', 'A', 'U', 1.40)]   // GGACAAAAAGUCAC 2×3
    [TestCase(4, 4, 'C', 'G', 'G', 'C', 'G', 'G', 'G', 'U', 1.10)]   // CGAAGGAAAAACUACGG generic, GG −1.2 (default file 1.30)
    [TestCase(5, 3, 'A', 'U', 'A', 'U', 'U', 'U', 'U', 'A', 4.20)]   // AUAAAUAAAAAAUAAUU generic, UU −0.7 (default file 4.30)
    public void InteriorLoop_MatchesViennaRnaEvalIntLoop(
        int n1, int n2, char o5, char o3, char i5, char i3,
        char mm5Outer, char mm3Outer, char mm5Inner, char mm3Inner, double expected)
    {
        double e = CalculateInternalLoopEnergy(n1, n2, o5, o3, i5, i3, mm5Outer, mm3Outer, mm5Inner, mm3Inner);
        Assert.That(e, Is.EqualTo(expected).Within(1e-9));
    }

    // VP-002: interior-loop symmetry — reading the same loop from the other helix gives the same energy.
    [Test]
    public void InteriorLoop_Int21_SameLoopFromBothSides()
    {
        // CACAAAAAGGAG read as (outer C-G, inner C-G, 1×2) and its 180° rotation (outer G-C ... 2×1).
        double forward = CalculateInternalLoopEnergy(1, 2, 'C', 'G', 'C', 'G', 'A', 'A', 'A', 'G');
        double rotated = CalculateInternalLoopEnergy(2, 1, 'G', 'C', 'G', 'C', 'G', 'A', 'A', 'A');
        Assert.That(rotated, Is.EqualTo(forward).Within(1e-9));
    }

    // VP-003: degree-2 special cases route to the stack / bulge models.
    [Test]
    public void InteriorLoop_StackAndBulgeCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CalculateInternalLoopEnergy(0, 0, 'G', 'C', 'G', 'C', 'A', 'A', 'A', 'A'),
                Is.EqualTo(-3.26).Within(1e-9), "n1 = n2 = 0 is the GG/CC stack (NNDB −3.26)");
            Assert.That(CalculateInternalLoopEnergy(0, 3, 'G', 'C', 'G', 'C', 'A', 'A', 'A', 'A'),
                Is.EqualTo(3.2).Within(1e-9), "0×3 is a 3-nt bulge (NNDB initiation 3.2, GC ends)");
            Assert.That(() => CalculateInternalLoopEnergy(-1, 2, 'G', 'C', 'G', 'C', 'A', 'A', 'A', 'A'),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    // VP-004: exterior loop, d2 — single 5'/3' dangles and the terminal mismatch are read with the helix end
    // seen from the loop (ViennaRNA dangle5[type(i,j)][s(i−1)], dangle3[type][s(j+1)], mismatchExt).
    [TestCase("GGGAAACCC", -1.12)]    // no neighbour
    [TestCase("AGGGAAACCC", -1.32)]   // 5' A on G-C: −0.2
    [TestCase("GGGAAACCCA", -2.82)]   // 3' A on C(-G): −1.7
    [TestCase("AGGGAAACCCA", -2.62)]  // terminal mismatch tm(C·A / A·G) = −1.5
    public void ExteriorLoop_DanglesAndTerminalMismatch_MatchViennaRna(string seq, double expected)
    {
        Assert.That(CalculateMinimumFreeEnergy(seq), Is.EqualTo(expected).Within(1e-9));
    }

    // VP-005: MFE and MFE structure equal ViennaRNA RNA.fold (d2, NNDB-precision parameters) — random
    // sequences (seed 20261125) covering interior loops, bulges, multiloops and exterior dangles.
    private static readonly (string Seq, string Structure, double Energy, double ViennaDefault)[] FoldPanel =
    {
        ("UGGCCCUGGAUCUUUGUCGAGUAACGUGAUACAGUGGAAG", "...(((((..((..(((......))).))..))).))...", -6.96, -7.00),
        ("CCUCUCAUGGCAUCUGAGGCUCCCUUUGCUGUAAGGAUGUGGAGUGUCUCGAGUCGGCC", "((.(((..((((.((...((..((((......))))..))..))))))..)))..))..", -13.89, -13.70),
        ("AUGUUUGUGCUGAAGGAGCCGUAGCAGCCUUGAAAAGGCUCGCCGUUGCUUAAGGAGCUCACCUAAACCGUGUUGCACUGCAUAUUAUA", "((((..((((....(((((....((((((((...)))))).))....)))..(((......)))...)).....)))).))))......", -20.69, -20.40),
        ("GUAGUGGGCCAGCAUAAUCUUGGGCGGUACGGGGAAUCUGAGUGACCCAUUACCCUAAGUGCAAAU", "((((((((.((.((...((((.(......).))))...))..)).)))))))).............", -16.79, -16.90),
        ("GAAAUUGCGCGCUAAUAUGCUACUGGGACCUGACUUAGGCCGUUCACUCGGCCCGUUGUACCUUAGGUACGAAACACAUCGUCCGUU", "..........((......)).....((((........(((((......)))))..(((((((...)))))))........))))...", -21.67, -21.70),
        ("ACGACGUGGGCAAUCAGCAGUGAUAGACCCACGAUGUAUGUUAACUUUGCGAUAGGUGAUGACGCACG", "....(((((((.((((....)))).).))))))..((.(((((((((......))))..))))).)).", -18.37, -18.20),
        ("GGCGUGGUAAGAUACGCGUGUUGCCCGCACUUAACUAGC", ".(((.(((((.((....)).))))))))...........", -11.15, -11.20),
        ("UUUAACUUAUGACCCGCCUGUAGGAACAUACGCGCCACAUACGACGGUACUCGACGCUCACUU", ".......((((..((.......))..)))).(((.(...(((....)))...).)))......", -7.50, -7.40),
        ("GUUUGGUGGUCGCUACGCAUCGGGCCUAUUUCCAGCUGC", "((((((((((....)).))))))))..............", -6.87, -6.80),
        ("AGAGUAAUCCAAGCGCGUCAGACCACUUUAAA", "(((((..((...........))..)))))...", -0.64, -0.70),
    };

    [Test]
    public void Mfe_RandomPanel_MatchesViennaRnaFold()
    {
        Assert.Multiple(() =>
        {
            foreach (var (seq, structure, energy, viennaDefault) in FoldPanel)
            {
                var mfe = CalculateMfeStructure(seq);
                Assert.That(CalculateMinimumFreeEnergy(seq), Is.EqualTo(energy).Within(1e-9), seq);
                Assert.That(mfe.FreeEnergy, Is.EqualTo(energy).Within(1e-9), seq);
                Assert.That(mfe.DotBracket, Is.EqualTo(structure), seq);
                Assert.That(EvaluateStructureEnergy(seq, structure), Is.EqualTo(energy).Within(1e-9), seq);
                Assert.That(energy, Is.EqualTo(viennaDefault).Within(0.35),
                    "the residual vs ViennaRNA's default (rounded) parameter file is parameter rounding only");
            }
        });
    }

    // VP-006: the fold is self-consistent — the energy of the MFE structure re-evaluated loop by loop equals the
    // DP optimum (ViennaRNA RNA.eval_structure counterpart), over random sequences (seed 20261125).
    [Test]
    public void Mfe_EqualsEnergyOfItsOwnStructure_RandomSequences()
    {
        var rng = new Random(20261125);
        Assert.Multiple(() =>
        {
            for (int t = 0; t < 40; t++)
            {
                string seq = GenerateRandomRna(rng.Next(20, 121), rng);
                var mfe = CalculateMfeStructure(seq);
                Assert.That(EvaluateStructureEnergy(seq, mfe.DotBracket), Is.EqualTo(mfe.FreeEnergy).Within(1e-9), seq);
            }
        });
    }

    // VP-007: McCaskill partition function on the same energy model — ensemble free energy equals ViennaRNA
    // fc.pf() (d2, NNDB-precision parameters) up to the gas-constant convention (1.987 vs 1.98717 cal/mol/K).
    [TestCase("CGUCUAGCAUAGUGCUCGACGAUGAUCGUGCGACGUGC", -10.1674)]
    [TestCase("GACUCAGUCGAUGCACAUGCUAGUCAGCGCAUCAGUCG", -9.9254)]
    public void PartitionFunction_EnsembleEnergy_MatchesViennaRnaPf(string seq, double viennaEnsemble)
    {
        var r = CalculateUnpairedProbabilities(seq);
        Assert.That(r.EnsembleFreeEnergy, Is.EqualTo(viennaEnsemble).Within(0.005));
    }

    // VP-008: NNDB special C bulge needs a paired C ADJACENT to the bulged C (same strand).
    [Test]
    public void BulgeSpecialC_RequiresAdjacentPairedC()
    {
        // Outer G-C, inner A-U, bulged C on the 5' strand between G and A: its neighbours are G and A.
        double sideKnown = CalculateBulgeLoopEnergy(1, 'C', 'G', 'C', 'A', 'U', bulgeOn5PrimeStrand: true);
        double sideOther = CalculateBulgeLoopEnergy(1, 'C', 'G', 'C', 'A', 'U', bulgeOn5PrimeStrand: false);
        Assert.Multiple(() =>
        {
            Assert.That(sideKnown, Is.EqualTo(3.8 - 2.35).Within(1e-9), "init 3.8 + GA/CU stack −2.35, no special C");
            Assert.That(sideOther, Is.EqualTo(3.8 - 2.35 - 0.9).Within(1e-9), "3'-strand bulge next to the paired C: −0.9");
        });
    }
}
