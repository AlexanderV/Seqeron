// RNA-ENERGY-001 — Turner 2004 free-energy terms (dangles, terminal mismatches, bulges, multibranch loops,
// coaxial stacking) and the shared Turner2004Parameters members.
// Evidence: docs/Evidence/RNA-ENERGY-001-Evidence.md (§8, B12 review 2026-11)
// TestSpec: tests/TestSpecs/RNA-ENERGY-001.md (EN-* rows)
// Sources: NNDB Turner 2004 (rna.urmc.rochester.edu/NNDB/rna_2004/: dangling ends, terminal mismatches, bulge
//          loops, coaxial stacking §8.2–8.4 + table turner_2004_coax, multibranch loops §9.1/9.3 + table
//          turner_2004_multi_branch); ViennaRNA 2.7.2 params/default.c (mismatchExt37 / mismatchM37 / dangle5_37 /
//          dangle3_37 'N' rows = max over A/C/G/U; ML_closing37 = 930, ML_intern37 = −90, ML_BASE37 = 0),
//          sequences/alphabet.c (T ≡ U, case-insensitive, other characters → N); Mathews et al. 2004 PNAS 101:7287.
// Oracle: ViennaRNA 2.7.2 Python RNA (default parameters, and the NNDB-precision parameter file of
//         RNA-STRUCT-001 Evidence §12.2 where NNDB stacking precision matters).

using Seqeron.Genomics.Core;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_Energy_Tests
{
    #region EN-1 Unknown bases (N / IUPAC / gap) in exterior and multibranch helix-end terms

    // ViennaRNA 2.7.2 default.c: the 'N' entries are the maximum over A/C/G/U.
    // GetTerminalMismatchEnergy(c5, c3, x, y) = mismatchExt37[pair(c3, c5)][y][x].
    [TestCase('C', 'G', 'N', 'N', -0.8, TestName = "EN-1 TM C·G NN = mismatchExt37[GC][N][N] -80")]
    [TestCase('G', 'C', 'A', 'N', -1.1, TestName = "EN-1 TM G·C AN = mismatchExt37[CG][N][A] -110")]
    [TestCase('A', 'U', 'R', 'G', -0.6, TestName = "EN-1 TM A·U RG (IUPAC R = unknown) = mismatchExt37[UA][G][N] -60")]
    [TestCase('U', 'G', '-', 'C', -0.6, TestName = "EN-1 TM U·G -C (gap = unknown) = mismatchExt37[GU][C][N] -60")]
    public void TerminalMismatch_UnknownBase_IsViennaMaxOverBases(char c5, char c3, char x, char y, double expected)
    {
        Assert.That(GetTerminalMismatchEnergy(c5, c3, x, y), Is.EqualTo(expected).Within(1e-9));
    }

    // GetDanglingEndEnergy(c5, c3, d, 3') = dangle3_37[pair(c3, c5)][d]; 5' = dangle5_37[pair(c3, c5)][d].
    [TestCase('C', 'G', 'N', true, -0.8, TestName = "EN-1 dangle3 on C·G, N = dangle3_37[GC][N] -80")]
    [TestCase('G', 'C', 'N', false, -0.1, TestName = "EN-1 dangle5 on G·C, N = dangle5_37[CG][N] -10")]
    [TestCase('A', 'U', 'Y', true, -0.5, TestName = "EN-1 dangle3 on A·U, Y = dangle3_37[UA][N] -50")]
    public void DanglingEnd_UnknownBase_IsViennaMaxOverBases(char c5, char c3, char d, bool is3Prime, double expected)
    {
        Assert.That(GetDanglingEndEnergy(c5, c3, d, is3Prime), Is.EqualTo(expected).Within(1e-9));
    }

    [Test]
    public void TerminalMismatchAndDangle_DnaTAndLowerCase_AreReadAsRna()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetDanglingEndEnergy('A', 'T', 'A', true), Is.EqualTo(-0.8).Within(1e-9), "dangle3_37[UA][A]");
            Assert.That(GetDanglingEndEnergy('g', 'c', 't', true), Is.EqualTo(-0.6).Within(1e-9), "dangle3_37[CG][U]");
            Assert.That(GetTerminalMismatchEnergy('a', 't', 'g', 't'), Is.EqualTo(-1.0).Within(1e-9), "AGUU = mismatchExt37[UA][U][G]");
            Assert.That(GetTerminalMismatchEnergy('G', 'T', 'G', 'G'), Is.EqualTo(-0.8).Within(1e-9), "GGGU = mismatchExt37[UG][G][G]");
        });
    }

    [Test]
    public void TerminalMismatchAndDangle_NonCanonicalClosingPair_ReturnZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetTerminalMismatchEnergy('A', 'A', 'G', 'G'), Is.EqualTo(0.0));
            Assert.That(GetDanglingEndEnergy('A', 'A', 'A', true), Is.EqualTo(0.0));
            Assert.That(GetDanglingEndEnergy('N', 'G', 'A', false), Is.EqualTo(0.0));
        });
    }

    // Engine impact: unknown neighbours of an exterior helix end. RNA.fold, NNDB-precision file:
    // NGGGAAACCCN -1.92 (default file -2.00); before the fix the N mismatch scored 0 (-1.12).
    [TestCase("NGGGAAACCCN", ".(((...))).", -1.92)]
    [TestCase("GGGAAACCCN", "(((...))).", -1.92)]
    [TestCase("NNGGGGAAACCCCNN", "..((((...))))..", -5.18)]
    public void Mfe_UnknownNeighbours_MatchViennaRnaFold(string seq, string mfeStructure, double viennaMfe)
    {
        var r = CalculateMfeStructure(seq);
        Assert.Multiple(() =>
        {
            Assert.That(r.DotBracket, Is.EqualTo(mfeStructure));
            Assert.That(r.FreeEnergy, Is.EqualTo(viennaMfe).Within(1e-9));
            Assert.That(CalculateMinimumFreeEnergy(seq), Is.EqualTo(viennaMfe).Within(1e-9));
        });
    }

    #endregion

    #region EN-2 Bulge loops: special-C strand, DNA T, size validation

    [Test]
    public void Bulge_SideUnknown_AppliesSpecialCOnlyWhenBothStrandsQualify()
    {
        // Outer G·C, inner A·U, bulged C: a paired C is adjacent only on the 3' strand (C of G·C).
        // Strand unknown → the least stabilising reading = ViennaRNA eval_int_loop (no special C): 3.8 − 2.35.
        Assert.Multiple(() =>
        {
            Assert.That(CalculateBulgeLoopEnergy(1, 'C', 'G', 'C', 'A', 'U'), Is.EqualTo(1.45).Within(1e-9));
            Assert.That(CalculateBulgeLoopEnergy(1, 'C', 'G', 'C', 'A', 'U', bulgeOn5PrimeStrand: false),
                Is.EqualTo(0.55).Within(1e-9));
            // C·G outer, G·C inner: a paired C is adjacent on both strands → bonus whatever the strand.
            Assert.That(CalculateBulgeLoopEnergy(1, 'C', 'C', 'G', 'G', 'C'), Is.EqualTo(0.54).Within(1e-9), "3.8 − 2.36 − 0.9");
            // C·G outer, C·G inner: the 3' strand neighbours are both G → no bonus without the strand.
            Assert.That(CalculateBulgeLoopEnergy(1, 'C', 'C', 'G', 'C', 'G'), Is.EqualTo(0.54).Within(1e-9), "3.8 − 3.26");
            Assert.That(CalculateBulgeLoopEnergy(1, 'C', 'C', 'G', 'C', 'G', bulgeOn5PrimeStrand: true),
                Is.EqualTo(-0.36).Within(1e-9), "3.8 − 3.26 − 0.9");
        });
    }

    [Test]
    public void Bulge_DnaTAndLowerCase_AreReadAsRna()
    {
        Assert.Multiple(() =>
        {
            // n = 1: stack GA/CU −2.35 (was not found for "GA/CT").
            Assert.That(CalculateBulgeLoopEnergy(1, 'a', 'G', 'C', 'A', 'T', bulgeOn5PrimeStrand: true), Is.EqualTo(1.45).Within(1e-9));
            // n = 2: A·T closure takes the AU/GU penalty: 2.8 + 0.45 (ViennaRNA eval_int_loop, NNDB-precision file 325).
            Assert.That(CalculateBulgeLoopEnergy(2, 'A', 'A', 'T', 'G', 'C'), Is.EqualTo(3.25).Within(1e-9));
        });
    }

    [TestCase(0)]
    [TestCase(-1)]
    public void Bulge_SizeBelowOne_Throws(int size)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateBulgeLoopEnergy(size, 'A', 'G', 'C', 'G', 'C'));
    }

    // NNDB extrapolation beyond the table: init(6) + 1.75·R·T·ln(n/6) (R = 1.987, T = 310.15).
    [Test]
    public void Bulge_LargerThan30_UsesNndbExtrapolation()
    {
        double expected = Math.Round(4.4 + 1.75 * 1.987 * 310.15 / 1000.0 * Math.Log(40 / 6.0), 2); // 6.45
        Assert.That(CalculateBulgeLoopEnergy(40, 'A', 'G', 'C', 'G', 'C'), Is.EqualTo(expected).Within(1e-9));
    }

    #endregion

    #region EN-3 Multibranch loops: linear (engine) model and NNDB efn2 model

    // ViennaRNA 2.7.2 eval_structure_verbose, dangles = 2, NNDB-precision file:
    // GAGGGAAAACCCAGGGAAAACCCAC (.(((....))).(((....))).) "Multi loop" = 250 = 660 + three terminal mismatches.
    [Test]
    public void Multibranch_Linear_PlusTerminalMismatches_MatchesViennaRnaD2()
    {
        const string s = "GAGGGAAAACCCAGGGAAAACCCAC";
        double stacking = GetTerminalMismatchEnergy(s[0], s[24], s[1], s[23])      // closing pair, inside
                        + GetTerminalMismatchEnergy(s[11], s[2], s[12], s[1])      // branch 3–12
                        + GetTerminalMismatchEnergy(s[22], s[13], s[23], s[12]);   // branch 14–23
        Assert.Multiple(() =>
        {
            Assert.That(stacking, Is.EqualTo(-4.1).Within(1e-9));
            Assert.That(CalculateMultibranchLoopEnergy(3, 3, stackingEnergy: stacking), Is.EqualTo(2.5).Within(1e-9));
        });
    }

    // NNDB efn2: 9.25 + 0.91·min(2, mean asymmetry) − 0.63·h (+3.14 strain for 3-way loops with < 2 unpaired).
    [TestCase(new[] { 2, 0, 1, 5 }, 0.0, 8.55, TestName = "EN-3 efn2 4-way, mean asymmetry 2.5 capped at 2 (NNDB §9.3: 8.6)")]
    [TestCase(new[] { 2, 0, 1, 5 }, -5.9, 2.65, TestName = "EN-3 efn2 NNDB §9.3 total with stacking -5.9 (NNDB: 2.7)")]
    [TestCase(new[] { 2, 2, 2 }, 0.0, 7.36, TestName = "EN-3 efn2 symmetric 3-way, asymmetry 0")]
    [TestCase(new[] { 3, 0, 0 }, 0.0, 9.18, TestName = "EN-3 efn2 3-way, asymmetry (3+3+0)/3 = 2, no strain")]
    [TestCase(new[] { 0, 0, 0 }, 0.0, 10.5, TestName = "EN-3 efn2 3-way, 0 unpaired: strain 3.14")]
    public void Multibranch_Efn2_MatchesNndbFormula(int[] gaps, double stacking, double expected)
    {
        Assert.That(CalculateMultibranchLoopEnergy(gaps, stacking), Is.EqualTo(expected).Within(1e-9));
    }

    [Test]
    public void Multibranch_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMultibranchLoopEnergy(2, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMultibranchLoopEnergy(3, -1));
            Assert.Throws<ArgumentException>(() => CalculateMultibranchLoopEnergy(3, 1, hasStrain: true));
            Assert.Throws<ArgumentException>(() => CalculateMultibranchLoopEnergy(new[] { 1, 1 }));
            Assert.Throws<ArgumentException>(() => CalculateMultibranchLoopEnergy(new[] { 1, -1, 0 }));
            Assert.Throws<ArgumentNullException>(() => CalculateMultibranchLoopEnergy((IReadOnlyList<int>)null!));
        });
    }

    #endregion

    #region EN-4 Coaxial stacking (NNDB §8)

    [TestCase('A', 'U', 'G', 'G', -2.9, TestName = "EN-4 mismatch coax NNDB §8.4: AU + GG mismatch -0.8 -2.1 = -2.9")]
    [TestCase('G', 'C', 'A', 'U', -4.0, TestName = "EN-4 mismatch coax: AU-capable mismatch tm(GAUC) -1.5 -2.1 -0.4")]
    [TestCase('C', 'G', 'G', 'U', -3.8, TestName = "EN-4 mismatch coax: GU-capable mismatch tm(CGUG) -1.5 -2.1 -0.2")]
    [TestCase('G', 'U', 'A', 'A', -2.4, TestName = "EN-4 mismatch coax: GU closing, AA mismatch: no bonus")]
    [TestCase('g', 'c', 'a', 't', -4.0, TestName = "EN-4 mismatch coax: lower case / T")]
    [TestCase('A', 'A', 'G', 'G', 0.0, TestName = "EN-4 mismatch coax: non-canonical closing pair = 0")]
    public void MismatchCoaxialStacking_MatchesNndbRule(char c5, char c3, char x, char y, double expected)
    {
        Assert.That(CalculateMismatchCoaxialStacking(c5, c3, x, y), Is.EqualTo(expected).Within(1e-9));
    }

    [Test]
    public void FlushCoaxialStacking_NndbExampleAndDnaT()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CalculateFlushCoaxialStacking('A', 'U', 'C', 'G'), Is.EqualTo(-2.24).Within(1e-9), "NNDB §8.4");
            Assert.That(CalculateFlushCoaxialStacking('a', 't', 'c', 'g'), Is.EqualTo(-2.24).Within(1e-9));
            Assert.That(CalculateFlushCoaxialStacking('G', 'U', 'U', 'G'), Is.EqualTo(1.29).Within(1e-9), "GU/UG stack");
        });
    }

    #endregion

    #region EN-6 CalculateFreeEnergy (structure evaluation, RNAeval counterpart)

    // ViennaRNA 2.7.2 RNA.fold_compound(s).eval_structure(db), md.dangles = 2, NNDB-precision parameter file.
    [TestCase("GGGAAACCC", "(((...)))", -1.12, TestName = "EN-6 eval hairpin GGGAAACCC")]
    [TestCase("GGGAUAAAUCCC", "((((....))))", -3.42, TestName = "EN-6 eval NNDB hairpin example 1 (inner AU)")]
    [TestCase("GGGAAACCCAGGGAAACCC", "(((...))).(((...)))", -4.14, TestName = "EN-6 eval two exterior stems")]
    [TestCase("GAGGGAAAACCCAGGGAAAACCCAC", "(.(((....))).(((....))).)", -1.54, TestName = "EN-6 eval 3-way multiloop d2")]
    [TestCase("NGGGAAACCCN", ".(((...))).", -1.92, TestName = "EN-6 eval unknown exterior neighbours")]
    [TestCase("gggaaaccc", "(((...)))", -1.12, TestName = "EN-6 eval lower case")]
    public void CalculateFreeEnergy_MatchesViennaRnaEvalStructure(string seq, string db, double expected)
    {
        Assert.That(CalculateFreeEnergy(seq, db), Is.EqualTo(expected).Within(1e-9));
    }

    [TestCase("GGGAAACCCAGGGAAACCCAAAGCGCAAAGCGCU")]
    [TestCase("ACGUACGUAGCUAGCUAGCGAUCGAUCGAUCGAUGC")]
    public void CalculateFreeEnergy_OfMfeStructure_EqualsMfe(string seq)
    {
        var mfe = CalculateMfeStructure(seq);
        Assert.That(CalculateFreeEnergy(seq, mfe.DotBracket), Is.EqualTo(mfe.FreeEnergy).Within(1e-9));
    }

    [Test]
    public void CalculateFreeEnergy_InvalidStructure_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => CalculateFreeEnergy("GGGAAACCC", "[[[...]]]"), "pseudoknot brackets");
            Assert.Throws<ArgumentException>(() => CalculateFreeEnergy("GGGAAACCC", "(((...))"), "length");
            Assert.Throws<ArgumentException>(() => CalculateFreeEnergy("GGGAAACCC", "((....)))"), "unbalanced");
            Assert.Throws<ArgumentException>(() => CalculateFreeEnergy("AAAAAAAAA", "(((...)))"), "non-canonical pair");
            Assert.Throws<ArgumentNullException>(() => CalculateFreeEnergy(null!, "."));
        });
    }

    #endregion

    #region EN-5 Turner2004Parameters

    [Test]
    public void TerminalPenaltyPair_IsCaseInsensitiveAndReadsTAsU()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Turner2004Parameters.IsTerminalPenaltyPair('a', 't'), Is.True);
            Assert.That(Turner2004Parameters.IsTerminalPenaltyPair('T', 'g'), Is.True);
            Assert.That(Turner2004Parameters.IsTerminalPenaltyPair('G', 'C'), Is.False);
            Assert.That(Turner2004Parameters.IsTerminalPenaltyPair('N', 'U'), Is.False);
        });
    }

    // Table sizes = NNDB Turner 2004 (16 WC + 20 GU stacks, hairpin initiation 3..30, 96 terminal mismatches).
    [Test]
    public void Tables_HaveNndbEntryCounts()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Turner2004Parameters.StackingEnergies, Has.Count.EqualTo(36));
            Assert.That(Turner2004Parameters.HairpinLoopInitiation.Keys, Is.EquivalentTo(Enumerable.Range(3, 28)));
            Assert.That(Turner2004Parameters.TerminalMismatchEnergies, Has.Count.EqualTo(96));
            Assert.That(Turner2004Parameters.TerminalAuGuPenalty, Is.EqualTo(0.45));
        });
    }

    #endregion
}
