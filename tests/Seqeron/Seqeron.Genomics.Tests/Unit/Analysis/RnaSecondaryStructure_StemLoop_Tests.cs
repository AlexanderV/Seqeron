// RNA-STEMLOOP-001 — Stem-loop (hairpin) detection and its Turner 2004 energies
// Evidence: docs/Evidence/RNA-STEMLOOP-001-Evidence.md (§11, B12 review 2026-11)
// TestSpec: tests/TestSpecs/RNA-STEMLOOP-001.md (SL-E*, SL-U* rows)
// Sources: NNDB Turner 2004 hairpin loops (rna.urmc.rochester.edu/NNDB/turner04/hairpin.html, worked
//          examples §5.3; special hairpin loop table), NNDB helix rules (per AU/GU end 0.45);
//          ViennaRNA 2.7.2 eval/eval_hairpin.c (vrna_E_hairpin returns Hexaloop_E), params/default.c
//          (mismatchH37 'N' row/column = max over A/C/G/U), sequences/alphabet.c (T ≡ U, other → N);
//          Zuker & Stiegler (1981) Nucleic Acids Res 9:133 (V(i,j) recursion, degree-2 loops).
// Oracle: ViennaRNA 2.7.2 Python RNA.fold_compound(...).eval_structure / eval_hp_loop / mfe with the
//         NNDB-precision parameter file of RNA-STRUCT-001 Evidence §12.2 (file b); default-file values in Evidence.

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_StemLoop_Tests
{
    private static StemLoop Single(string seq, string dotBracket, int maxInternalLoopSize = 0) =>
        FindStemLoops(seq, minStemLength: 1, minLoopSize: 3, maxLoopSize: 30, allowWobble: true, maxInternalLoopSize)
            .Single(sl => sl.Start == 0 && sl.End == seq.Length - 1 && sl.DotBracketNotation == dotBracket);

    #region SL-E1 NNDB worked examples (Turner 2004 hairpin page §5.3)

    // ΔG°37 = helix (stacks + AU/GU end penalties) + hairpin loop, exactly as NNDB adds them.
    [TestCase("CACAAAAAAAUGUG", "((((......))))", -1.41, TestName = "SL-E1 6-nt loop, no special terms: NNDB -1.4")]
    [TestCase("CACAGAAAGUGUG", "((((.....))))", -1.91, TestName = "SL-E1 5-nt loop, GG first mismatch: NNDB -1.9")]
    [TestCase("CACCCGAGGGUG", "((((....))))", -4.11, TestName = "SL-E1 special tetraloop CcgagG: NNDB -4.1")]
    [TestCase("CACACCCCCCUGUG", "((((......))))", 2.09, TestName = "SL-E1 6-nt all-C loop: NNDB 2.1")]
    [TestCase("CGGGGAAAGUCCG", "((((.....))))", -4.80, TestName = "SL-E1 special GU closure: NNDB -4.8")]
    public void FindStemLoops_NndbWorkedExamples_TotalFreeEnergy(string seq, string db, double expected)
    {
        var sl = Single(seq, db);
        Assert.That(sl.TotalFreeEnergy, Is.EqualTo(expected).Within(0.005));
    }

    #endregion

    #region SL-E2 A-U-closed special hexaloop: measured total already includes the closure penalty

    // Before the fix FindStemLoops added the +0.45 inner AU end on top of the measured hexaloop
    // (GGGACAGUGUUCCC −6.62, CGACAGUACUCG −1.46). ViennaRNA eval_structure (file b): −7.07 / −1.91.
    [TestCase("GGGACAGUGUUCCC", "((((......))))", -7.07, -8.42)]
    [TestCase("CGACAGUACUCG", "(((......)))", -1.91, -4.26)]
    public void FindStemLoops_AuClosedSpecialHexaloop_MatchesViennaEval(string seq, string db, double total, double helix)
    {
        var sl = Single(seq, db);
        Assert.Multiple(() =>
        {
            Assert.That(sl.TotalFreeEnergy, Is.EqualTo(total).Within(0.005), "ViennaRNA eval_structure");
            Assert.That(sl.Stem.FreeEnergy, Is.EqualTo(helix).Within(0.005), "NNDB helix (CalculateStemEnergy)");
            Assert.That(sl.TotalFreeEnergy, Is.EqualTo(EvaluateStructureEnergy(seq, db)).Within(1e-9),
                "TotalFreeEnergy = shared structure evaluator on the isolated motif");
        });
    }

    #endregion

    #region SL-E3 Unknown bases (N) and DNA T in the hairpin loop

    // ViennaRNA encodes any non-ACGU symbol as N; mismatchH37 N entries = max over A/C/G/U.
    // Before: N first/last mismatch scored 0 (GGGAAANCCC −0.92, GGGNAAANCCC −0.82, GGGCNNNNGCCC −4.34).
    [TestCase("GGGAAANCCC", "(((....)))", -2.02)]
    [TestCase("GGGNAAANCCC", "(((.....)))", -1.32)]
    [TestCase("GGGCNNNNGCCC", "((((....))))", -5.14)]
    public void FindStemLoops_UnknownBaseMismatch_MatchesViennaEval(string seq, string db, double expected)
    {
        Assert.That(Single(seq, db).TotalFreeEnergy, Is.EqualTo(expected).Within(0.005));
    }

    [Test]
    public void CalculateHairpinLoopEnergy_UnknownAndDnaBases_MatchViennaEvalHpLoop()
    {
        Assert.Multiple(() =>
        {
            // ViennaRNA eval_hp_loop("GAAANC") = 450 dcal/mol.
            Assert.That(CalculateHairpinLoopEnergy("AAAN", 'G', 'C'), Is.EqualTo(4.50).Within(0.005));
            // Special tetraloop CUUCGG 3.7 (NNDB; eval_hp_loop 370) — DNA / lower-case spelling must hit it too
            // (before: "TTCG" scored as a generic loop, 5.6).
            Assert.That(CalculateHairpinLoopEnergy("UUCG", 'C', 'G'), Is.EqualTo(3.7).Within(0.005));
            Assert.That(CalculateHairpinLoopEnergy("ttcg", 'c', 'g'), Is.EqualTo(3.7).Within(0.005));
        });
    }

    [Test]
    public void FindStemLoops_DnaAndLowerCase_SameAsRna()
    {
        var rna = FindStemLoops("GGGACAGUGUUCCC", 2, 3, 10).ToList();
        var dna = FindStemLoops("gggacagtgttccc", 2, 3, 10).ToList();
        Assert.That(dna.Select(s => (s.Start, s.End, s.DotBracketNotation, s.TotalFreeEnergy, s.Stem.FreeEnergy)),
            Is.EqualTo(rna.Select(s => (s.Start, s.End, s.DotBracketNotation, s.TotalFreeEnergy, s.Stem.FreeEnergy))));
    }

    #endregion

    #region SL-E4 CalculateStemEnergy reads the bases (case, T, pair class)

    [Test]
    public void CalculateStemEnergy_LowerCaseDnaAndMislabelledType_ScoredFromBases()
    {
        static BasePair P(int i, int j, char a, char b, BasePairType t) => new(i, j, a, b, t);
        Assert.Multiple(() =>
        {
            // NNDB GG/CC = −3.26 (before: lower-case keys missed the table → 0).
            Assert.That(CalculateStemEnergy("", [P(0, 9, 'g', 'c', BasePairType.WatsonCrick), P(1, 8, 'g', 'c', BasePairType.WatsonCrick)]),
                Is.EqualTo(-3.26).Within(0.005));
            // NNDB AA/UU −0.93 + 2 × AU end 0.45 = −0.03 (before: A·T → no stack, no AU penalty → 0).
            Assert.That(CalculateStemEnergy("", [P(0, 9, 'A', 'T', BasePairType.WatsonCrick), P(1, 8, 'A', 'T', BasePairType.WatsonCrick)]),
                Is.EqualTo(-0.03).Within(0.005));
            // G·U mislabelled WatsonCrick still gets its GU end penalty: GC/UG... CG/GU stack + 0.45.
            Assert.That(CalculateStemEnergy("", [P(0, 9, 'C', 'G', BasePairType.WatsonCrick), P(1, 8, 'G', 'U', BasePairType.WatsonCrick)]),
                Is.EqualTo(CalculateStemEnergy("", [P(0, 9, 'C', 'G', BasePairType.WatsonCrick), P(1, 8, 'G', 'U', BasePairType.Wobble)])));
        });
    }

    #endregion

    #region SL-U Stems with bulges / internal loops (opt-in maxInternalLoopSize)

    // GGG A CAG GAAA CUG CCC: a 1-nt bulge splits two 3-bp helices. ViennaRNA RNA.fold MFE (file b) is
    // exactly this stem-loop: (((.(((....)))))) −7.23 (default file −7.30).
    [Test]
    public void FindStemLoops_WithBulges_ReturnsOptimalUnbranchedStem_MatchingViennaMfe()
    {
        const string seq = "GGGACAGGAAACUGCCC";
        var perfect = FindStemLoops(seq, 3, 4, 4).Single();
        var bulged = FindStemLoops(seq, 3, 4, 4, allowWobble: true, maxInternalLoopSize: 1).Single();
        Assert.Multiple(() =>
        {
            Assert.That(perfect.DotBracketNotation, Is.EqualTo("(((....)))"));
            Assert.That(perfect.TotalFreeEnergy, Is.EqualTo(-1.09).Within(0.005));
            Assert.That(bulged.DotBracketNotation, Is.EqualTo("(((.(((....))))))"));
            Assert.That((bulged.Start, bulged.End, bulged.Stem.Length), Is.EqualTo((0, 16, 6)));
            Assert.That(bulged.TotalFreeEnergy, Is.EqualTo(-7.23).Within(0.005));
            Assert.That(bulged.Loop.Sequence, Is.EqualTo("GAAA"));
            Assert.That(bulged.Stem.FreeEnergy + CalculateHairpinLoopEnergy("GAAA", 'G', 'C'),
                Is.EqualTo(bulged.TotalFreeEnergy).Within(0.005));
        });
    }

    // maxInternalLoopSize = 0 (default) keeps the perfect-helix enumeration bit for bit.
    [Test]
    public void FindStemLoops_MaxInternalLoopSizeZero_IsBackwardCompatible()
    {
        const string seq = "GGGAGCUCAGAUAUCCGGGACAGUGUUCCCAUGCAUGGAAAC";
        var a = FindStemLoops(seq, 2, 3, 12, true).ToList();
        var b = FindStemLoops(seq, 2, 3, 12, true, maxInternalLoopSize: 0).ToList();
        Assert.That(b.Select(s => (s.Start, s.End, s.DotBracketNotation, s.TotalFreeEnergy, s.Stem.FreeEnergy)),
            Is.EqualTo(a.Select(s => (s.Start, s.End, s.DotBracketNotation, s.TotalFreeEnergy, s.Stem.FreeEnergy))));
    }

    // Oracle values: for every hairpin loop, min over unbranched chains of AU/GU end + Σ ViennaRNA
    // eval_int_loop + eval_hp_loop (file b); the reported structure re-evaluates to the same value.
    [TestCase("ACGUACUGACGUCAGCAUGCUAGCAGUCAG", 4)]
    [TestCase("GCAUGCUACGUAGCAUCGAUCAGCUAGCUAGC", 8)]
    public void FindStemLoops_Unbranched_EveryCandidateEqualsItsStructureEnergy(string seq, int maxInternal)
    {
        foreach (var sl in FindStemLoops(seq, 1, 3, 10, true, maxInternal))
        {
            string sub = seq.Substring(sl.Start, sl.End - sl.Start + 1);
            Assert.That(sl.TotalFreeEnergy, Is.EqualTo(EvaluateStructureEnergy(sub, sl.DotBracketNotation)).Within(1e-9));
            Assert.That(sl.Loop.Size, Is.InRange(3, 10));
            Assert.That(sl.Stem.BasePairs.All(bp => CanPair(bp.Base1, bp.Base2)), Is.True);
        }
    }

    // Exact oracle rows for GCAUGCUACGUAGCAUCGAUCAGCUAGCUAGC, maxInternalLoopSize = 8 (75/75 hairpins agree;
    // Evidence §11.4): loop [start..end], motif span, dot-bracket, ΔG = optimum over chains = eval_structure.
    [TestCase(14, 21, 4, 31, "((((..((((........))))..))))", -8.86)]
    [TestCase(10, 15, 0, 27, "((..(((.((......))...)))..))", -3.48)]
    [TestCase(7, 10, 2, 15, "(((((....)))))", -3.01)]
    [TestCase(3, 9, 0, 13, "(((.......).))", 4.12)]
    [TestCase(25, 28, 22, 31, "(((....)))", -0.25)]
    public void FindStemLoops_Unbranched_MatchesViennaOptimumOverChains(
        int loopStart, int loopEnd, int start, int end, string db, double expected)
    {
        var sl = FindStemLoops("GCAUGCUACGUAGCAUCGAUCAGCUAGCUAGC", 1, 3, 10, true, maxInternalLoopSize: 8)
            .Single(s => s.Loop.Start == loopStart && s.Loop.End == loopEnd);
        Assert.Multiple(() =>
        {
            Assert.That((sl.Start, sl.End, sl.DotBracketNotation), Is.EqualTo((start, end, db)));
            Assert.That(sl.TotalFreeEnergy, Is.EqualTo(expected).Within(0.005));
        });
    }

    [Test]
    public void FindStemLoops_NegativeMaxInternalLoopSize_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => FindStemLoops("GGGAAACCC", maxInternalLoopSize: -1));
    }

    #endregion
}
