// MIRNA-PAIR-001 — MiRNA-Target Pairing Analysis
// Evidence: docs/Evidence/MIRNA-PAIR-001-Evidence.md
// TestSpec: tests/TestSpecs/MIRNA-PAIR-001.md
// Source: Agarwal V et al. (2015) eLife 4:e05005 (A-U/C-G); Crick FHC (1966) J Mol Biol 19(2):548-555 (G-U wobble);
//         Lewis BP et al. (2005) Cell 120(1):15-20 (reverse-complement seed); NNDB Turner 2004 duplex model
//         (rna.urmc.rochester.edu/NNDB/rna_2004: helices, dangling ends, terminal mismatches, internal loops);
//         ViennaRNA 2.8.0 eval_structure cross-check; Biopython 1.88 reverse_complement_rna.

using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

[TestFixture]
public class MiRnaAnalyzer_AlignMiRnaToTarget_Tests
{
    #region CanPair

    // M1 — Watson-Crick pairs A-U, U-A, G-C, C-G all pair. Evidence: Agarwal et al. (2015) "A pairs with U, C pairs with G".
    [Test]
    public void CanPair_WatsonCrickPairs_ReturnsTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CanPair('A', 'U'), Is.True, "A-U is a Watson-Crick pair");
            Assert.That(CanPair('U', 'A'), Is.True, "U-A is a Watson-Crick pair");
            Assert.That(CanPair('G', 'C'), Is.True, "G-C is a Watson-Crick pair");
            Assert.That(CanPair('C', 'G'), Is.True, "C-G is a Watson-Crick pair");
        });
    }

    // M2 — G:U wobble pairs. Evidence: Crick (1966) wobble hypothesis; G-U is the standard RNA wobble pair.
    [Test]
    public void CanPair_WobblePairs_ReturnsTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CanPair('G', 'U'), Is.True, "G-U is the standard RNA wobble pair (Crick 1966)");
            Assert.That(CanPair('U', 'G'), Is.True, "U-G is the standard RNA wobble pair (Crick 1966)");
        });
    }

    // M3 — Non-pairing combinations. Evidence: Agarwal et al. — A pairs only with U, C only with G.
    [Test]
    public void CanPair_NonPairs_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CanPair('A', 'A'), Is.False, "A-A does not pair");
            Assert.That(CanPair('A', 'C'), Is.False, "A pairs only with U, not C");
            Assert.That(CanPair('A', 'G'), Is.False, "A pairs only with U, not G");
            Assert.That(CanPair('C', 'U'), Is.False, "C pairs only with G, not U");
            Assert.That(CanPair('C', 'C'), Is.False, "C-C does not pair");
            Assert.That(CanPair('G', 'G'), Is.False, "G-G does not pair");
        });
    }

    // M4 — Case-insensitive and DNA T handled (T treated as U for pairing).
    // Evidence: T↔U RNA normalisation contract (MiRNA_Target_Pairing.md §3.1, §6.1; TestSpec M4 "A-T → true").
    [Test]
    public void CanPair_LowercaseAndDnaT_ReturnsTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CanPair('a', 'u'), Is.True, "Lowercase a-u must pair (case-insensitive)");
            Assert.That(CanPair('g', 'c'), Is.True, "Lowercase g-c must pair (case-insensitive)");
            // DNA T is normalised to RNA U, so A-T must pair like A-U (contract §3.1, §6.1).
            Assert.That(CanPair('A', 'T'), Is.True, "DNA A-T must pair (T normalised to U)");
            Assert.That(CanPair('T', 'A'), Is.True, "DNA T-A must pair (T normalised to U)");
            Assert.That(CanPair('G', 'T'), Is.True, "G-T is a wobble (T normalised to U → G-U)");
            // Mixed case + DNA T together.
            Assert.That(CanPair('a', 't'), Is.True, "Lowercase a-t must pair (case + T→U)");
        });
    }

    // M4b — IsWobblePair also honours the DNA T→U contract: G-T / T-G are wobble.
    // Evidence: contract §3.1 "T treated via U where relevant"; Crick (1966) G-U wobble.
    [Test]
    public void IsWobblePair_DnaT_TreatedAsU()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IsWobblePair('G', 'T'), Is.True, "G-T is a wobble (T→U)");
            Assert.That(IsWobblePair('T', 'G'), Is.True, "T-G is a wobble (T→U)");
            Assert.That(IsWobblePair('A', 'T'), Is.False, "A-T (A-U) is Watson-Crick, not wobble");
        });
    }

    #endregion

    #region IsWobblePair

    // M5 — Only G:U / U:G are wobble. Evidence: PMC4870184 distinguishes G:U wobble from Watson-Crick.
    [Test]
    public void IsWobblePair_GU_ReturnsTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IsWobblePair('G', 'U'), Is.True, "G-U is a wobble pair");
            Assert.That(IsWobblePair('U', 'G'), Is.True, "U-G is a wobble pair");
        });
    }

    // M6 — Watson-Crick pairs are NOT wobble. Evidence: PMC4870184 — solid lines (WC) vs wobble (G:U).
    [Test]
    public void IsWobblePair_WatsonCrick_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IsWobblePair('A', 'U'), Is.False, "A-U is Watson-Crick, not wobble");
            Assert.That(IsWobblePair('G', 'C'), Is.False, "G-C is Watson-Crick, not wobble");
            Assert.That(IsWobblePair('A', 'A'), Is.False, "A-A is not a pair at all");
        });
    }

    #endregion

    #region GetReverseComplement

    // M7 — let-7a-5p seed (pos 2-8) GAGGUAG reverse complement = CUACCUC.
    // Evidence: Lewis (2005) targets are reverse complement of seed; A-U/G-C complementarity (Agarwal).
    [Test]
    public void GetReverseComplement_Let7aSeed_ReturnsExpected()
    {
        // GAGGUAG: complement G->C,A->U,G->C,G->C,U->A,A->U,G->C = CUCCAUC; reversed = CUACCUC
        string rc = GetReverseComplement("GAGGUAG");

        Assert.That(rc, Is.EqualTo("CUACCUC"), "Reverse complement of let-7a seed GAGGUAG must be CUACCUC");
    }

    // M8 — DNA input: T complements to A; output uses RNA U. Evidence: A-U/G-C; T->U convention.
    [Test]
    public void GetReverseComplement_DnaInput_ReturnsRnaReverseComplement()
    {
        // ACGT: complement A->U,C->G,G->C,T->A = UGCA; reversed = ACGU
        string rc = GetReverseComplement("ACGT");

        Assert.That(rc, Is.EqualTo("ACGU"), "RC of DNA ACGT must be ACGU (T complements to A, RNA U used)");
    }

    // M9 — Empty input returns empty (defensive contract).
    [Test]
    public void GetReverseComplement_Empty_ReturnsEmpty()
    {
        Assert.That(GetReverseComplement(""), Is.Empty, "RC of empty sequence is empty");
    }

    // M8b — IUPAC codes complemented, other symbols pass through (canonical GetRnaComplementBase).
    // Oracle: Biopython 1.88 Bio.Seq.reverse_complement_rna("ACGR") = "YCGU",
    // reverse_complement_rna("NNRYKMBDHVSW") = "WSBDHVKMRYNN", reverse_complement_rna("acgu-X") = "X-acgu"
    // (the canonical helper upper-cases recognised bases: "X-ACGU").
    [Test]
    public void GetReverseComplement_IupacAndUnknown_MatchesBiopython()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GetReverseComplement("ACGR"), Is.EqualTo("YCGU"), "R complements to Y (was N before)");
            Assert.That(GetReverseComplement("NNRYKMBDHVSW"), Is.EqualTo("WSBDHVKMRYNN"));
            Assert.That(GetReverseComplement("acgu-X"), Is.EqualTo("X-ACGU"), "gap and unknown symbol pass through");
        });
    }

    // C1 — Double reverse complement restores the RNA sequence (involution). Evidence: INV-04.
    [Test]
    public void GetReverseComplement_DoubleApplication_RestoresRna()
    {
        const string seq = "GAGGUAGUAGG";

        string twice = GetReverseComplement(GetReverseComplement(seq));

        Assert.That(twice, Is.EqualTo(seq), "RC applied twice must restore the original RNA sequence");
    }

    #endregion

    #region AlignMiRnaToTarget

    // M10 — Perfect Watson-Crick complement: every position pairs canonically.
    // Evidence: Watson-Crick A-U (Agarwal). Antiparallel: miRNA A[i] pairs target U.
    [Test]
    public void AlignMiRnaToTarget_PerfectComplement_AllWatsonCrickMatches()
    {
        var duplex = AlignMiRnaToTarget("AAAA", "UUUU");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.Matches, Is.EqualTo(4), "All 4 A-U positions are Watson-Crick matches");
            Assert.That(duplex.Mismatches, Is.EqualTo(0), "No mismatches for a perfect complement");
            Assert.That(duplex.GUWobbles, Is.EqualTo(0), "No wobbles for an A-U duplex");
            Assert.That(duplex.AlignmentString, Is.EqualTo("||||"), "All positions marked '|' (Watson-Crick)");
        });
    }

    // M11 — G:U duplex counted as wobbles, not matches. Evidence: Crick (1966); G:U != Watson-Crick (PMC4870184).
    [Test]
    public void AlignMiRnaToTarget_GUWobbleDuplex_CountedAsWobbles()
    {
        var duplex = AlignMiRnaToTarget("GGGG", "UUUU");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.GUWobbles, Is.EqualTo(4), "All 4 G-U positions are wobble pairs");
            Assert.That(duplex.Matches, Is.EqualTo(0), "G:U must NOT be counted as Watson-Crick matches");
            Assert.That(duplex.Mismatches, Is.EqualTo(0), "G:U pairs are valid pairs, not mismatches");
            Assert.That(duplex.AlignmentString, Is.EqualTo("::::"), "All positions marked ':' (wobble)");
        });
    }

    // M12 — Non-pairing duplex: A vs A cannot pair. Evidence: A pairs only with U (Agarwal).
    [Test]
    public void AlignMiRnaToTarget_NonPairingDuplex_AllMismatches()
    {
        var duplex = AlignMiRnaToTarget("AAAA", "AAAA");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.Mismatches, Is.EqualTo(4), "A-A cannot pair → 4 mismatches");
            Assert.That(duplex.Matches, Is.EqualTo(0), "No Watson-Crick matches");
            Assert.That(duplex.GUWobbles, Is.EqualTo(0), "No wobbles");
            Assert.That(duplex.AlignmentString, Is.EqualTo("    "), "All positions are spaces (mismatch)");
        });
    }

    // M13 — Empty miRNA returns empty duplex (defensive contract).
    [Test]
    public void AlignMiRnaToTarget_EmptyMiRna_ReturnsEmptyDuplex()
    {
        var duplex = AlignMiRnaToTarget("", "AAAA");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.AlignmentString, Is.Empty, "Empty miRNA yields empty alignment");
            Assert.That(duplex.Matches, Is.EqualTo(0), "No matches for empty input");
            Assert.That(duplex.Mismatches, Is.EqualTo(0), "No mismatches for empty input");
            Assert.That(duplex.GUWobbles, Is.EqualTo(0), "No wobbles for empty input");
        });
    }

    // M14 — Empty target returns empty duplex (defensive contract).
    [Test]
    public void AlignMiRnaToTarget_EmptyTarget_ReturnsEmptyDuplex()
    {
        var duplex = AlignMiRnaToTarget("AAAA", "");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.AlignmentString, Is.Empty, "Empty target yields empty alignment");
            Assert.That(duplex.Matches, Is.EqualTo(0), "No matches for empty target");
        });
    }

    // M15 — Count invariant on a genuinely mixed duplex exercising all three classes.
    // miRNA AGGU vs target AUCG, antiparallel (miRNA[i] pairs target[3-i]):
    //   i0: A-G → mismatch ' ';  i1: G-C → Watson-Crick '|';  i2: G-U → wobble ':';  i3: U-A → Watson-Crick '|'.
    // Exact counts derive from the pairing rules (Agarwal A-U/C-G; Crick G-U); INV-05 sum = overlap.
    [Test]
    public void AlignMiRnaToTarget_CountInvariant_SumEqualsOverlapLength()
    {
        const string mirna = "AGGU";
        const string target = "AUCG";
        var duplex = AlignMiRnaToTarget(mirna, target);
        int overlap = System.Math.Min(mirna.Length, target.Length);

        Assert.Multiple(() =>
        {
            Assert.That(duplex.Matches, Is.EqualTo(2), "G-C and U-A are the two Watson-Crick matches");
            Assert.That(duplex.GUWobbles, Is.EqualTo(1), "G-U is the single wobble");
            Assert.That(duplex.Mismatches, Is.EqualTo(1), "A-G is the single mismatch");
            Assert.That(duplex.AlignmentString, Is.EqualTo(" |:|"), "Alignment symbols per position");
            Assert.That(duplex.Matches + duplex.Mismatches + duplex.GUWobbles, Is.EqualTo(overlap),
                "Every overlap position is classified exactly once (INV-05)");
            Assert.That(duplex.Gaps, Is.EqualTo(0), "Aligner is ungapped (INV-05)");
            Assert.That(duplex.AlignmentString, Has.Length.EqualTo(overlap), "Alignment string spans the overlap");
        });
    }

    // S4 — Unequal lengths align over the shorter overlap. Evidence: INV-05.
    [Test]
    public void AlignMiRnaToTarget_UnequalLengths_AlignsOverOverlap()
    {
        var duplex = AlignMiRnaToTarget("AAAAAA", "UUUU");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.AlignmentString, Has.Length.EqualTo(4), "Overlap = min(6,4) = 4 positions");
            Assert.That(duplex.Matches, Is.EqualTo(4), "First 4 A-U positions pair (Watson-Crick)");
        });
    }

    // S3 — DNA input normalised to RNA before pairing (T→U). Evidence: T↔U normalisation contract.
    [Test]
    public void AlignMiRnaToTarget_DnaInput_NormalisedAndMatched()
    {
        var duplex = AlignMiRnaToTarget("AAAA", "TTTT");

        Assert.That(duplex.Matches, Is.EqualTo(4), "DNA T normalised to U → 4 A-U Watson-Crick matches");
    }

    // S1 — Fully Watson-Crick self-complementary duplex: exact NNDB Turner 2004 sum.
    // GCGCGCGC is its own reverse complement, so both strands are identical (NNDB symmetry term applies):
    // 4.09 (initiation) + 0.43 (symmetry) + 4 × (5'GC/3'CG −3.42) + 3 × (5'CG/3'GC −2.36) = −16.24.
    // ViennaRNA 2.8 eval_structure of GCGCGCGC&GCGCGCGC "(((((((())))))))" = −16.70 (0.1-rounded stacks, no symmetry term).
    [Test]
    public void AlignMiRnaToTarget_FullyPairedDuplex_FreeEnergyIsNndbSum()
    {
        string mirna = "GCGCGCGC";
        string target = GetReverseComplement(mirna);

        var duplex = AlignMiRnaToTarget(mirna, target);

        Assert.Multiple(() =>
        {
            Assert.That(duplex.Matches, Is.EqualTo(8), "Fully Watson-Crick paired duplex (8 pairs)");
            Assert.That(duplex.FreeEnergy, Is.EqualTo(-16.24).Within(1e-9),
                "NNDB Turner 2004: init 4.09 + symmetry 0.43 + 4×(−3.42) + 3×(−2.36)");
        });
    }

    // S2 — All-mismatch duplex: no base pair ⇒ no duplex ⇒ ΔG = 0 (ViennaRNA eval of the open structure = 0).
    [Test]
    public void AlignMiRnaToTarget_AllMismatchDuplex_FreeEnergyZero()
    {
        var duplex = AlignMiRnaToTarget("AAAA", "AAAA");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.Mismatches, Is.EqualTo(4), "No pairs form");
            Assert.That(duplex.FreeEnergy, Is.EqualTo(0.0), "No base pair → no duplex → ΔG = 0");
        });
    }

    // E1 — NNDB Turner 2004 Watson-Crick helix page, worked example "Self complementary duplex"
    // 5'AGCGCU3' / 3'UCGCGA5': 4.09 + 2×0.45 + 0.43 − 2.08 − 3.42 − 2.36 − 3.42 − 2.08 = −7.94 kcal/mol.
    [Test]
    public void AlignMiRnaToTarget_NndbSelfComplementaryExample_Minus7_94()
    {
        var duplex = AlignMiRnaToTarget("AGCGCU", "AGCGCU");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.AlignmentString, Is.EqualTo("||||||"));
            Assert.That(duplex.FreeEnergy, Is.EqualTo(-7.94).Within(1e-9), "NNDB worked example −7.94 kcal/mol");
        });
    }

    // E2 — NNDB worked example "Non-self complementary duplex" 5'GCACG3' / 3'CGUGC5':
    // 4.09 − 3.42 − 2.11 − 2.24 − 2.36 = −6.04 kcal/mol (no AU end, no symmetry term).
    [Test]
    public void AlignMiRnaToTarget_NndbNonSelfComplementaryExample_Minus6_04()
    {
        var duplex = AlignMiRnaToTarget("GCACG", "CGUGC");

        Assert.That(duplex.FreeEnergy, Is.EqualTo(-6.04).Within(1e-9), "NNDB worked example −6.04 kcal/mol");
    }

    // E3 — Dangling end: the E2 duplex with a 3'-dangling C on the miRNA (3' of the terminal G·C pair).
    // NNDB dangling-end page / table: 3' dangling C adjacent to GC = −0.4 → −6.04 − 0.4 = −6.44.
    // ViennaRNA 2.8 eval_structure GCACGC&CGUGC "(((((.)))))" = −6.40.
    [Test]
    public void AlignMiRnaToTarget_ThreePrimeDanglingEnd_NndbValue()
    {
        var duplex = AlignMiRnaToTarget("GCACGC", "CGUGC");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.AlignmentString, Is.EqualTo("|||||"), "overlap 5; miRNA nt 6 overhangs");
            Assert.That(duplex.FreeEnergy, Is.EqualTo(-6.44).Within(1e-9), "−6.04 + 3' dangling C on GC (−0.4)");
        });
    }

    // E4 — Terminal mismatch: NNDB terminal-mismatch page example extends the E1 helix by a
    // "UA followed by GA" mismatch (−1.1). Here the two strands differ (AGCGCUG vs AAGCGCU), so the
    // 0.43 symmetry term of E1 does not apply: −7.94 − 0.43 − 1.1 = −9.47 (NNDB's −9.0 keeps the symmetry term).
    // ViennaRNA 2.8 eval_structure = −9.40.
    [Test]
    public void AlignMiRnaToTarget_TerminalMismatch_NndbValue()
    {
        var duplex = AlignMiRnaToTarget("AGCGCUG", "AAGCGCU");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.AlignmentString, Is.EqualTo("|||||| "), "G·A at the 3' end of the miRNA does not pair");
            Assert.That(duplex.FreeEnergy, Is.EqualTo(-9.47).Within(1e-9), "−7.94 − 0.43 + terminal mismatch UA/GA (−1.1)");
        });
    }

    // E5 — A short A·U duplex is NOT stabilising: the intermolecular initiation (+4.09) outweighs three
    // 5'AA/3'UU stacks (−0.93 each) and two terminal AU penalties (+0.45 each): 4.09 − 2.79 + 0.90 = +2.20.
    // ViennaRNA 2.8 eval_structure AAAA&UUUU "(((())))" = +2.40 (init 4.10, stack −0.90, AU 0.50).
    // With two extra miRNA nt the 3'-dangling A on the A·U end adds −0.8 (NNDB dangle table) → +1.40 (Vienna +1.60).
    [Test]
    public void AlignMiRnaToTarget_ShortAuDuplex_PositiveNndbEnergy()
    {
        var duplex = AlignMiRnaToTarget("AAAA", "UUUU");
        var overhang = AlignMiRnaToTarget("AAAAAA", "UUUU");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.FreeEnergy, Is.EqualTo(2.20).Within(1e-9), "4.09 + 3×(−0.93) + 2×0.45");
            Assert.That(overhang.FreeEnergy, Is.EqualTo(1.40).Within(1e-9), "+2.20 + 3' dangling A on AU (−0.8)");
        });
    }

    // E6 — Cross-check against ViennaRNA 2.8.0: fold_compound(miRNA + "&" + target).eval_structure(s), where
    // s pairs every in-register Watson-Crick/G·U position of the reported alignment (dangles = 2).
    // Residual = Vienna's 0.1-rounded parameters (stacks, init 4.10, terminal AU 0.50); max |Δ| 0.46 over
    // 51 fixtures (0.43 of it = NNDB symmetry term on GCGCGCGC·GCGCGCGC, which eval_structure omits), 0.45 on
    // fixtures with 2×2 loops (was 2.2 before the canonical kernel gained the NNDB int11/int21/int22 tables and
    // vrna_E_internal rules — B12 RNA-STRUCT-001 F3; MIRNA-PAIR-001 Evidence §Duplex energy).
    // Exact Seqeron values (NNDB parameters) are locked; Vienna is asserted within 0.35 kcal/mol.
    // Loop terms of the two generic-loop fixtures (NNDB internal_loop_dg.txt: initiation(8) 2.3, (10) 2.5,
    // (14) 2.9; AU/GU closure 0.7 per pair; first mismatch GA −1.0, other mismatches 0; each value equals
    // ViennaRNA eval_int_loop):
    //   4×4: C·G/U·A 2.3 + 0.7 = 3.0;  G·U/G·U 2.3 + 2×0.7 + GA −1.0 = 2.7 (2.9 before B12: inner-pair
    //        mismatch read as AG −0.8)  ⇒ total 2.55 → 2.35 (Vienna 2.4).
    //   5×5/7×7: 7×7 A·U/C·G 2.9 + 0.7 + GA −1.0 = 2.6 (was 2.8, same mismatch orientation fix);
    //        5×5 G·C/U·G 2.5 + 0.7 = 3.2  ⇒ total 3.00 → 2.80 (Vienna 2.9).
    [TestCase("UGAGGUAGUAGGUUGUAUAGUU", "AACUAUACAACCUACUACCUCA", -35.90, -35.6, TestName = "Vienna_Let7a_PerfectComplement")]
    [TestCase("UGAGGUAGUAGGUUGUAUAGUU", "GGGGGCUACCUCA", -11.67, -11.4, TestName = "Vienna_Let7a_SiteWindow_3x3Loop")]
    [TestCase("UGAGGUAGUAGGUUGUAUAGUU", "CCAACUAUACAGCCUACUACCUCA", -36.09, -35.9, TestName = "Vienna_Let7a_GuAndOverhang")]
    [TestCase("UAGCUUAUCAGACUGAUGUUGA", "UCAACAUCAGUCUGAUAAGCUA", -35.79, -35.7, TestName = "Vienna_MiR21_PerfectComplement")]
    [TestCase("AUGCGAUCC", "GCGGAUGGCAU", -10.50, -10.5, TestName = "Vienna_1x1Loop")]
    [TestCase("UAUCUAUAUAAGCAGGGGA", "UCCCCUGCUAGUAUAGAUA", -26.91, -27.0, TestName = "Vienna_1x1Loop_GU")]
    [TestCase("GUGUGAAUCG", "CGAUUCAUAC", -11.71, -11.8, TestName = "Vienna_GuInHelix")]
    [TestCase("CGCAGUAUGCCAAGACU", "GUCGUCCUUGGCCUACUGCGCA", 2.35, 2.4, TestName = "Vienna_4x4Loop")]
    [TestCase("AGUAUGUCCAACGGCGAGCUUUA", "UGCUAAAGCUCGCCGUUGGACAUACUGUG", 2.80, 2.9, TestName = "Vienna_5x5And7x7Loops")]
    public void AlignMiRnaToTarget_FreeEnergy_MatchesViennaEval(string mirna, string target, double nndb, double vienna)
    {
        var duplex = AlignMiRnaToTarget(mirna, target);

        Assert.Multiple(() =>
        {
            Assert.That(duplex.FreeEnergy, Is.EqualTo(nndb).Within(1e-9), "Turner 2004 (NNDB) duplex energy");
            Assert.That(duplex.FreeEnergy, Is.EqualTo(vienna).Within(0.35), "ViennaRNA 2.8 eval_structure, 0.1-rounding residual");
        });
    }

    #endregion

    #region FoldMiRnaDuplex (RNAduplex model)

    // F1 — Perfect complement: the MFE duplex is the full 22-bp helix and its energy equals the
    // ungapped evaluation (same Turner 2004 structure). ViennaRNA 2.8 RNA.duplexfold: same structure, −35.6.
    [Test]
    public void FoldMiRnaDuplex_PerfectComplement_FullHelixSameEnergyAsAlignment()
    {
        const string let7a = "UGAGGUAGUAGGUUGUAUAGUU";
        string target = GetReverseComplement(let7a);

        var fold = FoldMiRnaDuplex(let7a, target);

        Assert.That(fold, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(fold!.Value.Structure, Is.EqualTo(new string('(', 22) + "&" + new string(')', 22)));
            Assert.That(fold.Value.BasePairs, Has.Count.EqualTo(22));
            Assert.That(fold.Value.BasePairs[0], Is.EqualTo((0, 21)), "miRNA nt 1 pairs the target 3' end");
            Assert.That(fold.Value.FreeEnergy, Is.EqualTo(AlignMiRnaToTarget(let7a, target).FreeEnergy).Within(1e-9));
            Assert.That(fold.Value.FreeEnergy, Is.EqualTo(-35.90).Within(1e-9));
        });
    }

    // F2 — Bulged target: an extra A in the let-7a site. The ungapped register loses the 3' half
    // (AlignMiRnaToTarget −11.11), the duplex fold recovers it with a 1-nt bulge.
    // Ungapped A12 register: loops 2×2 (int22 5′GGUU/3′CACA → +0.5; generic formula before B12 gave +0.8),
    // 5×5 (2.5 + 2×0.7 + GA −1.0 = 2.9), 1×1 (int11 1.9) ⇒ −11.11 (was −10.81); Vienna eval_structure −11.0.
    // 1×2 loop fixture: NNDB turner_2004_int21_dg.txt table 5'GXU/3'CYAA, X = G, Y = A → 1.9 (generic formula
    // before B12 gave 2.4) ⇒ −28.50 (was −28.00); Vienna duplexfold −28.2 (same pairs; int21 also 1.9 —
    // residual = 0.1-rounded stacks/init/AU).
    // ViennaRNA 2.8 RNA.duplexfold("UGAGGUAGUAGGUUGUAUAGUU", <target>) = same pairs, −31.8 (no NNDB bulge degeneracy term).
    [TestCase("AACUAUACAACACUACUACCUCA", "((((((((((((((((((((((&))))))))))).)))))))))))", -32.10, -31.8, -11.11, TestName = "Fold_OneNtBulge_A12")]
    [TestCase("AACUAUACAACCAUACUACCUCA", "((((((((((((((((((((((&)))))))))))).))))))))))", -32.10, -31.8, -7.13, TestName = "Fold_OneNtBulge_A13")]
    [TestCase("AACUAUACAAAACUACUACCUCA", "(((((((((((.((((((((((&))))))))))..)))))))))))", -28.50, -28.2, -11.34, TestName = "Fold_1x2InternalLoop")]
    public void FoldMiRnaDuplex_BulgedSite_RecoversHelixLikeVienna(string target, string structure, double nndb, double vienna, double ungapped)
    {
        const string let7a = "UGAGGUAGUAGGUUGUAUAGUU";

        var fold = FoldMiRnaDuplex(let7a, target);

        Assert.That(fold, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(fold!.Value.Structure, Is.EqualTo(structure), "same base pairs as ViennaRNA duplexfold");
            Assert.That(fold.Value.FreeEnergy, Is.EqualTo(nndb).Within(1e-9), "Turner 2004 (NNDB) energy");
            Assert.That(fold.Value.FreeEnergy, Is.EqualTo(vienna).Within(0.35), "ViennaRNA 2.8 duplexfold");
            Assert.That(AlignMiRnaToTarget(let7a, target).FreeEnergy, Is.EqualTo(ungapped).Within(1e-9),
                "the fixed-register ungapped duplex is far less stable");
        });
    }

    // F3 — Cross-check against ViennaRNA 2.8.0 RNA.duplexfold on generated fixtures: identical base-pair
    // set; energies within Vienna's 0.1-kcal parameter rounding. Coordinates are 0-based inclusive.
    [TestCase("GUGUGAAUCG", "CGAUUCAUAC", "((((((((((&))))))))))", 0, 9, 0, 9, -11.71, -11.8, TestName = "FoldVienna_GuHelix")]
    [TestCase("AUGCGAUCC", "GCGGAUGGCAU", "((((.((((&)))).))))", 0, 8, 2, 10, -10.50, -10.5, TestName = "FoldVienna_1x1Loop")]
    [TestCase("CAUUUCUGGAUGGCCA", "GACUGGCCAUUCAGAUAUG", "(((.((((((((((((&)))))))))))).)))", 0, 15, 3, 18, -21.87, -21.9, TestName = "FoldVienna_1x1LoopNearEnd")]
    [TestCase("AAGCAGGCACGACAGU", "ACUGUCAUGACCGCUUCAC", "((((.((((.((((((&)))))).)).))))))", 0, 15, 0, 15, -14.91, -14.8, TestName = "FoldVienna_Bulges")]
    [TestCase("CGGGCAUAUG", "GCCUAUGCAGG", "((((&))))", 1, 4, 0, 3, -6.75, -6.7, TestName = "FoldVienna_SubDuplex")]
    public void FoldMiRnaDuplex_MatchesViennaDuplexfold(string mirna, string target, string structure,
        int mStart, int mEnd, int tStart, int tEnd, double nndb, double vienna)
    {
        var fold = FoldMiRnaDuplex(mirna, target);

        Assert.That(fold, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(fold!.Value.Structure, Is.EqualTo(structure));
            Assert.That((fold.Value.MiRnaStart, fold.Value.MiRnaEnd, fold.Value.TargetStart, fold.Value.TargetEnd),
                Is.EqualTo((mStart, mEnd, tStart, tEnd)));
            Assert.That(fold.Value.FreeEnergy, Is.EqualTo(nndb).Within(1e-9));
            Assert.That(fold.Value.FreeEnergy, Is.EqualTo(vienna).Within(0.35));
        });
    }

    // F4 — The MFE duplex is never less stable than the fixed-register duplex (whose loops are all ≤ 30 nt
    // for ≤ 23-nt inputs, so it is one of the structures the recursion minimises over).
    [TestCase("UGAGGUAGUAGGUUGUAUAGUU", "GGGGGCUACCUCA")]
    [TestCase("UGAGGUAGUAGGUUGUAUAGUU", "ACUAUACAACGAUACUACCUCAG")]
    [TestCase("CGCAGUAUGCCAAGACU", "GUCGUCCUUGGCCUACUGCGCA")]
    [TestCase("AAAAAA", "UUUU")]
    public void FoldMiRnaDuplex_NeverAboveUngappedDuplex(string mirna, string target)
    {
        var align = AlignMiRnaToTarget(mirna, target);
        var fold = FoldMiRnaDuplex(mirna, target);

        Assert.That(fold, Is.Not.Null);
        Assert.That(fold!.Value.FreeEnergy, Is.LessThanOrEqualTo(align.FreeEnergy + 1e-9));
    }

    // F6 — NNDB GU-pair note: 5'GGUC/3'CUGG replaces its three stacks by −4.12 kcal/mol. Duplex
    // CAGGUCAC·GUGGUCUG (pairs C·G A·U G·C G·U U·G C·G A·U C·G): 4.09 − 2.11 (CA/GU) − 2.08 (AG/UC)
    // − 4.12 (GGUC special) − 2.11 (CA/GU) − 2.24 (AC/UG) = −8.57, in both the ungapped evaluation and
    // the duplex fold. ViennaRNA 2.8 has no special GGUC term (duplexfold −6.1, same structure).
    [Test]
    public void GgucSpecialStack_AppliedByAlignmentAndFold()
    {
        var align = AlignMiRnaToTarget("CAGGUCAC", "GUGGUCUG");
        var fold = FoldMiRnaDuplex("CAGGUCAC", "GUGGUCUG");

        Assert.Multiple(() =>
        {
            Assert.That(align.AlignmentString, Is.EqualTo("|||::|||"));
            Assert.That(align.FreeEnergy, Is.EqualTo(-8.57).Within(1e-9));
            Assert.That(fold!.Value.Structure, Is.EqualTo("((((((((&))))))))"));
            Assert.That(fold.Value.FreeEnergy, Is.EqualTo(-8.57).Within(1e-9));
        });
    }

    // F5 — maxLoopSize 0 allows only stacked pairs (one contiguous helix): the bulged let-7a site then
    // keeps only one helix. Null/empty input and inputs without any possible pair give null; a negative
    // loop bound throws.
    [Test]
    public void FoldMiRnaDuplex_LoopBoundAndDegenerateInputs()
    {
        var helixOnly = FoldMiRnaDuplex("UGAGGUAGUAGGUUGUAUAGUU", "AACUAUACAACACUACUACCUCA", maxLoopSize: 0);

        Assert.Multiple(() =>
        {
            Assert.That(helixOnly, Is.Not.Null);
            Assert.That(helixOnly!.Value.Structure, Does.Not.Contain("."), "no loop of any size is allowed");
            Assert.That(FoldMiRnaDuplex("", "UUUU"), Is.Null);
            Assert.That(FoldMiRnaDuplex("AAAA", null!), Is.Null);
            Assert.That(FoldMiRnaDuplex("AAAA", "AAAA"), Is.Null, "A·A never pairs");
            Assert.Throws<ArgumentOutOfRangeException>(() => FoldMiRnaDuplex("AAAA", "UUUU", -1));
        });
    }

    #endregion

    #region Property

    // C2 — Wobble pairs are a subset of pairable pairs (INV-03): IsWobblePair ⇒ CanPair, over all RNA base pairs.
    [Test]
    public void WobblePairs_AreSubsetOfPairs_Property()
    {
        const string bases = "ACGU";
        foreach (char b1 in bases)
        {
            foreach (char b2 in bases)
            {
                if (IsWobblePair(b1, b2))
                {
                    Assert.That(CanPair(b1, b2), Is.True,
                        $"Every wobble pair must also be a valid pair: {b1}-{b2}");
                }
            }
        }
    }

    #endregion
}
