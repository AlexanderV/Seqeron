// MIRNA-PAIR-001 — MiRNA-Target Pairing Analysis
// Evidence: docs/Evidence/MIRNA-PAIR-001-Evidence.md
// TestSpec: tests/TestSpecs/MIRNA-PAIR-001.md
// Source: Agarwal V et al. (2015) eLife 4:e05005 (A-U/C-G); Crick FHC (1966) J Mol Biol 19(2):548-555 (G-U wobble);
//         Lewis BP et al. (2005) Cell 120(1):15-20 (reverse-complement seed); NNDB Turner 2004 (stacking sign).

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

    // S1 — Fully Watson-Crick paired duplex has non-positive (stabilising) free energy.
    // Evidence: NNDB Turner 2004 nearest-neighbor stacking energies are negative for paired stacks (sign only — see ASSUMPTION 1).
    [Test]
    public void AlignMiRnaToTarget_FullyPairedDuplex_FreeEnergyNonPositive()
    {
        // GCGCGCGC is a palindrome; its reverse complement self-pairs perfectly (8 WC pairs).
        string mirna = "GCGCGCGC";
        string target = GetReverseComplement(mirna);

        var duplex = AlignMiRnaToTarget(mirna, target);

        Assert.Multiple(() =>
        {
            Assert.That(duplex.Matches, Is.EqualTo(8), "Fully Watson-Crick paired duplex (8 pairs)");
            Assert.That(duplex.FreeEnergy, Is.LessThanOrEqualTo(0.0),
                "Sum of negative Turner stacking energies must be stabilising (ΔG ≤ 0)");
        });
    }

    // S2 — All-mismatch duplex has no paired stacks → free energy is not stabilising (≥ 0).
    // Evidence: with no consecutive paired positions the stacking sum is 0 (INV-06).
    [Test]
    public void AlignMiRnaToTarget_AllMismatchDuplex_FreeEnergyNonNegative()
    {
        var duplex = AlignMiRnaToTarget("AAAA", "AAAA");

        Assert.Multiple(() =>
        {
            Assert.That(duplex.Mismatches, Is.EqualTo(4), "No pairs form");
            Assert.That(duplex.FreeEnergy, Is.GreaterThanOrEqualTo(0.0),
                "No paired stacks → stacking sum 0; duplex is not stabilising");
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
