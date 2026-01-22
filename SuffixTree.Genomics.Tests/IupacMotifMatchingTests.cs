using NUnit.Framework;
using SuffixTree.Genomics;
using System;
using System.Linq;
using System.Threading;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Comprehensive tests for IUPAC degenerate motif matching (PAT-IUPAC-001).
/// Evidence: IUPAC-IUB 1970, NC-IUB 1984, Wikipedia Nucleic acid notation, Bioinformatics.org IUPAC codes
/// </summary>
[TestFixture]
[Category("PAT-IUPAC-001")]
[Description("Test Unit PAT-IUPAC-001: IUPAC Degenerate Motif Matching")]
public class IupacMotifMatchingTests
{
    #region M1-M4: Standard Base Matching (Self-Match)

    [Test]
    [Description("M1: Standard base A matches itself (IUPAC 1970)")]
    public void MatchesIupac_A_MatchesA_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'A'), Is.True);
    }

    [Test]
    [Description("M2: Standard base C matches itself (IUPAC 1970)")]
    public void MatchesIupac_C_MatchesC_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'C'), Is.True);
    }

    [Test]
    [Description("M3: Standard base G matches itself (IUPAC 1970)")]
    public void MatchesIupac_G_MatchesG_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'G'), Is.True);
    }

    [Test]
    [Description("M4: Standard base T matches itself (IUPAC 1970)")]
    public void MatchesIupac_T_MatchesT_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'T'), Is.True);
    }

    [Test]
    [Description("M5: Standard base A does not match T (IUPAC 1970)")]
    public void MatchesIupac_A_DoesNotMatchT_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'T'), Is.False);
    }

    #endregion

    #region M6-M9: N (Any Base) Tests

    [Test]
    [Description("M6: N matches A (IUPAC 1970 - any base)")]
    public void MatchesIupac_N_MatchesA_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'N'), Is.True);
    }

    [Test]
    [Description("M7: N matches C (IUPAC 1970 - any base)")]
    public void MatchesIupac_N_MatchesC_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'N'), Is.True);
    }

    [Test]
    [Description("M8: N matches G (IUPAC 1970 - any base)")]
    public void MatchesIupac_N_MatchesG_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'N'), Is.True);
    }

    [Test]
    [Description("M9: N matches T (IUPAC 1970 - any base)")]
    public void MatchesIupac_N_MatchesT_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'N'), Is.True);
    }

    [Test]
    [Description("INV-8: N matches all four standard bases")]
    public void MatchesIupac_N_MatchesAllBases_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('A', 'N'), Is.True, "N should match A");
            Assert.That(IupacHelper.MatchesIupac('C', 'N'), Is.True, "N should match C");
            Assert.That(IupacHelper.MatchesIupac('G', 'N'), Is.True, "N should match G");
            Assert.That(IupacHelper.MatchesIupac('T', 'N'), Is.True, "N should match T");
        });
    }

    #endregion

    #region M10-M17: R (Purine) and Y (Pyrimidine) Tests

    [Test]
    [Description("M10: R (purine) matches A (Wikipedia)")]
    public void MatchesIupac_R_MatchesA_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'R'), Is.True);
    }

    [Test]
    [Description("M11: R (purine) matches G (Wikipedia)")]
    public void MatchesIupac_R_MatchesG_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'R'), Is.True);
    }

    [Test]
    [Description("M12: R (purine) does NOT match C (Wikipedia)")]
    public void MatchesIupac_R_DoesNotMatchC_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'R'), Is.False);
    }

    [Test]
    [Description("M13: R (purine) does NOT match T (Wikipedia)")]
    public void MatchesIupac_R_DoesNotMatchT_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'R'), Is.False);
    }

    [Test]
    [Description("M14: Y (pyrimidine) matches C (Wikipedia)")]
    public void MatchesIupac_Y_MatchesC_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'Y'), Is.True);
    }

    [Test]
    [Description("M15: Y (pyrimidine) matches T (Wikipedia)")]
    public void MatchesIupac_Y_MatchesT_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'Y'), Is.True);
    }

    [Test]
    [Description("M16: Y (pyrimidine) does NOT match A (Wikipedia)")]
    public void MatchesIupac_Y_DoesNotMatchA_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'Y'), Is.False);
    }

    [Test]
    [Description("M17: Y (pyrimidine) does NOT match G (Wikipedia)")]
    public void MatchesIupac_Y_DoesNotMatchG_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'Y'), Is.False);
    }

    #endregion

    #region M18-M25: S (Strong) and W (Weak) Tests

    [Test]
    [Description("M18: S (strong) matches G (Wikipedia)")]
    public void MatchesIupac_S_MatchesG_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'S'), Is.True);
    }

    [Test]
    [Description("M19: S (strong) matches C (Wikipedia)")]
    public void MatchesIupac_S_MatchesC_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'S'), Is.True);
    }

    [Test]
    [Description("M20: S (strong) does NOT match A (Wikipedia)")]
    public void MatchesIupac_S_DoesNotMatchA_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'S'), Is.False);
    }

    [Test]
    [Description("M21: S (strong) does NOT match T (Wikipedia)")]
    public void MatchesIupac_S_DoesNotMatchT_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'S'), Is.False);
    }

    [Test]
    [Description("M22: W (weak) matches A (Wikipedia)")]
    public void MatchesIupac_W_MatchesA_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'W'), Is.True);
    }

    [Test]
    [Description("M23: W (weak) matches T (Wikipedia)")]
    public void MatchesIupac_W_MatchesT_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'W'), Is.True);
    }

    [Test]
    [Description("M24: W (weak) does NOT match G (Wikipedia)")]
    public void MatchesIupac_W_DoesNotMatchG_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'W'), Is.False);
    }

    [Test]
    [Description("M25: W (weak) does NOT match C (Wikipedia)")]
    public void MatchesIupac_W_DoesNotMatchC_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'W'), Is.False);
    }

    #endregion

    #region M26-M33: K (Keto) and M (Amino) Tests

    [Test]
    [Description("M26: K (keto) matches G (Wikipedia)")]
    public void MatchesIupac_K_MatchesG_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'K'), Is.True);
    }

    [Test]
    [Description("M27: K (keto) matches T (Wikipedia)")]
    public void MatchesIupac_K_MatchesT_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'K'), Is.True);
    }

    [Test]
    [Description("M28: K (keto) does NOT match A (Wikipedia)")]
    public void MatchesIupac_K_DoesNotMatchA_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'K'), Is.False);
    }

    [Test]
    [Description("M29: K (keto) does NOT match C (Wikipedia)")]
    public void MatchesIupac_K_DoesNotMatchC_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'K'), Is.False);
    }

    [Test]
    [Description("M30: M (amino) matches A (Wikipedia)")]
    public void MatchesIupac_M_MatchesA_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'M'), Is.True);
    }

    [Test]
    [Description("M31: M (amino) matches C (Wikipedia)")]
    public void MatchesIupac_M_MatchesC_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'M'), Is.True);
    }

    [Test]
    [Description("M32: M (amino) does NOT match G (Wikipedia)")]
    public void MatchesIupac_M_DoesNotMatchG_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'M'), Is.False);
    }

    [Test]
    [Description("M33: M (amino) does NOT match T (Wikipedia)")]
    public void MatchesIupac_M_DoesNotMatchT_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'M'), Is.False);
    }

    #endregion

    #region M34-M49: Three-Base Ambiguity Codes (B, D, H, V)

    [Test]
    [Description("M34: B (not A) matches C (Wikipedia)")]
    public void MatchesIupac_B_MatchesC_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'B'), Is.True);
    }

    [Test]
    [Description("M35: B (not A) matches G (Wikipedia)")]
    public void MatchesIupac_B_MatchesG_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'B'), Is.True);
    }

    [Test]
    [Description("M36: B (not A) matches T (Wikipedia)")]
    public void MatchesIupac_B_MatchesT_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'B'), Is.True);
    }

    [Test]
    [Description("M37: B (not A) does NOT match A (Wikipedia)")]
    public void MatchesIupac_B_DoesNotMatchA_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'B'), Is.False);
    }

    [Test]
    [Description("M38: D (not C) matches A (Wikipedia)")]
    public void MatchesIupac_D_MatchesA_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'D'), Is.True);
    }

    [Test]
    [Description("M39: D (not C) matches G (Wikipedia)")]
    public void MatchesIupac_D_MatchesG_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'D'), Is.True);
    }

    [Test]
    [Description("M40: D (not C) matches T (Wikipedia)")]
    public void MatchesIupac_D_MatchesT_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'D'), Is.True);
    }

    [Test]
    [Description("M41: D (not C) does NOT match C (Wikipedia)")]
    public void MatchesIupac_D_DoesNotMatchC_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'D'), Is.False);
    }

    [Test]
    [Description("M42: H (not G) matches A (Wikipedia)")]
    public void MatchesIupac_H_MatchesA_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'H'), Is.True);
    }

    [Test]
    [Description("M43: H (not G) matches C (Wikipedia)")]
    public void MatchesIupac_H_MatchesC_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'H'), Is.True);
    }

    [Test]
    [Description("M44: H (not G) matches T (Wikipedia)")]
    public void MatchesIupac_H_MatchesT_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'H'), Is.True);
    }

    [Test]
    [Description("M45: H (not G) does NOT match G (Wikipedia)")]
    public void MatchesIupac_H_DoesNotMatchG_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'H'), Is.False);
    }

    [Test]
    [Description("M46: V (not T) matches A (Wikipedia)")]
    public void MatchesIupac_V_MatchesA_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('A', 'V'), Is.True);
    }

    [Test]
    [Description("M47: V (not T) matches C (Wikipedia)")]
    public void MatchesIupac_V_MatchesC_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('C', 'V'), Is.True);
    }

    [Test]
    [Description("M48: V (not T) matches G (Wikipedia)")]
    public void MatchesIupac_V_MatchesG_ReturnsTrue()
    {
        Assert.That(IupacHelper.MatchesIupac('G', 'V'), Is.True);
    }

    [Test]
    [Description("M49: V (not T) does NOT match T (Wikipedia)")]
    public void MatchesIupac_V_DoesNotMatchT_ReturnsFalse()
    {
        Assert.That(IupacHelper.MatchesIupac('T', 'V'), Is.False);
    }

    #endregion

    #region Comprehensive IUPAC Code Invariant Tests

    [Test]
    [Description("INV-3: R (purine) matches exactly A and G")]
    public void MatchesIupac_R_MatchesExactlyAG_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('A', 'R'), Is.True, "R should match A (purine)");
            Assert.That(IupacHelper.MatchesIupac('G', 'R'), Is.True, "R should match G (purine)");
            Assert.That(IupacHelper.MatchesIupac('C', 'R'), Is.False, "R should NOT match C");
            Assert.That(IupacHelper.MatchesIupac('T', 'R'), Is.False, "R should NOT match T");
        });
    }

    [Test]
    [Description("INV-3: Y (pyrimidine) matches exactly C and T")]
    public void MatchesIupac_Y_MatchesExactlyCT_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('C', 'Y'), Is.True, "Y should match C (pyrimidine)");
            Assert.That(IupacHelper.MatchesIupac('T', 'Y'), Is.True, "Y should match T (pyrimidine)");
            Assert.That(IupacHelper.MatchesIupac('A', 'Y'), Is.False, "Y should NOT match A");
            Assert.That(IupacHelper.MatchesIupac('G', 'Y'), Is.False, "Y should NOT match G");
        });
    }

    [Test]
    [Description("INV-3: S (strong) matches exactly G and C")]
    public void MatchesIupac_S_MatchesExactlyGC_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('G', 'S'), Is.True, "S should match G (strong)");
            Assert.That(IupacHelper.MatchesIupac('C', 'S'), Is.True, "S should match C (strong)");
            Assert.That(IupacHelper.MatchesIupac('A', 'S'), Is.False, "S should NOT match A");
            Assert.That(IupacHelper.MatchesIupac('T', 'S'), Is.False, "S should NOT match T");
        });
    }

    [Test]
    [Description("INV-3: W (weak) matches exactly A and T")]
    public void MatchesIupac_W_MatchesExactlyAT_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('A', 'W'), Is.True, "W should match A (weak)");
            Assert.That(IupacHelper.MatchesIupac('T', 'W'), Is.True, "W should match T (weak)");
            Assert.That(IupacHelper.MatchesIupac('G', 'W'), Is.False, "W should NOT match G");
            Assert.That(IupacHelper.MatchesIupac('C', 'W'), Is.False, "W should NOT match C");
        });
    }

    [Test]
    [Description("INV-3: K (keto) matches exactly G and T")]
    public void MatchesIupac_K_MatchesExactlyGT_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('G', 'K'), Is.True, "K should match G (keto)");
            Assert.That(IupacHelper.MatchesIupac('T', 'K'), Is.True, "K should match T (keto)");
            Assert.That(IupacHelper.MatchesIupac('A', 'K'), Is.False, "K should NOT match A");
            Assert.That(IupacHelper.MatchesIupac('C', 'K'), Is.False, "K should NOT match C");
        });
    }

    [Test]
    [Description("INV-3: M (amino) matches exactly A and C")]
    public void MatchesIupac_M_MatchesExactlyAC_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('A', 'M'), Is.True, "M should match A (amino)");
            Assert.That(IupacHelper.MatchesIupac('C', 'M'), Is.True, "M should match C (amino)");
            Assert.That(IupacHelper.MatchesIupac('G', 'M'), Is.False, "M should NOT match G");
            Assert.That(IupacHelper.MatchesIupac('T', 'M'), Is.False, "M should NOT match T");
        });
    }

    [Test]
    [Description("INV-3/4: B (not A) matches exactly C, G, T and excludes A")]
    public void MatchesIupac_B_MatchesExactlyCGT_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('C', 'B'), Is.True, "B should match C");
            Assert.That(IupacHelper.MatchesIupac('G', 'B'), Is.True, "B should match G");
            Assert.That(IupacHelper.MatchesIupac('T', 'B'), Is.True, "B should match T");
            Assert.That(IupacHelper.MatchesIupac('A', 'B'), Is.False, "B should NOT match A (not A)");
        });
    }

    [Test]
    [Description("INV-3/4: D (not C) matches exactly A, G, T and excludes C")]
    public void MatchesIupac_D_MatchesExactlyAGT_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('A', 'D'), Is.True, "D should match A");
            Assert.That(IupacHelper.MatchesIupac('G', 'D'), Is.True, "D should match G");
            Assert.That(IupacHelper.MatchesIupac('T', 'D'), Is.True, "D should match T");
            Assert.That(IupacHelper.MatchesIupac('C', 'D'), Is.False, "D should NOT match C (not C)");
        });
    }

    [Test]
    [Description("INV-3/4: H (not G) matches exactly A, C, T and excludes G")]
    public void MatchesIupac_H_MatchesExactlyACT_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('A', 'H'), Is.True, "H should match A");
            Assert.That(IupacHelper.MatchesIupac('C', 'H'), Is.True, "H should match C");
            Assert.That(IupacHelper.MatchesIupac('T', 'H'), Is.True, "H should match T");
            Assert.That(IupacHelper.MatchesIupac('G', 'H'), Is.False, "H should NOT match G (not G)");
        });
    }

    [Test]
    [Description("INV-3/4: V (not T) matches exactly A, C, G and excludes T")]
    public void MatchesIupac_V_MatchesExactlyACG_Invariant()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('A', 'V'), Is.True, "V should match A");
            Assert.That(IupacHelper.MatchesIupac('C', 'V'), Is.True, "V should match C");
            Assert.That(IupacHelper.MatchesIupac('G', 'V'), Is.True, "V should match G");
            Assert.That(IupacHelper.MatchesIupac('T', 'V'), Is.False, "V should NOT match T (not T)");
        });
    }

    #endregion

    #region M50-M54: FindDegenerateMotif Pattern Matching Tests

    [Test]
    [Description("M50: Purine R in pattern matches A and G in sequence (Wikipedia)")]
    public void FindDegenerateMotif_PurineR_MatchesAandG()
    {
        var sequence = new DnaSequence("ATGC");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "R").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(2), "R should match A at 0 and G at 2");
            Assert.That(matches.Select(m => m.Position), Is.EquivalentTo(new[] { 0, 2 }));
        });
    }

    [Test]
    [Description("M51: Pyrimidine Y in pattern matches C and T in sequence (Wikipedia)")]
    public void FindDegenerateMotif_PyrimidineY_MatchesCandT()
    {
        var sequence = new DnaSequence("ATGC");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "Y").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(2), "Y should match T at 1 and C at 3");
            Assert.That(matches.Select(m => m.Position), Is.EquivalentTo(new[] { 1, 3 }));
        });
    }

    [Test]
    [Description("M52: Any N in pattern matches all bases (Wikipedia)")]
    public void FindDegenerateMotif_AnyN_MatchesAllBases()
    {
        var sequence = new DnaSequence("ACGT");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "N").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(4), "N should match all 4 positions");
            Assert.That(matches.Select(m => m.Position), Is.EquivalentTo(new[] { 0, 1, 2, 3 }));
        });
    }

    [Test]
    [Description("M53: Mixed pattern RTG matches ATG and GTG (Wikipedia)")]
    public void FindDegenerateMotif_MixedPatternRTG_MatchesATGandGTG()
    {
        var sequence = new DnaSequence("ATGCGTGC");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "RTG").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(2), "RTG should match ATG at 0 and GTG at 4");
            Assert.That(matches.Select(m => m.Position), Is.EquivalentTo(new[] { 0, 4 }));
            Assert.That(matches[0].MatchedSequence, Is.EqualTo("ATG"));
            Assert.That(matches[1].MatchedSequence, Is.EqualTo("GTG"));
        });
    }

    [Test]
    [Description("M54: E-box pattern CANNTG matches biological motif")]
    public void FindDegenerateMotif_EBoxCANNTG_MatchesBiologicalMotif()
    {
        var sequence = new DnaSequence("CAGCTG");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "CANNTG").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(1), "CANNTG should match CAGCTG");
            Assert.That(matches[0].Position, Is.EqualTo(0));
            Assert.That(matches[0].MatchedSequence, Is.EqualTo("CAGCTG"));
        });
    }

    #endregion

    #region M55-M59: Edge Cases for FindDegenerateMotif

    [Test]
    [Description("M55: No match returns empty collection")]
    public void FindDegenerateMotif_NoMatch_ReturnsEmpty()
    {
        var sequence = new DnaSequence("AAAA");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "GGG").ToList();

        Assert.That(matches, Is.Empty);
    }

    [Test]
    [Description("M56: Empty pattern returns empty collection")]
    public void FindDegenerateMotif_EmptyPattern_ReturnsEmpty()
    {
        var sequence = new DnaSequence("ACGT");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "").ToList();

        Assert.That(matches, Is.Empty);
    }

    [Test]
    [Description("M58: Pattern longer than sequence returns empty")]
    public void FindDegenerateMotif_PatternLongerThanSequence_ReturnsEmpty()
    {
        var sequence = new DnaSequence("AC");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "ACGT").ToList();

        Assert.That(matches, Is.Empty);
    }

    [Test]
    [Description("M59: Null sequence throws ArgumentNullException")]
    public void FindDegenerateMotif_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.FindDegenerateMotif(null!, "ATG").ToList());
    }

    #endregion

    #region M60-M62: Result Verification Tests

    [Test]
    [Description("M60: Case insensitive pattern matching")]
    public void FindDegenerateMotif_CaseInsensitive_MatchesPattern()
    {
        var sequence = new DnaSequence("ATGC");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "atgc").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(1));
            Assert.That(matches[0].Position, Is.EqualTo(0));
        });
    }

    [Test]
    [Description("M61: Result contains correct matched sequence")]
    public void FindDegenerateMotif_ResultContainsMatchedSequence()
    {
        var sequence = new DnaSequence("CAGCTG");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "CANNTG").ToList();

        Assert.That(matches[0].MatchedSequence, Is.EqualTo("CAGCTG"));
    }

    [Test]
    [Description("M62/INV-5: Result positions are within valid range")]
    public void FindDegenerateMotif_ResultPositionsValid_Invariant()
    {
        var sequence = new DnaSequence("ACGTACGTACGT");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "NNN").ToList();

        int maxValidPosition = sequence.Length - 3;
        Assert.That(matches.All(m => m.Position >= 0 && m.Position <= maxValidPosition), Is.True,
            $"All positions should be in range [0, {maxValidPosition}]");
    }

    #endregion

    #region S1-S11: SHOULD Tests (Important Edge Cases)

    [Test]
    [Description("S1: W (weak) pattern finds AT alternation")]
    public void FindDegenerateMotif_WeakW_FindsATAlternation()
    {
        var sequence = new DnaSequence("ATATAT");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "WWW").ToList();

        Assert.That(matches.Count, Is.EqualTo(4), "WWW should match 4 overlapping positions");
    }

    [Test]
    [Description("S2: S (strong) pattern finds GC regions")]
    public void FindDegenerateMotif_StrongS_FindsGCRegions()
    {
        var sequence = new DnaSequence("GCGCGC");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "SSS").ToList();

        Assert.That(matches.Count, Is.EqualTo(4), "SSS should match 4 overlapping positions");
    }

    [Test]
    [Description("S3: K (keto) pattern finds GT alternation")]
    public void FindDegenerateMotif_KetoK_FindsGTAlternation()
    {
        var sequence = new DnaSequence("GTGTGT");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "KK").ToList();

        Assert.That(matches.Count, Is.EqualTo(5), "KK should match 5 overlapping positions");
    }

    [Test]
    [Description("S4: M (amino) pattern finds AC alternation")]
    public void FindDegenerateMotif_AminoM_FindsACAlternation()
    {
        var sequence = new DnaSequence("ACACAC");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "MM").ToList();

        Assert.That(matches.Count, Is.EqualTo(5), "MM should match 5 overlapping positions");
    }

    [Test]
    [Description("S5: Pattern at end of sequence is found")]
    public void FindDegenerateMotif_PatternAtEnd_Found()
    {
        var sequence = new DnaSequence("GCGATG");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "ATG").ToList();

        Assert.That(matches.Select(m => m.Position), Does.Contain(3));
    }

    [Test]
    [Description("S8: All B (not A) positions match CGT sequence")]
    public void FindDegenerateMotif_AllB_MatchesCGT()
    {
        var sequence = new DnaSequence("CGT");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "BBB").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(1));
            Assert.That(matches[0].MatchedSequence, Is.EqualTo("CGT"));
        });
    }

    [Test]
    [Description("S9: All D (not C) positions match AGT sequence")]
    public void FindDegenerateMotif_AllD_MatchesAGT()
    {
        var sequence = new DnaSequence("AGT");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "DDD").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(1));
            Assert.That(matches[0].MatchedSequence, Is.EqualTo("AGT"));
        });
    }

    [Test]
    [Description("S10: All H (not G) positions match ACT sequence")]
    public void FindDegenerateMotif_AllH_MatchesACT()
    {
        var sequence = new DnaSequence("ACT");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "HHH").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(1));
            Assert.That(matches[0].MatchedSequence, Is.EqualTo("ACT"));
        });
    }

    [Test]
    [Description("S11: All V (not T) positions match ACG sequence")]
    public void FindDegenerateMotif_AllV_MatchesACG()
    {
        var sequence = new DnaSequence("ACG");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "VVV").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(1));
            Assert.That(matches[0].MatchedSequence, Is.EqualTo("ACG"));
        });
    }

    [Test]
    [Description("S12: Cancellation token is respected")]
    public void FindDegenerateMotif_WithCancellationToken_CompletesNormally()
    {
        var sequence = new DnaSequence("ACGTACGTACGT");
        using var cts = new CancellationTokenSource();

        var matches = MotifFinder.FindDegenerateMotif(sequence, "ACGT", cts.Token).ToList();

        Assert.That(matches, Is.Not.Empty);
    }

    #endregion

    #region Unknown Character Handling (ASSUMPTION)

    [Test]
    [Description("S7/A2: Unknown character in pattern matches exactly (ASSUMPTION)")]
    public void MatchesIupac_UnknownCode_MatchesExactly()
    {
        // According to implementation, unknown codes fall back to exact match
        Assert.Multiple(() =>
        {
            Assert.That(IupacHelper.MatchesIupac('X', 'X'), Is.True, "Unknown code should match itself");
            Assert.That(IupacHelper.MatchesIupac('A', 'X'), Is.False, "Unknown code should not match other chars");
        });
    }

    #endregion
}
