// PRIMER-TM-001 — DNA hairpin folding + secondary-structure (hairpin) Tm (opt-in)
// Evidence: docs/Evidence/PRIMER-TM-001-HAIRPIN-Evidence.md
// TestSpec: tests/TestSpecs/PRIMER-TM-001-HAIRPIN.md
// Sources:
//   SantaLucia J, Hicks D (2004). Annu Rev Biophys 33:415-440 — Table 1 NN stacks,
//     Table 4 hairpin-loop ΔG°37 by size, Eqs 7-11 (hairpin model + unimolecular Tm).
//   SantaLucia J (1998). PNAS 95(4):1460-65 — unified NN ΔH°/ΔS° (stem stacks reused).
//   Vallone PM, Benight AS (1999). Biochemistry — hairpin Tm concentration-independence.

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the opt-in DNA hairpin folder + unimolecular hairpin Tm added under
/// PRIMER-TM-001. Expected ΔH°/ΔS°/ΔG°37/Tm are hand-derived from the SantaLucia
/// (1998) Table 1 NN stem stacks and the SantaLucia &amp; Hicks (2004) Table 4
/// hairpin-loop-initiation increments, independent of the implementation.
/// </summary>
[TestFixture]
public class PrimerDesigner_HairpinTm_Tests
{
    private const double Tol = 1e-9;

    // Canonical hairpin: 5'-GGGC + TTTT(loop) + GCCC-3', a 4-bp stem GGGC/GCCC closing a 4-nt loop.
    // Stem NN steps over the 5' arm GGGC: GG(-8.0,-19.9) GG(-8.0,-19.9) GC(-9.8,-24.4).
    //   Stem ΔH° = -25.8 ; Stem ΔS° = -64.2  (NO bimolecular init term for a unimolecular hairpin).
    // Loop of 4: ΔG°37 = 3.5 (Table 4) ; loop ΔH° = 0 ; loop ΔS° = -3.5*1000/310.15 = -11.28486216346929.
    //   Total ΔH° = -25.8 ; Total ΔS° = -75.48486216346927.
    //   Total ΔG°37 = -25.8 - 310.15*(-75.48486216346927)/1000 = -2.3883700000000054.
    //   Hairpin Tm = (-25.8*1000)/(-75.48486216346927) - 273.15 = 68.6403836682880 °C.
    private const string CanonicalHairpin = "GGGCTTTTGCCC";
    private const double CanonicalDeltaH = -25.8;
    private const double CanonicalDeltaS = -75.48486216346927;
    private const double CanonicalDeltaG37 = -2.3883700000000054;
    private const double CanonicalTm = 68.64038366828805;

    #region FindMostStableHairpin — thermodynamics (M1, M2, M4, M7)

    // M1 — exact ΔH°/ΔS° of the hand-derived canonical hairpin.
    [Test]
    public void FindMostStableHairpin_CanonicalHairpin_MatchesHandDerivedEnthalpyEntropy()
    {
        var hp = PrimerDesigner.FindMostStableHairpin(CanonicalHairpin);

        Assert.That(hp, Is.Not.Null, "GGGCTTTTGCCC forms a 4-bp stem closing a 4-nt loop.");
        Assert.Multiple(() =>
        {
            Assert.That(hp!.Value.DeltaH, Is.EqualTo(CanonicalDeltaH).Within(Tol),
                "ΔH° = Σ stem NN stacks GG+GG+GC (-25.8); loop ΔH° = 0 (Table 4); no bimolecular init.");
            Assert.That(hp.Value.DeltaS, Is.EqualTo(CanonicalDeltaS).Within(Tol),
                "ΔS° = stem(-64.2) + loop(-3.5*1000/310.15) per SantaLucia&Hicks 2004 Table 4 footnote a.");
        });
    }

    // M2 — exact ΔG°37 of the canonical hairpin (a wrong stem or loop table fails this).
    [Test]
    public void FindMostStableHairpin_CanonicalHairpin_MatchesHandDerivedDeltaG37()
    {
        var hp = PrimerDesigner.FindMostStableHairpin(CanonicalHairpin);

        Assert.That(hp!.Value.DeltaG37, Is.EqualTo(CanonicalDeltaG37).Within(Tol),
            "ΔG°37 = ΔH° - 310.15*ΔS°/1000 (= stem stacks + loop ΔG°37 3.5) per SantaLucia&Hicks 2004 Eq.10.");
    }

    // M4 — the folder FINDS the correct hairpin (full 4-bp stem, 4-nt loop, spanning the oligo),
    // not a worse partial structure.
    [Test]
    public void FindMostStableHairpin_CanonicalHairpin_FindsFullStemAndLoop()
    {
        var hp = PrimerDesigner.FindMostStableHairpin(CanonicalHairpin);

        Assert.Multiple(() =>
        {
            Assert.That(hp!.Value.StemLength, Is.EqualTo(4),
                "The most stable hairpin uses the full 4-bp stem GGGC/GCCC.");
            Assert.That(hp.Value.LoopSize, Is.EqualTo(4),
                "The 4-nt TTTT loop is closed by the stem.");
            Assert.That(hp.Value.StemStart, Is.EqualTo(0),
                "Stem 5' arm starts at index 0.");
            Assert.That(hp.Value.StemEnd, Is.EqualTo(11),
                "Stem 3' arm ends at the last index (11).");
        });
    }

    // M7 — the loop ΔS° increment follows ΔS° = -ΔG°37*1000/310.15 (a sign/T error must fail).
    // Total ΔS° minus the stem ΔS° (-64.2) must equal the loop increment.
    [Test]
    public void FindMostStableHairpin_LoopEntropy_FollowsTable4Rule()
    {
        const double stemDeltaS = -64.2;                       // GG+GG+GC entropy sum
        const double expectedLoopDeltaS = -3.5 * 1000.0 / 310.15;

        var hp = PrimerDesigner.FindMostStableHairpin(CanonicalHairpin);
        double loopDeltaS = hp!.Value.DeltaS - stemDeltaS;

        Assert.That(loopDeltaS, Is.EqualTo(expectedLoopDeltaS).Within(Tol),
            "Loop ΔS° = -ΔG°37*1000/310.15 (-11.2849); loop is destabilising so ΔS° < 0 (Table 4 footnote a).");
    }

    #endregion

    #region CalculateHairpinMeltingTemperature — unimolecular Tm (M3, M8)

    // M3 — exact unimolecular hairpin Tm (Eq.11).
    [Test]
    public void CalculateHairpinMeltingTemperature_CanonicalHairpin_MatchesHandDerivedTm()
    {
        double tm = PrimerDesigner.CalculateHairpinMeltingTemperature(CanonicalHairpin);

        Assert.That(tm, Is.EqualTo(CanonicalTm).Within(1e-7),
            "Tm = ΔH°*1000/ΔS° - 273.15 = 68.6404 °C (SantaLucia&Hicks 2004 Eq.11).");
    }

    // M8 — concentration independence: the hairpin Tm equals ΔH°*1000/ΔS° - 273.15 with NO
    // R*ln(C_T/x) strand-concentration term. We recompute that bare formula from the returned
    // ΔH°/ΔS° and require equality (the method exposes no concentration parameter at all).
    [Test]
    public void CalculateHairpinMeltingTemperature_IsConcentrationIndependent()
    {
        var hp = PrimerDesigner.FindMostStableHairpin(CanonicalHairpin);
        double expected = (hp!.Value.DeltaH * 1000.0) / hp.Value.DeltaS - 273.15;

        double tm = PrimerDesigner.CalculateHairpinMeltingTemperature(CanonicalHairpin);

        Assert.That(tm, Is.EqualTo(expected).Within(Tol),
            "Intramolecular hairpin Tm is concentration-independent: ΔH°*1000/ΔS° - 273.15, no C_T term.");
    }

    #endregion

    #region Folding correctness — non-hairpin & loop size (M5, M6, M9, M10)

    // M5 — poly-A cannot form a hairpin (no Watson-Crick stem) → null.
    [Test]
    public void FindMostStableHairpin_PolyA_ReturnsNull()
    {
        var hp = PrimerDesigner.FindMostStableHairpin("AAAAAAAAAAAA");

        Assert.That(hp, Is.Null, "A homopolymer has no complementary stem, so no hairpin exists.");
    }

    // M6 — poly-A hairpin Tm is NaN (no structure).
    [Test]
    public void CalculateHairpinMeltingTemperature_PolyA_ReturnsNaN()
    {
        double tm = PrimerDesigner.CalculateHairpinMeltingTemperature("AAAAAAAAAAAA");

        Assert.That(double.IsNaN(tm), Is.True, "No hairpin → NaN Tm.");
    }

    // M9 — a 5-nt loop uses the Table 4 size-5 increment (ΔG°37 = 3.3), not the size-4 value.
    // Stem GGGC (same -25.8/-64.2); loop of 5 ΔG°37 = 3.3 → loop ΔS° = -3.3*1000/310.15.
    //   Total ΔS° = -64.2 - 10.639787... = -74.84001289698533.
    //   Total ΔG°37 = -25.8 - 310.15*(-74.84001289698533)/1000 = -2.5883700000000054.
    [Test]
    public void FindMostStableHairpin_FiveNtLoop_UsesTable4Size5Increment()
    {
        var hp = PrimerDesigner.FindMostStableHairpin("GGGCAAAAAGCCC");

        Assert.Multiple(() =>
        {
            Assert.That(hp!.Value.LoopSize, Is.EqualTo(5),
                "The 5-nt AAAAA loop is closed by the 4-bp GGGC/GCCC stem.");
            Assert.That(hp.Value.DeltaG37, Is.EqualTo(-2.5883700000000054).Within(Tol),
                "Loop-of-5 ΔG°37 = 3.3 (Table 4), distinct from the loop-of-4 value 3.5.");
        });
    }

    // M10 — a loop shorter than 3 nt is never returned (sterically prohibited). GCGC could only
    // close a 0-nt loop, so no valid hairpin exists → null.
    [Test]
    public void FindMostStableHairpin_NoLoopOfAtLeastThree_ReturnsNull()
    {
        var hp = PrimerDesigner.FindMostStableHairpin("GCGC");

        Assert.That(hp, Is.Null,
            "GCGC can only close a 0-nt loop; loops < 3 nt are sterically prohibited (SantaLucia&Hicks 2004).");
    }

    #endregion

    #region Invalid input & guards (S1, S2, S3)

    // S1 — null / empty input → null and NaN.
    [Test]
    public void FindMostStableHairpin_NullOrEmpty_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindMostStableHairpin(null!), Is.Null, "null input → null.");
            Assert.That(PrimerDesigner.FindMostStableHairpin(""), Is.Null, "empty input → null.");
            Assert.That(double.IsNaN(PrimerDesigner.CalculateHairpinMeltingTemperature("")), Is.True,
                "empty input → NaN Tm.");
        });
    }

    // S2 — a non-ACGT base makes the structure not computable → null.
    [Test]
    public void FindMostStableHairpin_NonAcgtBase_ReturnsNull()
    {
        var hp = PrimerDesigner.FindMostStableHairpin("GGGCNNNNGCCC");

        Assert.That(hp, Is.Null, "Degenerate/non-ACGT bases are rejected (strict alphabet).");
    }

    // S3 — minStemLength < 2 (no NN stack possible) → null.
    [Test]
    public void FindMostStableHairpin_MinStemLengthBelowTwo_ReturnsNull()
    {
        var hp = PrimerDesigner.FindMostStableHairpin(CanonicalHairpin, minStemLength: 1);

        Assert.That(hp, Is.Null, "A stem needs at least one NN stack (≥ 2 bp); minStemLength < 2 → null.");
    }

    #endregion

    #region Extrapolation & opt-in bonus (S4, C1, C2)

    // S4 — Jacobson-Stockmayer: a loop larger than the largest tabulated size (30) is more
    // destabilising (larger loop ΔG°37) than a smaller loop with the same stem, so the longer-loop
    // hairpin has a HIGHER ΔG°37. We compare two oligos with the same 4-bp stem but loops of 6 vs 35.
    [Test]
    public void FindMostStableHairpin_LongerLoop_IsLessStable()
    {
        var smallLoop = PrimerDesigner.FindMostStableHairpin("GGGC" + new string('A', 6) + "GCCC");
        var bigLoop = PrimerDesigner.FindMostStableHairpin("GGGC" + new string('A', 35) + "GCCC");

        Assert.Multiple(() =>
        {
            Assert.That(smallLoop!.Value.LoopSize, Is.EqualTo(6), "6-nt loop closed by the GGGC stem.");
            Assert.That(bigLoop!.Value.LoopSize, Is.EqualTo(35), "35-nt loop (beyond Table 4) closed by the stem.");
            Assert.That(bigLoop.Value.DeltaG37, Is.GreaterThan(smallLoop.Value.DeltaG37),
                "Jacobson-Stockmayer: a longer loop is more destabilising (higher ΔG°37), Eq.7.");
        });
    }

    // C1 — the opt-in terminal-mismatch / special-loop bonus increment shifts ΔG°37 by exactly that
    // amount (and shifts ΔS° by -bonus*1000/310.15). A +1.0 bonus must raise ΔG°37 by 1.0.
    [Test]
    public void FindMostStableHairpin_LoopBonus_ShiftsDeltaG37ByBonus()
    {
        var baseline = PrimerDesigner.FindMostStableHairpin(CanonicalHairpin);
        var withBonus = PrimerDesigner.FindMostStableHairpin(CanonicalHairpin, loopBonusDeltaG37: 1.0);

        Assert.That(withBonus!.Value.DeltaG37 - baseline!.Value.DeltaG37, Is.EqualTo(1.0).Within(Tol),
            "A +1.0 kcal/mol caller-supplied loop bonus raises the hairpin ΔG°37 by exactly 1.0.");
    }

    // C2 — when two stems compete, the folder returns the more stable (lower ΔG°37) one.
    // CCCCAAAAGGGG: the only complementary stem is CCCC/GGGG (4 bp) closing AAAA (4 nt). A weaker
    // partial structure would have a higher ΔG°37; the MFE result is the full stem.
    [Test]
    public void FindMostStableHairpin_CompetingStems_PicksMostStable()
    {
        var hp = PrimerDesigner.FindMostStableHairpin("CCCCAAAAGGGG");

        Assert.Multiple(() =>
        {
            Assert.That(hp!.Value.StemLength, Is.EqualTo(4),
                "The MFE hairpin uses the full 4-bp CCCC/GGGG stem, not a shorter sub-stem.");
            Assert.That(hp.Value.LoopSize, Is.EqualTo(4), "AAAA forms the 4-nt loop.");
            Assert.That(hp.Value.DeltaG37, Is.LessThan(0.0),
                "The full-stem hairpin is stable (negative ΔG°37).");
        });
    }

    #endregion

    #region Varying stem length, minStemLength selectivity, palindrome, very long (E1–E4)

    // E1 — a DIFFERENT stem length (3 bp, not 4). GCCAAAGGC = stem GCC/GGC (3 bp) closing AAA (3 nt).
    // 5'-arm GCC NN steps: GC(-9.8,-24.4) + CC(-8.0,-19.9) → stem ΔH°=-17.8, ΔS°=-44.3.
    // Loop of 3 ΔG°37 = 3.5 (Table 4) → loop ΔS° = -3.5*1000/310.15 = -11.28486216346929.
    //   Total ΔS° = -55.58486216346929; ΔG°37 = -17.8 - 310.15*(-55.58486216346929)/1000 = -0.5603550000000013;
    //   Tm = -17.8*1000/-55.58486216346929 - 273.15 = 47.08107204353689 °C.
    // Hand-derived from SantaLucia&Hicks 2004 Table 1 (GC, CC stacks) + Table 4 (loop-3 = 3.5).
    [Test]
    public void FindMostStableHairpin_ThreeBasePairStem_MatchesHandDerived()
    {
        var hp = PrimerDesigner.FindMostStableHairpin("GCCAAAGGC");

        Assert.That(hp, Is.Not.Null, "GCCAAAGGC forms a 3-bp stem closing a 3-nt loop.");
        var h = hp!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(h.StemLength, Is.EqualTo(3), "Full 3-bp GCC/GGC stem.");
            Assert.That(h.LoopSize, Is.EqualTo(3), "3-nt AAA loop.");
            Assert.That(h.DeltaH, Is.EqualTo(-17.8).Within(Tol), "Stem ΔH° = GC(-9.8) + CC(-8.0).");
            Assert.That(h.DeltaS, Is.EqualTo(-55.58486216346929).Within(Tol),
                "Stem ΔS° (-44.3) + loop ΔS° (-3.5*1000/310.15).");
            Assert.That(h.DeltaG37, Is.EqualTo(-0.5603550000000013).Within(Tol), "ΔG°37 hand-derived.");
            Assert.That(PrimerDesigner.CalculateHairpinMeltingTemperature("GCCAAAGGC"),
                Is.EqualTo(47.08107204353689).Within(1e-7), "Unimolecular Tm = ΔH°*1000/ΔS° - 273.15.");
        });
    }

    // E2 — minStemLength is selective: GCCAAAGGC has only a 3-bp stem, so minStemLength=4 → null,
    // while the default (2) returns the 3-bp hairpin. Locks the stem-length gate beyond the <2 guard.
    [Test]
    public void FindMostStableHairpin_MinStemLength_FiltersStemsBelowThreshold()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindMostStableHairpin("GCCAAAGGC", minStemLength: 4), Is.Null,
                "No 4-bp stem exists; minStemLength=4 → null.");
            Assert.That(PrimerDesigner.FindMostStableHairpin("GCCAAAGGC", minStemLength: 3)!.Value.StemLength,
                Is.EqualTo(3), "minStemLength=3 admits the 3-bp stem.");
        });
    }

    // E3 — a self-complementary palindrome with NO interior loop forms no hairpin. GGGGCCCC pairs
    // fully (4-bp stem) but closes a 0-nt loop, which is sterically prohibited (< 3) → null.
    [Test]
    public void FindMostStableHairpin_PalindromeNoLoop_ReturnsNull()
    {
        var hp = PrimerDesigner.FindMostStableHairpin("GGGGCCCC");

        Assert.That(hp, Is.Null,
            "A perfect palindrome leaves no ≥3-nt loop (0-nt loop is sterically prohibited) → null.");
    }

    // E4 — a long oligo with a 6-bp stem closing a 10-nt loop. 5'-arm GGGGGG = 5 GG stacks
    // (ΔH°=5*-8.0=-40, ΔS°=5*-19.9=-99.5); loop-10 ΔG°37 = 4.6 (Table 4) → loop ΔS°=-4.6*1000/310.15.
    //   Total ΔS° = -114.33153312913106; ΔG°37 = -40 - 310.15*(-114.33153312913106)/1000 = -4.540075000000002.
    [Test]
    public void FindMostStableHairpin_LongStemAndLoop_MatchesHandDerived()
    {
        var hp = PrimerDesigner.FindMostStableHairpin("GGGGGG" + new string('A', 10) + "CCCCCC");

        Assert.That(hp, Is.Not.Null);
        var h = hp!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(h.StemLength, Is.EqualTo(6), "Full 6-bp GGGGGG/CCCCCC stem.");
            Assert.That(h.LoopSize, Is.EqualTo(10), "10-nt loop.");
            Assert.That(h.DeltaH, Is.EqualTo(-40.0).Within(Tol), "Five GG stacks, ΔH° = 5*-8.0.");
            Assert.That(h.DeltaS, Is.EqualTo(-114.33153312913106).Within(Tol),
                "Stem ΔS° (-99.5) + loop-10 ΔS° (-4.6*1000/310.15).");
            Assert.That(h.DeltaG37, Is.EqualTo(-4.540075000000002).Within(Tol), "ΔG°37 hand-derived (Table 4 loop-10 = 4.6).");
        });
    }

    #endregion
}
