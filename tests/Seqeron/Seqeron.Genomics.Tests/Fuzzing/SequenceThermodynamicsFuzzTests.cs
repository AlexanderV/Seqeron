namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Statistics-area DNA nearest-neighbor thermodynamics unit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (KeyNotFoundException, IndexOutOfRange on
/// the L−1 dinucleotide window, DivideByZero / NaN / ±Infinity result, …). Every
/// input must yield EITHER a well-defined, theory-correct value, OR a *documented,
/// intentional* outcome. A raw runtime exception, a NaN ΔG/ΔH/ΔS, or a hang on
/// garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-THERMO-001 — DNA Duplex Thermodynamics, nearest-neighbor (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 129.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — empty / null (L=0 ⇒ no dinucleotide step, the
///          (L−1) lower boundary that must NOT crash or go negative-length),
///          single base (L=1 ⇒ 0 NN steps, the absent-step boundary), all-AT
///          (weakest stacking ⇒ least-negative ΔG), all-GC (strongest stacking ⇒
///          most-negative ΔG). Off-by-one in the L−1 step sum, DivideByZero on
///          empty/single, and KeyNotFound on an out-of-table dinucleotide are the
///          target failure modes.
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE = граничні значення).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The thermodynamics contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry (SequenceStatistics.cs, src/.../Seqeron.Genomics.Analysis):
///   • SequenceStatistics.CalculateThermodynamics(string, double naConc, double primerConc)
///         → ThermodynamicProperties(DeltaH, DeltaS, DeltaG, MeltingTemperature)
///
/// Documented behaviour (DNA_Thermodynamics.md, Test Unit ID SEQ-THERMO-001):
///   • §2.2 / §4.1 (core model): for a duplex of length N ≥ 2,
///         ΔH° = Σ ΔH°(NN step over the N−1 overlapping dinucleotides)
///               + ΔH°(init, 5′ terminus) + ΔH°(init, 3′ terminus)
///         ΔS° = Σ ΔS°(NN step) + ΔS°(init, 5′) + ΔS°(init, 3′)
///                + 0.368·(N−1)·ln[Na⁺]                     (salt correction, method 5)
///         ΔG°₃₇ = ΔH° − 310.15·ΔS°/1000                    (INV-02)
///         Tm = (1000·ΔH°)/(ΔS° + R·ln(C_T/4)) − 273.15, R = 1.987   (INV-03)
///   • §2.4 INV-01: initiation contributes ONCE at EACH terminus (two init terms),
///     chosen by whether that terminal base is G·C (+0.1, −2.8) or A·T (+2.3, +4.1).
///   • §2.4 INV-06 / §6.1: empty or length-1 input returns (0,0,0,0) — NN undefined
///     for length < 2 (no dinucleotide); guarded BEFORE any (L−1) indexing.
///   • §3.3 / §6.1: only A/C/G/T recognized; an unrecognized dinucleotide adds 0 via
///     TryGetValue — NO KeyNotFound, no exception for valid-typed input.
///   • §2.4 INV-04: NN table is Watson-Crick symmetric (AA=TT, CA=TG, GT=AC, …).
///   • §2.4 INV-05: deterministic and case-insensitive (input ToUpperInvariant'd).
///
/// The NN table, initiation terms and salt/Gibbs/Tm constants pinned below are an
/// INDEPENDENT oracle copy of the documented constants (DNA_Thermodynamics.md §4.2),
/// so the tests assert the real BUSINESS formula, not merely non-throwing.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceThermodynamicsFuzzTests
{
    #region Helpers — oracle NN table + constants (independent copy of DNA_Thermodynamics.md §4.2)

    private const double Tolerance = 1e-9;

    // Default API parameters (DNA_Thermodynamics.md §3.1).
    private const double DefaultNa = 0.05;          // 50 mM Na⁺
    private const double DefaultPrimer = 0.00000025; // 250 nM total strand C_T

    private const double SaltCoeff = 0.368;          // method 5 salt coefficient
    private const double GasConstant = 1.987;        // R, cal/(mol·K)
    private const double RefTempK = 310.15;          // 37 °C
    private const double KelvinOffset = 273.15;
    private const double NonSelfFactor = 4.0;        // F

    /// <summary>NN ΔH° (kcal/mol) / ΔS° (cal/(mol·K)) — oracle copy of the documented §4.2 table,
    /// with both Watson-Crick complements listed (INV-04).</summary>
    private static readonly IReadOnlyDictionary<string, (double dH, double dS)> Nn =
        new Dictionary<string, (double, double)>
        {
            { "AA", (-7.9, -22.2) }, { "TT", (-7.9, -22.2) },
            { "AT", (-7.2, -20.4) },
            { "TA", (-7.2, -21.3) },
            { "CA", (-8.5, -22.7) }, { "TG", (-8.5, -22.7) },
            { "GT", (-8.4, -22.4) }, { "AC", (-8.4, -22.4) },
            { "CT", (-7.8, -21.0) }, { "AG", (-7.8, -21.0) },
            { "GA", (-8.2, -22.2) }, { "TC", (-8.2, -22.2) },
            { "CG", (-10.6, -27.2) },
            { "GC", (-9.8, -24.4) },
            { "GG", (-8.0, -19.9) }, { "CC", (-8.0, -19.9) }
        };

    // Two-terminus initiation by G·C vs A·T (DNA_Thermodynamics.md §4.2, INV-01).
    private static (double dH, double dS) Init(char terminal) =>
        terminal is 'G' or 'C' ? (0.1, -2.8) : (2.3, 4.1);

    /// <summary>Independent oracle that mirrors the documented contract exactly:
    /// two-terminus init + Σ NN over the N−1 overlapping dinucleotides + salt correction,
    /// then ΔG and Tm. Returns the API's rounding (ΔH/ΔS/ΔG to 2 dp, Tm to 1 dp).</summary>
    private static (double dH, double dS, double dG, double tm) Oracle(
        string seq, double na = DefaultNa, double primer = DefaultPrimer)
    {
        if (string.IsNullOrEmpty(seq) || seq.Length < 2)
            return (0, 0, 0, 0);

        string s = seq.ToUpperInvariant();
        double dH = 0, dS = 0;

        var (h5, s5) = Init(s[0]);
        var (h3, s3) = Init(s[^1]);
        dH += h5 + h3;
        dS += s5 + s3;

        for (int i = 0; i < s.Length - 1; i++)
            if (Nn.TryGetValue(s.Substring(i, 2), out var p)) { dH += p.dH; dS += p.dS; }

        dS += SaltCoeff * (s.Length - 1) * Math.Log(na);

        double dG = dH - (RefTempK * dS / 1000.0);
        double tm = (dH * 1000) / (dS + GasConstant * Math.Log(primer / NonSelfFactor)) - KelvinOffset;

        return (Math.Round(dH, 2), Math.Round(dS, 2), Math.Round(dG, 2), Math.Round(tm, 1));
    }

    /// <summary>The universal well-formedness contract for ANY input: ΔH, ΔS and ΔG are finite
    /// numbers (never NaN / ±Infinity). For length &lt; 2 the documented sentinel is all-zero;
    /// for a real duplex (N ≥ 2 with recognized steps) ΔH and ΔS are strictly negative (every
    /// table value and the GC init ΔH aside is negative, stacking is exothermic). We only assert
    /// finiteness universally here — the sign/magnitude assertions live in the dedicated tests
    /// where the duplex shape is known.</summary>
    private static void AssertFinite(SequenceStatistics.ThermodynamicProperties t)
    {
        double.IsFinite(t.DeltaH).Should().BeTrue("ΔH must be finite — no NaN/±Infinity (§8 fuzzing contract)");
        double.IsFinite(t.DeltaS).Should().BeTrue("ΔS must be finite — no NaN/±Infinity (§8 fuzzing contract)");
        double.IsFinite(t.DeltaG).Should().BeTrue("ΔG must be finite — no NaN/±Infinity (§8 fuzzing contract)");
        double.IsFinite(t.MeltingTemperature).Should().BeTrue("Tm must be finite — no NaN/±Infinity");
    }

    private static void ShouldMatchOracle(
        SequenceStatistics.ThermodynamicProperties t, string seq,
        double na = DefaultNa, double primer = DefaultPrimer)
    {
        var (dH, dS, dG, tm) = Oracle(seq, na, primer);
        t.DeltaH.Should().BeApproximately(dH, Tolerance, $"ΔH for '{seq}' = Σ NN + 2·init (§2.2)");
        t.DeltaS.Should().BeApproximately(dS, Tolerance, $"ΔS for '{seq}' = Σ NN + 2·init + salt (§2.2)");
        t.DeltaG.Should().BeApproximately(dG, Tolerance, $"ΔG₃₇ for '{seq}' = ΔH − 310.15·ΔS/1000 (INV-02)");
        t.MeltingTemperature.Should().BeApproximately(tm, Tolerance, $"Tm for '{seq}' (INV-03)");
    }

    /// <summary>Random A/C/G/T string of the given length (a well-formed DNA duplex).</summary>
    private static string RandomDna(Random rng, int length)
    {
        const string bases = "ACGT";
        var c = new char[length];
        for (int i = 0; i < length; i++) c[i] = bases[rng.Next(bases.Length)];
        return new string(c);
    }

    /// <summary>Random arbitrary BMP code points (control chars, null byte, lone surrogate
    /// halves, unicode) — fuzz fodder whose dinucleotides are absent from the NN table.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var c = new char[length];
        for (int i = 0; i < length; i++) c[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(c);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-THERMO-001 — nearest-neighbor ΔG/ΔH/ΔS : fuzz targets (BE)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-computed exact ΔH/ΔS/ΔG/Tm

    /// <summary>
    /// Positive baseline (not a boundary): the documented worked example "GCGC" must reproduce
    /// EXACTLY, including the two-terminus initiation, the three NN steps and the salt term.
    ///   • init (G end + C end): ΔH = 2×0.1 = 0.2;  ΔS = 2×(−2.8) = −5.6
    ///   • NN GC,CG,GC:          ΔH = −30.2;        ΔS = −76.0
    ///   • ΔH° = −30.0;  salt ΔS = 0.368·3·ln(0.05) ⇒ ΔS° = −84.91
    ///   • ΔG°₃₇ = −30.0 − 310.15·(−84.91)/1000 = −3.67;  Tm = −18.6 °C
    /// Confirms the suite asserts the real BUSINESS formula, not merely non-throwing.
    /// — DNA_Thermodynamics.md §7.1 worked example.
    /// </summary>
    [Test]
    public void Thermo_DocumentedWorkedExample_GCGC_MatchesHandComputedExactly()
    {
        var t = SequenceStatistics.CalculateThermodynamics("GCGC");

        t.DeltaH.Should().BeApproximately(-30.0, 1e-9, "0.2 (2·init G·C) − 30.2 (GC+CG+GC) = −30.0 (§7.1)");
        t.DeltaS.Should().BeApproximately(-84.91, 1e-9, "−5.6 − 76.0 + 0.368·3·ln(0.05) = −84.91 (§7.1)");
        t.DeltaG.Should().BeApproximately(-3.67, 1e-9, "−30.0 − 310.15·(−84.91)/1000 = −3.67 (§7.1)");
        t.MeltingTemperature.Should().BeApproximately(-18.6, 1e-9, "Tm = −18.6 °C (§7.1)");
    }

    /// <summary>
    /// Positive baseline (INV-01): initiation is added at BOTH termini, not once. A shortest
    /// duplex "GC" carries exactly two init terms plus a single NN step:
    ///   • init G + init C: ΔH = 0.2, ΔS = −5.6;  NN "GC": ΔH = −9.8, ΔS = −24.4
    ///   • ΔH° = −9.6;  salt ΔS = 0.368·1·ln(0.05) ⇒ ΔS° = −31.10
    /// Guards against the historical single-terminus regression (§5.4 deviation #2).
    /// — DNA_Thermodynamics.md §2.4 INV-01, §4.2.
    /// </summary>
    [Test]
    public void Thermo_ShortestDuplex_GC_HasTwoInitiationTerms()
    {
        var t = SequenceStatistics.CalculateThermodynamics("GC");

        t.DeltaH.Should().BeApproximately(-9.6, 1e-9, "−9.8 (NN GC) + 0.2 (two G·C init) = −9.6 (INV-01)");
        t.DeltaS.Should().BeApproximately(-31.10, 1e-9, "−24.4 + 2·(−2.8) + 0.368·1·ln(0.05) = −31.10 (INV-01)");
        ShouldMatchOracle(t, "GC");
    }

    /// <summary>
    /// Positive baseline (INV-05): the model is case-insensitive — lowercase input yields the
    /// IDENTICAL ΔH/ΔS/ΔG/Tm as upper-case (input ToUpperInvariant'd before lookup).
    /// — DNA_Thermodynamics.md §2.4 INV-05, §3.3.
    /// </summary>
    [Test]
    public void Thermo_LowercaseInput_EqualsUppercase()
    {
        var lower = SequenceStatistics.CalculateThermodynamics("gcgatatcgc");
        var upper = SequenceStatistics.CalculateThermodynamics("GCGATATCGC");

        lower.Should().Be(upper, "INV-05: thermodynamics are case-insensitive");
    }

    /// <summary>
    /// Positive baseline (INV-04): the NN table is Watson-Crick symmetric, so a sequence and the
    /// reverse complement of its complement strand read 5′→3′ share the same NN set. Concretely
    /// AA and TT, CA and TG etc. carry identical step values, so two sequences differing only by
    /// such complement substitutions (with matching termini) yield identical ΔH/ΔS. Here "AACC"
    /// and "GGTT" are reverse complements: same NN multiset {AA/TT, AC/GT, CC/GG} and same A·T /
    /// G·C terminus pattern ⇒ identical thermodynamics. — DNA_Thermodynamics.md §2.4 INV-04.
    /// </summary>
    [Test]
    public void Thermo_WatsonCrickSymmetry_ReverseComplementMatches()
    {
        var fwd = SequenceStatistics.CalculateThermodynamics("AACC");
        var rc = SequenceStatistics.CalculateThermodynamics("GGTT");

        rc.Should().Be(fwd, "INV-04: reverse-complement shares the WC-symmetric NN set and terminus pattern");
    }

    #endregion

    #region BE — Boundary: empty / null (length 0 ⇒ (0,0,0,0), no L−1 crash)

    /// <summary>
    /// BE: the empty string is the L=0 lower size boundary — there is NO dinucleotide and the
    /// (L−1) step loop must not iterate to index −1. The documented result is the (0,0,0,0)
    /// sentinel, reached via the `length &lt; 2` guard BEFORE any indexing — crucially NOT an
    /// IndexOutOfRange on s[0]/s[^1], not a DivideByZero, and not a NaN.
    /// — DNA_Thermodynamics.md §6.1 / INV-06.
    /// </summary>
    [Test]
    public void Thermo_EmptyString_IsAllZero_NoCrash()
    {
        var act = () => SequenceStatistics.CalculateThermodynamics(string.Empty);
        act.Should().NotThrow("L=0 is guarded before any (L−1) indexing (INV-06)");

        var t = act();
        t.Should().Be(new SequenceStatistics.ThermodynamicProperties(0, 0, 0, 0),
            "no dinucleotide step ⇒ documented (0,0,0,0) sentinel (INV-06)");
        AssertFinite(t);
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit) — (0,0,0,0), no
    /// NullReferenceException. — DNA_Thermodynamics.md §3.3 (null/empty ⇒ 0), §6.1.
    /// </summary>
    [Test]
    public void Thermo_Null_IsAllZero_NoThrow()
    {
        var act = () => SequenceStatistics.CalculateThermodynamics(null!);
        act.Should().NotThrow("null is 'no sequence', not an error");
        act().Should().Be(new SequenceStatistics.ThermodynamicProperties(0, 0, 0, 0));
    }

    #endregion

    #region BE — Boundary: single base (length 1 ⇒ 0 NN steps, no absent-step crash)

    /// <summary>
    /// BE: a single base is the L=1 boundary — there are ZERO nearest-neighbor steps (the step
    /// loop runs 0..(L−1−1) = empty) and the NN model is undefined. The documented result is the
    /// (0,0,0,0) sentinel via the `length &lt; 2` guard: no init term is applied (the duplex
    /// concept needs two paired ends), no s.Substring(i,2) on a length-1 string, no DivideByZero
    /// from ΔS=0 in the Tm denominator. Verified for every base, upper and lower case.
    /// — DNA_Thermodynamics.md §6.1 / INV-06.
    /// </summary>
    [TestCase("A")]
    [TestCase("C")]
    [TestCase("G")]
    [TestCase("T")]
    [TestCase("a")]
    [TestCase("g")]
    public void Thermo_SingleBase_IsAllZero_NoCrash(string seq)
    {
        var act = () => SequenceStatistics.CalculateThermodynamics(seq);
        act.Should().NotThrow($"L=1 '{seq}' has 0 NN steps and is guarded (INV-06)");

        var t = act();
        t.Should().Be(new SequenceStatistics.ThermodynamicProperties(0, 0, 0, 0),
            $"single base '{seq}' ⇒ no dinucleotide ⇒ (0,0,0,0) (INV-06)");
        AssertFinite(t);
    }

    /// <summary>
    /// BE: a single UNRECOGNIZED character is the L=1 boundary where the only base is also out of
    /// the alphabet — still the (0,0,0,0) sentinel, with no KeyNotFound and no crash, because the
    /// length guard fires before any table lookup. — DNA_Thermodynamics.md §6.1, §3.3.
    /// </summary>
    [TestCase("N")]
    [TestCase("-")]
    [TestCase("*")]
    [TestCase("1")]
    public void Thermo_SingleUnrecognizedChar_IsAllZero_NoThrow(string seq)
    {
        var act = () => SequenceStatistics.CalculateThermodynamics(seq);
        act.Should().NotThrow($"'{seq}' is length-1 and guarded before any lookup");
        act().Should().Be(new SequenceStatistics.ThermodynamicProperties(0, 0, 0, 0));
    }

    #endregion

    #region BE — Boundary: all-AT (weakest stacking, least-negative ΔG)

    /// <summary>
    /// BE: an all-AT duplex is the GC-content lower boundary — the weakest-stacking extreme. The
    /// documented A/T-rich thermodynamics must match the oracle EXACTLY over the AA/TA/AT NN steps
    /// plus two A·T init terms. Hand pins for the canonical forms:
    ///   • "ATAT": init A+T (ΔH 4.6, ΔS 8.2) + AT,TA,AT (ΔH −21.6, ΔS −62.1) + salt
    ///             ⇒ ΔH −17.0, ΔS −57.21, ΔG +0.74, Tm −84.6
    ///   • "AAAA": init A+A + AA,AA,AA ⇒ ΔH −19.1, ΔS −61.71, ΔG +0.04
    /// Both are checked against the independent oracle for arbitrary lengths.
    /// — DNA_Thermodynamics.md §2.2 / §4.2 (AA, AT, TA steps; A·T init).
    /// </summary>
    [Test]
    public void Thermo_AllAT_MatchesDocumentedATRichThermodynamics()
    {
        var atat = SequenceStatistics.CalculateThermodynamics("ATAT");
        atat.DeltaH.Should().BeApproximately(-17.0, 1e-9, "4.6 (2·A·T init) − 21.6 (AT+TA+AT) = −17.0");
        atat.DeltaS.Should().BeApproximately(-57.21, 1e-9, "8.2 − 62.1 + 0.368·3·ln(0.05) = −57.21");
        atat.DeltaG.Should().BeApproximately(0.74, 1e-9, "−17.0 − 310.15·(−57.21)/1000 = +0.74");
        ShouldMatchOracle(atat, "ATAT");

        var aaaa = SequenceStatistics.CalculateThermodynamics("AAAA");
        aaaa.DeltaH.Should().BeApproximately(-19.1, 1e-9, "4.6 (2·A·T init) − 23.7 (AA·3) = −19.1");
        ShouldMatchOracle(aaaa, "AAAA");

        foreach (string seq in new[] { "AT", "TA", "AAAAAAAA", "ATATATAT", "TATATA", "AATTAATT" })
            ShouldMatchOracle(SequenceStatistics.CalculateThermodynamics(seq), seq);
    }

    #endregion

    #region BE — Boundary: all-GC (strongest stacking, most-negative ΔG)

    /// <summary>
    /// BE: an all-GC duplex is the GC-content upper boundary — the strongest-stacking extreme. The
    /// documented G/C-rich thermodynamics must match the oracle EXACTLY over the GG/CG/GC/CC NN
    /// steps plus two G·C init terms. Hand pins for the canonical forms:
    ///   • "GCGC": ΔH −30.0, ΔS −84.91, ΔG −3.67 (the §7.1 worked example)
    ///   • "GGGG": init G+G + GG·3 ⇒ ΔH −23.8, ΔG −2.52
    /// Both checked against the independent oracle for arbitrary lengths.
    /// — DNA_Thermodynamics.md §2.2 / §4.2 (GG, GC, CG steps; G·C init).
    /// </summary>
    [Test]
    public void Thermo_AllGC_MatchesDocumentedGCRichThermodynamics()
    {
        var gggg = SequenceStatistics.CalculateThermodynamics("GGGG");
        gggg.DeltaH.Should().BeApproximately(-23.8, 1e-9, "0.2 (2·G·C init) − 24.0 (GG·3) = −23.8");
        gggg.DeltaG.Should().BeApproximately(-2.52, 1e-9, "−23.8 − 310.15·ΔS/1000 = −2.52");
        ShouldMatchOracle(gggg, "GGGG");

        foreach (string seq in new[] { "GC", "CG", "GGGGGGGG", "GCGCGCGC", "CGCGCG", "GGCCGGCC" })
            ShouldMatchOracle(SequenceStatistics.CalculateThermodynamics(seq), seq);
    }

    /// <summary>
    /// BE (the headline physical contrast): a GC-rich duplex is MORE stable than an AT-rich duplex
    /// of EQUAL length — strictly more-negative ΔG and ΔH (stronger stacking) and a higher Tm.
    /// This is the cross-boundary invariant the two GC-content extremes exist to expose; a sign
    /// flip here would mean the AT/GC step values were swapped or the init term mis-assigned.
    /// — DNA_Thermodynamics.md §2.1 (GC pairs stack more strongly), §4.2.
    /// </summary>
    [Test]
    public void Thermo_GCRich_MoreStableThan_ATRich_OfEqualLength()
    {
        foreach (int n in new[] { 4, 8, 16, 40 })
        {
            string gc = string.Concat(System.Linq.Enumerable.Repeat("GC", n / 2));
            string at = string.Concat(System.Linq.Enumerable.Repeat("AT", n / 2));

            var tgc = SequenceStatistics.CalculateThermodynamics(gc);
            var tat = SequenceStatistics.CalculateThermodynamics(at);

            tgc.DeltaG.Should().BeLessThan(tat.DeltaG,
                $"GC-rich (n={n}) is more stable: ΔG more negative than AT-rich");
            tgc.DeltaH.Should().BeLessThan(tat.DeltaH,
                $"GC stacking is more exothermic: ΔH more negative (n={n})");
            tgc.MeltingTemperature.Should().BeGreaterThan(tat.MeltingTemperature,
                $"GC-rich duplex melts at a higher Tm (n={n})");
        }
    }

    #endregion

    #region BE / MC — out-of-table dinucleotides: contribute 0, never KeyNotFound

    /// <summary>
    /// BE/contract: a dinucleotide ABSENT from the NN table contributes 0 via TryGetValue — no
    /// KeyNotFound. A length-≥2 sequence of all-unrecognized bases therefore has only the two
    /// init terms (computed from the actual terminal chars) plus the salt term, exactly as the
    /// oracle predicts, with finite output. — DNA_Thermodynamics.md §3.3, §6.1, §5.2 (TryGetValue).
    /// </summary>
    [TestCase("NN")]
    [TestCase("NNNN")]
    [TestCase("----")]
    [TestCase("XYZW")]
    public void Thermo_OutOfTableDinucleotides_ContributeZero_NoKeyNotFound(string seq)
    {
        var act = () => SequenceStatistics.CalculateThermodynamics(seq);
        act.Should().NotThrow($"'{seq}' dinucleotides are not in the NN table — TryGetValue, no KeyNotFound");

        var t = act();
        ShouldMatchOracle(t, seq);
        AssertFinite(t);
    }

    /// <summary>
    /// MC: a recognized duplex with interior junk — the junk-spanning dinucleotides are simply
    /// absent from the table and add 0, so the result equals the oracle (which models the exact
    /// same TryGetValue skipping over the L−1 windows). Pins that the L−1 step loop never throws
    /// on a window straddling an unknown base. — DNA_Thermodynamics.md §3.3 / §6.1.
    /// </summary>
    [Test]
    public void Thermo_RecognizedDuplexWithInteriorJunk_MatchesOracle()
    {
        foreach (string seq in new[] { "GC-GC", "AT*AT", "GNG", "ACGTNACGT", "G C G" })
        {
            var t = SequenceStatistics.CalculateThermodynamics(seq);
            ShouldMatchOracle(t, seq);
            AssertFinite(t);
        }
    }

    #endregion

    #region BE — very long duplex: linear, finite, terminates

    /// <summary>
    /// BE: a very long homopolymer duplex must remain finite and terminate quickly (O(n), one
    /// pass over n−1 dinucleotides) — no overflow in the double accumulation, no hang. The exact
    /// closed form is asserted against the oracle for n up to 1,000,000.
    /// — DNA_Thermodynamics.md §4.3 (O(n), O(1)).
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Thermo_VeryLongHomopolymer_IsFinite_MatchesOracle()
    {
        foreach (int n in new[] { 1000, 100_000, 1_000_000 })
        {
            string gc = new string('G', n);
            var t = SequenceStatistics.CalculateThermodynamics(gc);
            AssertFinite(t);

            var (dH, _, _, _) = Oracle(gc);
            // ΔH grows linearly: ~ (n−1)·(−8.0 GG) + 0.2 init; pin to the oracle exactly.
            t.DeltaH.Should().BeApproximately(dH, 1e-3, $"ΔH for poly-G n={n} matches oracle (linear, no overflow)");
        }
    }

    #endregion

    #region BE / MC — random garbage: never throws, always finite, matches oracle

    /// <summary>
    /// MC/BE: a large batch of arbitrary BMP strings (control chars, the null byte, lone surrogate
    /// halves, unicode, occasionally seeded with real bases) of length 0..200 must NEVER throw and
    /// ALWAYS yield finite ΔH/ΔS/ΔG/Tm that match the independent oracle. Core fuzz guarantee: no
    /// KeyNotFound on an out-of-table dinucleotide, no IndexOutOfRange on the L−1 window, no NaN.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Thermo_RandomGarbageStrings_NeverThrow_MatchOracle()
    {
        var rng = new Random(129_001);

        for (int iteration = 0; iteration < 3000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);

            SequenceStatistics.ThermodynamicProperties t = default;
            var act = () => t = SequenceStatistics.CalculateThermodynamics(input);
            act.Should().NotThrow($"garbage (len {len}) must never crash thermodynamics");

            ShouldMatchOracle(t, input);
            AssertFinite(t);
        }
    }

    /// <summary>
    /// MC/BE: randomly built genuine DNA duplexes (A/C/G/T, length 2..300, varied GC content and
    /// salt) must equal the independent oracle EXACTLY over many shapes, and every recognized
    /// duplex must have strictly-negative ΔH and ΔS (all NN steps are exothermic / entropy-losing;
    /// the small +0.1 G·C and +2.3 A·T init ΔH cannot outweigh the ≥1 negative step). Confirms the
    /// L−1 step sum and two-terminus init are correct across the full random space.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Thermo_RandomDnaDuplexes_MatchOracle_AndNegativeEnthalpyEntropy()
    {
        var rng = new Random(129_002);

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(2, 300);
            string seq = RandomDna(rng, len);
            double na = 0.01 + rng.NextDouble() * 0.99; // 10 mM .. ~1 M, always > 0

            var t = SequenceStatistics.CalculateThermodynamics(seq, na, DefaultPrimer);

            ShouldMatchOracle(t, seq, na, DefaultPrimer);
            AssertFinite(t);

            t.DeltaH.Should().BeLessThan(0, "every recognized DNA duplex has exothermic (negative) ΔH");
            t.DeltaS.Should().BeLessThan(0, "every recognized DNA duplex loses entropy on duplex formation (negative ΔS)");
        }
    }

    #endregion
}
