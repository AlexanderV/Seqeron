namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the MolTools <b>salt-corrected nearest-neighbour melting temperature</b>
/// (PRIMER-NNTM-001) — the opt-in SantaLucia/Hicks (2004) NN Tm with a published monovalent
/// (Owczarzy 2004) or divalent (Owczarzy 2008) salt correction. The plain composition Tm
/// (PRIMER-TM-001, row 21) is covered separately in MolToolsFuzzTests; THIS file targets the
/// salt-correction and concentration parameters and the bimolecular Tm equation, NOT the basic
/// Wallace/Marmur-Doty model.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Malformed / boundary / out-of-domain inputs must NEVER hang, throw an *unhandled* runtime
/// exception, or emit out-of-contract output. Concretely, every input must resolve to EITHER a
/// well-defined, theory-correct result (a finite, physically-plausible Tm), OR a *documented*
/// sentinel / validation exception. For a thermodynamic Tm the central hazards are non-finite
/// leakage — a NaN/±Inf escaping where a finite Tm is contracted, or a finite Tm where the doc
/// says the result is not computable — and undisciplined ln-of-non-positive blow-ups.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PRIMER-NNTM-001 — salt-corrected nearest-neighbour melting temperature
/// Checklist: docs/checklists/03_FUZZING.md, row 240.
/// Algorithm doc: docs/algorithms/MolTools/NearestNeighbor_Salt_Corrected_Tm.md.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the 1-bp primer (NN needs ≥2 bases → not computable),
///          the ZERO-salt boundary (sodiumMolar = 0 → ln(0) = −∞ hazard), the zero-strand-
///          concentration boundary (ln(0) in the Tm equation), and a very long primer (no
///          overflow / non-finite / perf blow-up).
///   • INJ = Injection — non-DNA / special / unicode characters and an all-N primer (no
///          A/C/G/T dinucleotide step), the null reference, and NaN-valued concentration
///          parameters injected into the numeric API.
/// — docs/checklists/03_FUZZING.md §Description (BE; INJ = injection of special chars / null
///   bytes / unicode); row 240 targets:
///   "1-bp, all-N, non-DNA chars, zero salt, negative concentration, very long".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The salt-corrected NN-Tm contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: PrimerDesigner.CalculateMeltingTemperatureNN(string primer,
///   double strandConcentrationMolar = 0.5e-6, double sodiumMolar = 0.05,
///   double magnesiumMolar = 0.0, double dntpMolar = 0.0,
///   SaltCorrectionMode saltMode = Owczarzy2004Monovalent)  (PrimerDesigner.cs).
/// Supporting entry: PrimerDesigner.CalculateNearestNeighborThermodynamics(string) →
///   (ΔH°, ΔS°, IsSelfComplementary)? — null for a non-computable sequence.
///
/// The Tm equation (doc §2.2, SantaLucia &amp; Hicks 2004 Eq. 3):
///   Tm = ΔH°·1000 / (ΔS° + R·ln(C_T / x)) − 273.15,   R = 1.9872, x ∈ {1, 4}.
/// The monovalent (Owczarzy 2004) correction is a 1/Tm form in ln[Na⁺]; the divalent
/// (Owczarzy 2008) correction is a 1/Tm form in ln[Mg²⁺]. EVERY one of these — the Tm
/// equation's R·ln(C_T/x) and both salt corrections — evaluates a logarithm whose argument is
/// a concentration; a non-positive argument makes ln undefined (−∞ for 0, NaN for a negative).
///
/// SEQUENCE domain (doc §3.3, INV-06, §6.1): empty / null / length-1 / any non-ACGT character
/// makes the NN dinucleotide sum fail → CalculateNearestNeighborThermodynamics returns null and
/// CalculateMeltingTemperatureNN returns double.NaN. This is a DOCUMENTED sentinel, not a throw —
/// pinned here for 1-bp, all-N, non-DNA / unicode, and null.
///
/// PARAMETER domain (doc §3.1 constraints; §3.3): strandConcentrationMolar &gt; 0, sodiumMolar &gt; 0,
/// magnesiumMolar ≥ 0, dntpMolar ≥ 0. A non-positive C_T or [Na⁺] (zero or negative — and zero
/// salt is the KEY boundary, since ln(0) = −∞ would otherwise leak a non-physical ≈ −273.15 °C or
/// a silent NaN), a negative [Mg²⁺], or a negative dNTP each throw ArgumentOutOfRangeException.
/// NaN-valued concentrations are likewise out of domain (NaN &gt; 0 is false) and throw. So the
/// total contract is a clean trichotomy: finite theory-correct Tm | NaN sentinel (bad sequence) |
/// ArgumentOutOfRangeException (bad parameter) — never an undisciplined NaN/Inf escaping a Tm that
/// the doc contracts as finite. (Source fix made for row 240: the parameter-domain guards were
/// added to CalculateMeltingTemperatureNN; previously zero/negative salt or concentration leaked a
/// non-physical finite Tm or NaN. The doc §3.3 records the new throws.)
///
/// THEORY invariants pinned on valid input (doc §2.4):
///   INV-03 — over the CALIBRATED salt range ([Na⁺] 0.05–1.1 M, doc §6.2) lower [Na⁺] ⇒ lower
///            Owczarzy-2004 Tm, i.e. raising [Na⁺] within range does NOT lower Tm (Owczarzy 2004,
///            quadratic-correction sign). Asserted only on the published domain — NOT extrapolated
///            outside it (the campaign rule: no idealised monotonicity beyond the parameters'
///            validity).
///   INV-05 — the divalent mode with [Mg²⁺] = 0 reduces to the monovalent 2004 mode (Method-7
///            fallback) — so the two modes agree exactly at Mg²⁺ = 0.
/// Every valid-input test additionally asserts the Tm is FINITE (not NaN, not ±Infinity).
///
/// CalculateMeltingTemperatureNN is a pure function (no iterator), so every probe calls it
/// directly. The very-long-sequence probe is pinned with [CancelAfter] so any perf blow-up fails
/// as a timeout rather than hanging.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PrimerNnTmFuzzTests
{
    private const PrimerDesigner.SaltCorrectionMode None = PrimerDesigner.SaltCorrectionMode.None;
    private const PrimerDesigner.SaltCorrectionMode SantaLucia = PrimerDesigner.SaltCorrectionMode.SantaLuciaEntropy;
    private const PrimerDesigner.SaltCorrectionMode Owczarzy04 = PrimerDesigner.SaltCorrectionMode.Owczarzy2004Monovalent;
    private const PrimerDesigner.SaltCorrectionMode Owczarzy08 = PrimerDesigner.SaltCorrectionMode.Owczarzy2008Divalent;

    private static readonly PrimerDesigner.SaltCorrectionMode[] AllModes =
        { None, SantaLucia, Owczarzy04, Owczarzy08 };

    #region Helpers

    /// <summary>Deterministic ACGT generator — seed fixed locally so fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PRIMER-NNTM-001 — salt-corrected NN melting temperature : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PRIMER-NNTM-001 — BE: 1-bp sequence (NN needs ≥ 2 bases → not computable)

    /// <summary>
    /// BE (1-bp): the NN model sums ΔH°/ΔS° over adjacent dinucleotide STEPS, so it is undefined
    /// for a single base — there is no step. The doc (§3.3, INV-06, §6.1) contracts a length-&lt;2
    /// input as NOT computable: CalculateNearestNeighborThermodynamics returns null and
    /// CalculateMeltingTemperatureNN returns double.NaN — a DOCUMENTED sentinel, NOT a throw and
    /// NOT a fabricated finite Tm. Pinned for every base and every salt mode.
    /// </summary>
    [Test]
    public void SingleBasePrimer_NotComputable_ReturnsNaN_NoThrow()
    {
        foreach (var b in new[] { "A", "C", "G", "T" })
        {
            PrimerDesigner.CalculateNearestNeighborThermodynamics(b)
                .Should().BeNull($"a single base '{b}' has no NN dinucleotide step (INV-06)");

            foreach (var mode in AllModes)
            {
                double tm = PrimerDesigner.CalculateMeltingTemperatureNN(b, saltMode: mode);
                double.IsNaN(tm).Should().BeTrue(
                    $"a 1-bp primer '{b}' is not computable → NaN sentinel under mode {mode}, never a fabricated Tm");
            }
        }
    }

    #endregion

    #region PRIMER-NNTM-001 — INJ: all-N primer (no A/C/G/T dinucleotide step)

    /// <summary>
    /// INJ (all-N): 'N' is the IUPAC "any base" code, but the unified NN table keys are only the
    /// 16 ACGT dinucleotides — an 'N' produces no parameter, so the dinucleotide lookup fails and
    /// the sequence is NOT computable (doc INV-06, §6.1). An all-N primer therefore returns null
    /// thermo and a NaN Tm across every salt mode — never a crash, never a spurious finite Tm.
    /// </summary>
    [Test]
    public void AllNPrimer_NotComputable_ReturnsNaN_NoThrow()
    {
        var allN = new string('N', 20);

        PrimerDesigner.CalculateNearestNeighborThermodynamics(allN)
            .Should().BeNull("an all-N primer has no ACGT dinucleotide step (INV-06)");

        foreach (var mode in AllModes)
        {
            double tm = PrimerDesigner.CalculateMeltingTemperatureNN(allN, saltMode: mode);
            double.IsNaN(tm).Should().BeTrue($"an all-N primer is not computable → NaN under mode {mode}");
        }
    }

    #endregion

    #region PRIMER-NNTM-001 — INJ: non-DNA / special / unicode characters and null

    /// <summary>
    /// INJ (non-DNA / special / unicode / null): any character outside A/C/G/T breaks a
    /// dinucleotide-table lookup, so a primer containing junk — symbols, digits, whitespace,
    /// unicode, a null byte — is NOT computable and returns a NaN Tm (doc §3.3, INV-06). The null
    /// reference and the empty string are caught by the IsNullOrEmpty guard and likewise return NaN
    /// (NOT a NullReferenceException). Junk interleaved into an otherwise-valid primer ALSO breaks
    /// the lookup at the junk step → NaN. Pinned on the default Owczarzy-2004 mode and (for the
    /// interleaved case) all modes; every result is NaN, never a finite Tm, never a crash.
    /// </summary>
    [Test]
    public void NonDnaSpecialUnicodeAndNull_NotComputable_ReturnsNaN_NoThrow()
    {
        string?[] junk = BuildJunkInputs();

        foreach (var s in junk)
        {
            var act = () => PrimerDesigner.CalculateMeltingTemperatureNN(s!, saltMode: Owczarzy04);
            act.Should().NotThrow($"junk input {Describe(s)} must be handled by the NaN sentinel, not crash");

            double tm = PrimerDesigner.CalculateMeltingTemperatureNN(s!, saltMode: Owczarzy04);
            double.IsNaN(tm).Should().BeTrue($"non-ACGT / null / empty input {Describe(s)} is not computable → NaN");

            // The lowercase guard upper-cases first, so a *clean lowercase* primer is the control:
            // it is computable; only genuinely non-ACGT content yields NaN.
        }

        // Junk interleaved into a real primer breaks the lookup at the junk step → NaN, every mode.
        foreach (var mode in AllModes)
        {
            double tm = PrimerDesigner.CalculateMeltingTemperatureNN("ATGC#ATGC", saltMode: mode);
            double.IsNaN(tm).Should().BeTrue($"a '#' embedded in a real primer breaks the NN step → NaN under {mode}");
        }
    }

    /// <summary>
    /// Malformed / out-of-alphabet primer strings: null, empty, an unknown letter, embedded
    /// whitespace / hyphen / digit, the RNA base U (absent from the DNA NN table), an embedded
    /// null byte, and multi-byte unicode (a Latin diacritic, Greek letters, an emoji surrogate
    /// pair). Built programmatically so the source bytes never depend on literal encoding.
    /// </summary>
    private static string?[] BuildJunkInputs() => new[]
    {
        null,
        "",
        "ATGX",
        "AT GC",
        "AT-GC",
        "AT5GC",
        "ATGCU",                       // U is RNA, not in the DNA NN table
        "AT\0GC",                      // embedded null byte
        "ATGÄC",                  // U+00C4 (Latin capital A with diaeresis)
        "αβγδ",     // Greek letters alpha-beta-gamma-delta
        "atg\U0001F600c"               // emoji surrogate pair injected into a lowercase primer
    };

    private static string Describe(string? s) =>
        s is null ? "<null>" : s.Length == 0 ? "<empty>" : $"\"{s}\"";

    #endregion

    #region PRIMER-NNTM-001 — BE/INJ: ZERO salt (sodiumMolar = 0 → ln(0) hazard) — KEY

    /// <summary>
    /// BE/INJ (ZERO salt — KEY): sodiumMolar = 0 is the central ln(0) = −∞ hazard. The Owczarzy
    /// 2004/2008 corrections take ln[Na⁺] (and the SantaLucia entropy mode takes 0.368·(N/2)·ln[Na⁺]);
    /// at [Na⁺] = 0 each is −∞, which would leak a non-physical Tm (≈ −273.15 °C) or a silent NaN. The
    /// doc (§3.1: sodiumMolar &gt; 0; §3.3) requires this be rejected: a non-positive [Na⁺] throws
    /// ArgumentOutOfRangeException. Pinned on EVERY salt mode (the guard precedes mode selection), with
    /// a VALID sequence so this isolates the salt-domain boundary from the sequence sentinel — the
    /// result must be a documented throw, NOT a NaN and NOT a finite −273.15 °C.
    /// </summary>
    [Test]
    public void ZeroSodium_OutOfDomain_ThrowsArgumentOutOfRange_AllModes()
    {
        const string primer = "ATGCATGCATGC"; // valid, computable

        foreach (var mode in AllModes)
        {
            var act = () => PrimerDesigner.CalculateMeltingTemperatureNN(
                primer, sodiumMolar: 0.0, saltMode: mode);
            act.Should().Throw<ArgumentOutOfRangeException>(
                    $"zero [Na⁺] makes ln[Na⁺] = −∞ — rejected as out of domain under mode {mode}, never a leaked NaN/−273.15")
                .Which.ParamName.Should().Be("sodiumMolar");
        }
    }

    /// <summary>
    /// BE: a NEGATIVE [Na⁺] is doubly out of domain (negative concentration is unphysical AND
    /// ln(negative) = NaN). The doc constraint sodiumMolar &gt; 0 rejects it with
    /// ArgumentOutOfRangeException — NOT a NaN leak — across every salt mode.
    /// </summary>
    [Test]
    public void NegativeSodium_OutOfDomain_ThrowsArgumentOutOfRange_AllModes()
    {
        const string primer = "ATGCATGCATGC";

        foreach (var mode in AllModes)
        {
            var act = () => PrimerDesigner.CalculateMeltingTemperatureNN(
                primer, sodiumMolar: -0.05, saltMode: mode);
            act.Should().Throw<ArgumentOutOfRangeException>(
                    $"negative [Na⁺] is rejected under mode {mode} (ln(negative) = NaN otherwise)")
                .Which.ParamName.Should().Be("sodiumMolar");
        }
    }

    /// <summary>
    /// INJ: a NaN-valued [Na⁺] is out of domain — NaN &gt; 0 is false, so the constraint
    /// sodiumMolar &gt; 0 rejects it with ArgumentOutOfRangeException rather than propagating a NaN Tm.
    /// </summary>
    [Test]
    public void NaNSodium_OutOfDomain_ThrowsArgumentOutOfRange()
    {
        var act = () => PrimerDesigner.CalculateMeltingTemperatureNN(
            "ATGCATGCATGC", sodiumMolar: double.NaN, saltMode: Owczarzy04);
        act.Should().Throw<ArgumentOutOfRangeException>("a NaN [Na⁺] is out of the >0 domain")
            .Which.ParamName.Should().Be("sodiumMolar");
    }

    #endregion

    #region PRIMER-NNTM-001 — BE: zero / negative / NaN strand concentration (ln(C_T/x) hazard)

    /// <summary>
    /// BE (zero / negative / NaN C_T): the Tm equation divides by (ΔS° + R·ln(C_T/x)), so a
    /// non-positive C_T makes ln(C_T/x) = −∞ (zero) or NaN (negative), leaking a non-physical Tm.
    /// The doc constraint strandConcentrationMolar &gt; 0 (§3.1) rejects ≤ 0 and NaN with
    /// ArgumentOutOfRangeException — pinned for 0, a negative, and NaN, with a valid sequence so the
    /// concentration boundary is isolated from the sequence sentinel.
    /// </summary>
    [Test]
    public void NonPositiveOrNaNStrandConcentration_OutOfDomain_ThrowsArgumentOutOfRange()
    {
        const string primer = "ATGCATGCATGC";

        foreach (var ct in new[] { 0.0, -0.5e-6, double.NaN })
        {
            var act = () => PrimerDesigner.CalculateMeltingTemperatureNN(
                primer, strandConcentrationMolar: ct, saltMode: None);
            act.Should().Throw<ArgumentOutOfRangeException>(
                    $"a non-positive/NaN strand concentration ({ct}) makes ln(C_T/x) undefined → rejected")
                .Which.ParamName.Should().Be("strandConcentrationMolar");
        }
    }

    #endregion

    #region PRIMER-NNTM-001 — BE: negative magnesium / dNTP (divalent mode)

    /// <summary>
    /// BE: the divalent mode takes ln[Mg²⁺]; a NEGATIVE [Mg²⁺] is unphysical and would make
    /// ln(negative) = NaN. The doc constraint magnesiumMolar ≥ 0 (§3.1) rejects a negative value
    /// with ArgumentOutOfRangeException. The dNTP chelation term likewise requires dntpMolar ≥ 0.
    /// Pinned for a negative Mg²⁺ and a negative dNTP under the divalent mode (Mg²⁺ = 0 is VALID and
    /// is exercised by the INV-05 equivalence test, NOT here).
    /// </summary>
    [Test]
    public void NegativeMagnesiumOrDntp_OutOfDomain_ThrowsArgumentOutOfRange()
    {
        const string primer = "ATGCATGCATGC";

        var negMg = () => PrimerDesigner.CalculateMeltingTemperatureNN(
            primer, magnesiumMolar: -0.003, saltMode: Owczarzy08);
        negMg.Should().Throw<ArgumentOutOfRangeException>("a negative [Mg²⁺] is out of the ≥0 domain")
            .Which.ParamName.Should().Be("magnesiumMolar");

        var negDntp = () => PrimerDesigner.CalculateMeltingTemperatureNN(
            primer, magnesiumMolar: 0.003, dntpMolar: -0.0008, saltMode: Owczarzy08);
        negDntp.Should().Throw<ArgumentOutOfRangeException>("a negative dNTP concentration is out of the ≥0 domain")
            .Which.ParamName.Should().Be("dntpMolar");
    }

    #endregion

    #region PRIMER-NNTM-001 — Positive sanity: valid input → finite, physically-plausible Tm (all modes)

    /// <summary>
    /// Positive sanity: a normal primer under in-domain parameters yields a FINITE Tm (not NaN, not
    /// ±Infinity) in every salt mode, and the value is physically plausible — above absolute zero
    /// (&gt; −273.15 °C) and below the boiling point of water at 1 atm (&lt; 100 °C) for these short
    /// oligos at realistic buffers. This pins that the equation and each correction produce a real
    /// number for the canonical operating point, so the NaN assertions elsewhere are meaningful
    /// (a code that always returned NaN would fail HERE).
    /// </summary>
    [Test]
    public void ValidPrimer_AllModes_FiniteAndPhysicallyPlausibleTm()
    {
        foreach (var seed in new[] { 240_001, 240_002, 240_003 })
        {
            string primer = RandomDna(length: 22, seed: seed);

            foreach (var mode in AllModes)
            {
                double mg = mode == Owczarzy08 ? 0.003 : 0.0;
                double dntp = mode == Owczarzy08 ? 0.0008 : 0.0;

                double tm = PrimerDesigner.CalculateMeltingTemperatureNN(
                    primer, strandConcentrationMolar: 0.5e-6, sodiumMolar: 0.05,
                    magnesiumMolar: mg, dntpMolar: dntp, saltMode: mode);

                double.IsNaN(tm).Should().BeFalse($"{primer} under {mode} must yield a finite Tm, not NaN");
                double.IsInfinity(tm).Should().BeFalse($"{primer} under {mode} must yield a finite Tm, not ±Inf");
                tm.Should().BeGreaterThan(-273.15, $"a real Tm is above absolute zero ({mode})");
                tm.Should().BeLessThan(100.0, $"a ~22-nt oligo Tm at realistic buffer stays below 100 °C ({mode})");
            }
        }
    }

    #endregion

    #region PRIMER-NNTM-001 — Theory INV-03: higher [Na⁺] does not lower Owczarzy-2004 Tm (calibrated range)

    /// <summary>
    /// Theory (INV-03, doc §2.4 / §6.2): over the CALIBRATED monovalent range ([Na⁺] 0.05–1.1 M)
    /// the Owczarzy-2004 quadratic correction makes a LOWER [Na⁺] give a LOWER Tm — equivalently,
    /// raising [Na⁺] within range does NOT lower Tm (sodium screens the phosphate backbone and
    /// stabilises the duplex). Asserted ONLY across the published domain (0.05 → 0.5 → 1.0 M) — NOT
    /// extrapolated outside it, per the campaign rule against idealised monotonicity beyond the
    /// parameters' validity. Each Tm is also pinned FINITE. Several sequences are checked so the
    /// monotone trend is not a single-sequence accident.
    /// </summary>
    [Test]
    public void Owczarzy2004_HigherSodiumDoesNotLowerTm_OverCalibratedRange()
    {
        double[] saltSteps = { 0.05, 0.5, 1.0 }; // all within the documented 0.05–1.1 M range

        foreach (var seed in new[] { 240_101, 240_102, 240_103, 240_104 })
        {
            string primer = RandomDna(length: 20, seed: seed);

            double prev = double.NegativeInfinity;
            foreach (var na in saltSteps)
            {
                double tm = PrimerDesigner.CalculateMeltingTemperatureNN(
                    primer, sodiumMolar: na, saltMode: Owczarzy04);

                double.IsNaN(tm).Should().BeFalse($"{primer} @ [Na⁺]={na} must be finite");
                double.IsInfinity(tm).Should().BeFalse($"{primer} @ [Na⁺]={na} must be finite");
                tm.Should().BeGreaterThanOrEqualTo(prev - 1e-9,
                    $"INV-03: raising [Na⁺] to {na} M within the calibrated range must not LOWER Tm for {primer}");
                prev = tm;
            }
        }
    }

    #endregion

    #region PRIMER-NNTM-001 — Theory INV-05: divalent mode at [Mg²⁺]=0 ≡ monovalent 2004 mode

    /// <summary>
    /// Theory (INV-05, doc §2.4): the Owczarzy-2008 divalent mode with [Mg²⁺] = 0 (and no dNTP)
    /// falls back to the Owczarzy-2004 monovalent form — so the two modes must agree EXACTLY at
    /// Mg²⁺ = 0. This is the documented Method-7 fallback and is the boundary between the divalent
    /// and monovalent regimes; pinning equality here proves the divalent path's zero-Mg edge is the
    /// monovalent path, not a discontinuity. Both values are also finite.
    /// </summary>
    [Test]
    public void Owczarzy2008_ZeroMagnesium_EqualsOwczarzy2004Monovalent()
    {
        foreach (var seed in new[] { 240_201, 240_202, 240_203 })
        {
            string primer = RandomDna(length: 24, seed: seed);

            double mono = PrimerDesigner.CalculateMeltingTemperatureNN(
                primer, sodiumMolar: 0.05, saltMode: Owczarzy04);
            double diZeroMg = PrimerDesigner.CalculateMeltingTemperatureNN(
                primer, sodiumMolar: 0.05, magnesiumMolar: 0.0, dntpMolar: 0.0, saltMode: Owczarzy08);

            double.IsNaN(mono).Should().BeFalse();
            double.IsNaN(diZeroMg).Should().BeFalse();
            diZeroMg.Should().BeApproximately(mono, 1e-9,
                $"INV-05: divalent mode at [Mg²⁺]=0 reduces to the monovalent 2004 mode for {primer}");
        }
    }

    #endregion

    #region PRIMER-NNTM-001 — BE: very long primer (no overflow / non-finite / hang)

    /// <summary>
    /// BE (very long primer): a long primer drives the O(n) NN sum and the O(1) salt correction over
    /// a large length. The result must stay FINITE (the ΔH°/ΔS° sums grow linearly with no integer
    /// overflow — they are doubles — and the Tm equation/corrections stay well-defined), and the call
    /// must COMPLETE promptly (the algorithm is linear; pinned with [CancelAfter] so any blow-up fails
    /// as a timeout, not a hang). Exercised across all salt modes on a 5000-nt fixed-seed primer; each
    /// Tm is finite and above absolute zero. NOTE: the doc (§6.2) flags that the salt corrections are
    /// calibrated for short oligos — this test pins NUMERICAL discipline (finite, no overflow, no
    /// hang) on the long-primer path, NOT the accuracy of the extrapolated value.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void VeryLongPrimer_AllModes_FiniteNoOverflowNoHang()
    {
        string longPrimer = RandomDna(length: 5000, seed: 240_301);

        foreach (var mode in AllModes)
        {
            double mg = mode == Owczarzy08 ? 0.003 : 0.0;

            double tm = PrimerDesigner.CalculateMeltingTemperatureNN(
                longPrimer, sodiumMolar: 0.05, magnesiumMolar: mg, saltMode: mode);

            double.IsNaN(tm).Should().BeFalse($"a 5000-nt primer Tm under {mode} must be finite, not NaN");
            double.IsInfinity(tm).Should().BeFalse($"a 5000-nt primer Tm under {mode} must be finite, not ±Inf");
            tm.Should().BeGreaterThan(-273.15, $"a real Tm is above absolute zero ({mode})");
        }

        // The raw thermodynamics over the long primer are also finite (no overflow in the linear sum).
        var thermo = PrimerDesigner.CalculateNearestNeighborThermodynamics(longPrimer);
        thermo.Should().NotBeNull("a clean 5000-nt ACGT primer is computable");
        double.IsFinite(thermo!.Value.DeltaH).Should().BeTrue("ΔH° over 5000 nt is finite");
        double.IsFinite(thermo.Value.DeltaS).Should().BeTrue("ΔS° over 5000 nt is finite");
    }

    #endregion
}
