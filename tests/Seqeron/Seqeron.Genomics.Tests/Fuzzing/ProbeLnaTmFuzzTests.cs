namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the MolTools <b>LNA-adjusted nearest-neighbour melting temperature</b>
/// (PROBE-LNATM-001) — the opt-in McTigue, Peterson &amp; Kahn (2004) LNA-DNA NN increments
/// added to the SantaLucia (1998) DNA NN stack for a DNA hybridisation probe carrying one or
/// more <b>internal</b> LNA (locked nucleic acid) substitutions. The plain composition Tm
/// (PRIMER-TM-001, row 21) and the salt-corrected perfect-match NN Tm (PRIMER-NNTM-001,
/// row 240) are covered separately; THIS file targets the LNA-position mask and the McTigue
/// increment application, NOT the underlying DNA NN model.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Malformed / boundary / out-of-domain inputs must NEVER hang, throw an *unhandled* runtime
/// exception (an IndexOutOfRange from a bad LNA position, a NaN leaking where a finite Tm is
/// contracted), or emit out-of-contract output. Every input must resolve to EITHER a
/// well-defined, theory-correct result (a finite Tm, OR the documented null/NaN sentinel),
/// OR a *documented* validation exception (ArgumentNullException on a null position set). For
/// the LNA Tm the central hazards are (a) an LNA mask whose positions index out of range or
/// land on a terminal base — which McTigue (2004) does NOT parameterise — leaking an
/// IndexOutOfRange or a fabricated finite Tm instead of the documented not-computable sentinel,
/// and (b) a non-finite leak on an empty / sub-2-base / non-ACGT probe.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROBE-LNATM-001 — LNA-adjusted nearest-neighbour melting temperature
/// Checklist: docs/checklists/03_FUZZING.md, row 243.
/// Algorithm doc: docs/algorithms/MolTools/LNA_Adjusted_Nearest_Neighbor_Tm.md (PROBE-DESIGN-001).
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — an INVALID LNA mask (positions out of range, negative, or
///          landing on a TERMINAL base 0 / length−1 → McTigue has no terminal parameter →
///          documented not-computable), a mask LONGER than the probe (more positions than
///          bases), and non-DNA / unicode junk characters in the probe (the DNA NN lookup
///          fails → not computable).
///   • BE = Boundary Exploitation — the EMPTY probe and the 1-bp probe (the NN model needs ≥ 2
///          bases → not computable), the EMPTY mask (must reproduce the perfect-match NN Tm
///          exactly, INV-01), and the smallest internal-LNA-capable probe (length 3, only
///          index 1 is internal).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 243 targets:
///   "invalid LNA mask, mask longer than probe, empty probe, non-DNA".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The LNA-adjusted NN-Tm contract under test (LNA_Adjusted_Nearest_Neighbor_Tm.md)
/// ───────────────────────────────────────────────────────────────────────────
/// API entries (PrimerDesigner.cs):
///   • CalculateNearestNeighborThermodynamicsLna(string sequence,
///       IReadOnlyCollection&lt;int&gt; lnaPositions) → (ΔH°, ΔS°, IsSelfComplementary)? —
///       null when the duplex is not computable.
///   • CalculateMeltingTemperatureNNLna(string sequence, IReadOnlyCollection&lt;int&gt;
///       lnaPositions, double strandConcentrationMolar = 0.5e-6, double sodiumMolar = 0.05,
///       double magnesiumMolar = 0, double dntpMolar = 0,
///       SaltCorrectionMode saltMode = Owczarzy2004Monovalent) → Tm °C, or double.NaN when
///       not computable.
///
/// The model (doc §2.2): the base DNA NN ΔH°/ΔS° (SantaLucia 1998 unified) is computed first,
/// then for each nearest-neighbour STEP (i, i+1) that contains an LNA base the verbatim McTigue
/// (2004) increment ΔΔH°/ΔΔS° is ADDED (5'-locked increment if base i is locked, 3'-locked
/// increment if base i+1 is locked). The bimolecular Tm equation
///   Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15,  R = 1.9872, x ∈ {1, 4}
/// is then applied identically to the perfect-match NN Tm.
///
/// DOMAIN / VALIDATION (doc §3.3, §6.1; source CalculateNearestNeighborThermodynamicsLna):
///   • null lnaPositions → ArgumentNullException (ThrowIfNull) — the ONLY throw.
///   • null / empty / &lt; 2 bases / any non-ACGT character → the base DNA NN returns null, so the
///     LNA result is null / NaN (DOCUMENTED sentinel, NOT a throw).
///   • ANY LNA position ≤ 0 OR ≥ length−1 — i.e. a TERMINAL base (index 0 or last) or an
///     OUT-OF-RANGE index (negative, or ≥ length) — makes the result not computable (null / NaN),
///     because McTigue (2004) did NOT parameterise a terminal LNA (ASM-01, INV-03). This is the
///     central MC/BE probe: a bad mask must yield the sentinel, NEVER an IndexOutOfRangeException
///     and NEVER a fabricated finite Tm.
///   • Order and DUPLICATE positions are tolerated (set semantics, §6.1).
///
/// ───────────────────────────────────────────────────────────────────────────
/// THEORY invariants pinned on valid input (doc §2.4) — and what we deliberately DO NOT assert
/// ───────────────────────────────────────────────────────────────────────────
///   • INV-01 — an EMPTY LNA-position set reproduces the perfect-match NN Tm EXACTLY (no
///     increment added). Pinned by equality vs CalculateMeltingTemperatureNN.
///   • A valid internal-LNA probe yields a FINITE Tm (not NaN, not ±Infinity).
///   • The increment IS actually applied: a valid internal LNA changes the Tm vs the bare-DNA
///     Tm (DIFFERS) — we assert it differs WITHOUT claiming a sign in general.
///   • LNA-C is the reliably-STABILISING case: all 8 C-locked McTigue steps are net-stabilising,
///     so substituting an internal LNA-C does NOT LOWER the Tm vs bare DNA. We MAY (and do)
///     assert this single, published-true monotone fact.
///   • Worked example (doc §7.1): CCATT(L)GCTACC, LNA@4, C=1e-4, Na=1, mode None → 63.5276 °C;
///     all-DNA → 59.6923 °C (a +3.84 °C internal-LNA stabilisation, here T-locked).
///
/// CRITICAL THERMODYNAMIC REALITY — what we MUST NOT assert (campaign rule):
///   The McTigue (2004) increments are the published per-step ΔΔH°/ΔΔS° and are MIXED-SIGN — NOT
///   uniformly stabilising (e.g. the G_L-A step has ΔΔH° = +3.162, ΔΔS° = +10.5; several A/G/T-
///   locked steps RAISE entropy more than enthalpy, LOWERING Tm). A single LNA-A — and LNA-G in
///   some steps — can LOWER Tm, and a fully internally-locked A/T-rich duplex can melt BELOW the
///   bare DNA. That is the model behaving CORRECTLY (verbatim from MELTING 5
///   McTigue2004lockedmn.xml). Therefore this file does NOT assert any "adding any LNA never
///   lowers Tm" uniform-stabilisation monotonicity (a FALSE claim, recently removed from the
///   property suite). The ONLY sign-bearing monotone assertion is the LNA-C-does-not-lower-Tm
///   case above, which the published parameters DO guarantee.
///
/// Both methods are pure functions (no iterator), so every probe calls them directly. The very-
/// long-probe probe is pinned with [CancelAfter] so any perf blow-up fails as a timeout, not a
/// hang.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProbeLnaTmFuzzTests
{
    private const PrimerDesigner.SaltCorrectionMode None = PrimerDesigner.SaltCorrectionMode.None;
    private const PrimerDesigner.SaltCorrectionMode SantaLucia = PrimerDesigner.SaltCorrectionMode.SantaLuciaEntropy;
    private const PrimerDesigner.SaltCorrectionMode Owczarzy04 = PrimerDesigner.SaltCorrectionMode.Owczarzy2004Monovalent;
    private const PrimerDesigner.SaltCorrectionMode Owczarzy08 = PrimerDesigner.SaltCorrectionMode.Owczarzy2008Divalent;

    private static readonly PrimerDesigner.SaltCorrectionMode[] AllModes =
        { None, SantaLucia, Owczarzy04, Owczarzy08 };

    // Reference state of the doc §7.1 worked example, so the increment-application
    // assertions are grounded in the published numbers rather than incidental defaults.
    private const double RefConc = 1e-4;
    private const double RefNa = 1.0;

    #region Helpers

    /// <summary>Deterministic ACGT generator — seed fixed LOCALLY so fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    private static double LnaTm(string seq, IEnumerable<int> mask,
        PrimerDesigner.SaltCorrectionMode mode = None) =>
        PrimerDesigner.CalculateMeltingTemperatureNNLna(
            seq, mask.ToArray(), strandConcentrationMolar: RefConc, sodiumMolar: RefNa, saltMode: mode);

    private static double DnaTm(string seq, PrimerDesigner.SaltCorrectionMode mode = None) =>
        PrimerDesigner.CalculateMeltingTemperatureNN(
            seq, strandConcentrationMolar: RefConc, sodiumMolar: RefNa, saltMode: mode);

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROBE-LNATM-001 — LNA-adjusted NN melting temperature : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROBE-LNATM-001 — MC: invalid LNA mask (terminal / out-of-range positions → not computable)

    /// <summary>
    /// MC (invalid mask — TERMINAL position): McTigue (2004) parameterises INTERNAL LNA only, so a
    /// position on the first base (index 0) or the last base (index length−1) is NOT computable
    /// (doc INV-03, ASM-01; source guard <c>pos &lt;= 0 || pos &gt;= seq.Length - 1</c>). The contract is
    /// the null / NaN sentinel, NEVER an IndexOutOfRangeException and NEVER a fabricated finite Tm.
    /// Pinned on both terminal ends, on a clean computable probe, for every salt mode.
    /// </summary>
    [Test]
    public void TerminalLnaPosition_NotComputable_NullAndNaN_NoThrow()
    {
        const string seq = "CCATTGCTACC"; // length 11, doc worked example; internal = 1..9
        foreach (int terminal in new[] { 0, seq.Length - 1 })
        {
            PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { terminal })
                .Should().BeNull($"a TERMINAL LNA at index {terminal} has no McTigue parameter (INV-03)");

            foreach (var mode in AllModes)
            {
                double tm = LnaTm(seq, new[] { terminal }, mode);
                double.IsNaN(tm).Should().BeTrue(
                    $"terminal LNA @{terminal} → NaN sentinel under {mode}, never a fabricated Tm");
            }
        }
    }

    /// <summary>
    /// MC (invalid mask — OUT-OF-RANGE / NEGATIVE position): a position &lt; 0 or ≥ length is plainly
    /// out of range. The same internal-only guard rejects it (any <c>pos &lt;= 0</c> or
    /// <c>pos &gt;= length−1</c>) → not computable. The KEY fuzz assertion: an out-of-range index must
    /// NEVER index past the sequence (no IndexOutOfRangeException from the <c>seq.Substring(i, 2)</c>
    /// step loop or the locked-set lookup) — it returns the null / NaN sentinel instead.
    /// </summary>
    [Test]
    public void OutOfRangeOrNegativeLnaPosition_NotComputable_NoIndexOutOfRange()
    {
        const string seq = "GCGTATACGC"; // length 10, internal = 1..8
        foreach (int bad in new[] { -1, -100, seq.Length, seq.Length + 5, int.MaxValue, int.MinValue })
        {
            Action thermo = () => PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { bad });
            thermo.Should().NotThrow($"out-of-range LNA index {bad} is a documented null sentinel, not a crash");

            PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { bad })
                .Should().BeNull($"out-of-range LNA index {bad} → not computable");

            double tm = 0;
            Action tmCall = () => tm = LnaTm(seq, new[] { bad });
            tmCall.Should().NotThrow($"out-of-range index {bad} must not leak IndexOutOfRange");
            double.IsNaN(tm).Should().BeTrue($"out-of-range LNA index {bad} → NaN sentinel");
        }
    }

    /// <summary>
    /// MC (invalid mask — a VALID internal position mixed with a bad one): ANY single bad position
    /// poisons the whole mask (the source returns null on the first offending index). So a mask that
    /// is mostly internal but contains one terminal / out-of-range index is wholly not computable —
    /// it does NOT silently drop the bad index and compute a partial Tm.
    /// </summary>
    [Test]
    public void MaskWithOneBadPositionAmongValid_WholeResultNotComputable()
    {
        const string seq = "ACGTACGTACGT"; // length 12, internal = 1..10
        var masks = new[]
        {
            new[] { 2, 5, 0 },              // a terminal index
            new[] { 3, 11 },                // last index (length−1)
            new[] { 1, 2, -1 },             // a negative index
            new[] { 4, 99 },                // far out of range
        };
        foreach (var mask in masks)
        {
            PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, mask)
                .Should().BeNull($"mask [{string.Join(",", mask)}] contains a non-internal index → not computable");
            double.IsNaN(LnaTm(seq, mask)).Should().BeTrue(
                $"mask [{string.Join(",", mask)}] → NaN sentinel (no partial Tm)");
        }
    }

    #endregion

    #region PROBE-LNATM-001 — MC/BE: mask longer than the probe (more positions than bases)

    /// <summary>
    /// MC/BE (mask LONGER than the probe): a mask with MORE positions than the probe has bases is
    /// degenerate. Whatever the positions are, the internal-only guard must catch any index that is
    /// ≥ length−1 (a mask longer than the probe necessarily contains such indices, or repeats), so a
    /// "lock every base 0..length−1" mask is not computable (it includes the terminal indices), and a
    /// mask densely packed with out-of-range indices likewise yields the null / NaN sentinel — never
    /// an IndexOutOfRangeException from over-long indexing.
    /// </summary>
    [Test]
    public void MaskLongerThanProbe_NoCrash_NotComputableOrFinite()
    {
        const string seq = "ATGCATGC"; // length 8, internal = 1..6
        // "Lock every base" (length 8 positions) — includes terminals 0 and 7 → not computable.
        var lockAll = Enumerable.Range(0, seq.Length).ToArray();
        PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, lockAll)
            .Should().BeNull("locking every base includes the terminal indices → not computable");
        double.IsNaN(LnaTm(seq, lockAll)).Should().BeTrue("lock-all mask → NaN sentinel");

        // A mask far longer than the probe, filled with out-of-range indices.
        var overlong = Enumerable.Range(0, seq.Length * 4).ToArray(); // 32 indices for an 8-mer
        Action a = () => LnaTm(seq, overlong);
        a.Should().NotThrow("an over-long mask must never IndexOutOfRange");
        double.IsNaN(LnaTm(seq, overlong)).Should().BeTrue("over-long mask hits terminal/out-of-range → NaN");

        // An over-long mask that is over-long ONLY by DUPLICATION of valid INTERNAL indices is
        // tolerated (set semantics, §6.1): more entries than bases, but every entry is internal.
        var dupInternal = new[] { 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6 }; // 12 entries for an 8-mer, all internal
        double tmDup = LnaTm(seq, dupInternal);
        double.IsNaN(tmDup).Should().BeFalse("duplicate INTERNAL indices are tolerated → finite Tm");
        double.IsInfinity(tmDup).Should().BeFalse("duplicate-internal mask Tm must be finite");
        // Duplicates collapse to the set — same Tm as the de-duplicated internal mask.
        tmDup.Should().BeApproximately(LnaTm(seq, new[] { 1, 2, 3, 4, 5, 6 }), 1e-9,
            "duplicate positions collapse to set semantics (§6.1)");
    }

    #endregion

    #region PROBE-LNATM-001 — BE: empty probe, 1-bp probe (NN needs ≥ 2 bases → not computable)

    /// <summary>
    /// BE (empty / 1-bp probe): the NN model sums ΔH°/ΔS° over adjacent dinucleotide STEPS, so it is
    /// undefined for &lt; 2 bases — the base DNA NN returns null, so the LNA result is null / NaN (doc
    /// §3.3, §6.1). A DOCUMENTED sentinel, NOT a throw (even with an empty mask). Pinned for "", null,
    /// and every single base, for every salt mode.
    /// </summary>
    [Test]
    public void EmptyOrSingleBaseProbe_NotComputable_NullAndNaN_NoThrow()
    {
        var degenerate = new[] { "", null, "A", "C", "G", "T" };
        foreach (var seq in degenerate)
        {
            PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq!, Array.Empty<int>())
                .Should().BeNull($"a &lt;2-base probe '{seq ?? "null"}' has no NN step → null");

            foreach (var mode in AllModes)
            {
                double tm = PrimerDesigner.CalculateMeltingTemperatureNNLna(
                    seq!, Array.Empty<int>(), strandConcentrationMolar: RefConc, sodiumMolar: RefNa, saltMode: mode);
                double.IsNaN(tm).Should().BeTrue(
                    $"&lt;2-base probe '{seq ?? "null"}' → NaN sentinel under {mode}, never a throw or fabricated Tm");
            }
        }
    }

    /// <summary>
    /// BE (empty probe with a NON-empty mask): even a populated mask on an empty / sub-2-base probe
    /// must short-circuit to the sentinel — the base NN fails first, so no index is ever taken and no
    /// IndexOutOfRange leaks from the impossible-to-satisfy mask.
    /// </summary>
    [Test]
    public void EmptyProbeWithNonEmptyMask_NotComputable_NoCrash()
    {
        foreach (var seq in new[] { "", "A" })
        {
            Action a = () => LnaTm(seq, new[] { 1, 2, 3 });
            a.Should().NotThrow($"an LNA mask on the degenerate probe '{seq}' must not crash");
            double.IsNaN(LnaTm(seq, new[] { 1, 2, 3 })).Should().BeTrue(
                $"degenerate probe '{seq}' with a mask → NaN sentinel");
            PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 1, 2, 3 })
                .Should().BeNull($"degenerate probe '{seq}' with a mask → null");
        }
    }

    #endregion

    #region PROBE-LNATM-001 — MC: non-DNA / unicode characters in the probe → not computable

    /// <summary>
    /// MC (non-DNA probe): the base DNA NN lookup is over A/C/G/T dinucleotide steps, so a probe with
    /// ANY non-ACGT character (IUPAC 'N', symbols, digits, whitespace, unicode) makes the NN sum fail
    /// → the LNA result is null / NaN (doc §6.1, edge "non-ACGT base → null/NaN"). A DOCUMENTED
    /// sentinel, NOT a crash, on BOTH a clean mask and an internal mask. Tested with the LNA position
    /// kept internal so the ONLY reason for the sentinel is the bad alphabet.
    /// </summary>
    [Test]
    public void NonDnaProbe_NotComputable_NullAndNaN_NoThrow()
    {
        var junk = new[]
        {
            "ACGTNACGT",        // IUPAC 'N'
            "ACGT-ACGT",        // gap
            "ACGT ACGT",        // whitespace
            "ACGT5ACGT",        // digit
            "ACGTéCGT",    // unicode (é)
            "ACGT\0CGT",        // null byte
            "acgt#cgt",         // symbol + lower-case
            "RYSWKMBDHVN",      // all-IUPAC-ambiguity
        };
        foreach (var seq in junk)
        {
            PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 4 })
                .Should().BeNull($"non-ACGT probe '{seq}' has no valid NN step → null");
            foreach (var mode in AllModes)
            {
                double tm = LnaTm(seq, new[] { 4 }, mode);
                double.IsNaN(tm).Should().BeTrue(
                    $"non-ACGT probe '{seq}' → NaN sentinel under {mode}, never a crash or fabricated Tm");
            }
        }

        // Fuzz sweep: junk char injected at every internal position of an otherwise-clean probe.
        const string clean = "ACGTACGTAC";
        for (int i = 1; i < clean.Length - 1; i++)
        {
            string mutated = string.Concat(clean.AsSpan(0, i), "?", clean.AsSpan(i + 1));
            double tm = LnaTm(mutated, new[] { Math.Min(i, clean.Length - 2) });
            double.IsNaN(tm).Should().BeTrue($"junk at index {i} of '{clean}' → NaN, no crash");
        }
    }

    #endregion

    #region PROBE-LNATM-001 — INV: only documented throw is null lnaPositions

    /// <summary>
    /// The ONLY documented throw is a NULL position collection (source ThrowIfNull; doc §3.3
    /// "null lnaPositions → ArgumentNullException"). Pinned on both API entries. Everything else —
    /// bad probe, bad mask — is the null / NaN sentinel, never a throw.
    /// </summary>
    [Test]
    public void NullLnaPositions_Throws_ArgumentNullException()
    {
        Action thermo = () => PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("ACGTACGT", null!);
        thermo.Should().Throw<ArgumentNullException>("a null position collection is the documented throw");

        Action tm = () => PrimerDesigner.CalculateMeltingTemperatureNNLna("ACGTACGT", null!);
        tm.Should().Throw<ArgumentNullException>("a null position collection is the documented throw");
    }

    #endregion

    #region PROBE-LNATM-001 — INV-01: empty mask reproduces the perfect-match NN Tm exactly

    /// <summary>
    /// INV-01 (doc §2.4): an EMPTY LNA-position set adds no increment, so it reproduces the
    /// perfect-match NN Tm EXACTLY. Pinned across random probes and every salt mode against the
    /// non-LNA CalculateMeltingTemperatureNN. (Both finite; bit-for-bit equal.)
    /// </summary>
    [Test]
    public void EmptyMask_ReproducesPerfectMatchNnTm_Exactly()
    {
        for (int seed = 0; seed < 12; seed++)
        {
            string seq = RandomDna(8 + (seed % 14), seed * 31 + 7);
            foreach (var mode in AllModes)
            {
                double lna = LnaTm(seq, Array.Empty<int>(), mode);
                double dna = DnaTm(seq, mode);
                double.IsNaN(lna).Should().BeFalse($"empty-mask LNA Tm of '{seq}' must be finite ({mode})");
                lna.Should().Be(dna, $"empty LNA mask reproduces the perfect-match NN Tm exactly (INV-01, {mode}) for '{seq}'");
            }
        }
    }

    #endregion

    #region PROBE-LNATM-001 — BE: smallest internal-LNA probe, and a valid internal LNA gives a finite Tm

    /// <summary>
    /// BE (smallest internal-LNA-capable probe): a 3-mer has exactly one INTERNAL index (1); 0 and 2
    /// are terminal. So index 1 → finite Tm; 0 and 2 → not computable. This pins the boundary of the
    /// internal/terminal partition at the smallest probe where it is non-trivial.
    /// </summary>
    [Test]
    public void ThreeMer_OnlyMiddleIndexIsInternal()
    {
        const string seq = "GCG"; // length 3
        double internalTm = LnaTm(seq, new[] { 1 });
        double.IsNaN(internalTm).Should().BeFalse("index 1 is the only internal LNA position of a 3-mer → finite Tm");
        double.IsInfinity(internalTm).Should().BeFalse("a valid internal-LNA Tm must be finite");

        foreach (int terminal in new[] { 0, 2 })
            double.IsNaN(LnaTm(seq, new[] { terminal })).Should().BeTrue(
                $"index {terminal} is terminal on a 3-mer → not computable");
    }

    /// <summary>
    /// A valid internal-LNA probe yields a FINITE Tm, and the McTigue increment is ACTUALLY applied:
    /// the LNA Tm DIFFERS from the bare-DNA Tm (the increment is non-zero for these steps). We assert
    /// it DIFFERS — we do NOT claim a sign, because the published increments are mixed-sign and a
    /// single internal LNA can RAISE or LOWER Tm depending on the step (campaign rule: no false
    /// uniform-stabilisation monotonicity). Swept over random internal positions of random probes.
    /// </summary>
    [Test]
    public void ValidInternalLna_FiniteTm_AndIncrementActuallyApplied()
    {
        for (int seed = 0; seed < 20; seed++)
        {
            string seq = RandomDna(10 + (seed % 10), seed * 53 + 3);
            var rng = new Random(seed * 7 + 1);
            int pos = 1 + rng.Next(seq.Length - 2); // a random INTERNAL index
            double lna = LnaTm(seq, new[] { pos });
            double dna = DnaTm(seq);

            double.IsNaN(lna).Should().BeFalse($"internal LNA @{pos} of '{seq}' → finite Tm");
            double.IsInfinity(lna).Should().BeFalse($"internal LNA @{pos} of '{seq}' → finite Tm");
            // Increment IS applied: the McTigue increment is non-zero for every step, so the LNA Tm
            // differs measurably from the bare-DNA Tm. (Sign is intentionally NOT asserted.)
            lna.Should().NotBe(dna,
                $"the McTigue increment must actually shift the Tm vs bare DNA ('{seq}' LNA@{pos})");
        }
    }

    #endregion

    #region PROBE-LNATM-001 — INV (LNA-C): the reliably-stabilising case does NOT lower Tm

    /// <summary>
    /// The ONLY sign-bearing monotone assertion the published parameters allow: LNA-C is reliably
    /// STABILISING — all eight C-locked McTigue steps (both 5'- and 3'-locked CX / XC) are net-
    /// stabilising — so substituting an internal LNA-C does NOT LOWER the Tm vs the all-DNA duplex.
    /// We pin <c>LNA-C Tm ≥ bare-DNA Tm</c> on probes whose chosen internal position is a 'C'. This is
    /// deliberately the LNA-C case ONLY: we do NOT generalise to LNA-A/G/T, whose mixed-sign
    /// increments can LOWER Tm (the campaign rule). The doc §7.1 worked example is a T-locked
    /// stabilisation (+3.84 °C) and is also pinned here as a numeric ground truth.
    /// </summary>
    [Test]
    public void InternalLnaC_DoesNotLowerTm_AndWorkedExample()
    {
        // 1) The §7.1 worked example: CCATT(L)GCTACC, LNA@4 (a 'T'), C=1e-4, Na=1, mode None.
        const string worked = "CCATTGCTACC";
        double lnaWorked = LnaTm(worked, new[] { 4 });
        double dnaWorked = DnaTm(worked);
        lnaWorked.Should().BeApproximately(63.527594, 1e-3, "doc §7.1 LNA-adjusted Tm");
        dnaWorked.Should().BeApproximately(59.692264, 1e-3, "doc §7.1 all-DNA Tm");
        (lnaWorked - dnaWorked).Should().BeApproximately(63.527594 - 59.692264, 1e-3,
            "doc §7.1: the single internal LNA raises Tm by +3.84 °C");

        // 2) LNA-C reliable stabilisation: for several probes, lock an internal 'C' and assert the
        //    Tm does NOT drop below bare DNA.
        for (int seed = 0; seed < 25; seed++)
        {
            string seq = RandomDna(10 + (seed % 12), seed * 41 + 9);
            // Find an internal index that holds a 'C'.
            int cPos = -1;
            for (int i = 1; i < seq.Length - 1; i++)
                if (seq[i] == 'C') { cPos = i; break; }
            if (cPos < 0) continue; // no internal C in this random probe; skip

            double lna = LnaTm(seq, new[] { cPos });
            double dna = DnaTm(seq);
            double.IsNaN(lna).Should().BeFalse($"internal LNA-C @{cPos} of '{seq}' → finite");
            lna.Should().BeGreaterThanOrEqualTo(dna - 1e-9,
                $"an internal LNA-C does NOT lower Tm (all 8 C-locked McTigue steps stabilise) — '{seq}' C@{cPos}");
        }
    }

    #endregion

    #region PROBE-LNATM-001 — BE: very long probe (no overflow / non-finite / hang)

    /// <summary>
    /// BE (very long probe): a long probe with several internal LNA positions must produce a FINITE
    /// Tm promptly — no overflow, no Inf/NaN leak, no perf blow-up. Pinned with [CancelAfter] so a
    /// hang fails as a timeout rather than wedging the suite.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void VeryLongProbe_FiniteTm_NoHang()
    {
        string seq = RandomDna(500, 20250626);
        var mask = Enumerable.Range(1, seq.Length - 2).Where(i => i % 37 == 0).ToArray(); // sparse internal LNAs
        foreach (var mode in AllModes)
        {
            double tm = LnaTm(seq, mask, mode);
            double.IsNaN(tm).Should().BeFalse($"long-probe LNA Tm must be finite ({mode})");
            double.IsInfinity(tm).Should().BeFalse($"long-probe LNA Tm must be finite ({mode})");
        }
    }

    #endregion
}
