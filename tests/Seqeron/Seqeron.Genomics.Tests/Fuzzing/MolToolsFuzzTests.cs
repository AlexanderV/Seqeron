using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the MolTools area — CRISPR PAM (protospacer adjacent motif)
/// site finding (CRISPR-PAM-001), CRISPR guide RNA design (CRISPR-GUIDE-001),
/// CRISPR off-target analysis (CRISPR-OFF-001), and primer melting-temperature
/// calculation (PRIMER-TM-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds malformed, out-of-domain and boundary inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output (spurious or out-of-bounds
/// sites), and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, ArgumentOutOfRangeException leaking from internal
/// indexing). Every input must resolve to EITHER a well-defined, theory-correct
/// result, OR a *documented, intentional* validation exception
/// (ArgumentException / ArgumentNullException). A raw runtime exception, a hang,
/// or a blow-up of bogus "PAM sites" on garbage input is a bug, not a passing
/// test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CRISPR-PAM-001 — PAM site finding
/// Checklist: docs/checklists/03_FUZZING.md, row 18.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — an invalid CRISPR system selector (an undefined
///          CrisprSystemType enum value, i.e. a "PAM specification" the API does
///          not know) and non-DNA junk characters in the sequence.
///   • BE = Boundary Exploitation — the empty sequence and the degenerate
///          "PAM/guide span longer than the sequence" case (no scan window fits).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 18 targets:
///   "Invalid PAM sequences, non-DNA characters, empty seq, PAM longer than seq".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The PAM-site-finding contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A PAM is a short motif (e.g. SpCas9 NGG, where N = any base) immediately
/// adjacent to a CRISPR protospacer; Cas nucleases require it to bind and cleave
/// (PAM_Site_Detection.md §2.1). The detector scans BOTH strands for the
/// system-specific PAM motif (specified with IUPAC ambiguity codes), extracts the
/// guide-length target on the side the system dictates (PAM-after-target for Cas9,
/// PAM-before-target for Cas12a), and yields a PamSite ONLY when that target fits
/// within sequence bounds (PAM_Site_Detection.md §2.2, §2.4 INV-01..INV-03).
///
/// API entry: CrisprDesigner.FindPamSites(...)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs lines 40–60)
/// in two overloads — typed DnaSequence and raw string:
///   • FindPamSites(DnaSequence, CrisprSystemType = SpCas9): null → ArgumentNullException
///     (ThrowIfNull, line 44). The DnaSequence itself is a *validated* type whose
///     constructor rejects any non-A/C/G/T character with ArgumentException, so a
///     non-DNA sequence can never even reach this overload.
///   • FindPamSites(string, CrisprSystemType = SpCas9): null OR empty → the empty
///     result via `yield break` (lines 55–56); otherwise the input is upper-cased
///     and scanned RAW (no A/C/G/T validation). This is the surface that can be fed
///     malformed content directly, so it is the primary MC target here.
///
/// IMPORTANT — what "invalid PAM" means for THIS API (the central design fact the
/// checklist row probes): the PAM motif is NOT a free-form caller string. It is
/// selected by the CrisprSystemType enum and resolved by GetSystem(...) into one of
/// a fixed set of IUPAC-valid motifs (NGG, NAG, NNGRRT, TTTV, TTCN — §4.2). So an
/// "invalid PAM specification" is an UNDEFINED enum value, which GetSystem rejects
/// with ArgumentException (CrisprDesigner.cs line 27). Because every built-in PAM
/// motif contains only valid IUPAC codes, the IUPAC matcher's own
/// ArgumentOutOfRangeException guard (IupacHelper.MatchesIupac, the `_ => throw`
/// arm) is never reached from the PAM side — the motif is always well-formed; it is
/// the *sequence* base that is tested against the (valid) PAM code. We pin both
/// facets: an undefined system throws, and the IUPAC matcher rejects truly invalid
/// codes (so the validity of the fixed motifs is load-bearing, not incidental).
///
/// Non-DNA characters in the SEQUENCE (MC) — pinned on both surfaces:
///   • Typed surface: rejected up front by the DnaSequence constructor
///     (ArgumentException, "Invalid nucleotide") — junk never reaches the scanner.
///   • Raw-string surface: tolerated, NOT crashing. A junk character is tested as
///     the `nucleotide` argument of MatchesIupac(seqChar, pamChar); it matches no
///     IUPAC code, so it can never satisfy a PAM position and never invents a site.
///     The reverse-strand pass complements via GetComplementBase, whose final arm
///     passes any non-IUPAC char THROUGH unchanged (SequenceExtensions.cs line 156),
///     so no exception and no out-of-range indexing — at worst, fewer or no sites.
///     We pin that junk in the sequence yields no bogus sites and no crash, and that
///     interleaving junk into an otherwise-valid PAM window suppresses that site.
///
/// Empty sequence (BE): the raw-string overload short-circuits null/empty to the
/// empty enumerable (yield break); the typed overload over an empty DnaSequence has
/// a forward-scan bound `i ≤ len − pamLen` that is negative, so the loop never runs.
/// No PAM sites, no division, no indexing past the end (PAM_Site_Detection.md §6.1).
///
/// PAM(+guide span) longer than the sequence (BE): when the sequence is shorter than
/// the PAM pattern, the forward and reverse scan bounds `i ≤ len − pamLen` are
/// negative, so neither loop body runs and no Substring is taken past the end. Even
/// when the PAM motif itself fits but the guide-length target cannot (sequence
/// shorter than PAM+guide), the in-loop `targetStart >= 0 && targetEnd < seq.Length`
/// guard suppresses the site (INV-02). Either way: an empty result, never an
/// IndexOutOfRangeException (PAM_Site_Detection.md §6.1, "Sequence shorter than PAM
/// plus guide → no valid sites").
///
/// Documented invariants pinned on every positive site (PAM_Site_Detection.md §2.4):
/// INV-01 the matched PAM satisfies the system motif under IUPAC matching;
/// INV-02 the extracted target fits within sequence bounds; INV-03 both strands are
/// searched. The raw-string overload's scan body lives inside a `yield` iterator, so
/// every test forces enumeration (`.ToList()`) — the documented short-circuit and
/// any hang surface only on enumeration, and any hang would manifest as a
/// non-terminating materialization.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CRISPR-GUIDE-001 — guide RNA design
/// Checklist: docs/checklists/03_FUZZING.md, row 19.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — non-DNA junk in the guide string fed to the raw
///          standalone evaluator, and (the dual surface) non-DNA fed to the typed
///          region-design path, which is rejected at DnaSequence construction.
///   • BE = Boundary Exploitation — a sequence with NO PAM sites (no guides), a
///          sequence too short for the system's guide to fit (guide-span > seq),
///          and the degenerate zero-length guide (empty string).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 19 targets:
///   "No PAM sites in seq, guide length > seq, guide length 0, non-DNA".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The guide-RNA-design contract under test (Guide_RNA_Design.md)
/// ───────────────────────────────────────────────────────────────────────────
/// Guide design has TWO public surfaces (Guide_RNA_Design.md §5.1):
///   • CrisprDesigner.DesignGuideRnas(DnaSequence, regionStart, regionEnd, system,
///     parameters) — extracts PAM-adjacent candidate guides inside a target region,
///     scores each, and YIELDS only candidates whose Score >= MinScore (INV-02).
///     It is a typed surface: the DnaSequence constructor already rejects non-DNA,
///     so guide design only ever scans clean A/C/G/T. It null-guards the sequence
///     (ArgumentNullException) and range-guards the region
///     (ArgumentOutOfRangeException when regionStart/End fall outside the sequence)
///     — §3.3, §6.1. The scan reuses FindPamSitesCore, so "no PAM sites" and
///     "guide-span longer than the sequence" inherit the PAM scanner's bounds
///     discipline: no NGG ⇒ no candidates; a sequence too short to hold the
///     protospacer ⇒ the in-loop `targetStart >= 0 && targetEnd < seq.Length` guard
///     suppresses every site ⇒ an empty result, never an IndexOutOfRangeException.
///   • CrisprDesigner.EvaluateGuideRna(string guideSequence, system, parameters) —
///     scores a standalone guide string from its composition (GC, seed GC, poly-T,
///     self-complementarity, restriction sites; §2.2). It does NOT enforce a guide
///     length (§3.3, §6.1: non-20-nt guides are accepted and scored naturally) but
///     it DOES reject the degenerate empty/null guide with ArgumentNullException
///     (source guard, CrisprDesigner.cs line 205). It scores the raw string without
///     A/C/G/T validation, so non-DNA is the MC target for THIS surface.
///
/// THE FOUR ROW-19 FUZZ TARGETS, mapped to the theory-correct contract:
///   • No PAM sites in seq (BE): an all-A region has no NGG, so DesignGuideRnas
///     yields the empty result — never a crash, never a fabricated guide.
///   • Guide length > seq (BE): a sequence shorter than the 20-nt protospacer (even
///     when an NGG fits) yields no candidate — the target would start at a negative
///     index and the bounds guard drops it (mirrors CRISPR-PAM-001's
///     "PAM fits but guide does not"). Empty result, no IndexOutOfRangeException.
///   • Guide length 0 (BE, KEY): a zero-length (empty) guide is meaningless and is
///     the div-by-zero hazard — CalculateSelfComplementarity divides by
///     `length * length` and CalculateGcContent would face an empty span. The source
///     short-circuits BOTH: the empty guard throws ArgumentNullException BEFORE any
///     scoring, and CalculateGcContent null/empty-guards to 0. So guideLen 0 is a
///     DOCUMENTED validation throw, not a DivideByZero/IndexOutOfRange crash.
///   • Non-DNA (MC): the typed design surface rejects it at DnaSequence construction
///     (ArgumentException); the raw evaluator TOLERATES it — junk is excluded from GC
///     (CalculateGcFraction counts only A/C/G/T/U, returns 0 on no-valid), and the
///     reverse-complement used by self-complementarity passes non-IUPAC chars through
///     unchanged (GetComplementBase fall-through). So a junk guide scores without
///     crashing and never divides by zero (length > 0).
///
/// Documented invariant pinned on every produced candidate (Guide_RNA_Design.md
/// §2.4): INV-01 the score is clamped to `>= 0` (Math.Max(0, score)); INV-02
/// DesignGuideRnas yields only Score >= MinScore. DesignGuideRnas is a yield
/// iterator, so every test forces enumeration (`.ToList()`); the documented
/// short-circuit and any hang surface only on enumeration.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CRISPR-OFF-001 — off-target analysis
/// Checklist: docs/checklists/03_FUZZING.md, row 20.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the zero mismatch tolerance (`maxMismatches = 0`,
///          the lower end of the documented 0..5 range) and the degenerate empty
///          guide (zero-length guide string).
///   • MC = Malformed Content — a guide of all N's (the IUPAC "any base" wildcard
///          fed as a literal guide), probing whether 'N' is a wildcard or an
///          ordinary non-A/C/G/T character in off-target matching.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 20 targets:
///   "Zero mismatch tolerance, empty guide, guide of all N's".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The off-target-analysis contract under test (Off_Target_Analysis.md)
/// ───────────────────────────────────────────────────────────────────────────
/// Off-target analysis scans a genome for PAM-supported, guide-length targets that
/// differ from the guide by a bounded, NON-ZERO number of mismatches — exact matches
/// are deliberately excluded so the on-target (and any perfect duplicate) is never
/// reported as an off-target (Off_Target_Analysis.md §2.2, INV-01). The mismatch
/// count is Hamming-style: CountMismatches compares the guide and the PAM-adjacent
/// target POSITION-BY-POSITION with a plain `!=` (CrisprDesigner.cs lines 393–405),
/// so there is NO IUPAC/wildcard semantics — 'N' is an ordinary character that
/// differs from A/C/G/T. A site is yielded only when `0 < mismatches <= maxMismatches`
/// (line 359, INV-02).
///
/// API entry: CrisprDesigner.FindOffTargets(string guide, DnaSequence genome,
///   int maxMismatches = 3, CrisprSystemType = SpCas9)
///   (CrisprDesigner.cs lines 331–373) — a yield iterator. Its documented guards
///   (Off_Target_Analysis.md §3.3, §6.1):
///   • null/empty guide → ArgumentNullException (line 337–338);
///   • null genome → ArgumentNullException (line 339);
///   • maxMismatches < 0 or > 5 → ArgumentOutOfRangeException (line 340–341);
///   • guide length != system.GuideLength → ArgumentException (line 345–348).
///   Because the body is a `yield` iterator, every guard fires only on enumeration,
///   so each test forces materialization (`.ToList()`).
///
/// THE THREE ROW-20 FUZZ TARGETS, mapped to the theory-correct contract:
///   • Zero mismatch tolerance (BE, KEY): `maxMismatches = 0` is in-range (the guard
///     rejects only < 0), so it does NOT throw. But the yield condition
///     `mismatches > 0 && mismatches <= 0` is unsatisfiable, so NO site is ever
///     emitted — not even an exact on-target (which is excluded by `mismatches > 0`
///     anyway, INV-01). Theory-correct result: the EMPTY enumerable. A genome that
///     contains a perfect copy of the guide therefore yields zero off-targets at
///     tolerance 0, and the same genome with a 1-mismatch site yields that site only
///     at tolerance >= 1 — never a crash, never a div-by-zero.
///   • Empty guide (BE): a zero-length (or null) guide is degenerate and is rejected
///     up front with ArgumentNullException (line 337–338) — BEFORE any length check,
///     PAM scan, or CountMismatches loop, so there is no IndexOutOfRange and no
///     division by a zero seed/guide length. Pinned for both "" and null, on
///     enumeration. (An empty-but-length-matching guide is impossible: GuideLength is
///     always >= 20, so "" can never pass the length check even if it reached it.)
///   • Guide of all N's (MC, KEY): 'N' is the central ambiguity question. This API
///     does NOT treat 'N' as a wildcard — CountMismatches is a literal `!=` compare.
///     A 20-nt all-N guide PASSES the SpCas9 length check (length 20 == GuideLength),
///     then mismatches EVERY base of any A/C/G/T target, so the per-site mismatch
///     count is the full guide length (20) — far above the 0..5 cap — and NO site is
///     ever yielded. So all-N does NOT "match everything" and does NOT blow up into
///     an off-target at every position (no hang): it is the OPPOSITE — it matches
///     nothing, because every position is a guaranteed mismatch. Theory-correct
///     result: the EMPTY enumerable, promptly, with no exception. (A wrong-length
///     all-N guide is instead rejected with ArgumentException by the length check,
///     pinned separately.)
///
/// Documented invariants pinned on every produced off-target (Off_Target_Analysis.md
/// §2.4): INV-01 only `mismatches > 0` sites are returned (exact matches excluded);
/// INV-02 returned mismatch counts are bounded by the requested `maxMismatches`.
/// FindOffTargets is a yield iterator, so every test forces enumeration (`.ToList()`).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PRIMER-TM-001 — primer melting-temperature calculation
/// Checklist: docs/checklists/03_FUZZING.md, row 21.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the empty sequence, the single base (1-bp), and a
///          long 100+-bp primer (upper size boundary; overflow / non-finite hazard).
///   • INJ = Injection — non-DNA / special / unicode characters in the primer string,
///          an all-N primer (no A/C/G/T at all), and the `null` reference.
/// — docs/checklists/03_FUZZING.md §Description (BE; INJ = injection of special chars /
///   null bytes / unicode); row 21 targets:
///   "Empty seq, 1-bp, 100+ bp, all-N, non-DNA chars, null".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The melting-temperature contract under test (Melting_Temperature.md, Primer_Design.md)
/// ───────────────────────────────────────────────────────────────────────────
/// Tm is the temperature at which half of a DNA duplex is dissociated. The MolTools
/// primer surface PrimerDesigner.CalculateMeltingTemperature(string) (PRIMER-TM-001;
/// Primer_Design.md §7.2 names it the PRIMER-DESIGN-001 prerequisite) selects ONE of two
/// empirical composition models by the count of VALID bases (PrimerDesigner.cs 197–219):
///   • Wallace rule for SHORT oligos with `validLength < 14`:
///         Tm = 2·(A+T) + 4·(G+C)            — Melting_Temperature.md §2.2 INV-01.
///   • Marmur-Doty GC formula for `validLength >= 14`:
///         Tm = max(0, 64.9 + 41·(GC − 16.4)/validLength)
///                                            — Melting_Temperature.md §2.2 INV-02.
/// The constants (2, 4; 64.9, 41, 16.4; the 14-nt threshold) live in
/// ThermoConstants (WallaceAtContribution, WallaceGcContribution, WallaceMaxLength,
/// MarmurDotyBase/GcCoefficient/GcOffset). Output is °C, a finite double; the method
/// throws NO exception for any string input — it GUARDS every degenerate case by
/// returning 0 (Melting_Temperature.md §3.3, §6.1). The KEY hazards the row probes are
/// the division by `validLength` in Marmur-Doty (a zero-valid-base input would divide by
/// zero → NaN) and overflow / non-finite leakage on a long primer.
///
/// THE SIX ROW-21 FUZZ TARGETS, mapped to the theory-correct contract:
///   • Empty seq (BE): `string.IsNullOrEmpty` short-circuits to Tm = 0 BEFORE any count
///     or division (line 199–200) — a DEFINED degenerate boundary, never a crash or NaN
///     (Melting_Temperature.md §6.1 "Empty / length-1 → 0").
///   • 1-bp (BE): a single base has validLength 1 < 14, so the Wallace rule applies and
///     yields a TINY, EXACT Tm: "A"/"T" → 2 °C, "G"/"C" → 4 °C (INV-01). No division,
///     no NaN — the lower non-empty size boundary is well-defined.
///   • 100+ bp (BE): a long primer has validLength >= 14, so Marmur-Doty applies; the
///     result is a LARGE but FINITE Tm (no overflow, no Inf/NaN). A 100-nt all-GC primer
///     pins the exact Marmur-Doty value 64.9 + 41·(100 − 16.4)/100 = 99.176 °C, proving
///     the long-primer path is finite and formula-correct, not overflowing.
///   • All-N (INJ/KEY div-by-zero hazard): 'N' is neither A/T nor G/C, so for an all-N
///     primer at + gc = validLength = 0. The source GUARDS this with an explicit
///     `validLength == 0 → return 0` (line 208–209) BEFORE the Marmur-Doty division, so
///     the denominator is never zero — the theory-correct result is 0 °C, NOT a NaN /
///     DivideByZero. This is the central INJ probe: a model whose denominator is the
///     valid-base count must not divide by zero when no base is valid.
///   • Non-DNA chars (INJ): special characters, digits, spaces, and unicode are simply
///     not counted as A/T or G/C — they are IGNORED (the method's own contract: "Only
///     standard DNA bases (A, C, G, T) are recognized; all other characters are
///     ignored", PrimerDesigner.cs 194–195). So junk never crashes and never inflates
///     Tm; a primer of pure junk has validLength 0 → Tm 0, and junk interleaved into a
///     real primer contributes nothing to the count (the Tm equals that of the cleaned
///     primer). Lower-case is upper-cased first (line 202), so case is irrelevant.
///   • Null (INJ): a null reference is caught by `string.IsNullOrEmpty(primer)` and
///     returns 0 — an ArgumentNullException is NOT thrown here, and crucially there is
///     NEVER a NullReferenceException (the guard runs before any dereference). The
///     salt-corrected sibling CalculateMeltingTemperatureWithSalt(string, double) shares
///     the same null/empty short-circuit (line 229–230), so null is safe on both.
///
/// Documented invariants pinned (Melting_Temperature.md §2.4): INV-01 Wallace
/// Tm = 2·(A+T) + 4·(G+C); INV-02 Marmur-Doty Tm = 64.9 + 41·(GC − 16.4)/N; INV-05
/// empty / length-1 (and, here, all-non-DNA) input → 0. CalculateMeltingTemperature is
/// a pure function (no iterator), so every probe calls it directly. Every test asserts
/// the result is FINITE (not NaN, not ±Infinity) in addition to its exact value.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MolToolsFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
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
    //  CRISPR-PAM-001 — PAM site finding : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CRISPR-PAM-001 — PAM site finding

    #region MC — Invalid PAM specification (undefined CRISPR system)

    /// <summary>
    /// MC: an "invalid PAM" for this API is an undefined CRISPR system selector —
    /// the PAM motif is chosen by the CrisprSystemType enum, not a free-form string
    /// (CrisprDesigner.cs §PAM Definitions). An out-of-range enum value resolves to
    /// the `_ => throw new ArgumentException` arm of GetSystem (line 27), so the
    /// detector REJECTS the bad specification rather than scanning with a garbage
    /// motif or crashing. Pinned on both surfaces; the raw-string overload resolves
    /// the system inside its iterator, so we force enumeration to surface the throw.
    /// </summary>
    [Test]
    public void FindPamSites_UndefinedCrisprSystem_ThrowsArgumentException()
    {
        const CrisprSystemType bogus = (CrisprSystemType)9999;

        var typed = () => CrisprDesigner.FindPamSites(new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG"), bogus).ToList();
        var raw = () => CrisprDesigner.FindPamSites("ACGTACGTACGTACGTACGTACGTAGG", bogus).ToList();

        typed.Should().Throw<ArgumentException>(
            "an undefined CrisprSystemType is an invalid PAM specification; GetSystem rejects it rather than scanning a garbage motif");
        raw.Should().Throw<ArgumentException>(
            "the raw-string overload resolves the same GetSystem mapping and rejects the undefined system on enumeration");
    }

    /// <summary>
    /// MC: the fixed PAM motifs (NGG, NAG, NNGRRT, TTTV, TTCN) are load-bearing —
    /// they are valid IUPAC patterns, which is WHY scanning never trips the IUPAC
    /// matcher's own validity guard. This pins that guard directly: a genuinely
    /// invalid IUPAC code IS rejected by IupacHelper.MatchesIupac (the `_ => throw`
    /// arm), proving the no-crash guarantee for PAM scanning rests on the motifs
    /// being well-formed, not on the matcher silently swallowing junk.
    /// </summary>
    [Test]
    public void MatchesIupac_InvalidPamCode_ThrowsArgumentOutOfRange()
    {
        var act = () => IupacHelper.MatchesIupac('A', 'Z');

        act.Should().Throw<ArgumentOutOfRangeException>(
            "'Z' is not one of the 15 IUPAC codes; the matcher rejects an invalid PAM code rather than matching arbitrarily");
    }

    #endregion

    #region MC — Non-DNA characters in the sequence

    /// <summary>
    /// MC: non-DNA junk in the sequence fed to the TYPED overload is rejected up
    /// front — the DnaSequence constructor validates A/C/G/T and throws
    /// ArgumentException on the first offending character (DnaSequence.cs
    /// ValidateSequence), so junk never reaches the PAM scanner at all.
    /// </summary>
    [Test]
    public void FindPamSites_NonDnaSequence_TypedOverload_RejectedAtConstruction()
    {
        var act = () => new DnaSequence("ACGT$#@!XYZ123");

        act.Should().Throw<ArgumentException>(
            "the validated DnaSequence type rejects any non-A/C/G/T character, so the typed PAM scanner only ever sees clean DNA");
    }

    /// <summary>
    /// MC: non-DNA junk fed to the RAW-STRING overload must NOT crash and must NOT
    /// invent sites. Each sequence character is tested as the `nucleotide` argument
    /// of MatchesIupac(seqChar, pamChar); a junk char matches no IUPAC code, so it
    /// can never satisfy a PAM position. The reverse-strand pass complements via
    /// GetComplementBase, whose fall-through arm passes non-IUPAC chars through
    /// unchanged — no exception, no out-of-range indexing. Pure-junk input yields no
    /// PAM sites. We force enumeration so the in-iterator scan actually runs.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPamSites_NonDnaSequence_RawOverload_NoSitesNoCrash()
    {
        var act = () => CrisprDesigner.FindPamSites("$#@!XYZ123 ￿qwerty", CrisprSystemType.SpCas9).ToList();

        act.Should().NotThrow(
            "non-DNA characters are tested against the PAM via IUPAC matching and complemented via a pass-through arm; neither path indexes out of range");

        CrisprDesigner.FindPamSites("$#@!XYZ123 ￿qwerty", CrisprSystemType.SpCas9)
            .Should().BeEmpty(
                "a character that is not A/C/G/T matches no PAM position, so pure junk can never produce a spurious PAM site");
    }

    /// <summary>
    /// MC: junk INTERLEAVED into an otherwise-valid PAM window must suppress that
    /// site, not crash. "…AGG" at the right offset is a real forward SpCas9 site;
    /// replacing the final 'G' with '#' breaks the NGG match (the '#' satisfies no
    /// IUPAC code), so no forward site is reported there — proving junk neither
    /// fabricates nor silently "rounds to" a valid PAM.
    /// </summary>
    [Test]
    public void FindPamSites_JunkInsidePamWindow_SuppressesThatSite()
    {
        // 20 A's (guide) + "AG#" — the NGG window at index 20 is broken by '#'.
        string broken = new string('A', 20) + "AG#";

        var sites = CrisprDesigner.FindPamSites(broken, CrisprSystemType.SpCas9).ToList();

        sites.Should().NotContain(s => s.IsForwardStrand && s.Position == 20,
            "the '#' in the third PAM position fails IUPAC matching, so the forward NGG site at index 20 is correctly NOT reported");
        sites.Should().OnlyContain(s => s.PamSequence.All(c => c == 'A' || c == 'C' || c == 'G' || c == 'T'),
            "any site that IS reported still has a clean A/C/G/T PAM string — junk never leaks into a reported PAM");
    }

    #endregion

    #region BE — Empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The raw-string overload
    /// short-circuits null/empty to the empty enumerable (yield break, lines 55–56);
    /// the typed overload over an empty DnaSequence has a negative forward-scan bound
    /// `i ≤ 0 − pamLen`, so the loop never runs. Neither path divides, indexes past
    /// the end, or hangs (PAM_Site_Detection.md §6.1). Pinned for both surfaces plus
    /// the raw null input.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPamSites_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var typed = () => CrisprDesigner.FindPamSites(new DnaSequence(string.Empty), CrisprSystemType.SpCas9).ToList();
        var rawEmpty = () => CrisprDesigner.FindPamSites(string.Empty, CrisprSystemType.SpCas9).ToList();
        var rawNull = () => CrisprDesigner.FindPamSites((string)null!, CrisprSystemType.SpCas9).ToList();

        typed.Should().NotThrow("an empty sequence has no scan window; the forward bound is negative so the loop never runs");
        rawEmpty.Should().NotThrow("the raw-string overload short-circuits empty input to an empty result");
        rawNull.Should().NotThrow("the raw-string overload treats null input as empty, not as an error");

        CrisprDesigner.FindPamSites(new DnaSequence(string.Empty), CrisprSystemType.SpCas9).Should().BeEmpty();
        CrisprDesigner.FindPamSites(string.Empty, CrisprSystemType.SpCas9).Should().BeEmpty();
        CrisprDesigner.FindPamSites((string)null!, CrisprSystemType.SpCas9).Should().BeEmpty();
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no typed input". The typed
    /// overload guards it with ArgumentNullException (ThrowIfNull, line 44), raised
    /// eagerly at the call — never a NullReferenceException dereferencing
    /// `sequence.Sequence` (PAM_Site_Detection.md §6.1).
    /// </summary>
    [Test]
    public void FindPamSites_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => CrisprDesigner.FindPamSites((DnaSequence)null!, CrisprSystemType.SpCas9);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload null-guards its sequence; null is rejected, never dereferenced into a NullReferenceException");
    }

    #endregion

    #region BE — PAM (and guide span) longer than the sequence

    /// <summary>
    /// BE: a sequence SHORTER than the PAM pattern itself is the degenerate "PAM
    /// longer than seq" case. SaCas9's PAM is 6 nt (NNGRRT); on a 3-nt sequence the
    /// scan bounds `i ≤ len − pamLen` are negative on both strands, so neither loop
    /// runs and no Substring is taken past the end — an empty result, never an
    /// IndexOutOfRangeException (PAM_Site_Detection.md §6.1). Pinned on both surfaces.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPamSites_PamLongerThanSequence_IsEmptyAndDoesNotThrow()
    {
        var typed = () => CrisprDesigner.FindPamSites(new DnaSequence("ACG"), CrisprSystemType.SaCas9).ToList();
        var raw = () => CrisprDesigner.FindPamSites("ACG", CrisprSystemType.SaCas9).ToList();

        typed.Should().NotThrow(
            "the 6-nt NNGRRT PAM cannot fit in a 3-nt sequence; the scan bound is negative so no Substring is taken past the end");
        raw.Should().NotThrow("the raw-string overload is equally guarded against indexing past the sequence end");

        CrisprDesigner.FindPamSites(new DnaSequence("ACG"), CrisprSystemType.SaCas9).Should().BeEmpty(
            "a PAM longer than the whole sequence yields no sites, not a crash");
        CrisprDesigner.FindPamSites("ACG", CrisprSystemType.SaCas9).Should().BeEmpty();
    }

    /// <summary>
    /// BE: the PAM motif FITS but the guide-length target does not — the sequence is
    /// long enough to hold an NGG but too short for the 20-nt protospacer upstream of
    /// it. The in-loop bounds guard `targetStart >= 0 && targetEnd < seq.Length`
    /// (line 87) suppresses the site (INV-02) rather than slicing a negative-start
    /// Substring. "AGG" alone matches NGG at index 0 but targetStart = 0 − 20 = −20,
    /// so no site is yielded. Pinned: no crash, empty result.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPamSites_PamFitsButGuideDoesNot_IsEmptyAndDoesNotThrow()
    {
        var typed = () => CrisprDesigner.FindPamSites(new DnaSequence("AGG"), CrisprSystemType.SpCas9).ToList();
        var raw = () => CrisprDesigner.FindPamSites("AGG", CrisprSystemType.SpCas9).ToList();

        typed.Should().NotThrow(
            "the NGG matches but the 20-nt target would start at -20; the bounds guard suppresses it instead of slicing past the start");
        raw.Should().NotThrow();

        CrisprDesigner.FindPamSites(new DnaSequence("AGG"), CrisprSystemType.SpCas9).Should().BeEmpty(
            "no protospacer fits upstream of an NGG at index 0, so the site is correctly not yielded");
        CrisprDesigner.FindPamSites("AGG", CrisprSystemType.SpCas9).Should().BeEmpty();
    }

    #endregion

    #region Positive sanity — known PAM sites are found at the correct positions

    /// <summary>
    /// Positive sanity: alongside the degenerate probes, a textbook SpCas9 forward
    /// site must be detected at the CORRECT position with the CORRECT PAM and target,
    /// so the boundary hardening never silently breaks the core function. The
    /// sequence is 20 A's (the protospacer) followed by "AGG": the only GG pair sits
    /// at indices 21–22, so NGG matches at PAM position 20 (N = seq[20] = 'A'). The
    /// forward site therefore has Position 20, PamSequence "AGG", and a 20-nt all-A
    /// target starting at index 0 (INV-01, INV-02).
    /// </summary>
    [Test]
    public void FindPamSites_KnownForwardNggSite_DetectedAtCorrectPosition()
    {
        string seq = new string('A', 20) + "AGG"; // length 23

        var forward = CrisprDesigner.FindPamSites(seq, CrisprSystemType.SpCas9)
            .Where(s => s.IsForwardStrand)
            .ToList();

        var site = forward.Should().ContainSingle(s => s.Position == 20).Subject;
        site.PamSequence.Should().Be("AGG", "the NGG window at index 20 is 'AGG' (N = the upstream 'A')");
        site.PamSequence[1].Should().Be('G');
        site.PamSequence[2].Should().Be('G', "INV-01: the matched PAM satisfies the NGG motif under IUPAC matching");
        site.TargetSequence.Should().Be(new string('A', 20), "the 20-nt protospacer immediately upstream of the PAM is all-A");
        site.TargetStart.Should().Be(0, "INV-02: the target fits within bounds, starting at the first base");
        site.System.PamSequence.Should().Be("NGG");
    }

    /// <summary>
    /// Positive sanity / INV-03: both strands are searched. A forward NGG on one
    /// strand appears as a reverse-complement CCN motif on the other; a construct
    /// carrying a valid reverse-strand site must therefore produce at least one site
    /// with IsForwardStrand == false. "CCT" + 20 A's contains "CC" at the 5' end,
    /// which is the reverse complement of an "GG" PAM read on the opposite strand, so
    /// the reverse-strand scan must report a non-forward site. Pinned: the detector
    /// does not silently search only the forward strand.
    /// </summary>
    [Test]
    public void FindPamSites_ReverseStrandSite_IsAlsoSearched()
    {
        // Reverse complement of ("CCT" + 20×A) is (20×T + "AGG"): a forward NGG on the
        // reverse-complement strand → a reverse-strand site on the original.
        string seq = "CCT" + new string('A', 20);

        var sites = CrisprDesigner.FindPamSites(seq, CrisprSystemType.SpCas9).ToList();

        sites.Should().Contain(s => !s.IsForwardStrand,
            "INV-03: the reverse strand is scanned, so a reverse-complement (CCN) PAM yields a non-forward site");
        sites.Should().OnlyContain(s =>
                s.TargetStart >= 0 &&
                s.PamSequence.Length == 3,
            "every yielded site is well-formed: a 3-nt PAM and a non-negative target start");
    }

    /// <summary>
    /// Positive sanity: a Cas12a (Cpf1) site exercises the PAM-BEFORE-target branch
    /// (the opposite orientation from Cas9), so the boundary work does not break the
    /// alternate code path. Cas12a PAM is TTTV (V = A/C/G), guide length 23, PAM 5'
    /// of the target. "TTTA" + 23 valid bases places a TTTV at index 0 with a 23-nt
    /// target immediately downstream. We pin that at least one Cas12a site is found
    /// with the PAM upstream of (before) its target start.
    /// </summary>
    [Test]
    public void FindPamSites_Cas12aPamBeforeTarget_DetectedWithCorrectOrientation()
    {
        string seq = "TTTA" + new string('C', 23); // TTTV at 0, 23-nt target downstream

        var forward = CrisprDesigner.FindPamSites(seq, CrisprSystemType.Cas12a)
            .Where(s => s.IsForwardStrand)
            .ToList();

        var site = forward.Should().ContainSingle(s => s.Position == 0).Subject;
        site.PamSequence.Should().Be("TTTA", "TTTV matches 'TTTA' (V = A) at the 5' end");
        site.System.PamAfterTarget.Should().BeFalse("Cas12a is a PAM-before-target system");
        site.TargetStart.Should().Be(4, "the 23-nt target lies immediately DOWNSTREAM of the 4-nt PAM");
        site.TargetSequence.Should().Be(new string('C', 23), "the protospacer is the 23 C's following the PAM");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// produce only well-formed sites — no out-of-bounds targets, no malformed PAM
    /// strings, no hang — so the degenerate-boundary guards do not corrupt the scan
    /// on ordinary input. Every site must satisfy INV-01 (PAM matches NGG under
    /// IUPAC) and INV-02 (target within bounds) regardless of random content.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void FindPamSites_RandomSequence_ProducesOnlyWellFormedSites()
    {
        string seq = RandomDna(2000, seed: 18_001);

        var sites = CrisprDesigner.FindPamSites(seq, CrisprSystemType.SpCas9).ToList();

        sites.Should().NotBeEmpty("a 2000-nt random sequence is overwhelmingly likely to contain NGG sites with room for a guide");
        // INV-01 is strand-aware: a forward site's PamSequence matches NGG (pos1=G, pos2=G);
        // a reverse site's PamSequence is reverse-complemented back to forward orientation, so it
        // reads as CCN (pos0=C, pos1=C) — the reverse complement of an NGG read on the other strand.
        sites.Should().OnlyContain(s =>
                s.PamSequence.Length == 3 &&
                (s.IsForwardStrand
                    ? (s.PamSequence[1] == 'G' && s.PamSequence[2] == 'G')
                    : (s.PamSequence[0] == 'C' && s.PamSequence[1] == 'C')) &&  // INV-01
                s.TargetSequence.Length == 20 &&                                // guide length
                s.TargetStart >= 0,                                             // INV-02
            "every site on random input is well-formed: an NGG PAM (forward) or its CCN reverse-complement (reverse) with a 20-nt in-bounds target");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CRISPR-GUIDE-001 — guide RNA design : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CRISPR-GUIDE-001 — guide RNA design

    #region BE — No PAM sites in the sequence

    /// <summary>
    /// BE (no PAM sites in seq): an all-A region contains no NGG, so the PAM scan
    /// underlying DesignGuideRnas finds nothing and the iterator yields the empty
    /// result — never a crash and never a fabricated guide. The region is valid
    /// (in-bounds), so the no-PAM behaviour is isolated from the range guards.
    /// DesignGuideRnas is a yield iterator, so we force enumeration to run the scan.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignGuideRnas_NoPamSites_IsEmptyAndDoesNotThrow()
    {
        var seq = new DnaSequence(new string('A', 60)); // no NGG anywhere

        var act = () => CrisprDesigner.DesignGuideRnas(seq, regionStart: 10, regionEnd: 50,
            CrisprSystemType.SpCas9).ToList();

        act.Should().NotThrow(
            "a sequence with no NGG PAM has no candidate guides; the scan yields nothing instead of crashing");

        CrisprDesigner.DesignGuideRnas(seq, 10, 50, CrisprSystemType.SpCas9)
            .Should().BeEmpty("no PAM site means no PAM-adjacent guide can be extracted");
    }

    #endregion

    #region BE — Guide span longer than the sequence

    /// <summary>
    /// BE (guide length > seq): the system's protospacer (20 nt for SpCas9) cannot
    /// fit even though an NGG may match — the sequence is shorter than PAM+guide. The
    /// PAM scanner's in-loop `targetStart >= 0 && targetEnd &lt; seq.Length` guard
    /// suppresses every site (no negative-start Substring), so DesignGuideRnas yields
    /// the empty result rather than an IndexOutOfRangeException. "AAGG" holds an NGG at
    /// index 1 but the 20-nt target would start at 1 − 20 = −19, so nothing is yielded.
    /// The region spans the whole 4-nt sequence (in-bounds), isolating the guide-fit
    /// boundary from the range guards.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignGuideRnas_GuideSpanLongerThanSequence_IsEmptyAndDoesNotThrow()
    {
        var seq = new DnaSequence("AAGG"); // NGG fits, 20-nt protospacer cannot

        var act = () => CrisprDesigner.DesignGuideRnas(seq, regionStart: 0, regionEnd: 3,
            CrisprSystemType.SpCas9).ToList();

        act.Should().NotThrow(
            "the 20-nt protospacer cannot fit upstream of the PAM in a 4-nt sequence; the bounds guard suppresses the site instead of slicing past the start");

        CrisprDesigner.DesignGuideRnas(seq, 0, 3, CrisprSystemType.SpCas9)
            .Should().BeEmpty("no protospacer fits, so no guide candidate is produced");
    }

    /// <summary>
    /// BE: the region-bounds guards are the upper boundary of "region outside the
    /// sequence". DesignGuideRnas range-checks regionStart and regionEnd against the
    /// sequence length and throws ArgumentOutOfRangeException — never an internal
    /// IndexOutOfRange while scanning (Guide_RNA_Design.md §3.3, §6.1). Pinned for a
    /// region start past the end and a region end past the end.
    /// </summary>
    [Test]
    public void DesignGuideRnas_RegionOutsideSequence_ThrowsArgumentOutOfRange()
    {
        var seq = new DnaSequence(new string('A', 20) + "AGG"); // length 23

        var startPastEnd = () => CrisprDesigner.DesignGuideRnas(seq, regionStart: 100, regionEnd: 100).ToList();
        var endPastEnd = () => CrisprDesigner.DesignGuideRnas(seq, regionStart: 0, regionEnd: 100).ToList();

        startPastEnd.Should().Throw<ArgumentOutOfRangeException>(
            "a region start beyond the sequence is rejected eagerly, not indexed into");
        endPastEnd.Should().Throw<ArgumentOutOfRangeException>(
            "a region end beyond the sequence is rejected eagerly, not indexed into");
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no typed input". DesignGuideRnas
    /// guards it with ArgumentNullException (ThrowIfNull) — never a NullReferenceException
    /// dereferencing the sequence (Guide_RNA_Design.md §3.3, §6.1).
    /// </summary>
    [Test]
    public void DesignGuideRnas_NullSequence_ThrowsArgumentNullException()
    {
        var act = () => CrisprDesigner.DesignGuideRnas((DnaSequence)null!, 0, 0).ToList();

        act.Should().Throw<ArgumentNullException>(
            "the design surface null-guards its sequence; null is rejected, never dereferenced");
    }

    #endregion

    #region BE — Guide length 0 (degenerate empty guide) — KEY div-by-zero hazard

    /// <summary>
    /// BE (guide length 0 — KEY): a zero-length (empty) guide is meaningless and is the
    /// division-by-zero hazard — CalculateSelfComplementarity divides by
    /// <c>length * length</c> and the GC computation would face an empty span. The
    /// source short-circuits this BEFORE any scoring: EvaluateGuideRna throws
    /// ArgumentNullException on a null/empty guide (CrisprDesigner.cs line 205), so the
    /// degenerate guide is a DOCUMENTED validation throw — NOT a DivideByZeroException,
    /// IndexOutOfRangeException, or a NaN score (Guide_RNA_Design.md §6.1). Pinned for
    /// both the empty string and null.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_ZeroLengthGuide_ThrowsArgumentNullException()
    {
        var empty = () => CrisprDesigner.EvaluateGuideRna(string.Empty);
        var nullGuide = () => CrisprDesigner.EvaluateGuideRna((string)null!);

        empty.Should().Throw<ArgumentNullException>(
            "a zero-length guide is degenerate; the source rejects it up front rather than dividing by length² or scoring an empty span");
        nullGuide.Should().Throw<ArgumentNullException>(
            "a null guide is rejected by the same guard, never dereferenced");
    }

    #endregion

    #region MC — Non-DNA characters in the guide / sequence

    /// <summary>
    /// MC (non-DNA, typed design surface): non-DNA junk fed to DesignGuideRnas is
    /// rejected at DnaSequence construction — the validated type throws ArgumentException
    /// on the first non-A/C/G/T character, so the guide-design scanner only ever sees
    /// clean DNA (Guide_RNA_Design.md §3.1; DnaSequence ValidateSequence).
    /// </summary>
    [Test]
    public void DesignGuideRnas_NonDnaSequence_RejectedAtConstruction()
    {
        var act = () => new DnaSequence("ACGT$#@!XYZ123");

        act.Should().Throw<ArgumentException>(
            "the validated DnaSequence type rejects non-DNA, so guide design never scans junk");
    }

    /// <summary>
    /// MC (non-DNA, raw evaluator): the standalone EvaluateGuideRna scores a raw string
    /// without A/C/G/T validation, so it must TOLERATE junk without crashing. Non-DNA
    /// characters are excluded from GC (CalculateGcFraction counts only A/C/G/T/U and
    /// returns 0 when none are valid) and the reverse complement used by
    /// self-complementarity passes non-IUPAC chars through unchanged (GetComplementBase
    /// fall-through) — no exception, no out-of-range indexing, and no division by zero
    /// because the guide is non-empty (length > 0). We pin: a pure-junk guide scores
    /// without throwing, with a finite, clamped (>= 0) score (INV-01) and GC of 0.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void EvaluateGuideRna_NonDnaGuide_ScoresWithoutCrash()
    {
        var act = () => CrisprDesigner.EvaluateGuideRna("$#@!XYZ123 ￿qwerty");

        act.Should().NotThrow(
            "non-DNA characters are excluded from GC and pass through the complement unchanged; neither path indexes out of range or divides by zero");

        var candidate = CrisprDesigner.EvaluateGuideRna("$#@!XYZ123 ￿qwerty");
        candidate.Score.Should().BeGreaterThanOrEqualTo(0,
            "INV-01: the score is clamped to >= 0 even on garbage input");
        candidate.Score.Should().NotBe(double.NaN, "a junk guide must not produce a NaN score");
        candidate.GcContent.Should().Be(0, "no character is A/C/G/T, so GC content is 0, not a divide-by-zero");
    }

    #endregion

    #region Positive sanity — a clear NGG yields a correct-length guide adjacent to the PAM

    /// <summary>
    /// Positive sanity: alongside the degenerate probes, a textbook SpCas9 construct must
    /// yield a real guide of the CORRECT length (20 nt) immediately adjacent to (upstream
    /// of) a clear NGG PAM, so the boundary hardening never silently breaks the core
    /// function. The sequence is 20 distinct-ish bases (the protospacer) + "AGG": the NGG
    /// matches at PAM index 20, and the 20-nt forward protospacer starts at index 0. We
    /// target a region covering the Cas9 cut site (3 bp upstream of the PAM) and pin that a
    /// candidate exists with a 20-nt forward guide equal to that protospacer, INV-01
    /// (Score >= 0) and INV-02 (Score >= MinScore) holding. MinScore is lowered so the
    /// candidate is not filtered out purely on heuristic composition.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignGuideRnas_ClearNggPam_YieldsCorrectLengthGuideAdjacentToPam()
    {
        // 20-nt protospacer with mixed bases (avoids extreme GC/poly-T heuristics) + NGG PAM.
        const string protospacer = "ACGTACGTACGTACGTACGT"; // 20 nt, GC = 50%
        var seq = new DnaSequence(protospacer + "AGG"); // PAM "AGG" (NGG) at index 20
        var loose = GuideRnaParameters.Default with { MinScore = 0 };

        // Cas9 cuts 3 bp upstream of the PAM (index 17); region must cover that cut site.
        var candidates = CrisprDesigner.DesignGuideRnas(seq, regionStart: 0, regionEnd: 22,
            CrisprSystemType.SpCas9, loose).ToList();

        var guide = candidates.Should().ContainSingle(c => c.IsForwardStrand && c.Sequence == protospacer).Subject;
        guide.Sequence.Length.Should().Be(20, "the SpCas9 protospacer is 20 nt");
        guide.Position.Should().Be(0, "the guide is extracted immediately upstream of the PAM, starting at index 0");
        guide.Score.Should().BeGreaterThanOrEqualTo(0, "INV-01: the score is clamped to >= 0");
        candidates.Should().OnlyContain(c => c.Score >= loose.MinScore,
            "INV-02: DesignGuideRnas yields only candidates meeting MinScore");
    }

    /// <summary>
    /// Positive sanity (standalone evaluator): a clean, well-formed 20-nt guide must be
    /// accepted and scored with finite, in-range metrics — GC in [0, 100], a clamped
    /// score (INV-01), and the system's 20-nt length echoed back — so the degenerate-guide
    /// guard does not break ordinary evaluation. A fixed-seed random guide keeps this
    /// deterministic.
    /// </summary>
    [Test]
    public void EvaluateGuideRna_CleanGuide_ProducesWellFormedScore()
    {
        string guide = RandomDna(20, seed: 19_001);

        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        candidate.Sequence.Should().Be(guide, "the evaluator scores the guide as given (upper-cased)");
        candidate.GcContent.Should().BeInRange(0, 100, "GC content is a percentage");
        candidate.Score.Should().BeInRange(0, 100, "INV-01: the heuristic score is clamped to [0, 100]");
        candidate.System.GuideLength.Should().Be(20, "SpCas9's documented guide length is 20 nt");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CRISPR-OFF-001 — off-target analysis : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CRISPR-OFF-001 — off-target analysis

    #region BE — Zero mismatch tolerance (maxMismatches = 0)

    /// <summary>
    /// BE (zero mismatch tolerance — KEY): <c>maxMismatches = 0</c> is the lower end of
    /// the documented 0..5 range, so it is IN-range and must NOT throw (the guard rejects
    /// only <c>&lt; 0</c>, CrisprDesigner.cs line 340). But the yield condition
    /// <c>mismatches &gt; 0 &amp;&amp; mismatches &lt;= 0</c> is unsatisfiable, so NO site
    /// is ever emitted — not even an EXACT on-target, which is excluded by the
    /// <c>mismatches &gt; 0</c> filter anyway (INV-01). Here the genome contains a PERFECT
    /// copy of the guide (the would-be on-target) AND a 1-mismatch variant; at tolerance 0
    /// the result must be EMPTY (exact excluded, 1-mismatch excluded), proving zero
    /// tolerance yields no off-targets rather than crashing or dividing by zero.
    /// FindOffTargets is a yield iterator, so we force enumeration.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOffTargets_ZeroMismatchTolerance_YieldsNothing()
    {
        const string guide = "ACGTACGTACGTACGTACGT"; // 20 nt
        // Two SpCas9 forward sites: a PERFECT copy of the guide (would-be on-target) and a
        // 1-mismatch variant (first base T instead of A), each followed by its own NGG PAM.
        const string oneMismatch = "TCGTACGTACGTACGTACGT"; // differs from guide at index 0
        var genome = new DnaSequence(guide + "AGG" + oneMismatch + "AGG");

        var act = () => CrisprDesigner.FindOffTargets(guide, genome, maxMismatches: 0,
            CrisprSystemType.SpCas9).ToList();

        act.Should().NotThrow(
            "maxMismatches = 0 is in the documented 0..5 range, so it is accepted, not rejected");

        CrisprDesigner.FindOffTargets(guide, genome, maxMismatches: 0, CrisprSystemType.SpCas9)
            .Should().BeEmpty(
                "at zero tolerance the yield condition 0 < mm <= 0 is unsatisfiable: the exact on-target is excluded (INV-01) and the 1-mismatch site exceeds the cap, so nothing is returned");
    }

    /// <summary>
    /// BE: the negative tolerance is the boundary just BELOW zero — distinct from the
    /// in-range zero. <c>maxMismatches = -1</c> is rejected eagerly with
    /// ArgumentOutOfRangeException (CrisprDesigner.cs line 340–341), as is the upper
    /// out-of-range value <c>6</c> (the guard rejects <c>&gt; 5</c>). Pinned on
    /// enumeration so the in-iterator guard actually fires.
    /// </summary>
    [Test]
    public void FindOffTargets_OutOfRangeMismatchTolerance_ThrowsArgumentOutOfRange()
    {
        var genome = new DnaSequence(new string('A', 20) + "AGG");

        var negative = () => CrisprDesigner.FindOffTargets("ACGTACGTACGTACGTACGT", genome, -1).ToList();
        var tooLarge = () => CrisprDesigner.FindOffTargets("ACGTACGTACGTACGTACGT", genome, 6).ToList();

        negative.Should().Throw<ArgumentOutOfRangeException>(
            "a negative mismatch tolerance is below the documented 0..5 range and is rejected, not silently clamped");
        tooLarge.Should().Throw<ArgumentOutOfRangeException>(
            "a tolerance above 5 is above the documented range and is rejected");
    }

    #endregion

    #region BE — Empty guide (degenerate zero-length guide)

    /// <summary>
    /// BE (empty guide): a zero-length (or null) guide is degenerate and is rejected up
    /// front with ArgumentNullException (CrisprDesigner.cs line 337–338) — BEFORE the
    /// length check, the PAM scan, or any CountMismatches loop, so there is no
    /// IndexOutOfRange and no division by a zero guide/seed length
    /// (Off_Target_Analysis.md §3.3, §6.1). The throw lives inside the yield iterator,
    /// so enumeration is forced. Pinned for both the empty string and null.
    /// </summary>
    [Test]
    public void FindOffTargets_EmptyGuide_ThrowsArgumentNullException()
    {
        var genome = new DnaSequence(new string('A', 20) + "AGG");

        var empty = () => CrisprDesigner.FindOffTargets(string.Empty, genome, 3, CrisprSystemType.SpCas9).ToList();
        var nullGuide = () => CrisprDesigner.FindOffTargets((string)null!, genome, 3, CrisprSystemType.SpCas9).ToList();

        empty.Should().Throw<ArgumentNullException>(
            "an empty guide is degenerate; it is rejected before any scoring or scan, never dividing by a zero guide length");
        nullGuide.Should().Throw<ArgumentNullException>(
            "a null guide is rejected by the same guard, never dereferenced");
    }

    #endregion

    #region MC — Guide of all N's (wildcard vs literal)

    /// <summary>
    /// MC (all-N guide — KEY): 'N' is the IUPAC "any base" code, so the central question
    /// is whether off-target matching treats it as a WILDCARD. This API does NOT:
    /// CountMismatches is a literal position-by-position <c>!=</c> compare
    /// (CrisprDesigner.cs lines 393–405) with no IUPAC semantics. A 20-nt all-N guide
    /// PASSES the SpCas9 length check (length 20 == GuideLength), then MISMATCHES every
    /// A/C/G/T base of any PAM-adjacent target, so each site's mismatch count is the full
    /// guide length (20) — far above the 0..5 cap — and NO site is ever yielded. So all-N
    /// is the OPPOSITE of a wildcard blow-up: instead of matching every position (a hang /
    /// off-target-everywhere hazard) it matches NOTHING, because every position is a
    /// guaranteed mismatch. Theory-correct result on a real, PAM-rich genome: the EMPTY
    /// enumerable, promptly (CancelAfter guards against any hang), with no exception.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOffTargets_AllNGuide_MatchesNothingAndDoesNotHang()
    {
        const string allN = "NNNNNNNNNNNNNNNNNNNN"; // 20 N's == SpCas9 guide length
        // A PAM-rich genome so the scan actually has many targets to compare against.
        var genome = new DnaSequence(RandomDna(500, seed: 20_001));

        var act = () => CrisprDesigner.FindOffTargets(allN, genome, maxMismatches: 5,
            CrisprSystemType.SpCas9).ToList();

        act.Should().NotThrow(
            "'N' is compared literally (not as a wildcard); it differs from every A/C/G/T base but never indexes out of range or hangs");

        CrisprDesigner.FindOffTargets(allN, genome, maxMismatches: 5, CrisprSystemType.SpCas9)
            .Should().BeEmpty(
                "an all-N guide mismatches every target base, so each candidate has 20 mismatches — above the cap of 5 — and no off-target is ever yielded; 'N' is not a wildcard");
    }

    /// <summary>
    /// MC/BE: a WRONG-length all-N guide is the dual boundary — it is rejected by the
    /// guide-length check with ArgumentException BEFORE any matching
    /// (Off_Target_Analysis.md §3.3, §6.1), pinning that the length contract is enforced
    /// regardless of the (junk) content of the guide.
    /// </summary>
    [Test]
    public void FindOffTargets_WrongLengthAllNGuide_ThrowsArgumentException()
    {
        var genome = new DnaSequence(new string('A', 20) + "AGG");

        // 10 N's: not equal to SpCas9's 20-nt guide length.
        var act = () => CrisprDesigner.FindOffTargets("NNNNNNNNNN", genome, 3, CrisprSystemType.SpCas9).ToList();

        act.Should().Throw<ArgumentException>(
            "a guide whose length does not match the system's guide length is rejected, even when its content is all-N");
    }

    #endregion

    #region MC — Null genome

    /// <summary>
    /// MC/INJ: a null genome is the boundary of "no search space". FindOffTargets
    /// null-guards the genome with ArgumentNullException (ThrowIfNull, line 339) — never a
    /// NullReferenceException dereferencing <c>genome.Sequence</c>
    /// (Off_Target_Analysis.md §3.3). The guard is inside the yield iterator, so
    /// enumeration is forced.
    /// </summary>
    [Test]
    public void FindOffTargets_NullGenome_ThrowsArgumentNullException()
    {
        var act = () => CrisprDesigner.FindOffTargets("ACGTACGTACGTACGTACGT", (DnaSequence)null!, 3,
            CrisprSystemType.SpCas9).ToList();

        act.Should().Throw<ArgumentNullException>(
            "the genome is null-guarded; null is rejected, never dereferenced into a NullReferenceException");
    }

    #endregion

    #region Positive sanity — a known 1-mismatch off-target is found at tolerance >= 1 but NOT at 0

    /// <summary>
    /// Positive sanity (the contract pivot): a genome carrying a guide-length target that
    /// differs from the guide by EXACTLY ONE base, immediately upstream of a real NGG PAM,
    /// must be reported as an off-target with Mismatches == 1 at tolerance >= 1, but must
    /// NOT be reported at tolerance 0 — so the boundary hardening never silently disables
    /// the core detection. The genome is the 1-mismatch protospacer followed by "AGG"
    /// (NGG at index 20). At tolerance 1: one forward site with Mismatches == 1 and the
    /// mismatch recorded at index 0 (INV-02, bounded by the cap). At tolerance 0: empty
    /// (INV-01 / unsatisfiable yield condition). This is the exact contrast the row probes.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOffTargets_KnownOneMismatchSite_FoundAtToleranceOneNotZero()
    {
        const string guide = "ACGTACGTACGTACGTACGT";       // 20 nt
        const string oneMismatch = "TCGTACGTACGTACGTACGT";  // differs from guide ONLY at index 0
        var genome = new DnaSequence(oneMismatch + "AGG");   // 1-mismatch target + NGG PAM at index 20

        var atOne = CrisprDesigner.FindOffTargets(guide, genome, maxMismatches: 1, CrisprSystemType.SpCas9)
            .Where(o => o.IsForwardStrand)
            .ToList();

        var site = atOne.Should().ContainSingle(o => o.Position == 20).Subject;
        site.Mismatches.Should().Be(1, "the off-target differs from the guide by exactly one base");
        site.Sequence.Should().Be(oneMismatch, "the reported off-target target is the PAM-adjacent protospacer");
        site.MismatchPositions.Should().ContainSingle().Which.Should().Be(0,
            "the single mismatch is at index 0, where the guide's 'A' meets the target's 'T'");

        CrisprDesigner.FindOffTargets(guide, genome, maxMismatches: 0, CrisprSystemType.SpCas9)
            .Should().BeEmpty(
                "at zero tolerance the same 1-mismatch site is excluded: 0 < 1 <= 0 is false");
    }

    /// <summary>
    /// Positive sanity (specificity score): CalculateSpecificityScore returns 100 when no
    /// off-targets exist (Off_Target_Analysis.md §6.1) and a value strictly in
    /// <c>[0, 100]</c> otherwise (INV-03), with the SCORE STRICTLY LOWER once a real
    /// off-target is present — so the boundary work does not corrupt the score path. A
    /// guide with NO near-match in an all-A genome (no NGG even) scores 100; the same
    /// guide against a genome carrying a 1-mismatch site scores below 100.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void CalculateSpecificityScore_ReflectsPresenceOfOffTargets()
    {
        const string guide = "ACGTACGTACGTACGTACGT";
        const string oneMismatch = "TCGTACGTACGTACGTACGT";

        var noSites = new DnaSequence(new string('A', 60)); // no NGG -> no off-targets
        var withSite = new DnaSequence(oneMismatch + "AGG"); // one 1-mismatch off-target

        double clean = CrisprDesigner.CalculateSpecificityScore(guide, noSites, CrisprSystemType.SpCas9);
        double withOff = CrisprDesigner.CalculateSpecificityScore(guide, withSite, CrisprSystemType.SpCas9);

        clean.Should().Be(100.0, "no off-targets means maximum specificity per the documented edge case");
        withOff.Should().BeInRange(0, 100, "INV-03: the specificity score is clamped to [0, 100]");
        withOff.Should().BeLessThan(100.0, "a real off-target reduces specificity below the no-off-target maximum");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PRIMER-TM-001 — primer melting temperature : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PRIMER-TM-001 — primer melting temperature

    #region BE — Empty sequence

    /// <summary>
    /// BE (empty seq): the empty string is the lower size boundary. The source
    /// short-circuits `string.IsNullOrEmpty` to Tm = 0 BEFORE any base count or
    /// Marmur-Doty division (PrimerDesigner.cs 199–200), so the degenerate zero-length
    /// primer is a DEFINED boundary value, never an IndexOutOfRange, DivideByZero, or
    /// NaN (Melting_Temperature.md §6.1 "Empty / length-1 → 0"). The salt-corrected
    /// sibling shares the same guard.
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_EmptySequence_ReturnsZeroFinite()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperature(string.Empty);

        tm.Should().Be(0, "an empty primer is short-circuited to a defined Tm of 0 before any count or division");
        double.IsNaN(tm).Should().BeFalse("the empty-input guard returns 0, never a NaN");
        double.IsInfinity(tm).Should().BeFalse("the empty-input guard returns a finite 0");

        PrimerDesigner.CalculateMeltingTemperatureWithSalt(string.Empty)
            .Should().Be(0, "the salt-corrected sibling shares the same empty-input short-circuit");
    }

    #endregion

    #region BE — Single base (1-bp)

    /// <summary>
    /// BE (1-bp): a single base is the lower NON-empty size boundary. validLength is 1,
    /// which is below the 14-nt Wallace threshold, so the Wallace rule applies and yields
    /// a TINY, EXACT Tm — an A·T base contributes 2 °C and a G·C base contributes 4 °C
    /// (Melting_Temperature.md §2.2 INV-01: Tm = 2·(A+T) + 4·(G+C)). No division by
    /// length, so no NaN. We pin all four single bases against the documented formula.
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_SingleBase_ReturnsWallaceTm()
    {
        // Wallace: A/T contribute 2, G/C contribute 4 (ThermoConstants 2/4 contributions).
        PrimerDesigner.CalculateMeltingTemperature("A").Should().Be(2, "one A·T base: Wallace Tm = 2·1 + 4·0 = 2 °C");
        PrimerDesigner.CalculateMeltingTemperature("T").Should().Be(2, "one A·T base: Wallace Tm = 2 °C");
        PrimerDesigner.CalculateMeltingTemperature("G").Should().Be(4, "one G·C base: Wallace Tm = 2·0 + 4·1 = 4 °C");
        PrimerDesigner.CalculateMeltingTemperature("C").Should().Be(4, "one G·C base: Wallace Tm = 4 °C");

        double tm = PrimerDesigner.CalculateMeltingTemperature("G");
        double.IsNaN(tm).Should().BeFalse("a 1-bp Wallace Tm is finite, with no length division");
        double.IsInfinity(tm).Should().BeFalse();
    }

    #endregion

    #region BE — Long primer (100+ bp) — overflow / non-finite hazard

    /// <summary>
    /// BE (100+ bp): a long primer is the UPPER size boundary and the overflow / non-finite
    /// leakage hazard. validLength >= 14 selects Marmur-Doty
    /// (Tm = 64.9 + 41·(GC − 16.4)/N, Melting_Temperature.md §2.2 INV-02). A 100-nt all-GC
    /// primer pins the exact value 64.9 + 41·(100 − 16.4)/100 = 99.176 °C — a LARGE but
    /// FINITE Tm, proving the long-primer path neither overflows nor leaks Inf/NaN. A
    /// 100-nt all-A·T primer (GC = 0) yields 64.9 + 41·(0 − 16.4)/100 = 58.176 °C, also
    /// finite and positive — the documented `max(0, …)` clamp is never even needed in the
    /// valid >= 14 regime, confirming the formula stays well-behaved at the upper boundary.
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_LongPrimer_IsFiniteMarmurDotyTm()
    {
        string allGc = new string('G', 100);
        string allAt = new string('A', 100);

        double gcTm = PrimerDesigner.CalculateMeltingTemperature(allGc);
        double atTm = PrimerDesigner.CalculateMeltingTemperature(allAt);

        // Marmur-Doty: 64.9 + 41*(GC - 16.4)/N.
        gcTm.Should().BeApproximately(64.9 + 41.0 * (100 - 16.4) / 100, 1e-9,
            "a 100-nt all-GC primer follows Marmur-Doty: 64.9 + 41·(100 − 16.4)/100 = 99.176 °C");
        atTm.Should().BeApproximately(64.9 + 41.0 * (0 - 16.4) / 100, 1e-9,
            "a 100-nt all-A·T primer follows Marmur-Doty: 64.9 + 41·(0 − 16.4)/100 = 58.176 °C");

        double.IsNaN(gcTm).Should().BeFalse("a long-primer Tm must be finite, not NaN");
        double.IsInfinity(gcTm).Should().BeFalse("a long-primer Tm must not overflow to ±Infinity");
        double.IsNaN(atTm).Should().BeFalse();
        double.IsInfinity(atTm).Should().BeFalse();
        gcTm.Should().BeGreaterThan(atTm, "G·C pairs are more stable, so the all-GC primer melts higher");
    }

    /// <summary>
    /// BE: the Wallace↔Marmur-Doty switchover at validLength 14 is the model-selection
    /// boundary. A 13-nt all-A primer (validLength 13 &lt; 14) takes the Wallace path
    /// (Tm = 2·13 = 26 °C); a 14-nt all-A primer (validLength 14) crosses into
    /// Marmur-Doty (Tm = 64.9 + 41·(0 − 16.4)/14 ≈ 16.871 °C). Pinning both sides proves
    /// the threshold (ThermoConstants.WallaceMaxLength = 14) is applied exactly and that
    /// neither side produces a non-finite value.
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_WallaceMarmurDotyThreshold_SwitchesAt14()
    {
        double wallace13 = PrimerDesigner.CalculateMeltingTemperature(new string('A', 13));
        double marmur14 = PrimerDesigner.CalculateMeltingTemperature(new string('A', 14));

        wallace13.Should().Be(2 * 13, "13 valid bases is below the 14-nt threshold: Wallace Tm = 2·13 = 26 °C");
        marmur14.Should().BeApproximately(64.9 + 41.0 * (0 - 16.4) / 14, 1e-9,
            "14 valid bases crosses into Marmur-Doty: 64.9 + 41·(0 − 16.4)/14 ≈ 16.871 °C");

        double.IsNaN(marmur14).Should().BeFalse("the Marmur-Doty side of the threshold is finite");
        double.IsInfinity(marmur14).Should().BeFalse();
    }

    #endregion

    #region INJ — All-N primer (KEY div-by-zero hazard)

    /// <summary>
    /// INJ (all-N — KEY div-by-zero hazard): 'N' is neither A/T nor G/C, so for an all-N
    /// primer at + gc = validLength = 0. Marmur-Doty divides by validLength, so a
    /// zero-valid-base input is the direct division-by-zero / NaN hazard. The source
    /// GUARDS it with an explicit `validLength == 0 → return 0` (PrimerDesigner.cs
    /// 208–209) BEFORE the division, so the theory-correct result is a finite 0 °C, NOT a
    /// NaN or DivideByZeroException. This holds for any non-DNA fill, so we also pin an
    /// all-'X' primer and a long all-N primer (still 0, no overflow).
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_AllNPrimer_ReturnsZeroNoDivideByZero()
    {
        double tmN = PrimerDesigner.CalculateMeltingTemperature(new string('N', 20));
        double tmX = PrimerDesigner.CalculateMeltingTemperature(new string('X', 30));
        double tmLongN = PrimerDesigner.CalculateMeltingTemperature(new string('N', 150));

        tmN.Should().Be(0, "an all-N primer has zero valid bases; the validLength==0 guard returns 0 before the Marmur-Doty division");
        tmX.Should().Be(0, "any non-A/C/G/T fill has zero valid bases and is guarded to 0");
        tmLongN.Should().Be(0, "a long all-N primer is still all-invalid: 0, never an overflow");

        double.IsNaN(tmN).Should().BeFalse("the zero-valid-base guard prevents a divide-by-zero NaN");
        double.IsInfinity(tmN).Should().BeFalse();
    }

    #endregion

    #region INJ — Non-DNA / special / unicode characters

    /// <summary>
    /// INJ (non-DNA chars): special characters, digits, spaces, and unicode are NOT
    /// counted as A/T or G/C — the method's contract states "Only standard DNA bases
    /// (A, C, G, T) are recognized; all other characters are ignored" (PrimerDesigner.cs
    /// 194–195). So a primer of PURE junk has validLength 0 and is guarded to 0 (no crash,
    /// no NaN), and junk INTERLEAVED into a real primer contributes nothing — the Tm equals
    /// that of the same primer with the junk stripped out. We pin both: pure junk → 0, and
    /// "AC$#GT 1￿" (4 valid bases A,C,G,T) → the Wallace Tm of "ACGT" (2 A·T + 2 G·C = 12 °C).
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_NonDnaChars_IgnoredNotCounted()
    {
        double junkOnly = PrimerDesigner.CalculateMeltingTemperature("$#@!123 \t￿qwz");
        junkOnly.Should().Be(0, "no character is A/C/G/T, so validLength is 0 and the guard returns 0, never a crash or NaN");
        double.IsNaN(junkOnly).Should().BeFalse();
        double.IsInfinity(junkOnly).Should().BeFalse();

        // "AC$#GT 1￿" has exactly the valid bases A,C,G,T (2 A·T + 2 G·C); junk is ignored.
        double interleaved = PrimerDesigner.CalculateMeltingTemperature("AC$#GT 1￿");
        double cleaned = PrimerDesigner.CalculateMeltingTemperature("ACGT");
        interleaved.Should().Be(cleaned,
            "injected special/unicode characters are ignored, so the Tm equals that of the cleaned 'ACGT' primer");
        interleaved.Should().Be(2 * 2 + 4 * 2, "Wallace Tm of A,C,G,T = 2·(A+T) + 4·(G+C) = 2·2 + 4·2 = 12 °C");
    }

    /// <summary>
    /// INJ: lower-case input is upper-cased before counting (PrimerDesigner.cs line 202),
    /// so a lower-case primer must give exactly the same Tm as its upper-case form — case
    /// is not a hidden injection vector that drops bases into the "ignored" bucket.
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_LowerCaseInput_SameAsUpperCase()
    {
        double lower = PrimerDesigner.CalculateMeltingTemperature("acgtacgtacgtacgt");   // 16 nt
        double upper = PrimerDesigner.CalculateMeltingTemperature("ACGTACGTACGTACGT");

        lower.Should().Be(upper, "input is upper-cased before counting, so case does not change which bases are valid");
        double.IsNaN(lower).Should().BeFalse();
    }

    #endregion

    #region INJ — Null reference

    /// <summary>
    /// INJ (null): a null reference is caught by `string.IsNullOrEmpty(primer)` and
    /// returns 0 — there is NEVER a NullReferenceException (the guard runs before any
    /// dereference; PrimerDesigner.cs 199–200). This API documents 0, not an
    /// ArgumentNullException, for null/empty (Melting_Temperature.md §3.3). The
    /// salt-corrected sibling shares the same null guard (line 229–230).
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_NullPrimer_ReturnsZeroNoNullReference()
    {
        var act = () => PrimerDesigner.CalculateMeltingTemperature((string)null!);

        act.Should().NotThrow<NullReferenceException>(
            "the null/empty guard runs before any dereference, so null can never raise a NullReferenceException");

        PrimerDesigner.CalculateMeltingTemperature((string)null!)
            .Should().Be(0, "null is treated as empty and returns a defined Tm of 0");
        PrimerDesigner.CalculateMeltingTemperatureWithSalt((string)null!)
            .Should().Be(0, "the salt-corrected sibling shares the same null short-circuit");
    }

    #endregion

    #region Positive sanity — a known primer Tm matches the documented formula

    /// <summary>
    /// Positive sanity: alongside the degenerate probes, a textbook 20-nt primer must
    /// produce the EXACT documented Tm so the boundary hardening never silently breaks the
    /// core calculation. "ACGTACGTACGTACGTACGT" has 10 A·T and 10 G·C (validLength 20 >= 14),
    /// so Marmur-Doty applies: Tm = 64.9 + 41·(10 − 16.4)/20 = 51.78 °C
    /// (Melting_Temperature.md §2.2 INV-02). A short 10-nt mixed primer with 4 A·T and 6 G·C
    /// (validLength 10 &lt; 14) takes Wallace: Tm = 2·4 + 4·6 = 32 °C (INV-01). Both pin
    /// the model selection AND the exact value, with a finite result.
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_KnownPrimers_MatchDocumentedFormula()
    {
        // 20-nt, GC = 50% (10 G·C, 10 A·T) → Marmur-Doty.
        double tm20 = PrimerDesigner.CalculateMeltingTemperature("ACGTACGTACGTACGTACGT");
        tm20.Should().BeApproximately(64.9 + 41.0 * (10 - 16.4) / 20, 1e-9,
            "Marmur-Doty Tm of a 20-nt GC-50% primer = 64.9 + 41·(10 − 16.4)/20 = 51.78 °C");

        // 10-nt with 4 A·T + 6 G·C → Wallace (validLength 10 < 14).
        double tm10 = PrimerDesigner.CalculateMeltingTemperature("ATATGCGCGC"); // A,T,A,T,G,C,G,C,G,C → AT=4, GC=6
        tm10.Should().Be(2 * 4 + 4 * 6, "Wallace Tm of a 10-nt primer with 4 A·T + 6 G·C = 2·4 + 4·6 = 32 °C");

        double.IsNaN(tm20).Should().BeFalse("a known primer Tm must be finite");
        double.IsInfinity(tm20).Should().BeFalse();
    }

    #endregion

    #endregion
}
