using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests.Fuzzing;

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
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PRIMER-DESIGN-001 — primer (pair) design
/// Checklist: docs/checklists/03_FUZZING.md, row 22.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — a template too short to hold ANY MinLength primer
///          (no candidate window fits), and the degenerate GC range 0–0 (MinGcContent =
///          MaxGcContent = 0, only all-A·T primers can qualify).
///   • MC = Malformed Content — a CONTRADICTORY (inverted) Tm range (MinTm > MaxTm), an
///          unsatisfiable constraint whose feasible set is empty.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 22 targets:
///   "Seq shorter than min primer, GC range 0-0, Tm range inverted".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The primer-design contract under test (Primer_Design.md)
/// ───────────────────────────────────────────────────────────────────────────
/// PrimerDesigner.DesignPrimers(DnaSequence template, int targetStart, int targetEnd,
///   PrimerParameters? parameters) (PrimerDesigner.cs lines 40–115) enumerates candidate
/// primers in the flanking regions (up to 200 bp upstream of targetStart for the forward
/// primer, up to 200 bp downstream of targetEnd for the reverse), scores each, and selects
/// the highest-scoring VALID forward and reverse candidate independently, then checks pair
/// compatibility (Primer_Design.md §2.2, §4.1). A candidate is VALID only when it satisfies
/// EVERY screen — length in [MinLength, MaxLength], GC% in [MinGcContent, MaxGcContent], Tm in
/// [MinTm, MaxTm], homopolymer/dinucleotide limits, no hairpin, and (if enabled) 3' stability
/// (EvaluatePrimer, lines 138–169: a candidate is valid iff its issue list is empty).
///
/// The ONLY documented throw is the target-region guard: targetStart &lt; 0,
/// targetEnd &gt;= template.Length, or targetStart &gt;= targetEnd → ArgumentException
/// (lines 48–49; §3.3, §6.1 "Invalid target region"). The PARAMETER RANGES themselves are
/// NOT validated — a degenerate or contradictory range is not an error but simply shrinks
/// (or empties) the feasible candidate set. INV-01 (Primer_Design.md §2.4) governs the
/// degenerate outcome: when either side has NO valid candidate, DesignPrimers returns an
/// INVALID PrimerPairResult (Forward = Reverse = null, IsValid = false, ProductSize = 0;
/// lines 93–100) — never a crash, hang, div-by-zero, or invalid primer.
///
/// THE THREE ROW-22 FUZZ TARGETS, mapped to the theory-correct contract:
///   • Seq shorter than min primer (BE): the template is too short to hold even one
///     MinLength-bp primer in either flank. The forward loop bound `start + len &lt;= targetStart`
///     (len &gt;= MinLength) and the reverse loop bound `end - len &gt;= targetEnd` are never
///     satisfiable, so NO candidate is generated, both best candidates are null, and INV-01
///     returns the invalid PrimerPairResult — no Substring is ever taken past the end, no
///     IndexOutOfRangeException. (The target region itself is kept VALID so this isolates the
///     "no primer fits" boundary from the region guard.)
///   • GC range 0–0 (BE, KEY): MinGcContent = MaxGcContent = 0 means the GC screen
///     `gc &lt; 0 || gc &gt; 0` passes ONLY for a 0% GC (all-A·T) primer. There is NO division by
///     the range width (the check is a pair of comparisons, not a normalisation), and
///     CalculateGcContent null/empty-guards to 0, so there is no div-by-zero. An all-A·T flank
///     can therefore still yield an all-A·T forward primer (GC exactly 0, inside [0,0]); a
///     GC-bearing flank yields none. Either way: a well-formed result, never a crash.
///   • Tm range inverted (MC, KEY): MinTm &gt; MaxTm is a CONTRADICTION whose feasible set is
///     empty. The source does NOT validate the range, so it does NOT throw — instead the Tm
///     screen `tm &lt; MinTm || tm &gt; MaxTm` is satisfied by EVERY finite Tm (any value is either
///     below the high MinTm or above the low MaxTm), so EVERY candidate is rejected, both sides
///     are null, and INV-01 returns the invalid PrimerPairResult. A contradictory constraint
///     thus yields NO primers promptly — it never loops forever and never returns a primer that
///     violates the (impossible) range.
///
/// Documented invariant pinned (Primer_Design.md §2.4): INV-01 DesignPrimers returns
/// IsValid = false (with null candidates) when either side has no valid candidate; INV-03
/// ProductSize = reverse.Position + reverse.Sequence.Length − forward.Position on a valid pair.
/// DesignPrimers is NOT an iterator (it returns a materialised PrimerPairResult), so every
/// probe calls it directly; the positive-sanity test additionally pins that a real pair's
/// chosen primers fall INSIDE the requested GC and Tm ranges.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PRIMER-STRUCT-001 — primer secondary structure (hairpin / dimer detection)
/// Checklist: docs/checklists/03_FUZZING.md, row 23.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — an EXTREMELY SHORT primer (1–3 bp, below the
///          minimum stem+loop+stem span), an EXTREMELY LONG primer (the >100-bp
///          suffix-tree path; overflow / non-finite / hang hazard), a PALINDROMIC
///          (self-reverse-complementary) primer (the maximal self-structure positive
///          signal), and an ALL-G homopolymer (no Watson–Crick self-complement).
/// — docs/checklists/03_FUZZING.md §Description (BE); row 23 targets:
///   "Extremely short, extremely long, palindromic primer, all-G".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The secondary-structure contract under test (Primer_Design.md §1, §2.5)
/// ───────────────────────────────────────────────────────────────────────────
/// Primer secondary structure (PRIMER-STRUCT-001; Primer_Design.md §2.5 names it the
/// "Hairpin and dimer detection (used in evaluation)" prerequisite) is detected by two
/// pure BOOLEAN predicates in PrimerDesigner (PrimerDesigner.cs):
///   • HasHairpinPotential(string, minStemLength = 4, minLoopLength = 3) → bool
///     (lines 307–323) — INTRAMOLECULAR self-complementarity: a stem-loop. It scans for
///     a minStemLength-bp fragment whose reverse-complement-pairing partner occurs at
///     least minLoopLength bases downstream (so a stem can fold back over a loop). It
///     GUARDS the degenerate small case explicitly: `IsNullOrEmpty(seq) ||
///     seq.Length < minStemLength*2 + minLoopLength` (i.e. < 4·2+3 = 11) → false
///     (lines 309–310). For seq.Length < 100 it uses the simple O(n²) scanner; for
///     >= 100 it uses the suffix-tree O(n) path (lines 316–322) — BOTH paths are exercised
///     here (short palindrome vs long stem-loop) so neither indexes out of range or hangs.
///   • HasPrimerDimer(string p1, string p2, minComplementarity = 4) → bool
///     (lines 398–419) — INTERMOLECULAR 3'-end complementarity: it reverse-complements p2,
///     compares the (≤ 8-base) 3' ends position-by-position, and returns true when at
///     least minComplementarity positions pair. It null/empty-guards BOTH primers to
///     false (lines 400–401).
/// Both predicates return a plain bool — there is NO ΔG, NO division, NO NaN/Inf surface;
/// the only failure modes a boundary input could provoke are an out-of-range Substring /
/// AsSpan slice, a suffix-tree blow-up on a long primer, or a wrong boolean. Neither method
/// throws for ANY string input (the empty/short guards run before any slicing).
///
/// THE FOUR ROW-23 FUZZ TARGETS, mapped to the theory-correct contract:
///   • Extremely short (BE): a 1–3 bp primer is far below the 11-base stem+loop+stem
///     minimum, so HasHairpinPotential short-circuits to FALSE (the guard, line 309) — too
///     short to fold a stem, NEVER an IndexOutOfRange on the Substring(i, minStemLength)
///     slices. HasPrimerDimer on two 1–3 bp primers is likewise FALSE: the 3' overlap is
///     shorter than minComplementarity (4), so the count can never reach the threshold.
///     LOW/zero self-structure, no crash.
///   • Extremely long (BE): a long primer drives the >= 100-bp SUFFIX-TREE hairpin path; the
///     result is a FINITE boolean computed promptly (no overflow, no hang — pinned with
///     [CancelAfter]). A long random primer returns a well-defined true/false, and a long
///     explicit stem-loop (60×G · loop · 60×C, length 130) returns TRUE via the suffix-tree
///     path — proving the long-primer branch is correct, not merely non-crashing.
///   • Palindromic primer (BE, KEY positive signal): a self-reverse-complementary primer
///     has MAXIMAL hairpin potential — every base can fold back onto its partner. A clean
///     inverted-repeat stem-loop ("GGGGAAACCCC", an 11-base stem-loop) and a true palindrome
///     whose reverse complement equals itself ("GGAATTCCGGAATTCC") both yield
///     HasHairpinPotential == TRUE. This is the load-bearing positive assertion: the
///     detector MUST flag the self-structure, not silently miss it.
///   • All-G homopolymer (BE): G does NOT Watson–Crick pair with G, so a poly-G primer has
///     NO intramolecular self-complement and HasHairpinPotential == FALSE on BOTH the short
///     O(n²) path (30×G) and the long suffix-tree path (150×G) — LOW/zero self-structure,
///     never a crash. (Biochemically distinct, and pinned as a deliberate note: two SEPARATE
///     poly-G primers DO form an inter-molecular dimer, because reverse-complementing the
///     second poly-G yields poly-C which pairs the first — so HasPrimerDimer(G…, G…) is TRUE.
///     That is correct and is NOT "self"-structure; the all-G self-structure signal is the
///     hairpin, which is correctly FALSE.)
///
/// Documented contract pinned (Primer_Design.md §2.5, §3.3): the structure screens are pure
/// heuristic booleans with explicit small-input guards; the positive-sanity test pins that an
/// ordinary primer flows through HasHairpinPotential AND HasPrimerDimer without throwing and
/// that a designed inverted repeat is detected. Both methods are pure (no iterator), so every
/// probe calls them directly.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROBE-DESIGN-001 — hybridization probe design
/// Checklist: docs/checklists/03_FUZZING.md, row 24.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — a target SHORTER than the minimum probe length
///          (no candidate window fits), and the EMPTY target (the lower size
///          boundary). Both must yield no probes, never an IndexOutOfRangeException.
///   • MC = Malformed Content — a CONTRADICTORY (inverted) Tm range (MinTm > MaxTm),
///          an unsatisfiable constraint whose feasible Tm set is empty.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 24 targets:
///   "Seq shorter than min probe, Tm range inverted, empty seq".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The probe-design contract under test (Hybridization_Probe_Design.md)
/// ───────────────────────────────────────────────────────────────────────────
/// ProbeDesigner.DesignProbes(string targetSequence, ProbeParameters? parameters = null,
///   int maxProbes = 10) (ProbeDesigner.cs lines 127–146) is a YIELD iterator that scans the
/// target for every length-`[MinLength..MaxLength]` window, scores each candidate by a fixed
/// additive penalty model (GC, Tm, homopolymer, self-complementarity, secondary structure,
/// repeats, terminal G/C — EvaluateProbeWithGc, lines 253–327), keeps only RAW-SCORE-POSITIVE
/// candidates (`score <= 0 → null`, line 315; INV-01), ranks them by score, and yields the top
/// `maxProbes`. Each Probe carries its Tm, GcContent, Score, and a Warnings list.
///
/// CRITICAL DESIGN FACT the checklist row probes — UNLIKE PrimerDesigner.DesignPrimers,
/// ProbeDesigner does NOT hard-reject a candidate when a screen fails. A failed screen only
/// DEDUCTS a penalty and APPENDS a warning; the candidate survives as long as its residual
/// score stays positive. So a contradictory Tm range does NOT empty the result — it makes the
/// Tm screen fail for EVERY candidate (a fixed −0.3 penalty plus a "Tm … outside range"
/// warning), but probes can still be returned. This is the theory-correct contract for THIS
/// surface; the test asserts it precisely (no throw, no hang, and every returned probe carries
/// the Tm warning) rather than the primer-design "empty result" — a contradictory constraint
/// must neither crash nor loop forever, and the source must NOT silently drop the warning.
///
/// The ONLY documented short-circuit (Hybridization_Probe_Design.md §3.3, §6.1): the target is
/// null, empty, OR shorter than `param.MinLength` → `yield break` (lines 134–135). There is NO
/// throw for any string target and NO validation of the parameter ranges — a degenerate or
/// contradictory range simply shrinks (Tm-warns) or empties the candidate set.
///
/// THE THREE ROW-24 FUZZ TARGETS, mapped to the theory-correct contract:
///   • Seq shorter than min probe (BE): a target shorter than `MinLength` (Microarray default
///     MinLength = 50) trips the length short-circuit `targetSequence.Length < param.MinLength`
///     (line 134) → `yield break`. The candidate-window loop `start <= n − length` (line 222)
///     would in any case never run when `length > n` (the `length <= n` guard, line 220), so no
///     Substring is taken past the end. Empty result, never an IndexOutOfRangeException.
///   • Tm range inverted (MC, KEY): MinTm > MaxTm is a CONTRADICTION whose feasible Tm set is
///     empty. The source does NOT validate the range, so it does NOT throw; instead the Tm
///     screen `tm < MinTm || tm > MaxTm` is satisfied by EVERY finite Tm, so EVERY candidate
///     incurs the −0.3 Tm penalty and records the "Tm … outside range" warning. Candidates with
///     enough residual score still return — so the theory-correct result is NOT empty here; it
///     is "no throw, no hang, and every returned probe carries the Tm-out-of-range warning".
///   • Empty seq (BE): `string.IsNullOrEmpty(targetSequence)` short-circuits to `yield break`
///     (line 134) BEFORE any GC prefix sum or windowing — a DEFINED degenerate boundary, never
///     a crash. Pinned for both "" and null target.
///
/// Documented invariants pinned (Hybridization_Probe_Design.md §2.4): INV-01 only raw-score-
/// positive candidates are retained (Score > 0 on every returned probe); INV-02 GcContent is a
/// fraction of length (0 <= GcContent <= 1). DesignProbes is a yield iterator, so every probe
/// forces enumeration (`.ToList()`); the positive-sanity test additionally pins that a real
/// probe's GC and Tm fall INSIDE the requested ranges with no warnings.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROBE-VALID-001 — probe validation
/// Checklist: docs/checklists/03_FUZZING.md, row 25.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — a probe with degenerate IUPAC 'N' bases (the
///          ambiguity wildcard fed as a literal probe base), and a probe of pure
///          non-DNA junk (symbols / digits / unicode), each fed to the un-validated
///          raw-string validation surface.
///   • INJ = Injection — special characters, null bytes, and unicode interleaved
///          into a probe string, plus the `null` reference itself.
/// — docs/checklists/03_FUZZING.md §Description (MC; INJ = injection of special chars /
///   null bytes / unicode); row 25 targets:
///   "Non-DNA probe, extremely short, null, probe with N's".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The probe-validation contract under test (Probe_Validation.md)
/// ───────────────────────────────────────────────────────────────────────────
/// ProbeDesigner.ValidateProbe(string probeSequence, IEnumerable&lt;string&gt;
/// referenceSequences, int maxMismatches = 3, double selfComplementarityThreshold = 0.3)
/// (ProbeDesigner.cs lines 491–564) assesses a candidate hybridization probe against a
/// set of reference sequences and returns a ProbeValidation record (IsValid,
/// SpecificityScore, OffTargetHits, SelfComplementarity, HasSecondaryStructure, Issues).
/// It is NOT an iterator — it returns a materialised record, so every probe calls it
/// directly. It uppercases the probe (line 500), counts substitution-tolerant
/// fixed-length approximate hits across every reference (FindApproximateMatches,
/// lines 789–804), maps the hit count to a specificity in [0,1] (0 hits → 0.0, 1 → 1.0,
/// N → 1/N; INV-01), measures self-complementarity as a fraction of positions matching
/// the reverse complement (CalculateSelfComplementarity, lines 721–733; INV-02), and
/// flags hairpin potential (HasSecondaryStructurePotential, lines 735–761).
///
/// CRITICAL DESIGN FACT the checklist row probes — ValidateProbe does NOT validate the
/// probe ALPHABET. Unlike the typed DnaSequence surface, the raw probe string is scanned
/// AS-IS (only upper-cased): there is NO A/C/G/T screen and NO minimum-length screen. The
/// ONLY documented throws are the two null guards (ArgumentNullException on a null probe
/// OR a null reference collection, lines 497–498). Every other degenerate input is HANDLED,
/// not rejected: the empty probe is a structured invalid result ("Empty probe sequence"
/// issue, SpecificityScore 0.0, OffTargetHits 0; lines 504–513) rather than a throw.
///
/// THE FOUR ROW-25 FUZZ TARGETS, mapped to the theory-correct contract:
///   • Non-DNA probe (MC/INJ): junk (symbols, digits, unicode, null bytes) is upper-cased
///     and scanned literally. FindApproximateMatches compares characters with a plain `!=`,
///     so junk simply mismatches A/C/G/T references (no crash, no inflated hits). The
///     reverse complement used by self-complementarity passes any non-IUPAC char THROUGH
///     unchanged (GetComplementBase fall-through arm, SequenceExtensions.cs line 156), so
///     a junk base equals its own "complement" and is counted as self-complementary — a
///     defined number, never an exception. Because the probe length is &gt; 0,
///     CalculateSelfComplementarity never divides by zero → NO NaN. The result is a
///     well-formed ProbeValidation with SpecificityScore and SelfComplementarity both
///     finite and in [0,1], never an IndexOutOfRange/NullReference.
///   • Extremely short probe (BE/INJ): a 1–3 bp probe has NO minimum-length screen, so it
///     is scanned, NOT rejected. FindApproximateMatches' bound `i &lt;= text.Length −
///     pattern.Length` simply runs over each reference; a 1-nt probe matches many
///     positions (high OffTargetHits, low specificity) but never crashes. The hairpin
///     screen's loop bound `i &lt;= len − stemLen*2 − 3` is negative for a tiny probe, so it
///     short-circuits to false. Self-complementarity divides by the (non-zero) length →
///     finite, no NaN. So an extremely short probe yields a finite, in-range result — the
///     theory-correct "handled, scored, never crash" outcome (the only length-based reject
///     is the EMPTY probe, pinned separately).
///   • Null (INJ, KEY): a null probe is the documented ArgumentNullException
///     (ThrowIfNull, line 497) — raised eagerly, NEVER a NullReferenceException from
///     dereferencing `probeSequence.ToUpperInvariant()`. A null reference COLLECTION is
///     the sibling throw (line 498). Both pinned; a null reference STRING inside a non-null
///     collection is a separate hazard pinned as NotThrow-or-documented.
///   • Probe with N's (MC, KEY): 'N' is the IUPAC "any base" code. ValidateProbe does NOT
///     treat 'N' as a wildcard in matching — FindApproximateMatches is a literal `!=`, so
///     an 'N' in the probe mismatches an A/C/G/T reference base (counts toward mismatches),
///     and 'N' is NOT counted as G/C. In self-complementarity 'N' complements to 'N'
///     (GetComplementBase, SequenceExtensions.cs line 155), so an all-N probe is maximally
///     self-complementary (every position N==N) — a DEFINED 1.0, never a NaN. An all-N
///     probe is therefore handled and scored (no A/C/G/T, finite specificity and
///     self-complementarity), never crashing and never producing NaN in any numeric field.
///
/// Documented invariants pinned on every result (Probe_Validation.md §2.4): INV-01
/// `0.0 &lt;= SpecificityScore &lt;= 1.0`; INV-02 `0.0 &lt;= SelfComplementarity &lt;= 1.0`;
/// INV-03 `OffTargetHits &gt;= 0`; INV-04 `OffTargetHits == 1` ⇒ `SpecificityScore == 1.0`.
/// ValidateProbe is a pure function (no iterator), so every probe calls it directly; the
/// positive-sanity test pins that a clean, unique probe validates as IsValid == true with
/// no issues, SpecificityScore 1.0, and OffTargetHits 1 (INV-04).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: RESTR-FIND-001 — restriction site finding
/// Checklist: docs/checklists/03_FUZZING.md, row 26.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — an UNKNOWN enzyme name (a name absent from the
///          built-in catalog) and non-DNA junk characters in the sequence.
///   • BE = Boundary Exploitation — the EMPTY enzyme list (the multi-enzyme overload
///          called with zero names) and the EMPTY sequence (lower size boundary).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 26 targets:
///   "Empty enzyme list, unknown enzyme names, empty sequence, non-DNA".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The restriction-site-finding contract under test (Restriction_Site_Detection.md)
/// ───────────────────────────────────────────────────────────────────────────
/// Restriction site detection slides a recognition-sequence window across BOTH
/// strands of a DNA sequence, matches each window character-by-character under IUPAC
/// ambiguity rules, and yields a RestrictionSite per match with the enzyme, strand,
/// cut position and matched bases (Restriction_Site_Detection.md §2.2, §4.1). The
/// enzyme is NOT a free-form motif — it is resolved by NAME from a fixed built-in
/// catalog (RestrictionAnalyzer.cs §Built-in Enzyme Database), so "an invalid PAM"
/// here is an UNKNOWN enzyme NAME. EcoRI = GAATTC (CutPositionForward 1); the catalog
/// also carries IUPAC-degenerate motifs (HincII = GTYRAC, SfiI = GGCCNNNNNGGCC), so
/// the matcher's IUPAC path is load-bearing (§5.2, §6.1).
///
/// API surfaces under test (RestrictionAnalyzer.cs lines 121–230):
///   • FindSites(DnaSequence, string enzymeName): null sequence → ArgumentNullException
///     (ThrowIfNull, line 129); null/empty enzymeName → ArgumentNullException (line
///     130–131); UNKNOWN enzymeName → ArgumentException "Unknown enzyme" (line 133–134,
///     raised EAGERLY because this overload is a plain method that returns
///     FindSitesCore's enumerable). The DnaSequence itself is a validated type whose
///     constructor rejects any non-A/C/G/T base (ArgumentException), so non-DNA never
///     reaches the scanner on this surface.
///   • FindSites(string sequence, string enzymeName): null/empty sequence → empty
///     result via `yield break` (line 144–145, INV-04); UNKNOWN enzymeName →
///     ArgumentException, but because the body is a `yield` ITERATOR the throw fires
///     only on ENUMERATION (line 147–148), so every probe forces `.ToList()`. The raw
///     sequence is upper-cased and scanned RAW (no A/C/G/T validation) — this is the MC
///     surface for non-DNA: a junk sequence char is tested as the `nucleotide` argument
///     of MatchesIupac(seqChar, patChar); it satisfies no enzyme-pattern code, so it
///     never matches a site and never crashes (IupacHelper.cs §MatchesIupac; the `_ =>
///     throw` arm is reached only by an invalid PATTERN code, and every catalog motif is
///     IUPAC-valid, so it is never reached from the sequence side).
///   • FindSites(DnaSequence, params string[] enzymeNames): null sequence →
///     ArgumentNullException (line 209). An EMPTY enzymeNames array iterates zero times,
///     so the result is the EMPTY enumerable — NOT a throw. This is the "empty enzyme
///     list" BE target: no enzymes ⇒ no sites, never a crash. (A NULL enzymeNames array
///     is a separate hazard — the `foreach` over null would NullReference — and is NOT
///     the row target; the row probes the empty list, which is the well-defined no-op.)
///
/// THE FOUR ROW-26 FUZZ TARGETS, mapped to the theory-correct contract:
///   • Empty enzyme list (BE): the params overload with `new string[0]` (or no names)
///     iterates no enzymes and yields the EMPTY result — never a crash, never a
///     fabricated site. Pinned on enumeration.
///   • Unknown enzyme names (MC, KEY): a name absent from the catalog is REJECTED with
///     ArgumentException "Unknown enzyme" — a DELIBERATE, documented validation throw
///     (Restriction_Site_Detection.md §3.3, §6.1), NOT a raw KeyNotFoundException
///     leaking from the dictionary. GetEnzyme uses TryGetValue (line 76–77), so the
///     lookup miss is mapped to a typed ArgumentException by the `?? throw` (line 133,
///     147), never a KeyNotFound. Pinned on BOTH surfaces: the typed overload throws
///     eagerly; the raw-string iterator throws on enumeration.
///   • Empty sequence (BE): the raw-string overload short-circuits null/empty to the
///     empty result (yield break, INV-04); the typed overload over an empty DnaSequence
///     has a negative forward-scan bound `i <= 0 − patternLen`, so neither the forward
///     nor the reverse loop runs and no Substring is taken past the end. Empty result,
///     never an IndexOutOfRangeException (§6.1). Pinned for both surfaces plus raw null.
///   • Non-DNA (MC): junk fed to the TYPED surface is rejected at DnaSequence
///     construction (ArgumentException); junk fed to the RAW-string surface is TOLERATED
///     — a non-A/C/G/T char matches no enzyme pattern position (and the reverse-strand
///     pass complements via GetReverseComplementString whose fall-through passes
///     non-IUPAC chars through unchanged), so pure junk yields no sites and never
///     crashes, and junk interleaved into a real GAATTC window suppresses that site.
///
/// Documented invariants pinned on every produced site (Restriction_Site_Detection.md
/// §2.4): INV-01 `0 <= Position <= len − recognitionLength`; INV-02
/// `RecognizedSequence.Length == Enzyme.RecognitionSequence.Length`; INV-03 enzyme-name
/// lookup is case-insensitive; INV-04 empty raw-string input yields no sites. The
/// raw-string and params overloads are yield iterators, so every probe forces
/// enumeration (`.ToList()`); the positive-sanity test pins that EcoRI's GAATTC site is
/// found at the CORRECT forward position with the CORRECT cut position and matched bases.
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

    // ═══════════════════════════════════════════════════════════════════
    //  PRIMER-DESIGN-001 — primer design : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PRIMER-DESIGN-001 — primer design

    #region BE — Template shorter than the minimum primer

    /// <summary>
    /// BE (seq shorter than min primer): the template is far too short to hold even one
    /// MinLength-bp primer in either flank. With a valid target region on a 4-nt template,
    /// the forward candidate loop bound <c>start + len &lt;= targetStart</c> (len >= MinLength,
    /// default 18) and the reverse bound <c>end - len >= targetEnd</c> are never satisfiable,
    /// so NO candidate is generated — no Substring is ever sliced past the template end. Both
    /// best candidates are null, so DesignPrimers returns the INVALID PrimerPairResult (INV-01,
    /// PrimerDesigner.cs 93–100): Forward/Reverse null, IsValid false, ProductSize 0 — an empty
    /// feasible set, never an IndexOutOfRangeException (Primer_Design.md §6.1).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignPrimers_TemplateShorterThanMinPrimer_ReturnsInvalidNoCrash()
    {
        // 4-nt template; target region [1,2) is VALID (1<2, 2<4) so this isolates the
        // "no primer fits" boundary from the target-region guard.
        var template = new DnaSequence("ACGT");

        var act = () => PrimerDesigner.DesignPrimers(template, targetStart: 1, targetEnd: 2);

        act.Should().NotThrow(
            "no MinLength-bp primer can fit in a 4-nt template; the candidate loops never run instead of slicing past the end");

        var result = PrimerDesigner.DesignPrimers(template, 1, 2);
        result.IsValid.Should().BeFalse("INV-01: no valid candidate on either side yields an invalid pair");
        result.Forward.Should().BeNull("no forward primer fits, so the forward candidate is null");
        result.Reverse.Should().BeNull("no reverse primer fits, so the reverse candidate is null");
        result.ProductSize.Should().Be(0, "an empty pair reports zero product size");
    }

    /// <summary>
    /// BE: the target-region guard is the documented throw boundary. Unlike a degenerate
    /// PARAMETER range (which is silently empty), an INVALID target region is rejected eagerly
    /// with ArgumentException (PrimerDesigner.cs 48–49; Primer_Design.md §6.1). Pinned for the
    /// three guard conditions: targetStart &lt; 0, targetEnd &gt;= length, and targetStart &gt;= targetEnd.
    /// </summary>
    [Test]
    public void DesignPrimers_InvalidTargetRegion_ThrowsArgumentException()
    {
        var template = new DnaSequence(RandomDna(100, seed: 22_001));

        var negativeStart = () => PrimerDesigner.DesignPrimers(template, targetStart: -1, targetEnd: 50);
        var endPastEnd = () => PrimerDesigner.DesignPrimers(template, targetStart: 10, targetEnd: 100);
        var startGeEnd = () => PrimerDesigner.DesignPrimers(template, targetStart: 50, targetEnd: 50);

        negativeStart.Should().Throw<ArgumentException>("a negative target start is rejected by the region guard");
        endPastEnd.Should().Throw<ArgumentException>("a target end at/past the template length is rejected");
        startGeEnd.Should().Throw<ArgumentException>("a target start not strictly before the end is rejected");
    }

    #endregion

    #region BE — GC range 0-0 (only all-A·T primers can qualify)

    /// <summary>
    /// BE (GC range 0–0 — KEY): MinGcContent = MaxGcContent = 0 means the GC screen
    /// <c>gc &lt; 0 || gc &gt; 0</c> passes ONLY for a 0% GC (all-A·T) primer. The check is a pair of
    /// comparisons, NOT a normalisation by the (zero-width) range, and CalculateGcContent
    /// null/empty-guards to 0, so there is NO division by zero. An all-A·T template therefore
    /// can still produce a valid all-A·T forward primer whose GC is EXACTLY 0 (inside [0,0]).
    /// We pin: the call does not throw, and any chosen primer's GC is exactly 0. The template is
    /// pure A's flanking the target so an all-A primer is available; Tm/structure screens are
    /// widened so the GC screen is the only discriminator under test.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignPrimers_GcRangeZeroZero_OnlyAllAtPrimersQualifyNoDivByZero()
    {
        // 250 A's: every candidate window is all-A (GC = 0%). Wide Tm/length/structure screens
        // so the GC range 0-0 is the discriminating constraint.
        var template = new DnaSequence(new string('A', 250));
        var param = PrimerDesigner.DefaultParameters with
        {
            MinGcContent = 0,
            MaxGcContent = 0,
            MinTm = 0,
            MaxTm = 200,
            MaxHomopolymer = 1000,   // an all-A primer is one long homopolymer; do not reject on it
            MaxDinucleotideRepeats = 1000,
            Check3PrimeStability = false
        };

        var act = () => PrimerDesigner.DesignPrimers(template, targetStart: 120, targetEnd: 130, param);

        act.Should().NotThrow(
            "a zero-width GC range is a pair of comparisons, not a division; an all-A·T primer satisfies gc<0||gc>0 being false");

        var result = PrimerDesigner.DesignPrimers(template, 120, 130, param);
        result.IsValid.Should().BeTrue("an all-A·T template offers all-A primers whose 0% GC sits inside the [0,0] range");
        result.Forward!.GcContent.Should().Be(0, "the only qualifying GC under a [0,0] range is exactly 0%");
        result.Reverse!.GcContent.Should().Be(0, "the reverse primer's GC is likewise exactly 0% (its reverse complement of all-A is all-T)");
        result.Forward.MeltingTemperature.Should().NotBe(double.NaN, "the Tm of a qualifying primer is finite, never NaN");
    }

    /// <summary>
    /// BE (GC range 0–0, dual): a GC-BEARING flank can satisfy NO primer under a [0,0] GC range —
    /// every window contains at least one G/C, so its GC% exceeds 0 and fails <c>gc &gt; 0</c>. The
    /// theory-correct result is therefore the empty feasible set: an invalid PrimerPairResult,
    /// no crash and no div-by-zero. A random (GC-rich-enough) template makes every candidate fail.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignPrimers_GcRangeZeroZero_GcBearingTemplateYieldsNoPair()
    {
        var template = new DnaSequence(RandomDna(250, seed: 22_002));
        var param = PrimerDesigner.DefaultParameters with
        {
            MinGcContent = 0,
            MaxGcContent = 0,
            MinTm = 0,
            MaxTm = 200,
            Check3PrimeStability = false
        };

        var act = () => PrimerDesigner.DesignPrimers(template, targetStart: 120, targetEnd: 130, param);

        act.Should().NotThrow("a [0,0] GC range over a GC-bearing template empties the feasible set rather than crashing");

        PrimerDesigner.DesignPrimers(template, 120, 130, param)
            .IsValid.Should().BeFalse("every candidate carries a G/C, so its GC% > 0 fails the [0,0] range — no valid pair");
    }

    #endregion

    #region MC — Tm range inverted (contradictory constraint)

    /// <summary>
    /// MC (Tm range inverted — KEY): MinTm > MaxTm is a CONTRADICTORY constraint whose feasible
    /// set is empty. The source does NOT validate the parameter range, so it does NOT throw —
    /// instead the Tm screen <c>tm &lt; MinTm || tm &gt; MaxTm</c> is satisfied by EVERY finite Tm
    /// (any value is either below the high MinTm or above the low MaxTm), so EVERY candidate is
    /// rejected and both best candidates are null. DesignPrimers returns the INVALID
    /// PrimerPairResult (INV-01) — promptly, with no infinite loop and no primer that violates the
    /// impossible range (Primer_Design.md §2.4, §6.1). The template is rich enough that the ONLY
    /// reason no primer qualifies is the inverted Tm range (all other screens are widened).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignPrimers_TmRangeInverted_ReturnsInvalidNeverHangs()
    {
        var template = new DnaSequence(RandomDna(250, seed: 22_003));
        var inverted = PrimerDesigner.DefaultParameters with
        {
            MinTm = 90,   // MinTm > MaxTm: a contradiction
            MaxTm = 10,
            MinGcContent = 0,
            MaxGcContent = 100,
            MaxHomopolymer = 1000,
            MaxDinucleotideRepeats = 1000,
            Check3PrimeStability = false
        };

        var act = () => PrimerDesigner.DesignPrimers(template, targetStart: 120, targetEnd: 130, inverted);

        act.Should().NotThrow(
            "an inverted Tm range is not validated as an error; it simply makes every candidate fail the Tm screen");

        var result = PrimerDesigner.DesignPrimers(template, 120, 130, inverted);
        result.IsValid.Should().BeFalse("no finite Tm can satisfy tm >= 90 AND tm <= 10, so no candidate is valid (INV-01)");
        result.Forward.Should().BeNull("the contradictory Tm range leaves no valid forward primer");
        result.Reverse.Should().BeNull("the contradictory Tm range leaves no valid reverse primer");
    }

    #endregion

    #region Positive sanity — a reasonable template yields a valid pair inside the requested ranges

    /// <summary>
    /// Positive sanity: alongside the degenerate probes, a reasonable template with DEFAULT
    /// parameters must yield a VALID primer pair whose chosen forward and reverse primers fall
    /// INSIDE the requested length, GC and Tm ranges — so the boundary hardening never silently
    /// breaks the core design. A long fixed-seed random template gives both flanks ample,
    /// balanced candidate windows; the default screens (18–25 bp, 40–60% GC, 57–63 °C Tm) are the
    /// contract under test. We pin INV-01 (a valid result has non-null candidates), that each
    /// chosen primer satisfies its own EvaluatePrimer (IsValid, no issues), and INV-03
    /// (ProductSize = reverse.Position + reverse.Length − forward.Position).
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void DesignPrimers_ReasonableTemplate_YieldsValidPairInsideRequestedRanges()
    {
        var param = PrimerDesigner.DefaultParameters;
        // A 600-nt balanced random template with a central target region leaves ~200 bp on each
        // flank for candidate enumeration — deterministic via a fixed seed.
        var template = new DnaSequence(RandomDna(600, seed: 22_000));

        var result = PrimerDesigner.DesignPrimers(template, targetStart: 300, targetEnd: 320, param);

        result.Forward.Should().NotBeNull("a balanced 600-nt template should offer a valid forward candidate");
        result.Reverse.Should().NotBeNull("a balanced 600-nt template should offer a valid reverse candidate");

        var fwd = result.Forward!;
        var rev = result.Reverse!;

        // Each chosen primer is itself valid and inside the requested ranges.
        foreach (var primer in new[] { fwd, rev })
        {
            primer.IsValid.Should().BeTrue("DesignPrimers only selects from valid (issue-free) candidates");
            primer.Issues.Should().BeEmpty("a valid candidate carries no constraint violations");
            primer.Length.Should().BeInRange(param.MinLength, param.MaxLength,
                "the chosen primer's length is inside the requested [MinLength, MaxLength] range");
            primer.GcContent.Should().BeInRange(param.MinGcContent, param.MaxGcContent,
                "the chosen primer's GC% is inside the requested [MinGcContent, MaxGcContent] range");
            primer.MeltingTemperature.Should().BeInRange(param.MinTm, param.MaxTm,
                "the chosen primer's Tm is inside the requested [MinTm, MaxTm] range");
            double.IsNaN(primer.MeltingTemperature).Should().BeFalse("a chosen primer's Tm is finite");
        }

        // INV-03: product size is computed directly from the chosen candidates.
        result.ProductSize.Should().Be(rev.Position + rev.Sequence.Length - fwd.Position,
            "INV-03: ProductSize = reverse.Position + reverse.Length − forward.Position");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PRIMER-STRUCT-001 — primer secondary structure : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PRIMER-STRUCT-001 — primer secondary structure

    #region BE — Extremely short primer (below the stem+loop+stem minimum)

    /// <summary>
    /// BE (extremely short — KEY no-crash boundary): a 1–3 bp primer is far below the
    /// 11-base stem+loop+stem minimum, so HasHairpinPotential short-circuits to FALSE via
    /// its `seq.Length &lt; minStemLength*2 + minLoopLength` guard (PrimerDesigner.cs line
    /// 309) BEFORE any `Substring(i, 4)` slice is taken — too short to fold a stem, and
    /// never an IndexOutOfRangeException from slicing a 4-base fragment out of a 1–3-base
    /// string. Pinned for "", "A", "AC", and "ACG" plus an explicit no-throw guard.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void HasHairpinPotential_ExtremelyShortPrimer_NoStructureNoCrash()
    {
        var act = () => PrimerDesigner.HasHairpinPotential("ACG");

        act.Should().NotThrow(
            "a 3-bp primer is shorter than the 4-base stem fragment; the length guard returns before any Substring slice");

        foreach (var s in new[] { string.Empty, "A", "AC", "ACG" })
            PrimerDesigner.HasHairpinPotential(s).Should().BeFalse(
                $"'{s}' is below the 11-base stem+loop+stem minimum, so no hairpin can form (zero self-structure)");
    }

    /// <summary>
    /// BE (extremely short, dimer surface): two 1–3 bp primers have a 3' overlap shorter
    /// than the minComplementarity default (4), so the per-position complementary count can
    /// never reach the threshold — HasPrimerDimer returns FALSE without crashing. The
    /// `checkLength = min(8, min(len1, len2))` clamp (PrimerDesigner.cs line 407) keeps the
    /// 3'-end Substring in-bounds even for a single base. Pinned: no dimer, no IndexOutOfRange.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void HasPrimerDimer_ExtremelyShortPrimers_NoDimerNoCrash()
    {
        var act = () => PrimerDesigner.HasPrimerDimer("A", "A");

        act.Should().NotThrow(
            "the checkLength clamp keeps the 3'-end slice in-bounds for a single-base primer");

        PrimerDesigner.HasPrimerDimer("A", "A").Should().BeFalse(
            "a 1-base 3' overlap cannot reach the minComplementarity of 4, so no dimer is reported");
        PrimerDesigner.HasPrimerDimer("ACG", "ACG").Should().BeFalse(
            "a 3-base 3' overlap is still below the 4-base complementarity threshold");
    }

    #endregion

    #region BE — Palindromic primer (maximal self-structure — KEY positive signal)

    /// <summary>
    /// BE (palindromic — KEY positive signal): a self-reverse-complementary primer has
    /// MAXIMAL hairpin potential. A clean inverted-repeat stem-loop "GGGGAAACCCC" (4-G stem,
    /// 3-A loop, 4-C stem = revcomp of the 5' stem; exactly the 11-base minimum) and a true
    /// palindrome whose reverse complement equals itself ("GGAATTCCGGAATTCC") both MUST be
    /// flagged: HasHairpinPotential == TRUE. This is the load-bearing assertion of the unit —
    /// the detector must DETECT the self-structure, never silently miss it. Both run the short
    /// O(n²) path (&lt; 100 bp).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void HasHairpinPotential_PalindromicPrimer_DetectsHairpin()
    {
        // Reverse complement of "GGAATTCCGGAATTCC" is itself — a true palindrome.
        DnaSequence.GetReverseComplementString("GGAATTCCGGAATTCC")
            .Should().Be("GGAATTCCGGAATTCC", "the test fixture is a genuine self-reverse-complementary palindrome");

        PrimerDesigner.HasHairpinPotential("GGGGAAACCCC").Should().BeTrue(
            "an inverted-repeat stem-loop (GGGG…loop…CCCC) is the canonical hairpin and MUST be detected");
        PrimerDesigner.HasHairpinPotential("GGAATTCCGGAATTCC").Should().BeTrue(
            "a self-reverse-complementary palindrome has maximal self-complementarity, so a hairpin is detected");
    }

    #endregion

    #region BE — All-G homopolymer (no Watson–Crick self-complement)

    /// <summary>
    /// BE (all-G homopolymer): G does NOT Watson–Crick pair with G, so a poly-G primer has
    /// NO intramolecular self-complement and HasHairpinPotential == FALSE — LOW/zero
    /// self-structure, never a crash. Pinned on BOTH the short O(n²) path (30×G) AND the long
    /// suffix-tree path (150×G, &gt;= 100 bp), so the homopolymer exercises both branches without
    /// fabricating a hairpin or hanging.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void HasHairpinPotential_AllGHomopolymer_NoHairpinOnBothPaths()
    {
        string g30 = new string('G', 30);    // short O(n²) path
        string g150 = new string('G', 150);  // long suffix-tree path

        PrimerDesigner.HasHairpinPotential(g30).Should().BeFalse(
            "G–G does not pair, so a 30-base poly-G primer has no self-complementary stem (no hairpin)");
        PrimerDesigner.HasHairpinPotential(g150).Should().BeFalse(
            "the >=100-bp suffix-tree path agrees: a poly-G homopolymer forms no intramolecular hairpin");
    }

    /// <summary>
    /// BE (all-G, dimer surface — deliberate biochemical contrast): two SEPARATE poly-G
    /// primers DO form an INTER-molecular dimer, because reverse-complementing the second
    /// poly-G yields poly-C, which is fully complementary to the first poly-G's 3' end. So
    /// HasPrimerDimer(G…, G…) is TRUE — and this is correct, NOT a bug: it is duplex
    /// (G:C) pairing between two molecules, not the absent intramolecular self-structure. We
    /// pin it to document that "all-G ⇒ no SELF structure" is a hairpin fact, while the dimer
    /// predicate correctly sees the G:C complementarity across two strands.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void HasPrimerDimer_TwoAllGPrimers_DetectsGcDimerAcrossStrands()
    {
        string g30 = new string('G', 30);

        PrimerDesigner.HasPrimerDimer(g30, g30).Should().BeTrue(
            "revcomp(poly-G) is poly-C, which pairs the first poly-G's 3' end — a correct inter-molecular G:C dimer, distinct from a (absent) hairpin");
    }

    #endregion

    #region BE — Extremely long primer (suffix-tree path; finite, no hang)

    /// <summary>
    /// BE (extremely long): a long primer drives the &gt;= 100-bp SUFFIX-TREE hairpin path,
    /// which must return a FINITE boolean PROMPTLY — no overflow, no Inf/NaN (there is no
    /// floating-point surface), and no hang (pinned with [CancelAfter]). A 5000-bp fixed-seed
    /// random primer returns a well-defined true/false without throwing; an explicit long
    /// stem-loop (60×G · 10×A loop · 60×C, length 130 — a self-complementary inverted repeat)
    /// returns TRUE via the suffix-tree branch, proving the long-primer path is correct, not
    /// merely non-crashing.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void HasHairpinPotential_ExtremelyLongPrimer_FiniteResultNoHang()
    {
        string longRandom = RandomDna(5000, seed: 23_001);
        string longStemLoop = new string('G', 60) + new string('A', 10) + new string('C', 60); // length 130

        var act = () => PrimerDesigner.HasHairpinPotential(longRandom);
        act.Should().NotThrow(
            "a 5000-bp primer runs the suffix-tree path; it must terminate with a boolean, never overflow or hang");

        PrimerDesigner.HasHairpinPotential(longStemLoop).Should().BeTrue(
            "a long explicit inverted-repeat stem-loop (60×G…60×C) is a hairpin and the >=100-bp suffix-tree path detects it");
    }

    #endregion

    #region Positive sanity — an ordinary primer flows through both structure screens

    /// <summary>
    /// Positive sanity: an ordinary, well-behaved primer must flow through BOTH structure
    /// predicates without throwing, returning plain booleans — so the boundary hardening never
    /// breaks the core screens. A balanced 20-nt primer with no designed inverted repeat is not
    /// flagged as a hairpin and does not self-dimer with a non-complementary partner, while the
    /// designed inverted repeat IS flagged — pinning that the detector discriminates structure
    /// from non-structure rather than answering a constant.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void StructureScreens_OrdinaryPrimer_ReturnWellDefinedBooleans()
    {
        // A poly-A 20-mer has no inverted repeat (A does not pair A): no hairpin. For the
        // dimer screen we pair it with a poly-C partner: revcomp(poly-C) is poly-G, and A:G
        // does not pair, so the 3' ends are non-complementary and no dimer is reported.
        string ordinary = new string('A', 20);    // 20-nt, no self-complement
        const string hairpinPrimer = "GGGGAAACCCC"; // designed stem-loop

        var hairpin = () => PrimerDesigner.HasHairpinPotential(ordinary);
        var dimer = () => PrimerDesigner.HasPrimerDimer(ordinary, new string('C', 20));

        hairpin.Should().NotThrow("an ordinary 20-nt primer is screened for hairpins without throwing");
        dimer.Should().NotThrow("an ordinary 20-nt primer is screened for dimers without throwing");

        PrimerDesigner.HasHairpinPotential(ordinary).Should().BeFalse(
            "a poly-A primer has no self-complementary stem, so it forms no hairpin");
        PrimerDesigner.HasPrimerDimer(ordinary, new string('C', 20)).Should().BeFalse(
            "a poly-A vs poly-C pair has non-complementary 3' ends (A:G after reverse-complementing poly-C), so no dimer");
        PrimerDesigner.HasHairpinPotential(hairpinPrimer).Should().BeTrue(
            "the designed inverted repeat IS detected — the screen discriminates structure from non-structure");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROBE-DESIGN-001 — hybridization probe design : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROBE-DESIGN-001 — probe design

    #region BE — Target shorter than the minimum probe length

    /// <summary>
    /// BE (seq shorter than min probe): a target SHORTER than the configured MinLength
    /// (Microarray default MinLength = 50) cannot hold even one probe window. The
    /// length short-circuit `targetSequence.Length &lt; param.MinLength` (ProbeDesigner.cs
    /// line 134) fires `yield break` BEFORE the GC prefix sums or the candidate-window
    /// loop, so no Substring is ever taken past the end — an empty result, never an
    /// IndexOutOfRangeException (Hybridization_Probe_Design.md §6.1 "Sequence shorter
    /// than MinLength → no probes"). DesignProbes is a yield iterator, so enumeration is
    /// forced. Pinned with the default params and with an explicit short window.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignProbes_TargetShorterThanMinLength_IsEmptyAndDoesNotThrow()
    {
        // 30-nt target, but Microarray MinLength is 50 — no window fits.
        string shortTarget = RandomDna(30, seed: 24_001);

        var act = () => ProbeDesigner.DesignProbes(shortTarget).ToList();

        act.Should().NotThrow(
            "a target shorter than MinLength trips the length short-circuit before any windowing; no Substring is taken past the end");

        ProbeDesigner.DesignProbes(shortTarget).Should().BeEmpty(
            "no probe window of MinLength can fit in a sub-MinLength target, so no probe is produced");
    }

    /// <summary>
    /// BE: the boundary is exact — a target ONE base shorter than MinLength yields no
    /// probe, while a target of EXACTLY MinLength can yield one. Using a custom param
    /// set (MinLength = MaxLength = 20, qPCR-style) isolates the length boundary from
    /// the other screens: a 19-nt target is empty; a 20-nt target with in-range
    /// composition produces at least one 20-nt probe.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignProbes_TargetAtMinLengthBoundary_EmptyBelowNonEmptyAt()
    {
        var param = new ProbeDesigner.ProbeParameters(
            MinLength: 20, MaxLength: 20,
            MinTm: 0, MaxTm: 200,           // permissive Tm so the length boundary is isolated
            MinGc: 0.0, MaxGc: 1.0,         // permissive GC
            MaxHomopolymer: 20,
            AvoidSecondaryStructure: false,
            MaxSelfComplementarity: 1.0);

        // A balanced 20-mer (10 GC) so it is a clean candidate; the 19-mer is one short.
        const string at20 = "ACGTACGTACGTACGTACGT"; // length 20
        string below19 = at20.Substring(0, 19);     // length 19

        ProbeDesigner.DesignProbes(below19, param).Should().BeEmpty(
            "a target one base shorter than MinLength has no candidate window");
        ProbeDesigner.DesignProbes(at20, param).Should().NotBeEmpty(
            "a target of exactly MinLength holds one window, which under permissive screens is a valid probe");
    }

    #endregion

    #region MC — Inverted (contradictory) Tm range

    /// <summary>
    /// MC (Tm range inverted, KEY): MinTm &gt; MaxTm is a CONTRADICTION whose feasible Tm
    /// set is empty. UNLIKE PrimerDesigner, ProbeDesigner does NOT hard-reject a
    /// candidate when a screen fails — the Tm screen `tm &lt; MinTm || tm &gt; MaxTm`
    /// (ProbeDesigner.cs line 267) is satisfied by EVERY finite Tm under an inverted
    /// range, so EVERY candidate incurs the fixed −0.3 Tm penalty and records a
    /// "Tm … outside range" warning, yet survives if its residual score stays positive.
    /// So the theory-correct contract is NOT an empty result: it is "no throw, no hang,
    /// and every returned probe carries the Tm-out-of-range warning". The range is never
    /// validated and never throws. Pinned on a real 60-nt target with the Microarray
    /// length/GC defaults but an inverted Tm window (MinTm 95 &gt; MaxTm 5).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignProbes_InvertedTmRange_NoThrowEveryProbeWarnsTmOutOfRange()
    {
        var param = ProbeDesigner.Defaults.Microarray with { MinTm = 95.0, MaxTm = 5.0 };
        string target = RandomDna(120, seed: 24_002);

        var act = () => ProbeDesigner.DesignProbes(target, param).ToList();

        act.Should().NotThrow(
            "an inverted Tm range is a contradictory constraint; the source does not validate it, so it neither throws nor hangs");

        var probes = ProbeDesigner.DesignProbes(target, param).ToList();
        probes.Should().OnlyContain(
            p => p.Warnings.Any(w => w.Contains("Tm") && w.Contains("outside range")),
            "an inverted Tm range fails the Tm screen for EVERY candidate, so every returned probe records the Tm-out-of-range warning rather than being silently accepted");
        probes.Should().OnlyContain(p => p.Score > 0,
            "INV-01: only raw-score-positive candidates are retained, even when the Tm screen always fails");
    }

    #endregion

    #region BE — Empty / null target sequence

    /// <summary>
    /// BE (empty seq): the empty target is the lower size boundary.
    /// `string.IsNullOrEmpty(targetSequence)` short-circuits to `yield break`
    /// (ProbeDesigner.cs line 134) BEFORE any GC prefix sum or windowing — a defined
    /// degenerate boundary, never a crash, division, or out-of-range index
    /// (Hybridization_Probe_Design.md §6.1 "Null or empty target → no probes"). Pinned
    /// for both the empty string and the null reference (the same IsNullOrEmpty guard
    /// catches both, so null yields the empty result rather than a NullReferenceException).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DesignProbes_EmptyOrNullTarget_IsEmptyAndDoesNotThrow()
    {
        var empty = () => ProbeDesigner.DesignProbes(string.Empty).ToList();
        var nullTarget = () => ProbeDesigner.DesignProbes((string)null!).ToList();

        empty.Should().NotThrow("an empty target short-circuits to an empty result before any windowing");
        nullTarget.Should().NotThrow(
            "IsNullOrEmpty catches the null target, so it yields the empty result rather than dereferencing into a NullReferenceException");

        ProbeDesigner.DesignProbes(string.Empty).Should().BeEmpty();
        ProbeDesigner.DesignProbes((string)null!).Should().BeEmpty();
    }

    #endregion

    #region Positive sanity — a reasonable target yields an in-range valid probe

    /// <summary>
    /// Positive sanity: alongside the degenerate probes, a reasonable target must yield at
    /// least one VALID probe whose GC and Tm fall INSIDE the requested ranges — so the
    /// boundary hardening never silently breaks the core function. The parameter ranges are
    /// chosen MUTUALLY SATISFIABLE for the salt-adjusted Tm model: for a 50–60-nt probe the
    /// source computes Tm = 81.5 + 16.6·log10(0.05) + 41·GC − 600/length, so a 0.40–0.60-GC
    /// window maps to roughly 64–74 °C — hence the GC window [0.40, 0.60] is paired with a Tm
    /// window [60, 80] (the Microarray Tm default of 75–85 is, by contrast, unreachable at
    /// ≤ 0.60 GC for this formula, which is exactly why the inverted-range test asserts on
    /// warnings rather than the default ranges). The target is a 300-nt fixed-seed sequence;
    /// AT LEAST ONE returned probe must satisfy INV-01 (Score &gt; 0), INV-02
    /// (0 ≤ GcContent ≤ 1), have GC within [MinGc, MaxGc] AND Tm within [MinTm, MaxTm], and —
    /// being fully in-range — carry NO GC or Tm out-of-range warning.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void DesignProbes_ReasonableTarget_YieldsInRangeValidProbe()
    {
        // GC and Tm windows are mutually satisfiable for the 50–60-nt salt-adjusted Tm model.
        var param = new ProbeDesigner.ProbeParameters(
            MinLength: 50, MaxLength: 60,
            MinTm: 60, MaxTm: 80,
            MinGc: 0.40, MaxGc: 0.60,
            MaxHomopolymer: 5,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.3);
        string target = RandomDna(300, seed: 24_003);

        var probes = ProbeDesigner.DesignProbes(target, param, maxProbes: 20).ToList();

        probes.Should().NotBeEmpty(
            "a 300-nt random target overwhelmingly contains a 50–60-nt window meeting the (mutually satisfiable) screens");

        // INV-01 / INV-02 hold for EVERY returned probe, and each is a well-formed in-bounds window.
        probes.Should().OnlyContain(p =>
                p.Score > 0 &&                                  // INV-01
                p.GcContent >= 0.0 && p.GcContent <= 1.0 &&     // INV-02
                p.Sequence.Length >= param.MinLength && p.Sequence.Length <= param.MaxLength &&
                p.Start >= 0 && p.End == p.Start + p.Sequence.Length - 1 && p.End < target.Length,
            "every returned probe is raw-score-positive (INV-01), has GC as a fraction (INV-02), and is an in-bounds window of the configured length");

        // At least one returned probe is FULLY in-range and therefore warning-free on GC and Tm.
        var inRange = probes.Where(p =>
                p.GcContent >= param.MinGc && p.GcContent <= param.MaxGc &&
                p.Tm >= param.MinTm && p.Tm <= param.MaxTm)
            .ToList();

        inRange.Should().NotBeEmpty(
            "with mutually satisfiable GC and Tm windows, a 300-nt random target yields at least one probe whose GC and Tm both fall inside the requested ranges");

        var best = inRange.First();
        best.Warnings.Should().NotContain(w => w.Contains("GC content") && w.Contains("outside range"),
            "a fully in-range probe records no GC-out-of-range warning");
        best.Warnings.Should().NotContain(w => w.Contains("Tm") && w.Contains("outside range"),
            "a fully in-range probe records no Tm-out-of-range warning");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROBE-VALID-001 — probe validation : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROBE-VALID-001 — probe validation

    #region MC/INJ — Non-DNA probe (junk / special chars / unicode)

    /// <summary>
    /// MC/INJ: a probe of pure non-DNA junk (symbols, digits, spaces, unicode, null
    /// bytes) must NOT crash and must NOT produce NaN. ValidateProbe does not validate
    /// the alphabet — the junk probe is upper-cased and scanned literally against the
    /// (clean A/C/G/T) reference. FindApproximateMatches compares with a plain `!=`, so
    /// junk simply mismatches every reference base; the reverse complement used by
    /// self-complementarity passes non-IUPAC chars through unchanged (GetComplementBase
    /// fall-through). Because the probe length is &gt; 0, the self-complementarity division
    /// is never by zero — the result is a well-formed ProbeValidation with finite,
    /// in-range SpecificityScore and SelfComplementarity (INV-01, INV-02, INV-03).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ValidateProbe_NonDnaJunkProbe_HandledWithFiniteInRangeResult()
    {
        string junk = "$#@!12 3￿\0xyzé❤";
        var references = new[] { RandomDna(500, seed: 25_001) };

        var act = () => ProbeDesigner.ValidateProbe(junk, references);
        act.Should().NotThrow(
            "a non-DNA probe is scanned literally, not validated; junk mismatches references and never indexes out of range");

        var result = ProbeDesigner.ValidateProbe(junk, references);

        result.SpecificityScore.Should().BeInRange(0.0, 1.0,
            "INV-01: the specificity score is mapped to [0,1] regardless of the probe alphabet");
        result.SelfComplementarity.Should().BeInRange(0.0, 1.0,
            "INV-02: self-complementarity is a fraction of positions in [0,1] — finite, never NaN");
        double.IsNaN(result.SpecificityScore).Should().BeFalse("a non-empty junk probe never divides by zero");
        double.IsNaN(result.SelfComplementarity).Should().BeFalse("self-complementarity divides by the non-zero probe length");
        result.OffTargetHits.Should().BeGreaterThanOrEqualTo(0, "INV-03: hit counts are accumulated, never negative");
    }

    /// <summary>
    /// INJ: special / unicode / null-byte characters INTERLEAVED into an otherwise valid
    /// probe must not crash and must not inflate the hit count — an injected junk base
    /// fails the literal `!=` comparison against any reference base, so it can only ever
    /// REDUCE matches, never fabricate an off-target. The result stays finite and in range.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ValidateProbe_JunkInterleavedIntoProbe_NoCrashNoInflatedHits()
    {
        string probe = "ACGTACGT\0❤#ACGTACGT";
        var references = new[] { "ACGTACGTACGTACGTACGT", RandomDna(300, seed: 25_002) };

        var act = () => ProbeDesigner.ValidateProbe(probe, references);
        act.Should().NotThrow(
            "injected null bytes / unicode are compared literally and never index out of range");

        var result = ProbeDesigner.ValidateProbe(probe, references);
        result.SpecificityScore.Should().BeInRange(0.0, 1.0);
        result.SelfComplementarity.Should().BeInRange(0.0, 1.0);
        result.OffTargetHits.Should().BeGreaterThanOrEqualTo(0,
            "INV-03: injected junk can only suppress a match, never invent an off-target");
    }

    #endregion

    #region BE/INJ — Extremely short probe

    /// <summary>
    /// BE: an extremely short (1–3 bp) probe is NOT rejected — ValidateProbe enforces no
    /// minimum length (only the EMPTY probe is special-cased). A tiny probe is scanned:
    /// FindApproximateMatches' bound `i &lt;= text.Length − pattern.Length` runs over each
    /// reference (a 1-nt probe matches many positions), the hairpin screen's negative loop
    /// bound short-circuits to false, and self-complementarity divides by the non-zero
    /// length. The result is therefore finite and in-range — handled and scored, never a
    /// crash and never a NaN (only the empty probe yields the structured invalid result).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ValidateProbe_ExtremelyShortProbe_ScoredNotRejected()
    {
        var references = new[] { RandomDna(400, seed: 25_003) };

        foreach (var probe in new[] { "A", "AC", "ACG" })
        {
            var act = () => ProbeDesigner.ValidateProbe(probe, references);
            act.Should().NotThrow($"a {probe.Length}-nt probe has no minimum-length screen and is scanned, not rejected");

            var result = ProbeDesigner.ValidateProbe(probe, references);
            result.SpecificityScore.Should().BeInRange(0.0, 1.0, "INV-01 holds for any non-empty probe length");
            result.SelfComplementarity.Should().BeInRange(0.0, 1.0, "INV-02: a 1–3 nt probe never divides by zero");
            double.IsNaN(result.SelfComplementarity).Should().BeFalse("non-zero length ⇒ no NaN");
            result.HasSecondaryStructure.Should().BeFalse(
                "a probe shorter than 2·stem+loop (11 nt) cannot form a hairpin; the screen short-circuits");
        }
    }

    #endregion

    #region INJ — Null probe and null references

    /// <summary>
    /// INJ (KEY): a null probe is the documented ArgumentNullException (ThrowIfNull,
    /// line 497) — raised eagerly, NEVER a NullReferenceException from dereferencing
    /// `probeSequence.ToUpperInvariant()`. The sibling null-reference-collection guard
    /// (line 498) is pinned alongside it.
    /// </summary>
    [Test]
    public void ValidateProbe_NullProbeOrNullReferences_ThrowsArgumentNullException()
    {
        var references = new[] { "ACGTACGTACGT" };

        var nullProbe = () => ProbeDesigner.ValidateProbe(null!, references);
        var nullRefs = () => ProbeDesigner.ValidateProbe("ACGTACGT", null!);

        nullProbe.Should().Throw<ArgumentNullException>(
            "a null probe is null-guarded eagerly, never dereferenced into a NullReferenceException");
        nullRefs.Should().Throw<ArgumentNullException>(
            "a null reference collection is null-guarded the same way");
    }

    /// <summary>
    /// BE: the empty probe is the lower length boundary — it is a DEFINED degenerate input
    /// returning a structured invalid result ("Empty probe sequence" issue, SpecificityScore
    /// 0.0, OffTargetHits 0), NOT a throw and NOT a crash (Probe_Validation.md §6.1).
    /// </summary>
    [Test]
    public void ValidateProbe_EmptyProbe_ReturnsStructuredInvalidResult()
    {
        var references = new[] { "ACGTACGTACGT" };

        var act = () => ProbeDesigner.ValidateProbe(string.Empty, references);
        act.Should().NotThrow("the empty probe is a defined degenerate case, not an exception");

        var result = ProbeDesigner.ValidateProbe(string.Empty, references);
        result.IsValid.Should().BeFalse("an empty probe cannot hybridize specifically");
        result.SpecificityScore.Should().Be(0.0, "the empty-probe structured result fixes specificity at 0.0");
        result.OffTargetHits.Should().Be(0);
        result.Issues.Should().Contain("Empty probe sequence");
    }

    #endregion

    #region MC — Probe with N's (degenerate IUPAC bases)

    /// <summary>
    /// MC (KEY): a probe of all-N (the IUPAC "any base" code) must be HANDLED, not crash,
    /// and must produce NO NaN. ValidateProbe does NOT treat 'N' as a matching wildcard —
    /// FindApproximateMatches is a literal `!=`, so an 'N' mismatches an A/C/G/T reference
    /// base and 'N' is never counted as G/C. In self-complementarity 'N' complements to 'N'
    /// (GetComplementBase, line 155), so an all-N probe is maximally self-complementary
    /// (every position N==N → 1.0) — a DEFINED value, never a NaN. The result is finite and
    /// in-range on every numeric field.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ValidateProbe_AllNProbe_HandledNoNaN()
    {
        string probe = new string('N', 20);
        var references = new[] { RandomDna(500, seed: 25_004) };

        var act = () => ProbeDesigner.ValidateProbe(probe, references);
        act.Should().NotThrow("'N' is not a wildcard here; it is compared literally and complemented to 'N' — no crash");

        var result = ProbeDesigner.ValidateProbe(probe, references);

        double.IsNaN(result.SpecificityScore).Should().BeFalse("a non-empty all-N probe never divides by zero");
        double.IsNaN(result.SelfComplementarity).Should().BeFalse("N complements to N; self-complementarity is a finite fraction");
        result.SpecificityScore.Should().BeInRange(0.0, 1.0, "INV-01");
        result.SelfComplementarity.Should().BeInRange(0.0, 1.0, "INV-02");
        result.SelfComplementarity.Should().Be(1.0,
            "every position of an all-N probe equals its own complement (N↔N), so it is maximally self-complementary");
        result.OffTargetHits.Should().BeGreaterThanOrEqualTo(0, "INV-03");
    }

    /// <summary>
    /// MC: a probe with a FEW interspersed N's against a reference that contains that exact
    /// probe verbatim must still match within the mismatch tolerance — the N's are ordinary
    /// (non-wildcard) characters that match themselves. No crash, finite result, hits &gt;= 1.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ValidateProbe_FewNsMatchingReference_HandledFinite()
    {
        string probe = "ACGTNCGTACGTNCGTACGT"; // 20 nt with two N's
        var references = new[] { "TTTT" + probe + "TTTT" }; // probe present verbatim

        var result = ProbeDesigner.ValidateProbe(probe, references);

        double.IsNaN(result.SpecificityScore).Should().BeFalse();
        double.IsNaN(result.SelfComplementarity).Should().BeFalse();
        result.OffTargetHits.Should().BeGreaterThanOrEqualTo(1,
            "the probe occurs verbatim (its N's match themselves), so at least one approximate hit is counted");
        result.SpecificityScore.Should().BeInRange(0.0, 1.0, "INV-01");
    }

    #endregion

    #region Positive sanity — a clean unique probe validates

    /// <summary>
    /// Positive sanity: a clean A/C/G/T probe that occurs EXACTLY ONCE in the references and
    /// carries no self-complementarity or hairpin issues must validate as IsValid == true with
    /// SpecificityScore 1.0, OffTargetHits 1, and NO issues — so the fuzz hardening never
    /// silently breaks the core success path. INV-04 (OffTargetHits == 1 ⇒ SpecificityScore
    /// == 1.0) is pinned directly.
    /// </summary>
    [Test]
    public void ValidateProbe_CleanUniqueProbe_IsValidWithUniqueSpecificity()
    {
        // A non-self-complementary probe (self-complementarity 0.0) present exactly once.
        string probe = "GTACGGATCCATGCTAACGT"; // 20 nt, reverse complement differs at every position
        string reference = "TTTTTTTTTT" + probe + "TTTTTTTTTT";
        var references = new[] { reference };

        var result = ProbeDesigner.ValidateProbe(probe, references, maxMismatches: 0);

        result.OffTargetHits.Should().Be(1, "the probe occurs exactly once at mismatch tolerance 0");
        result.SpecificityScore.Should().Be(1.0, "INV-04: a single hit maps to specificity 1.0");
        result.IsValid.Should().BeTrue("a unique, low-self-complementarity probe is a valid probe");
        result.Issues.Should().BeEmpty("a clean unique probe records no validation issues");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  RESTR-FIND-001 — restriction site finding : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region RESTR-FIND-001 — restriction site finding

    #region BE — Empty enzyme list

    /// <summary>
    /// BE: the multi-enzyme overload called with an EMPTY enzyme list is the lower
    /// boundary of "how many enzymes to scan". The `foreach (var name in enzymeNames)`
    /// loop (RestrictionAnalyzer.cs line 211) iterates zero times, so the result is the
    /// EMPTY enumerable — never a crash and never a fabricated site. Because the overload
    /// is a `yield` iterator, the body runs only on enumeration, so we force `.ToList()`.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindSites_EmptyEnzymeList_IsEmptyAndDoesNotThrow()
    {
        var sequence = new DnaSequence("AAAGAATTCAAAGGATCCAAA");

        var act = () => RestrictionAnalyzer.FindSites(sequence, Array.Empty<string>()).ToList();

        act.Should().NotThrow(
            "an empty enzyme list scans no enzymes; the per-enzyme loop never runs, so there is nothing to throw");
        RestrictionAnalyzer.FindSites(sequence, Array.Empty<string>())
            .Should().BeEmpty("no enzymes requested ⇒ no restriction sites, never a spurious cut");
    }

    #endregion

    #region MC — Unknown enzyme names

    /// <summary>
    /// MC (KEY): an enzyme NAME absent from the built-in catalog is the central
    /// "invalid specification" probe. The lookup uses TryGetValue (GetEnzyme, line 76–77)
    /// and the `?? throw` maps a miss to a DELIBERATE, typed ArgumentException
    /// "Unknown enzyme" (line 133–134) — it must NEVER let a raw KeyNotFoundException
    /// leak from the dictionary. The typed overload returns FindSitesCore's enumerable
    /// directly, so the throw fires EAGERLY at the call; we also pin that it is NOT a
    /// KeyNotFoundException (Restriction_Site_Detection.md §3.3, §6.1).
    /// </summary>
    [Test]
    public void FindSites_UnknownEnzymeName_TypedOverload_ThrowsArgumentException()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");

        var act = () => RestrictionAnalyzer.FindSites(sequence, "NotARealEnzyme");

        act.Should().Throw<ArgumentException>(
                "an unknown enzyme name is rejected deliberately rather than scanned with a garbage motif")
            .And.Message.Should().Contain("Unknown enzyme",
                "the rejection is the intentional 'Unknown enzyme' validation, not an internal lookup miss");
        act.Should().NotThrow<KeyNotFoundException>(
            "a raw KeyNotFoundException must never escape the dictionary lookup; GetEnzyme uses TryGetValue");
    }

    /// <summary>
    /// MC (KEY): the raw-string overload resolves the SAME catalog lookup, but inside a
    /// `yield` iterator (line 147–148), so an unknown enzyme name surfaces its
    /// ArgumentException only on ENUMERATION. We force `.ToList()` so the documented throw
    /// actually fires, and pin that it is the intentional "Unknown enzyme" rejection and
    /// not a leaked KeyNotFoundException.
    /// </summary>
    [Test]
    public void FindSites_UnknownEnzymeName_RawOverload_ThrowsArgumentExceptionOnEnumeration()
    {
        var act = () => RestrictionAnalyzer.FindSites("AAAGAATTCAAA", "Bogus123").ToList();

        act.Should().Throw<ArgumentException>(
                "the raw-string overload rejects an unknown enzyme name on enumeration")
            .And.Message.Should().Contain("Unknown enzyme");
        act.Should().NotThrow<KeyNotFoundException>(
            "the unknown-name rejection is a typed ArgumentException, never a raw KeyNotFoundException");
    }

    /// <summary>
    /// MC/BE: a null OR empty enzyme name is the degenerate "no enzyme specified" input.
    /// The typed overload guards it with ArgumentNullException (line 130–131) BEFORE the
    /// catalog lookup, so an empty name is a documented validation throw — never a
    /// silent empty scan and never a KeyNotFoundException.
    /// </summary>
    [Test]
    public void FindSites_NullOrEmptyEnzymeName_TypedOverload_ThrowsArgumentNullException()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");

        var nullName = () => RestrictionAnalyzer.FindSites(sequence, (string)null!);
        var emptyName = () => RestrictionAnalyzer.FindSites(sequence, string.Empty);

        nullName.Should().Throw<ArgumentNullException>(
            "a null enzyme name is rejected up front, not dereferenced into a lookup");
        emptyName.Should().Throw<ArgumentNullException>(
            "an empty enzyme name is the same degenerate 'no enzyme' input and is rejected up front");
    }

    #endregion

    #region MC — Non-DNA characters in the sequence

    /// <summary>
    /// MC: non-DNA junk fed to the TYPED overload is rejected up front — the DnaSequence
    /// constructor validates A/C/G/T and throws ArgumentException on the first offending
    /// character, so junk never reaches the restriction scanner at all.
    /// </summary>
    [Test]
    public void FindSites_NonDnaSequence_TypedOverload_RejectedAtConstruction()
    {
        var act = () => new DnaSequence("GAATTC$#@!XYZ123");

        act.Should().Throw<ArgumentException>(
            "the validated DnaSequence type rejects any non-A/C/G/T character, so the typed scanner only ever sees clean DNA");
    }

    /// <summary>
    /// MC: non-DNA junk fed to the RAW-string overload must NOT crash and must NOT invent
    /// sites. Each sequence character is tested as the `nucleotide` argument of
    /// MatchesIupac(seqChar, patChar); a junk char satisfies no enzyme-pattern code, so it
    /// can never complete a recognition match. The reverse-strand pass complements via
    /// GetReverseComplementString, whose fall-through passes non-IUPAC chars through
    /// unchanged — no exception, no out-of-range indexing. Pure-junk input yields no
    /// sites. We force enumeration so the in-iterator scan actually runs.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindSites_NonDnaSequence_RawOverload_NoSitesNoCrash()
    {
        const string junk = "$#@!XYZ123 ￿qwerty";

        var act = () => RestrictionAnalyzer.FindSites(junk, "EcoRI").ToList();

        act.Should().NotThrow(
            "non-DNA characters are tested against the enzyme pattern via IUPAC matching and complemented via a pass-through arm; neither path indexes out of range");
        RestrictionAnalyzer.FindSites(junk, "EcoRI")
            .Should().BeEmpty(
                "a character that is not A/C/G/T matches no recognition position, so pure junk can never produce a spurious restriction site");
    }

    /// <summary>
    /// MC: junk INTERLEAVED into an otherwise-valid GAATTC window must suppress that site,
    /// not crash. "GAATTC" at offset 3 is a real EcoRI site; replacing the final 'C' with
    /// '#' breaks the recognition match (the '#' satisfies no pattern code), so no forward
    /// EcoRI site is reported there — proving junk neither fabricates nor silently
    /// "rounds to" a valid recognition site.
    /// </summary>
    [Test]
    public void FindSites_JunkInsideRecognitionWindow_SuppressesThatSite()
    {
        // "AAA" + "GAATT#" + "AAA": the EcoRI window at index 3 is broken by '#'.
        const string broken = "AAAGAATT#AAA";

        var sites = RestrictionAnalyzer.FindSites(broken, "EcoRI").ToList();

        sites.Should().NotContain(s => s.IsForwardStrand && s.Position == 3,
            "the '#' in the sixth recognition position fails IUPAC matching, so the forward EcoRI site at index 3 is correctly NOT reported");
        sites.Should().OnlyContain(s => s.RecognizedSequence.All(c => c == 'A' || c == 'C' || c == 'G' || c == 'T'),
            "any site that IS reported still has a clean A/C/G/T recognition string — junk never leaks into a reported site");
    }

    #endregion

    #region BE — Empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The raw-string overload
    /// short-circuits null/empty to the empty enumerable (yield break, INV-04); the typed
    /// overload over an empty DnaSequence has a negative forward-scan bound
    /// `i <= 0 − patternLen`, so neither the forward nor the reverse loop runs. Neither
    /// path divides, indexes past the end, or hangs (Restriction_Site_Detection.md §6.1).
    /// Pinned for both surfaces plus the raw null input. NOTE: the raw-string overload's
    /// empty short-circuit precedes the enzyme lookup, so an empty sequence with an UNKNOWN
    /// enzyme is still the empty result, not the unknown-enzyme throw.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindSites_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var typed = () => RestrictionAnalyzer.FindSites(new DnaSequence(string.Empty), "EcoRI").ToList();
        var rawEmpty = () => RestrictionAnalyzer.FindSites(string.Empty, "EcoRI").ToList();
        var rawNull = () => RestrictionAnalyzer.FindSites((string)null!, "EcoRI").ToList();

        typed.Should().NotThrow("an empty sequence has no scan window; the forward bound is negative so the loop never runs");
        rawEmpty.Should().NotThrow("the raw-string overload short-circuits empty input to an empty result");
        rawNull.Should().NotThrow("the raw-string overload treats null input as empty, not as an error");

        RestrictionAnalyzer.FindSites(new DnaSequence(string.Empty), "EcoRI").Should().BeEmpty();
        RestrictionAnalyzer.FindSites(string.Empty, "EcoRI").Should().BeEmpty();
        RestrictionAnalyzer.FindSites((string)null!, "EcoRI").Should().BeEmpty();
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no typed input". The typed overload
    /// guards it with ArgumentNullException (ThrowIfNull, line 129), raised eagerly at the
    /// call — never a NullReferenceException dereferencing `sequence.Sequence`.
    /// </summary>
    [Test]
    public void FindSites_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => RestrictionAnalyzer.FindSites((DnaSequence)null!, "EcoRI");

        act.Should().Throw<ArgumentNullException>(
            "the typed overload null-guards its sequence; null is rejected, never dereferenced into a NullReferenceException");
    }

    #endregion

    #region BE — Recognition site longer than the sequence

    /// <summary>
    /// BE: a sequence SHORTER than the recognition site is the degenerate "site longer
    /// than seq" case. NotI's site is 8 nt (GCGGCCGC); on a 3-nt sequence the scan bounds
    /// `i <= len − patternLen` are negative on both strands, so neither loop runs and no
    /// Substring is taken past the end — an empty result, never an IndexOutOfRangeException
    /// (Restriction_Site_Detection.md §6.1). Pinned on both surfaces.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindSites_SiteLongerThanSequence_IsEmptyAndDoesNotThrow()
    {
        var typed = () => RestrictionAnalyzer.FindSites(new DnaSequence("ACG"), "NotI").ToList();
        var raw = () => RestrictionAnalyzer.FindSites("ACG", "NotI").ToList();

        typed.Should().NotThrow(
            "the 8-nt NotI site cannot fit in a 3-nt sequence; the scan bound is negative so no Substring is taken past the end");
        raw.Should().NotThrow("the raw-string overload is equally guarded against indexing past the sequence end");

        RestrictionAnalyzer.FindSites(new DnaSequence("ACG"), "NotI").Should().BeEmpty(
            "a recognition site longer than the whole sequence yields no sites, not a crash");
        RestrictionAnalyzer.FindSites("ACG", "NotI").Should().BeEmpty();
    }

    #endregion

    #region Positive sanity — known EcoRI sites are found at the correct positions

    /// <summary>
    /// Positive sanity: alongside the degenerate probes, a textbook EcoRI forward site
    /// must be detected at the CORRECT position with the CORRECT cut position and matched
    /// bases, so the boundary hardening never silently breaks the core function. EcoRI =
    /// GAATTC, CutPositionForward = 1; in "AAAGAATTCAAA" the GAATTC window sits at index 3,
    /// so the forward site has Position 3, RecognizedSequence "GAATTC", and CutPosition
    /// 3 + 1 = 4 (INV-01, INV-02).
    /// </summary>
    [Test]
    public void FindSites_KnownEcoRISite_DetectedAtCorrectPosition()
    {
        const string seq = "AAAGAATTCAAA"; // GAATTC at index 3

        var forward = RestrictionAnalyzer.FindSites(seq, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        var site = forward.Should().ContainSingle(s => s.Position == 3).Subject;
        site.RecognizedSequence.Should().Be("GAATTC", "INV-02: the matched bases are exactly the EcoRI recognition sequence");
        site.RecognizedSequence.Length.Should().Be(site.Enzyme.RecognitionLength, "INV-02: matched length == recognition length");
        site.CutPosition.Should().Be(4, "EcoRI cuts after the leading G: CutPosition = Position 3 + CutPositionForward 1");
        site.Position.Should().BeInRange(0, seq.Length - 6, "INV-01: the site start is a valid window start");
        site.Enzyme.Name.Should().Be("EcoRI");
    }

    /// <summary>
    /// Positive sanity / INV-03: enzyme-name lookup is case-insensitive (the catalog uses
    /// StringComparer.OrdinalIgnoreCase). The same GAATTC site must be found whether the
    /// caller writes "EcoRI", "ecori", or "ECORI" — so a fuzzed-case name resolves to the
    /// real enzyme rather than falling through to the unknown-enzyme throw.
    /// </summary>
    [Test]
    public void FindSites_EnzymeNameLookup_IsCaseInsensitive()
    {
        const string seq = "AAAGAATTCAAA";

        var lower = RestrictionAnalyzer.FindSites(seq, "ecori").Where(s => s.IsForwardStrand).ToList();
        var upper = RestrictionAnalyzer.FindSites(seq, "ECORI").Where(s => s.IsForwardStrand).ToList();

        lower.Should().ContainSingle(s => s.Position == 3,
            "INV-03: 'ecori' resolves to EcoRI case-insensitively and finds the GAATTC site");
        upper.Should().ContainSingle(s => s.Position == 3,
            "INV-03: 'ECORI' resolves to EcoRI case-insensitively and finds the GAATTC site");
    }

    /// <summary>
    /// Positive sanity: an IUPAC-degenerate enzyme exercises the ambiguity-matching path so
    /// the boundary work does not break it. HincII = GTYRAC (Y = C/T, R = A/G). The window
    /// "GTCGAC" satisfies GTYRAC (Y matches C, R matches G), so a forward HincII site is
    /// found — proving the degenerate matcher is reached and correct, not bypassed.
    /// </summary>
    [Test]
    public void FindSites_DegenerateIupacEnzyme_MatchesViaAmbiguityCodes()
    {
        const string seq = "AAAGTCGACAAA"; // GTCGAC at index 3 satisfies HincII GTYRAC

        var forward = RestrictionAnalyzer.FindSites(seq, "HincII")
            .Where(s => s.IsForwardStrand)
            .ToList();

        var site = forward.Should().ContainSingle(s => s.Position == 3).Subject;
        site.RecognizedSequence.Should().Be("GTCGAC",
            "GTCGAC satisfies the degenerate GTYRAC motif (Y=C, R=G) under IUPAC matching");
        site.Enzyme.RecognitionSequence.Should().Be("GTYRAC");
    }

    #endregion

    #region RESTR-DIGEST-001 — restriction digest

    // ─────────────────────────────────────────────────────────────────────────
    // Unit: RESTR-DIGEST-001 — restriction digest
    // Checklist: docs/checklists/03_FUZZING.md, row 27.
    // Fuzz strategies exercised for THIS unit:
    //   • BE = Boundary Exploitation — a sequence with NO cut sites (the whole
    //          molecule survives as one fragment), a CIRCULAR molecule with ZERO
    //          cut sites (the wrap-around branch must not crash on the uncut
    //          circle), a 100+-ENZYME digest (a quadratic-hang / mass-loss hazard),
    //          and the EMPTY sequence (lower size boundary).
    // — docs/checklists/03_FUZZING.md §Description (strategy codes); row 27 targets:
    //   "No cut sites, circular with 0 sites, 100+ enzymes, empty seq".
    //
    // ─────────────────────────────────────────────────────────────────────────
    // The restriction-digest contract under test (Restriction_Digest_Simulation.md)
    // ─────────────────────────────────────────────────────────────────────────
    // A restriction digest cleaves a DNA molecule at enzyme cut positions and emits
    // the fragments bounded by those cuts (Restriction_Digest_Simulation.md §2.1).
    // The source surface is RestrictionAnalyzer.Digest (RestrictionAnalyzer.cs lines
    // 259–363) in two overloads, plus the private circular path DigestCircular
    // (lines 365–437):
    //   • Digest(DnaSequence, params string[] enzymeNames) — the LINEAR digest. It
    //     null-guards the sequence (ArgumentNullException, line 261) and requires at
    //     least one enzyme (ArgumentException "At least one enzyme is required",
    //     line 262–263). It collects DISTINCT FORWARD-strand cut positions only
    //     (line 273, to avoid double-counting palindromic sites), then partitions the
    //     molecule at boundaries [0, c1, …, ck, L] into the half-open intervals
    //     [0,c1),[c1,c2),…,[ck,L) — k cuts ⇒ k+1 fragments (§2.2, INV-01). Because
    //     it is a `yield` iterator, every probe forces enumeration (`.ToList()`).
    //   • Digest(DnaSequence, MoleculeTopology, params string[]) — dispatches to the
    //     linear path for MoleculeTopology.Linear and to DigestCircular for
    //     MoleculeTopology.Circular. Same null/enzyme guards (lines 355–357).
    //   • DigestCircular models a closed plasmid: k distinct forward-strand cut
    //     positions yield exactly k fragments (NOT k+1) because the molecule has no
    //     free ends — the last-cut→origin→first-cut piece joins into ONE
    //     origin-spanning fragment (§2.4 remark; "cut a circle once → one linear
    //     fragment"). A circle with ZERO cuts yields ONE full-length uncut circular
    //     fragment (line 387–397) — the wrap-around boundary that a naive
    //     implementation crashes on.
    //
    // THE KEY INVARIANT pinned on EVERY digest below (Restriction_Digest_Simulation.md
    // §2.4): INV-02 MASS CONSERVATION — the sum of fragment lengths equals the
    // original sequence length, because adjacent cut boundaries partition the
    // sequence with no overlap and no gap. This is the load-bearing fuzz assertion:
    // boundary input must never lose, duplicate, or invent sequence mass. INV-01
    // (k cuts ⇒ k+1 linear / k circular fragments) is pinned alongside.
    //
    // THE FOUR ROW-27 FUZZ TARGETS, mapped to the theory-correct contract:
    //   • No cut sites (BE): a sequence the enzyme does not cut (e.g. an all-A
    //     sequence digested with EcoRI=GAATTC) has zero forward-strand cuts, so the
    //     explicit special case (line 281–292) yields a SINGLE fragment equal to the
    //     WHOLE sequence (Length == sequence length, LeftEnzyme == RightEnzyme ==
    //     null) — never an empty result, never a crash (§6.1 "No cut sites found").
    //   • Circular with 0 sites (BE, KEY wrap-around bug): a CIRCULAR molecule the
    //     enzyme does not cut takes the DigestCircular zero-cut branch (line 387) and
    //     yields ONE full-length circular fragment — the wrap-around join logic
    //     (`seq.Substring(start, n - start) + seq.Substring(0, nextCut)`) is NEVER
    //     reached with a degenerate index, so no IndexOutOfRange / no crash on the
    //     uncut circle. Fragment length == sequence length (mass conserved).
    //   • 100+ enzymes (BE, hang hazard): digesting with 100+ enzyme names must
    //     COMPLETE promptly (no quadratic blow-up, pinned with [CancelAfter]) and the
    //     fragments must STILL conserve total length (INV-02) regardless of how many
    //     enzymes cut. Unknown names in the list contribute no sites and must not
    //     crash; known cutters contribute their sites and the partition stays exact.
    //   • Empty seq (BE): an empty DnaSequence has length 0 and no possible cut, so
    //     the linear path takes the no-cut branch and yields ONE fragment of length 0
    //     (the whole — empty — sequence); the circular path likewise yields one
    //     length-0 fragment. Either way: a single empty fragment, never a crash and
    //     never a Substring past the end. (A DnaSequence("") is the empty molecule;
    //     DnaSequence rejects only non-A/C/G/T, not emptiness.)
    //
    // Documented invariants pinned (Restriction_Digest_Simulation.md §2.4): INV-01
    // k cuts ⇒ k+1 linear / k circular fragments; INV-02 Σ fragment lengths ==
    // sequence length. Digest is a yield iterator, so every probe forces enumeration;
    // the positive-sanity test pins that a sequence with a KNOWN number of EcoRI sites
    // yields the right fragment count AND the fragment lengths sum to the sequence
    // length.

    /// <summary>
    /// BE — No cut sites: a sequence the enzyme cannot cut must survive as a SINGLE
    /// whole fragment. An all-A sequence has no EcoRI (GAATTC) site, so the digest
    /// takes the no-cut special case and yields exactly one fragment whose sequence
    /// and length equal the WHOLE input (LeftEnzyme/RightEnzyme null) — never an empty
    /// result, never a crash. Mass conservation (INV-02) is trivially the whole
    /// sequence. (Restriction_Digest_Simulation.md §6.1 "No cut sites found".)
    /// </summary>
    [Test]
    public void Digest_LinearNoCutSites_YieldsSingleWholeFragment()
    {
        const string seq = "AAAAAAAAAAAAAAAAAAAA"; // 20 A's — no GAATTC anywhere
        var dna = new DnaSequence(seq);

        var fragments = () => RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();

        var list = fragments.Should().NotThrow(
            "a sequence with no recognition site must digest without crashing").Subject;

        list.Should().ContainSingle("no cut sites ⇒ the whole molecule is one fragment (INV-01, k=0 ⇒ 1)");
        var only = list[0];
        only.Sequence.Should().Be(seq, "the single fragment IS the whole undigested sequence");
        only.Length.Should().Be(seq.Length, "INV-02: the lone fragment carries the full sequence mass");
        only.LeftEnzyme.Should().BeNull("no enzyme cuts before position 0");
        only.RightEnzyme.Should().BeNull("no enzyme cuts after the final boundary");
        list.Sum(f => f.Length).Should().Be(seq.Length, "INV-02: Σ fragment lengths == sequence length");
    }

    /// <summary>
    /// BE (KEY wrap-around) — Circular with 0 sites: a CIRCULAR molecule the enzyme
    /// cannot cut must take the uncut-circle branch and yield ONE full-length circular
    /// fragment, WITHOUT touching the origin-spanning Substring join (the common
    /// boundary bug). An all-A plasmid digested with EcoRI has zero cuts ⇒ exactly one
    /// fragment whose length == sequence length. No IndexOutOfRange on the wrap-around.
    /// (Restriction_Digest_Simulation.md §2.4 "zero cut sites yields a single
    /// full-length uncut circular fragment".)
    /// </summary>
    [Test]
    public void Digest_CircularZeroCutSites_YieldsSingleCircularFragment_NoWrapAroundCrash()
    {
        const string seq = "AAAAAAAAAAAAAAAAAAAACCCCCCCCCC"; // 30 nt, no GAATTC
        var dna = new DnaSequence(seq);

        var fragments = () => RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, "EcoRI").ToList();

        var list = fragments.Should().NotThrow(
            "the uncut-circle branch must not crash on the wrap-around logic").Subject;

        list.Should().ContainSingle("a circle with 0 cuts ⇒ exactly ONE full-length circular fragment");
        list[0].Length.Should().Be(seq.Length, "INV-02: the uncut circular fragment carries the whole sequence mass");
        list[0].Sequence.Should().Be(seq, "the lone circular fragment IS the whole undigested sequence");
        list.Sum(f => f.Length).Should().Be(seq.Length, "INV-02: Σ fragment lengths == sequence length");
    }

    /// <summary>
    /// BE (KEY circular fragment count) — Circular with k sites: a circle cut at k
    /// distinct forward-strand sites yields exactly k fragments (NOT k+1), because the
    /// ends join. "GAATTC...GAATTC..." with two EcoRI sites ⇒ two circular fragments,
    /// and Σ lengths still equals the sequence length (the origin-spanning fragment
    /// joins last-cut→end with start→first-cut). This pins that the wrap-around branch
    /// is not only crash-free but length-exact. (§2.4 "k cut sites → k fragments".)
    /// </summary>
    [Test]
    public void Digest_CircularWithTwoSites_YieldsTwoFragments_MassConserved()
    {
        // Two EcoRI (GAATTC) sites, well separated, on a clean A/C/G/T sequence.
        const string seq = "AAAAGAATTCAAAAAAAAGAATTCAAAA"; // sites at index 4 and 18
        var dna = new DnaSequence(seq);

        var linearCuts = RestrictionAnalyzer.FindSites(dna, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .Select(s => s.CutPosition)
            .Distinct()
            .Count();
        linearCuts.Should().Be(2, "the fixture is constructed to hold exactly two distinct EcoRI cut positions");

        var fragments = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, "EcoRI").ToList();

        fragments.Should().HaveCount(2, "INV-01 (circular): k=2 distinct cut sites ⇒ k=2 fragments, not k+1");
        fragments.Sum(f => f.Length).Should().Be(seq.Length,
            "INV-02: Σ circular fragment lengths == sequence length (the origin-spanning fragment joins the ends)");
        fragments.Should().OnlyContain(f => f.Length > 0, "no zero-length fragment is emitted");
        fragments.Select(f => f.Sequence.Length).Sum().Should().Be(seq.Length,
            "the emitted fragment SEQUENCES also reconstruct the full sequence mass");
    }

    /// <summary>
    /// BE (hang hazard) — 100+ enzymes: digesting with 100+ enzyme names must COMPLETE
    /// promptly (no quadratic blow-up) and the fragments must STILL conserve total
    /// length (INV-02), no matter how many enzymes cut. The list mixes the full
    /// built-in catalog (repeated to exceed 100 names) so many real cutters fire; the
    /// partition must stay exact. Pinned with [CancelAfter] so a hang fails as a
    /// timeout rather than wedging the suite.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void Digest_OneHundredPlusEnzymes_CompletesAndConservesMass()
    {
        var rng = new Random(27_001);
        var bases = "ACGT";
        var sb = new System.Text.StringBuilder(500);
        for (int i = 0; i < 500; i++)
            sb.Append(bases[rng.Next(bases.Length)]);
        var dna = new DnaSequence(sb.ToString());

        // Build a 100+-name enzyme list from the built-in catalog (repeat to exceed 100).
        var catalog = RestrictionAnalyzer.Enzymes.Values.Select(e => e.Name).ToList();
        catalog.Should().NotBeEmpty("the built-in catalog must supply enzyme names for this probe");
        var enzymeNames = new List<string>();
        while (enzymeNames.Count <= 110)
            enzymeNames.AddRange(catalog);
        enzymeNames.Count.Should().BeGreaterThan(100, "the digest is fuzzed with 100+ enzyme names");

        var enzymeArray = enzymeNames.ToArray();
        var fragments = () => RestrictionAnalyzer.Digest(dna, enzymeArray).ToList();

        var list = fragments.Should().NotThrow(
            "100+ enzymes must complete without crashing or hanging").Subject;

        list.Should().NotBeEmpty("a non-empty sequence always yields at least one fragment");
        list.Sum(f => f.Length).Should().Be(dna.Length,
            "INV-02: Σ fragment lengths == sequence length no matter how many enzymes cut");
        list.Should().OnlyContain(f => f.Length > 0, "no zero-length fragment is ever emitted");

        // The same 100+-name digest on the CIRCULAR topology must also conserve mass.
        var circular = RestrictionAnalyzer
            .Digest(dna, MoleculeTopology.Circular, enzymeNames.ToArray()).ToList();
        circular.Sum(f => f.Length).Should().Be(dna.Length,
            "INV-02 (circular): Σ fragment lengths == sequence length under a 100+-enzyme digest");
    }

    /// <summary>
    /// BE — Empty seq: an empty molecule has no possible cut, so both the linear and
    /// circular digests take their no-cut branch and yield a SINGLE length-0 fragment
    /// (the whole — empty — sequence). No crash, no Substring past the end, and mass
    /// conservation holds trivially (Σ lengths == 0). DnaSequence("") is the empty
    /// molecule (DnaSequence rejects only non-A/C/G/T, not emptiness).
    /// </summary>
    [Test]
    public void Digest_EmptySequence_YieldsSingleEmptyFragment()
    {
        var dna = new DnaSequence(string.Empty);
        dna.Length.Should().Be(0, "the empty molecule has length 0");

        var linear = () => RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();
        var circular = () => RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, "EcoRI").ToList();

        var lin = linear.Should().NotThrow("the empty sequence must digest without crashing").Subject;
        lin.Should().ContainSingle("no cut sites on an empty molecule ⇒ one fragment");
        lin[0].Length.Should().Be(0, "the single fragment of an empty molecule has length 0");
        lin[0].Sequence.Should().BeEmpty("the single fragment IS the empty sequence");
        lin.Sum(f => f.Length).Should().Be(0, "INV-02: Σ fragment lengths == sequence length (0)");

        var circ = circular.Should().NotThrow(
            "the empty circular molecule must not crash on the wrap-around branch").Subject;
        circ.Should().ContainSingle("a circle with 0 cuts ⇒ one fragment");
        circ[0].Length.Should().Be(0, "the empty circular fragment has length 0");
        circ.Sum(f => f.Length).Should().Be(0, "INV-02: Σ fragment lengths == sequence length (0)");
    }

    /// <summary>
    /// Positive sanity: a sequence with a KNOWN number of EcoRI sites must yield the
    /// right fragment COUNT and the fragment lengths must SUM to the sequence length,
    /// so the boundary hardening never silently breaks the core partition. Three
    /// well-separated GAATTC sites ⇒ 3 distinct forward-strand cuts ⇒ exactly 4 linear
    /// fragments (INV-01, k+1), and Σ lengths == sequence length (INV-02). The same
    /// molecule digested as a CIRCLE yields exactly 3 fragments (INV-01 circular, k).
    /// </summary>
    [Test]
    public void Digest_KnownEcoRISites_YieldCorrectFragmentCountAndConserveMass()
    {
        // Three EcoRI (GAATTC) sites, well separated.
        const string seq = "AAAAGAATTCAAAAAAAAGAATTCAAAAAAAAGAATTCAAAA";
        var dna = new DnaSequence(seq);

        var cuts = RestrictionAnalyzer.FindSites(dna, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .Select(s => s.CutPosition)
            .Distinct()
            .Count();
        cuts.Should().Be(3, "the fixture holds exactly three distinct EcoRI cut positions");

        var linear = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();
        linear.Should().HaveCount(4, "INV-01 (linear): k=3 cuts ⇒ k+1=4 fragments");
        linear.Sum(f => f.Length).Should().Be(seq.Length, "INV-02: Σ linear fragment lengths == sequence length");
        linear.First().LeftEnzyme.Should().BeNull("the first fragment has no enzyme to its left");
        linear.Last().RightEnzyme.Should().BeNull("the last fragment has no enzyme to its right");

        var circular = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, "EcoRI").ToList();
        circular.Should().HaveCount(3, "INV-01 (circular): k=3 cuts ⇒ k=3 fragments (the ends join)");
        circular.Sum(f => f.Length).Should().Be(seq.Length, "INV-02: Σ circular fragment lengths == sequence length");
    }

    #endregion

    #endregion
}
