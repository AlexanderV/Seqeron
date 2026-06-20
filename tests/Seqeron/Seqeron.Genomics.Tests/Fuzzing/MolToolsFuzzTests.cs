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
}
