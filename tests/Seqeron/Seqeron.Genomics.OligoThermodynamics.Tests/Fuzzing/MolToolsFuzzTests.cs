// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
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
/// ThermodynamicPrimerDesigner.DesignPrimers(DnaSequence template, int targetStart, int targetEnd,
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
/// ThermodynamicProbeDesigner.DesignProbes(string targetSequence, ProbeParameters? parameters = null,
///   int maxProbes = 10) (ProbeDesigner.cs lines 127–146) is a YIELD iterator that scans the
/// target for every length-`[MinLength..MaxLength]` window, scores each candidate by a fixed
/// additive penalty model (GC, Tm, homopolymer, self-complementarity, secondary structure,
/// repeats, terminal G/C — EvaluateProbeWithGc, lines 253–327), keeps only RAW-SCORE-POSITIVE
/// candidates (`score <= 0 → null`, line 315; INV-01), ranks them by score, and yields the top
/// `maxProbes`. Each Probe carries its Tm, GcContent, Score, and a Warnings list.
///
/// CRITICAL DESIGN FACT the checklist row probes — UNLIKE ThermodynamicPrimerDesigner.DesignPrimers,
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
/// ThermodynamicProbeDesigner.ValidateProbe(string probeSequence, IEnumerable&lt;string&gt;
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



    // ═══════════════════════════════════════════════════════════════════
    //  PRIMER-DESIGN-001 — primer design : fuzz targets
    // ═══════════════════════════════════════════════════════════════════



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

        var act = () => ThermodynamicPrimerDesigner.DesignPrimers(template, targetStart: 1, targetEnd: 2);

        act.Should().NotThrow(
            "no MinLength-bp primer can fit in a 4-nt template; the candidate loops never run instead of slicing past the end");

        var result = ThermodynamicPrimerDesigner.DesignPrimers(template, 1, 2);
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

        var negativeStart = () => ThermodynamicPrimerDesigner.DesignPrimers(template, targetStart: -1, targetEnd: 50);
        var endPastEnd = () => ThermodynamicPrimerDesigner.DesignPrimers(template, targetStart: 10, targetEnd: 100);
        var startGeEnd = () => ThermodynamicPrimerDesigner.DesignPrimers(template, targetStart: 50, targetEnd: 50);

        negativeStart.Should().Throw<ArgumentException>("a negative target start is rejected by the region guard");
        endPastEnd.Should().Throw<ArgumentException>("a target end at/past the template length is rejected");
        startGeEnd.Should().Throw<ArgumentException>("a target start not strictly before the end is rejected");
    }



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
        var param = ThermodynamicPrimerDesigner.DefaultParameters with
        {
            MinGcContent = 0,
            MaxGcContent = 0,
            MinTm = 0,
            MaxTm = 200,
            MaxHomopolymer = 1000,   // an all-A primer is one long homopolymer; do not reject on it
            MaxDinucleotideRepeats = 1000,
            Check3PrimeStability = false
        };

        var act = () => ThermodynamicPrimerDesigner.DesignPrimers(template, targetStart: 120, targetEnd: 130, param);

        act.Should().NotThrow(
            "a zero-width GC range is a pair of comparisons, not a division; an all-A·T primer satisfies gc<0||gc>0 being false");

        var result = ThermodynamicPrimerDesigner.DesignPrimers(template, 120, 130, param);
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
        var param = ThermodynamicPrimerDesigner.DefaultParameters with
        {
            MinGcContent = 0,
            MaxGcContent = 0,
            MinTm = 0,
            MaxTm = 200,
            Check3PrimeStability = false
        };

        var act = () => ThermodynamicPrimerDesigner.DesignPrimers(template, targetStart: 120, targetEnd: 130, param);

        act.Should().NotThrow("a [0,0] GC range over a GC-bearing template empties the feasible set rather than crashing");

        ThermodynamicPrimerDesigner.DesignPrimers(template, 120, 130, param)
            .IsValid.Should().BeFalse("every candidate carries a G/C, so its GC% > 0 fails the [0,0] range — no valid pair");
    }



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
        var inverted = ThermodynamicPrimerDesigner.DefaultParameters with
        {
            MinTm = 90,   // MinTm > MaxTm: a contradiction
            MaxTm = 10,
            MinGcContent = 0,
            MaxGcContent = 100,
            MaxHomopolymer = 1000,
            MaxDinucleotideRepeats = 1000,
            Check3PrimeStability = false
        };

        var act = () => ThermodynamicPrimerDesigner.DesignPrimers(template, targetStart: 120, targetEnd: 130, inverted);

        act.Should().NotThrow(
            "an inverted Tm range is not validated as an error; it simply makes every candidate fail the Tm screen");

        var result = ThermodynamicPrimerDesigner.DesignPrimers(template, 120, 130, inverted);
        result.IsValid.Should().BeFalse("no finite Tm can satisfy tm >= 90 AND tm <= 10, so no candidate is valid (INV-01)");
        result.Forward.Should().BeNull("the contradictory Tm range leaves no valid forward primer");
        result.Reverse.Should().BeNull("the contradictory Tm range leaves no valid reverse primer");
    }



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
        var param = ThermodynamicPrimerDesigner.DefaultParameters;
        // A 600-nt balanced random template with a central target region leaves ~200 bp on each
        // flank for candidate enumeration — deterministic via a fixed seed.
        var template = new DnaSequence(RandomDna(600, seed: 22_000));

        var result = ThermodynamicPrimerDesigner.DesignPrimers(template, targetStart: 300, targetEnd: 320, param);

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



    // ═══════════════════════════════════════════════════════════════════
    //  PROBE-DESIGN-001 — hybridization probe design : fuzz targets
    // ═══════════════════════════════════════════════════════════════════



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

        var act = () => ThermodynamicProbeDesigner.DesignProbes(shortTarget).ToList();

        act.Should().NotThrow(
            "a target shorter than MinLength trips the length short-circuit before any windowing; no Substring is taken past the end");

        ThermodynamicProbeDesigner.DesignProbes(shortTarget).Should().BeEmpty(
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
        var param = new ThermodynamicProbeDesigner.ProbeParameters(
            MinLength: 20, MaxLength: 20,
            MinTm: 0, MaxTm: 200,           // permissive Tm so the length boundary is isolated
            MinGc: 0.0, MaxGc: 1.0,         // permissive GC
            MaxHomopolymer: 20,
            AvoidSecondaryStructure: false,
            MaxSelfComplementarity: 1.0);

        // A balanced 20-mer (10 GC) so it is a clean candidate; the 19-mer is one short.
        const string at20 = "ACGTACGTACGTACGTACGT"; // length 20
        string below19 = at20.Substring(0, 19);     // length 19

        ThermodynamicProbeDesigner.DesignProbes(below19, param).Should().BeEmpty(
            "a target one base shorter than MinLength has no candidate window");
        ThermodynamicProbeDesigner.DesignProbes(at20, param).Should().NotBeEmpty(
            "a target of exactly MinLength holds one window, which under permissive screens is a valid probe");
    }



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
        var param = ThermodynamicProbeDesigner.Defaults.Microarray with { MinTm = 95.0, MaxTm = 5.0 };
        string target = RandomDna(120, seed: 24_002);

        var act = () => ThermodynamicProbeDesigner.DesignProbes(target, param).ToList();

        act.Should().NotThrow(
            "an inverted Tm range is a contradictory constraint; the source does not validate it, so it neither throws nor hangs");

        var probes = ThermodynamicProbeDesigner.DesignProbes(target, param).ToList();
        probes.Should().OnlyContain(
            p => p.Warnings.Any(w => w.Contains("Tm") && w.Contains("outside range")),
            "an inverted Tm range fails the Tm screen for EVERY candidate, so every returned probe records the Tm-out-of-range warning rather than being silently accepted");
        probes.Should().OnlyContain(p => p.Score > 0,
            "INV-01: only raw-score-positive candidates are retained, even when the Tm screen always fails");
    }



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
        var empty = () => ThermodynamicProbeDesigner.DesignProbes(string.Empty).ToList();
        var nullTarget = () => ThermodynamicProbeDesigner.DesignProbes((string)null!).ToList();

        empty.Should().NotThrow("an empty target short-circuits to an empty result before any windowing");
        nullTarget.Should().NotThrow(
            "IsNullOrEmpty catches the null target, so it yields the empty result rather than dereferencing into a NullReferenceException");

        ThermodynamicProbeDesigner.DesignProbes(string.Empty).Should().BeEmpty();
        ThermodynamicProbeDesigner.DesignProbes((string)null!).Should().BeEmpty();
    }



    /// <summary>
    /// Positive sanity: alongside the degenerate probes, a reasonable target must yield at
    /// least one VALID probe whose GC and Tm fall INSIDE the requested ranges — so the
    /// boundary hardening never silently breaks the core function. The parameter ranges are
    /// chosen MUTUALLY SATISFIABLE for the salt-adjusted Tm model: for a 50–60-nt probe the
    /// source computes Tm = 81.5 + 16.6·log10(0.05) + 41·GC − 600/length, so a 0.40–0.60-GC
    /// window maps to roughly 64–74 °C — hence the GC window [0.40, 0.60] is paired with a Tm
    /// window [60, 80] (the former Microarray Tm default of 75–85 was unreachable at ≤ 0.60 GC
    /// for this formula; since audit round 3 A3-10 the Microarray preset uses OligoArray's
    /// 1 M / 1 µM nearest-neighbour scale with a reachable 82–90 window — this test keeps the
    /// Primer3 50 mM conditions of a plain ProbeParameters). The target is a 300-nt fixed-seed sequence;
    /// AT LEAST ONE returned probe must satisfy INV-01 (Score &gt; 0), INV-02
    /// (0 ≤ GcContent ≤ 1), have GC within [MinGc, MaxGc] AND Tm within [MinTm, MaxTm], and —
    /// being fully in-range — carry NO GC or Tm out-of-range warning.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void DesignProbes_ReasonableTarget_YieldsInRangeValidProbe()
    {
        // GC and Tm windows are mutually satisfiable for the 50–60-nt salt-adjusted Tm model.
        var param = new ThermodynamicProbeDesigner.ProbeParameters(
            MinLength: 50, MaxLength: 60,
            MinTm: 60, MaxTm: 80,
            MinGc: 0.40, MaxGc: 0.60,
            MaxHomopolymer: 5,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.3);
        string target = RandomDna(300, seed: 24_003);

        var probes = ThermodynamicProbeDesigner.DesignProbes(target, param, maxProbes: 20).ToList();

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



    // ═══════════════════════════════════════════════════════════════════
    //  PROBE-VALID-001 — probe validation : fuzz targets
    // ═══════════════════════════════════════════════════════════════════



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

        var act = () => ThermodynamicProbeDesigner.ValidateProbe(junk, references);
        act.Should().NotThrow(
            "a non-DNA probe is scanned literally, not validated; junk mismatches references and never indexes out of range");

        var result = ThermodynamicProbeDesigner.ValidateProbe(junk, references);

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

        var act = () => ThermodynamicProbeDesigner.ValidateProbe(probe, references);
        act.Should().NotThrow(
            "injected null bytes / unicode are compared literally and never index out of range");

        var result = ThermodynamicProbeDesigner.ValidateProbe(probe, references);
        result.SpecificityScore.Should().BeInRange(0.0, 1.0);
        result.SelfComplementarity.Should().BeInRange(0.0, 1.0);
        result.OffTargetHits.Should().BeGreaterThanOrEqualTo(0,
            "INV-03: injected junk can only suppress a match, never invent an off-target");
    }



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
            var act = () => ThermodynamicProbeDesigner.ValidateProbe(probe, references);
            act.Should().NotThrow($"a {probe.Length}-nt probe has no minimum-length screen and is scanned, not rejected");

            var result = ThermodynamicProbeDesigner.ValidateProbe(probe, references);
            result.SpecificityScore.Should().BeInRange(0.0, 1.0, "INV-01 holds for any non-empty probe length");
            result.SelfComplementarity.Should().BeInRange(0.0, 1.0, "INV-02: a 1–3 nt probe never divides by zero");
            double.IsNaN(result.SelfComplementarity).Should().BeFalse("non-zero length ⇒ no NaN");
            result.HasSecondaryStructure.Should().BeFalse(
                "a probe shorter than 2·stem+loop (11 nt) cannot form a hairpin; the screen short-circuits");
        }
    }



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

        var nullProbe = () => ThermodynamicProbeDesigner.ValidateProbe(null!, references);
        var nullRefs = () => ThermodynamicProbeDesigner.ValidateProbe("ACGTACGT", null!);

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

        var act = () => ThermodynamicProbeDesigner.ValidateProbe(string.Empty, references);
        act.Should().NotThrow("the empty probe is a defined degenerate case, not an exception");

        var result = ThermodynamicProbeDesigner.ValidateProbe(string.Empty, references);
        result.IsValid.Should().BeFalse("an empty probe cannot hybridize specifically");
        result.SpecificityScore.Should().Be(0.0, "the empty-probe structured result fixes specificity at 0.0");
        result.OffTargetHits.Should().Be(0);
        result.Issues.Should().Contain("Empty probe sequence");
    }



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

        var act = () => ThermodynamicProbeDesigner.ValidateProbe(probe, references);
        act.Should().NotThrow("'N' is not a wildcard here; it is compared literally and complemented to 'N' — no crash");

        var result = ThermodynamicProbeDesigner.ValidateProbe(probe, references);

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

        var result = ThermodynamicProbeDesigner.ValidateProbe(probe, references);

        double.IsNaN(result.SpecificityScore).Should().BeFalse();
        double.IsNaN(result.SelfComplementarity).Should().BeFalse();
        result.OffTargetHits.Should().BeGreaterThanOrEqualTo(1,
            "the probe occurs verbatim (its N's match themselves), so at least one approximate hit is counted");
        result.SpecificityScore.Should().BeInRange(0.0, 1.0, "INV-01");
    }



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

        var result = ThermodynamicProbeDesigner.ValidateProbe(probe, references, maxMismatches: 0);

        result.OffTargetHits.Should().Be(1, "the probe occurs exactly once at mismatch tolerance 0");
        result.SpecificityScore.Should().Be(1.0, "INV-04: a single hit maps to specificity 1.0");
        result.IsValid.Should().BeTrue("a unique, low-self-complementarity probe is a valid probe");
        result.Issues.Should().BeEmpty("a clean unique probe records no validation issues");
    }


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
}
