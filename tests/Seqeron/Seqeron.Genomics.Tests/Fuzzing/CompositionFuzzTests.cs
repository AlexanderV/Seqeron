using System;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Composition area — GC content.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, OverflowException, …). Every input must result in
/// EITHER a well-defined, theory-correct value, OR a *documented, intentional*
/// validation exception (ArgumentException / ArgumentNullException). A raw
/// runtime exception or a hang on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-GC-001 — GC content (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 1.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — empty string, single char, extremely long.
///   • INJ = Injection — non-ACGT characters, null, unicode (combining marks,
///           astral/surrogate-pair code points, null byte '\0').
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The GC% contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// GC% = (G + C) / (A + T + G + C) × 100, case-insensitive.
///   — docs/algorithms/Statistics/GC_Content_Profile.md §2.2 [Wikipedia GC-content];
///     docs/algorithms/Sequence_Composition/Sequence_Composition.md §2.2.
///
/// API entry: DnaSequence.GcContent()
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs), backed by
///   SequenceExtensions.CalculateGcContent
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs).
///
/// Documented input handling along the public DnaSequence path
/// (DnaSequence.cs lines 22–33, 112–124):
///   • null or empty string  → an empty sequence; GcContent() == 0
///     (string.IsNullOrEmpty short-circuit). This is a *defined result*, NOT an
///     exception: the public surface treats null as "no sequence".
///   • input is case-folded with ToUpperInvariant, then validated; so lowercase
///     a/c/g/t round-trips to the same GC% as uppercase.
///   • ANY character that is not A/C/G/T after upper-casing (digits, whitespace,
///     N/IUPAC ambiguity codes, U, '\0', unicode letters, combining marks,
///     astral code points) → a *documented, intentional* ArgumentException from
///     ValidateSequence. This is the contract's validation gate, not a crash.
///   • a valid, extremely long sequence → computed without overflow/hang;
///     result stays in [0, 100].
///
/// The backing span primitive SequenceExtensions.CalculateGcContent is *lenient*
/// by separate documented design: it excludes non-A/T/G/C/U symbols from BOTH
/// numerator and denominator and returns 0 when no valid base is present
/// (SequenceExtensions.cs lines 17–58, matching Biopython gc_fraction "remove").
/// We pin that contract too, so the boundary between the strict public API and
/// the lenient primitive is explicit and cannot silently drift.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-COMP-001 — DNA complement (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 2.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — empty string, single char, null.
///   • INJ = Injection — non-DNA characters, mixed case, unicode (accented
///           Latin, Greek, combining marks, full-width look-alikes, and
///           astral/surrogate-pair code points).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// The DNA-complement contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// DNA complement maps A↔T and C↔G (Watson–Crick base pairing).
///   — docs/algorithms/Sequence_Composition/RNA_Complement.md §2.1–2.2 (the DNA
///     sibling GetComplementBase, explicitly cited there as SEQ-COMP-001);
///     docs/mcp/tools/sequence/complement_base.md.
///
/// SEQ-COMP-001 has TWO documented surfaces with DIFFERENT, intentional contracts.
/// Fuzzing must pin both, and pin the boundary between them so neither can drift:
///
/// (1) The STRICT public path — DnaSequence.Complement()
///     (src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs lines 54–62),
///     which validates its input at construction (DnaSequence ctor +
///     ValidateSequence, lines 22–33 / 112–124):
///       • null or empty string  → an empty sequence; Complement() is the empty
///         sequence; no exception (string.IsNullOrEmpty short-circuit). A defined
///         result, NOT an error.
///       • input is case-folded with ToUpperInvariant before validation, so
///         lowercase / mixed case a-c-g-t is accepted and complements identically
///         to uppercase.
///       • ANY character that is not A/C/G/T after upper-casing (digits,
///         whitespace, N/IUPAC ambiguity codes, U, unicode letters, combining
///         marks, astral code points) → a *documented, intentional*
///         ArgumentException from ValidateSequence. The validation gate, not a
///         crash and not a silent mis-complement.
///
/// (2) The LENIENT char/span primitive — SequenceExtensions.GetComplementBase
///     (SequenceExtensions.cs lines 137–157) and the span overload
///     ReadOnlySpan&lt;char&gt;.TryGetComplement (lines 204–215). By separate
///     documented design (RNA_Complement.md §3.3, the DNA sibling) this surface:
///       • is IUPAC-complete: A↔T, C↔G, U→A, and the eleven ambiguity codes
///         (R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N);
///       • is case-insensitive and always emits UPPERCASE for recognized symbols;
///       • passes ANY non-IUPAC character (gaps, digits, whitespace, '\0',
///         unicode letters, surrogate halves) THROUGH UNCHANGED and NEVER throws;
///       • on an empty span succeeds and writes nothing — no exception, no hang.
///     We pin this so the lenient/strict boundary is explicit: the same garbage
///     that the public DnaSequence path REJECTS, the primitive must carry through
///     without crashing.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-REVCOMP-001 — reverse complement (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 3.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — empty string, single char, null.
///   • INJ = Injection — non-DNA characters, unicode (accented Latin, Greek,
///           combining marks, full-width look-alikes, and astral/surrogate-pair
///           code points), embedded null byte.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// The reverse-complement contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Reverse complement = reverse ∘ complement: complement each base (A↔T, C↔G,
/// Watson–Crick pairing) and read the strand 5'→3', i.e. in reverse order. It is
/// an INVOLUTION: revcomp(revcomp(x)) == x for any sequence.
///   — docs/algorithms/Sequence_Composition/RNA_Complement.md §"DNA span helpers"
///     (the DNA reverse-complement is composed from GetComplementBase, cited there
///     as the DNA sibling); Biopython Bio.Seq.reverse_complement worked examples
///     (RNA_Complement.md ref 4).
///
/// SEQ-REVCOMP-001 has THREE documented surfaces with DIFFERENT, intentional
/// contracts. Fuzzing pins all three, and the boundary between them, so none can
/// silently drift:
///
/// (1) The STRICT public path — DnaSequence.ReverseComplement()
///     (DnaSequence.cs lines 68–76), which validates its input at construction
///     (DnaSequence ctor + ValidateSequence, lines 22–33 / 112–124):
///       • null or empty string  → an empty sequence; ReverseComplement() is the
///         empty sequence; no exception (string.IsNullOrEmpty short-circuit). A
///         defined result, NOT an error.
///       • input is case-folded with ToUpperInvariant before validation, so a
///         single lowercase base and mixed case a-c-g-t are accepted and
///         reverse-complement identically to uppercase.
///       • a single base maps to its complement (A→T, C→G, G→C, T→A): with one
///         base, reverse is a no-op, so revcomp is just the complement.
///       • ANY character that is not A/C/G/T after upper-casing (digits,
///         whitespace, N/IUPAC ambiguity codes, U, unicode letters, combining
///         marks, astral code points, '\0') → a *documented, intentional*
///         ArgumentException from ValidateSequence. The validation gate, not a
///         crash and not a silent mis-complement.
///       • the result is re-wrapped as a DnaSequence; since the complement of
///         valid A/C/G/T is again valid A/C/G/T, that re-validation never throws.
///         This pins that the involution holds: revcomp(revcomp(x)) == x.
///
/// (2) The LENIENT static string helper — DnaSequence.GetReverseComplementString
///     (DnaSequence.cs lines 149–160). By design it does NOT validate: it maps
///     through GetComplementBase, so it is IUPAC-complete, always emits UPPERCASE
///     for recognized symbols, passes ANY non-IUPAC character (gap, digit,
///     whitespace, '\0', unicode letter, surrogate half) THROUGH UNCHANGED, and
///     NEVER throws. On null/empty it returns the input verbatim — no exception.
///
/// (3) The LENIENT span primitive — ReadOnlySpan&lt;char&gt;.TryGetReverseComplement
///     (SequenceExtensions.cs lines 220–231). Same lenient char mapping; on an
///     empty span it succeeds and writes nothing — no exception, no hang. We pin
///     this so the strict/lenient boundary is explicit: the same garbage the
///     public path REJECTS, both lenient surfaces must carry through without
///     crashing and without shifting the reverse-complement of the valid bases.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-VALID-001 — sequence validation (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 4.
/// Fuzz strategies exercised for THIS unit:
///   • RB  = Random Bytes — fixed-seed sweeps of random chars (full BMP code-point
///           range, including control chars and surrogate halves) and random
///           ASCII strings; the validators must never throw and must classify
///           valid iff every char is in the accepted alphabet.
///   • INJ = Injection — non-ASCII letters, null bytes, mixed-case, unicode
///           (combining marks, full-width look-alikes, astral/surrogate-pair code
///           points), control characters (\0, \t, \n, \r, BEL, DEL, ESC).
///   • BE  = Boundary Exploitation — empty string, single char, extremely long.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The sequence-validation contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// SEQ-VALID-001 is the sequence-validation entry point: a per-character
/// set-membership scan against an alphabet, normalized to uppercase, returning a
/// boolean (never throwing) — Sequence_Validation.md §2.2, §4.1, §5. Strict mode:
/// DNA accepts ONLY {A,C,G,T}; RNA accepts ONLY {A,C,G,U}. IUPAC ambiguity codes
/// (N,R,Y,S,W,K,M,B,D,H,V) and the gap '-' are REJECTED (Sequence_Validation.md
/// §5.2–5.4, INV-01..03). The EMPTY sequence is valid by vacuous truth — no
/// invalid character is present (Sequence_Validation.md §3.3, §6.1).
///
/// SEQ-VALID-001 has THREE documented surfaces. Fuzzing pins all three and the
/// boundary between them so none can silently drift:
///
/// (1) SequenceExtensions.IsValidDna(ReadOnlySpan&lt;char&gt;)
///     (SequenceExtensions.cs lines 302–311): a TOTAL predicate — for ANY input
///     it returns true/false and NEVER throws. Case-insensitive
///     (char.ToUpperInvariant per char). true iff every char ∈ {A,C,G,T};
///     empty → true (vacuous truth). Because it folds char-by-char, surrogate
///     halves, null bytes, control chars and astral code points are simply "not
///     A/C/G/T" → false, never a crash and never an encoding surprise.
///
/// (2) SequenceExtensions.IsValidRna(ReadOnlySpan&lt;char&gt;)
///     (SequenceExtensions.cs lines 317–326): identical contract over {A,C,G,U}.
///     The DNA/RNA asymmetry is documented (Sequence_Validation.md §5.2 table):
///     "ACGT" is valid DNA but INVALID RNA; "ACGU" is valid RNA but INVALID DNA.
///     We pin that asymmetry so neither alphabet can drift into the other.
///
/// (3) DnaSequence.TryCreate(string, out DnaSequence?)
///     (DnaSequence.cs lines 129–141): factory validation. Returns true with a
///     materialized sequence when the DnaSequence ctor accepts the input; returns
///     false with null ONLY when the ctor raises the documented ArgumentException
///     (Sequence_Validation.md §3.3, §5.1). null/empty input → true with an empty
///     sequence (the ctor's IsNullOrEmpty short-circuit, DnaSequence.cs lines
///     24–28) — TryCreate does NOT treat "no input" as a failure. TryCreate only
///     catches ArgumentException, so any *other* exception type leaking from the
///     validation path would surface here — fuzzing pins that no such leak occurs
///     on random bytes, control chars, null bytes, unicode or huge input.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-COMPLEX-001 — sequence (linguistic) complexity (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 5.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — empty string, single char, all-same
///           nucleotide (homopolymer), extremely long.
///   • RB  = Random Bytes — fixed-seed random sequences over {A,C,G,T} (and, on
///           the lenient raw-string surface, random BMP code points / arbitrary
///           chars) asserting no unhandled throw, a finite result, and the
///           theory-correct ordering homopolymer ≤ diverse.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The sequence-complexity contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// SEQ-COMPLEX-001 is LINGUISTIC complexity (LC): the summation-form ratio of
/// observed distinct subwords to the maximum possible distinct subwords, summed
/// over word lengths 1..min(maxWordLength, N):
///     LC = Σ_i V_i / Σ_i V_max,i,   V_max,i = min(4^i, N − i + 1)
/// (DNA alphabet K = 4). LC = 1 is maximum complexity; low values flag
/// repeats / low-complexity DNA. A homopolymer observes exactly ONE distinct word
/// per length → MINIMUM LC; a maximally diverse sequence approaches MAXIMUM.
///   — docs/algorithms/Sequence_Composition/Linguistic_Complexity.md §2.2, §2.4
///     (INV-01: 0 ≤ LC ≤ 1 for DNA; INV-02: empty → 0), §6.1 (edge cases).
///
/// API entry: SequenceComplexity.CalculateLinguisticComplexity(...)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs).
///
/// SEQ-COMPLEX-001 has TWO documented surfaces with DIFFERENT validation
/// contracts. Fuzzing pins both, and the boundary between them, so neither drifts:
///
/// (1) The TYPED overload — CalculateLinguisticComplexity(DnaSequence, int)
///     (SequenceComplexity.cs lines 22–28). The DnaSequence argument has ALREADY
///     passed the strict ctor validation gate, so only A/C/G/T (upper-cased) ever
///     reach the metric. Documented validation (Linguistic_Complexity.md §3.3):
///       • null DnaSequence            → ArgumentNullException (explicit guard);
///       • maxWordLength < 1           → ArgumentOutOfRangeException;
///       • empty sequence              → 0 (INV-02 short-circuit, no division);
///       • single base                 → a positive value (one distinct 1-mer);
///       • homopolymer                 → low LC (one word per length);
///       • valid extremely long input  → finite LC in [0, 1], no overflow/hang.
///     Because non-ACGT cannot reach this overload, the result is always in the
///     DNA-bounded [0, 1] interval (INV-01).
///
/// (2) The RAW-STRING overload — CalculateLinguisticComplexity(string, int)
///     (SequenceComplexity.cs lines 33–37). LENIENT by documented design: it
///     short-circuits null/empty to 0, upper-cases, and does NOT validate the
///     alphabet (Linguistic_Complexity.md §5.2, §6.1). So it must NEVER throw on
///     arbitrary chars (digits, gaps, unicode, null byte, surrogate halves) — it
///     just counts whatever subwords appear. The denominator stays the DNA 4^i,
///     so non-ACGT input MAY exceed 1 (documented in §5.3/§6.1) — a *defined*
///     consequence, not a crash. maxWordLength < 1 is NOT validated here; with the
///     min(maxWordLength, N) cap the inner loop simply never runs and the result
///     is 0. We pin that this surface stays finite and non-throwing on pure
///     random-byte garbage, so the strict/lenient boundary is explicit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-ENTROPY-001 — Shannon entropy of base composition (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 6.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — empty string, single symbol, all-same
///           nucleotide (homopolymer), extremely long.
///   • RB  = Random Bytes — fixed-seed random sequences over {A,C,G,T} (and, on
///           the lenient raw-string surface, arbitrary BMP code points) asserting
///           no unhandled throw, NO NaN/Inf leak, a finite result inside the
///           documented bounds, and the theory extremes all-same → 0 and uniform
///           over k symbols → log2(k).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The Shannon-entropy contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// SEQ-ENTROPY-001 is base-composition Shannon entropy, H = −Σ pᵢ log₂ pᵢ, in
/// BITS (log base 2). For DNA over the 4-symbol alphabet {A,T,G,C} the value lies
/// in [0, 2]: an all-same homopolymer is the MINIMUM (H = 0); the uniform 25%/base
/// distribution is the MAXIMUM (H = log2(4) = 2). More generally, a uniform
/// distribution over k observed symbols has H = log2(k). The p·log(p) term is
/// handled as 0 when p = 0 by SKIPPING zero-count bases (count>0 guard), so no
/// −0·log2(0) = NaN ever leaks, and a 0/0 division (no counted base present) is
/// avoided by an explicit total==0 → 0 guard.
///   — docs/algorithms/Sequence_Composition/Shannon_Entropy.md §2.2 (core model),
///     §2.4 (INV-01 range [0,2]; INV-02 empty → 0; INV-03 homopolymer → 0),
///     §4.2 (reference table: uniform → 2.0, 50/50 → 1.0, single → 0.0),
///     §5.2/§5.3 (counts only A/T/G/C; ignores N/ambiguity; 0 if none present),
///     §6.1 (edge cases). Sources: Shannon (1948), Cover &amp; Thomas (1991).
///
/// API entry: SequenceComplexity.CalculateShannonEntropy(...)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs).
///
/// SEQ-ENTROPY-001 has TWO documented surfaces with DIFFERENT validation
/// contracts. Fuzzing pins both, and the boundary between them, so neither drifts:
///
/// (1) The TYPED overload — CalculateShannonEntropy(DnaSequence)
///     (SequenceComplexity.cs lines 78–82). The DnaSequence argument has ALREADY
///     passed the strict ctor validation gate, so only A/C/G/T (upper-cased) ever
///     reach the metric. Documented validation (Shannon_Entropy.md §3.3):
///       • null DnaSequence            → ArgumentNullException (explicit guard,
///                                        never a NullReferenceException);
///       • empty sequence              → 0 (INV-02 short-circuit, no division);
///       • single base / homopolymer   → exactly 0 (INV-03, the p=1 term is 0);
///       • uniform over k bases         → exactly log2(k); 25%/base → 2.0;
///       • valid extremely long input  → finite H in [0, 2], no overflow/hang.
///     Because non-ACGT cannot reach this overload, the result is always in the
///     DNA-bounded [0, 2] interval (INV-01).
///
/// (2) The RAW-STRING overload — CalculateShannonEntropy(string)
///     (SequenceComplexity.cs lines 87–91). LENIENT by documented design: it
///     short-circuits null/empty to 0, upper-cases, and counts ONLY A/T/G/C —
///     every other character (N, ambiguity codes, U, digits, gaps, whitespace,
///     null byte, unicode, surrogate halves) is IGNORED, not validated
///     (Shannon_Entropy.md §5.2, §5.3). So it must NEVER throw on arbitrary chars,
///     and because only ≤4 DNA symbols are ever counted the result stays a FINITE
///     value in [0, 2] — ignored garbage does not shift the entropy of the A/T/G/C
///     bases that ARE present. When NO A/T/G/C base is present at all, the
///     total==0 guard returns the defined 0 rather than a 0/0 NaN. We pin that
///     this surface stays finite, non-throwing, and NaN-free on pure random-byte
///     garbage, so the strict/lenient boundary is explicit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-GCSKEW-001 — GC skew (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 7.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — empty string, single base, no G or C
///           (G+C=0, the DivideByZero boundary), alternating GC (skew 0),
///           and an extremely long sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The GC-skew contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// SEQ-GCSKEW-001 is strand-asymmetry GC skew, defined as
///     GC skew = (G − C) / (G + C),
/// counting only guanine and cytosine. The value lies in the CLOSED interval
/// [−1, 1]: an all-G run is the MAXIMUM (+1), an all-C run is the MINIMUM (−1),
/// and any sequence with equal G and C counts (e.g. alternating GCGC…) is exactly
/// 0. THE critical boundary is G + C = 0 (empty input, or a sequence with no G and
/// no C at all): the contract returns the DEFINED value 0, NOT a DivideByZero
/// crash and NOT a NaN — the implementation guards the zero denominator
/// (GcSkewCalculator.cs lines 38–45: `total > 0 ? (G−C)/total : 0`).
///   — docs/algorithms/Sequence_Composition/GC_Skew.md §2.2 (core model),
///     §2.4 (INV-01 range [−1,1]; INV-02 empty / no-G-or-C → 0), §6.1 (edge cases).
///     Sources: Lobry (1996), Grigoriev (1998), Wikipedia "GC skew".
///
/// API entry: GcSkewCalculator.CalculateGcSkew(...)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs).
///
/// SEQ-GCSKEW-001 has TWO documented surfaces with DIFFERENT validation contracts.
/// Fuzzing pins both, and the boundary between them, so neither drifts:
///
/// (1) The TYPED overload — CalculateGcSkew(DnaSequence)
///     (GcSkewCalculator.cs lines 21–25). The DnaSequence argument has ALREADY
///     passed the strict ctor validation gate, so only A/C/G/T (upper-cased) ever
///     reach the metric. Documented validation (GC_Skew.md §3.3):
///       • null DnaSequence            → ArgumentNullException (explicit guard,
///                                        never a NullReferenceException);
///       • empty sequence              → 0 (INV-02 zero-denominator guard);
///       • single G/C                  → +1 / −1; single A/T → 0 (no G/C present);
///       • alternating GC (equal G,C)  → exactly 0;
///       • valid extremely long input  → finite skew in [−1, 1], no overflow/hang.
///     Because non-ACGT cannot reach this overload, the result is always in the
///     DNA-bounded [−1, 1] interval (INV-01) and never NaN.
///
/// (2) The RAW-STRING overload — CalculateGcSkew(string)
///     (GcSkewCalculator.cs lines 30–36). LENIENT by documented design: it
///     short-circuits null/empty to 0, upper-cases, and counts ONLY G and C —
///     every other character (A, T, N, ambiguity codes, U, digits, gaps,
///     whitespace, null byte, unicode, surrogate halves) is IGNORED, not validated
///     (GC_Skew.md §5.2). So it must NEVER throw on arbitrary chars, and because
///     only G and C are ever counted the result stays a FINITE value in [−1, 1] —
///     ignored garbage does not shift the skew of the G/C bases that ARE present.
///     When NO G and NO C is present at all, the total==0 guard returns the defined
///     0 rather than a 0/0 NaN (the DivideByZero boundary). We pin that this surface
///     stays finite, non-throwing, and NaN-free on pure random-byte garbage, so the
///     strict/lenient boundary is explicit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-ATSKEW-001 — AT skew (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 210.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — balanced AT (equal A and T ⇒ skew 0),
///           all-A (skew +1), no AT present (A+T=0, the DivideByZero boundary),
///           and the "window edge" boundary — AT skew computed over caller
///           substring windows at window = 1, window &gt; sequence length, and
///           the partial trailing window.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The AT-skew contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// SEQ-ATSKEW-001 is strand A/T compositional asymmetry, defined as
///     AT skew = (A − T) / (A + T),
/// counting ONLY adenine and thymine (Charneski et al. 2011; Lobry 1996). The value
/// lies in the CLOSED interval [−1, 1]: an all-A run is the MAXIMUM (+1), an all-T
/// run is the MINIMUM (−1), and any sequence with equal A and T counts (e.g.
/// alternating ATAT…, or balanced ⇒ A = T) is exactly 0 (INV-02). THE critical
/// boundary is A + T = 0 (empty input, or a sequence with no A and no T at all,
/// e.g. "GGCC"): the contract returns the DEFINED value 0, NOT a DivideByZero crash
/// and NOT a NaN — the implementation guards the zero denominator
/// (GcSkewCalculator.cs line 205: `total > 0 ? (double)(A−T)/total : 0`).
///   — docs/algorithms/Extended_GC_Skew_Analysis/AT_Skew.md §2.2 (core model),
///     §2.4 (INV-01 range [−1,1]; INV-02 A=T ⇒ 0; INV-03 A+T=0 ⇒ 0, no exception/NaN;
///     INV-04 case-insensitive; INV-05 non-A/T ignored), §6.1 (edge cases),
///     §7.1 (worked examples: "AAAT" → 0.5, "AAATGGGCCC" → 0.5, pure A → +1,
///     pure T → −1). Sources: Lobry (1996), Charneski et al. (2011),
///     Biopython GC_skew zero-denominator convention.
///
/// "WINDOW EDGE" — the checklist names a window-edge boundary, but AT_Skew.md §5.3
/// (Not implemented) and §6.2 are EXPLICIT that this unit computes a single GLOBAL
/// scalar: windowed / cumulative AT-skew profiles are out of scope (users rely on
/// the GC-skew windowed methods for localization). There is therefore no native
/// AT-skew windowing API to fuzz. We honour the boundary the only spec-faithful way:
/// AT skew computed over CALLER-supplied substring windows — window length 1, a
/// window LONGER than the sequence (clamped by the caller), and the partial trailing
/// window — must equal the global formula applied to that exact slice and must stay
/// in [−1, 1] with no crash/NaN. This pins the window-edge boundary without inventing
/// a windowing semantic the unit does not document.
///
/// API entry: GcSkewCalculator.CalculateAtSkew(...)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs).
///
/// SEQ-ATSKEW-001 has TWO documented surfaces with DIFFERENT validation contracts.
/// Fuzzing pins both, and the boundary between them, so neither drifts:
///
/// (1) The TYPED overload — CalculateAtSkew(DnaSequence)
///     (GcSkewCalculator.cs lines 174–178). The DnaSequence argument has ALREADY
///     passed the strict ctor validation gate, so only A/C/G/T (upper-cased) ever
///     reach the metric. Documented validation (AT_Skew.md §3.3):
///       • null DnaSequence            → ArgumentNullException (explicit guard,
///                                        never a NullReferenceException);
///       • empty sequence              → 0 (INV-03 zero-denominator guard);
///       • single A/T                  → +1 / −1; single G/C → 0 (no A/T present);
///       • balanced / alternating AT   → exactly 0 (INV-02);
///       • valid extremely long input  → finite skew in [−1, 1], no overflow/hang.
///     Because non-ACGT cannot reach this overload, the result is always in the
///     DNA-bounded [−1, 1] interval (INV-01) and never NaN.
///
/// (2) The RAW-STRING overload — CalculateAtSkew(string)
///     (GcSkewCalculator.cs lines 189–195). LENIENT by documented design: it
///     short-circuits null/empty to 0, upper-cases, and counts ONLY A and T —
///     every other character (G, C, N, ambiguity codes, U, digits, gaps, whitespace,
///     null byte, unicode, surrogate halves) is IGNORED, not validated
///     (AT_Skew.md §3.3, INV-05). So it must NEVER throw on arbitrary chars, and
///     because only A and T are ever counted the result stays a FINITE value in
///     [−1, 1] — ignored garbage does not shift the skew of the A/T bases that ARE
///     present. When NO A and NO T is present at all, the total==0 guard returns the
///     defined 0 rather than a 0/0 NaN (the DivideByZero boundary). We pin that this
///     surface stays finite, non-throwing, and NaN-free on pure random-byte garbage,
///     so the strict/lenient boundary is explicit. NOTE: there is NO T↔U conversion
///     (AT_Skew.md §6.2): 'U' is ignored, not treated as T — we pin that too.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-RNACOMP-001 — RNA-specific per-base complement (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 212.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the empty span / empty destination boundary
///           of the whole-sequence composition, and the single-character minimum
///           (one base per recognized symbol).
///   • MC  = Malformed Content — out-of-alphabet (non-RNA) symbols, the "T instead
///           of U" case (T treated as U), and lowercase input (case-folding).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The RNA-complement contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// SEQ-RNACOMP-001 is the Watson–Crick complement of a SINGLE nucleotide IN THE
/// RNA ALPHABET — it emits U rather than T. The mapping is the Biopython
/// `ambiguous_rna_complement` lookup, IUPAC-complete:
///     A→U, U→A, C→G, G→C,  T→A,
///     R→Y, Y→R, S→S, W→W, K→M, M→K, B→V, V→B, D→H, H→D, N→N.
///   — docs/algorithms/Sequence_Composition/RNA_Complement.md §2.2, §4.2 (table),
///     Evidence SEQ-RNACOMP-001-Evidence.md (Biopython verbatim dicts).
///
/// THE CRUX of this unit (the three checklist boundaries) is how the edges are
/// handled, derived from the DOC/primary source, NOT from the code:
///   • "T instead of U" — T/t is TREATED AS U: T → A. Biopython builds the table
///     with `ambiguous_rna_complement["T"] = ...["U"]`, and `["U"] == "A"`
///     (RNA_Complement.md §2.2, INV-03, §6.1; Evidence pt. "T treated as U").
///     This is the RNA-vs-DNA distinction: the DNA sibling GetComplementBase maps
///     A→T, but RNA maps A→U and a stray T is silently read as a U.
///   • "non-RNA base" (MC) — characters OUTSIDE the IUPAC nucleotide set (gaps
///     '-'/'.', digits, 'Z', whitespace, the null byte, unicode letters, surrogate
///     halves) PASS THROUGH UNCHANGED, including their ORIGINAL CASE, and NEVER
///     throw (RNA_Complement.md §3.3, §6.1). This is a TOTAL function over char:
///     no exception is thrown for any input.
///   • "lowercase" — input is CASE-INSENSITIVE; every RECOGNIZED symbol is
///     normalized to UPPERCASE on output (repo convention §5.4, mirroring the DNA
///     sibling SEQ-COMP-001). So 'a' → 'U', 't' → 'A'. Unrecognized characters keep
///     their original case (so lowercase 'z' stays 'z').
///   • empty — there is no per-CHAR empty; the BE empty boundary is honoured on the
///     whole-sequence span composition (RNA complement of the empty span writes
///     nothing, never throws — RNA_Complement.md §6.2: whole-sequence complement is
///     composed by the caller).
///
/// API entry: SequenceExtensions.GetRnaComplementBase(char)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs lines
///    175–194). It is the single-char primitive; we additionally pin a
///    caller-composed whole-sequence RNA complement (the §6.2 documented usage) so
///    the empty-input BE boundary and the per-base fractions/identities of a
///    hand-checkable worked example (RNA_Complement.md §7.1) are exercised end to
///    end. Worked example (§7.1): RNA-complementing "ACGTUacgtuXYZxyz" under this
///    repo's uppercase convention yields "UGCAAUGCAAXRZxRz" — recognized bases
///    uppercased, X/Z/x/z (non-IUPAC) passed through verbatim.
///
/// Discipline: every expected value below is derived from the Biopython
/// `ambiguous_rna_complement` table and the §2.2/§6.1 T→A and case rules — NOT read
/// off the switch arms in the code. A test that would still pass against a DNA-style
/// A→T implementation, or against one that dropped the T→A rule, is invalid.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-GC-ANALYSIS-001 — comprehensive GC analysis (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 233.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — empty input (the zero-denominator / no-window
///           boundary), all-GC (GC=100, GcSkew bounds), all-AT (GC=0), non-ACGT
///           characters (ignored in both num &amp; denom), and a very long sequence
///           (O(n) scalars, no overflow / hang).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes);
///   docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The comprehensive-GC-analysis contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// SEQ-GC-ANALYSIS-001 is a single-pass AGGREGATION that bundles exact closed-form
/// composition statistics into a `GcAnalysisResult` record. For base counts
/// G, C, A, T over a sequence of length n
/// (Comprehensive_GC_Analysis.md §2.2, §3.2):
///   • OverallGcContent = (G+C)/(A+T+G+C)×100   — a PERCENTAGE in [0, 100] (INV-03),
///     repository/Brock convention (Biopython gc_fraction ×100) [3].
///   • OverallGcSkew    = (G−C)/(G+C), in [−1, +1] (INV-01); G+C=0 → 0 [1][2][5].
///   • OverallAtSkew    = (A−T)/(A+T), in [−1, +1] (INV-02); A+T=0 → 0 [6].
///   • WindowedGcSkew / WindowedGcContent — per-window values over every FULL window
///     of length `windowSize` advancing by `stepSize`; count = ⌊(n−w)/step⌋+1 when
///     n ≥ w, else 0 (INV-05) [5].
///   • GcSkewVariance / GcContentVariance — POPULATION variance σ²=Σ(xᵢ−μ)²/N (÷N,
///     not Bessel ÷N−1; the windows are the full population) of the windowed values,
///     both ≥ 0 (INV-04); 0 when there are no windows [7].
///   • SequenceLength = n.
/// Counting is CASE-INSENSITIVE; ONLY A/C/G/T contribute — every other symbol (N,
/// IUPAC ambiguity, U, digits, gaps, whitespace, null byte, unicode) is IGNORED in
/// BOTH numerator and denominator (Comprehensive_GC_Analysis.md §3.3, §6.2; matches
/// Biopython GC_skew ignoring ambiguous bases) [5].
///   — docs/algorithms/Extended_GC_Skew_Analysis/Comprehensive_GC_Analysis.md
///     §2.2 (model), §2.4 (INV-01..05), §3 (contract), §6.1 (edge cases),
///     §7.1 (worked examples). Sources: Lobry (1996), Grigoriev (1998),
///     Madigan &amp; Martinko / Brock, Charneski et al. (2011), Biopython, Cuemath.
///
/// SEQ-GC-ANALYSIS-001 has TWO documented surfaces with DIFFERENT validation
/// contracts. Fuzzing pins both, and the boundary between them, so neither drifts:
///
/// (1) The TYPED entry — AnalyzeGcContent(DnaSequence, windowSize, stepSize)
///     (GcSkewCalculator.cs lines 310–317). The DnaSequence argument has ALREADY
///     passed the strict ctor validation gate, so only A/C/G/T (upper-cased) ever
///     reach the metrics. Documented validation (Comprehensive_GC_Analysis.md §3.3,
///     §6.1): a null DnaSequence → ArgumentNullException (explicit guard, never a
///     NullReferenceException). Because non-ACGT cannot reach this surface, every
///     scalar lands in its DNA-bounded interval and never NaN.
///
/// (2) The LENIENT string entry — AnalyzeGcContent(string, windowSize, stepSize)
///     (GcSkewCalculator.cs lines 325–334). By design it short-circuits null/empty
///     to the ZERO result (all scalars 0, empty windowed lists, SequenceLength 0 —
///     §3.3, §6.1), upper-cases, and counts ONLY A/C/G/T — every other character is
///     IGNORED, not validated. So it must NEVER throw on arbitrary chars, every
///     scalar stays FINITE and in range, and ignored garbage does not shift the
///     statistics of the A/C/G/T bases that ARE present. When NO counted base is
///     present at all the zero-denominator guards return the DEFINED 0 rather than a
///     0/0 NaN (the DivideByZero boundary). We pin that this surface stays finite,
///     non-throwing and NaN-free on pure random-byte garbage.
///
/// Worked examples reproduced as positive-sanity tests (independently hand-derived
/// from §2.2, NOT read off the code — Comprehensive_GC_Analysis.md §7.1):
///   • "GGGCCAT" (G=3,C=2,A=1,T=1,n=7): GcContent = 5/7×100 = 71.42857142857143,
///     GcSkew = (3−2)/5 = 0.2, AtSkew = (1−1)/2 = 0.0.
///   • "GGCC", window 2 step 2 → windows GG (skew +1, GC% 100) and CC (skew −1,
///     GC% 100): GcSkewVariance = ((1−0)²+(−1−0)²)/2 = 1.0; GcContentVariance = 0.0.
/// A test that would still pass against an implementation that used the sample
/// variance (÷N−1), dropped the zero-denominator guard, or counted non-ACGT bases is
/// invalid.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class CompositionFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>Deterministic RNG — seed fixed so generated fuzz inputs are reproducible.</summary>
    private static readonly Random Rng = new(20260620);

    /// <summary>Generates a random valid DNA string of the given length over {A,C,G,T}.</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Generates a random string of arbitrary BMP code points (0x0000–0xFFFF),
    /// deliberately spanning control characters, the null byte, and lone surrogate
    /// halves — pure random-byte (RB) fuzz fodder for the validators.
    /// </summary>
    private static string RandomBmpChars(int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)Rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-GC-001 — GC content : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region BE — Boundary: empty string

    /// <summary>
    /// BE: the empty string is the lower size boundary. The public DnaSequence
    /// path defines this as an empty sequence whose GC% is 0 (no division by
    /// zero, no exception) — DnaSequence.cs lines 24–28; GC% definition has an
    /// empty denominator, resolved to 0 by the repository zero-division
    /// convention (GC_Content_Profile.md INV-05).
    /// </summary>
    [Test]
    public void GcContent_EmptyString_IsZeroAndDoesNotThrow()
    {
        var act = () =>
        {
            var seq = new DnaSequence(string.Empty);
            seq.Length.Should().Be(0);
            seq.GcContent().Should().Be(0.0,
                because: "an empty sequence has no bases; GC% is defined as 0 by the zero-division convention");
        };

        act.Should().NotThrow("the empty string is a defined boundary input, not an error");
    }

    #endregion

    #region BE — Boundary: single character

    /// <summary>
    /// BE: a one-base sequence is the minimal non-empty input. GC% is the binary
    /// extreme — exactly 100 for a single G/C, exactly 0 for a single A/T — with
    /// no rounding drift. Lowercase is accepted (case-folded) and yields the same
    /// value. Verified over all four bases in both cases.
    /// </summary>
    [TestCase('G', 100.0)]
    [TestCase('C', 100.0)]
    [TestCase('A', 0.0)]
    [TestCase('T', 0.0)]
    [TestCase('g', 100.0)]
    [TestCase('c', 100.0)]
    [TestCase('a', 0.0)]
    [TestCase('t', 0.0)]
    public void GcContent_SingleCharacter_IsBinaryExtreme(char baseChar, double expectedGc)
    {
        var seq = new DnaSequence(baseChar.ToString());

        seq.GcContent().Should().BeApproximately(expectedGc, Tolerance,
            because: $"a single '{baseChar}' is {(expectedGc == 100.0 ? "a GC base → 100%" : "an AT base → 0%")}");
    }

    #endregion

    #region INJ — Injection: non-ACGT characters

    /// <summary>
    /// INJ: characters that are not A/C/G/T after upper-casing must be rejected
    /// by the public DnaSequence path with the *documented, intentional*
    /// ArgumentException (DnaSequence.ValidateSequence, lines 112–124). This is
    /// the validation contract — NOT a silent miscount and NOT a raw runtime
    /// exception. Covers digits, whitespace, punctuation, the ambiguity code N,
    /// the RNA base U (DNA does not accept U), and an embedded null byte.
    /// </summary>
    [TestCase("N", TestName = "GcContent_NonAcgt_AmbiguityCodeN_Throws")]
    [TestCase("ACGTN", TestName = "GcContent_NonAcgt_TrailingN_Throws")]
    [TestCase("ACGU", TestName = "GcContent_NonAcgt_RnaBaseU_Throws")]
    [TestCase("ACGT123", TestName = "GcContent_NonAcgt_Digits_Throws")]
    [TestCase("AC GT", TestName = "GcContent_NonAcgt_Whitespace_Throws")]
    [TestCase("ACGT-", TestName = "GcContent_NonAcgt_Punctuation_Throws")]
    [TestCase("ACG\0T", TestName = "GcContent_NonAcgt_EmbeddedNullByte_Throws")]
    public void GcContent_NonAcgtCharacters_ThrowDocumentedArgumentException(string input)
    {
        var act = () => _ = new DnaSequence(input);

        act.Should().Throw<ArgumentException>(
            "non-ACGT input is rejected at construction by the documented validation gate, " +
            "not miscounted and not crashed");
    }

    #endregion

    #region INJ / BE — Injection: null

    /// <summary>
    /// INJ/BE: a null reference is the boundary of "no input". The public
    /// DnaSequence path defines null as an empty sequence (string.IsNullOrEmpty
    /// short-circuit, DnaSequence.cs lines 24–28) — so it must NOT throw
    /// NullReferenceException and GcContent() must be the defined 0. This pins
    /// that null is handled gracefully rather than crashing.
    /// </summary>
    [Test]
    public void GcContent_NullSequence_IsTreatedAsEmptyAndDoesNotThrow()
    {
        var act = () =>
        {
            var seq = new DnaSequence(null!);
            seq.Length.Should().Be(0);
            seq.GcContent().Should().Be(0.0);
        };

        act.Should().NotThrow<NullReferenceException>(
            "null must be handled by the documented IsNullOrEmpty gate, never dereferenced");
        act.Should().NotThrow(
            "null is a defined 'empty sequence' input on the public path, not an error");
    }

    #endregion

    #region INJ — Injection: unicode

    /// <summary>
    /// INJ: unicode injection — non-ASCII letters, combining diacritics,
    /// full-width look-alikes, and astral/surrogate-pair code points. None of
    /// these are A/C/G/T, so the public DnaSequence path must reject every one
    /// with the documented ArgumentException — never an
    /// IndexOutOfRange/encoding surprise from surrogate handling. The astral case
    /// (😀, a surrogate pair) specifically guards char-by-char validation against
    /// crashing on the high/low surrogate halves.
    /// </summary>
    [TestCase("ÀCGT", TestName = "GcContent_Unicode_AccentedLatin_Throws")]
    [TestCase("ACGTα", TestName = "GcContent_Unicode_GreekLetter_Throws")]
    [TestCase("ÁCGT", TestName = "GcContent_Unicode_CombiningAcute_Throws")]
    [TestCase("ＡＣＧＴ", TestName = "GcContent_Unicode_FullWidthLatin_Throws")]
    [TestCase("ACG😀T", TestName = "GcContent_Unicode_AstralSurrogatePair_Throws")]
    public void GcContent_UnicodeCharacters_ThrowDocumentedArgumentException(string input)
    {
        var act = () => _ = new DnaSequence(input);

        act.Should().Throw<ArgumentException>(
            "unicode characters are not valid nucleotides; the validation gate must reject them " +
            "via ArgumentException, including surrogate-pair (astral) code points");
    }

    #endregion

    #region BE — Boundary: extremely long

    /// <summary>
    /// BE/OVF: an extremely long valid sequence (1,000,000 bases) must compute
    /// without overflow, hang, or precision blow-up, and the result must stay in
    /// the closed range [0, 100]. The denominator is an int count; this guards
    /// that the (double)gc / valid * 100 arithmetic does not overflow or drift
    /// out of range at scale. A known-composition long input pins the exact
    /// value too.
    /// </summary>
    [Test]
    public void GcContent_ExtremelyLong_StaysInRangeAndDoesNotHang()
    {
        const int length = 1_000_000;

        // Known composition: exactly half GC — "AG" repeated is 50% A (AT base)
        // and 50% G (GC base), so GC% is precisely 50.
        var halfGc = new DnaSequence(string.Concat(Enumerable.Repeat("AG", length / 2)));
        halfGc.Length.Should().Be(length);
        halfGc.GcContent().Should().BeApproximately(50.0, Tolerance,
            because: "a sequence that is exactly half G/C has GC% = 50, at any length");

        // Fixed-seed random long sequence: result must be a valid percentage.
        var random = new DnaSequence(RandomDna(length));
        random.GcContent().Should().BeInRange(0.0, 100.0,
            because: "GC% is a proportion times 100; it can never escape [0, 100], even at scale");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-GC-001 — backing primitive SequenceExtensions.CalculateGcContent
    //  (documented as *lenient*: excludes non-A/T/G/C/U from num & denom)
    // ═══════════════════════════════════════════════════════════════════

    #region BE/INJ — backing span primitive contract

    /// <summary>
    /// The lenient backing primitive must, by its documented contract
    /// (SequenceExtensions.cs lines 17–58), return 0 — not NaN, not a crash — for
    /// the empty span and for input that contains no valid nucleotide at all
    /// (avoiding 0/0). Pins the zero-division convention at the primitive layer.
    /// </summary>
    [TestCase("", TestName = "CalculateGcContent_Empty_IsZero")]
    [TestCase("N", TestName = "CalculateGcContent_OnlyAmbiguity_IsZero")]
    [TestCase("12345", TestName = "CalculateGcContent_OnlyDigits_IsZero")]
    [TestCase("----", TestName = "CalculateGcContent_OnlyPunctuation_IsZero")]
    public void CalculateGcContent_NoValidNucleotide_IsZero(string input)
    {
        double result = input.AsSpan().CalculateGcContent();

        result.Should().Be(0.0,
            because: "the lenient primitive excludes invalid symbols and returns 0 when the denominator is empty");
    }

    /// <summary>
    /// The lenient primitive must exclude non-A/T/G/C/U symbols from BOTH the
    /// numerator and the denominator (Biopython "remove" mode). So "GC" with any
    /// amount of injected garbage (N, digits, whitespace, null byte) interspersed
    /// must still read as 100% GC, and the canonical Biopython case
    /// gc_fraction("ACTGN") = 0.50 → 50% must hold. This pins that injection does
    /// not silently corrupt the count at the primitive layer.
    /// </summary>
    [TestCase("ACTGN", 50.0, TestName = "CalculateGcContent_BiopythonAcgtnRemoveMode_Is50")]
    [TestCase("G N C ", 100.0, TestName = "CalculateGcContent_GarbageAroundGc_Is100")]
    [TestCase("AT12", 0.0, TestName = "CalculateGcContent_DigitsAroundAt_Is0")]
    [TestCase("g\0c", 100.0, TestName = "CalculateGcContent_NullByteBetweenGc_Is100")]
    [TestCase("acgu", 50.0, TestName = "CalculateGcContent_LowercaseWithU_Is50")]
    public void CalculateGcContent_ExcludesInvalidSymbols(string input, double expected)
    {
        double result = input.AsSpan().CalculateGcContent();

        result.Should().BeApproximately(expected, Tolerance,
            because: "non-A/T/G/C/U symbols are excluded from both numerator and denominator (remove mode)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-COMP-001 — DNA complement
    //  Strict public path: DnaSequence.Complement() (validates at construction)
    //  Lenient primitive:  SequenceExtensions.GetComplementBase / TryGetComplement
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-COMP-001 — DNA complement

    #region INJ — Injection: non-DNA characters (strict path rejects)

    /// <summary>
    /// INJ: characters that are not A/C/G/T after upper-casing must be rejected by
    /// the strict public DnaSequence path with the *documented, intentional*
    /// ArgumentException (DnaSequence ctor → ValidateSequence, lines 22–33 /
    /// 112–124) BEFORE any complement is taken — never a silent mis-complement and
    /// never a raw runtime exception. Covers digits, whitespace, punctuation/gap,
    /// the ambiguity code N, the RNA base U (DNA does not accept U), and an
    /// embedded null byte.
    /// </summary>
    [TestCase("N", TestName = "Complement_NonDna_AmbiguityCodeN_Throws")]
    [TestCase("ACGTN", TestName = "Complement_NonDna_TrailingN_Throws")]
    [TestCase("ACGU", TestName = "Complement_NonDna_RnaBaseU_Throws")]
    [TestCase("ACGT123", TestName = "Complement_NonDna_Digits_Throws")]
    [TestCase("AC GT", TestName = "Complement_NonDna_Whitespace_Throws")]
    [TestCase("ACGT-", TestName = "Complement_NonDna_GapDash_Throws")]
    [TestCase("ACG\0T", TestName = "Complement_NonDna_EmbeddedNullByte_Throws")]
    public void Complement_NonDnaCharacters_ThrowDocumentedArgumentException(string input)
    {
        var act = () => _ = new DnaSequence(input).Complement();

        act.Should().Throw<ArgumentException>(
            "non-A/C/G/T input is rejected at construction by the documented validation gate, " +
            "so the complement is never computed on garbage and never crashes");
    }

    #endregion

    #region BE — Boundary: empty string (strict path)

    /// <summary>
    /// BE: the empty string is the lower size boundary. The strict DnaSequence
    /// path defines it as an empty sequence (DnaSequence.cs lines 24–28); its
    /// complement is therefore the empty sequence — no division, no indexing, no
    /// exception. Complement of nothing is nothing.
    /// </summary>
    [Test]
    public void Complement_EmptyString_IsEmptyAndDoesNotThrow()
    {
        var act = () =>
        {
            var complement = new DnaSequence(string.Empty).Complement();
            complement.Length.Should().Be(0);
            complement.Sequence.Should().BeEmpty(
                because: "the complement of an empty sequence is the empty sequence");
        };

        act.Should().NotThrow("the empty string is a defined boundary input, not an error");
    }

    #endregion

    #region INJ / BE — Injection: null (strict path treats as empty)

    /// <summary>
    /// INJ/BE: a null reference is the boundary of "no input". The strict
    /// DnaSequence path defines null as an empty sequence (string.IsNullOrEmpty
    /// short-circuit, DnaSequence.cs lines 24–28), so Complement() must NOT throw
    /// NullReferenceException and must yield the empty sequence. Pins that null is
    /// handled gracefully rather than dereferenced.
    /// </summary>
    [Test]
    public void Complement_NullSequence_IsTreatedAsEmptyAndDoesNotThrow()
    {
        var act = () =>
        {
            var complement = new DnaSequence(null!).Complement();
            complement.Length.Should().Be(0);
            complement.Sequence.Should().BeEmpty();
        };

        act.Should().NotThrow<NullReferenceException>(
            "null must be handled by the documented IsNullOrEmpty gate, never dereferenced");
        act.Should().NotThrow(
            "null is a defined 'empty sequence' input on the public path, not an error");
    }

    #endregion

    #region INJ — Injection: mixed case (strict path accepts, case-folded)

    /// <summary>
    /// INJ: mixed and lower case input is accepted by the strict path — it is
    /// upper-cased before validation (ToUpperInvariant, DnaSequence.cs line 30) —
    /// and must complement IDENTICALLY to the uppercase form, always emitting
    /// uppercase A/T/G/C. This guards that case-folding neither rejects valid DNA
    /// nor corrupts the A↔T / C↔G mapping.
    /// </summary>
    [TestCase("acgt", "ACGT", TestName = "Complement_MixedCase_AllLower_FoldsAndComplements")]
    [TestCase("AcGt", "ACGT", TestName = "Complement_MixedCase_Alternating_FoldsAndComplements")]
    [TestCase("aCgT", "ACGT", TestName = "Complement_MixedCase_AlternatingInverse_FoldsAndComplements")]
    public void Complement_MixedCase_FoldsToUppercaseAndComplements(string input, string upper)
    {
        // The complement of A/C/G/T is T/G/C/A; computed once on the canonical
        // uppercase form, it must equal the complement of the mixed-case input.
        string expected = new DnaSequence(upper).Complement().Sequence;

        var complement = new DnaSequence(input).Complement();

        complement.Sequence.Should().Be(expected,
            because: "input is case-folded before complementing; case must not change the result");
        complement.Sequence.Should().MatchRegex("^[ACGT]*$",
            because: "the complement always emits uppercase canonical bases");
    }

    #endregion

    #region INJ — Injection: unicode (strict path rejects)

    /// <summary>
    /// INJ: unicode injection — accented Latin, Greek letters, combining
    /// diacritics, full-width look-alikes, and astral/surrogate-pair code points.
    /// None are A/C/G/T, so the strict DnaSequence path must reject every one with
    /// the documented ArgumentException — never an IndexOutOfRange/encoding
    /// surprise from surrogate handling. The astral case (😀, a surrogate pair)
    /// specifically guards char-by-char validation against crashing on the
    /// high/low surrogate halves before the complement is ever taken.
    /// </summary>
    [TestCase("ÀCGT", TestName = "Complement_Unicode_AccentedLatin_Throws")]
    [TestCase("ACGTα", TestName = "Complement_Unicode_GreekLetter_Throws")]
    [TestCase("ÁCGT", TestName = "Complement_Unicode_CombiningAcute_Throws")]
    [TestCase("ＡＣＧＴ", TestName = "Complement_Unicode_FullWidthLatin_Throws")]
    [TestCase("ACG😀T", TestName = "Complement_Unicode_AstralSurrogatePair_Throws")]
    public void Complement_UnicodeCharacters_ThrowDocumentedArgumentException(string input)
    {
        var act = () => _ = new DnaSequence(input).Complement();

        act.Should().Throw<ArgumentException>(
            "unicode characters are not valid nucleotides; the validation gate must reject them " +
            "via ArgumentException, including surrogate-pair (astral) code points");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Lenient primitive: GetComplementBase / span TryGetComplement
    //  (IUPAC-complete; non-IUPAC passes through unchanged; never throws)
    // ───────────────────────────────────────────────────────────────────

    #region INJ — lenient primitive: non-DNA / unicode pass through, never throw

    /// <summary>
    /// The lenient char primitive must NEVER throw on any char and, by its
    /// documented contract (RNA_Complement.md §3.3, the DNA sibling), must pass
    /// any non-IUPAC character through UNCHANGED. Recognized IUPAC symbols are
    /// complemented and emitted uppercase (A↔T, C↔G, U→A, N↔N, R↔Y). This pins
    /// that injection (digits, gap, whitespace, null byte, unicode letters,
    /// surrogate halves) cannot crash or silently corrupt the recognized mapping.
    /// </summary>
    [TestCase('A', 'T', TestName = "GetComplementBase_A_IsT")]
    [TestCase('C', 'G', TestName = "GetComplementBase_C_IsG")]
    [TestCase('U', 'A', TestName = "GetComplementBase_U_IsA")]
    [TestCase('N', 'N', TestName = "GetComplementBase_AmbiguityN_IsN")]
    [TestCase('R', 'Y', TestName = "GetComplementBase_AmbiguityR_IsY")]
    [TestCase('1', '1', TestName = "GetComplementBase_Digit_PassesThrough")]
    [TestCase('-', '-', TestName = "GetComplementBase_GapDash_PassesThrough")]
    [TestCase(' ', ' ', TestName = "GetComplementBase_Whitespace_PassesThrough")]
    [TestCase('\0', '\0', TestName = "GetComplementBase_NullByte_PassesThrough")]
    [TestCase('α', 'α', TestName = "GetComplementBase_GreekLetter_PassesThrough")]
    [TestCase('Z', 'Z', TestName = "GetComplementBase_NonIupacLetter_PassesThrough")]
    [TestCase('\uD83D', '\uD83D', TestName = "GetComplementBase_HighSurrogateHalf_PassesThrough")]
    [TestCase('\uDE00', '\uDE00', TestName = "GetComplementBase_LowSurrogateHalf_PassesThrough")]
    public void GetComplementBase_AnyChar_NeverThrowsAndPassesNonIupacThrough(char input, char expected)
    {
        char result = '￿';
        var act = () => result = SequenceExtensions.GetComplementBase(input);

        act.Should().NotThrow("the lenient primitive is total over char and never throws");
        result.Should().Be(expected,
            because: "recognized IUPAC symbols are complemented (uppercase); everything else passes through unchanged");
    }

    /// <summary>
    /// The lenient span complement must carry a mix of recognized and non-IUPAC
    /// garbage through char-by-char without throwing: recognized bases are
    /// complemented and uppercased, every other character (gap, digit, null byte)
    /// is preserved verbatim. This pins that injected garbage interspersed with
    /// real bases neither crashes nor shifts the complement of the valid bases.
    /// </summary>
    [Test]
    public void TryGetComplement_GarbageInterspersed_ComplementsBasesAndPreservesGarbage()
    {
        const string input = "aC-G1N\0T";       // mixed case, gap, digit, N, null byte
        Span<char> destination = new char[input.Length];

        bool ok = input.AsSpan().TryGetComplement(destination);

        ok.Should().BeTrue("the destination is exactly the source length");
        new string(destination).Should().Be("TG-C1N\0A",
            because: "A→T C→G G→C T→A N→N (uppercase), and non-IUPAC '-','1','\\0' pass through unchanged");
    }

    #endregion

    #region BE — lenient primitive: empty span

    /// <summary>
    /// BE: the empty span is the lower size boundary for the lenient primitive.
    /// TryGetComplement must succeed (destination length ≥ source length holds
    /// trivially) and write nothing — no exception, no hang, no out-of-range.
    /// </summary>
    [Test]
    public void TryGetComplement_EmptySpan_SucceedsAndWritesNothing()
    {
        var act = () =>
        {
            bool ok = ReadOnlySpan<char>.Empty.TryGetComplement(Span<char>.Empty);
            ok.Should().BeTrue("an empty complement always fits an empty destination");
        };

        act.Should().NotThrow("the empty span is a defined boundary, not an error");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-REVCOMP-001 — reverse complement
    //  Strict public path: DnaSequence.ReverseComplement() (validates at ctor)
    //  Lenient string:     DnaSequence.GetReverseComplementString (never throws)
    //  Lenient primitive:  ReadOnlySpan<char>.TryGetReverseComplement
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-REVCOMP-001 — reverse complement

    #region INJ — Injection: non-DNA characters (strict path rejects)

    /// <summary>
    /// INJ: characters that are not A/C/G/T after upper-casing must be rejected by
    /// the strict public DnaSequence path with the *documented, intentional*
    /// ArgumentException (DnaSequence ctor → ValidateSequence, lines 22–33 /
    /// 112–124) BEFORE any reverse-complement is taken — never a silent
    /// mis-complement and never a raw runtime exception. Covers digits, whitespace,
    /// punctuation/gap, the ambiguity code N, the RNA base U (DNA does not accept
    /// U), and an embedded null byte.
    /// </summary>
    [TestCase("N", TestName = "ReverseComplement_NonDna_AmbiguityCodeN_Throws")]
    [TestCase("ACGTN", TestName = "ReverseComplement_NonDna_TrailingN_Throws")]
    [TestCase("ACGU", TestName = "ReverseComplement_NonDna_RnaBaseU_Throws")]
    [TestCase("ACGT123", TestName = "ReverseComplement_NonDna_Digits_Throws")]
    [TestCase("AC GT", TestName = "ReverseComplement_NonDna_Whitespace_Throws")]
    [TestCase("ACGT-", TestName = "ReverseComplement_NonDna_GapDash_Throws")]
    [TestCase("ACG\0T", TestName = "ReverseComplement_NonDna_EmbeddedNullByte_Throws")]
    public void ReverseComplement_NonDnaCharacters_ThrowDocumentedArgumentException(string input)
    {
        var act = () => _ = new DnaSequence(input).ReverseComplement();

        act.Should().Throw<ArgumentException>(
            "non-A/C/G/T input is rejected at construction by the documented validation gate, " +
            "so the reverse complement is never computed on garbage and never crashes");
    }

    #endregion

    #region BE — Boundary: empty string (strict path)

    /// <summary>
    /// BE: the empty string is the lower size boundary. The strict DnaSequence path
    /// defines it as an empty sequence (DnaSequence.cs lines 24–28); its reverse
    /// complement is therefore the empty sequence — no division, no indexing, no
    /// exception. Reverse complement of nothing is nothing.
    /// </summary>
    [Test]
    public void ReverseComplement_EmptyString_IsEmptyAndDoesNotThrow()
    {
        var act = () =>
        {
            var revComp = new DnaSequence(string.Empty).ReverseComplement();
            revComp.Length.Should().Be(0);
            revComp.Sequence.Should().BeEmpty(
                because: "the reverse complement of an empty sequence is the empty sequence");
        };

        act.Should().NotThrow("the empty string is a defined boundary input, not an error");
    }

    #endregion

    #region INJ / BE — Injection: null (strict path treats as empty)

    /// <summary>
    /// INJ/BE: a null reference is the boundary of "no input". The strict
    /// DnaSequence path defines null as an empty sequence (string.IsNullOrEmpty
    /// short-circuit, DnaSequence.cs lines 24–28), so ReverseComplement() must NOT
    /// throw NullReferenceException and must yield the empty sequence. Pins that
    /// null is handled gracefully rather than dereferenced.
    /// </summary>
    [Test]
    public void ReverseComplement_NullSequence_IsTreatedAsEmptyAndDoesNotThrow()
    {
        var act = () =>
        {
            var revComp = new DnaSequence(null!).ReverseComplement();
            revComp.Length.Should().Be(0);
            revComp.Sequence.Should().BeEmpty();
        };

        act.Should().NotThrow<NullReferenceException>(
            "null must be handled by the documented IsNullOrEmpty gate, never dereferenced");
        act.Should().NotThrow(
            "null is a defined 'empty sequence' input on the public path, not an error");
    }

    #endregion

    #region BE — Boundary: single character (strict path)

    /// <summary>
    /// BE: a one-base sequence is the minimal non-empty input. With a single base
    /// the reverse is a no-op, so the reverse complement is exactly the complement:
    /// A→T, C→G, G→C, T→A. Lowercase is accepted (case-folded) and yields the same
    /// uppercase result. Verified over all four bases in both cases. This pins that
    /// the length-1 boundary neither off-by-ones the index nor skips the
    /// complement.
    /// </summary>
    [TestCase('A', "T")]
    [TestCase('C', "G")]
    [TestCase('G', "C")]
    [TestCase('T', "A")]
    [TestCase('a', "T")]
    [TestCase('c', "G")]
    [TestCase('g', "C")]
    [TestCase('t', "A")]
    public void ReverseComplement_SingleCharacter_IsItsComplement(char baseChar, string expected)
    {
        var revComp = new DnaSequence(baseChar.ToString()).ReverseComplement();

        revComp.Sequence.Should().Be(expected,
            because: $"with one base the reverse is a no-op, so revcomp('{baseChar}') is its complement '{expected}'");
    }

    #endregion

    #region INJ — Injection: unicode (strict path rejects)

    /// <summary>
    /// INJ: unicode injection — accented Latin, Greek letters, combining
    /// diacritics, full-width look-alikes, and astral/surrogate-pair code points.
    /// None are A/C/G/T, so the strict DnaSequence path must reject every one with
    /// the documented ArgumentException — never an IndexOutOfRange/encoding surprise
    /// from surrogate handling. The astral case (😀, a surrogate pair) specifically
    /// guards char-by-char validation against crashing on the high/low surrogate
    /// halves before the reverse complement is ever taken.
    /// </summary>
    [TestCase("ÀCGT", TestName = "ReverseComplement_Unicode_AccentedLatin_Throws")]
    [TestCase("ACGTα", TestName = "ReverseComplement_Unicode_GreekLetter_Throws")]
    [TestCase("ÁCGT", TestName = "ReverseComplement_Unicode_CombiningAcute_Throws")]
    [TestCase("ＡＣＧＴ", TestName = "ReverseComplement_Unicode_FullWidthLatin_Throws")]
    [TestCase("ACG😀T", TestName = "ReverseComplement_Unicode_AstralSurrogatePair_Throws")]
    public void ReverseComplement_UnicodeCharacters_ThrowDocumentedArgumentException(string input)
    {
        var act = () => _ = new DnaSequence(input).ReverseComplement();

        act.Should().Throw<ArgumentException>(
            "unicode characters are not valid nucleotides; the validation gate must reject them " +
            "via ArgumentException, including surrogate-pair (astral) code points");
    }

    #endregion

    #region Robustness — involution on fuzzed-but-valid inputs (strict path)

    /// <summary>
    /// Robustness: reverse complement is an INVOLUTION — revcomp(revcomp(x)) == x
    /// for any valid sequence. Asserted over deterministic fuzzed valid DNA across
    /// sizes (including the single-base and odd/even-length boundaries) and over a
    /// known fixed case. This is the theory-correct contract: applying the
    /// transform twice must return the exact original, with no drift, truncation,
    /// or off-by-one from the reverse indexing.
    /// </summary>
    [TestCase(1, TestName = "ReverseComplement_Involution_Len1")]
    [TestCase(2, TestName = "ReverseComplement_Involution_Len2")]
    [TestCase(7, TestName = "ReverseComplement_Involution_Len7_Odd")]
    [TestCase(64, TestName = "ReverseComplement_Involution_Len64")]
    [TestCase(1000, TestName = "ReverseComplement_Involution_Len1000")]
    public void ReverseComplement_AppliedTwice_IsIdentity(int length)
    {
        var original = new DnaSequence(RandomDna(length));

        var doubleRevComp = original.ReverseComplement().ReverseComplement();

        doubleRevComp.Sequence.Should().Be(original.Sequence,
            because: "reverse complement is an involution: applying it twice returns the original");
    }

    /// <summary>
    /// Robustness: a known, fully-pinned reverse-complement case. revcomp("ATGC")
    /// = complement "TACG" read 5'→3' (reversed) = "GCAT". Pins both the A↔T/C↔G
    /// mapping AND the reversal direction together, so neither can drift
    /// independently.
    /// </summary>
    [Test]
    public void ReverseComplement_KnownCase_IsComplementReversed()
    {
        var revComp = new DnaSequence("ATGC").ReverseComplement();

        revComp.Sequence.Should().Be("GCAT",
            because: "complement of ATGC is TACG; read 5'→3' (reversed) it is GCAT");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Lenient surfaces: GetReverseComplementString / span TryGetReverseComplement
    //  (IUPAC-complete; non-IUPAC passes through unchanged; never throw)
    // ───────────────────────────────────────────────────────────────────

    #region BE — lenient string helper: null and empty pass through

    /// <summary>
    /// BE: the lenient static string helper returns null/empty input verbatim
    /// (DnaSequence.cs lines 151–152) — no NullReferenceException, no exception, no
    /// hang. Pins that the "no input" boundary is a defined pass-through, not a
    /// crash.
    /// </summary>
    [Test]
    public void GetReverseComplementString_NullAndEmpty_ReturnInputAndDoNotThrow()
    {
        string? nullResult = null;
        var actNull = () => nullResult = DnaSequence.GetReverseComplementString(null!);
        actNull.Should().NotThrow<NullReferenceException>(
            "the lenient helper short-circuits null via IsNullOrEmpty, never dereferencing it");
        actNull.Should().NotThrow();
        nullResult.Should().BeNull("null is returned verbatim");

        DnaSequence.GetReverseComplementString(string.Empty).Should().BeEmpty(
            because: "the reverse complement of the empty string is the empty string");
    }

    #endregion

    #region INJ — lenient string helper: non-DNA / unicode pass through, never throw

    /// <summary>
    /// INJ: the lenient string helper NEVER throws and passes any non-IUPAC
    /// character through UNCHANGED while complementing recognized IUPAC bases
    /// (uppercased) and reversing the whole thing. So garbage interspersed with
    /// real bases neither crashes nor shifts the reverse complement of the valid
    /// bases. For "aC-G1N\0T": complement char-by-char is T,G,-,C,1,N,\0,A; reading
    /// that 5'→3' (reversed) gives "A\0N1C-GT". This pins that injection cannot
    /// corrupt the mapping or the reversal on this surface, including the embedded
    /// null byte. A unicode letter (Greek α) and a lone surrogate half likewise
    /// pass through unchanged without an encoding crash.
    /// </summary>
    [TestCase("aC-G1N\0T", "A\0N1C-GT", TestName = "GetReverseComplementString_GarbageInterspersed_PreservesAndReverses")]
    [TestCase("Gαc", "GαC", TestName = "GetReverseComplementString_GreekLetter_PassesThrough")]
    [TestCase("G\uD83Dc", "G\uD83DC", TestName = "GetReverseComplementString_LoneHighSurrogate_PassesThrough")]
    public void GetReverseComplementString_NonDnaAndUnicode_PassThroughAndNeverThrow(string input, string expected)
    {
        string result = "￿";
        var act = () => result = DnaSequence.GetReverseComplementString(input);

        act.Should().NotThrow("the lenient string helper is total over char and never throws");
        result.Should().Be(expected,
            because: "recognized bases are complemented (uppercase) and the string reversed; " +
                     "non-IUPAC characters pass through unchanged");
    }

    #endregion

    #region BE — lenient span primitive: empty span

    /// <summary>
    /// BE: the empty span is the lower size boundary for the lenient span primitive.
    /// TryGetReverseComplement must succeed (destination length ≥ source length
    /// holds trivially) and write nothing — no exception, no hang, no out-of-range.
    /// </summary>
    [Test]
    public void TryGetReverseComplement_EmptySpan_SucceedsAndWritesNothing()
    {
        var act = () =>
        {
            bool ok = ReadOnlySpan<char>.Empty.TryGetReverseComplement(Span<char>.Empty);
            ok.Should().BeTrue("an empty reverse complement always fits an empty destination");
        };

        act.Should().NotThrow("the empty span is a defined boundary, not an error");
    }

    /// <summary>
    /// INJ: the lenient span primitive carries a mix of recognized and non-IUPAC
    /// garbage through char-by-char without throwing: recognized bases are
    /// complemented (uppercased), every other character (gap, digit, null byte) is
    /// preserved verbatim, and the whole result is reversed. Pins that injected
    /// garbage neither crashes nor shifts the reverse complement of the valid bases
    /// on the span surface, mirroring the string helper.
    /// </summary>
    [Test]
    public void TryGetReverseComplement_GarbageInterspersed_ComplementsBasesReversedAndPreservesGarbage()
    {
        const string input = "aC-G1N\0T";   // mixed case, gap, digit, N, null byte
        Span<char> destination = new char[input.Length];

        bool ok = input.AsSpan().TryGetReverseComplement(destination);

        ok.Should().BeTrue("the destination is exactly the source length");
        new string(destination).Should().Be("A\0N1C-GT",
            because: "each base is complemented (uppercase) and the sequence reversed; " +
                     "non-IUPAC '-','1','\\0' pass through unchanged");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-VALID-001 — sequence validation
    //  Total predicates:  SequenceExtensions.IsValidDna / IsValidRna (never throw)
    //  Factory:           DnaSequence.TryCreate (false+null only on ArgumentException)
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-VALID-001 — sequence validation

    #region BE — Boundary: empty string

    /// <summary>
    /// BE: the empty span is the lower size boundary. Both validators return true
    /// by VACUOUS TRUTH — no invalid character is present (Sequence_Validation.md
    /// §3.3, §6.1; INV-03). No division, no indexing, no exception.
    /// </summary>
    [Test]
    public void IsValid_EmptyString_IsTrueForBothAlphabets()
    {
        var act = () =>
        {
            ReadOnlySpan<char>.Empty.IsValidDna().Should().BeTrue(
                "an empty sequence has no invalid characters (vacuous truth)");
            ReadOnlySpan<char>.Empty.IsValidRna().Should().BeTrue(
                "an empty sequence has no invalid characters (vacuous truth)");
        };

        act.Should().NotThrow("the empty span is a defined boundary input, not an error");
    }

    /// <summary>
    /// BE: TryCreate on null/empty input is NOT a validation failure — the
    /// DnaSequence ctor short-circuits "no input" to an empty sequence
    /// (DnaSequence.cs lines 24–28), so TryCreate returns true with a non-null,
    /// length-0 sequence and never throws NullReferenceException
    /// (Sequence_Validation.md §3.3).
    /// </summary>
    [TestCase("", TestName = "TryCreate_EmptyString_SucceedsWithEmptySequence")]
    [TestCase(null, TestName = "TryCreate_Null_SucceedsWithEmptySequence")]
    public void TryCreate_NullOrEmpty_SucceedsWithEmptySequence(string? input)
    {
        bool created = false;
        DnaSequence? result = null;
        var act = () => created = DnaSequence.TryCreate(input!, out result);

        act.Should().NotThrow<NullReferenceException>(
            "null must be handled by the documented IsNullOrEmpty gate, never dereferenced");
        act.Should().NotThrow();
        created.Should().BeTrue("'no input' is a defined empty sequence, not a validation failure");
        result.Should().NotBeNull();
        result!.Length.Should().Be(0);
    }

    #endregion

    #region BE — Boundary: single character

    /// <summary>
    /// BE: a one-character input is the minimal non-empty case. Each unambiguous
    /// base is valid for exactly the alphabet that contains it, case-insensitively
    /// (char.ToUpperInvariant per char). This pins the documented DNA/RNA asymmetry
    /// at the smallest scale: A/C/G are valid for both; T is DNA-only; U is RNA-only
    /// (Sequence_Validation.md §5.2 table).
    /// </summary>
    [TestCase('A', true, true)]
    [TestCase('C', true, true)]
    [TestCase('G', true, true)]
    [TestCase('T', true, false)]
    [TestCase('U', false, true)]
    [TestCase('a', true, true)]
    [TestCase('t', true, false)]
    [TestCase('u', false, true)]
    public void IsValid_SingleCharacter_MatchesItsAlphabet(char c, bool validDna, bool validRna)
    {
        ReadOnlySpan<char> span = stackalloc char[] { c };

        span.IsValidDna().Should().Be(validDna,
            because: $"'{c}' is {(validDna ? "" : "not ")}a DNA base (A/C/G/T, case-insensitive)");
        span.IsValidRna().Should().Be(validRna,
            because: $"'{c}' is {(validRna ? "" : "not ")}an RNA base (A/C/G/U, case-insensitive)");
    }

    #endregion

    #region INJ — Injection: mixed case (validators fold to uppercase)

    /// <summary>
    /// INJ: mixed and lower case must be accepted because validation folds each
    /// char with char.ToUpperInvariant before the membership test
    /// (Sequence_Validation.md §3.3, "case-insensitive"). Lowercase/mixed a-c-g-t
    /// is valid DNA and TryCreate materializes it; a-c-g-u is valid RNA. Case must
    /// neither reject valid bases nor flip the classification.
    /// </summary>
    [TestCase("acgt", true, false, TestName = "IsValid_MixedCase_LowerAcgt_DnaOnly")]
    [TestCase("AcGt", true, false, TestName = "IsValid_MixedCase_AlternatingAcgt_DnaOnly")]
    [TestCase("acgu", false, true, TestName = "IsValid_MixedCase_LowerAcgu_RnaOnly")]
    [TestCase("aCgU", false, true, TestName = "IsValid_MixedCase_AlternatingAcgu_RnaOnly")]
    public void IsValid_MixedCase_FoldsBeforeMembershipTest(string input, bool validDna, bool validRna)
    {
        input.AsSpan().IsValidDna().Should().Be(validDna,
            because: "characters are upper-cased before the A/C/G/T membership test");
        input.AsSpan().IsValidRna().Should().Be(validRna,
            because: "characters are upper-cased before the A/C/G/U membership test");
    }

    [TestCase("acgt", TestName = "TryCreate_MixedCase_LowerAcgt_SucceedsUppercased")]
    [TestCase("AcGt", TestName = "TryCreate_MixedCase_AlternatingAcgt_SucceedsUppercased")]
    public void TryCreate_MixedCaseDna_SucceedsAndStoresUppercase(string input)
    {
        bool created = DnaSequence.TryCreate(input, out var result);

        created.Should().BeTrue("lowercase/mixed-case A/C/G/T is valid DNA after case folding");
        result.Should().NotBeNull();
        result!.Sequence.Should().Be(input.ToUpperInvariant(),
            because: "the ctor normalizes accepted input to uppercase (DnaSequence.cs line 30)");
    }

    #endregion

    #region INJ — Injection: non-ASCII / unicode (rejected, never throw)

    /// <summary>
    /// INJ: non-ASCII letters, combining diacritics, full-width look-alikes, Greek
    /// letters, and astral/surrogate-pair code points are NOT in either alphabet.
    /// Both total predicates must return false WITHOUT throwing — char-by-char
    /// folding makes a surrogate half simply "not a base", never an encoding crash
    /// or IndexOutOfRange. The astral case (😀) specifically guards the high/low
    /// surrogate halves of a pair.
    /// </summary>
    [TestCase("ÀCGT", TestName = "IsValid_Unicode_AccentedLatin_RejectedNoThrow")]
    [TestCase("ACGTα", TestName = "IsValid_Unicode_GreekLetter_RejectedNoThrow")]
    [TestCase("ÁCGT", TestName = "IsValid_Unicode_CombiningAcute_RejectedNoThrow")]
    [TestCase("ＡＣＧＴ", TestName = "IsValid_Unicode_FullWidthLatin_RejectedNoThrow")]
    [TestCase("ACG😀T", TestName = "IsValid_Unicode_AstralSurrogatePair_RejectedNoThrow")]
    public void IsValid_UnicodeCharacters_AreRejectedAndNeverThrow(string input)
    {
        bool dna = true, rna = true;
        var act = () =>
        {
            dna = input.AsSpan().IsValidDna();
            rna = input.AsSpan().IsValidRna();
        };

        act.Should().NotThrow(
            "the validators are total over char, including surrogate halves — never an encoding crash");
        dna.Should().BeFalse("unicode characters are not A/C/G/T");
        rna.Should().BeFalse("unicode characters are not A/C/G/U");
    }

    /// <summary>
    /// INJ: the same unicode injection routed through the factory must yield
    /// false+null (the ctor's documented ArgumentException, caught by TryCreate),
    /// NOT a leaked exception. Pins that astral/surrogate input does not escape the
    /// ArgumentException-only catch in TryCreate (DnaSequence.cs lines 129–141).
    /// </summary>
    [TestCase("ÀCGT", TestName = "TryCreate_Unicode_AccentedLatin_FalseNull")]
    [TestCase("ACGTα", TestName = "TryCreate_Unicode_GreekLetter_FalseNull")]
    [TestCase("ＡＣＧＴ", TestName = "TryCreate_Unicode_FullWidthLatin_FalseNull")]
    [TestCase("ACG😀T", TestName = "TryCreate_Unicode_AstralSurrogatePair_FalseNull")]
    public void TryCreate_UnicodeCharacters_ReturnFalseAndNullWithoutLeakingException(string input)
    {
        bool created = true;
        DnaSequence? result = new DnaSequence("A");
        var act = () => created = DnaSequence.TryCreate(input, out result);

        act.Should().NotThrow(
            "TryCreate converts the documented ArgumentException to false+null; no other exception leaks");
        created.Should().BeFalse("unicode input is invalid DNA");
        result.Should().BeNull("a failed validation yields a null result");
    }

    #endregion

    #region INJ — Injection: null bytes and control characters (rejected, never throw)

    /// <summary>
    /// INJ: the null byte and the ASCII control characters (TAB, LF, CR, BEL, ESC,
    /// DEL) are outside both alphabets. The validators must reject them and never
    /// throw — these are classic crash triggers for naive char handling. Both an
    /// isolated control char and one embedded between real bases are covered.
    /// </summary>
    [TestCase("\0", TestName = "IsValid_Control_NullByteAlone_Rejected")]
    [TestCase("ACG\0T", TestName = "IsValid_Control_EmbeddedNullByte_Rejected")]
    [TestCase("AC\tGT", TestName = "IsValid_Control_Tab_Rejected")]
    [TestCase("AC\nGT", TestName = "IsValid_Control_LineFeed_Rejected")]
    [TestCase("AC\rGT", TestName = "IsValid_Control_CarriageReturn_Rejected")]
    [TestCase("AC\aGT", TestName = "IsValid_Control_Bell_Rejected")]
    [TestCase("AC\u001BGT", TestName = "IsValid_Control_Escape_Rejected")]
    [TestCase("AC\u007FGT", TestName = "IsValid_Control_Delete_Rejected")]
    public void IsValid_ControlCharacters_AreRejectedAndNeverThrow(string input)
    {
        bool dna = true, rna = true;
        var act = () =>
        {
            dna = input.AsSpan().IsValidDna();
            rna = input.AsSpan().IsValidRna();
        };

        act.Should().NotThrow("control characters and null bytes must not crash the per-char scan");
        dna.Should().BeFalse("control characters are not A/C/G/T");
        rna.Should().BeFalse("control characters are not A/C/G/U");
    }

    /// <summary>
    /// INJ: the factory must reject control/null-byte input via false+null without
    /// leaking any non-ArgumentException. The embedded null byte specifically guards
    /// against C-string truncation surprises in the validation message path.
    /// </summary>
    [TestCase("ACG\0T", TestName = "TryCreate_Control_EmbeddedNullByte_FalseNull")]
    [TestCase("AC\tGT", TestName = "TryCreate_Control_Tab_FalseNull")]
    [TestCase("AC\u001BGT", TestName = "TryCreate_Control_Escape_FalseNull")]
    public void TryCreate_ControlCharacters_ReturnFalseAndNullWithoutLeakingException(string input)
    {
        bool created = true;
        DnaSequence? result = new DnaSequence("A");
        var act = () => created = DnaSequence.TryCreate(input, out result);

        act.Should().NotThrow("TryCreate must surface invalid control input as false+null, not a leak");
        created.Should().BeFalse();
        result.Should().BeNull();
    }

    #endregion

    #region INJ — Injection: IUPAC ambiguity codes and gap (strict mode rejects)

    /// <summary>
    /// INJ: strict mode rejects IUPAC ambiguity codes (N,R,Y,S,W,K,M,B,D,H,V) and
    /// the gap '-', even though the IUPAC standard defines them
    /// (Sequence_Validation.md §5.2–5.4, INV-01..02). A single ambiguity code makes
    /// the whole sequence invalid for BOTH alphabets; this pins the strict-mode
    /// deviation so it cannot silently widen.
    /// </summary>
    [TestCase("N", TestName = "IsValid_Iupac_AnyN_Rejected")]
    [TestCase("ACGTN", TestName = "IsValid_Iupac_TrailingN_Rejected")]
    [TestCase("R", TestName = "IsValid_Iupac_PurineR_Rejected")]
    [TestCase("Y", TestName = "IsValid_Iupac_PyrimidineY_Rejected")]
    [TestCase("ACGT-", TestName = "IsValid_Iupac_GapDash_Rejected")]
    public void IsValid_IupacAmbiguityAndGap_AreRejectedInStrictMode(string input)
    {
        input.AsSpan().IsValidDna().Should().BeFalse(
            "strict DNA validation does not accept IUPAC ambiguity codes or the gap");
        input.AsSpan().IsValidRna().Should().BeFalse(
            "strict RNA validation does not accept IUPAC ambiguity codes or the gap");
    }

    #endregion

    #region INJ — DNA/RNA alphabet asymmetry (T vs U cannot cross over)

    /// <summary>
    /// INJ: the documented DNA/RNA asymmetry must hold exactly (Sequence_Validation.md
    /// §5.2 table). "ACGT" is valid DNA but INVALID RNA (T ∉ RNA); "ACGU" is valid
    /// RNA but INVALID DNA (U ∉ DNA). A sequence mixing T and U is invalid for both.
    /// This pins that neither alphabet leaks into the other.
    /// </summary>
    [TestCase("ACGT", true, false, TestName = "IsValid_Asymmetry_Acgt_DnaOnly")]
    [TestCase("ACGU", false, true, TestName = "IsValid_Asymmetry_Acgu_RnaOnly")]
    [TestCase("ACGTU", false, false, TestName = "IsValid_Asymmetry_MixedTandU_NeitherAlphabet")]
    public void IsValid_DnaRnaAlphabets_DoNotCrossOver(string input, bool validDna, bool validRna)
    {
        input.AsSpan().IsValidDna().Should().Be(validDna,
            because: "T is a DNA base; U is not — the alphabets are disjoint on T/U");
        input.AsSpan().IsValidRna().Should().Be(validRna,
            because: "U is an RNA base; T is not — the alphabets are disjoint on T/U");
    }

    #endregion

    #region BE — Boundary: extremely long (no hang, classifies consistently)

    /// <summary>
    /// BE/OVF: an extremely long valid sequence (1,000,000 bases) must validate
    /// without hang or overflow — the scan is O(n), O(1) space. A long valid input
    /// classifies true; flipping a single buried character to garbage must flip the
    /// result to false (the scan reaches it), proving the predicate does not bail
    /// early or short-circuit incorrectly at scale.
    /// </summary>
    [Test]
    public void IsValid_ExtremelyLong_DoesNotHangAndClassifiesConsistently()
    {
        const int length = 1_000_000;

        var longValid = RandomDna(length);
        bool dnaValid = true;
        var act = () => dnaValid = longValid.AsSpan().IsValidDna();
        act.Should().NotThrow("a long valid sequence must not overflow or hang");
        dnaValid.Should().BeTrue("a million A/C/G/T characters are all valid DNA");

        // Inject one invalid character deep in the interior: the full scan must
        // still reach it and return false.
        var corrupted = longValid.ToCharArray();
        corrupted[length / 2] = 'N';
        new string(corrupted).AsSpan().IsValidDna().Should().BeFalse(
            "a single buried invalid character must be detected even at scale");

        DnaSequence.TryCreate(longValid, out var created).Should().BeTrue(
            "the long valid sequence materializes through the factory without leaking");
        created!.Length.Should().Be(length);
    }

    #endregion

    #region RB — Random-byte sweeps (never throw; classify valid iff all-in-alphabet)

    /// <summary>
    /// RB: a fixed-seed sweep of random BMP code points — deliberately including
    /// control characters, null bytes and lone surrogate halves — must NEVER throw
    /// from either predicate, and the classification must be EXACTLY equivalent to
    /// the independent oracle "every char ∈ alphabet (case-folded)". Random garbage
    /// is overwhelmingly invalid; the point is total, crash-free, consistent
    /// classification rather than any particular verdict.
    /// </summary>
    [Test]
    public void IsValid_RandomBmpBytes_NeverThrowAndMatchMembershipOracle()
    {
        const string dnaAlphabet = "ACGT";
        const string rnaAlphabet = "ACGU";

        for (int trial = 0; trial < 2000; trial++)
        {
            string input = RandomBmpChars(Rng.Next(0, 33));

            bool dna = false, rna = false;
            var act = () =>
            {
                dna = input.AsSpan().IsValidDna();
                rna = input.AsSpan().IsValidRna();
            };
            act.Should().NotThrow(
                $"the validators are total over any char sequence; offending input: {Describe(input)}");

            bool oracleDna = input.All(ch => dnaAlphabet.Contains(char.ToUpperInvariant(ch)));
            bool oracleRna = input.All(ch => rnaAlphabet.Contains(char.ToUpperInvariant(ch)));

            dna.Should().Be(oracleDna,
                because: $"IsValidDna must equal the membership oracle; offending input: {Describe(input)}");
            rna.Should().Be(oracleRna,
                because: $"IsValidRna must equal the membership oracle; offending input: {Describe(input)}");
        }
    }

    /// <summary>
    /// RB: a fixed-seed sweep over the printable-ASCII range (0x20–0x7E) — letters,
    /// digits, punctuation — must never throw and must classify exactly per the
    /// membership oracle, AND IsValidDna must agree with TryCreate's success flag
    /// for the very same string. This pins that the total predicate and the factory
    /// never disagree about validity on random ASCII, and that TryCreate never
    /// leaks a non-ArgumentException.
    /// </summary>
    [Test]
    public void IsValid_RandomAscii_AgreesWithTryCreateAndNeverThrows()
    {
        const string dnaAlphabet = "ACGT";

        for (int trial = 0; trial < 2000; trial++)
        {
            int len = Rng.Next(0, 17);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = (char)Rng.Next(0x20, 0x7F); // printable ASCII
            string input = new string(chars);

            bool predicate = false;
            DnaSequence? result = null;
            bool created = false;
            var act = () =>
            {
                predicate = input.AsSpan().IsValidDna();
                created = DnaSequence.TryCreate(input, out result);
            };
            act.Should().NotThrow(
                $"neither the predicate nor the factory may throw on random ASCII; input: {Describe(input)}");

            bool oracle = input.All(ch => dnaAlphabet.Contains(char.ToUpperInvariant(ch)));
            predicate.Should().Be(oracle,
                because: $"IsValidDna must equal the membership oracle; input: {Describe(input)}");
            created.Should().Be(oracle,
                because: $"TryCreate success must equal validity; input: {Describe(input)}");
            (result is null).Should().Be(!oracle,
                because: $"result is non-null iff validation succeeded; input: {Describe(input)}");
        }
    }

    /// <summary>Renders a fuzz string with escaped non-printables so failures are diagnosable.</summary>
    private static string Describe(string s)
    {
        var sb = new System.Text.StringBuilder("\"");
        foreach (char c in s)
            sb.Append(c is >= ' ' and < '\u007F' ? c.ToString() : $"\\u{(int)c:X4}");
        return sb.Append('"').ToString();
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-COMPLEX-001 — sequence (linguistic) complexity
    //  Typed overload:      CalculateLinguisticComplexity(DnaSequence, int)
    //                       (only validated A/C/G/T reach it; result ∈ [0,1])
    //  Raw-string overload: CalculateLinguisticComplexity(string, int)
    //                       (lenient: never throws; DNA denominator; may exceed 1)
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-COMPLEX-001 — sequence (linguistic) complexity

    #region BE — Boundary: empty string

    /// <summary>
    /// BE: the empty string is the lower size boundary. The typed DnaSequence path
    /// (the ctor short-circuits "" to an empty sequence) and the lenient raw-string
    /// overload (IsNullOrEmpty → 0) must BOTH return exactly 0 — INV-02, no
    /// accumulation, no division by zero, no exception
    /// (Linguistic_Complexity.md §2.4 INV-02, §6.1).
    /// </summary>
    [Test]
    public void LinguisticComplexity_EmptyString_IsZeroAndDoesNotThrow()
    {
        var act = () =>
        {
            SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(string.Empty))
                .Should().Be(0.0, because: "an empty sequence has no subword vocabulary (INV-02)");
            SequenceComplexity.CalculateLinguisticComplexity(string.Empty)
                .Should().Be(0.0, because: "the lenient overload short-circuits empty input to 0");
        };

        act.Should().NotThrow("the empty string is a defined boundary input, not an error");
    }

    #endregion

    #region BE — Boundary: single character

    /// <summary>
    /// BE: a one-base sequence is the minimal non-empty input. The only word length
    /// in range is 1, with exactly one distinct 1-mer observed out of
    /// V_max,1 = min(4, 1) = 1 possible, so LC is exactly 1.0 — a positive,
    /// well-defined value (Linguistic_Complexity.md §6.1 "single nucleotide →
    /// positive value"). Verified over all four bases and a lowercase form
    /// (case-folded), with no off-by-one on the length-1 boundary.
    /// </summary>
    [TestCase('A')]
    [TestCase('C')]
    [TestCase('G')]
    [TestCase('T')]
    [TestCase('g')]
    public void LinguisticComplexity_SingleCharacter_IsPositiveAndInRange(char baseChar)
    {
        double lc = SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(baseChar.ToString()));

        lc.Should().BeInRange(0.0, 1.0,
            because: "linguistic complexity is DNA-bounded to [0,1] (INV-01)");
        lc.Should().BeApproximately(1.0, Tolerance,
            because: "one base has one distinct 1-mer out of one possible → LC = 1 (single 1-mer)");
    }

    #endregion

    #region BE — Boundary: all-same nucleotide (homopolymer → minimum complexity)

    /// <summary>
    /// BE/theory: a homopolymer (all-same nucleotide) is the MINIMUM-complexity
    /// extreme. At every word length only ONE distinct word is observed, while the
    /// denominator min(4^i, N−i+1) grows with available positions, so LC collapses
    /// toward 0 as the sequence lengthens (Linguistic_Complexity.md §6.1
    /// "homopolymer → low value"). For each tested length LC must stay finite,
    /// inside [0,1], and a longer homopolymer must be no MORE complex than a shorter
    /// one — the monotone signature of a low-complexity tract.
    /// </summary>
    [TestCase('A')]
    [TestCase('C')]
    [TestCase('G')]
    [TestCase('T')]
    public void LinguisticComplexity_Homopolymer_IsMinimalAndDecreasesWithLength(char baseChar)
    {
        double lcShort = SequenceComplexity.CalculateLinguisticComplexity(
            new DnaSequence(new string(baseChar, 10)));
        double lcLong = SequenceComplexity.CalculateLinguisticComplexity(
            new DnaSequence(new string(baseChar, 200)));

        lcShort.Should().BeInRange(0.0, 1.0, because: "LC is DNA-bounded to [0,1] (INV-01)");
        lcLong.Should().BeInRange(0.0, 1.0, because: "LC is DNA-bounded to [0,1] (INV-01)");
        lcLong.Should().BeLessThanOrEqualTo(lcShort,
            because: "a homopolymer is minimum complexity; a longer one is no more complex than a shorter one");
        lcLong.Should().BeLessThan(0.25,
            because: "a long homopolymer observes one word per length → LC near the minimum");
    }

    /// <summary>
    /// Theory ordering: a homopolymer must be ≤ a maximally diverse sequence of the
    /// SAME length. "AAAA…" (one distinct word per length) is the minimum; a
    /// permutation-rich sequence approaches the maximum. This pins the core
    /// complexity contract — homopolymer ≤ diverse — that fuzzing must never let
    /// invert (Linguistic_Complexity.md §6.1).
    /// </summary>
    [Test]
    public void LinguisticComplexity_Homopolymer_IsNoGreaterThanDiverse()
    {
        const int length = 256;
        double homopolymer = SequenceComplexity.CalculateLinguisticComplexity(
            new DnaSequence(new string('A', length)));
        double diverse = SequenceComplexity.CalculateLinguisticComplexity(
            new DnaSequence(RandomDna(length)));

        homopolymer.Should().BeLessThanOrEqualTo(diverse,
            because: "a homopolymer is minimum complexity; a diverse sequence of equal length is at least as complex");
    }

    #endregion

    #region BE — Boundary: extremely long (no hang, stays in range)

    /// <summary>
    /// BE/OVF: an extremely long valid sequence (200,000 bases) must compute without
    /// hang, overflow, or precision blow-up, and the result must stay finite and in
    /// the closed DNA range [0, 1] (INV-01). The accumulators are <c>long</c>; this
    /// guards that the observed/possible summation and the final ratio do not
    /// overflow or drift out of range at scale. Both a fixed-seed random input and a
    /// long homopolymer (the low-complexity extreme) are exercised.
    /// </summary>
    [Test]
    public void LinguisticComplexity_ExtremelyLong_StaysInRangeAndDoesNotHang()
    {
        const int length = 200_000;

        double random = SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(RandomDna(length)));
        double.IsFinite(random).Should().BeTrue(
            because: "LC must be a finite number even at scale");
        random.Should().BeInRange(0.0, 1.0,
            because: "LC is DNA-bounded to [0,1]; it cannot escape the range at scale (INV-01)");

        double homopolymer = SequenceComplexity.CalculateLinguisticComplexity(
            new DnaSequence(new string('G', length)));
        homopolymer.Should().BeInRange(0.0, 1.0,
            because: "even a huge homopolymer stays in [0,1] without overflow");
        homopolymer.Should().BeLessThanOrEqualTo(random,
            because: "a long homopolymer is minimum complexity; it is no more complex than a long random " +
                     "sequence (at this scale both LC values can collapse to the same near-zero ratio, so " +
                     "the theory-correct relation is ≤, matching the order-robust RB sweep assertion)");
    }

    #endregion

    #region BE — typed overload: validation gate (null / maxWordLength)

    /// <summary>
    /// BE: the typed overload's documented validation gate
    /// (Linguistic_Complexity.md §3.3, §6.1): a null DnaSequence is an
    /// ArgumentNullException (explicit guard, never a NullReferenceException), and
    /// maxWordLength &lt; 1 is an ArgumentOutOfRangeException — both *intentional*
    /// validation exceptions, not crashes.
    /// </summary>
    [Test]
    public void LinguisticComplexity_Typed_NullSequence_ThrowsArgumentNullException()
    {
        var act = () => SequenceComplexity.CalculateLinguisticComplexity((DnaSequence)null!);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload guards null explicitly; this is documented validation, not a crash");
    }

    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void LinguisticComplexity_Typed_MaxWordLengthBelowOne_ThrowsArgumentOutOfRange(int maxWordLength)
    {
        var act = () => SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence("ACGTACGT"), maxWordLength);

        act.Should().Throw<ArgumentOutOfRangeException>(
            "maxWordLength < 1 is rejected by the documented validation gate, not a silent miscompute");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Lenient raw-string overload: never throws; DNA denominator; may exceed 1
    // ───────────────────────────────────────────────────────────────────

    #region INJ/BE — lenient raw-string: null / arbitrary chars never throw

    /// <summary>
    /// BE: the lenient raw-string overload returns 0 for null/empty input
    /// (IsNullOrEmpty short-circuit) — no NullReferenceException, no division, no
    /// hang (Linguistic_Complexity.md §3.3). maxWordLength &lt; 1 is NOT validated
    /// here, but the min(maxWordLength, N) loop cap means the inner loop never runs
    /// and the result is a defined 0 — a graceful no-op, not a crash.
    /// </summary>
    [Test]
    public void LinguisticComplexity_RawString_NullEmptyAndZeroWordLength_AreZeroAndDoNotThrow()
    {
        var act = () =>
        {
            SequenceComplexity.CalculateLinguisticComplexity((string)null!)
                .Should().Be(0.0, because: "null is short-circuited to 0, never dereferenced");
            SequenceComplexity.CalculateLinguisticComplexity(string.Empty)
                .Should().Be(0.0, because: "empty input is short-circuited to 0");
            SequenceComplexity.CalculateLinguisticComplexity("ACGTACGT", 0)
                .Should().Be(0.0, because: "the min(maxWordLength,N) cap makes maxWordLength=0 a defined 0, not a throw");
        };

        act.Should().NotThrow("the lenient raw-string overload is total and never throws on these inputs");
    }

    /// <summary>
    /// INJ: the lenient raw-string overload does NOT validate the alphabet, so it
    /// must NEVER throw on arbitrary characters — digits, gaps, whitespace, the
    /// ambiguity code N, the RNA base U, an embedded null byte, unicode letters,
    /// full-width look-alikes, an astral surrogate pair, even a lone surrogate half.
    /// The result must always be a FINITE, non-negative number. Because the
    /// denominator stays the DNA 4^i, a value &gt; 1 is a *documented* consequence
    /// for non-ACGT input (Linguistic_Complexity.md §5.3, §6.1), not a defect — so
    /// we assert finiteness and ≥ 0, never an unhandled crash.
    /// </summary>
    [TestCase("ACGTN", TestName = "LinguisticComplexity_RawString_AmbiguityN_NeverThrows")]
    [TestCase("ACGU", TestName = "LinguisticComplexity_RawString_RnaBaseU_NeverThrows")]
    [TestCase("ACGT1234", TestName = "LinguisticComplexity_RawString_Digits_NeverThrows")]
    [TestCase("AC GT-AC", TestName = "LinguisticComplexity_RawString_WhitespaceAndGap_NeverThrows")]
    [TestCase("ACG\0TACG", TestName = "LinguisticComplexity_RawString_EmbeddedNullByte_NeverThrows")]
    [TestCase("ACGTαβγ", TestName = "LinguisticComplexity_RawString_GreekLetters_NeverThrows")]
    [TestCase("ＡＣＧＴ", TestName = "LinguisticComplexity_RawString_FullWidthLatin_NeverThrows")]
    [TestCase("ACG😀TACG", TestName = "LinguisticComplexity_RawString_AstralSurrogatePair_NeverThrows")]
    [TestCase("ACG\uD83DTAC", TestName = "LinguisticComplexity_RawString_LoneHighSurrogate_NeverThrows")]
    public void LinguisticComplexity_RawString_ArbitraryChars_AreFiniteAndNeverThrow(string input)
    {
        double lc = double.NaN;
        var act = () => lc = SequenceComplexity.CalculateLinguisticComplexity(input);

        act.Should().NotThrow(
            "the lenient raw-string overload does not validate the alphabet and must never crash on garbage");
        double.IsFinite(lc).Should().BeTrue(
            because: "LC must be a finite number even on non-ACGT input");
        lc.Should().BeGreaterThanOrEqualTo(0.0,
            because: "the LC ratio is non-negative; non-ACGT may exceed 1 (DNA denominator) but never goes negative");
    }

    #endregion

    #region RB — Random-byte sweeps (never throw; finite; homopolymer ≤ diverse)

    /// <summary>
    /// RB: fixed-seed random-byte sweeps over the strict typed path. Every random
    /// VALID DNA string of varied length must yield a finite LC inside [0, 1]
    /// (INV-01) without throwing, and a homopolymer of the SAME length must be no
    /// more complex than the random one — the theory-correct ordering must hold on
    /// every draw, never inverted.
    /// </summary>
    [Test]
    public void LinguisticComplexity_RandomValidDna_IsFiniteInRangeAndDominatesHomopolymer()
    {
        for (int trial = 0; trial < 200; trial++)
        {
            int length = Rng.Next(1, 300);
            string randomSeq = RandomDna(length);

            double lc = double.NaN;
            var act = () => lc = SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(randomSeq));
            act.Should().NotThrow($"valid DNA must never throw; length={length}");

            double.IsFinite(lc).Should().BeTrue(
                because: $"LC must be finite; length={length}");
            lc.Should().BeInRange(0.0, 1.0,
                because: $"LC is DNA-bounded to [0,1] (INV-01); length={length}");

            double homopolymer = SequenceComplexity.CalculateLinguisticComplexity(
                new DnaSequence(new string(randomSeq[0], length)));
            homopolymer.Should().BeLessThanOrEqualTo(lc + Tolerance,
                because: $"a homopolymer is minimum complexity; it cannot exceed the random sequence; length={length}");
        }
    }

    /// <summary>
    /// RB: fixed-seed random-byte sweeps over the LENIENT raw-string path with
    /// arbitrary BMP code points (control chars, null bytes, lone surrogate halves —
    /// pure garbage). The overload must NEVER throw — no IndexOutOfRange/encoding
    /// surprise from surrogate halves — and must always return a finite, non-negative
    /// number. This pins that the strict/lenient boundary holds: the same random
    /// bytes the typed path would reject at construction, the lenient overload
    /// carries through to a defined finite score without crashing.
    /// </summary>
    [Test]
    public void LinguisticComplexity_RawString_RandomBmpBytes_AreFiniteAndNeverThrow()
    {
        for (int trial = 0; trial < 200; trial++)
        {
            int length = Rng.Next(0, 64);
            string garbage = RandomBmpChars(length);

            double lc = double.NaN;
            var act = () => lc = SequenceComplexity.CalculateLinguisticComplexity(garbage);
            act.Should().NotThrow(
                $"the lenient overload must never throw on random bytes; input: {Describe(garbage)}");

            double.IsFinite(lc).Should().BeTrue(
                because: $"LC must be finite on random bytes; input: {Describe(garbage)}");
            lc.Should().BeGreaterThanOrEqualTo(0.0,
                because: $"the LC ratio is non-negative on any input; input: {Describe(garbage)}");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-ENTROPY-001 — Shannon entropy of base composition
    //  Typed overload:      CalculateShannonEntropy(DnaSequence)
    //                       (only validated A/C/G/T reach it; throws on null)
    //  Raw-string overload: CalculateShannonEntropy(string)
    //                       (lenient: never throws; counts only A/T/G/C)
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-ENTROPY-001 — Shannon entropy

    #region BE — Boundary: empty string

    /// <summary>
    /// BE: the empty string is the lower size boundary. Both surfaces define it as
    /// entropy 0 (INV-02): the typed path goes through the ctor's IsNullOrEmpty
    /// short-circuit to an empty sequence whose core returns 0; the raw-string
    /// overload short-circuits empty input to 0 directly. No frequency computation,
    /// no division by zero, no log(0), no exception
    /// (Shannon_Entropy.md §2.4 INV-02, §6.1 "Empty sequence → 0.0").
    /// </summary>
    [Test]
    public void ShannonEntropy_EmptyString_IsZeroAndDoesNotThrow()
    {
        var act = () =>
        {
            SequenceComplexity.CalculateShannonEntropy(new DnaSequence(string.Empty))
                .Should().Be(0.0, because: "an empty sequence has no information content (INV-02)");
            SequenceComplexity.CalculateShannonEntropy(string.Empty)
                .Should().Be(0.0, because: "the lenient overload short-circuits empty input to 0");
        };

        act.Should().NotThrow("the empty string is a defined boundary input, not an error");
    }

    #endregion

    #region BE — Boundary: single character

    /// <summary>
    /// BE: a one-base sequence is the minimal non-empty input. A single symbol has
    /// probability 1.0, so its only term is −1·log2(1) = 0 — the MINIMUM entropy.
    /// This is the p=1 edge of the −p·log2(p) sum and pins that a length-1 input
    /// neither off-by-ones nor leaks a NaN. Verified over all four bases and a
    /// lowercase form (case-folded by the typed ctor / the raw-string overload).
    /// </summary>
    [TestCase('A')]
    [TestCase('C')]
    [TestCase('G')]
    [TestCase('T')]
    [TestCase('g')]
    public void ShannonEntropy_SingleCharacter_IsZero(char baseChar)
    {
        double typed = SequenceComplexity.CalculateShannonEntropy(new DnaSequence(baseChar.ToString()));
        double raw = SequenceComplexity.CalculateShannonEntropy(baseChar.ToString());

        typed.Should().BeApproximately(0.0, Tolerance,
            because: "one symbol has probability 1; −1·log2(1) = 0 (minimum entropy)");
        raw.Should().BeApproximately(0.0, Tolerance,
            because: "the raw-string overload counts the single A/T/G/C base identically");
    }

    #endregion

    #region BE — Boundary: all-same nucleotide (homopolymer → entropy 0)

    /// <summary>
    /// BE/theory: a homopolymer (all-same nucleotide) is the MINIMUM-entropy extreme
    /// (INV-03). One symbol has probability 1.0 and every other base has probability
    /// 0; only the p=1 term contributes and it is −1·log2(1) = 0, while the p=0 bases
    /// are SKIPPED (count>0 guard) so no −0·log2(0) = NaN can leak. Entropy must be
    /// EXACTLY 0 at every length, not merely "small"
    /// (Shannon_Entropy.md §2.4 INV-03, §6.1 "Homopolymer AAAA → 0.0").
    /// </summary>
    [TestCase('A', 1)]
    [TestCase('C', 2)]
    [TestCase('G', 50)]
    [TestCase('T', 1000)]
    public void ShannonEntropy_Homopolymer_IsExactlyZeroAndFinite(char baseChar, int length)
    {
        double entropy = SequenceComplexity.CalculateShannonEntropy(new DnaSequence(new string(baseChar, length)));

        double.IsNaN(entropy).Should().BeFalse(
            because: "the zero-probability bases are skipped, so no −0·log2(0) NaN can leak");
        entropy.Should().Be(0.0,
            because: "a homopolymer has one symbol at probability 1 → entropy is exactly 0 (INV-03)");
    }

    #endregion

    #region Theory — uniform distribution → maximum entropy log2(k)

    /// <summary>
    /// Theory (the MAXIMUM-entropy contract): a sequence using k of the four DNA
    /// bases in EQUAL proportion has entropy exactly log2(k) — the maximum for that
    /// support (Shannon_Entropy.md §2.2, §4.2 reference table). k=4 uniform "ACGT"
    /// → 2.0 bits (the DNA maximum, INV-01 upper bound); k=2 "AT"/"GC" 50/50 → 1.0
    /// bit; k=1 → 0.0. This pins the upper edge of the contract so it can never
    /// drift, and confirms log base 2 (bits), not natural or base-10.
    /// </summary>
    [TestCase("ACGT", 2.0, TestName = "ShannonEntropy_UniformAcgt_Is2Bits")]
    [TestCase("ACGTACGT", 2.0, TestName = "ShannonEntropy_UniformAcgtRepeated_Is2Bits")]
    [TestCase("AT", 1.0, TestName = "ShannonEntropy_HalfHalfAt_Is1Bit")]
    [TestCase("GC", 1.0, TestName = "ShannonEntropy_HalfHalfGc_Is1Bit")]
    [TestCase("AATT", 1.0, TestName = "ShannonEntropy_HalfHalfAatt_Is1Bit")]
    [TestCase("AAAACCGT", 1.75, TestName = "ShannonEntropy_DyadicComposition_Is1Point75Bits")]
    public void ShannonEntropy_UniformOverKSymbols_IsLog2K(string sequence, double expected)
    {
        double typed = SequenceComplexity.CalculateShannonEntropy(new DnaSequence(sequence));
        double raw = SequenceComplexity.CalculateShannonEntropy(sequence);

        typed.Should().BeApproximately(expected, Tolerance,
            because: $"uniform over k symbols gives H = log2(k); '{sequence}' → {expected} bits");
        raw.Should().BeApproximately(expected, Tolerance,
            because: "the raw-string overload computes the identical base-2 entropy");
    }

    #endregion

    #region BE — Boundary: extremely long (no hang, finite, in [0,2])

    /// <summary>
    /// BE/OVF: an extremely long valid sequence (1,000,000 bases) must compute
    /// without hang, overflow, or precision blow-up, and the result must stay finite
    /// and inside the closed DNA range [0, 2] (INV-01). A known-composition input
    /// ("ACGT" tiled = exactly uniform) pins the exact 2.0 maximum at scale; a
    /// fixed-seed random input pins finiteness and range. The counts are int; this
    /// guards that the (double)count/total ratio and the −p·log2(p) sum neither
    /// overflow nor drift out of range when N is large.
    /// </summary>
    [Test]
    public void ShannonEntropy_ExtremelyLong_StaysInRangeAndDoesNotHang()
    {
        const int length = 1_000_000;

        // Exactly uniform composition: "ACGT" tiled is 25% of each base → 2.0 bits.
        var uniform = new DnaSequence(string.Concat(Enumerable.Repeat("ACGT", length / 4)));
        uniform.Length.Should().Be(length);
        SequenceComplexity.CalculateShannonEntropy(uniform).Should().BeApproximately(2.0, Tolerance,
            because: "a sequence that is exactly 25% of each base has entropy 2 bits, at any length");

        // Fixed-seed random long sequence: result must be finite and in [0, 2].
        double random = SequenceComplexity.CalculateShannonEntropy(new DnaSequence(RandomDna(length)));
        double.IsFinite(random).Should().BeTrue(because: "entropy must be a finite number even at scale");
        random.Should().BeInRange(0.0, 2.0,
            because: "DNA Shannon entropy is bounded to [0, 2] bits; it cannot escape at scale (INV-01)");
    }

    #endregion

    #region INJ — typed overload: null sequence throws ArgumentNullException

    /// <summary>
    /// INJ/BE: a null DnaSequence is the boundary of "no input" on the typed path.
    /// It is guarded explicitly with ArgumentNullException.ThrowIfNull — a
    /// *documented, intentional* validation exception, NOT a raw
    /// NullReferenceException (Shannon_Entropy.md §3.3). Pins that null is rejected
    /// at the gate rather than dereferenced.
    /// </summary>
    [Test]
    public void ShannonEntropy_Typed_NullSequence_ThrowsArgumentNullException()
    {
        var act = () => SequenceComplexity.CalculateShannonEntropy((DnaSequence)null!);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload guards null explicitly; this is documented validation, not a crash");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Lenient raw-string overload: never throws; counts only A/T/G/C;
    //  non-ACGT symbols are ignored (no NaN/Inf, no crash on garbage)
    // ───────────────────────────────────────────────────────────────────

    #region INJ/BE — lenient raw-string: null / non-nucleotide chars never throw

    /// <summary>
    /// BE: the lenient raw-string overload returns 0 for null/empty input
    /// (IsNullOrEmpty short-circuit) — no NullReferenceException, no division, no
    /// hang (Shannon_Entropy.md §3.3, §5.2). Pins that "no input" is a defined 0,
    /// not a crash, on the surface documented to accept arbitrary text.
    /// </summary>
    [Test]
    public void ShannonEntropy_RawString_NullAndEmpty_AreZeroAndDoNotThrow()
    {
        var act = () =>
        {
            SequenceComplexity.CalculateShannonEntropy((string)null!)
                .Should().Be(0.0, because: "null is short-circuited to 0, never dereferenced");
            SequenceComplexity.CalculateShannonEntropy(string.Empty)
                .Should().Be(0.0, because: "empty input is short-circuited to 0");
        };

        act.Should().NotThrow("the lenient raw-string overload is total and never throws on null/empty");
    }

    /// <summary>
    /// INJ: the lenient raw-string overload ignores every non-A/T/G/C symbol
    /// (Shannon_Entropy.md §5.2, §5.3: "counts only A/T/G/C and ignores
    /// non-standard bases such as N or other ambiguity codes"). So injected garbage
    /// — ambiguity code N, the RNA base U, digits, gaps, whitespace, an embedded
    /// null byte, unicode letters — does NOT change the entropy of the A/T/G/C bases
    /// that ARE present, and never throws. Each case interleaves garbage into a
    /// uniform ACGT core, which must therefore still read as exactly 2.0 bits.
    /// </summary>
    [TestCase("ACGTN", TestName = "ShannonEntropy_RawString_AmbiguityN_Ignored_Is2Bits")]
    [TestCase("ACGTU", TestName = "ShannonEntropy_RawString_RnaBaseU_Ignored_Is2Bits")]
    [TestCase("A1C2G3T4", TestName = "ShannonEntropy_RawString_Digits_Ignored_Is2Bits")]
    [TestCase("A-C-G-T", TestName = "ShannonEntropy_RawString_Gaps_Ignored_Is2Bits")]
    [TestCase("A C G T", TestName = "ShannonEntropy_RawString_Whitespace_Ignored_Is2Bits")]
    [TestCase("A\0C\0G\0T", TestName = "ShannonEntropy_RawString_NullBytes_Ignored_Is2Bits")]
    [TestCase("AαCβGγT", TestName = "ShannonEntropy_RawString_GreekLetters_Ignored_Is2Bits")]
    public void ShannonEntropy_RawString_NonNucleotideChars_AreIgnoredAndNeverThrow(string input)
    {
        double entropy = double.NaN;
        var act = () => entropy = SequenceComplexity.CalculateShannonEntropy(input);

        act.Should().NotThrow(
            "the lenient raw-string overload does not validate the alphabet and must never crash on garbage");
        double.IsFinite(entropy).Should().BeTrue(
            because: "non-ACGT symbols are ignored, so no NaN/Inf can leak from the entropy sum");
        entropy.Should().BeApproximately(2.0, Tolerance,
            because: "the four A/T/G/C bases are uniform; ignored non-nucleotide chars do not shift their entropy");
    }

    /// <summary>
    /// INJ/BE: a raw string containing NO A/T/G/C base at all — pure non-nucleotide
    /// garbage — must read as exactly 0, not NaN and not a crash. The core's
    /// total==0 guard prevents the 0/0 division that would otherwise produce NaN
    /// (Shannon_Entropy.md §5.2 "returns 0.0 if no counted DNA bases are present
    /// after filtering"). This is the critical log(0)/divide-by-zero boundary.
    /// </summary>
    [TestCase("N", TestName = "ShannonEntropy_RawString_OnlyAmbiguityN_IsZero")]
    [TestCase("NNNN", TestName = "ShannonEntropy_RawString_OnlyAmbiguityNNNN_IsZero")]
    [TestCase("12345", TestName = "ShannonEntropy_RawString_OnlyDigits_IsZero")]
    [TestCase("------", TestName = "ShannonEntropy_RawString_OnlyGaps_IsZero")]
    [TestCase("   ", TestName = "ShannonEntropy_RawString_OnlyWhitespace_IsZero")]
    [TestCase("\0\0\0", TestName = "ShannonEntropy_RawString_OnlyNullBytes_IsZero")]
    [TestCase("αβγδ", TestName = "ShannonEntropy_RawString_OnlyUnicode_IsZero")]
    [TestCase("😀😀", TestName = "ShannonEntropy_RawString_OnlyAstral_IsZero")]
    public void ShannonEntropy_RawString_NoNucleotidePresent_IsZeroNotNaN(string input)
    {
        double entropy = double.NaN;
        var act = () => entropy = SequenceComplexity.CalculateShannonEntropy(input);

        act.Should().NotThrow(
            "input with no A/T/G/C must not crash, even from surrogate halves of an astral pair");
        double.IsNaN(entropy).Should().BeFalse(
            because: "the total==0 guard prevents the 0/0 division that would otherwise leak NaN");
        entropy.Should().Be(0.0,
            because: "no counted DNA base is present after filtering → entropy is the defined 0 (§5.2)");
    }

    #endregion

    #region RB — Random-byte sweeps (never throw; finite; in [0,2]; all-same→0; uniform→log2 k)

    /// <summary>
    /// RB: fixed-seed random-byte sweeps over the strict typed path. Every random
    /// VALID DNA string of varied length must yield a finite entropy inside [0, 2]
    /// bits (INV-01) without throwing, and a homopolymer of the SAME first base must
    /// be exactly 0 — the theory-correct minimum that fuzzing must never let drift
    /// (Shannon_Entropy.md §2.4 INV-01/INV-03).
    /// </summary>
    [Test]
    public void ShannonEntropy_RandomValidDna_IsFiniteInRangeAndHomopolymerIsZero()
    {
        for (int trial = 0; trial < 300; trial++)
        {
            int length = Rng.Next(1, 300);
            string randomSeq = RandomDna(length);

            double entropy = double.NaN;
            var act = () => entropy = SequenceComplexity.CalculateShannonEntropy(new DnaSequence(randomSeq));
            act.Should().NotThrow($"valid DNA must never throw; length={length}");

            double.IsFinite(entropy).Should().BeTrue(because: $"entropy must be finite; length={length}");
            entropy.Should().BeInRange(0.0, 2.0,
                because: $"DNA entropy is bounded to [0,2] bits (INV-01); length={length}");

            double homopolymer = SequenceComplexity.CalculateShannonEntropy(
                new DnaSequence(new string(randomSeq[0], length)));
            homopolymer.Should().Be(0.0,
                because: $"a homopolymer has one symbol at probability 1 → entropy exactly 0 (INV-03); length={length}");
        }
    }

    /// <summary>
    /// RB/theory: a fixed-seed sweep of EXACTLY-uniform sequences over a random
    /// support of k ∈ {1,2,3,4} distinct bases. With each chosen base repeated an
    /// equal number of times, the entropy must equal log2(k) to within tolerance —
    /// the maximum-entropy contract across the whole support spectrum, asserted on
    /// every draw so neither the log base nor the probability normalization can
    /// drift (Shannon_Entropy.md §2.2).
    /// </summary>
    [Test]
    public void ShannonEntropy_RandomUniformSupport_EqualsLog2OfSupportSize()
    {
        const string bases = "ACGT";

        for (int trial = 0; trial < 200; trial++)
        {
            int k = Rng.Next(1, 5);                 // 1..4 distinct bases
            int perBase = Rng.Next(1, 20);          // equal count of each → uniform
            var chosen = new System.Collections.Generic.HashSet<char>();
            while (chosen.Count < k)
                chosen.Add(bases[Rng.Next(bases.Length)]);

            var sb = new System.Text.StringBuilder();
            foreach (char b in chosen)
                sb.Append(b, perBase);
            string uniform = sb.ToString();

            double entropy = SequenceComplexity.CalculateShannonEntropy(new DnaSequence(uniform));

            entropy.Should().BeApproximately(Math.Log2(k), Tolerance,
                because: $"a uniform distribution over k={k} symbols has entropy log2(k)={Math.Log2(k):F4} bits");
        }
    }

    /// <summary>
    /// RB: fixed-seed random-byte sweeps over the LENIENT raw-string path with
    /// arbitrary BMP code points (control chars, null bytes, lone surrogate halves —
    /// pure garbage). The overload must NEVER throw — no IndexOutOfRange/encoding
    /// surprise from surrogate halves — and must always return a FINITE, non-negative
    /// entropy in [0, 2] (only A/T/G/C are ever counted, so the support is at most 4
    /// and the value cannot exceed log2(4) = 2). Pins the strict/lenient boundary:
    /// the same random bytes the typed path would reject at construction, the lenient
    /// overload carries through to a defined finite entropy without crashing or
    /// leaking NaN.
    /// </summary>
    [Test]
    public void ShannonEntropy_RawString_RandomBmpBytes_AreFiniteInRangeAndNeverThrow()
    {
        for (int trial = 0; trial < 300; trial++)
        {
            int length = Rng.Next(0, 64);
            string garbage = RandomBmpChars(length);

            double entropy = double.NaN;
            var act = () => entropy = SequenceComplexity.CalculateShannonEntropy(garbage);
            act.Should().NotThrow(
                $"the lenient overload must never throw on random bytes; input: {Describe(garbage)}");

            double.IsFinite(entropy).Should().BeTrue(
                because: $"only A/T/G/C are counted, so no NaN/Inf can leak; input: {Describe(garbage)}");
            entropy.Should().BeInRange(0.0, 2.0,
                because: $"entropy over a ≤4-symbol DNA support is bounded to [0,2]; input: {Describe(garbage)}");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-GCSKEW-001 — GC skew : (G − C) / (G + C), range [−1, 1]
    //  Typed overload:  CalculateGcSkew(DnaSequence) — strict, ArgumentNull on null
    //  Raw-string:      CalculateGcSkew(string)      — lenient, counts only G/C
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-GCSKEW-001 — GC skew

    #region BE — Boundary: empty string

    /// <summary>
    /// BE: the empty string is the lower size boundary AND the canonical
    /// DivideByZero boundary (G + C = 0). The typed path defines the empty
    /// DnaSequence (DnaSequence.cs lines 24–28) and CalculateGcSkew must return the
    /// defined 0 via the zero-denominator guard (GcSkewCalculator.cs lines 44) —
    /// no DivideByZeroException, no NaN, no crash. The lenient raw-string overload
    /// returns 0 for "" and null by its IsNullOrEmpty short-circuit (lines 32–33).
    /// </summary>
    [Test]
    public void GcSkew_EmptyAndNull_IsZeroAndDoesNotThrow()
    {
        var typed = () =>
        {
            double skew = GcSkewCalculator.CalculateGcSkew(new DnaSequence(string.Empty));
            skew.Should().Be(0.0,
                because: "an empty sequence has no G or C; skew is defined as 0 by the zero-denominator guard");
        };
        typed.Should().NotThrow("the empty sequence is a defined boundary, not a DivideByZero error");

        var rawEmpty = () => GcSkewCalculator.CalculateGcSkew(string.Empty).Should().Be(0.0);
        var rawNull = () => GcSkewCalculator.CalculateGcSkew((string)null!).Should().Be(0.0);
        rawEmpty.Should().NotThrow("the lenient overload short-circuits empty input to 0");
        rawNull.Should().NotThrow<NullReferenceException>(
            "null is handled by the IsNullOrEmpty gate, never dereferenced");
        rawNull.Should().NotThrow("null is a defined 'no input' boundary returning 0, not an error");
    }

    /// <summary>
    /// BE: a null DnaSequence is the explicit-guard boundary on the strict typed
    /// path. CalculateGcSkew(DnaSequence) must throw the *documented, intentional*
    /// ArgumentNullException (GcSkewCalculator.cs line 23), never a raw
    /// NullReferenceException from dereferencing <c>.Sequence</c>.
    /// </summary>
    [Test]
    public void GcSkew_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => GcSkewCalculator.CalculateGcSkew((DnaSequence)null!);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload guards null with the documented ArgumentNullException, not a crash");
    }

    #endregion

    #region BE — Boundary: single base

    /// <summary>
    /// BE: a one-base sequence is the minimal non-empty input and the binary
    /// extreme of the skew. A single 'G' is the maximum +1 (G=1, C=0 → 1/1);
    /// a single 'C' is the minimum −1; a single 'A' or 'T' has NO G/C, so the
    /// G+C=0 guard returns the defined 0 (NOT NaN). Lowercase is accepted on the
    /// typed path (case-folded at construction) and yields the same value. Verified
    /// over all four bases in both cases on the strict typed overload.
    /// </summary>
    [TestCase('G', +1.0)]
    [TestCase('C', -1.0)]
    [TestCase('A', 0.0)]
    [TestCase('T', 0.0)]
    [TestCase('g', +1.0)]
    [TestCase('c', -1.0)]
    [TestCase('a', 0.0)]
    [TestCase('t', 0.0)]
    public void GcSkew_SingleBase_IsBinaryExtremeOrZero(char baseChar, double expected)
    {
        double skew = GcSkewCalculator.CalculateGcSkew(new DnaSequence(baseChar.ToString()));

        skew.Should().BeApproximately(expected, Tolerance,
            because: expected switch
            {
                +1.0 => $"a single '{baseChar}' is all-G → maximum skew +1",
                -1.0 => $"a single '{baseChar}' is all-C → minimum skew −1",
                _ => $"a single '{baseChar}' has no G/C, so the G+C=0 guard returns the defined 0, not NaN"
            });
    }

    #endregion

    #region BE — Boundary: no G or C (the DivideByZero boundary)

    /// <summary>
    /// BE: the critical G + C = 0 boundary. A sequence with NO guanine and NO
    /// cytosine (all-A, all-T, AT-only) makes the denominator zero. The typed path
    /// over valid A/T-only DNA, and the lenient raw-string path over A/T plus
    /// ignored garbage (N, digits, gap, null byte, unicode), must BOTH return the
    /// defined 0 — never DivideByZeroException, never NaN. This is the headline
    /// fuzz target for SEQ-GCSKEW-001 (GC_Skew.md INV-02 / §6.1).
    /// </summary>
    [TestCase("A", TestName = "GcSkew_NoGorC_SingleA_IsZero")]
    [TestCase("AAAA", TestName = "GcSkew_NoGorC_AllA_IsZero")]
    [TestCase("TTTT", TestName = "GcSkew_NoGorC_AllT_IsZero")]
    [TestCase("ATATAT", TestName = "GcSkew_NoGorC_AtAlternating_IsZero")]
    public void GcSkew_TypedNoGorC_IsZeroAndDoesNotThrow(string atOnly)
    {
        double skew = double.NaN;
        var act = () => skew = GcSkewCalculator.CalculateGcSkew(new DnaSequence(atOnly));

        act.Should().NotThrow("no G/C means G+C=0, the guarded boundary — not a DivideByZero crash");
        double.IsNaN(skew).Should().BeFalse("the zero-denominator guard returns 0, never a 0/0 NaN");
        skew.Should().Be(0.0, because: "with no G and no C the skew is defined as 0");
    }

    /// <summary>
    /// BE: the lenient raw-string path at the G+C=0 boundary. Input that contains
    /// NO G and NO C after upper-casing — pure A/T, ambiguity codes, digits, gaps,
    /// an embedded null byte, or unicode — counts no G/C, so the total==0 guard
    /// returns the defined 0 (GcSkewCalculator.cs line 44). The overload must NEVER
    /// throw and never leak NaN, even on garbage that the typed path would reject.
    /// </summary>
    [TestCase("AAAA", TestName = "GcSkew_RawNoGorC_AllA_IsZero")]
    [TestCase("ATAT", TestName = "GcSkew_RawNoGorC_AtAlternating_IsZero")]
    [TestCase("NNNN", TestName = "GcSkew_RawNoGorC_OnlyAmbiguity_IsZero")]
    [TestCase("12345", TestName = "GcSkew_RawNoGorC_OnlyDigits_IsZero")]
    [TestCase("----", TestName = "GcSkew_RawNoGorC_OnlyGaps_IsZero")]
    [TestCase("a\0t", TestName = "GcSkew_RawNoGorC_NullByteBetweenAt_IsZero")]
    [TestCase("ααα", TestName = "GcSkew_RawNoGorC_Unicode_IsZero")]
    public void GcSkew_RawStringNoGorC_IsZeroAndDoesNotThrow(string input)
    {
        double skew = double.NaN;
        var act = () => skew = GcSkewCalculator.CalculateGcSkew(input);

        act.Should().NotThrow("the lenient overload counts only G/C and never throws on non-G/C garbage");
        double.IsNaN(skew).Should().BeFalse("the total==0 guard returns 0, never a 0/0 NaN");
        skew.Should().Be(0.0, because: "no G and no C is present, so the skew is the defined 0");
    }

    #endregion

    #region BE — Boundary: alternating GC and the all-G / all-C extremes

    /// <summary>
    /// BE: the symmetry boundary. A sequence with EQUAL G and C counts — an
    /// alternating GCGC… (and its CGCG… phase) — has numerator G−C = 0, so the skew
    /// is exactly 0 regardless of length. The all-G run is the maximum +1, the
    /// all-C run is the minimum −1. Together these pin the full closed range
    /// [−1, 1] and the zero-crossing at G = C. Verified on the strict typed path
    /// (valid G/C-only DNA) across lengths.
    /// </summary>
    [TestCase("GC", 0.0, TestName = "GcSkew_Alternating_GC_IsZero")]
    [TestCase("CG", 0.0, TestName = "GcSkew_Alternating_CG_IsZero")]
    [TestCase("GCGCGCGC", 0.0, TestName = "GcSkew_Alternating_GCGCGCGC_IsZero")]
    [TestCase("GGGG", +1.0, TestName = "GcSkew_AllG_IsPlusOne")]
    [TestCase("CCCC", -1.0, TestName = "GcSkew_AllC_IsMinusOne")]
    [TestCase("GGGC", +0.5, TestName = "GcSkew_ThreeGOneC_IsHalf")]
    [TestCase("GCCC", -0.5, TestName = "GcSkew_OneGThreeC_IsMinusHalf")]
    public void GcSkew_GcComposition_MatchesFormulaAcrossRange(string input, double expected)
    {
        double skew = GcSkewCalculator.CalculateGcSkew(new DnaSequence(input));

        skew.Should().BeApproximately(expected, Tolerance,
            because: $"(G−C)/(G+C) for \"{input}\" is {expected}; the skew stays in the closed range [−1, 1]");
    }

    #endregion

    #region BE — Boundary: extremely long

    /// <summary>
    /// BE/OVF: an extremely long valid sequence (1,000,000 bases) must compute the
    /// skew without overflow, hang, or precision blow-up, and the result must stay
    /// in the closed range [−1, 1] and be finite (no NaN/Inf). A known-composition
    /// long input pins the exact value: "AG" repeated has G but no C, so the skew is
    /// exactly +1 at any length; "GC" repeated has equal G and C, so the skew is
    /// exactly 0. A fixed-seed random long sequence pins the range invariant at scale.
    /// </summary>
    [Test]
    public void GcSkew_ExtremelyLong_StaysInRangeAndDoesNotHang()
    {
        const int length = 1_000_000;

        // "AG" repeated: G present, C absent → skew = (G−0)/(G+0) = +1, at any length.
        var allGskew = new DnaSequence(string.Concat(Enumerable.Repeat("AG", length / 2)));
        allGskew.Length.Should().Be(length);
        GcSkewCalculator.CalculateGcSkew(allGskew).Should().BeApproximately(+1.0, Tolerance,
            because: "a long sequence with G but no C has skew +1, regardless of length");

        // "GC" repeated: equal G and C → skew = 0, at any length.
        var balanced = new DnaSequence(string.Concat(Enumerable.Repeat("GC", length / 2)));
        GcSkewCalculator.CalculateGcSkew(balanced).Should().BeApproximately(0.0, Tolerance,
            because: "a long sequence with equal G and C has skew 0, regardless of length");

        // Fixed-seed random long sequence: a locally-seeded RNG so other fixtures are
        // not perturbed by the shared static Rng. The skew must be a finite value in
        // the closed range [−1, 1].
        var local = new Random(7351);
        var chars = new char[length];
        const string bases = "ACGT";
        for (int i = 0; i < length; i++)
            chars[i] = bases[local.Next(bases.Length)];

        double randomSkew = GcSkewCalculator.CalculateGcSkew(new DnaSequence(new string(chars)));
        double.IsFinite(randomSkew).Should().BeTrue("the skew can never be NaN/Inf, even at scale");
        randomSkew.Should().BeInRange(-1.0, 1.0,
            because: "|G−C| ≤ G+C, so the skew can never escape [−1, 1], even at scale");
    }

    #endregion

    #region BE/RB — range invariant on fuzzed-but-valid inputs (both surfaces)

    /// <summary>
    /// BE/RB: the [−1, 1] range invariant (INV-01) and the NaN-free / non-throwing
    /// contract must hold for ANY input on both surfaces. Locally-seeded (so the
    /// shared static Rng is untouched) sweeps feed (a) valid random DNA to the strict
    /// typed path and (b) arbitrary BMP code points — control chars, null bytes, lone
    /// surrogate halves — to the lenient raw-string path. Every result must be finite
    /// and inside [−1, 1]; neither surface may throw on its allowed input. This pins
    /// the strict/lenient boundary: the same random bytes the typed path rejects at
    /// construction, the lenient overload carries through to a defined finite skew.
    /// </summary>
    [Test]
    public void GcSkew_RandomSweeps_AreFiniteInRangeAndNeverThrow()
    {
        var local = new Random(99173);
        const string bases = "ACGT";

        for (int trial = 0; trial < 300; trial++)
        {
            // (a) Strict typed path over valid random DNA.
            int validLen = local.Next(0, 64);
            var validChars = new char[validLen];
            for (int i = 0; i < validLen; i++)
                validChars[i] = bases[local.Next(bases.Length)];
            string validDna = new(validChars);

            double typedSkew = double.NaN;
            var typedAct = () => typedSkew = GcSkewCalculator.CalculateGcSkew(new DnaSequence(validDna));
            typedAct.Should().NotThrow($"valid DNA must never make the typed path throw; input: \"{validDna}\"");
            double.IsFinite(typedSkew).Should().BeTrue($"skew is never NaN/Inf; input: \"{validDna}\"");
            typedSkew.Should().BeInRange(-1.0, 1.0, $"skew is bounded to [−1,1]; input: \"{validDna}\"");

            // (b) Lenient raw-string path over arbitrary BMP garbage.
            int garbageLen = local.Next(0, 64);
            var garbageChars = new char[garbageLen];
            for (int i = 0; i < garbageLen; i++)
                garbageChars[i] = (char)local.Next(0x0000, 0x10000);
            string garbage = new(garbageChars);

            double rawSkew = double.NaN;
            var rawAct = () => rawSkew = GcSkewCalculator.CalculateGcSkew(garbage);
            rawAct.Should().NotThrow("the lenient overload must never throw on random bytes");
            double.IsFinite(rawSkew).Should().BeTrue(
                "only G/C are counted, so no NaN/Inf can leak from the lenient overload");
            rawSkew.Should().BeInRange(-1.0, 1.0,
                "the lenient skew is bounded to [−1,1] even on pure garbage");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-ATSKEW-001 — AT skew : (A − T) / (A + T), range [−1, 1]
    //  Typed overload:  CalculateAtSkew(DnaSequence) — strict, ArgumentNull on null
    //  Raw-string:      CalculateAtSkew(string)      — lenient, counts only A/T
    //  AT_Skew.md §2.2/§2.4/§3.3/§6.1/§7.1; ADVANCED_TESTING_CHECKLIST.md §8.
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-ATSKEW-001 — AT skew

    #region Positive sanity — worked examples from the doc (§7.1)

    /// <summary>
    /// Positive sanity: the exact, hand-checkable worked examples from AT_Skew.md
    /// §7.1 and §6.1, derived INDEPENDENTLY from the formula (A−T)/(A+T), NOT read
    /// off the implementation. "AAAT" → A=3,T=1 → 2/4 = 0.5; the numerical
    /// walk-through "AAATGGGCCC" → A=3,T=1 (G/C ignored, INV-05) → 0.5; pure A → +1;
    /// pure T → −1; balanced 2A/2T → 0 (INV-02); a mixed count "AATTTT" →
    /// (2−4)/6 = −1/3. Verified on the strict typed path over valid A/C/G/T DNA.
    /// </summary>
    [TestCase("AAAT", 0.5, TestName = "AtSkew_Doc_AAAT_IsHalf")]
    [TestCase("AAATGGGCCC", 0.5, TestName = "AtSkew_Doc_AAATGGGCCC_GcIgnored_IsHalf")]
    [TestCase("AAAA", +1.0, TestName = "AtSkew_Doc_PureA_IsPlusOne")]
    [TestCase("TTTT", -1.0, TestName = "AtSkew_Doc_PureT_IsMinusOne")]
    [TestCase("AATT", 0.0, TestName = "AtSkew_Doc_Balanced2A2T_IsZero")]
    [TestCase("AATTTT", -1.0 / 3.0, TestName = "AtSkew_Doc_TwoAFourT_IsMinusThird")]
    public void AtSkew_TypedWorkedExamples_MatchFormula(string input, double expected)
    {
        double skew = GcSkewCalculator.CalculateAtSkew(new DnaSequence(input));

        skew.Should().BeApproximately(expected, Tolerance,
            because: $"(A−T)/(A+T) for \"{input}\" is {expected} per AT_Skew.md §7.1/§6.1; "
                   + "G/C are ignored (INV-05) and the skew stays in the closed range [−1, 1]");
    }

    #endregion

    #region BE — Boundary: empty and null

    /// <summary>
    /// BE: the empty string is the lower size boundary AND a face of the
    /// DivideByZero boundary (A + T = 0). The typed path defines the empty
    /// DnaSequence (DnaSequence.cs lines 24–28) and CalculateAtSkew must return the
    /// defined 0 via the zero-denominator guard (GcSkewCalculator.cs line 205) —
    /// no DivideByZeroException, no NaN, no crash (AT_Skew.md INV-03 / §6.1). The
    /// lenient raw-string overload returns 0 for "" and null by its IsNullOrEmpty
    /// short-circuit (lines 191–192).
    /// </summary>
    [Test]
    public void AtSkew_EmptyAndNull_IsZeroAndDoesNotThrow()
    {
        var typed = () =>
        {
            double skew = GcSkewCalculator.CalculateAtSkew(new DnaSequence(string.Empty));
            skew.Should().Be(0.0,
                because: "an empty sequence has no A or T; skew is defined as 0 by the zero-denominator guard");
        };
        typed.Should().NotThrow("the empty sequence is a defined boundary, not a DivideByZero error");

        var rawEmpty = () => GcSkewCalculator.CalculateAtSkew(string.Empty).Should().Be(0.0);
        var rawNull = () => GcSkewCalculator.CalculateAtSkew((string)null!).Should().Be(0.0);
        rawEmpty.Should().NotThrow("the lenient overload short-circuits empty input to 0");
        rawNull.Should().NotThrow<NullReferenceException>(
            "null is handled by the IsNullOrEmpty gate, never dereferenced");
        rawNull.Should().NotThrow("null is a defined 'no input' boundary returning 0, not an error");
    }

    /// <summary>
    /// BE: a null DnaSequence is the explicit-guard boundary on the strict typed
    /// path. CalculateAtSkew(DnaSequence) must throw the *documented, intentional*
    /// ArgumentNullException (GcSkewCalculator.cs line 176, AT_Skew.md §3.3), never
    /// a raw NullReferenceException from dereferencing <c>.Sequence</c>.
    /// </summary>
    [Test]
    public void AtSkew_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => GcSkewCalculator.CalculateAtSkew((DnaSequence)null!);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload guards null with the documented ArgumentNullException, not a crash");
    }

    #endregion

    #region BE — Boundary: all-A and the single-base extremes

    /// <summary>
    /// BE: the checklist's "all-A" target and the binary extreme of the skew. A
    /// one-base 'A' is the maximum +1 (A=1,T=0 → 1/1); a single 'T' is the minimum
    /// −1; a single 'G' or 'C' has NO A/T, so the A+T=0 guard returns the defined 0
    /// (NOT NaN). Lowercase is accepted on the typed path (case-folded at
    /// construction, INV-04) and yields the same value. Verified over all four bases
    /// in both cases on the strict typed overload.
    /// </summary>
    [TestCase('A', +1.0)]
    [TestCase('T', -1.0)]
    [TestCase('G', 0.0)]
    [TestCase('C', 0.0)]
    [TestCase('a', +1.0)]
    [TestCase('t', -1.0)]
    [TestCase('g', 0.0)]
    [TestCase('c', 0.0)]
    public void AtSkew_SingleBase_IsBinaryExtremeOrZero(char baseChar, double expected)
    {
        double skew = GcSkewCalculator.CalculateAtSkew(new DnaSequence(baseChar.ToString()));

        skew.Should().BeApproximately(expected, Tolerance,
            because: expected switch
            {
                +1.0 => $"a single '{baseChar}' is all-A → maximum skew +1",
                -1.0 => $"a single '{baseChar}' is all-T → minimum skew −1",
                _ => $"a single '{baseChar}' has no A/T, so the A+T=0 guard returns the defined 0, not NaN"
            });
    }

    /// <summary>
    /// BE: the "all-A" boundary across lengths and the all-T mirror. An all-A run is
    /// the maximum +1 and an all-T run is the minimum −1 at ANY length (T=0 ⇒ +1,
    /// A=0 ⇒ −1, INV-01). Pins that the extreme value does not drift with length.
    /// </summary>
    [TestCase("A", +1.0, TestName = "AtSkew_AllA_SingleA_IsPlusOne")]
    [TestCase("AAAA", +1.0, TestName = "AtSkew_AllA_Length4_IsPlusOne")]
    [TestCase("AAAAAAAAAAAA", +1.0, TestName = "AtSkew_AllA_Length12_IsPlusOne")]
    [TestCase("TTTT", -1.0, TestName = "AtSkew_AllT_Length4_IsMinusOne")]
    public void AtSkew_AllAOrAllT_IsExtreme(string input, double expected)
    {
        double skew = GcSkewCalculator.CalculateAtSkew(new DnaSequence(input));

        skew.Should().BeApproximately(expected, Tolerance,
            because: $"an all-{input[0]} run has skew {expected} regardless of length");
    }

    #endregion

    #region BE — Boundary: balanced AT (A = T ⇒ skew 0)

    /// <summary>
    /// BE: the checklist's "balanced AT" target — the symmetry boundary. A sequence
    /// with EQUAL A and T counts (alternating ATAT… and its TATA… phase, or any
    /// balanced mix) has numerator A−T = 0, so the skew is exactly 0 regardless of
    /// length and regardless of ignored G/C content (INV-02, INV-05). Verified on the
    /// strict typed path over valid A/C/G/T DNA across lengths and phases.
    /// </summary>
    [TestCase("AT", TestName = "AtSkew_Balanced_AT_IsZero")]
    [TestCase("TA", TestName = "AtSkew_Balanced_TA_IsZero")]
    [TestCase("ATATATAT", TestName = "AtSkew_Balanced_ATATATAT_IsZero")]
    [TestCase("AATT", TestName = "AtSkew_Balanced_AATT_IsZero")]
    [TestCase("ATGCAT", TestName = "AtSkew_Balanced_WithGcIgnored_IsZero")]
    [TestCase("AAGGTTCC", TestName = "AtSkew_Balanced_TwoAtwoT_GcIgnored_IsZero")]
    public void AtSkew_BalancedAt_IsZero(string input)
    {
        double skew = GcSkewCalculator.CalculateAtSkew(new DnaSequence(input));

        skew.Should().Be(0.0,
            because: $"\"{input}\" has equal A and T (G/C ignored), so A−T=0 and the skew is exactly 0");
    }

    #endregion

    #region BE — Boundary: no A or T (the DivideByZero boundary, A + T = 0)

    /// <summary>
    /// BE: the checklist's "no AT" target — the critical A + T = 0 boundary. A
    /// sequence with NO adenine and NO thymine (all-G, all-C, GC-only) makes the
    /// denominator zero. The typed path over valid G/C-only DNA must return the
    /// defined 0 — never DivideByZeroException, never NaN. This is the headline fuzz
    /// target for the AT-skew DivideByZero face (AT_Skew.md INV-03 / §6.1).
    /// </summary>
    [TestCase("G", TestName = "AtSkew_NoAorT_SingleG_IsZero")]
    [TestCase("GGGG", TestName = "AtSkew_NoAorT_AllG_IsZero")]
    [TestCase("CCCC", TestName = "AtSkew_NoAorT_AllC_IsZero")]
    [TestCase("GGCC", TestName = "AtSkew_NoAorT_GcOnly_IsZero")]
    [TestCase("GCGCGC", TestName = "AtSkew_NoAorT_GcAlternating_IsZero")]
    public void AtSkew_TypedNoAorT_IsZeroAndDoesNotThrow(string gcOnly)
    {
        double skew = double.NaN;
        var act = () => skew = GcSkewCalculator.CalculateAtSkew(new DnaSequence(gcOnly));

        act.Should().NotThrow("no A/T means A+T=0, the guarded boundary — not a DivideByZero crash");
        double.IsNaN(skew).Should().BeFalse("the zero-denominator guard returns 0, never a 0/0 NaN");
        skew.Should().Be(0.0, because: "with no A and no T the skew is defined as 0");
    }

    /// <summary>
    /// BE: the lenient raw-string path at the A+T=0 boundary. Input that contains NO
    /// A and NO T after upper-casing — pure G/C, ambiguity codes, digits, gaps, an
    /// embedded null byte, unicode, or the RNA base U (NOT converted to T, AT_Skew.md
    /// §6.2) — counts no A/T, so the total==0 guard returns the defined 0
    /// (GcSkewCalculator.cs line 205). The overload must NEVER throw and never leak
    /// NaN, even on garbage the typed path would reject.
    /// </summary>
    [TestCase("GGCC", TestName = "AtSkew_RawNoAorT_GcOnly_IsZero")]
    [TestCase("GCGC", TestName = "AtSkew_RawNoAorT_GcAlternating_IsZero")]
    [TestCase("NNNN", TestName = "AtSkew_RawNoAorT_OnlyAmbiguity_IsZero")]
    [TestCase("12345", TestName = "AtSkew_RawNoAorT_OnlyDigits_IsZero")]
    [TestCase("----", TestName = "AtSkew_RawNoAorT_OnlyGaps_IsZero")]
    [TestCase("UUUU", TestName = "AtSkew_RawNoAorT_RnaUNotTreatedAsT_IsZero")]
    [TestCase("g\0c", TestName = "AtSkew_RawNoAorT_NullByteBetweenGc_IsZero")]
    [TestCase("ααα", TestName = "AtSkew_RawNoAorT_Unicode_IsZero")]
    public void AtSkew_RawStringNoAorT_IsZeroAndDoesNotThrow(string input)
    {
        double skew = double.NaN;
        var act = () => skew = GcSkewCalculator.CalculateAtSkew(input);

        act.Should().NotThrow("the lenient overload counts only A/T and never throws on non-A/T garbage");
        double.IsNaN(skew).Should().BeFalse("the total==0 guard returns 0, never a 0/0 NaN");
        skew.Should().Be(0.0, because: "no A and no T is present, so the skew is the defined 0");
    }

    /// <summary>
    /// INV-05 boundary: the lenient overload IGNORES G/C and other non-A/T symbols
    /// and computes the skew of ONLY the A/T bases present — so injected garbage
    /// interspersed with real A/T must not shift the value. "A-G_T1N" has one A and
    /// one T → balanced → 0; "aXgtT" (lowercase, INV-04) has one A and two T →
    /// (1−2)/3 = −1/3; the doc's own "AAATGGGCCC" reads 0.5 on this surface too.
    /// Pins that injection cannot silently corrupt the count.
    /// </summary>
    [TestCase("A-G_T1N", 0.0, TestName = "AtSkew_RawIgnoresGarbage_BalancedAt_IsZero")]
    [TestCase("aXgtT", -1.0 / 3.0, TestName = "AtSkew_RawIgnoresGarbage_OneATwoT_IsMinusThird")]
    [TestCase("AAATGGGCCC", 0.5, TestName = "AtSkew_RawIgnoresGc_DocExample_IsHalf")]
    [TestCase("A\0\0\0", +1.0, TestName = "AtSkew_RawIgnoresNullBytes_PureA_IsPlusOne")]
    public void AtSkew_RawStringIgnoresNonAtSymbols(string input, double expected)
    {
        double skew = GcSkewCalculator.CalculateAtSkew(input);

        skew.Should().BeApproximately(expected, Tolerance,
            because: $"only A/T are counted (INV-05); ignored symbols do not shift the skew of \"{input}\"");
    }

    #endregion

    #region BE — Boundary: window edge (caller substring windows; no native window API)

    /// <summary>
    /// BE — "window edge": AT_Skew.md §5.3/§6.2 are explicit that this unit is a
    /// single GLOBAL scalar with NO native windowing API. We honour the checklist's
    /// window-edge boundary the only spec-faithful way: AT skew computed over a
    /// CALLER substring window must equal the global formula on that exact slice and
    /// stay in [−1, 1] with no crash/NaN. This pins the three window-edge faces:
    ///   • window length 1 — each single-base window is +1/−1 (A/T) or 0 (G/C);
    ///   • window LONGER than the sequence (clamped by the caller to the whole seq)
    ///     — equals the global skew, no over-read / IndexOutOfRange;
    ///   • the partial trailing window — the last window shorter than the step still
    ///     computes correctly from its own A/T counts.
    /// Locally-seeded RNG (never the shared static Rng) so sibling fixtures are
    /// untouched. [CancelAfter] guards the sweep against a hang.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void AtSkew_CallerWindows_MatchGlobalFormulaAtWindowEdges()
    {
        var local = new Random(20260620);
        const string bases = "ACGT";

        for (int trial = 0; trial < 200; trial++)
        {
            int n = local.Next(1, 40);
            var chars = new char[n];
            for (int i = 0; i < n; i++)
                chars[i] = bases[local.Next(bases.Length)];
            string seq = new(chars);

            // window length sweeps the three edges: 1, every interior size, and
            // a window strictly LONGER than the sequence (clamped to the whole seq).
            foreach (int window in new[] { 1, Math.Max(1, n / 2), n, n + 5, n + 1000 })
            {
                for (int start = 0; start < n; start += window)
                {
                    int len = Math.Min(window, n - start);          // partial trailing window
                    string slice = seq.Substring(start, len);

                    // Independent reference: the §2.2 formula on this exact slice.
                    int a = slice.Count(c => c == 'A');
                    int t = slice.Count(c => c == 'T');
                    double expected = (a + t) > 0 ? (double)(a - t) / (a + t) : 0.0;

                    double skew = double.NaN;
                    var act = () => skew = GcSkewCalculator.CalculateAtSkew(slice);
                    act.Should().NotThrow(
                        $"a caller window must never crash the global statistic; "
                        + $"slice=\"{slice}\" window={window} start={start}");
                    double.IsNaN(skew).Should().BeFalse(
                        $"the A+T=0 guard prevents NaN even on a G/C-only window; slice=\"{slice}\"");
                    skew.Should().BeInRange(-1.0, 1.0,
                        $"|A−T| ≤ A+T, so a window skew stays in [−1,1]; slice=\"{slice}\"");
                    skew.Should().BeApproximately(expected, Tolerance,
                        $"AT skew on a caller window equals (A−T)/(A+T) of that slice; slice=\"{slice}\"");
                }
            }
        }
    }

    #endregion

    #region BE — Boundary: extremely long

    /// <summary>
    /// BE/OVF: an extremely long valid sequence (1,000,000 bases) must compute the
    /// skew without overflow, hang, or precision blow-up; the result must stay in
    /// [−1, 1] and be finite. Known-composition long inputs pin exact values: "AG"
    /// repeated has A but no T → +1 at any length; "AT" repeated has equal A and T →
    /// 0 at any length; "GC" repeated has no A/T → the A+T=0 guard returns 0. A
    /// fixed-seed random long sequence pins the range invariant at scale.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void AtSkew_ExtremelyLong_StaysInRangeAndDoesNotHang()
    {
        const int length = 1_000_000;

        // "AG" repeated: A present, T absent → skew = (A−0)/(A+0) = +1, at any length.
        var allAskew = new DnaSequence(string.Concat(Enumerable.Repeat("AG", length / 2)));
        allAskew.Length.Should().Be(length);
        GcSkewCalculator.CalculateAtSkew(allAskew).Should().BeApproximately(+1.0, Tolerance,
            because: "a long sequence with A but no T has skew +1, regardless of length");

        // "AT" repeated: equal A and T → skew = 0, at any length.
        var balanced = new DnaSequence(string.Concat(Enumerable.Repeat("AT", length / 2)));
        GcSkewCalculator.CalculateAtSkew(balanced).Should().BeApproximately(0.0, Tolerance,
            because: "a long sequence with equal A and T has skew 0, regardless of length");

        // "GC" repeated: no A/T → the A+T=0 guard returns the defined 0.
        var noAt = new DnaSequence(string.Concat(Enumerable.Repeat("GC", length / 2)));
        GcSkewCalculator.CalculateAtSkew(noAt).Should().Be(0.0,
            because: "a long sequence with no A and no T is the A+T=0 boundary, defined as 0");

        // Fixed-seed random long sequence (locally-seeded so other fixtures are not
        // perturbed by the shared static Rng): a finite value in the closed range.
        var local = new Random(48211);
        var chars = new char[length];
        const string bases = "ACGT";
        for (int i = 0; i < length; i++)
            chars[i] = bases[local.Next(bases.Length)];

        double randomSkew = GcSkewCalculator.CalculateAtSkew(new DnaSequence(new string(chars)));
        double.IsFinite(randomSkew).Should().BeTrue("the skew can never be NaN/Inf, even at scale");
        randomSkew.Should().BeInRange(-1.0, 1.0,
            because: "|A−T| ≤ A+T, so the skew can never escape [−1, 1], even at scale");
    }

    #endregion

    #region BE/RB — range invariant on fuzzed-but-valid inputs (both surfaces)

    /// <summary>
    /// BE/RB: the [−1, 1] range invariant (INV-01) and the NaN-free / non-throwing
    /// contract must hold for ANY input on both surfaces. Locally-seeded (so the
    /// shared static Rng is untouched) sweeps feed (a) valid random DNA to the strict
    /// typed path and (b) arbitrary BMP code points — control chars, null bytes, lone
    /// surrogate halves — to the lenient raw-string path. Every result must be finite
    /// and inside [−1, 1]; neither surface may throw on its allowed input. This pins
    /// the strict/lenient boundary: the same random bytes the typed path rejects at
    /// construction, the lenient overload carries through to a defined finite skew.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void AtSkew_RandomSweeps_AreFiniteInRangeAndNeverThrow()
    {
        var local = new Random(33107);
        const string bases = "ACGT";

        for (int trial = 0; trial < 300; trial++)
        {
            // (a) Strict typed path over valid random DNA.
            int validLen = local.Next(0, 64);
            var validChars = new char[validLen];
            for (int i = 0; i < validLen; i++)
                validChars[i] = bases[local.Next(bases.Length)];
            string validDna = new(validChars);

            double typedSkew = double.NaN;
            var typedAct = () => typedSkew = GcSkewCalculator.CalculateAtSkew(new DnaSequence(validDna));
            typedAct.Should().NotThrow($"valid DNA must never make the typed path throw; input: \"{validDna}\"");
            double.IsFinite(typedSkew).Should().BeTrue($"skew is never NaN/Inf; input: \"{validDna}\"");
            typedSkew.Should().BeInRange(-1.0, 1.0, $"skew is bounded to [−1,1]; input: \"{validDna}\"");

            // (b) Lenient raw-string path over arbitrary BMP garbage.
            int garbageLen = local.Next(0, 64);
            var garbageChars = new char[garbageLen];
            for (int i = 0; i < garbageLen; i++)
                garbageChars[i] = (char)local.Next(0x0000, 0x10000);
            string garbage = new(garbageChars);

            double rawSkew = double.NaN;
            var rawAct = () => rawSkew = GcSkewCalculator.CalculateAtSkew(garbage);
            rawAct.Should().NotThrow($"the lenient overload must never throw on random bytes; input: {Describe(garbage)}");
            double.IsFinite(rawSkew).Should().BeTrue(
                $"only A/T are counted, so no NaN/Inf can leak; input: {Describe(garbage)}");
            rawSkew.Should().BeInRange(-1.0, 1.0,
                $"the lenient skew is bounded to [−1,1] even on pure garbage; input: {Describe(garbage)}");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-REPLICATION-001 — replication origin/terminus prediction via the
    //  cumulative GC-skew minimum/maximum (Composition).
    //  Checklist: docs/checklists/03_FUZZING.md, row 211 (BE; "flat skew,
    //  single minimum, circular wrap").
    //
    //  THEORY (Replication_Origin_Prediction.md, Test Unit ID SEQ-REPLICATION-001):
    //    Cumulative skew Skew_0 = 0, Skew_{i+1} = Skew_i + s(Genome[i]) with
    //    s(G)=+1, s(C)=−1, s(A)=s(T)=other=0 (§2.2 [2]). There are n+1 prefix
    //    values Skew_0…Skew_n over 0-based prefix indices i ∈ [0, n] (position i =
    //    the boundary BEFORE base i, Rosalind BA1F). The predicted ORIGIN is the
    //    FIRST prefix index minimizing Skew_i; the predicted TERMINUS is the FIRST
    //    prefix index maximizing it (§2.4 INV-01/INV-02, ties → smallest index).
    //    OriginSkew ≤ 0 ≤ TerminusSkew because Skew_0 = 0 is always a prefix value
    //    (INV-03). IsSignificant ⇔ max > min (INV-05): a FLAT diagram (no net G/C
    //    asymmetry) carries no origin signal. A and T never move the diagram
    //    (INV-06). O(n) single pass (§4.3). Typed null → ArgumentNullException;
    //    null/empty string → zero prediction, not significant (§3.3, §6.1).
    //
    //  Fuzz strategy for THIS unit: BE = Boundary Exploitation.
    //    • FLAT SKEW — no G/C asymmetry (empty, no-G/C, equal-and-balanced G=C,
    //      A/T-only). The diagram never leaves 0, so origin = terminus = 0,
    //      OriginSkew = TerminusSkew = 0, IsSignificant = false. No crash, no NaN.
    //    • SINGLE CLEAR MINIMUM — a sequence with one unambiguous global minimum;
    //      the predicted origin index (and the symmetric maximum/terminus) are
    //      recovered EXACTLY, derived independently from the §2.2 prefix recurrence,
    //      including the Rosalind BA1F sample (first minimizer 53, §7.1).
    //    • CIRCULAR WRAP — the molecule is circular (§1, ASM-01/ASM-02): the true
    //      origin can sit ACROSS the linearization join. We pin the spec-faithful
    //      consequence (§6.2): the unit treats the input as a LINEAR string in
    //      genome coordinates, so ROTATING the sequence shifts the predicted origin
    //      by the rotation; locating the cumulative minimum near / across the end
    //      boundary must still recover the correct first-minimizer index with no
    //      over-read / IndexOutOfRange.
    //  — docs/checklists/03_FUZZING.md §Description (BE); ADVANCED_TESTING_CHECKLIST.md §8.
    //
    //  Surfaces (Replication_Origin_Prediction.md §5.1):
    //    • PredictReplicationOrigin(DnaSequence) — strict; ArgumentNullException on null.
    //    • PredictReplicationOrigin(string)      — lenient; null/empty → zero prediction;
    //                                              only G/C affect the skew.
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-REPLICATION-001 — replication origin prediction

    #region Helpers — independent reference computed straight from §2.2

    /// <summary>
    /// Independent reference implementation of the cumulative-skew origin/terminus
    /// derived DIRECTLY from Replication_Origin_Prediction.md §2.2 / §2.4 — NEVER
    /// read off GcSkewCalculator. Skew_0 = 0; each base adds +1 for G, −1 for C, 0
    /// otherwise; the origin is the FIRST (smallest) prefix index minimizing the
    /// running skew and the terminus the FIRST index maximizing it (ties → smallest
    /// index, via strict comparison). Returns (origin, terminus, minSkew, maxSkew).
    /// </summary>
    private static (int Origin, int Terminus, int MinSkew, int MaxSkew) ReferenceOrigin(string seq)
    {
        int cumulative = 0;                 // Skew_0
        int minSkew = 0, maxSkew = 0;       // Skew_0 = 0 is always a prefix value (INV-03)
        int minPos = 0, maxPos = 0;         // first-minimizer / first-maximizer (ties → smallest)

        for (int i = 0; i < seq.Length; i++)
        {
            char c = char.ToUpperInvariant(seq[i]);
            if (c == 'G') cumulative += 1;
            else if (c == 'C') cumulative -= 1;
            // A, T and any other symbol leave the running skew unchanged (INV-06).

            int prefixIndex = i + 1;        // Skew_{i+1} after consuming base i
            if (cumulative < minSkew) { minSkew = cumulative; minPos = prefixIndex; }
            if (cumulative > maxSkew) { maxSkew = cumulative; maxPos = prefixIndex; }
        }

        return (minPos, maxPos, minSkew, maxSkew);
    }

    #endregion

    #region Positive sanity — worked examples from the doc (§7.1)

    /// <summary>
    /// Positive sanity: the EXACT, hand-checkable worked examples from
    /// Replication_Origin_Prediction.md §7.1, with expected values derived
    /// INDEPENDENTLY from the §2.2 prefix recurrence (NOT read off the code).
    ///
    /// • "CCGGGG" → diagram 0,−1,−2,−1,0,+1,+2 → min −2 first at prefix index 2
    ///   (origin), max +2 at index 6 (terminus); amplitude 4 ⇒ IsSignificant.
    /// • Single base "G" → diagram 0,+1 → origin 0 (skew 0), terminus 1 (skew +1).
    /// • Tied extremum "CCGGCC" → diagram 0,−1,−2,−1,0,−1,−2 → min −2 FIRST at index
    ///   2 (origin); the max value 0 is first seen at Skew_0, so terminus = 0;
    ///   amplitude 2 ⇒ significant. Pins the first-index tie-break (INV-01/INV-02).
    /// </summary>
    [TestCase("CCGGGG", 2, 6, -2.0, +2.0, true, TestName = "Replication_Doc_CCGGGG_Origin2Terminus6")]
    [TestCase("G", 0, 1, 0.0, +1.0, true, TestName = "Replication_Doc_SingleG_Origin0Terminus1")]
    [TestCase("CCGGCC", 2, 0, -2.0, 0.0, true, TestName = "Replication_Doc_CCGGCC_TieBreakFirstIndex")]
    public void Replication_TypedWorkedExamples_MatchPrefixRecurrence(
        string input, int origin, int terminus, double originSkew, double terminusSkew, bool significant)
    {
        var pred = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(input));

        pred.PredictedOrigin.Should().Be(origin,
            because: $"the cumulative-skew minimum of \"{input}\" first occurs at prefix index {origin} (§2.2/§7.1)");
        pred.PredictedTerminus.Should().Be(terminus,
            because: $"the cumulative-skew maximum of \"{input}\" first occurs at prefix index {terminus} (§2.2/§7.1)");
        pred.OriginSkew.Should().BeApproximately(originSkew, Tolerance,
            because: "OriginSkew is the minimum cumulative value and is ≤ 0 (INV-03)");
        pred.TerminusSkew.Should().BeApproximately(terminusSkew, Tolerance,
            because: "TerminusSkew is the maximum cumulative value and is ≥ 0 (INV-03)");
        pred.IsSignificant.Should().Be(significant,
            because: "IsSignificant ⇔ max > min (INV-05)");
    }

    /// <summary>
    /// Positive sanity: the Rosalind BA1F sample genome (§7.1 / §5.3). The doc pins
    /// the published answer — first minimizer at prefix index 53, OriginSkew −4 — and
    /// the maximizing terminus is verified against the independent §2.2 reference,
    /// not the code. This is the canonical "single clear minimum" anchor.
    /// </summary>
    [Test]
    public void Replication_RosalindBA1FSample_OriginIs53()
    {
        const string genome =
            "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG";

        var pred = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(genome));
        var reference = ReferenceOrigin(genome);

        pred.PredictedOrigin.Should().Be(53,
            because: "Rosalind BA1F publishes the first minimizer of this genome as position 53 (§5.3/§7.1)");
        pred.OriginSkew.Should().BeApproximately(-4.0, Tolerance,
            because: "the doc's worked example pins OriginSkew == −4 for the BA1F sample (§7.1)");
        pred.PredictedTerminus.Should().Be(reference.Terminus,
            because: "the terminus is the first maximizer of the §2.2 prefix diagram (INV-02)");
        pred.TerminusSkew.Should().BeApproximately(reference.MaxSkew, Tolerance);
        pred.IsSignificant.Should().BeTrue("the BA1F genome has a non-flat diagram (max > min)");
    }

    #endregion

    #region BE — Boundary: flat skew (no net G/C asymmetry ⇒ degenerate origin)

    /// <summary>
    /// BE — "flat skew": when the diagram never leaves 0 there is NO origin signal.
    /// This happens for (a) A/T-only input (s(A)=s(T)=0, INV-06: G/C absent), and
    /// (b) sequences whose running skew returns to and stays bounded at 0 only when
    /// G and C are perfectly interleaved so the cumulative value never dips below or
    /// rises above 0 (e.g. "GCGCGC": 0,+1,0,+1,0,+1,0 — max +1 > min 0 IS a signal,
    /// so that is NOT flat). A genuinely FLAT diagram requires the cumulative value
    /// to stay EXACTLY 0 throughout — i.e. no G/C at all. For such input the contract
    /// (§6.1) is origin = terminus = 0, OriginSkew = TerminusSkew = 0,
    /// IsSignificant = FALSE — a defined, non-throwing degenerate result. Expected
    /// origin/terminus are derived from the §2.2 reference, not the code.
    /// </summary>
    [TestCase("AAAATTTT", TestName = "Replication_Flat_AtOnly_NoOriginSignal")]
    [TestCase("AAAA", TestName = "Replication_Flat_AllA_NoOriginSignal")]
    [TestCase("TTTTTTTT", TestName = "Replication_Flat_AllT_NoOriginSignal")]
    [TestCase("ATATATAT", TestName = "Replication_Flat_AlternatingAT_NoOriginSignal")]
    public void Replication_FlatDiagram_OriginAndTerminusAreZeroNotSignificant(string input)
    {
        var reference = ReferenceOrigin(input);
        reference.MinSkew.Should().Be(0, "by construction this input has no G/C, so the diagram stays at 0");
        reference.MaxSkew.Should().Be(0, "by construction this input has no G/C, so the diagram stays at 0");

        ReplicationOriginPrediction pred = default;
        var act = () => pred = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(input));

        act.Should().NotThrow("a flat diagram is a defined degenerate case, not an error");
        pred.PredictedOrigin.Should().Be(0, "a flat diagram pins the origin at Skew_0 (prefix index 0)");
        pred.PredictedTerminus.Should().Be(0, "a flat diagram pins the terminus at Skew_0 (prefix index 0)");
        pred.OriginSkew.Should().Be(0.0);
        pred.TerminusSkew.Should().Be(0.0);
        pred.IsSignificant.Should().BeFalse(
            "a flat diagram (max == min == 0) carries no origin signal (INV-05)");
    }

    /// <summary>
    /// BE — "flat skew" at the size floor: the EMPTY sequence. The typed path defines
    /// an empty DnaSequence (DnaSequence.cs lines 24–28); the contract (§6.1) is the
    /// zero prediction (origin = terminus = 0, skews 0, not significant) with NO
    /// division and NO exception. The lenient raw-string overload returns the same
    /// zero prediction for "" and for null via its IsNullOrEmpty short-circuit
    /// (GcSkewCalculator.cs lines 248–250).
    /// </summary>
    [Test]
    public void Replication_EmptyAndNullString_AreZeroPredictionAndDoNotThrow()
    {
        ReplicationOriginPrediction typed = default;
        var typedAct = () => typed = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(string.Empty));
        typedAct.Should().NotThrow("the empty sequence is a defined boundary, not an error");
        typed.Should().Be(new ReplicationOriginPrediction(0, 0, 0, 0, false),
            because: "an empty diagram has the single prefix value Skew_0 = 0 → zero prediction (§6.1)");

        var rawEmpty = () => GcSkewCalculator.PredictReplicationOrigin(string.Empty)
            .Should().Be(new ReplicationOriginPrediction(0, 0, 0, 0, false));
        var rawNull = () => GcSkewCalculator.PredictReplicationOrigin((string)null!)
            .Should().Be(new ReplicationOriginPrediction(0, 0, 0, 0, false));
        rawEmpty.Should().NotThrow("the lenient overload short-circuits empty input to the zero prediction");
        rawNull.Should().NotThrow<NullReferenceException>(
            "null is handled by the IsNullOrEmpty gate, never dereferenced");
        rawNull.Should().NotThrow("null/empty string is a defined 'no input' boundary, not an error");
    }

    /// <summary>
    /// BE: a null DnaSequence is the explicit-guard boundary on the strict typed path.
    /// PredictReplicationOrigin(DnaSequence) must throw the *documented, intentional*
    /// ArgumentNullException (GcSkewCalculator.cs line 238, §3.3), never a raw
    /// NullReferenceException from dereferencing <c>.Sequence</c>.
    /// </summary>
    [Test]
    public void Replication_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => GcSkewCalculator.PredictReplicationOrigin((DnaSequence)null!);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload guards null with the documented ArgumentNullException, not a crash");
    }

    #endregion

    #region BE — Boundary: single clear minimum (origin recovered exactly)

    /// <summary>
    /// BE — "single clear minimum": a sequence engineered to have ONE unambiguous
    /// global minimum (a falling C-run followed by a longer rising G-run). The origin
    /// is the prefix index at the bottom of the C-run and the terminus the index at
    /// the top of the G-run, BOTH derived independently from the §2.2 recurrence.
    ///
    /// • "CCCGGGGG": 0,−1,−2,−3,−2,−1,0,+1,+2 → origin 3 (skew −3), terminus 8 (+2).
    /// • "ATCCCGGAT" (A/T padding, INV-06): G/C-relevant subsequence C C C G G →
    ///   0,0,0,−1,−2,−3,−2,−1,−1,−1 → origin 6 (skew −3), terminus 0 (max stays 0).
    /// • "GGGCCCCCC": 0,+1,+2,+3,+2,+1,0,−1,−2,−3 → origin 9 (skew −3), terminus 3 (+3).
    /// The exact origin/terminus indices and skews are pinned against the reference.
    /// </summary>
    [TestCase("CCCGGGGG", TestName = "Replication_SingleMin_FallingThenRising")]
    [TestCase("ATCCCGGAT", TestName = "Replication_SingleMin_AtPaddingIgnored")]
    [TestCase("GGGCCCCCC", TestName = "Replication_SingleMin_RisingThenFalling")]
    [TestCase("AGCTAGCTGGGGCCCCATAT", TestName = "Replication_SingleMin_MixedGenomeLike")]
    public void Replication_SingleClearMinimum_RecoversOriginExactly(string input)
    {
        var reference = ReferenceOrigin(input);

        var pred = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(input));

        pred.PredictedOrigin.Should().Be(reference.Origin,
            because: $"the cumulative minimum of \"{input}\" first occurs at prefix index {reference.Origin} (§2.2 INV-01)");
        pred.PredictedTerminus.Should().Be(reference.Terminus,
            because: $"the cumulative maximum of \"{input}\" first occurs at prefix index {reference.Terminus} (§2.2 INV-02)");
        pred.OriginSkew.Should().BeApproximately(reference.MinSkew, Tolerance,
            because: "OriginSkew is the minimum cumulative value (INV-03: ≤ 0)");
        pred.TerminusSkew.Should().BeApproximately(reference.MaxSkew, Tolerance,
            because: "TerminusSkew is the maximum cumulative value (INV-03: ≥ 0)");
        pred.OriginSkew.Should().BeLessThanOrEqualTo(0.0, "Skew_0 = 0 is always a prefix value, so min ≤ 0 (INV-03)");
        pred.TerminusSkew.Should().BeGreaterThanOrEqualTo(0.0, "Skew_0 = 0 is always a prefix value, so max ≥ 0 (INV-03)");
        pred.IsSignificant.Should().Be(reference.MaxSkew > reference.MinSkew,
            because: "IsSignificant ⇔ max > min (INV-05)");
    }

    /// <summary>
    /// BE — "single clear minimum", lowercase / mixed case. The typed path upper-cases
    /// at construction (§3.3, INV: case-insensitive), so the origin/terminus of a
    /// lowercase sequence must equal those of its uppercase form. Pins that
    /// case-folding does not shift the predicted positions.
    /// </summary>
    [TestCase("cccggggg", "CCCGGGGG", TestName = "Replication_SingleMin_AllLower_MatchesUpper")]
    [TestCase("CcCgGgGg", "CCCGGGGG", TestName = "Replication_SingleMin_MixedCase_MatchesUpper")]
    public void Replication_MixedCase_MatchesUppercase(string lower, string upper)
    {
        var lowerPred = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(lower));
        var upperPred = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(upper));

        lowerPred.Should().Be(upperPred,
            because: "input is case-folded before the cumulative skew is computed; case must not shift the origin");
    }

    #endregion

    #region BE — Boundary: circular wrap (rotation shifts the predicted origin)

    /// <summary>
    /// BE — "circular wrap": a bacterial chromosome is CIRCULAR (§1), so the true
    /// origin can lie ACROSS the linearization join. The unit, by documented design
    /// (§6.2, ASM-02), treats its input as a LINEAR string in genome coordinates: a
    /// rotated input shifts/mirrors the predicted positions. We pin the spec-faithful
    /// consequence — ROTATING a sequence so the cumulative minimum sits at / near the
    /// END boundary still recovers the correct FIRST-minimizer index from the §2.2
    /// reference, with no over-read past the last base (no IndexOutOfRange).
    ///
    /// "GGGGCCCC" linear: 0,+1,+2,+3,+4,+3,+2,+1,0 → origin 0 (min 0), terminus 4.
    /// Rotating the falling C-run to the END, "CCCCGGGG": 0,−1,−2,−3,−4,−3,−2,−1,0 →
    /// origin 4 (the minimum now sits mid-sequence). Rotating so the minimum lands on
    /// the LAST prefix, "CCCCGGGGGGGG"… we instead test the wrap face directly: a
    /// sequence whose global minimum is the FINAL prefix index n (origin == n).
    /// All expected indices come from the §2.2 reference, not the code.
    /// </summary>
    [TestCase("GGGGCCCC", TestName = "Replication_Wrap_MaxAtMiddle_MinAtEnds")]
    [TestCase("CCCCGGGG", TestName = "Replication_Wrap_MinAtMiddle")]
    [TestCase("GGGGGCCCCCCC", TestName = "Replication_Wrap_MinOnFinalPrefix")]
    public void Replication_CircularRotation_RecoversFirstMinimizerWithNoOverRead(string input)
    {
        var reference = ReferenceOrigin(input);

        var pred = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(input));

        pred.PredictedOrigin.Should().Be(reference.Origin,
            because: $"the first minimizer of \"{input}\" is prefix index {reference.Origin} even when it sits near the end (§2.2 INV-01)");
        pred.PredictedTerminus.Should().Be(reference.Terminus,
            because: $"the first maximizer of \"{input}\" is prefix index {reference.Terminus} (§2.2 INV-02)");
        pred.PredictedOrigin.Should().BeInRange(0, input.Length,
            because: "prefix indices range over [0, n]; the scan must never over-read past the last base (INV-04)");
        pred.PredictedTerminus.Should().BeInRange(0, input.Length, "INV-04");
    }

    /// <summary>
    /// BE — "circular wrap", the rotation invariant made explicit. For ANY rotation
    /// of a genome the unit predicts the SAME biological origin position MODULO the
    /// rotation amount (ASM-02, §6.2): if the linear origin of the unrotated genome
    /// is o, then rotating left by k bases moves the cumulative minimum to (o − k)
    /// in the rotated coordinates (a circular shift of the diagram). We pin this on a
    /// fixed genome with a single clear minimum by rotating it through every offset
    /// and checking the predicted origin tracks the rotation, computed independently
    /// from the §2.2 reference on each rotated string. This pins that "circular wrap"
    /// is handled exactly as the linear definition demands — no special-casing, no
    /// crash when the minimum crosses the join.
    /// </summary>
    [Test]
    public void Replication_AllRotations_TrackTheCircularShift()
    {
        // A genome with one clear cumulative minimum; rotating it sweeps the minimum
        // through every position, including across the linearization join.
        const string genome = "GGGAGCTCCCCGGTTAACCGG";

        for (int k = 0; k < genome.Length; k++)
        {
            string rotated = genome.Substring(k) + genome.Substring(0, k);
            var reference = ReferenceOrigin(rotated);

            var pred = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(rotated));

            pred.PredictedOrigin.Should().Be(reference.Origin,
                because: $"rotation by {k} moves the cumulative minimum to its §2.2 first-minimizer index in the rotated frame");
            pred.PredictedTerminus.Should().Be(reference.Terminus,
                because: $"rotation by {k} moves the cumulative maximum to its §2.2 first-maximizer index in the rotated frame");
            pred.PredictedOrigin.Should().BeInRange(0, rotated.Length, "INV-04");
            pred.PredictedTerminus.Should().BeInRange(0, rotated.Length, "INV-04");
        }
    }

    #endregion

    #region BE — Boundary: extremely long

    /// <summary>
    /// BE/OVF: an extremely long valid sequence must compute in O(n) without hang or
    /// integer overflow, and the prediction must obey the invariants. A constructed
    /// long genome — a descending C-block then an ascending G-block — has its single
    /// global minimum at the block boundary; the exact origin/terminus indices are
    /// derived from the §2.2 reference at scale. A fixed-seed random long genome pins
    /// the structural invariants (indices in [0,n], OriginSkew ≤ 0 ≤ TerminusSkew).
    /// Locally-seeded RNG so sibling fixtures are untouched. [CancelAfter] guards hang.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void Replication_ExtremelyLong_ObeysInvariantsAndDoesNotHang()
    {
        const int block = 250_000;

        // C-block (descending to −block at prefix index `block`) then G-block (rising).
        // The single global minimum is at prefix index `block`; the maximum is the
        // final prefix index 2*block at skew 0 — wait: G-block only returns to 0, so
        // the maximum value is 0 first seen at Skew_0 → terminus 0. Derived below from
        // the reference, NOT hand-asserted, so the expectation cannot be mis-stated.
        string constructed = new string('C', block) + new string('G', block);
        var reference = ReferenceOrigin(constructed);

        ReplicationOriginPrediction pred = default;
        var act = () => pred = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(constructed));
        act.Should().NotThrow("an O(n) scalar fold must not hang or overflow at 500k bases");

        pred.PredictedOrigin.Should().Be(reference.Origin,
            because: "the single global minimum sits at the C/G block boundary (§2.2 reference)");
        pred.PredictedTerminus.Should().Be(reference.Terminus);
        pred.OriginSkew.Should().BeApproximately(reference.MinSkew, Tolerance);
        pred.TerminusSkew.Should().BeApproximately(reference.MaxSkew, Tolerance);
        pred.OriginSkew.Should().BeLessThanOrEqualTo(0.0, "INV-03");
        pred.TerminusSkew.Should().BeGreaterThanOrEqualTo(0.0, "INV-03");

        // Fixed-seed random long genome: structural invariants must hold at scale.
        var local = new Random(71193);
        const string bases = "ACGT";
        var chars = new char[2 * block];
        for (int i = 0; i < chars.Length; i++)
            chars[i] = bases[local.Next(bases.Length)];
        var randomPred = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(new string(chars)));

        randomPred.PredictedOrigin.Should().BeInRange(0, chars.Length, "INV-04");
        randomPred.PredictedTerminus.Should().BeInRange(0, chars.Length, "INV-04");
        randomPred.OriginSkew.Should().BeLessThanOrEqualTo(0.0, "INV-03");
        randomPred.TerminusSkew.Should().BeGreaterThanOrEqualTo(0.0, "INV-03");
        double.IsFinite(randomPred.OriginSkew).Should().BeTrue();
        double.IsFinite(randomPred.TerminusSkew).Should().BeTrue();
    }

    #endregion

    #region BE/RB — randomized boundary sweep (invariants + reference agreement)

    /// <summary>
    /// BE/RB: a locally-seeded (never the shared static Rng) randomized sweep over
    /// short genomes — including the empty, single-base, all-G, all-C, A/T-only and
    /// tied-extremum boundaries that arise by chance — asserting on EVERY input:
    ///   • no crash, no hang ([CancelAfter]), finite skews, no NaN/Inf;
    ///   • PredictedOrigin / PredictedTerminus equal the independent §2.2 reference's
    ///     first-minimizer / first-maximizer (the real algorithmic contract);
    ///   • OriginSkew ≤ 0 ≤ TerminusSkew (INV-03) and both indices in [0, n] (INV-04);
    ///   • IsSignificant ⇔ max > min (INV-05).
    /// Both surfaces are exercised: the strict typed path over valid DNA and the
    /// lenient raw-string path (which must, per INV-06 / §3.3, ignore non-G/C bytes).
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void Replication_RandomSweep_MatchesReferenceAndObeysInvariants()
    {
        var local = new Random(90218);
        const string bases = "ACGT";

        for (int trial = 0; trial < 400; trial++)
        {
            int n = local.Next(0, 50);
            var chars = new char[n];
            for (int i = 0; i < n; i++)
                chars[i] = bases[local.Next(bases.Length)];
            string seq = new(chars);

            var reference = ReferenceOrigin(seq);

            // (a) strict typed path over valid DNA
            ReplicationOriginPrediction typed = default;
            var typedAct = () => typed = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(seq));
            typedAct.Should().NotThrow($"valid DNA must never make the typed path throw; input: \"{seq}\"");

            typed.PredictedOrigin.Should().Be(reference.Origin,
                $"the typed origin must be the §2.2 first-minimizer; input: \"{seq}\"");
            typed.PredictedTerminus.Should().Be(reference.Terminus,
                $"the typed terminus must be the §2.2 first-maximizer; input: \"{seq}\"");
            typed.OriginSkew.Should().BeApproximately(reference.MinSkew, Tolerance, $"input: \"{seq}\"");
            typed.TerminusSkew.Should().BeApproximately(reference.MaxSkew, Tolerance, $"input: \"{seq}\"");
            typed.OriginSkew.Should().BeLessThanOrEqualTo(0.0, $"INV-03; input: \"{seq}\"");
            typed.TerminusSkew.Should().BeGreaterThanOrEqualTo(0.0, $"INV-03; input: \"{seq}\"");
            typed.PredictedOrigin.Should().BeInRange(0, n, $"INV-04; input: \"{seq}\"");
            typed.PredictedTerminus.Should().BeInRange(0, n, $"INV-04; input: \"{seq}\"");
            typed.IsSignificant.Should().Be(reference.MaxSkew > reference.MinSkew,
                $"IsSignificant ⇔ max > min (INV-05); input: \"{seq}\"");
            double.IsFinite(typed.OriginSkew).Should().BeTrue($"input: \"{seq}\"");
            double.IsFinite(typed.TerminusSkew).Should().BeTrue($"input: \"{seq}\"");

            // (b) lenient raw-string path — same input plus injected non-G/C garbage
            //     that INV-06 must ignore: it must agree with the reference computed on
            //     the same (upper-cased) string, since only G/C move the diagram.
            string rawInput = seq.Length > 0 && (trial % 3 == 0) ? seq.ToLowerInvariant() : seq;
            var rawReference = ReferenceOrigin(rawInput);
            ReplicationOriginPrediction raw = default;
            var rawAct = () => raw = GcSkewCalculator.PredictReplicationOrigin(rawInput);
            rawAct.Should().NotThrow($"the lenient overload must never throw; input: \"{rawInput}\"");
            raw.PredictedOrigin.Should().Be(rawReference.Origin, $"input: \"{rawInput}\"");
            raw.PredictedTerminus.Should().Be(rawReference.Terminus, $"input: \"{rawInput}\"");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-RNACOMP-001 — RNA-specific per-base complement
    //  Primitive: SequenceExtensions.GetRnaComplementBase(char)
    //   (IUPAC-complete; A→U emits the RNA alphabet; T treated as U → A;
    //    recognized bases uppercased; non-IUPAC passes through unchanged;
    //    TOTAL over char — never throws)
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-RNACOMP-001 — RNA complement

    #region Helpers — independent (doc-derived) RNA complement reference

    /// <summary>
    /// Independent reference for the RNA complement of a single char, transcribed
    /// directly from the Biopython <c>ambiguous_rna_complement</c> table and the
    /// §2.2/§6.1 rules of RNA_Complement.md — NOT from the switch arms under test:
    ///   • A→U, U→A, C→G, G→C (canonical RNA pairing; emits U, never T);
    ///   • T→A (T treated as U: ambiguous_rna_complement["T"]=["U"]="A");
    ///   • R↔Y, S→S, W→W, K↔M, B↔V, D↔H, N→N (IUPAC ambiguity);
    ///   • recognized symbols are UPPERCASED on output (repo convention §5.4);
    ///   • any other char passes through UNCHANGED, keeping its original case.
    /// </summary>
    private static char ExpectedRnaComplement(char nucleotide)
    {
        // Fold to uppercase only to classify; output case is fixed by the rule.
        char upper = char.ToUpperInvariant(nucleotide);
        return upper switch
        {
            'A' => 'U',
            'U' => 'A',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A', // T treated as U → A (RNA_Complement.md INV-03, §6.1)
            'R' => 'Y',
            'Y' => 'R',
            'S' => 'S',
            'W' => 'W',
            'K' => 'M',
            'M' => 'K',
            'B' => 'V',
            'V' => 'B',
            'D' => 'H',
            'H' => 'D',
            'N' => 'N',
            _ => nucleotide, // non-IUPAC: unchanged, original case
        };
    }

    /// <summary>
    /// Whole-sequence RNA complement composed per RNA_Complement.md §6.2 (the
    /// caller composes the single-char primitive over the string). Used only to
    /// exercise the BE empty boundary and the §7.1 worked example end to end.
    /// </summary>
    private static string RnaComplementString(string sequence) =>
        new(sequence.Select(SequenceExtensions.GetRnaComplementBase).ToArray());

    #endregion

    #region Positive sanity — §7.1 worked example & canonical table

    /// <summary>
    /// Positive sanity (§7.1 worked example): the RNA complement of the canonical
    /// Biopython string "ACGTUacgtuXYZxyz" must be "UGCAAUGCAAXRZxRz" under this
    /// repo's uppercase convention. Hand-derived from the table, independent of the
    /// code:
    ///   A→U C→G G→C T→A U→A         ⇒ UGCAA
    ///   a→U c→G g→C t→A u→A         ⇒ UGCAA  (lowercase folds, output uppercase)
    ///   X→X (non-IUPAC) Y→R Z→Z     ⇒ XRZ
    ///   x→x (non-IUPAC) y→R z→z     ⇒ xRz
    /// This pins ALL three checklist boundaries at once: T-as-U (T,t → A), non-RNA
    /// pass-through (X,Z,x,z), and lowercase case-folding (a,c,g,t,u,y).
    /// </summary>
    [Test]
    public void GetRnaComplementBase_WorkedExampleString_MatchesDocumentedUppercaseConvention()
    {
        const string input = "ACGTUacgtuXYZxyz";
        const string expected = "UGCAAUGCAAXRZxRz";

        RnaComplementString(input).Should().Be(expected,
            because: "RNA_Complement.md §7.1: recognized bases are RNA-complemented and uppercased " +
                     "(T/t treated as U → A); non-IUPAC X/Z/x/z pass through verbatim");
    }

    /// <summary>
    /// Positive sanity: the full canonical RNA complement table (RNA_Complement.md
    /// §4.2), pinned base by base from the Biopython dict — emitting U not T, with
    /// the reciprocal/self-complementary ambiguity codes. The crux assertions:
    ///   • A → U (NOT 'T' — this is what distinguishes RNA from the DNA sibling);
    ///   • T → A (T treated as U), never 'T' on output.
    /// </summary>
    [TestCase('A', 'U', TestName = "GetRnaComplementBase_A_IsU_NotT")]
    [TestCase('U', 'A', TestName = "GetRnaComplementBase_U_IsA")]
    [TestCase('C', 'G', TestName = "GetRnaComplementBase_C_IsG")]
    [TestCase('G', 'C', TestName = "GetRnaComplementBase_G_IsC")]
    [TestCase('R', 'Y', TestName = "GetRnaComplementBase_R_IsY")]
    [TestCase('Y', 'R', TestName = "GetRnaComplementBase_Y_IsR")]
    [TestCase('S', 'S', TestName = "GetRnaComplementBase_S_IsS_SelfComplement")]
    [TestCase('W', 'W', TestName = "GetRnaComplementBase_W_IsW_SelfComplement")]
    [TestCase('K', 'M', TestName = "GetRnaComplementBase_K_IsM")]
    [TestCase('M', 'K', TestName = "GetRnaComplementBase_M_IsK")]
    [TestCase('B', 'V', TestName = "GetRnaComplementBase_B_IsV")]
    [TestCase('V', 'B', TestName = "GetRnaComplementBase_V_IsB")]
    [TestCase('D', 'H', TestName = "GetRnaComplementBase_D_IsH")]
    [TestCase('H', 'D', TestName = "GetRnaComplementBase_H_IsD")]
    [TestCase('N', 'N', TestName = "GetRnaComplementBase_N_IsN_SelfComplement")]
    public void GetRnaComplementBase_CanonicalTable_MatchesBiopythonRnaDict(char input, char expected)
    {
        SequenceExtensions.GetRnaComplementBase(input).Should().Be(expected,
            because: "the RNA complement is the Biopython ambiguous_rna_complement lookup (RNA_Complement.md §4.2)");
    }

    /// <summary>
    /// Positive sanity (INV-02): NO recognized base maps to 'T' — the whole point of
    /// the RNA alphabet is it emits U, never T. We assert this over the entire
    /// recognized set so a DNA-style A→T regression cannot pass.
    /// </summary>
    [Test]
    public void GetRnaComplementBase_RecognizedBases_NeverEmitT()
    {
        foreach (char b in "AUCGTRYSWKMBVDHN")
            SequenceExtensions.GetRnaComplementBase(b).Should().NotBe('T',
                because: $"INV-02: the RNA alphabet emits U, never T (input '{b}')");
    }

    /// <summary>
    /// Positive sanity (INV-06): the complement is an INVOLUTION on the recognized
    /// RNA bases and the eleven ambiguity codes — complementing twice returns the
    /// (uppercased) original. T is special: T→A→U (T is treated as U, so the round
    /// trip lands on U, the canonical form), which we assert explicitly so the
    /// involution claim is honest about the T-as-U folding.
    /// </summary>
    [Test]
    public void GetRnaComplementBase_RecognizedSymbols_AreInvolution()
    {
        // The canonical RNA symbols (excluding T, which folds to U) round-trip exactly.
        foreach (char b in "AUCGRYSWKMBVDHN")
        {
            char once = SequenceExtensions.GetRnaComplementBase(b);
            SequenceExtensions.GetRnaComplementBase(once).Should().Be(b,
                because: $"INV-06: complement is an involution on the RNA alphabet (input '{b}')");
        }

        // T is treated as U: T → A → U (not back to T), the documented folding.
        char t1 = SequenceExtensions.GetRnaComplementBase('T');
        t1.Should().Be('A');
        SequenceExtensions.GetRnaComplementBase(t1).Should().Be('U',
            because: "INV-03/§6.1: T is treated as U, so T→A→U lands on the canonical U, not T");
    }

    #endregion

    #region MC — "T instead of U": T/t treated as U (→ A)

    /// <summary>
    /// MC (checklist crux "T instead of U"): a thymine in an RNA context is TREATED
    /// AS A URACIL — T and t both complement to 'A' (RNA_Complement.md §2.2/§6.1,
    /// INV-03; Biopython ambiguous_rna_complement["T"]=["U"]="A"). This is NOT the
    /// DNA behaviour (where A→T): here T is folded to U and complemented as U. We
    /// pin BOTH cases and uppercase output. A regression that rejected T, passed it
    /// through unchanged, or mapped it DNA-style would be caught here.
    /// </summary>
    [TestCase('T', TestName = "GetRnaComplementBase_UppercaseT_TreatedAsU_IsA")]
    [TestCase('t', TestName = "GetRnaComplementBase_LowercaseT_TreatedAsU_IsA")]
    public void GetRnaComplementBase_TIsTreatedAsU_ComplementsToA(char thymine)
    {
        SequenceExtensions.GetRnaComplementBase(thymine).Should().Be('A',
            because: "RNA_Complement.md INV-03/§6.1: a T (or t) is treated as a U, whose complement is A");
    }

    /// <summary>
    /// MC: the RNA-vs-DNA divergence at the 'A' base — the DNA sibling
    /// GetComplementBase(SEQ-COMP-001) maps A→T, but the RNA primitive maps A→U.
    /// Pinning both side by side guards that the two surfaces cannot drift into one
    /// another (RNA_Complement.md §6.1 "A (RNA vs DNA)" row).
    /// </summary>
    [Test]
    public void GetRnaComplementBase_A_DivergesFromDnaSibling_EmitsUNotT()
    {
        SequenceExtensions.GetRnaComplementBase('A').Should().Be('U',
            because: "RNA emits U for A");
        SequenceExtensions.GetComplementBase('A').Should().Be('T',
            because: "the DNA sibling emits T for A — the surfaces must stay distinct");
    }

    #endregion

    #region MC — "lowercase": case-insensitive input, uppercase output

    /// <summary>
    /// MC (checklist "lowercase"): input is CASE-INSENSITIVE and every RECOGNIZED
    /// symbol is normalized to UPPERCASE on output (RNA_Complement.md §3.3, §5.4,
    /// §6.1). So lowercase a/u/c/g/t and the lowercase ambiguity codes complement
    /// IDENTICALLY to their uppercase forms and ALWAYS emit an uppercase letter.
    /// This pins both case-folding (lower accepted) and the uppercase convention.
    /// </summary>
    [TestCase('a', 'U', TestName = "GetRnaComplementBase_LowerA_IsUpperU")]
    [TestCase('u', 'A', TestName = "GetRnaComplementBase_LowerU_IsUpperA")]
    [TestCase('c', 'G', TestName = "GetRnaComplementBase_LowerC_IsUpperG")]
    [TestCase('g', 'C', TestName = "GetRnaComplementBase_LowerG_IsUpperC")]
    [TestCase('r', 'Y', TestName = "GetRnaComplementBase_LowerR_IsUpperY")]
    [TestCase('n', 'N', TestName = "GetRnaComplementBase_LowerN_IsUpperN")]
    public void GetRnaComplementBase_LowercaseRecognized_FoldsAndEmitsUppercase(char input, char expected)
    {
        char result = SequenceExtensions.GetRnaComplementBase(input);

        result.Should().Be(expected,
            because: "RNA_Complement.md §3.3/§5.4: recognized lowercase symbols fold and emit uppercase complement");
        char.IsUpper(result).Should().BeTrue(
            because: "recognized RNA complements are always returned uppercase");
    }

    #endregion

    #region MC / BE — "non-RNA base": out-of-alphabet pass-through, total over char

    /// <summary>
    /// MC (checklist crux "non-RNA base"): characters OUTSIDE the IUPAC nucleotide
    /// set must PASS THROUGH UNCHANGED, keeping their ORIGINAL case, and the method
    /// must NEVER throw (RNA_Complement.md §3.3, §6.1 — a total function over char).
    /// Covers gaps, digits, the non-IUPAC letters X/Z (lower and upper), whitespace,
    /// the null byte, a unicode letter, and lone surrogate halves. Note 'X' and 'Z'
    /// are NOT in this repo's recognized RNA set, so they pass through verbatim
    /// (matching the §7.1 worked example where X/Z/x/z are untouched).
    /// </summary>
    [TestCase('-', '-', TestName = "GetRnaComplementBase_GapDash_PassesThrough")]
    [TestCase('.', '.', TestName = "GetRnaComplementBase_GapDot_PassesThrough")]
    [TestCase('5', '5', TestName = "GetRnaComplementBase_Digit_PassesThrough")]
    [TestCase('X', 'X', TestName = "GetRnaComplementBase_UpperX_PassesThrough")]
    [TestCase('x', 'x', TestName = "GetRnaComplementBase_LowerX_PassesThroughKeepsCase")]
    [TestCase('Z', 'Z', TestName = "GetRnaComplementBase_UpperZ_PassesThrough")]
    [TestCase('z', 'z', TestName = "GetRnaComplementBase_LowerZ_PassesThroughKeepsCase")]
    [TestCase(' ', ' ', TestName = "GetRnaComplementBase_Whitespace_PassesThrough")]
    [TestCase('\0', '\0', TestName = "GetRnaComplementBase_NullByte_PassesThrough")]
    [TestCase('α', 'α', TestName = "GetRnaComplementBase_GreekLetter_PassesThrough")]
    [TestCase('\uD83D', '\uD83D', TestName = "GetRnaComplementBase_HighSurrogateHalf_PassesThrough")]
    [TestCase('\uDE00', '\uDE00', TestName = "GetRnaComplementBase_LowSurrogateHalf_PassesThrough")]
    public void GetRnaComplementBase_NonRnaCharacter_PassesThroughUnchangedAndNeverThrows(char input, char expected)
    {
        char result = '￿';
        var act = () => result = SequenceExtensions.GetRnaComplementBase(input);

        act.Should().NotThrow("the RNA complement primitive is total over char and never throws");
        result.Should().Be(expected,
            because: "non-IUPAC characters pass through unchanged, keeping their original case (RNA_Complement.md §3.3/§6.1)");
    }

    /// <summary>
    /// BE (empty boundary): there is no empty char, so the BE empty boundary is
    /// honoured on the §6.2 caller-composed whole-sequence complement — the RNA
    /// complement of the empty string is the empty string, with no throw and no hang.
    /// </summary>
    [Test]
    public void GetRnaComplementBase_EmptySequenceComposition_IsEmptyAndDoesNotThrow()
    {
        var act = () => RnaComplementString(string.Empty).Should().BeEmpty(
            because: "the RNA complement of the empty sequence is the empty sequence (RNA_Complement.md §6.2)");

        act.Should().NotThrow("the empty sequence is a defined boundary, not an error");
    }

    #endregion

    #region MC / BE — randomized boundary sweep (total, matches doc-derived reference)

    /// <summary>
    /// MC/BE: a LOCALLY-seeded (never the shared static Rng) randomized sweep over
    /// arbitrary BMP code points — deliberately spanning recognized RNA bases, T/t,
    /// non-IUPAC garbage, control chars, the null byte and lone surrogate halves —
    /// asserting on EVERY char:
    ///   • the primitive NEVER throws (total function over char, RNA_Complement.md §3.3);
    ///   • it EXACTLY equals the independent doc-derived ExpectedRnaComplement (the
    ///     real algorithmic contract: A→U, T→A as U, uppercase recognized, non-IUPAC
    ///     pass-through). A code that drifted from the spec would diverge here.
    /// Mixed in are the explicit boundary chars so the sweep cannot miss them by
    /// chance. [CancelAfter] guards against any pathological hang.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void GetRnaComplementBase_RandomCharSweep_NeverThrowsAndMatchesDocReference()
    {
        var local = new Random(2120601);
        // Boundary chars guaranteed in the sweep regardless of RNG draws.
        char[] forced =
        {
            'A', 'U', 'C', 'G', 'T', 't', 'a', 'u', 'c', 'g',
            'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'V', 'D', 'H', 'N',
            'r', 'y', 's', 'w', 'k', 'm', 'b', 'v', 'd', 'h', 'n',
            'X', 'Z', 'x', 'z', '-', '.', ' ', '\0', '5', 'α',
            '\uD83D', '\uDE00', '￿',
        };

        foreach (char c in forced)
            AssertMatchesReference(c);

        for (int trial = 0; trial < 5000; trial++)
            AssertMatchesReference((char)local.Next(0x0000, 0x10000));

        return;

        static void AssertMatchesReference(char c)
        {
            char result = '￿';
            var act = () => result = SequenceExtensions.GetRnaComplementBase(c);
            act.Should().NotThrow($"the RNA complement primitive is total over char; input U+{(int)c:X4}");
            result.Should().Be(ExpectedRnaComplement(c),
                $"the RNA complement must equal the doc-derived reference; input U+{(int)c:X4}");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-GC-ANALYSIS-001 — comprehensive GC analysis
    //  Typed entry:  GcSkewCalculator.AnalyzeGcContent(DnaSequence,…) — strict ctor
    //  Lenient entry: GcSkewCalculator.AnalyzeGcContent(string,…) — counts only ACGT
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-GC-ANALYSIS-001 — comprehensive GC analysis

    #region Helpers — independent doc-derived reference

    /// <summary>
    /// Independent reference for the doc-derived overall scalars
    /// (Comprehensive_GC_Analysis.md §2.2). Counts ONLY A/C/G/T, case-insensitively,
    /// and applies the documented zero-denominator → 0 convention. This is computed
    /// from FIRST PRINCIPLES, NOT from the code under test, so a test comparing the
    /// algorithm output against it cannot rubber-stamp a wrong implementation.
    /// </summary>
    private static (double gcContent, double gcSkew, double atSkew) ExpectedOverallScalars(string input)
    {
        int g = 0, c = 0, a = 0, t = 0;
        foreach (char ch in input)
        {
            switch (char.ToUpperInvariant(ch))
            {
                case 'G': g++; break;
                case 'C': c++; break;
                case 'A': a++; break;
                case 'T': t++; break;
            }
        }

        int counted = g + c + a + t;
        double gcContent = counted > 0 ? (double)(g + c) / counted * 100.0 : 0.0;
        double gcSkew = (g + c) > 0 ? (double)(g - c) / (g + c) : 0.0;
        double atSkew = (a + t) > 0 ? (double)(a - t) / (a + t) : 0.0;
        return (gcContent, gcSkew, atSkew);
    }

    /// <summary>Asserts every numeric scalar of a result is a real, finite number (no NaN/Inf).</summary>
    private static void AssertAllScalarsFinite(GcAnalysisResult r)
    {
        double.IsNaN(r.OverallGcContent).Should().BeFalse("OverallGcContent must never be NaN");
        double.IsNaN(r.OverallGcSkew).Should().BeFalse("OverallGcSkew must never be NaN");
        double.IsNaN(r.OverallAtSkew).Should().BeFalse("OverallAtSkew must never be NaN");
        double.IsNaN(r.GcContentVariance).Should().BeFalse("GcContentVariance must never be NaN");
        double.IsNaN(r.GcSkewVariance).Should().BeFalse("GcSkewVariance must never be NaN");
        double.IsInfinity(r.OverallGcContent).Should().BeFalse("OverallGcContent must never be Infinity");
        double.IsInfinity(r.OverallGcSkew).Should().BeFalse("OverallGcSkew must never be Infinity");
        double.IsInfinity(r.OverallAtSkew).Should().BeFalse("OverallAtSkew must never be Infinity");
        double.IsInfinity(r.GcContentVariance).Should().BeFalse("GcContentVariance must never be Infinity");
        double.IsInfinity(r.GcSkewVariance).Should().BeFalse("GcSkewVariance must never be Infinity");
    }

    /// <summary>Asserts the documented invariant ranges (INV-01..05) of a result.</summary>
    private static void AssertInvariantRanges(GcAnalysisResult r)
    {
        r.OverallGcContent.Should().BeInRange(0.0, 100.0, "INV-03: GC% ∈ [0,100]");
        r.OverallGcSkew.Should().BeInRange(-1.0, 1.0, "INV-01: GC skew ∈ [−1,+1]");
        r.OverallAtSkew.Should().BeInRange(-1.0, 1.0, "INV-02: AT skew ∈ [−1,+1]");
        r.GcContentVariance.Should().BeGreaterThanOrEqualTo(0.0, "INV-04: a variance is ≥ 0");
        r.GcSkewVariance.Should().BeGreaterThanOrEqualTo(0.0, "INV-04: a variance is ≥ 0");
    }

    #endregion

    #region Positive sanity — worked examples (Comprehensive_GC_Analysis.md §7.1)

    /// <summary>
    /// Positive sanity: the §7.1 worked example "GGGCCAT" (G=3,C=2,A=1,T=1,n=7).
    /// Every overall scalar is hand-derived from §2.2 and pinned EXACTLY:
    /// GcContent = 5/7×100 = 71.42857142857143, GcSkew = (3−2)/5 = 0.2,
    /// AtSkew = (1−1)/2 = 0.0. Reproduced from the doc, not the code.
    /// </summary>
    [Test]
    public void AnalyzeGcContent_WorkedExample_Gggccat_MatchesDocScalars()
    {
        var result = GcSkewCalculator.AnalyzeGcContent(new DnaSequence("GGGCCAT"));

        result.SequenceLength.Should().Be(7);
        result.OverallGcContent.Should().BeApproximately(5.0 / 7.0 * 100.0, Tolerance,
            because: "Gc% = (G+C)/n×100 = 5/7×100 (Comprehensive_GC_Analysis.md §7.1)");
        result.OverallGcSkew.Should().BeApproximately(0.2, Tolerance,
            because: "GC skew = (G−C)/(G+C) = (3−2)/5 = 0.2 (§7.1)");
        result.OverallAtSkew.Should().BeApproximately(0.0, Tolerance,
            because: "AT skew = (A−T)/(A+T) = (1−1)/2 = 0 (§7.1)");
    }

    /// <summary>
    /// Positive sanity: the §7.1 variance worked example. "GGCC", window 2, step 2
    /// yields exactly two full windows — GG (skew +1, GC% 100) and CC (skew −1,
    /// GC% 100). The POPULATION variance (÷N) of the per-window skews is
    /// ((1−0)²+(−1−0)²)/2 = 1.0; of the per-window GC% it is 0.0 (both windows 100).
    /// This specifically pins the ÷N population variance: the sample variance (÷N−1)
    /// would give 2.0 for the skew, so a code that used Bessel's correction FAILS.
    /// </summary>
    [Test]
    public void AnalyzeGcContent_WorkedExample_Ggcc_PopulationVarianceMatchesDoc()
    {
        var result = GcSkewCalculator.AnalyzeGcContent(new DnaSequence("GGCC"), windowSize: 2, stepSize: 2);

        result.WindowedGcSkew.Should().HaveCount(2, "INV-05: ⌊(4−2)/2⌋+1 = 2 full windows");
        result.WindowedGcContent.Should().HaveCount(2);
        result.WindowedGcSkew.Select(p => p.GcSkew).Should().Equal(new[] { 1.0, -1.0 },
            because: "GG → (2−0)/2 = +1, CC → (0−2)/2 = −1");
        result.GcSkewVariance.Should().BeApproximately(1.0, Tolerance,
            because: "population variance ((1−0)²+(−1−0)²)/2 = 1.0; sample variance ÷N−1 would be 2.0 (§7.1)");
        result.GcContentVariance.Should().BeApproximately(0.0, Tolerance,
            because: "both windows are 100% GC, so their variance is exactly 0 (§7.1)");
    }

    #endregion

    #region BE — Boundary: empty input (zero result, no NaN / no divide-by-zero)

    /// <summary>
    /// BE: the empty input is the lower size boundary and the zero-denominator
    /// boundary all at once. Per §3.3 / §6.1 it must produce the DEFINED ZERO RESULT
    /// — all scalars 0, both windowed lists empty, SequenceLength 0 — NOT a
    /// DivideByZero crash and NOT a NaN. Covered on BOTH surfaces: the typed
    /// DnaSequence("") path and the lenient string ("" and null) path.
    /// </summary>
    [Test]
    public void AnalyzeGcContent_EmptyInput_IsZeroResultAndDoesNotThrow()
    {
        var act = () =>
        {
            var typed = GcSkewCalculator.AnalyzeGcContent(new DnaSequence(string.Empty));
            var emptyStr = GcSkewCalculator.AnalyzeGcContent(string.Empty);
            var nullStr = GcSkewCalculator.AnalyzeGcContent((string)null!);

            foreach (var r in new[] { typed, emptyStr, nullStr })
            {
                r.SequenceLength.Should().Be(0, "empty input has length 0");
                r.OverallGcContent.Should().Be(0.0, "no bases → GC% 0 (numerator 0)");
                r.OverallGcSkew.Should().Be(0.0, "G+C=0 → skew 0 (zero-division guard)");
                r.OverallAtSkew.Should().Be(0.0, "A+T=0 → skew 0 (zero-division guard)");
                r.GcContentVariance.Should().Be(0.0, "no windows → variance 0");
                r.GcSkewVariance.Should().Be(0.0, "no windows → variance 0");
                r.WindowedGcSkew.Should().BeEmpty("no full window fits an empty sequence");
                r.WindowedGcContent.Should().BeEmpty();
                AssertAllScalarsFinite(r);
            }
        };

        act.Should().NotThrow("empty input is a defined zero-result boundary, never a crash or NaN");
    }

    /// <summary>
    /// BE: a null DnaSequence is the documented ArgumentNullException boundary on the
    /// typed surface (§6.1, contract parity) — an INTENTIONAL validation throw, never
    /// a raw NullReferenceException.
    /// </summary>
    [Test]
    public void AnalyzeGcContent_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => GcSkewCalculator.AnalyzeGcContent((DnaSequence)null!);

        act.Should().Throw<ArgumentNullException>(
            "a null DnaSequence is rejected by the documented explicit guard, never dereferenced");
    }

    #endregion

    #region BE — Boundary: all-GC (GC=100, skew bounds)

    /// <summary>
    /// BE: a pure-G sequence (§6.1) is the GC-content maximum and the GC-skew maximum:
    /// GcContent = 100, GcSkew = +1; a pure-C sequence gives GcContent = 100,
    /// GcSkew = −1. With no A/T present the AT-skew zero-denominator guard returns the
    /// defined 0, not NaN. All scalars stay inside the invariant ranges.
    /// </summary>
    [TestCase("GGGGGGGG", 100.0, 1.0, TestName = "AnalyzeGcContent_AllG_Gc100_Skew+1")]
    [TestCase("CCCCCCCC", 100.0, -1.0, TestName = "AnalyzeGcContent_AllC_Gc100_Skew-1")]
    [TestCase("GCGCGCGC", 100.0, 0.0, TestName = "AnalyzeGcContent_AlternatingGC_Gc100_Skew0")]
    public void AnalyzeGcContent_AllGc_GcIs100(string input, double expectedGc, double expectedSkew)
    {
        var result = GcSkewCalculator.AnalyzeGcContent(new DnaSequence(input));

        result.OverallGcContent.Should().BeApproximately(expectedGc, Tolerance,
            because: "a sequence of only G/C is 100% GC");
        result.OverallGcSkew.Should().BeApproximately(expectedSkew, Tolerance,
            because: "GC skew = (G−C)/(G+C) for a pure / balanced G·C run");
        result.OverallAtSkew.Should().Be(0.0,
            because: "no A and no T present → AT-skew zero-denominator guard returns 0, not NaN");
        AssertInvariantRanges(result);
        AssertAllScalarsFinite(result);
    }

    #endregion

    #region BE — Boundary: all-AT (GC=0)

    /// <summary>
    /// BE: a pure-A / pure-T / balanced-AT sequence is the GC-content minimum:
    /// GcContent = 0 (numerator 0). With no G/C present the GC-skew zero-denominator
    /// guard returns the defined 0 (not NaN). The AT skew follows (A−T)/(A+T):
    /// all-A → +1, all-T → −1, balanced AT → 0.
    /// </summary>
    [TestCase("AAAAAAAA", 1.0, TestName = "AnalyzeGcContent_AllA_Gc0_AtSkew+1")]
    [TestCase("TTTTTTTT", -1.0, TestName = "AnalyzeGcContent_AllT_Gc0_AtSkew-1")]
    [TestCase("ATATATAT", 0.0, TestName = "AnalyzeGcContent_AlternatingAT_Gc0_AtSkew0")]
    public void AnalyzeGcContent_AllAt_GcIs0(string input, double expectedAtSkew)
    {
        var result = GcSkewCalculator.AnalyzeGcContent(new DnaSequence(input));

        result.OverallGcContent.Should().Be(0.0,
            because: "a sequence of only A/T has no G or C → GC% is exactly 0");
        result.OverallGcSkew.Should().Be(0.0,
            because: "no G and no C present → GC-skew zero-denominator guard returns 0, not NaN");
        result.OverallAtSkew.Should().BeApproximately(expectedAtSkew, Tolerance,
            because: "AT skew = (A−T)/(A+T) for a pure / balanced A·T run");
        AssertInvariantRanges(result);
        AssertAllScalarsFinite(result);
    }

    #endregion

    #region BE — Boundary: non-ACGT characters (ignored in num & denom)

    /// <summary>
    /// BE: non-A/C/G/T characters are IGNORED in BOTH numerator and denominator on
    /// the lenient string surface (§3.3, §6.2; Biopython GC_skew ignores ambiguous
    /// bases). So injecting N/ambiguity/U/digits/gaps/whitespace/null-byte/unicode
    /// among a known A/C/G/T core must NEVER throw and must leave every overall
    /// scalar EXACTLY equal to the scalar computed on the counted bases alone.
    /// Here the counted core of each input is "GC" → GcContent 100, GcSkew 0
    /// (G=C=1), AtSkew 0 (no A/T). The garbage must not shift any of them.
    /// </summary>
    [TestCase("GNCN", TestName = "AnalyzeGcContent_NonAcgt_AmbiguityN_Ignored")]
    [TestCase("G-C.", TestName = "AnalyzeGcContent_NonAcgt_Gaps_Ignored")]
    [TestCase("G C ", TestName = "AnalyzeGcContent_NonAcgt_Whitespace_Ignored")]
    [TestCase("G1C2", TestName = "AnalyzeGcContent_NonAcgt_Digits_Ignored")]
    [TestCase("GUCU", TestName = "AnalyzeGcContent_NonAcgt_RnaBaseU_Ignored")]
    [TestCase("G\0C", TestName = "AnalyzeGcContent_NonAcgt_NullByte_Ignored")]
    [TestCase("GαCβ", TestName = "AnalyzeGcContent_NonAcgt_Unicode_Ignored")]
    public void AnalyzeGcContent_NonAcgtCharacters_AreIgnoredNotCounted(string input)
    {
        GcAnalysisResult result = default!;
        var act = () => result = GcSkewCalculator.AnalyzeGcContent(input);

        act.Should().NotThrow("non-ACGT symbols are ignored, not validated, on the lenient string surface");

        var (gc, gcSkew, atSkew) = ExpectedOverallScalars(input);
        gc.Should().BeApproximately(100.0, Tolerance, "sanity: the counted core 'GC' is 100% GC");
        gcSkew.Should().BeApproximately(0.0, Tolerance, "sanity: counted core G=C → skew 0");

        result.OverallGcContent.Should().BeApproximately(gc, Tolerance,
            because: "ignored garbage must not shift GC% of the counted A/C/G/T bases");
        result.OverallGcSkew.Should().BeApproximately(gcSkew, Tolerance,
            because: "ignored garbage must not shift the GC skew of the counted bases");
        result.OverallAtSkew.Should().BeApproximately(atSkew, Tolerance,
            because: "ignored garbage must not shift the AT skew of the counted bases");
        AssertAllScalarsFinite(result);
        AssertInvariantRanges(result);
    }

    /// <summary>
    /// BE: an input that is ENTIRELY non-ACGT (no counted base at all) is the full
    /// zero-denominator boundary on the lenient surface: every scalar is the defined
    /// 0, no NaN, no crash — yet SequenceLength still reflects the raw input length
    /// (§3.3: only the METRICS ignore the symbols; the length is the analyzed length).
    /// </summary>
    [TestCase("NNNNN", TestName = "AnalyzeGcContent_AllNonAcgt_AllN")]
    [TestCase("-----", TestName = "AnalyzeGcContent_AllNonAcgt_AllGaps")]
    [TestCase("12345", TestName = "AnalyzeGcContent_AllNonAcgt_AllDigits")]
    public void AnalyzeGcContent_AllNonAcgt_IsZeroScalarsNoNaN(string input)
    {
        GcAnalysisResult result = default!;
        var act = () => result = GcSkewCalculator.AnalyzeGcContent(input);

        act.Should().NotThrow("an all-garbage input is ignored base-by-base, never a crash");
        result.OverallGcContent.Should().Be(0.0, "no counted base → GC% 0");
        result.OverallGcSkew.Should().Be(0.0, "no G/C → zero-division guard → 0, not NaN");
        result.OverallAtSkew.Should().Be(0.0, "no A/T → zero-division guard → 0, not NaN");
        result.SequenceLength.Should().Be(input.Length,
            because: "the analyzed length is the raw input length; only the metrics ignore non-ACGT");
        AssertAllScalarsFinite(result);
    }

    #endregion

    #region BE / OVF — Boundary: very long (O(n) scalars, no overflow / hang)

    /// <summary>
    /// BE/OVF: an extremely long valid sequence (1,000,000 bases) must compute the
    /// overall scalars without overflow, hang, or NaN, with every scalar inside its
    /// invariant range. A known-composition long input ("AG" repeated = 50% A,
    /// 50% G) pins exact values: GcContent 50, GcSkew +1 (all G/C is G), AtSkew +1
    /// (all A/T is A). [CancelAfter] guards the O(n) scalar path against any hang.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void AnalyzeGcContent_ExtremelyLong_StaysInRangeAndDoesNotHang()
    {
        const int length = 1_000_000;

        var halfGc = string.Concat(Enumerable.Repeat("AG", length / 2));
        var known = GcSkewCalculator.AnalyzeGcContent(new DnaSequence(halfGc));
        known.SequenceLength.Should().Be(length);
        known.OverallGcContent.Should().BeApproximately(50.0, Tolerance,
            because: "exactly half the bases are G → 50% GC at any length");
        known.OverallGcSkew.Should().BeApproximately(1.0, Tolerance,
            because: "the only G/C base is G → GC skew (G−C)/(G+C) = +1");
        known.OverallAtSkew.Should().BeApproximately(1.0, Tolerance,
            because: "the only A/T base is A → AT skew (A−T)/(A+T) = +1");
        AssertAllScalarsFinite(known);
        AssertInvariantRanges(known);

        var random = GcSkewCalculator.AnalyzeGcContent(new DnaSequence(RandomDna(length)));
        AssertAllScalarsFinite(random);
        AssertInvariantRanges(random);
    }

    #endregion

    #region BE — randomized boundary sweep (locally-seeded, never the shared Rng)

    /// <summary>
    /// BE: a LOCALLY-seeded (never the shared static Rng) randomized sweep. For each
    /// trial a random-length sequence is built either from pure {A,C,G,T} or from a
    /// soup of {A,C,G,T,N,-,U,space,digit} (forcing the non-ACGT ignore path), the
    /// analysis is run on the lenient string surface, and we assert on EVERY result:
    ///   • it NEVER throws and produces no NaN/Infinity in any scalar;
    ///   • every scalar lands in its documented invariant range (INV-01..04);
    ///   • the overall scalars EXACTLY equal the independent doc-derived reference
    ///     (ExpectedOverallScalars) — the real algorithmic contract, so a code that
    ///     miscounted, dropped the zero-division guard, or counted garbage diverges.
    /// [CancelAfter] guards against any pathological hang.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void AnalyzeGcContent_RandomBoundarySweep_FiniteInRangeAndMatchesReference()
    {
        var local = new Random(2330601);
        const string pure = "ACGT";
        const string soup = "ACGTN-U 5";

        for (int trial = 0; trial < 3000; trial++)
        {
            bool useSoup = local.Next(2) == 0;
            string alphabet = useSoup ? soup : pure;
            int len = local.Next(0, 40);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = alphabet[local.Next(alphabet.Length)];
            string input = new string(chars);

            GcAnalysisResult result = default!;
            var act = () => result = GcSkewCalculator.AnalyzeGcContent(input);
            act.Should().NotThrow($"the analysis must never throw on fuzz input \"{input}\"");

            AssertAllScalarsFinite(result);
            AssertInvariantRanges(result);

            var (gc, gcSkew, atSkew) = ExpectedOverallScalars(input);
            result.OverallGcContent.Should().BeApproximately(gc, Tolerance,
                because: $"GC% must match the doc-derived reference for \"{input}\"");
            result.OverallGcSkew.Should().BeApproximately(gcSkew, Tolerance,
                because: $"GC skew must match the doc-derived reference for \"{input}\"");
            result.OverallAtSkew.Should().BeApproximately(atSkew, Tolerance,
                because: $"AT skew must match the doc-derived reference for \"{input}\"");
        }
    }

    #endregion

    #endregion
}
