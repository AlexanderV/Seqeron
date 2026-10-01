# Validation Report: REP-DIRECT-001 — Direct Repeat Detection

- **Validated:** 2026-09-30 (review campaign 2026-09, batch B04: F11 + completeness audit WP2); first pass 2026-06-24 superseded
- **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindDirectRepeats(DnaSequence|string, minLength=5, maxLength=50, minSpacing=1)` — `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs:3660/3678` (core `FindDirectRepeatsCore` :3704, shared engine `EnumerateMaximalPairs` :3765; variants `FindReverseComplementRepeats` :3943, `FindApproximateDirectRepeats` :4108, `FindDegenerateRepeats` :4377, `FindSupermaximalRepeats` :4737)
- **Variants (audit WP2):** `FindReverseComplementRepeats` (:2432/2450), `FindApproximateDirectRepeats` (:2578/2596), `FindSupermaximalRepeats` (:2792/2803)
- **Variant (audit WP8, 2026-10-01):** `FindDegenerateRepeats(DnaSequence|string, minLength, maxDifferences, ApproximateRepeatDistance, reverseComplement, maxLength, minSpacing)` — Vmatch `-e k` / `-p -e k` / `-p -h k` / `-h k -allmax`
- **Stage A verdict:** FAIL → corrected (the 2026-06 description accepted "every (i, j, len) window"; the reporting convention is now sourced: maximal repeated pairs)
- **Stage B verdict:** FAIL → fixed (B04 F11); variants added and reference-verified (B04 audit WP2)
- **State:** FIXED

## Stage A — Description

### Sources opened
- MUMmer 4 `src/tigr/repeat-match.cc` (raw GitHub; compiled and run): "identifies maximal exact repeat regions"; a pair is
  reported only if left-maximal (`Data[i−1] ≠ Data[j−1]`) and right-maximal (full common prefix); `-f` = forward strand
  only; without `-f` reverse-complement pairs are printed with `r` (Start2 = 1-based last base of copy 2, kept when `k ≥ i`).
  `src/essaMEM/mummer.cpp` `-n`: "match only the characters a, c, g, or t".
- Gusfield 1997 §7.12 (maximal pairs / maximal repeats), §7.12.1 + Theorem 7.12.4 (supermaximal repeats);
  Abouelhoda, Kurtz & Ohlebusch 2004 (enhanced suffix arrays; lcp-interval traversal).
- Kurtz & Schleiermacher 1999 and Kurtz et al. 2001 (REPuter: forward / palindromic, exact / k-mismatch repeats; seeds
  of length ⌊ℓ/(k+1)⌋ extended by the maximum-error strategy); Vmatch 2.3.1 manual (Kurtz, ISC; `virtman.tex`
  Appendix A: palindromic match with `i ≤ j`, k-mismatch match `d_H ≤ k`, maximal = not contained, supermaximal repeat;
  wildcards always mismatch) and the Vmatch binaries.
- Wikipedia "Direct repeat" / "Repeated sequence (DNA)" (orientation terminology; unchanged from the first pass).

### Definition / conventions confirmed
- Direct repeat = same-orientation exact copies; reported as **maximal repeated pairs** `(i, j, L)`, `i < j`, 0-based,
  `Spacing = j − i − L`; `minLength ≤ L ≤ maxLength` and `Spacing ≥ minSpacing` are filters on maximal pairs (no
  truncation into sub-windows). Only A/C/G/T match; case-insensitive.
- Reverse-complement pairs: `S[i..i+L) = revcomp(S[k..k+L))`, `i ≤ k`, maximal outward and inward; forward-strand
  starts reported (repeat-match Start2 = k + L).
- k-mismatch repeats: Hamming distance ≤ k, maximal on the diagonal (default; k = 0 ≡ maximal pairs) or, optionally,
  not contained in any k-mismatch repeat on another diagonal (Vmatch's literal Appendix A reading, = `vmatch -h k -allmax`).
- Supermaximal repeat: maximal repeat not a substring of another maximal repeat; one record per string with all occurrences.

### Findings
- The first-pass description ("every window of every length", hand-counted C1 = 14 hits for `AAAAAATTTTAAAAAA`) did not
  match any reference tool: the maximal pairs are (0,10,6) (0,11,5) (0,12,4) (1,10,5) (2,10,4) (5 pairs). Corrected.

## Stage B — Implementation

### Code path
ACGT → 0..3, any other symbol a unique code (never matches); suffix array (prefix doubling) + Kasai LCP — shared
`SequenceComplexity.BuildSuffixArray` / `BuildLcpArray`; bottom-up lcp-interval traversal with per-(strand,)
left-character lists emitting exactly the maximal pairs (O(n log² n + z)); sorted output; eager validation.
Variants: the same engine on `S · # · revcomp(S)` (cross-strand pairs only) for reverse-complement pairs; exact maximal
pairs of length ⌊m/(k+1)⌋ as seeds + per-seed windows from the first k+1 mismatches each side (+ optional
cross-diagonal containment filter by binary search) for k-mismatch repeats; LCP local maxima with left-diverse
suffixes for supermaximal repeats.

### Defects found and fixed (B04 F11)
1. Every nested `(i, j, len)` window reported — O(L²) hits per repeat (`AAAAAATTTTAAAAAA`, 4–6 → 14 instead of 5;
   snapshot `ACGTACGTTTTTTTTTACGTACGT` → 20 instead of 4).
2. Negative `minSpacing` produced self-pairs (i, i).
3. N / non-ACGT runs reported as repeats.
4. O(r · n · (m + k)) per-window `Substring` + suffix-tree lookup.

### Cross-verification (all 0 mismatches; details in `docs/Evidence/REP-DIRECT-001-Evidence.md`)
| Method | Reference | Cases |
|---|---|---|
| `FindDirectRepeats` | `repeat-match -f` (compiled) | 2 000 × 1–200 bp, 200 × ≤ 5 kb (12.8 M pairs), 50 kb–1 Mb; re-run after the WP2 engine refactor: 2 000 (223 155 pairs) |
| `FindDirectRepeats` | brute force (N/IUPAC/lowercase, ±int limits) | 8 000 (373 960 pairs) |
| `FindReverseComplementRepeats` | `repeat-match` (no `-f`) | 3 000 + 100 × ≤ 5 kb (6.6 M pairs) with an `N` sentinel; without it the only differences are repeat-match's shared-`$` leaf loss (31 cases, classified) |
| `FindReverseComplementRepeats` | `vmatch -p`; brute force | 2 000 (186 444 pairs); 3 000 |
| `FindApproximateDirectRepeats` | brute force (both modes); `vmatch -h k -allmax` | 3 000; 3 040 (1.09 M repeats) |
| `FindSupermaximalRepeats` | `vmatch -supermax`; brute force | 3 060 (38 434 pairs); 3 000 |
| `FindDegenerateRepeats` (4 modes) | brute force of Vmatch App. A; Vmatch 2.3.1 built from source with the left-extension seed shortcut disabled; stock Vmatch | 6 000 (78 507 repeats): 0 / 0; + 800 × 100–1 500 bp (534 587 repeats) and 1 Mb (2 333): 0 vs shortcut-free Vmatch; stock Vmatch differs only for edit mode (71 cases), only through the shortcut |
| `FindDegenerateRepeats(…, reporting, vmatchCompatible)` (WP15: `BestPerSeed` = Vmatch default without `-allmax`; `vmatchCompatible` = stock shortcut + first-seed distance label) | stock `vmatch` (compatible) / source build with `VM_NOPRUNE` + `VM_NOPRUNE_H` (complete), 16 configurations, multiset rows | 6 000 × 8–50 bp (1 423 142 rows) + 2 000 Hamming k ≤ 4 (233 692) + 300 × 100–1 500 bp (2 338 723) + 1 Mb / 200 kb (19 166): 0, except the documented default-mode distance label (1 case, finding 3 of Evidence §WP15) |

Locked values: `AAAAAATTTTAAAAAA` 4–6 → the 5 pairs above; `ACGTACGTTTTTTTTTACGTACGT` min 4 → (0,16,8) (0,20,4)
(3,15,5) (7,12,4) at spacing ≥ 1; RC `AAAAAAAACGTTGCAACGTAAAA` min 3 → (6,6,6) (7,7,12) (15,15,4) (repeat-match
`7 12r 6`, `8 19r 12`, `16 19r 4`); k-mismatch (0,21,17, 2 mm) (vmatch -h 2); supermaximal `CAGCAG` @ 0,3,12.

### Tests
`RepeatFinder_DirectRepeat_Tests` (re-locked to repeat-match), `RepeatFinder_RepeatVariants_Tests` (new, 21),
`RepeatFinder_DegenerateRepeats_Tests` (WP8, 12: Vmatch-locked lists, brute force of the definition in all 4 modes),
differential / combinatorial / snapshot / fuzz / property tests (F11), heavy tier `RepDirectVariantsProperties`
(3 brute-force oracles), `RepDirectVariantsMetamorphicTests` (5 relations), `RepDirectVariantsFuzzTests` (5);
MCP `FindDirectRepeatsTests`.

## Verdict & follow-ups
- Stage A: FAIL → corrected. Stage B: FAIL → fixed. **State: FIXED**; variants implemented and reference-identical.
- MCP: `find_direct_repeats` delegates to `FindDirectRepeats`. ~~The three variants are C# API only~~ — resolved by B04 F49:
  `find_reverse_complement_repeats`, `find_approximate_direct_repeats`, `find_degenerate_repeats`, `find_supermaximal_repeats`
  (Analysis server, delegating; tool counts updated additively).
- ~~REPuter/Vmatch k-differences (edit-distance, `vmatch -e`) repeats are not provided (Hamming only).~~ Implemented in WP8 (`FindDegenerateRepeats`, B04 F47), incl. approximate palindromic repeats (`-p -h`, `-p -e`); Vmatch's left-extension seed shortcut (incomplete for edit matches) is documented and, since WP15 (B04 F59), reproducible on request (`vmatchCompatible`); Vmatch's default best-per-seed output is `DegenerateRepeatReporting.BestPerSeed` (F58). Tests `RepeatFinder_VmatchReporting_Tests` (10) + 2 MCP.
