# Evidence Artifact: RESTR-FIND-001

**Test Unit ID:** RESTR-FIND-001
**Algorithm:** Restriction Site Detection (enzyme table, `FindSites`, `FindAllSites`, `FindPatternOccurrences`, `CreateMap` site reporting)
**Date Collected:** 2026-10-14
**Reviewer note:** every number in this file was computed in the 2026-10-14 review session; nothing is carried over from an earlier pass.

---

## Reference Oracle

### Biopython 1.88 — `Bio.Restriction`

**Package:** `biopython` 1.88, installed in the validation environment
(`/usr/local/lib/python3.13/dist-packages/Bio/Restriction/`)
**Accessed:** 2026-10-14, by importing and running it (no network access needed)
**Authority rank:** 1 for coordinates/cut arithmetic (the reference implementation used throughout
this campaign), 1 for the enzyme table (it is a verbatim REBASE snapshot — see below)

Reference code read this session, with line numbers in
`/usr/local/lib/python3.13/dist-packages/Bio/Restriction/Restriction.py`:

| Member | Lines | What it establishes |
|---|---|---|
| `FormattedSeq.__init__` / `_table` | 146-169 | The sequence is upper-cased and prefixed with one space, making `finditer` offsets 1-based |
| `FormattedSeq.finditer` | 215-228 | Linear search over `data`; **circular** search over `data + data[1:size]`, so a site may span the origin |
| `AbstractCut._drop` | 542-568 | **Linear:** keep a cut only when `1 < cut_on_watson <= length` **and** `1 < cut_on_crick <= length`, with `cut_on_crick = cut_on_watson - cls.ovhg`. **Circular:** fold out-of-range cuts back by ±`length` |
| `Defined._modify` | 789-814 | Forward-strand cut `= location + fst5` (docstring worked example: EcoRI at 1-based 10 → 11) |
| `Defined._rev_modify` | 817-825 | Minus-strand site's **watson** cut `= location - fst3` |
| `Palindromic._search` | 967-978 | One result per match, then `_drop` |
| `NonPalindromic._search` | 996-1019 | Plus matches through `_modify`, minus matches through `_rev_modify`, then one shared `_drop` |
| `Restriction_Dictionary.rest_dict["EcoRI"]["compsite"]` | — | `'(?=(?P<EcoRI>GAATTC))'` — a **lookahead**, so overlapping sites are reported; `SfiI` is `'(?=(?P<SfiI>GGCC.....GGCC))'`, i.e. the site's `N` expands to `.` |

Derived facts used below:

- `ovhg = fst5 - (len(site) + fst3) = CutPositionForward - CutPositionReverse`.
- Biopython's 1-based watson cut equals this library's forward-strand `CutPosition + 1`
  (verified: `EcoRI.search(Seq("AAAGAATTCAAA")) == [5]`, Seqeron forward `CutPosition == 4`).
- `_drop`'s linear filter in 0-based terms is `1 <= cut <= length - 1` on **both** strands — i.e.
  the cut must fall on a phosphodiester bond that actually exists.
- `crick = watson - ovhg` holds for minus-strand sites too: `(p + L - rev) - (fwd - rev) = p + L - fwd`.

### REBASE — EMBOSS files version 404 (2024)

**Primary host:** `http://rebase.neb.com/` / `https://rebase.neb.com/` — **not reachable** from this
environment. Attempts made 2026-10-14: `rebase.neb.com` does not resolve over plain HTTP and the
agent proxy answered **403 to CONNECT** for `rebase.neb.com:443`
(recorded in the proxy's own `recentRelayFailures`). `en.wikipedia.org` is denied the same way, so the
Wikipedia enzyme tables the pre-existing TestSpec cites could not be re-fetched either.
**Substitute used:** the REBASE snapshot that Biopython ships. The generated file
`Bio/Restriction/Restriction_Dictionary.py` states its own provenance in its module docstring,
verbatim: *"Used REBASE emboss files version 404 (2024)."* It carries 1088 enzymes with
`site`, `fst5`, `fst3`, `ovhg`, `ovhgseq`, `suppl` and a REBASE id/URI per enzyme
(e.g. EcoRI id 993, `https://identifiers.org/rebase:993`). It carries **no organism and no
methylation field**, which is why those two columns needed the separate `allenz` source below.

**Pinned and corroborated (WP-SRC, 2026-10-14):** the installed file is
sha256 `2a79099295dbad6061ea67a11e053787c591fcb2eb10fc8c0f89ead908dfa02b` (562 312 B, 1088 enzymes)
and is **byte-identical** to the same path at the pinned tag
`https://raw.githubusercontent.com/biopython/biopython/biopython-188/Bio/Restriction/Restriction_Dictionary.py`
(200, 562 312 B, same sha256) — so the snapshot is reproducible without reaching REBASE.
The organism and methylation columns come from REBASE **`allenz.104`** (REBASE version 104,
Copyright Dr. Richard J. Roberts 2021, released Mar 31 2021) vendored in `crisprVerse/crisprBase`.
Full dataset list, hashes, mirror cross-check numbers and the reachability record are in
[`B08-Source-Provenance.md`](B08-Source-Provenance.md).

---

## 1. Enzyme-table audit — all 39 entries vs REBASE 404

Method: the 39 dictionary literals were parsed out of `RestrictionAnalyzer.cs` and each compared on
four fields against `Bio.Restriction`: recognition sequence, `CutPositionForward` vs `fst5`,
`CutPositionReverse` vs `len(site) + fst3`, and the blunt/5'/3' class vs the reference's own
`is_blunt()` / `is_5overhang()` / `is_3overhang()`.

**Result: 39 entries checked, 39 matched on all four fields, 0 wrong entries.**
(37 distinct enzymes: `MspI`/`HpaII` both `CCGG`, `DpnI`/`Sau3AI`/`MboI` all `GATC`.)

| Entry | Seqeron site / fwd / rev | REBASE 404 site / fst5 / len+fst3 | elucidate() | Overhang (Seqeron / REBASE) | Verdict |
|---|---|---|---|---|---|
| EcoRI | GAATTC / 1 / 5 | GAATTC / 1 / 5 | G^AATT_C | 5' / 5' | ✅ |
| BamHI | GGATCC / 1 / 5 | GGATCC / 1 / 5 | G^GATC_C | 5' / 5' | ✅ |
| HindIII | AAGCTT / 1 / 5 | AAGCTT / 1 / 5 | A^AGCT_T | 5' / 5' | ✅ |
| XhoI | CTCGAG / 1 / 5 | CTCGAG / 1 / 5 | C^TCGA_G | 5' / 5' | ✅ |
| SalI | GTCGAC / 1 / 5 | GTCGAC / 1 / 5 | G^TCGA_C | 5' / 5' | ✅ |
| PstI | CTGCAG / 5 / 1 | CTGCAG / 5 / 1 | C_TGCA^G | 3' / 3' | ✅ |
| SphI | GCATGC / 5 / 1 | GCATGC / 5 / 1 | G_CATG^C | 3' / 3' | ✅ |
| KpnI | GGTACC / 5 / 1 | GGTACC / 5 / 1 | G_GTAC^C | 3' / 3' | ✅ |
| SacI | GAGCTC / 5 / 1 | GAGCTC / 5 / 1 | G_AGCT^C | 3' / 3' | ✅ |
| XbaI | TCTAGA / 1 / 5 | TCTAGA / 1 / 5 | T^CTAG_A | 5' / 5' | ✅ |
| NcoI | CCATGG / 1 / 5 | CCATGG / 1 / 5 | C^CATG_G | 5' / 5' | ✅ |
| NdeI | CATATG / 2 / 4 | CATATG / 2 / 4 | CA^TA_TG | 5' / 5' | ✅ |
| NheI | GCTAGC / 1 / 5 | GCTAGC / 1 / 5 | G^CTAG_C | 5' / 5' | ✅ |
| SpeI | ACTAGT / 1 / 5 | ACTAGT / 1 / 5 | A^CTAG_T | 5' / 5' | ✅ |
| AvrII | CCTAGG / 1 / 5 | CCTAGG / 1 / 5 | C^CTAG_G | 5' / 5' | ✅ |
| ClaI | ATCGAT / 2 / 4 | ATCGAT / 2 / 4 | AT^CG_AT | 5' / 5' | ✅ |
| AgeI | ACCGGT / 1 / 5 | ACCGGT / 1 / 5 | A^CCGG_T | 5' / 5' | ✅ |
| NotI | GCGGCCGC / 2 / 6 | GCGGCCGC / 2 / 6 | GC^GGCC_GC | 5' / 5' | ✅ |
| SfiI | GGCCNNNNNGGCC / 8 / 5 | GGCCNNNNNGGCC / 8 / 5 | GGCCN_NNN^NGGCC | 3' / 3' | ✅ |
| PacI | TTAATTAA / 5 / 3 | TTAATTAA / 5 / 3 | TTA_AT^TAA | 3' / 3' | ✅ |
| AscI | GGCGCGCC / 2 / 6 | GGCGCGCC / 2 / 6 | GG^CGCG_CC | 5' / 5' | ✅ |
| FseI | GGCCGGCC / 6 / 2 | GGCCGGCC / 6 / 2 | GG_CCGG^CC | 3' / 3' | ✅ |
| SwaI | ATTTAAAT / 4 / 4 | ATTTAAAT / 4 / 4 | ATTT^_AAAT | blunt / blunt | ✅ |
| MspI | CCGG / 1 / 3 | CCGG / 1 / 3 | C^CG_G | 5' / 5' | ✅ |
| HpaII | CCGG / 1 / 3 | CCGG / 1 / 3 | C^CG_G | 5' / 5' | ✅ |
| TaqI | TCGA / 1 / 3 | TCGA / 1 / 3 | T^CG_A | 5' / 5' | ✅ |
| AluI | AGCT / 2 / 2 | AGCT / 2 / 2 | AG^_CT | blunt / blunt | ✅ |
| RsaI | GTAC / 2 / 2 | GTAC / 2 / 2 | GT^_AC | blunt / blunt | ✅ |
| HaeIII | GGCC / 2 / 2 | GGCC / 2 / 2 | GG^_CC | blunt / blunt | ✅ |
| DpnI | GATC / 2 / 2 | GATC / 2 / 2 | GA^_TC | blunt / blunt | ✅ |
| Sau3AI | GATC / 0 / 4 | GATC / 0 / 4 | N^GATC_N | 5' / 5' | ✅ |
| MboI | GATC / 0 / 4 | GATC / 0 / 4 | N^GATC_N | 5' / 5' | ✅ |
| EcoRV | GATATC / 3 / 3 | GATATC / 3 / 3 | GAT^_ATC | blunt / blunt | ✅ |
| SmaI | CCCGGG / 3 / 3 | CCCGGG / 3 / 3 | CCC^_GGG | blunt / blunt | ✅ |
| HincII | GTYRAC / 3 / 3 | GTYRAC / 3 / 3 | GTY^_RAC | blunt / blunt | ✅ |
| ScaI | AGTACT / 3 / 3 | AGTACT / 3 / 3 | AGT^_ACT | blunt / blunt | ✅ |
| StuI | AGGCCT / 3 / 3 | AGGCCT / 3 / 3 | AGG^_CCT | blunt / blunt | ✅ |
| BglII | AGATCT / 1 / 5 | AGATCT / 1 / 5 | A^GATC_T | 5' / 5' | ✅ |
| ApaI | GGGCCC / 5 / 1 | GGGCCC / 5 / 1 | G_GGCC^C | 3' / 3' | ✅ |

Verified since (WP-SRC, 2026-10-14): the `Organism` string of each entry — **39/39 agree at genus +
species** against REBASE `allenz` field `<3>`, which Biopython's snapshot does not carry but which is
vendored on an allowed host. REBASE's string is consistently Seqeron's plus a strain suffix (e.g.
`EcoRI` *Escherichia coli* vs REBASE *Escherichia coli RY13*). Checked against two independent
copies, which agree with each other on all 39. Source-derived, not executed; no algorithm reads it.
Datasets, hashes and the mirror cross-check are in
[`B08-Source-Provenance.md`](B08-Source-Provenance.md) §3-§4.

---

## 2. Site-search comparison against Biopython

### 2.1 Harness

A throwaway .NET console dumper (scratchpad only, not checked in) drives
`RestrictionAnalyzer.FindSites` and prints `(Position, CutPosition, RecognizedSequence)` per record
per strand; a Python driver runs `Bio.Restriction.<enzyme>.search(Seq(s), linear=...)` on the same
inputs and compares:

1. the **set of watson-strand cuts** Seqeron implies (`Position + fst5` for forward records,
   `Position + L - rev` for reverse records) against Biopython's result minus 1;
2. each forward record's `CutPosition == Position + fst5`;
3. each reverse record's `CutPosition == watson - ovhg` (the reference's own crick relation);
4. each record's `RecognizedSequence` against the window read on its own strand.

Case set: 21 designed sequences (including the self-overlapping `TTAATTAATTAA` and
`GGCCAAAAAGGCCAAAAAGGCC`, all four `GTYRAC` expansions, sites filling the whole molecule, and
origin-spanning layouts) plus **400 pseudo-random A/C/G/T sequences, length 4-120, Python
`random.Random(20261014)`**, each run against all 39 built-in enzymes, linear **and** circular; plus
the 20 REBASE non-palindromic in-site cutters (see §2.3) over 128 sequences, linear.

### 2.2 Counts

| Pass | Comparisons | Mismatches |
|---|---|---|
| **Before** (pre-fix algorithm transcribed from `RestrictionAnalyzer.cs` @12250f9b, linear only — the pre-fix code had no circular mode) | 18,706 | **144** |
| **After** (this unit, linear + circular + non-palindromic customs) | 34,906 | **0** |

492 further cells are excluded from the "after" figure — every (enzyme, sequence, topology) cell where
the molecule is **shorter** than the recognition site (219 linear + 219 circular built-in cells and 54
non-palindromic cells). The reference and this library agree on all the linear ones (both report no
site); the circular ones are the one deliberate divergence, §3. The "before" figure uses the same
guard, which is why its 18,706 linear cells and the after figure's 34,906 (linear + circular) are
consistent: 2·219 + 54 = 492.

The 144 pre-fix mismatches break down as:

| Kind | Count | Defect |
|---|---|---|
| `revrecog` — reverse record's `RecognizedSequence` was the forward-strand window | 121, over 16 distinct enzymes — of the **built-ins**, `HincII` and `SfiI` (the two degenerate palindromes, whose two strands read differently at a non-self-complementary instance); the rest are the non-palindromic customs of §2.3 | F19 |
| `revcrick` — reverse record's `CutPosition` used `Position + CutPositionReverse` | 15 (`CdiI` 8, `SimI` 7) | F19 |
| `watson` — a site was reported whose cut falls outside the linear molecule | 8 (`Sau3AI` 4, `MboI` 4 — the only two table entries with a cut at offset 0) | F17 |

### 2.3 Non-palindromic in-site cutters (for F19)

A `RestrictionEnzyme` can express only cuts **inside** the recognition window. Scanning all 1088
REBASE 404 enzymes for `0 <= fst5 <= L`, `0 <= L + fst3 <= L` and a pure-ACGT site gives 20
non-palindromic candidates; of those, exactly **two** have asymmetric cuts
(`fst5 + (L + fst3) != L`), which is the only condition under which the pre-fix reverse-strand
formula differs from the correct one:

| Enzyme | Site | fst5 | len+fst3 | `elucidate()` | `ovhg` |
|---|---|---|---|---|---|
| CdiI | CATCG | 4 | 4 | `CATC^_G` | 0 |
| SimI | GGGTC | 2 | 5 | `GG^GTC_N` | -3 |

Locked numbers:

- `CdiI.search(Seq("AAACATCGAAACGATGAAA")) == [8, 13]` → 0-based watson cuts 7 and 12. Seqeron now
  reports forward (pos 3, cut 7) and reverse (pos 11, cut 12, recognized `CATCG`); **before the fix the
  reverse record read (pos 11, cut 15, recognized `CGATG`)**.
- `SimI.search(Seq("AAAGGGTCAAAGACCCAAA")) == [6, 12]` → 0-based watson cuts 5 and 11. Seqeron now
  reports forward (pos 3, cut 5) and reverse (pos 11, cut 14 = 11 − (−3), recognized `GGGTC`);
  **before, the reverse record read (pos 11, cut 16, recognized `GACCC`)**.

### 2.4 Individually locked reference values (all computed 2026-10-14)

| Call | Biopython result (1-based) | Seqeron after the fix |
|---|---|---|
| `EcoRI.search(Seq("AAAGAATTCAAA"))` | `[5]` | fwd (3, 4), rev (3, 8) |
| `EcoRI.search(Seq("GAATTC"))` | `[2]` | fwd (0, 1), rev (0, 5) |
| `Sau3AI.search(Seq("GATCGATCGATC"))` | `[5]` | fwd (4, 4) only — the GATC at 0 and 8 are not cuttable |
| `Sau3AI.search(Seq("GATC"))`, `MboI.search(Seq("GATC"))` | `[]`, `[]` | no sites |
| `Sau3AI.search(Seq("AAAAAAAAGATC"))` | `[]` | no sites |
| `Sau3AI.search(Seq("GATCGATCGATC"), linear=False)` | `[1, 5, 9]` | fwd cuts 0, 4, 8 |
| `Sau3AI.search(Seq("AAAAAAAAGATC"), linear=False)` | `[9]` | fwd (8, 8), rev (8, 0) — the crick cut at 12 is the origin bond |
| `EcoRI.search(Seq("TTCAAAAAAAGAA"))` / `..., linear=False` | `[]` / `[12]` | none / fwd (10, 11) `GAATTC` |
| `BamHI.search(Seq("CCAAAAAAAAGGAT"))` / `..., linear=False` | `[]` / `[12]` | none / fwd (10, 11) |
| `PacI.search(Seq("TTAATTAATTAA"))` | `[6, 10]` | fwd (0, 5) and (4, 9) — overlapping sites both reported |
| `SfiI.search(Seq("GGCCAAAAAGGCCAAAAAGGCC"))` | `[9, 18]` | fwd (0, 8) and (9, 17) |
| `HincII.search(Seq("AAGTCAACAA"))` | `[6]` | fwd (2, 5) `GTCAAC`, rev (2, 5) `GTTGAC` (blunt, so both cuts coincide); **before the fix the reverse record read `GTCAAC`** |
| `ClaI.search(Seq("GATC"), linear=False)` | `[4]` | **no sites** (divergence, §3) |

---

## 3. Deliberate divergence from the oracle

**Circular molecule shorter than the recognition site.** `FormattedSeq.finditer` (L215-228) builds the
circular search string as `data + data[1:size]` unconditionally, so when `len(seq) < len(site)` the
molecule is wrapped more than once and a site is reported that cannot physically be present:
`ClaI.search(Seq("GATC"), linear=False) == [4]` for the 6-bp site `ATCGAT` on a 4-bp circle. This
library returns no sites when `Length < RecognitionLength`, on either topology. Locked by
`FindSites_Circular_MoleculeShorterThanSite_ReportsNothing`.

One reference quirk checked and found **not** to be a divergence: `_drop` computes
`cut_on_crick = cut_on_watson - ovhg` for minus-strand results too, using the plus-orientation `ovhg`.
That is arithmetically correct for minus-strand sites — `(p + L - rev) - (fwd - rev) = p + L - fwd` —
so the filter is sound on both strands and this library uses the same rule.

---

## 4. Not verified / open

| Item | Sources tried | Status |
|---|---|---|
| Each entry's `Organism` string | ~~`rebase.neb.com`, `en.wikipedia.org` (both proxy 403 to CONNECT), Biopython's snapshot (no organism field)~~ → **resolved** via REBASE `allenz.104` vendored in `crisprVerse/crisprBase` on `raw.githubusercontent.com` | **Verified 39/39** (genus+species) by WP-SRC, corroborated on a second mirror — see [`B08-Source-Provenance.md`](B08-Source-Provenance.md) §4c |
| REBASE release newer than 404 | same blocked hosts | Not attempted further; 404 (2024) is the newest snapshot available offline |
| `RestrictionEnzyme`'s overhang **sequence** (`GetOverhang`) against REBASE `ovhgseq` | Biopython exposes `ovhgseq` per enzyme (e.g. EcoRI `AATT`, `ovhg` −4) | Out of this unit — belongs to RESTR-FILTER-001; recorded in `docs/Validation/review-2026-09/B08.md` |
