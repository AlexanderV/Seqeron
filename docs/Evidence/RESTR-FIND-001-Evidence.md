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

### 2.5 Ambiguity codes in the searched sequence (B08 WP3 / audit item D3, F40)

**Reference rule** (Biopython 1.88, installed locally 2026-10-14): `FormattedSeq` keeps every letter of
the sequence (upper-cased; only whitespace and digits are removed — `_table`, `_remove_chars`), and each
enzyme's `compsite` compiles a site `N` to `.` and any other code to the class of its A/C/G/T expansion:
`SfiI` `(?=(?P<SfiI>GGCC.....GGCC))`, `HincII` `(?=(?P<HincII>GT[CT][AG]AC))`, `AccI`
`(?=(?P<AccI>GT[AC][GT]AC))`. So a sequence `N` (or `R`, `Y`, …) is matched by a site `N` and by
nothing else. The pre-fix code let a sequence ambiguity code match **nothing**, so `SfiI` on
`GGCCNNNNNGGCC` reported no site where `SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]`.

**Harness** (scratchpad, not checked in; same comparisons 1-4 as §2.1): 816 sequences — 16 designed
(`GGCCNNNNNGGCC`, flanked/lower-case/origin-spanning variants, `GAANTC`, `GTYRAC`/`GTNNAC`,
`GGCNNNNNNGGCC`, all-`N`), 300 random A/C/G/T with `N` at p = 0.15, 300 random over all 15 IUPAC codes
(ambiguity p = 0.3), 200 `SfiI`-enriched `GGCC`+5+`GGCC` with random 15-code spacers and flanks;
lengths 4-60, Python `random.Random(20261014)` — × 39 built-ins × linear and circular, through the
raw-string `FindSites(string, string, MoleculeTopology)`. Separately, `FindPatternOccurrences` against
each Biopython `compsite`'s forward group for **all 513** single-site REBASE 404 enzymes carrying an
ambiguity code (the one `|`-joined multi-site entry excluded), over 120 random sequences plus the
all-`N` sequence and the site itself.

| Pass | Comparisons | Mismatches before (@c36e74fa) | Mismatches after |
|---|---|---|---|
| `FindSites`, 39 built-ins × 816 sequences × 2 topologies (cells with Length < site excluded, §2.2) | 62,112 | **409** (all `SfiI` — the only built-in with an `N`; all of kind `watson`: a cut Biopython reports that Seqeron missed) | **0** |
| `FindPatternOccurrences`, 513 ambiguous REBASE sites | 62,586 | **475** (top: `CNNR` 76, `CNNG` 47, `CCNNNNNNNGG` 35, `GGNCC` 24, `CCNGG` 18) | **0** |

No mismatch in either "before" pass was in the other direction (Seqeron reporting a site Biopython
does not) — the old rule was strictly narrower; degenerate site codes other than `N` already behaved
like Biopython, which the "after" figures confirm across all 513 ambiguous sites.

**Fix.** One predicate (private to `RestrictionAnalyzer` here; since B08 round 5 A5-13 the internal `IupacPattern.MatchesCode`, shared with the methylation mask and the CRISPR PAM scan, which had each kept a copy), `SequenceBaseMatchesSiteCode(sequenceBase, siteCode)` =
`siteCode == 'N' ? IupacHelper.IsNucleotideCode(sequenceBase) : IupacHelper.MatchesIupac(sequenceBase, siteCode)`,
composed from the canonical Core helpers — no new Core API and no second ambiguity table. Note the
rule is **not** "codes intersect" (as `CodesIntersect` in the overhang-compatibility code is): under
intersection `GAANTC` would match `EcoRI` and `GTYRAC` would match `HincII`, and Biopython reports
neither (`EcoRI.search(Seq("AGAANTCA")) == []`, `HincII.search(Seq("AGTYRACA")) == []`).

Reachability: only the raw-string `FindSites` overloads and `FindPatternOccurrences` — `DnaSequence`
rejects non-ACGT, so `Digest`, `CreateMap`, `FindAllSites` and every MCP tool (all of which build a
`DnaSequence`) are unchanged.

Locked values (Biopython 1.88, 1-based → Seqeron 0-based top-strand cut): `SfiI` `GGCCNNNNNGGCC`
`[9]` → fwd (0, 8), rev (0, 5), linear and circular; `AGGCCNNNNNGGCCA` `[10]` → 9;
`ggccnnnnnggcc` `[9]` → 8; `NNNNNGGCCAGGCC` circular `[5]` → 4; `GGCCRYSWKGGCC` `[9]` → 8;
`EcoRI` `AGAANTCA`, `HincII` `AGTYRACA` / `AGTNNACA`, `SfiI` `GGCNNNNNNGGCC`, `EcoRI` `N`×16 → `[]`.

### 2.6 Single-sided distal (Type IIS) cutters (B08 WP11a / audit item D6 part 1, F54, 2026-10-15)

**Claim retracted.** This unit had recorded that a `RestrictionEnzyme` "carries one cut offset per
strand inside the recognition window, so Type IIS / Type IIG distal cutters cannot be expressed".
The record never enforced that: both offsets are signed `int`s with no range check, and the site
finder's formulas (§2.1 conventions) are Biopython's own `OneCut._modify` (`location + fst5`) and
`OneCut._rev_modify` (`location − fst3`), read with `inspect.getsource` from the installed 1.88.

**Offset mapping (stated once, used everywhere):** `CutPositionForward = fst5`,
`CutPositionReverse = size + fst3`, so `ovhg = CutPositionForward − CutPositionReverse` and REBASE
`SITE(n/m)` = `(L + n, L + m)`. Values from Biopython 1.88: BsaI GGTCTC `fst5 7, fst3 5` → (7, 11);
BsmBI CGTCTC (7, 11); BbsI GAAGAC `8, 6` → (8, 12); SapI GCTCTTC `8, 4` → (8, 11); FokI GGATG
`14, 13` → (14, 18); MlyI GAGTC `10, 5` → (10, 10) blunt; AcuI CTGAAG `22, 14` → (22, 20), 3' `ovhg 2`.

**Snapshot census (Biopython 1.88, 1088 enzymes):** 576 one-cut in-site, **153 one-cut distal**,
25 two-cut distal (`TwoCuts`: BaeI, BcgI, …), 334 without a defined cut. Among the 153: 16 have one
strand cut inside and one outside (e.g. BsmI GAATGC (7, 5)); exactly two cut upstream of the site
(TscAI/TspRI CASTG (7, −2), palindromic); four are palindromic with distal cuts (HauII, SgeI, TscAI,
TspRI). No one-cut enzyme has a negative `fst5`.

**Harness** (scratchpad, not checked in): a .NET dumper over `FindSites(DnaSequence, RestrictionEnzyme,
MoleculeTopology)`, `Digest(…, IEnumerable<RestrictionEnzyme>)` and `PartialDigest(…)`, driven by a
Python script that builds every OneCut enzyme as a custom record, plants 0–4 sites (either orientation,
some across the origin) into random ACGT sequences of length `L … 300`, and compares on both
topologies: the multiset of top-strand cuts (forward records' `CutPosition`; reverse records'
`CutPosition + fwd − rev` unless the window also matches forward) against `search(seq, linear)`
reduced mod `n`; that every circular `CutPosition` is in `[0, n)`; fragment start/length against
`catalyse` (linear) / the sorted cut list (circular, see RESTR-DIGEST-001 oracle caveat); both end
overhangs against the window between `cut_on_watson` and `cut_on_crick = cut_on_watson − ovhg`;
the partial-digest multiset against all cut-list intervals.

| run (20 / 60 seq. per distal enzyme; ¼ that per in-site) | distal cases | distal cuts | before (cases) | after |
|---|---|---|---|---|
| N = 20 | 6,120 | 14,331 | circular `CutPosition` negative: **836** | 0 |
| N = 60 | 18,360 | 45,193 | circular `CutPosition` negative: **2,575** | 0 (+5 documented TspRI/TscAI) |

Linear `FindSites` cut multisets already matched before the fix (0 mismatches): the linear formulas
were correct for distal offsets. The circular finder reduced with C#'s signed `%`, so a reverse-strand
site whose bottom-strand cut fell before the origin (FokI `CATCC` at 0: `0 + 5 − 14 = −9`) was
reported at a negative position. Fixed with a non-negative modulo. Locked by
`RestrictionAnalyzer_DistalCutter_Tests` (FokI circular → 33, top-strand 1-based 30 = Biopython
`[30]`; BsaI forward/reverse records; FokI linear drop at the molecule end, Biopython `[17]` / `[]`).

**Documented divergence (Biopython data defect).** Biopython 1.88 stores TscAI/TspRI with
`ovhg = 10` although `fst5 − (size + fst3) = 7 − (5 − 7) = 9` and `ovhgseq = "NNCASTGNN"` (9 nt);
these are the only 2 of 729 OneCut enzymes where the stored `ovhg` disagrees with `fst5`/`fst3`.
`_drop` uses the stored value, so a linear site whose bottom-strand cut is the bond after base 0 is
dropped: `TspRI.search(Seq("AAACAGTGAAAAAA")) == []`, where this library reports the site (top cut 10,
bottom cut 1). In the N = 60 run this is the only class of difference (5 cases), all of exactly that
shape. Locked by `FindSites_TspRI_BottomCutAfterBaseZero_IsKept_DocumentedDivergence`.

*(Both items below were closed by WP11b — §2.7 and RESTR-FILTER-001 A8.)* **Not in this unit:** two-sided Type IIB cutters (two cut pairs per site; Biopython `TwoCuts` stores the
second pair as `charac[2:4]` = `scd5`/`scd3` and `search` yields both — `BaeI.search` on
`A×20 ACAAAAGTACC T×20` = `[11, 44]`, `catalyse` 10/33/8 — only its `elucidate` answers "cut twice, not
yet implemented sorry"); enzyme-level overhang/compatibility of distal ends
(`GetOverhangSequence` → `null`, `AreCompatible` → `false`; Biopython classes them `Ambiguous` with
`ovhgseq` `NNNN`) — both handed to WP11b.

---

### 2.7 Two-sided (Type IIB) cutters (B08 WP11b / D6 part 2, F56, 2026-10-15)

**Reference read** (Biopython 1.88, installed locally, `Restriction.py`): `TwoCuts._modify` yields
`location + fst5` and `location + scd5`; `TwoCuts._rev_modify` yields `location − fst3` and
`location − scd3`; `_drop` filters the cut list element by element with one `ovhg`. Census: **25**
`TwoCuts` enzymes (AjuI, AlfI, AloI, ArsI, BaeI, BarI, BceSIV, BcgI, BdaI, BplI, BsaXI, Bsp24I, CjeI,
CjePI, CspCI, FalI, Hin4I, NgoAVIII, NmeDI, PpiI, PsrI, RdeGBIII, SdeOSI, TstI, UcoMSI; 7 palindromic);
for every one `fst5 − (size + fst3) == scd5 − (size + scd3) == ovhg`, so one stagger describes both cuts.

**Mapping:** `SecondCut = (scd5, size + scd3)`, e.g. BaeI `("ACNNNNGTAYC", −10, −15) { SecondCut = (23, 18) }`.

**Differential** (harness `wp11b/diff2.py` + `dumper` in the session scratchpad; all 25 enzymes × 60
planted random sequences × linear and circular = 3 000 cases; 7 227 Biopython cuts, 8 859 fragments,
17 718 end overhangs, 11 332 partial-digest fragments compared):

| | before (HEAD `3bb17a2d`, first pair only — no way to express the second) | after |
|---|---|---|
| `FindSites` top-strand cut multiset ≠ `search` | 2 029 cases | **0** |
| `Digest` fragments ≠ fragments of the sorted distinct cuts | 1 988 | **0** |
| end overhangs ≠ window between `cut_on_watson` and `cut_on_watson − ovhg` | — | **0** |
| `PartialDigest` multiset | 680 | **0** |

**Biopython-side defects, documented, not reproduced:** (1) `Palindromic._search` does not sort a
two-cut result list (`[p1+fst5, p1+scd5, p2+fst5, …]`), so `catalyse` slices it out of order: UcoMSI
on `A×20 GAGCTC AAA GAGCTC A×30` (65 nt) gives `search == [14, 32, 23, 41]` and fragments
`[13, 18, 0, 18, 25]` (an empty fragment; lengths sum to 74). 111 of the 3 000 cases; this library cuts
at the sorted distinct positions. (2) The circular one-pass fold leaves cuts outside `[1, n]` / unsorted
(66 cases; e.g. BaeI on a 41-nt circle `[26, 18]`); compared modulo `n`, they agree.

## 3. Deliberate divergences from the oracle

**Circular molecule shorter than the recognition site (revised by review 2026-09 B08 H6 → slice F100).**
`FormattedSeq.finditer` (L215-228) builds the circular search string as `data + data[1:size]`, whose
`data` carries a leading space, so the slice is the molecule's first `size - 1` bases. The paragraph
that stood here claimed such a reading "wraps the molecule more than once" and that no site can
physically be present, and the library returned no sites when `Length < RecognitionLength` on either
topology. **That was wrong on both counts.** On a closed circle the window reads past the origin and
back into the molecule, so it does exist; and `ATCGAT` on the 4-bp circle `GATC` starts at index 1
and needs only two bases past the origin — exactly Biopython's own one-pass extension. The length
guard is now scoped to linear topology, and the circular search repeats the molecule to `n + len - 1`
bases. All 9 enzyme/molecule rows that the round-4 72384-row differential had left unexplained now
agree with Biopython exactly (BsiWI/Pfl23II/PspLI/SplI on `ACGT` → cut 2; MslI/DrdV/RseI/SmiMI on
`ACATGA` → 0; AlfI on `TATTGCAAAA` → 4, 8), and a 63336-row re-run (754 REBASE enzymes × 84
molecules) shows 0 divergences at every length that can hold the site.

The remaining divergence is the other way round and is a **superset**: Biopython's `data[1:size]`
slice truncates to the molecule, so for `n < len` it can only reach start indices
`0 .. n + min(n, len-1) - len` and never tests the rest, although those windows are physically
present on the circle. 44 of the 63336 rows, all at `n < len`. Example:
`Bio.Restriction.SnaBI.search(Seq("ACGT"), linear=False) == []`, while `TACGTA` is the window at
start 3 of the circle and this library reports it with cut 2. EMBOSS `restrict -plasmid` reports
nothing at all for `n < len` (checked on all three molecules above), so it does not arbitrate; it
agrees with Biopython on an ordinary wrap (`ACGAAAAAAACGT` + BsiWI → EMBOSS `5prime = 11`,
Biopython `[12]` 1-based). Locked by
`FindSites_Circular_MoleculeShorterThanSite_FindsTheWrappedSite`,
`Digest_Circular_MoleculeShorterThanSite_MatchesBiopython` and
`FindSites_Circular_ShortMolecule_ReportsWindowsBiopythonTruncatesAway`.

**Non-IUPAC character under a site `N` (D3).** Biopython's `.` matches any letter `FormattedSeq`
keeps, so `SfiI.search(Seq("GGCCXXXXXGGCC")) == [9]`. This library lets a site `N` match only the 15
IUPAC nucleotide codes (`IupacHelper.IsNucleotideCode`), keeping the raw-string overloads' documented
contract that a non-IUPAC character never matches. The §2.5 differential is over the 15-code alphabet,
where the two agree exactly. Locked by `FindSites_NonIupacCharacterInSpacer_MatchesNothing`.

**TscAI/TspRI stored `ovhg` (F54).** Biopython keeps `ovhg = 10` against its own `fst5`/`fst3`/`ovhgseq`
(= 9), so its `_drop` rejects a linear site whose bottom-strand cut is the bond after base 0; this
library uses the offsets and keeps it (§2.6).

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

---

## 5. Methylation-site field — REBASE `allenz.104` field `<6>` (B08 D7 / F52, 2026-10-15)

**Source.** REBASE `allenz.104` (REBASE version 104, © R. J. Roberts 2021) as vendored in
`crisprVerse/crisprBase`: `inst/enzymes/allenz.txt` re-downloaded 2026-10-15 into a fresh directory,
sha256 `43b12432ed3c1a19e961d019e3a131336e8ee690a4568a68ddf55ed0e3869388` — **identical** to the pin in
[`B08-Source-Provenance.md`](B08-Source-Provenance.md) §3a; header `inst/enzymes/notation.txt` sha256
`3820fd9dc95c713ea706f2da9c5000975f6f6ab66125d6181448cd562128a64b` (identical). The corroborating
mirror (§3b there) agrees on `<6>` 39/39 for the built-ins (WP-SRC, §4c).

**What the field means** (format header, quoted): *"The site of methylation by the cognate methylase
when known is indicated X(Y) or X,X2(Y,Y2), where X is the base within the recognition sequence that
is modified. A negative number indicates the complementary strand, numbered from the 5' base of that
strand, and Y is the specific type of methylation involved: (6) = N6-methyladenosine, (5) =
5-methylcytosine, (4) = N4-methylcytosine."* In the data, several sites are comma-separated
(`3(6),-3(6)`, `3(5), 2(5)`) and `?` stands for X when only the type is known (`?(4)`). It is the
**protection** site of the enzyme's own R-M system — not whether Dam, Dcm or CpG methylation by another
MTase blocks the enzyme. No source opened in this campaign tabulates that sensitivity for the 39
built-ins (Biopython 1.88 `Bio.Restriction` has no methylation attribute; REBASE's sensitivity pages
are on the blocked host), so this field drives **no** blocking behaviour. (Superseded in part, 2026-10-21, F68: a separate opt-in blocking mode exists, driven by the *methylase* sites of EMBOSS's `Emethylsites.dat` rather than by an enzyme-sensitivity table; see §7. The per-enzyme table remains unavailable.)

**Method.** `extract_meth.py` (16 lines; parses `<n>`-tagged records, prints `<5>` and `<6>` for the 39
table names) was run with `python3 -I` on the downloaded file; a second script rewrote the 39
`_enzymes` initialisers with the extracted `<6>` string (asserting exactly 39 substitutions), and the
same output generated the 39 `[TestCase]`s of `EnzymeDatabase_MethylationField_MatchesRebaseAllenz104`.
Nothing was typed by hand. Extraction output, verbatim:

```
allenz.txt sha256 43b12432ed3c1a19e961d019e3a131336e8ee690a4568a68ddf55ed0e3869388
EcoRI	<5>G^AATTC	<6>3(6)
BamHI	<5>G^GATCC	<6>5(4)
HindIII	<5>A^AGCTT	<6>1(6)
XhoI	<5>C^TCGAG	<6>5(6)
SalI	<5>G^TCGAC	<6>5(6)
PstI	<5>CTGCA^G	<6>5(6)
SphI	<5>GCATG^C	<6>
KpnI	<5>GGTAC^C	<6>4(6)
SacI	<5>GAGCT^C	<6>4(5)
XbaI	<5>T^CTAGA	<6>6(6)
NcoI	<5>C^CATGG	<6>?(4)
NdeI	<5>CA^TATG	<6>4(6)
NheI	<5>G^CTAGC	<6>6(4)
SpeI	<5>A^CTAGT	<6>1(6)
AvrII	<5>C^CTAGG	<6>
ClaI	<5>AT^CGAT	<6>5(6)
AgeI	<5>A^CCGGT	<6>3(5)
NotI	<5>GC^GGCCGC	<6>2(4)
SfiI	<5>GGCCNNNN^NGGCC	<6>?(4)
PacI	<5>TTAAT^TAA	<6>
AscI	<5>GG^CGCGCC	<6>?(5)
FseI	<5>GGCCGG^CC	<6>?(5)
SwaI	<5>ATTT^AAAT	<6>?(6)
MspI	<5>C^CGG	<6>1(5)
HpaII	<5>C^CGG	<6>2(5)
TaqI	<5>T^CGA	<6>4(6)
AluI	<5>AG^CT	<6>3(5)
RsaI	<5>GT^AC	<6>4(4)
HaeIII	<5>GG^CC	<6>3(5)
DpnI	<5>GA^TC	<6>
Sau3AI	<5>^GATC	<6>4(5)
MboI	<5>^GATC	<6>2(6)
EcoRV	<5>GAT^ATC	<6>2(6)
SmaI	<5>CCC^GGG	<6>2(4)
HincII	<5>GTY^RAC	<6>5(6)
ScaI	<5>AGT^ACT	<6>5(4)
StuI	<5>AGG^CCT	<6>
BglII	<5>A^GATCT	<6>5(4)
ApaI	<5>GGGCC^C	<6>4(5)
```

**Summary:** 34 entries carry a site, 5 have an empty `<6>` (SphI, AvrII, PacI, DpnI, StuI); 29 known
positions + 5 `?` positions (NcoI, SfiI, AscI, FseI, SwaI); no built-in has a complementary-strand or
multi-site entry. By type: m6A 15, m5C 10, m4C 9 (incl. `?`). Consistency check (locked in the same
test): every known position lies inside the recognition sequence and names an `A` for `(6)` and a `C`
for `(4)`/`(5)` — 29/29.

**Representation.** `RestrictionEnzyme.RebaseMethylationField` (the string verbatim, validated at
construction; default `""` for user-defined enzymes) and `CognateMethylationSites`
(`CognateMethylationSite(int? Position, bool IsComplementaryStrand, MethylatedBase Type)`), the typed
parse. Parser cases for notation absent from the table (negative position, lists, `?` in a list, and
the undocumented `(99)` type that occurs 3× in the file — rejected) are locked by
`CognateMethylationSites_RebaseNotation_ComplementaryStrandAndLists` and
`RestrictionEnzyme_MalformedMethylationField_Throws`.

## 6. The vendored REBASE set — all enzymes with a defined cut (B08 D4 part 1 / F59, F60, 2026-10-15)

**Inputs (both re-hashed this session, identical to the pins in `B08-Source-Provenance.md`):**
Biopython 1.88 `Bio/Restriction/Restriction_Dictionary.py` sha256
`2a79099295dbad6061ea67a11e053787c591fcb2eb10fc8c0f89ead908dfa02b` (installed copy in the session venv);
REBASE `allenz.104` (`crisprVerse/crisprBase` `inst/enzymes/allenz.txt`) sha256
`43b12432ed3c1a19e961d019e3a131336e8ee690a4568a68ddf55ed0e3869388`.

**Generator:** `scripts/rebase/generate_rebase_enzymes.py` (`table` reads both files with
`ast.literal_eval` / a line parser and refuses any other sha256; `oracle` imports Biopython 1.88 and checks
the installed dictionary's sha256). Outputs: `src/.../MolTools/Resources/rebase_enzymes.tsv` (embedded,
754 rows, provenance + census in its `#` header) and
`tests/.../TestData/Restriction/biopython_rebase_search.tsv` (754 rows, `random.Random(20261015)`).

**Census (Biopython 1.88 class membership, computed with `issubclass` on `Bio.Restriction.AllEnzymes`):**

| Biopython class | Count | Vendored |
|---|---|---|
| `OneCut` + `Defined` (fully defined overhang) | 389 | yes |
| `OneCut` + `Ambiguous` (in-site with `N` overhang, or distal) | 340 | yes |
| `TwoCuts` + `Ambiguous` | 25 | yes (`SecondCut`) |
| `NoCut` + `NotDefined` (`fst5 = fst3 = scd5 = ovhg = None`) | 334 | **no** |
| **total** | **1088** | **754** |

The 729 OneCut split as 576 in-site + 153 distal (F54). Among the 754: 137 blunt, 388 5', 229 3' overhang
(Biopython `is_blunt`/`is_5overhang`/`is_3overhang`); recognition-string length (`.size`) 6 → 385,
4–8 → 676. **Join with allenz.104:** 747 of 754 names found; 7 absent (BscXI, CviKI_1, EcoO157SI, KroNI,
PaqCI, TagI, UpaP162I) → organism and methylation empty. 282 vendored enzymes carry a methylation field;
every one parses under the F52 grammar (no `(99)` type among them). For the 39 common enzymes the vendored
rows equal the curated table on site, both offsets and methylation field, and on genus + species of the
organism (`VendoredRows_OfTheCommonEnzymes_AgreeWithTheCuratedTable`); the loader substitutes the curated
instances, so their behaviour is byte-identical.

**Why the 334 are excluded:** REBASE gives their site but no cut position; `RestrictionEnzyme` must carry
integer offsets, and a digest/overhang/compatibility result built on an invented offset would be false.
They are not resolvable by name (`GetEnzyme("Dpi3084I") == null`). Representing "cut unknown" would need a
record that most members (digest, map, overhang) could not use — not done; recorded here as the decision.

**Licence notes (texts opened this session):**
- Biopython `Restriction_Dictionary.py` header: *"This code is part of the Biopython distribution and governed
  by its license. Please see the LICENSE file that should have been included as part of this package."*
  Biopython 1.88 `LICENSE.rst` (dist-info), Biopython License Agreement: *"Permission to use, copy, modify,
  and distribute this software and its documentation with or without modifications and for any purpose and
  without fee is hereby granted, provided that any copyright notices appear in all copies and that both
  those copyright notices and this permission notice appear in supporting documentation …"* — the vendored
  table's header names the source and licence.
- REBASE `allenz.104` header (`notation.txt`): *"REBASE, The Restriction Enzyme Database
  http://rebase.neb.com — Copyright (c) Dr. Richard J. Roberts, 2021. All rights reserved."* The file
  carries no redistribution grant of its own; it is redistributed inside crisprBase, whose `DESCRIPTION`
  says `License: MIT + file LICENSE` (`LICENSE`: *"Copyright (c) 2022, Genentech, Inc. Permission is hereby
  granted, free of charge …"*; fetched 2026-10-15, sha256 `cacd57a7…` / `023d7b6c…`). Only two factual
  columns (organism, methylation site) of 747 records are taken. **Open for the project owner:** whether
  REBASE's "All rights reserved" notice requires permission from REBASE/NEB for redistributing these two
  columns at this scale (the 39-entry use predates this WP). The site/cut data does not depend on it.

**Differential (by name, 754 enzymes × linear + circular; `RestrictionAnalyzer_RebaseSet_Tests`):**

| check | comparisons | mismatches |
|---|---|---|
| `FindSites` top-strand cuts (F60 rule) = `search(seq, linear=True)` | 754 | **0** |
| same, circular (Biopython cuts reduced into [1, n]) | 754 | **0** |
| `Digest` fragment lengths = gaps of sorted distinct `search` cuts, linear | 754 | **0** |
| same, circular | 754 | **0** |

**F60 — what a naive comparison shows.** Comparing *every* `FindSites` record with `search` gives
9 mismatching (enzyme, topology) cases in 5 enzymes: HauII TGGCCA (palindromic, offsets 17/15 — distal), SgeI CNNG
(palindromic, offsets 13/17 — distal), MspJI CNNR, LpnPI CCDG, SgrTI CCDS (non-palindromic, degenerate; both readings
can match one window). E.g. HauII on the oracle molecule: `FindSites` top cuts `[20, 32]`, `search` `[20]`.
Biopython counts one reading per window, forward first (`Palindromic._search` → `_modify` only;
`NonPalindromic` compsite `(?=(?P<MspJI>C..[AG]))|(?=(?P<MspJI_as>[CT]..G))`), and drops infeasible cuts
only afterwards — so for MspJI and SgeI on a linear molecule, a window whose forward reading cuts off the
end contributes nothing even though its reverse reading would cut inside. Setting aside reverse records
whose window reads forward (pattern match, not "has a forward record") gives 0 mismatches. This is the rule
`Digest`/`CreateMap`/`PartialDigest` already applied (F55/F58); `FindSites` itself reports readings by
design (documented in the class's convention comment and algorithm doc §5.2a). The WP11a figure "FindSites
0 mismatches" (§2.6) was measured with the same rule applied in the harness.

The TscAI/TspRI `ovhg` divergence (§2.6) did not arise on this corpus (no site at the 5' edge).


---

## 7. Opt-in methylation-blocked recognition — EMBOSS `restrict -methylation` (B08 round 2 E2 / F68, 2026-10-21)

**What was wrong.** §5 / the algorithm doc recorded "Not implemented: no source reachable tabulates enzyme
sensitivity". The proof was defeated by the reference tool: EMBOSS `restrict` has `-methylation` ("If this is set
then RE recognition sites will not match methylated bases") driven by `-mfile` (default `Emethylsites.dat`), which
lists the **methylase** sites, not enzyme sensitivities. No REBASE data is involved (the vendored columns are untouched;
`rebase.neb.com` was not opened).

**Sources opened (neutral User-Agent `Seqeron-review/1.0`):**

| Source | Used for |
|---|---|
| `raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/data/Emethylsites.dat` (11 lines; no commit id reachable — the GitHub API refused the session — so pinned by content: byte-identical to `/usr/share/EMBOSS/data/Emethylsites.dat` of Ubuntu `emboss 6.6.0+dfsg-12ubuntu2`) | the five methylase entries: `dam GmATC`, `dcm mCmCAGG`, `dcm mCmCTGG`, `MEcoKI mAmACNNNNNNGTGC`, `MEcoKI GCmACNNNNNNGTT`; header: "any methylated base is preceded by 'm'" |
| `raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/restrict.c` (1032 lines) | `methyl = ajAcdGetBoolean("methylation")` (L128) passed to `embPatRestrictMatch` (L189) |
| `raw.githubusercontent.com/kimrutherford/EMBOSS/master/nucleus/embpat.c` (5822 lines) | `patRestrictReadMethyl` (L1352: `M` -> `N` in the *Replace* string, plain chars incl. pattern `N` copied, so a pattern `N` also becomes `N`; entries pushed with `ajListPush` = front), `patRestrictMethylMod` (L1451: search each entry, overwrite with Replace, "Know out base on opposite strand too": `rp[slen-fpos-1] = 'N'`), `embPatRestrictMatch` (L4607: mask `substr`/`revstr` BEFORE the plasmid doubling, forward pass then reverse pass), `embPatRestrictScan` (L1595: `if(!(*q & v) || v==15) break;` — a sequence `N` never matches, whatever the site base) |
| Debian `restrict.1e` man page (`manpages.debian.org/trixie/emboss/restrict.1e.en.html`) | option text only (the fetch was not needed for semantics) |
| **Executed**: Ubuntu `emboss 6.6.0+dfsg-12ubuntu2` `restrict`, with an enzyme table generated from Biopython 1.88 `Bio.Restriction` (729 enzymes: name, site, length, cut1 = `fst5`, cut2 = `fst5 - ovhg`) because `rebaseextract` needs REBASE files | the oracle values below |

**Semantics confirmed by execution (not assumed):**

1. *Window rule.* `restrict -methylation` on `AAAAGGATCCAAAATGATCAAAAAGATCAAAAACCAGGAAAACCTGGAAAAATTT`, enzymes
   MboI,Sau3AI,BamHI,BclI,EcoRII: without `-methylation` MboI/Sau3AI 5,15,24; BamHI 4; BclI 14; EcoRII 33,42 (0-based
   window starts); **with it: no hit at all** — BamHI included (GGATCC contains the masked A).
2. *Mask probe.* A one-window enzyme `PROBE NN` is reported exactly where the window holds no masked base, so it reads the
   mask out directly. Findings: `AAGATCAAGATCAA` masks one A per GATC (probe windows starting at
   `0,1,4,5,6,7,10,11,12` stay unmasked); `TGAACCGACCAGTGCTAACG` masks the two A's **and the six spacer bases** of
   `mAmACNNNNNNGTGC`; `GCACGGTACCGTTAACTT` masks `GCmAC` + spacer; `CAACCTCATTGTGCCTGGAGGG` (CCTGG shares its first C with
   the GTGC of an EcoKI site) masks both the EcoKI site and the CC — which is only possible if EcoKI is visited **before**
   dcm, i.e. in reverse file order.
3. *Asymmetry.* `GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA`: MboI 10,36 -> blocked; EcoRI 0,18 and SmaI 27 kept;
   TaqI 38,44 **kept** (the TCGA of `GATCGA` starts at the T, which EMBOSS never masks).
4. *Spacer.* `TTTTAACGGTACCGTGCTTTTGAATTCTTTTGCACGGTACCGTTTTTGGTACCTTT`: KpnI/Acc65I 7,35,47 -> **47 only**, EcoRI 21 kept.
5. *Circular (`-plasmid`).* `TAAAAAGAATTCAAAAAAAAAGATCAGAATTCAAAAGAAAAAGATC` (46 bp): EcoRI 6,26 kept; BglII 41 (AGATC|T across the origin) and
   MboI blocked. `ATTCAAAAAGAATTCAAAAAAGATCAAAAAAAAAAAAGA`: EcoRI 9 and 37 (GA|ATTC) kept, MboI blocked.
   `AGGAAAAGAATTCAAAAAAAAAGATCAAAAACC`: EMBOSS keeps EcoRII 31 (CC|AGG across the origin is not masked); the library blocks it (documented divergence).
6. *EMBOSS bug, avoided in the harness:* `methfile` is read once and never rewound, so with several sequences in one run
   only the **first** is masked. The oracle was therefore run one sequence per invocation.

**Differential (2026-10-21, harness outside the repo).** 40 random 500-bp sequences injected with Dam/Dcm/EcoKI motifs and enzyme
sites glued to them at random overlaps; every enzyme of the Biopython table that the library resolves (`GetEnzyme`).
Linear: mask via the `NN` probe 14 051 unmasked windows, **0 mismatching sequences**; sites kept/removed — 50 325
EMBOSS-and-library common sites, 26 439 removed by EMBOSS, **0 mismatching (enzyme, sequence)**. Circular: 58 750 common sites, 31 088 removed,
20 differences, **all** origin-adjacent and all "library removes more" (a methylase site spanning the origin). First
implementation attempts, kept as the audit trail: masking only the `m` bases (the spacer is masked too), file-order entry
visiting (reverse order), and a multi-sequence EMBOSS run (only the first sequence is masked) each produced mismatches that
the run above removed.

**Locked by tests:** `RestrictionAnalyzer_Methylation_Tests` (42 cases: the values above, default-off identity,
reverse strand, circular wrap, rotation invariance, digest / summary / map) and MCP `RestrictionMethylationToolsTests`.

**Limitations kept honest.** (i) A window rule, not an enzyme-sensitivity table: BamHI (Dam-insensitive at the bench) is blocked,
`DpnI` (needs methylation) and `MspI`/`HpaII` are not modelled. (ii) EMBOSS's one-A-per-GATC marking and masked EcoKI spacer are
reproduced, not corrected. (iii) CpG methylation is not in `Emethylsites.dat` and not offered. (iv) `PartialDigest` takes the set since B08 round 3 G3 (F83; Evidence RESTR-DIGEST-001 §12) and `FindAllSites(EnzymeSet)`, the raw-string
`FindSites` and `FindSites(…, enzymes, methylases)` since G4 (F84, below). `FindPatternOccurrences` does not (no reference: EMBOSS `fuzznuc` has no
methylation option); the `EnzymeSet` listing members never scan a sequence.

### 7.1 Every remaining scanning entry point (B08 round 3 G4 / F84, 2026-10-21)

Added: `FindAllSites(sequence, topology, EnzymeSet, methylases)`, raw-string `FindSites(string, enzymeName, topology, methylases)`,
`FindSites(sequence, topology, IEnumerable<RestrictionEnzyme>, methylases)`; MCP `find_all_restriction_sites` `methylation`. All reuse
`MethylationMask.Compute` and the `blocked` argument of `FindSitesCore` (nothing copied); `FindAllSites` computes the mask once per
molecule. Audit of every public member that reaches `FindSitesCore` / `CollectForwardCuts` / `DigestCore`: carried now = `FindSites` (all
seven forms), `FindAllSites`, `Digest`, `GetDigestSummary`, `CreateMap`, `PartialDigest`; the rest never reads a sequence
(`GetEnzymes*`, `GetBluntCutters`, `GetStickyCutters`, `FindCompatibleEnzymes`, `AreCompatible`, `FindIsoschizomers` ...) or scans a bare
pattern (`FindPatternOccurrences`, also used by `CodonOptimizer`): EMBOSS has no methylation mode for pattern scanning (`fuzznuc`,
`fuzzpro`, `patmatmotifs` list none; `restrict`, `remap`, `showseq`, `restover` do), so there is nothing to confirm a pattern mask against.

**Differential (executed `restrict`, 2026-10-21, harness outside the repo):** 40 seeded random 500-bp sequences with injected Dam / Dcm /
EcoKI / enzyme motifs, enzyme table = all 754 enzymes generated from Biopython 1.88, EMBOSS run once per sequence and mode (baseline and
`-methylation`), compared with `FindAllSites(seq, topology, EnzymeSet.Rebase)` and `...(…, MethylaseSets.Emboss)`. Sites present in both
baselines (the "common sites" of F68; the library additionally lists Type IIS windows EMBOSS does not and drops none): linear 39 058 common,
9 164 removed by EMBOSS, **0 mismatching sequences / 0 differing sites**; circular 39 537 common, 9 229 removed, 33 differing sites in 2
sequences, **all** library-removed / EMBOSS-kept and all caused by a GATC spanning the origin (the documented deliberate difference of F68).
(An earlier run that wrote one FASTA per sequence for several concurrent jobs produced spurious whole-sequence mismatches; the harness now
uses one file per job.) Locked values: the 39 common enzymes on L1, L2, L3 and a seeded 500-bp sequence under none / dam / dcm / emboss
(16 rows, EMBOSS = library, linear = circular); e.g. L2: none `DpnI:10,36;EcoRI:0,18;HpaII:28;MboI:10,36;MspI:28;Sau3AI:10,36;SmaI:27;TaqI:38,44`,
Dam / emboss `EcoRI:0,18;HpaII:28;MspI:28;SmaI:27;TaqI:38,44`. A raw-string case with sequence `N`s (`GGCCNNNNNGGCCAAGATCAAGAATTCAA`): MboI 15 removed, EcoRI 21 kept
as EMBOSS; SfiI 0 stays (Biopython's `N` rule, F40 - EMBOSS never matches a sequence `N`, a pre-existing difference).
**Locked by tests:** `RestrictionAnalyzer_MethylationEntryPoints_Tests` (G4 region, 24 cases) and MCP `RestrictionMethylationToolsTests`
(`FindAllRestrictionSites_*`).
