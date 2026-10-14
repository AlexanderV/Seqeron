# B08 — Source provenance and reachability

**Unit:** WP-SRC (source-provenance hardening), batch B08, 2026-10-14.
**Why this file exists:** B08's restriction findings (the 39/39 enzyme-table audit, the derived
ligation-compatibility rule, the overhang sequences) and two audit items (D7, D11) rest on
downloaded or vendored copies of REBASE. This file pins exactly which bytes were used, where an
independent copy can be had without network access to REBASE itself, and which claims are
executed, which are source-derived, and which cannot be verified in this environment at all.
It is referenced from `RESTR-FIND-001-Evidence.md`, `RESTR-FILTER-001-Evidence.md` and
`RESTR-DIGEST-001-Evidence.md`.

---

## 1. Host reachability

**Denied by organization egress policy** — the agent proxy answers **403 to CONNECT**
(`kind: connect_rejected`, `detail: "gateway answered 403 to CONNECT (policy denial or upstream
failure)"`). Verbatim from `curl -sS "$HTTPS_PROXY/__agentproxy/status"` → `recentRelayFailures`,
read 2026-10-14 (this unit did **not** re-attempt any denied host — see §6):

```
2026-10-14T17:57:47.955Z connect_rejected en.wikipedia.org:443
2026-10-14T17:57:48.248Z connect_rejected rebase.neb.com:443
2026-10-14T17:57:48.574Z connect_rejected www.addgene.org:443
2026-10-14T19:18:52.434Z connect_rejected rebase.neb.com:443
2026-10-14T19:18:52.700Z connect_rejected en.wikipedia.org:443
2026-10-14T19:18:53.016Z connect_rejected www.addgene.org:443
2026-10-14T19:18:56.700Z connect_rejected rebase.neb.com:443
2026-10-14T20:03:00.493Z connect_rejected www.neb.com:443
```

So: `rebase.neb.com`, `www.neb.com`, `en.wikipedia.org`, `www.addgene.org` are **unreachable as a
matter of policy, not of availability**. The most recent denial is 20:03Z, ~1.5 h before this unit
ran, which is why the policy is treated as still in force without re-triggering it.

**Allowed and used:** `raw.githubusercontent.com` (all dataset fetches below), and the GitHub code
*search* API through the repo's MCP GitHub tool (used only to *locate* candidate mirrors).

**Allowed but separately gated:** the GitHub *contents/commits* REST API. `gh api
repos/MartinGuehmann/diffdigester` and `mcp__github__get_file_contents` for
`crisprVerse/crisprBase` both return *"GitHub access to this repository is not enabled for this
session"*. This is a **session repository allowlist** (only `alexanderv/seqeron` is attached), not
an egress denial, and it is not worked around: every third-party file below was fetched from
`raw.githubusercontent.com` instead, which serves file content without repo attachment. The
consequence is that **commit dates and repository metadata for the third-party mirrors could not be
read**, so mirror recency is argued from file content (§4) rather than from git history.

---

## 2. Dataset 1 — the REBASE snapshot the findings actually rest on

| Property | Value |
|---|---|
| Installed package | `biopython` **1.88** (`Bio.__version__ == '1.88'`) |
| File | `/usr/local/lib/python3.13/dist-packages/Bio/Restriction/Restriction_Dictionary.py` |
| sha256 | `2a79099295dbad6061ea67a11e053787c591fcb2eb10fc8c0f89ead908dfa02b` |
| Size / lines | 562 312 B / 24 296 lines |
| Version header (verbatim, module docstring) | `Used REBASE emboss files version 404 (2024).` |
| Enzyme count | **1088** (`len(rest_dict)`); `len(typedict) == 51` |
| Per-enzyme fields | `charac, compsite, dna, freq, fst3, fst5, id, inact_temp, opt_temp, ovhg, ovhgseq, results, scd3, scd5, site, size, substrat, suppl, uri` |

REBASE's own version scheme is `<year-in-decade><month>`: `404` = **April 2024**, consistent with
`104` = April 2021 and `304` = April 2003 in the independently-headered files of §3.

**Reproducing it without network access:** `pip install biopython==1.88` (`pypi.org` and
`files.pythonhosted.org` are in the proxy's `noProxy` list, i.e. directly reachable) and hash the
file above.

**Independent corroboration from a second allowed location — exact match.**

| | |
|---|---|
| URL | `https://raw.githubusercontent.com/biopython/biopython/biopython-188/Bio/Restriction/Restriction_Dictionary.py` (pinned tag) |
| Result | **200, 562 312 B** |
| sha256 | `2a79099295dbad6061ea67a11e053787c591fcb2eb10fc8c0f89ead908dfa02b` |

The pinned-tag file is **byte-identical** to the installed file (`cmp` reports no differences; the
sha256 values are equal). No field-level diff was necessary — whole-file identity is strictly
stronger. The generated header matches verbatim, including the REBASE version line. The installed
`Restriction_Dictionary.py` is therefore exactly the file at tag `biopython-188`, and nothing in the
local environment has modified it.

**Important limitation of this dataset:** `rest_dict` carries **no organism field and no methylation
field** (see the field list above). It cannot corroborate audit items D11 (organism) or D7
(methylation) at all — hence §3.

---

## 3. Dataset 2 — REBASE `allenz`, for the organism and methylation columns

Two independent copies were obtained. Only the second carries a readable REBASE version header, so
it is the one the Evidence cites.

### 3a. Primary (versioned): REBASE `allenz.104` vendored in the Bioconductor package `crisprBase`

| | |
|---|---|
| Header URL | `https://raw.githubusercontent.com/crisprVerse/crisprBase/master/inst/enzymes/notation.txt` — **200, 4 931 B**, sha256 `3820fd9dc95c713ea706f2da9c5000975f6f6ab66125d6181448cd562128a64b` |
| Data URL | `https://raw.githubusercontent.com/crisprVerse/crisprBase/master/inst/enzymes/allenz.txt` — **200, 982 400 B**, sha256 `43b12432ed3c1a19e961d019e3a131336e8ee690a4568a68ddf55ed0e3869388` |
| Records | **5 322** enzymes (`grep -c '^<1>'`) |
| Format | REBASE `<n>`-tagged records; the header lives in `notation.txt`, the records in `allenz.txt` |

Header, verbatim (`notation.txt` lines 1-9):

```
REBASE version 104                                              allenz.104
 
    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    REBASE, The Restriction Enzyme Database   http://rebase.neb.com
    Copyright (c)  Dr. Richard J. Roberts, 2021.   All rights reserved.
    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 
Rich Roberts                                                    Mar 31 2021
```

This is the genuine REBASE `allenz` file header: named version (`104` = April 2021), named release
date, and Rich Roberts' copyright. `notation.txt` also documents the column semantics used below —
`<1>` name, `<2>` prototype, `<3>` microorganism, `<4>` source, `<5>` recognition sequence with
cleavage site, `<6>` methylation site, `<7>` commercial availability, `<8>` references.

### 3b. Secondary (unversioned): the `diffdigester` mirror the round-1 auditor used

| | |
|---|---|
| URL | `https://raw.githubusercontent.com/MartinGuehmann/diffdigester/master/allenz.commercial.txt` |
| Result | **200, 49 224 B** |
| sha256 | `c74e635890c34db395cf6974ffd2aa6edce94cd07b836f2db0e2d5d3798196cb` |
| Records | **622** enzymes, tab-separated, all 622 lines with the same 8-column shape |
| REBASE version header | **ABSENT** |
| Roberts copyright | **ABSENT** |

**Its self-declared provenance is nil.** The file is not a verbatim REBASE download: it begins
`database = ` + a backtick and ends with a backtick + `;`, i.e. it is the allenz payload embedded as
a source-code string literal, with REBASE's entire header block (version, release date, copyright,
column documentation, supplier key) stripped. `grep -iE 'rebase|roberts|copyright|version|allenz'`
matches **no header line** — only enzyme records. **The REBASE version of this file cannot be
determined from the file itself.**

Repository context (the one provenance signal that is available): the repo's `README.md`
(**200, 1 201 B**) identifies it as the data behind a peer-reviewed tool — *Martin Gühmann,
Stefanie Reuter, Ralf Mrowka, "Diffdigester.uni-jena.de: a tool for optimized selection of
restriction enzymes for plasmid identification in cloning procedures", Nucleic Acids Research, 2025,
gkaf418*, DOI `10.1093/nar/gkaf418`, software DOI `10.5281/zenodo.15207957`. That makes it a
credible research artifact rather than an anonymous dump, but it is still an **undated, unversioned,
preprocessed, commercial-only subset**.

### 3c. Further independent copies found (REBASE `withrefm`, same 8 columns, with full headers)

| URL | Status / bytes | sha256 | Header version | Records |
|---|---|---|---|---|
| `raw.githubusercontent.com/bebop/poly/main/io/rebase/data/rebase_test.txt` | 200, 53 648 B | `1976c7e805102c394ee1ede0ca474f7c8e75fa49d5c53c5d481284444af1ad40` | `REBASE version 104 / withrefm.104`, Roberts 2021, Mar 31 2021 | 92 |
| `raw.githubusercontent.com/kimrutherford/EMBOSS/master/test/data/withrefm` | 200, 25 336 B | `998af2fc0a89e1c1c620c5815e115625beed409f7bd32e5f46e11d2468a7c3c6` | `REBASE version 106 / withrefm.106`, Roberts 2001, May 31 2001 | 33 |
| `raw.githubusercontent.com/lauringlab/CodonShuffle/master/lib/EMBOSS-6.6.0/test/data/withrefm` | 200, 25 336 B | (identical to the row above) | same | 33 |
| `raw.githubusercontent.com/gitpan/BioPerl/master/t/data/rebase.withrefm` | 200, 10 405 B | `8ee016da79342e27936b92b8c07de64dd951e510a79e6b2ffa2b7ec91a07eb1b` | `REBASE version 304 / withrefm.304`, Roberts 2003, Mar 24 2003 | 12 |

These are **test fixtures** (12-92 enzymes), not full databases, so they are recorded as
corroboration that the REBASE header format and column semantics are as stated — not as a data
source. The EMBOSS-lineage and BioPerl-lineage copies carry the identical header block, which is
what establishes that §3b's file is missing one rather than that REBASE never had one.

---

## 4. Mirror cross-check numbers

### 4a. `diffdigester` (§3b) vs REBASE `allenz.104` (§3a) — the fields D11 and D7 depend on

All **622** enzymes of the diffdigester file are present in `allenz.104`'s 5 322 records.

| Column | Agreement |
|---|---|
| `<3>` **Microorganism** (D11) | **622 / 622** exact string match |
| `<6>` **Methylation site** (D7) | **622 / 622** exact string match |
| `<2>` Prototype | **622 / 622** exact string match |
| `<5>` Recognition sequence + cut | **621 / 622** — one difference, `AbaSI` |
| `<7>` Commercial availability | **590 / 622** — 32 differences |

**The one site difference** is `AbaSI`: diffdigester `CNNNNNNNNNNNNNNNNNNNNG`, `allenz.104`
`C(11/9)`. Biopython 404 agrees with `allenz.104` (`site='C'`, `fst5=12`, `fst3=9`). So this is
diffdigester's own expansion of an offset-notated site into explicit `N`s — plausible for a tool
that compiles sites to regexes. `AbaSI` is **not** one of the 39 built-ins and no B08 finding
touches it.

**The 32 supplier differences all point one way:** diffdigester carries an *extra* supplier code
that both `allenz.104` (2021) and Biopython 404 (2024) have dropped — e.g. `AccI` `BJKMNQRX` vs
`BJKMNQX`, `AccIII` `JKR` vs `JK`, `BanI` `NR` vs `N`, `EcoICRI` `IRV` vs `IV` (mostly `R` =
Promega, some `S` = Sigma). `allenz.104` and Biopython 404 agree with each other exactly on these.
Since supplier listings are only ever withdrawn over time, **the diffdigester file is older than
REBASE 104 (April 2021)**. Its exact version remains unknown because its header was stripped.

### 4b. `diffdigester` (§3b) vs Biopython 404 (§2) — sites and cut positions

620 enzymes in common (diffdigester-only: `CviKI-1`, which Biopython keys as `CviKI_1`, and
`EcoP15I`, which Biopython's `rest_dict` does not carry at all — naming/inclusion, not data,
differences).

| Comparison | Result |
|---|---|
| Recognition site string | **616 / 620** agree |
| Cut positions (top, bottom), where expressible | **605 / 608** agree (12 not expressible: 11 two-sided cutters, 1 with no cut information) |

**Every one of the 7 differences is accounted for, and none is a data conflict:**

- `SgeI`, `TscAI`, `TspRI` — Biopython **strips trailing `N`s** from the site string. Verified
  exactly: `allenz_site.rstrip('N') == biopython_site` for all three (`CASTGNN` → `CASTG`,
  `CNNGNNNNNNNNN` → `CNNG`). Their *top-strand* cut positions agree in both sources (7 and 7, 13 and
  13); only the derived bottom-strand position differs, because it is computed from the site length
  that Biopython shortened. These same three therefore account for 3 of the 4 site differences and
  all 3 cut differences.
- `AbaSI` — the fourth site difference, which is §4a's diffdigester-side expansion, not a Biopython
  convention.

### 4c. The 39 built-ins

All 39 entries of `RestrictionAnalyzer._enzymes` are present in **both** allenz copies, and the two
copies give **identical** organism and methylation strings for all 39.

- **Organism (D11): 39 / 39** agree at genus + species against both copies. REBASE's string is
  consistently the Seqeron string plus a strain suffix — e.g. `EcoRI` Seqeron *Escherichia coli* vs
  REBASE *Escherichia coli RY13*; `HindIII` *Haemophilus influenzae* vs *Haemophilus influenzae Rd*;
  `ApaI` *Acetobacter pasteurianus* vs *Acetobacter pasteurianus sub. pasteurianus*. The
  round-1 auditor's 39/39 claim is **independently confirmed**, now on two mirrors.
- **Methylation (D7): 39 / 39** identical between the two copies (e.g. `EcoRI` `3(6)`, `BamHI`
  `5(4)`, `AluI` `3(5)`; `AvrII`, `DpnI`, `PacI`, `SphI`, `StuI` empty). The premise of D7 — that
  REBASE field 6 is obtainable here — holds.

---

## 5. Verdict on the `allenz` mirror (audit items D7 and D11)

**YES — trustworthy enough to source D7 and D11 from, but only on the following basis, and the
primary citation must be `allenz.104`, not the diffdigester mirror.**

Reasoning, stated without softening in either direction:

1. **On its own the diffdigester mirror would not be sufficient.** It carries no REBASE version, no
   release date and no copyright line; it is a code-embedded, header-stripped, commercial-only
   subset; at least one site record is reformatted; and it is demonstrably *stale* relative to
   REBASE 104 (§4a). A finding whose only citation was this file would be uncheckable and
   silently out of date.
2. **It is nevertheless genuine REBASE content, not corrupted or fabricated.** Against a second,
   independently-vendored copy that *does* carry the verbatim REBASE header, it agrees 622/622 on
   organism, 622/622 on methylation and 622/622 on prototype, and 621/622 on the site field. Against
   a third, fully independent lineage (Biopython's EMBOSS-derived 404 snapshot) its sites and cut
   positions agree for every enzyme once Biopython's documented trailing-`N` convention is applied.
3. **Therefore the data is corroborated and the citation is upgraded.** D7 and D11 are **not**
   blocked-pending-network. They are source-derived from a versioned REBASE file obtained on an
   allowed host.

**What the Evidence must cite for D7 / D11:**

- Primary: **REBASE `allenz.104`** (REBASE version 104, Copyright Dr. Richard J. Roberts 2021,
  released Mar 31 2021), vendored in `crisprVerse/crisprBase`:
  `inst/enzymes/notation.txt` (sha256 `3820fd9dc95c713ea706f2da9c5000975f6f6ab66125d6181448cd562128a64b`)
  + `inst/enzymes/allenz.txt` (sha256 `43b12432ed3c1a19e961d019e3a131336e8ee690a4568a68ddf55ed0e3869388`,
  5 322 records), fields `<3>` organism and `<6>` methylation.
- Corroborating mirror: `MartinGuehmann/diffdigester` `allenz.commercial.txt`
  (sha256 `c74e635890c34db395cf6974ffd2aa6edce94cd07b836f2db0e2d5d3798196cb`, 622 records,
  **no version header**, older than REBASE 104), NAR 2025 `10.1093/nar/gkaf418`.
- Cross-check numbers: organism 622/622, methylation 622/622, prototype 622/622, site 621/622
  (`AbaSI` only), suppliers 590/622; and 39/39 for the built-ins on organism (genus+species) and on
  methylation.
- Both are **source-derived**, not executed: REBASE is a database, so these columns are read, not
  computed. No test asserts an organism or methylation string as a computed result.

---

## 6. No circumvention was attempted

This unit used only hosts that are *allowed*. It did **not** attempt, and does not recommend, any
proxy or VPN, alternate port, IP literal, text-extraction relay (`r.jina.ai` and similar), search-engine
cache, or any other route to `rebase.neb.com`, `www.neb.com`, `en.wikipedia.org` or
`www.addgene.org`. It did not re-request those hosts at all: the denial record in §1 is read from
the proxy's own status endpoint, which is why the timestamps there are from earlier units rather
than from this one. Per `/root/.ccr/README.md`, organization policy denials are reported, not
routed around; TLS verification was never disabled and `HTTPS_PROXY` was never unset.

---

## 7. Classification of every claim B08's restriction findings rest on

### 7a. Reference-executed (a reference implementation was run locally and its output compared)

Biopython 1.88 `Bio.Restriction` is installed, so these are executed, not read:

- the 39/39 enzyme-table audit on recognition sequence, `CutPositionForward` vs `fst5`,
  `CutPositionReverse` vs `len(site) + fst3`, and overhang class vs `is_blunt()` /
  `is_5overhang()` / `is_3overhang()`;
- the site-search comparisons (34 906 comparisons, 0 mismatches, both topologies);
- the overhang-sequence agreement with `ovhgseq` (39/39);
- the ligation-compatibility cross-check against Biopython's `%` operator.

### 7b. Source-derived (read from a pinned dataset; never presented as executed)

- Everything taken from the **REBASE 404 snapshot** as data rather than behaviour: enzyme counts,
  `suppl`, REBASE ids/URIs, and the population statistics quoted in B08.md and in the
  `RestrictionAnalyzer` comments.
- **Organism (D11)** and **methylation (D7)** from REBASE `allenz.104` (§5). No algorithm reads
  either column.

All numeric source-derived claims in the restriction code comments were re-verified against the
pinned snapshot by this unit, and **all of them hold**:

| Claim (location) | Recomputed |
|---|---|
| 1088 REBASE enzymes | **1088** ✅ |
| "REBASE has exactly two non-palindromic in-site cutters, `CdiI` `CATC^G` and `SimI` `GG^GTC_N`" (`RestrictionAnalyzer.cs:163`) | exactly **2**: `CdiI`, `SimI` ✅ |
| "389 REBASE enzymes whose `ovhgseq` is blunt or fully A/C/G/T, all 75 466 unordered pairs" (`:1124`) | **389**, and 389·388/2 = **75 466** ✅ |
| "asymmetric for 2386 of 283 881 REBASE pairs" (`:1127`) | the pair count matches the **754** enzymes with a non-`None` `ovhgseq`: 754·753/2 = **283 881** ✅ |
| 39 built-ins: 10 blunt, 21 5'-overhang, 8 3'-overhang (`:118`, `:130`) | **10 / 21 / 8** ✅ (sums to 39) |
| "39/39 agreement" against Biopython `.size` (`:92`) | **39/39** ✅ |
| D4: "1088 enzymes, **576** expressible in the current record, **252** with a fully defined overhang" (B08.md:516) | **576** expressible; of those, **252** have a non-blunt, fully A/C/G/T `ovhgseq` ✅ |

**No misattributed citation was found in the restriction code.** (The WP1-class defect — a comment
crediting `Tm_NN(DNA_NN3)` where upstream Azimuth calls the obsolete `Tm_staluc` — is already fixed:
`AzimuthRuleSet2.cs:308-313` now states that upstream calls `MeltingTemp.Tm_staluc` and documents
the exact equivalence to `Tm_NN`.)

### 7c. Unverifiable in this environment (citation cannot be opened by anyone in this campaign)

These claims are **not withdrawn** — each is either definitional or independently backed by an
executed reference — but **their cited source was never opened by this campaign**, and a reader must
not read them as checked. Annotated here rather than deleted.

| Where | Citation | Status |
|---|---|---|
| `tests/TestSpecs/RESTR-FILTER-001.md` rows 1-3 | `en.wikipedia.org` — *Sticky and blunt ends*, *Restriction enzyme*, *List of restriction enzyme cutting sites* (dated 2026-06-13, pre-campaign) | **Unverifiable here** (host denied). The *substance* is re-sourced: the blunt / 5' / 3' definitions are executed against Biopython `is_blunt()` / `is_5overhang()` / `is_3overhang()` (39/39), and the cutting-site table is superseded by REBASE 404 (§2) and `allenz.104` (§3a) |
| `tests/TestSpecs/RESTR-FILTER-001.md` rows 4-5 | `www.neb.com` product pages R0142 (KpnI), R3101 (EcoRI-HF) | **Unverifiable here** (host denied). The site / cut / overhang data for KpnI and EcoRI is re-sourced from REBASE 404 and `allenz.104`; the product-specific content (HF variant, buffers) is **not** re-sourced and nothing in B08 depends on it |
| `tests/TestSpecs/CRISPR-GUIDE-001.md` row "Addgene CRISPR Guide" | `www.addgene.org/guides/crispr/` | ~~**Unverifiable here** (host denied)~~ → **VERIFIED 2026-10-14 (WP3)** — seed range, spacer/scaffold, NGG/NNGRRT/TTTV PAMs confirmed; the Pol III `TTTT` claim is **not on the page**. See §8 |
| `tests/TestSpecs/CRISPR-GUIDE-001.md` rows "Wikipedia: Guide RNA", "Wikipedia: PAM" | `en.wikipedia.org` | **Unverifiable here** (host denied). PAM definitions were separately re-confirmed against CRISPOR's tables (`CRISPR-PAM-001-Evidence.md`) |
| `tests/TestSpecs/CRISPR-OFF-001.md` row "Wikipedia: Off-target genome editing" | `en.wikipedia.org` | **Unverifiable here** (host denied) |
| `RestrictionAnalyzer.cs:85`, `:105` | Wikipedia *Restriction enzyme* — "Type II sites are 4-8 nt" | **Unverifiable here.** Descriptive range only; no behaviour depends on it. The same remarks block already notes the interrupted-palindrome exception (`SfiI`, `size == 13`) measured against Biopython `.size` |
| `RestrictionAnalyzer.cs:115`, `:126`, `:1329` | Wikipedia *Sticky and blunt ends* | **Unverifiable here**, but each is paired in the same comment with an executed Biopython backing (`is_blunt()`, `Blunt._mod1`) |
| `RestrictionAnalyzer.cs:765`, `:983`, `:1497`, `:1559` | Addgene *Plasmids 101* / *Restriction Digest Protocol*, Univ. of Illinois MolBio | **Unverifiable here.** Both claims they support (a circle with *k* cuts gives *k* fragments; fragment lengths sum to the molecule length) are arithmetic identities that the digest tests assert directly, and `:765` additionally cites Biopython `catalyse(..., linear=False)` |
| ~~`CrisprDesigner.cs:258`~~ | ~~"Evidence: Addgene — seed sequence is 8-10 bases at 3' end; using upper bound (10bp)"~~ | **Superseded — VERIFIED 2026-10-14 (WP3, §8).** Original annotation: **Unverifiable here** (host denied): the 8-10 nt range cannot be checked. Not changed and not a misattribution — the comment states its own derivation (upper bound of the cited range), and that is the same convention `CRISPR-OFF-001.md:38` applies to the *off-target* seed (12 bp = upper bound of Hsu 2013's 8-12 bp, a reachable primary source). The guide-design seed (10) and the off-target seed (12) are different parameters resting on different sources; only the former's source is denied |
| `CRISPR-GUIDE-001-Evidence.md:43` | `bitbucket.org/valenlab/chopchop/raw/HEAD/chopchop.py` — **404** | Still 404; this is a *missing path*, not a policy denial. CHOPCHOP's self-complementarity measure was sourced from the reachable copy recorded in that Evidence file |
| `CRISPR-OFF-001-Evidence.md:29`, `:95` | `api.github.com/repos/maximilianh/crisporWebsite/contents/bin/src` — **403** | Still unavailable, and now explained: this is the **session repository allowlist** (§1), not egress policy. The directory cannot be listed from here; `raw.githubusercontent.com` serves individual files but cannot enumerate a directory |

### 7d. Items previously recorded as "not verified" that are now verified

| Previously | Now |
|---|---|
| `B08.md:145`, `RESTR-FIND-001-Evidence.md:110` and its §4 table: *"each entry's `Organism` string — the two sources that carry it are the blocked hosts and Biopython's snapshot has no organism field… Unverified."* | **Verified 39/39** against REBASE `allenz.104` on an allowed host, corroborated by a second mirror (§4c). The premise was correct that Biopython has no organism field; what was missed is that `allenz` is vendored by third-party packages on allowed hosts |
| D7: *"Methylation data claimed unobtainable"* | REBASE field `<6>` **is** obtainable here, and the two mirrors agree 622/622 (39/39 for the built-ins) |
| *"REBASE release newer than 404 — same blocked hosts"* | **Still true.** 404 (April 2024) remains the newest release reachable from here; `allenz.104` (April 2021) is *older*, and is used only for the two columns the 404 snapshot omits |

## 8. Addgene seed-range citation — verified (B08 WP3, 2026-10-14)

`www.addgene.org` now answers 200 through the agent proxy (the §1 denials are superseded for this
host). Fetched with `curl -sS -L` (TLS verified, proxy unchanged) on **2026-10-14**:

| Field | Value |
|---|---|
| URL | https://www.addgene.org/guides/crispr/ (HTTP 200, no redirect, 392 036 bytes) |
| SHA-256 of the HTML as fetched | `23226099b62de0567dd0a842e61743dc7db622b338e5c402fd3e3d85670fe7f0` |
| Section | *CRISPR Basics* → *Basic Mechanisms of CRISPR* |
| Exact sentence | "Once the Cas9-gRNA complex binds a putative DNA target, the seed sequence (8–10 bases at the 3′ end of the gRNA targeting sequence) will begin to anneal to the target DNA." |
| Following sentences | "If the seed and target DNA sequences match, the gRNA will continue to anneal to the target DNA in a 3′ to 5′ direction. … Mismatches between the target sequence in the 3′ seed sequence inhibit target cleavage, while mismatches toward the 5′ end distal to the PAM often permit target cleavage." |

**Verdict: the citation is CONFIRMED as written** — 8-10 bases, at the 3' end of the targeting
sequence; `CrisprDesigner.EvaluateGuideRna` uses the upper bound, 10. No value changed; the code
comment (`CrisprDesigner.cs`, seed-region block) now carries the quote, URL and fetch date. Scope
caveat recorded in the comment: the sentence is about **Cas9**; the PAM-before-target (Cas12a)
branch's "first 10 nt" mirrors it at the PAM-proximal end and is not itself stated on this page.

Other `CRISPR-GUIDE-001.md` claims attributed to this page, checked on the same fetch: the ~20-nt
spacer + Cas-binding scaffold ("a user-defined ∼20-nucleotide spacer that defines the genomic
target"), SpCas9 `NGG`, SaCas9 `NNGRRT` and Cas12a `TTTV` (PAM table and Cas12a section) — **present**.
The Pol III `TTTT` terminator — **absent** (no "TTTT", "Pol III", "U6" or "terminat" on the page):
that attribution remains unverified, and is the only Addgene-Guide row not closed here.
