---
type: source
title: "Evidence: ONCO-DRIVER-001 (driver-gene classification — 20/20 rule)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-DRIVER-001-Evidence.md
sources:
  - docs/Evidence/ONCO-DRIVER-001-Evidence.md
source_commit: f640eb404dd41ebb270208e6664505bba6c4cb8e
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-DRIVER-001

The validation-evidence artifact for test unit **ONCO-DRIVER-001** — **Driver Mutation Detection
(20/20 rule)**. This is the **twelfth ingested unit of the Oncology family** and one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct
mutation-pattern classification method is synthesized in its own concept,
[[driver-gene-classification-20-20-rule]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (all mutually consistent, one glyph difference):**
  - **Vogelstein B et al. (2013)** *Cancer Genome Landscapes*, *Science* 339(6127):1546–1558 (rank 1,
    **originating source** of the 20/20 rule; PMC full text behind CAPTCHA and the DOI HTTP 403 at fetch
    time, so the exact wording was taken from the three open-access secondary sources below that quote the
    primary text verbatim). Rule: an **oncogene** has **> 20%** of recorded mutations at **recurrent
    positions and missense (activating)**; a **tumor-suppressor gene** has **> 20%** of mutations
    **inactivating (truncating / loss of function)**. The rule is **lenient** — all well-documented cancer
    genes far surpass it. **IDH1** worked example: nearly all mutations at the identical amino acid
    (codon 132, R132H) → classified an **oncogene**.
  - **Tokheim & Karchin (2020)** *Somatic selection distinguishes oncogenes and tumor suppressor genes*
    (20/20+), *Bioinformatics* 36(6):1712 (rank 1, quotes the rule verbatim) — "OGs have >20% mutations
    causing **missense changes at recurrent positions** and TSGs have >20% mutations causing
    **inactivating changes**." Inactivating = protein-truncating: **nonsense and frameshift**.
  - **Schroeder et al. (2014)** *OncodriveROLE* (rank 1, reference implementation) — restates the rule
    and gives the **explicit truncating list**: "truncating mutations include mutations causing a
    **frameshift, a gained or lost stop codon** as well as mutations in **splice donor or acceptor
    sites**." Writes "≥20%" for TSGs (the one glyph difference vs the primary source's strict ">20%").
  - **Miller et al. (2017)** *Identification and analysis of mutational hotspots* (rank 1) — restates the
    rule, defines a **recurrent / hotspot position** as "at least two mutations of the same class" at an
    identical location (**same protein position observed ≥ 2 times**), confirms **IDH1 codon 132
    (R132H)**, and lists truncations as "nonsense mutations and frameshift insertions and deletions."

- **Documented corner cases / failure modes:**
  - **Passenger truncations in oncogenes** (Tokheim 2020): random truncating mutations drift up in
    frequency with no fitness impact and can mislead — the rule is a **heuristic, not a statistical
    test**, and can misclassify genes with few mutations.
  - **Low-recurrence drivers** (OncodriveROLE): the rule "fails to identify … the lowly recurrent ones";
    a gene meeting neither criterion is left **Ambiguous** (correct conservative behaviour).
  - **Single-amino-acid recurrence** (IDH1): when virtually all missense fall on one codon, the
    recurrent-missense fraction is ~1.0 ≫ 0.20 → a clear oncogene call.

- **Datasets (documented oracles):**
  - **IDH1** — 10 missense all at codon 132: truncating **0/10 = 0.00**, recurrent-missense **10/10 =
    1.00** → **Oncogene**.
  - **TSG archetype** — 5 nonsense + 2 frameshift + 1 splice = 8 truncating (distinct positions) + 2
    missense (distinct): truncating **8/10 = 0.80**, recurrent-missense **0.00** → **TumorSuppressor**.
  - **Boundary** — truncating fraction **exactly 0.20 → NOT a TSG** (strict `>`); **0.30 →
    TumorSuppressor**.

## Deviations and assumptions

- **ASSUMPTION 1 — dual-pass tie-break.** If a gene satisfies **both** the OG (>20% recurrent missense)
  and TSG (>20% truncating) criteria, the sources prescribe no single label. The unit classifies by the
  **larger fraction** and reports **Ambiguous on an exact tie**. Justification: well-documented genes
  "far surpass" one criterion, so a dual pass is atypical; choosing the dominant signal is the
  least-surprising deterministic resolution. The two archetypes are unaffected.
- **ASSUMPTION 2 — strict `> 0.20` threshold.** Vogelstein / Tokheim / Miller say ">20%"; OncodriveROLE
  writes "≥20%" for TSGs. The unit uses **strict `>` 0.20 for both**, matching the primary source and the
  verbatim Tokheim quote; a fraction of exactly 0.20 is not sufficient.
- **ASSUMPTION 3 — `ScoreDriverPotential` proxy.** The checklist names CADD/SIFT/PolyPhen, but those are
  externally trained models that cannot be retrieved/reproduced (forbidden to fabricate). The unit
  returns the **20/20-rule driver-signal fraction in [0,1] = max of the two criterion fractions** as a
  transparent, source-derived score, documenting that external pathogenicity scores are caller-supplied /
  not implemented.
- **Coverage recommendations:** MUST-test the IDH1 all-recurrent-missense → Oncogene (fraction 1.00),
  dispersed-truncating → TumorSuppressor (fraction > 0.20), the boundary (exactly 0.20 → not TSG), and
  `MatchCancerHotspots` flagging a (gene, position) in the caller-supplied hotspot set. SHOULD-test
  neither-criterion → Ambiguous and `IdentifyDriverMutations` returning a subset of input variants
  (driver ⊆ somatic). COULD-test that a singleton missense is not recurrent (recurrence needs the same
  position ≥ 2 times).

No source contradictions — the four references are mutually consistent; the sole `≥` vs `>` glyph
difference (OncodriveROLE) is resolved in favour of the primary Vogelstein / Tokheim strict ">20%".
