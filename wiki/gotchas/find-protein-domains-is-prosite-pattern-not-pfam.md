---
type: gotcha
title: "find_protein_domains matches exact PROSITE patterns + 3 bundled profiles — not a full Pfam/HMMER scan"
tags: [protein-features, gotcha]
mcp_tools:
  - find_protein_domains
sources:
  - docs/algorithms/ProteinMotif/Domain_Prediction.md
  - docs/algorithms/ProteinMotif/Profile_HMM_Domain_Detection.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# find_protein_domains matches exact PROSITE patterns + 3 bundled profiles — not a full Pfam/HMMER scan

**The trap.** `FindDomains` detects domains that have an **exact PROSITE PATTERN** (a regex-like
signature), plus a Plan7 profile-HMM path for which **only 3 CC0 profiles are bundled** (any other
family must be a caller-supplied `.hmm`). Its domain score is an **internal information-content
match-strength indicator, not a calibrated bit-score**, and standard HMMER MSV/bias prefilters are
not part of the reported score.

**Why it bites.** This is **not** a Pfam/InterPro-scale domain scanner. A protein whose domain has no
exact PROSITE PATTERN (most Pfam families) and no supplied `.hmm` returns **no hit** — absence here
means "not in the tiny built-in set", not "no domain". And the match-strength number is not
comparable to a HMMER bit-score or E-value, so ranking/ thresholding it like one is wrong.

**What to rely on instead.** Supply the relevant `.hmm` profiles, or run HMMER/InterProScan against
Pfam for genome-scale annotation; use the built-in PROSITE/CC0 path for the specific bundled
families. Full model: [[protein-domain-and-signal-peptide-prediction]]; simple fixed motifs are
[[common-protein-motifs]]. Research-grade, [[research-grade-limitations]].
