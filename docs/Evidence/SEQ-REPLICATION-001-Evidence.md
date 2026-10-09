# Evidence Artifact: SEQ-REPLICATION-001

**Test Unit ID:** SEQ-REPLICATION-001
**Algorithm:** Replication Origin Prediction (cumulative GC-skew minimum)
**Date Collected:** 2026-06-14

---

## Online Sources

### Rosalind — Minimum Skew Problem (BA1F)

**URL:** https://rosalind.info/problems/ba1f/
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference problem with an exact published worked example / canonical algorithm definition)

**Retrieved via:** WebFetch of the URL above.

**Key Extracted Points:**

1. **Skew definition:** "The skew of a DNA string Genome, denoted Skew(Genome), [is] the difference between the total number of occurrences of 'G' and 'C' in Genome." Computed by iterating each position: G contributes +1, C contributes −1, and A/T contribute 0, starting from 0 at position 0.
2. **Minimum Skew Problem:** "Find a position in a genome minimizing the skew. Given: A DNA string Genome. Return: All integer(s) i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)." So positions are 0-based prefix indices in [0, |Genome|].
3. **Sample input (verbatim):** `CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG` (length 100).
4. **Sample output (verbatim):** `53 97` — the positions minimizing the skew.

**Independent re-derivation (in this session):** running the per-nucleotide cumulative skew (G:+1, C:−1, A/T:0; Skew_0 = 0) over the sample input yields a global minimum value of −4 at prefix indices 53 and 97, reproducing the published output exactly.

---

### Grigoriev A (1998) — Analyzing genomes with cumulative skew diagrams

**URL:** https://academic.oup.com/nar/article/26/10/2286/1030593 (DOI: https://doi.org/10.1093/nar/26.10.2286)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary literature, Nucleic Acids Research)

**Retrieved via:** WebSearch ("Grigoriev 1998 Analyzing genomes with cumulative skew diagrams Nucleic Acids Research") then WebFetch of the Oxford Academic article page.

**Key Extracted Points:**

1. **Abstract (verbatim):** "A novel method of cumulative diagrams shows that the nucleotide composition of a microbial chromosome changes at two points separated by about a half of its length. These points coincide with sites of replication origin and terminus for all bacteria where such sites are known."
2. **Construction:** "a sum of (G−C)/(G+C) in adjacent windows from an arbitrary start to a given point in a sequence" (WebSearch snippet of the article, re-checked 2026-09-28; academic.oup.com blocked for WebFetch). With a one-base window each window's skew is +1 (G), −1 (C) or 0, so this is exactly the per-nucleotide #G−#C diagram of Rosalind BA1F; Seqeron computes the prediction from the canonical `CalculateCumulativeGcSkew` kernel at window 1.
3. **Origin/terminus location:** the cumulative GC-skew diagram reaches its global minimum at the replication origin and its global maximum near the terminus; the two extrema are separated by roughly half the chromosome length.
4. **Strand bias:** the leading strand contains more guanine than cytosine.

---

### Wikipedia — GC skew (used for its cited primaries)

**URL:** https://en.wikipedia.org/wiki/GC_skew
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia; relied on for the formula statement and its cited primaries Lobry 1996, Grigoriev 1998)

**Retrieved via:** WebSearch ("Lobry 1996 ... GC skew origin replication") then WebFetch of the article.

**Key Extracted Points:**

1. **Formula (verbatim):** "GC skew = (G − C)/(G + C)".
2. **Replication features (verbatim):** "the maximum value of the cumulative skew corresponds to the terminal, and the minimum value corresponds to the origin of replication."
3. **Strand bias (verbatim):** "the leading strand contains more guanine (G) and thymine (T), whereas the lagging strand contains more adenine (A) and cytosine (C)."
4. **Cited primaries:** Lobry, J. R. (1996) Mol Biol Evol 13:660–665; Grigoriev, A. (1998) Nucleic Acids Res 26:2286–2290.

---

### Lobry JR (1996) — Asymmetric substitution patterns in the two DNA strands of bacteria

**URL:** https://pubmed.ncbi.nlm.nih.gov/8676740/ (Mol Biol Evol 13(5):660–665)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary literature)

**Retrieved via:** WebSearch ("Lobry 1996 asymmetric substitution patterns ... GC skew origin replication"); PubMed record and the review summarizing it confirmed the points below.

**Key Extracted Points:**

1. Lobry first reported (1996) compositional asymmetry between the two strands in three bacteria (E. coli, B. subtilis, H. influenzae): departure from intrastrand A=T and C=G equifrequency.
2. GC/AT skews switch sign at the origin and terminus of replication; this sign switch is used to confirm/predict the origin of replication.

---

### Rosalind BA1F extra dataset (mirror)

**URL:** https://raw.githubusercontent.com/charlesreid1/go-rosalind/master/rosalind/data/minimum_skew.txt (Rosalind's own extra-dataset page needs a login; rosalind.info is blocked for curl)
**Accessed:** 2026-10-09 (curl; md5 7155ffe57770d00d51c215e9f0aaa5a9; stored verbatim as `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Rosalind/ba1f_extra_dataset.txt`)

**Key Extracted Points:** `Input` genome of 93 523 bp, `Output` `89969 89970 89971 90345 90346`. Python brute force (per-base walk, all argmin): min −184 at exactly those indices; max +154 at [20377, 20378, 20379]; Skew_n = −41.

---

### Lu J, Salzberg SL (2020) — SkewIT (Skew Index Test)

**URL:** https://doi.org/10.1371/journal.pcbi.1008439 (PLoS Comput Biol 16:e1008439; PMC7717575); code https://github.com/jenniferlu717/SkewIT
**Accessed:** 2026-10-09
**Authority rank:** 1 (peer-reviewed) + original reference implementation

**Retrieved via:** curl of `raw.githubusercontent.com/jenniferlu717/SkewIT/master/{README.md, src/skewi.py, src/gcskew.py, data/RefSeq97_Bacteria_GenusSkewIThresholds.txt}` (all HTTP 200). The paper itself (journals.plos.org, biorxiv, par.nsf.gov) returned proxy 403; WebSearch confirmed the bibliographic data.

**Key Extracted Points (from `skewi.py`):**

1. Windows `for i in range(0, len(seq), k)` over `seq[i:i+k]` (partial tail kept); each window scores `1/−1/0` by sign(`count("G") − count("C")`) (upper case only). Default k = 20000; `-f` is validated but not used in the computation.
2. `half_len = round(L/2)`, `curr_range = round(L*0.04)` (Python 3 round), `skew += skew[:L]`; `max_diff = max |sum(skew[i:t]) − sum(skew[t:i+L])|` over `i ∈ [0, L)`, `t ∈ [i+half−cr, i+half+cr)`.
3. `skewi = max_diff / len(seq) * k`, capped at 1.0; written only when `max_diff > 0`.
4. CLI filters: `--min-len` 500 000 default; header must contain "complete" (`--complete` default); "plasmid" excluded by default.
5. README: SkewI in [0, 1], higher = stronger GC-skew signal; per-genus thresholds = mean − 2 SD for genera with ≥ 10 genomes (RefSeq 97), listed in the data file (e.g. `g__Escherichia 934 0.8486 0.0688 0.7110`). No universal default cutoff.

**Reference run (skewi.py itself, Python 3 + Biopython 1.88, `--min-len 0`, `%0.10f` patched to `%r`):** `GGGGCCCC`×6+`GGGG` k4 → 0.23076923076923078; 58-mer `GCTAAAGACA…GTGAAT` k4 → 0.6206896551724138; 161-mer `GGGCCAATTG…CCTCCGCC` k4 → 0.6211180124223602; `G`×52 k4 → 0.2307…; `GGGGCCCC`×6 (12 windows) and `A`×52 → no output; `gggg`×13 → no output (lower case; Seqeron upper-cases → 0.2307…, documented deviation); BA1F extra dataset k1000 → 0.203158581311549, k20 → 0.034430033253851994; 100 kb synthetic genome k1000/k20 → 1.0, k20000 (5 windows) → no output.

**Threshold table and license (2026-10-09, A2-3):** `data/RefSeq97_Bacteria_GenusSkewIThresholds.txt` re-fetched (HTTP 200, 34 659 bytes, CRLF, sha256 `601779b2dcfe5f516b7554412c29b611dae12a4b45603a2cdfb18aa345b3743b`). Header `Genus\tNum_Genomes\tMean\tSTDEV \tThreshold`; 1 147 genus rows (`g__<Genus>`), 160 with a threshold (all genera with ≥ 10 genomes; the threshold column is empty below 10). Values: Escherichia 934 / 0.8486 / 0.0688 / 0.7110, Bordetella 0.2200, Mycobacterium 0.3959, Streptomyces 0.046, Azospirillum 0.319, min Synechococcus −0.222, max Borrelia 0.979. Thresholds equal mean − 2 SD within the 4-decimal rounding of the published columns (largest gap Azospirillum: 0.6174 − 2·0.1495 = 0.3184 vs 0.319). README: "for genera with >= 10 genomes) the SkewI threshold (2 standard deviations below mean)"; skewi.py default window "non-overlapping/adjacent windows of size 20kb". The repository `LICENSE` is **GNU GPL v3** (fetched, HTTP 200), so the table is not embedded in MIT-licensed Seqeron; callers load it with `ParseSkewIGenusThresholds`. A one-off parse of the full downloaded file gave 160 entries (Escherichia 0.711, Azospirillum 0.319, Streptomyces 0.046). The unit tests use 8 rows copied verbatim. Decision reference (skewi.py `-k 1000 --min-len 0`, re-run): BA1F extra 0.203158581311549 → below Escherichia/Mycobacterium/Bordetella, not below Streptomyces/Synechococcus; synthetic genome 1.0 → never below.

---

## Documented Corner Cases and Failure Modes

### From Rosalind BA1F

1. **Ties:** the problem asks for ALL positions minimizing the skew (e.g. `53 97`); an implementation returning a single position must define a deterministic tie-break (this unit returns the first/smallest minimizing index).
2. **Prefix indexing:** positions range over i ∈ [0, |Genome|], i.e. there are |Genome|+1 prefix values; Skew_0 = 0 before any base.

### From Grigoriev 1998 / Lobry 1996

3. **No asymmetry → no signal:** if a sequence has no net G/C strand bias the cumulative diagram is flat (amplitude 0) and the origin/terminus are not meaningfully resolved.

---

## Test Datasets

### Dataset: Rosalind BA1F sample

**Source:** Rosalind, Minimum Skew Problem (BA1F), https://rosalind.info/problems/ba1f/

| Parameter | Value |
|-----------|-------|
| Genome | `CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG` |
| Length | 100 |
| Minimum skew value | −4 |
| Positions of minimum skew | 53, 97 |
| First minimizing position (this unit's tie-break) | 53 |

### Dataset: Tiny worked examples (derived from the BA1F definition)

**Source:** Definition derivation from Rosalind BA1F / Grigoriev 1998.

| Sequence | Skew diagram (Skew_0..Skew_n) | Min value @ pos | Max value @ pos |
|----------|-------------------------------|-----------------|-----------------|
| `CCGGGG` | 0,−1,−2,−1,0,+1,+2 | −2 @ 2 | +2 @ 6 |
| `GGGCCC` | 0,+1,+2,+3,+2,+1,0 | 0 @ 0 | +3 @ 3 |
| `AATT`   | 0,0,0,0,0 | 0 @ 0 | 0 @ 0 |

### Dataset: circular-genome rotation (python brute force over Skew_0..Skew_{n−1})

| Sequence | Skew_n | circular minimizers | circular maximizers |
|----------|--------|---------------------|---------------------|
| 100 kb synthetic (`ACGTC`×6000+`AGGTC`×10000+`ACGTC`×4000) | 0 | [29997, 30000, 30001] | [79998, 79999] |
| same, rotated left 12345 | 0 | [17652, 17655, 17656] (−3531) → +12345 mod n = original | [67653, 67654] (+6470) |
| `CCGGG` | +1 | [2] | [0, 4] |
| `GGCCG` (= `CCGGG` rotated 3) | +1 | [0, 4] → mapped back {3, 2} ≠ {2} | [2] |
| `GGGCCC` / `G` | 0 / +1 | [0] / [0] | [3] / [0] |

### Dataset: windowed (Grigoriev) prediction — `numpy.cumsum(Bio.SeqUtils.GC_skew(seq, w)[:n//w])`, first argmin/argmax, centre = idx·w + w/2

| Sequence | w | origin (value) | terminus (value) |
|----------|---|----------------|------------------|
| 100 kb synthetic | 1000 | 29500 (−10.0) | 79500 (6.666666666666664) |
| BA1F extra | 1000 | 89500 (−0.35203934153706795) | 19500 (0.28817798587280524) |
| BA1F extra | 5000 | 87500 (−0.07055938354312184) | 17500 (0.056616604707736475) |
| BA1F sample | 10 | 45 (−0.3095238095238095) | 15 (0.5) |

---

## Assumptions

1. **ASSUMPTION: IsSignificant semantics.** No authoritative source defines a numeric "significance" cutoff for an origin call. The previous implementation used an invented threshold (`amplitude > count × 0.01`); that constant is untraceable and is removed. `IsSignificant` is redefined as the threshold-free, evidence-neutral predicate `max > min` (the diagram has non-zero amplitude, i.e. a detectable strand-composition asymmetry exists per Lobry 1996 / Grigoriev 1998). This is the weakest non-invented definition; callers needing a quantitative confidence measure should inspect the skew amplitude directly. (2026-10-09: the sourced quantitative measure is now available as `CalculateSkewIndex` — SkewIT SkewI; SkewIT gives per-genus thresholds only, so `IsSignificant` stays threshold-free. The per-genus rule is `IsSkewIBelowGenusThreshold` / `IsSkewIBelowThreshold`, with a caller-supplied table because SkewIT is GPL-3.0.)

---

## Recommendations for Test Coverage

1. **MUST Test:** BA1F sample genome returns PredictedOrigin = 53 (first minimizing prefix index). — Evidence: Rosalind BA1F sample output `53 97`.
2. **MUST Test:** Per-nucleotide skew uses G:+1, C:−1, A/T:0 with Skew_0 = 0 (small derived examples `CCGGGG`, `GGGCCC`). — Evidence: Rosalind BA1F definition.
3. **MUST Test:** PredictedTerminus = position of the global maximum. — Evidence: Grigoriev 1998 / Wikipedia (max = terminus).
4. **MUST Test:** First-occurrence tie-break when several positions share the extreme value. — Evidence: BA1F returns multiple minimizers.
5. **SHOULD Test:** Flat diagram (no G/C, or balanced) → origin = terminus = 0, IsSignificant = false. — Rationale: documented "no asymmetry" corner case.
6. **SHOULD Test:** Case-insensitive on the string overload; A/T bases do not move the diagram. — Rationale: counting convention.
7. **COULD Test:** Empty / null input returns zero prediction (string) / throws (DnaSequence). — Rationale: documented input handling.

---

## References

1. Rosalind. Minimum Skew Problem (BA1F). https://rosalind.info/problems/ba1f/
2. Grigoriev, A. (1998). Analyzing genomes with cumulative skew diagrams. Nucleic Acids Research 26(10):2286–2290. https://doi.org/10.1093/nar/26.10.2286
3. Lobry, J. R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. Molecular Biology and Evolution 13(5):660–665. https://pubmed.ncbi.nlm.nih.gov/8676740/
4. Wikipedia. GC skew. https://en.wikipedia.org/wiki/GC_skew (accessed 2026-06-14; used for cited primaries 2 and 3).

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-28** (review campaign B01): Grigoriev construction quoted (windowed (G−C)/(G+C) sum; window 1 = BA1F diagram). Reference cross-check added: python BA1F re-implementation + `numpy.cumsum(Bio.SeqUtils.GC_skew(seq, w))` on a 100 kb synthetic genome (`ACGTC`×6000 + `AGGTC`×10000 + `ACGTC`×4000; ori 30000 / ter 80000): per-base minimizers [29997, 30000, 30001, …] = −6000, maximizers [79998, 79999] = +4001; windowed (w=1000) min −10 in window ending 30000, max 20/3 in window ending 80000; rotation by 50000 (circular genome) → 79997/29998, which map back to 29997/79998. BA1F sample: Biopython `GC_skew(seq,1)` cumsum min −4 at [53, 97], max +2 at [16, 20, 21].
- **2026-10-09** (B01 finisher A1-2/A1-3, F16/F17): BA1F extra dataset (mirror) + all-minimizers API; circular option (positions mod n; rotation datasets); windowed Grigoriev prediction (Biopython cumsum datasets); SkewIT `skewi.py` opened and run (SkewI datasets, per-genus threshold source).
- **2026-10-09** (B01 finisher A2-1, F29): MCP `predict_replication_origin` `windowSize` / `skewIndexWindow`. Reference runs: Biopython 1.88 `numpy.cumsum(GC_skew(ba1f, w)[:n//w])` first argmin/argmax → centre: w10 45 (−0.3095238095238095) / 15 (0.5), w25 37 (−0.15384615384615385) / 62 (0.11888111888111885). SkewIT `skewi.py -k 4 --min-len 0` (headers with "complete", `%r`): BA1F sample → 0.16 (plus the k4 values above; 12-window `GGGGCCCC`×6 → no output).
- **2026-10-09** (B01 finisher A2-3, F31): SkewIT threshold table + LICENSE (GPL-3.0) fetched; per-genus decision API with a caller-supplied table; full-file parse 160/1 147; skewi.py `-k 1000` decision datasets.
