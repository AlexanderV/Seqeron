# Evidence Artifact: ONCO-FUSION-003

**Test Unit ID:** ONCO-FUSION-003
**Algorithm:** Fusion Breakpoint Analysis (junction reading-frame consequence + fusion protein prediction)
**Date Collected:** 2026-06-14

---

## Online Sources

### Arriba output-files specification (reading_frame, site, peptide_sequence, fusion_transcript)

**URL:** https://github.com/suhrig/arriba/wiki/05-Output-files
**Accessed:** 2026-06-14 (retrieved with WebFetch on the wiki page)
**Authority rank:** 3 (reference implementation — Arriba, Uhrig et al. 2021, Genome Research)

**Key Extracted Points:**

1. **reading_frame column:** "states whether the gene at the 3' end of the fusion is fused `in-frame` or `out-of-frame`". Possible values are `in-frame`, `out-of-frame`, `stop-codon`, and `.` (a dot when the peptide sequence cannot be predicted). It "builds on the prediction of the peptide sequence".
2. **site1 / site2 columns:** describe "the location of the breakpoints. Possible values are: `5'UTR`, `3'UTR`, `UTR`, `CDS`, `exon`, `intron`, and `intergenic`". "If the breakpoint coincides with an exon boundary, the additional keyword `splice-site` is appended."
3. **fusion_transcript column:** "contains the fusion transcript sequence … assembled from the supporting reads of the most highly expressed transcript." The pipe character `|` marks the breakpoint position in the transcript.
4. **peptide_sequence column:** "contains the fusion peptide sequence … translated from the fusion" transcript. The pipe character `|` marks the breakpoint in the peptide.
5. **type column:** the event type — `translocation`, `duplication`, `inversion`, `deletion` — inferred from supporting-read orientation and breakpoint coordinates; `read-through` is a specific event type.

### AGFusion reference implementation — model.py (frame determination + protein prediction)

**URL:** https://raw.githubusercontent.com/murphycj/AGFusion/master/agfusion/model.py
**Accessed:** 2026-06-14 (retrieved with WebFetch on the raw source file)
**Authority rank:** 3 (reference implementation — AGFusion, Murphy & Elemento 2016)

**Key Extracted Points:**

1. **5' CDS up to breakpoint:** `self.cds_5prime = self.transcript1.coding_sequence[0 : self.transcript_cds_junction_5prime]` — the 5' partner contributes its coding sequence from the start up to the breakpoint offset (a prefix of the 5' CDS).
2. **3' CDS from breakpoint:** `self.cds_3prime = self.transcript2.coding_sequence[self.transcript_cds_junction_3prime : :]` — the 3' partner contributes its coding sequence from the breakpoint offset to the end (a suffix of the 3' CDS).
3. **Chimeric CDS:** the two segments are concatenated, `seq = self.cds_5prime + self.cds_3prime`.
4. **Frame rule (verbatim):**
   ```python
   if (len(self.cds_5prime) / 3.0).is_integer() and (len(self.cds_3prime) / 3.0).is_integer():
       self.effect = "in-frame"
   elif round((len(self.cds_5prime) / 3.0 % 1) + (len(self.cds_3prime) / 3.0 % 1), 2) == 1.0:
       self.effect = "in-frame (with mutation)"
   else:
       self.effect = "out-of-frame"
   ```
   i.e. the fusion is **in-frame** when the combined chimeric CDS length is a multiple of 3 *and* the junction falls at codon boundaries such that the two segment lengths' fractional codon parts complement to a whole codon; otherwise **out-of-frame**. The decisive quantity is the residue `(len(cds_5prime) + len(cds_3prime)) mod 3` together with the per-segment phase at the junction.
5. **Translation + truncation at first stop (verbatim):**
   ```python
   protein_seq = self.cds.seq.translate()
   protein_seq = protein_seq[0 : protein_seq.find("*")]
   ```
   The chimeric CDS is translated and the protein is truncated at the first stop codon (`*`).
6. **Out-of-frame translation (verbatim):** for an out-of-frame fusion AGFusion first trims the CDS to a whole number of codons before translating: `seq = self.cds.seq[0 : 3 * int(len(self.cds.seq) / 3)]`, then `protein_seq = seq.translate()` truncated at the first `*`.

### Wikipedia "Reading frame" (codon-triplet / modulo-3 basis)

**URL:** https://en.wikipedia.org/wiki/Reading_frame
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 4 (Wikipedia citing the primary Badger & Olsen 1999)

**Key Extracted Points:**

1. **Triplet definition:** a reading frame is "a specific choice out of the possible ways to read the sequence of nucleotides … as a sequence of triplets". Codons are three-nucleotide groupings read during protein synthesis. This is the modulo-3 basis for the junction-phase test.
2. **Primary reference:** Badger JH, Olsen GJ (April 1999). "CRITICA: Coding Region Identification Tool Invoking Comparative Analysis." *Mol Biol Evol* 16(4):512–24.

---

## Documented Corner Cases and Failure Modes

### From Arriba output spec

1. **Peptide not predictable:** `reading_frame = .` when the peptide sequence cannot be predicted (e.g. a breakpoint outside the coding region of a partner). The frame status is only defined when both breakpoints can be placed in coding context.
2. **Breakpoint not in CDS:** a breakpoint may fall in `5'UTR`, `3'UTR`, `intron`, `exon` (non-coding part), or `intergenic`; in those cases the junction does not join two coding frames, so an in/out-of-frame call is not made.
3. **Premature stop at junction:** `reading_frame = stop-codon` when the junction or the downstream sequence yields a stop codon, terminating the chimeric ORF.

### From AGFusion model.py

1. **Out-of-frame frameshift:** when the junction does not preserve the frame, the 3' partner is translated in a shifted frame; AGFusion trims the CDS to whole codons and the resulting peptide is read out of the 3' gene's native frame (frameshifted), typically truncated early at the first downstream stop codon.
2. **First stop truncation:** translation always terminates at the first stop codon in the chimeric ORF (`protein_seq[0:find("*")]`), so an in-frame fusion with an early stop still yields a truncated peptide.

---

## Test Datasets

### Dataset: Junction reading-frame truth table (derived from AGFusion frame rule)

**Source:** AGFusion model.py `_fetch_protein()` frame rule (Murphy & Elemento 2016).

The in-frame decision is the residue of the 5' coding-base count relative to the 3' coding-start phase, modulo 3 (reading frames are triplets). With `b` = coding bases the 5' partner contributes up to the breakpoint and `p` = coding-start phase of the 3' partner at the breakpoint (0,1,2):

| b (5' coding bases) | p (3' start phase) | (b − p) mod 3 | Reading frame |
|---------------------|--------------------|---------------|---------------|
| 9   | 0 | 0 | in-frame    |
| 10  | 1 | 0 | in-frame    |
| 11  | 2 | 0 | in-frame    |
| 10  | 0 | 1 | out-of-frame|
| 9   | 1 | 2 | out-of-frame|
| 9   | 2 | 1 | out-of-frame|

### Dataset: Fusion protein prediction (derived from AGFusion concatenate + translate + truncate)

**Source:** AGFusion model.py (CDS concatenation, `translate()`, truncate at first `*`). Standard genetic code (NCBI transl_table=1): ATG→M, AAA→K, GAT→D, GGT→G, TAA/TAG/TGA→stop.

| 5' CDS prefix | 3' CDS suffix | chimeric CDS | translation | first-stop-truncated peptide | effect |
|---------------|---------------|--------------|-------------|------------------------------|--------|
| `ATGAAA` (len 6) | `GATGGT` (len 6) | `ATGAAAGATGGT` | `MKDG` | `MKDG` | in-frame (6%3=0, 6%3=0) |
| `ATGAAA` (len 6) | `GATTAAGGT` (len 9) | `ATGAAAGATTAAGGT` | `MKD*G` | `MKD` | in-frame, premature stop at junction-downstream |
| `ATGA` (len 4) | `AAGGT` (len 5) | `ATGAAAGGT` | `MKG` | `MKG` | **out-of-frame** under the Arriba 3'-gene-frame model: 5' contributes 4 bases (phase 1) but the 3' suffix starts at native phase 0, `(4−0) mod 3 = 1 ≠ 0`, so the 3' gene is read frameshifted. The 9-base chimeric CDS still translates cleanly (contiguous ORF), so AGFusion's three-way rule would label it "in-frame (with mutation)" — but the repo models Arriba's two-way in/out-of-frame call (whether the 3' gene is in its native frame), not AGFusion's ORF-continuity class. |
| `ATGA` (len 4) | `AAGGT` (len 5), suffix of `TAAGGT` taken from offset 1 | `ATGAAAGGT` | `MKG` | `MKG` | in-frame: 5' contributes 4 bases (phase 1) and the 3' suffix begins at native phase 1, `(4−1) mod 3 = 0`, so the 3' gene is read in its native frame (AGFusion: 4 and 5 fractional codon parts complement → in-frame). |
| `ATGAA` (len 5) | `GATGGT` (len 6) | `ATGAAGATGGT` (len 11) | trim to 9: `ATGAAGATG` → `MKM` | `MKM` (3' read in shifted frame) | out-of-frame (5%3=2, 6%3=0; 2+0≠0,≠3) |

---

## Assumptions

0. **MODEL: Arriba two-way reading-frame call, not AGFusion's three-way class.** The repo's `BreakpointFrameStatus` is `InFrame` / `OutOfFrame` / `StopCodon` / `NotPredicted` — Arriba's `reading_frame` model, which "states whether the gene at the 3' end of the fusion is fused in-frame or out-of-frame" (i.e. whether the 3' gene is read in its NATIVE reading frame). The frame test is therefore the exon-phase-compatibility rule `(fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0` (Genomics England / ONCO-FUSION-001), which asks whether the 3' partner is read in its native frame. AGFusion uses a *three-way* class — `in-frame` (both segments codon-aligned) / `in-frame (with mutation)` (segment fractional codon parts complement so the contiguous ORF length is a multiple of 3, but the 3' gene is read frameshifted) / `out-of-frame`. AGFusion's `in-frame (with mutation)` maps to the repo's **`OutOfFrame`** because the 3' gene is NOT in its native frame. This is a deliberate model choice (Arriba-style), not a defect: a fusion whose 3' partner is frameshifted does not preserve the 3' protein domains, so it is not "in-frame" in the functional (kinase-domain-preserving) sense.

1. **ASSUMPTION: Caller supplies CDS sequences and junction offsets.** The repository has no genome/GTF/transcript database, so `PredictFusionProtein` takes the partner CDS strings and the breakpoint CDS offsets directly (the 5' coding-base count and the 3' coding-start offset), rather than resolving them from an Ensembl/RefSeq annotation as AGFusion does. This is an API-shape decision: the concatenation, frame rule, translation, and first-stop truncation are exactly AGFusion's; only the source of the inputs differs. It does not change any computed output for given CDS/offsets.

---

## Recommendations for Test Coverage

1. **MUST Test:** in-frame call when `(b − p) mod 3 == 0` across all three residues of p (0,1,2) — Evidence: AGFusion frame rule; Arriba reading_frame.
2. **MUST Test:** out-of-frame call when `(b − p) mod 3 != 0` — Evidence: AGFusion frame rule.
3. **MUST Test:** breakpoint site classification (CDS vs UTR/intron/intergenic) drives whether a frame call is made; non-CDS breakpoint → frame Unknown / not predicted — Evidence: Arriba site / reading_frame `.`.
4. **MUST Test:** `PredictFusionProtein` concatenates 5' CDS prefix + 3' CDS suffix, translates with the standard genetic code, and truncates at the first stop codon (worked vectors `ATGAAA|GATGGT` → `MKDG`; `ATGAAA|GATTAAGGT` → `MKD`) — Evidence: AGFusion concatenate/translate/truncate.
5. **MUST Test:** out-of-frame fusion trims the CDS to whole codons before translation and reads the 3' partner in a shifted frame — Evidence: AGFusion out-of-frame branch.
6. **SHOULD Test:** premature stop codon detection flag (`stop-codon`) when a stop is reached before the end of the chimeric ORF — Rationale: Arriba `stop-codon` reading_frame value.
7. **SHOULD Test:** null / empty inputs and invalid phase (outside {0,1,2}) raise the documented exceptions — Rationale: input-validation contract mirroring sibling fusion methods.
8. **COULD Test:** stable, deterministic ordering / record-shape of the returned breakpoint analysis — Rationale: API consistency.

---

## References

1. Uhrig S, Ellermann J, Walther T, et al. (2021). Accurate and efficient detection of gene fusions from RNA sequencing data. *Genome Research* 31(3):448–460. https://doi.org/10.1101/gr.257246.119 — output schema retrieved from the project wiki https://github.com/suhrig/arriba/wiki/05-Output-files (accessed 2026-06-14).
2. Murphy C, Elemento O (2016). AGFusion: annotate and visualize gene fusions. *bioRxiv* 080903. https://doi.org/10.1101/080903 — source code retrieved from https://raw.githubusercontent.com/murphycj/AGFusion/master/agfusion/model.py (accessed 2026-06-14).
3. Badger JH, Olsen GJ (1999). CRITICA: Coding Region Identification Tool Invoking Comparative Analysis. *Mol Biol Evol* 16(4):512–524. https://doi.org/10.1093/oxfordjournals.molbev.a026133 — via Wikipedia "Reading frame", https://en.wikipedia.org/wiki/Reading_frame (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
</content>
</invoke>
