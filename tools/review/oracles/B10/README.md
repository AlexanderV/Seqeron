# B10 (Assembly) reference oracles

Python oracles used to cross-check B10 fixes (see docs/Validation/review-2026-09/B10.md and the unit Evidence files).
Run with `python3 -I <file>`.

| File | Unit | Reference transcribed |
|---|---|---|
| `scaffold_oracle.py` | ASSEMBLY-SCAFFOLD-001 | Huson, Reinert & Myers 2002 Alg. 5.1 greedy scaffolding; ABySS overlap collapse; AGP v2.1 gaps |
| `cover_oracle.py` | ASSEMBLY-COVER-001 | Bowtie 2 end-to-end ungapped placement → SAM → pysam `samtools depth -a` (`-Q 2`) |
| `emboss_merger_port.py` | ASSEMBLY-MERGE-001 | EMBOSS `merger.c` overlap conflict rule |
| `gaps_oracle.py` | ASSEMBLY-GAPS-001 | QUAST `qutils.split_by_ns` (verbatim) + `--min-contig` filter; seqtk `stk_gap`; numpy median |
| `musket_oracle.py` | ASSEMBLY-CORRECT-001 | Musket 1.1 `ParaCorrect` stages 1–3 (defects (a)/(b) corrected, see Evidence) |
| `trim_cutadapt_oracle.py` | ASSEMBLY-TRIM-001 | cutadapt 5.2 `quality_trim_index`: 3'-only default (-q C = (0,C)) vs both ends, 20000 random reads vs the running-sum formulation |
| `correct_default_oracle.py` | ASSEMBLY-CORRECT-001 | `musket_oracle.py` two-sided stage, canonical vs forward-only spectrum (default flip A9-1) |
| `completeness_oracle.py` | ASSEMBLY-COMPLETE-001 | Merqury `kmerHistToPloidyDepth` (filt), `spectra-cn.sh` completeness, `qv.sh` QV; GenomeScope 1.0 `uniroot` error rate + `total_len`; BUSCO `hmmer.py` percentages/notation; compleasm v0.2.6 pass filter + `record_1st_gene_label`; Biopython 1.88 local alignment |

The DBG (Rosalind BA3M) and consensus (Biopython 1.85 `dumb_consensus`) oracles are embedded in their Evidence files;
the STATS oracle is QUAST `quast_libs/N50.py`, imported unchanged; the TRIM oracle is pip `cutadapt==5.2`
`cutadapt.qualtrim.quality_trim_index`.
