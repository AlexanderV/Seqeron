# B10 (Assembly) reference oracles

Python oracles used to cross-check B10 fixes (see docs/Validation/review-2026-09/B10.md and the unit Evidence files).
Run with `python3 -I <file>`.

| File | Unit | Reference transcribed |
|---|---|---|
| `scaffold_oracle.py` | ASSEMBLY-SCAFFOLD-001 | Huson, Reinert & Myers 2002 Alg. 5.1 greedy scaffolding; ABySS overlap collapse; AGP v2.1 gaps |
| `cover_oracle.py` | ASSEMBLY-COVER-001 | Bowtie 2 end-to-end ungapped placement → SAM → pysam `samtools depth -a` (`-Q 2`) |
| `emboss_merger_port.py` | ASSEMBLY-MERGE-001 | EMBOSS `merger.c` overlap conflict rule |
| `musket_oracle.py` | ASSEMBLY-CORRECT-001 | Musket 1.1 `ParaCorrect` stages 1–3 (defects (a)/(b) corrected, see Evidence) |

The DBG (Rosalind BA3M) and consensus (Biopython 1.85 `dumb_consensus`) oracles are embedded in their Evidence files;
the STATS oracle is QUAST `quast_libs/N50.py`, imported unchanged; the TRIM oracle is pip `cutadapt==5.2`
`cutadapt.qualtrim.quality_trim_index`.
