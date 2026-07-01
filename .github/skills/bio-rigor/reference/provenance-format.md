# Provenance block — format & worked example

Every computed biological result MUST end with a **Provenance** block. It makes the result
reproducible and auditable: anyone can re-run the same tools / `Method ID`s with the same
parameters and get the same answer.

## Format

A fenced block titled `Provenance`, listing every tool (MCP) or `Method ID` (C# API) **in call
order**, each with its parameters and the value it produced. Then any cross-checks and caveats.

```
Provenance
1) <tool_or_MethodID>(param=value, …) → <key output>
2) <tool_or_MethodID>(param=value, …) → <key output>
…
Cross-check: <second independent tool/invariant> → agrees / disagrees
Envelope: <guarded units touched + LimitationPolicy mode, or "none guarded">
Caveat: <alpha / not-for-clinical-use if decision-relevant; else omit>
```

Rules:
- **In order.** The list is the actual call sequence, not a post-hoc summary.
- **Parameters explicit.** Record the exact inputs (pattern, k, matrix, mode, thresholds).
- **No manual steps.** If a value was not produced by a tool/`Method ID`, it does not belong here —
  and per rule 1 it should not have been computed at all.
- **Units & coordinate base** stated on outputs where they matter (°C, %, bp, 0-based).

## Worked example — cloning insert QC (from `docs/mcp/README.md`)

**Task:** for an insert `>seq1 GCGCGAATTCATGGATCCATAT`, report GC% and EcoRI/BamHI sites (0-based).

Result:

```
| id   | length | gc_percent | EcoRI_sites | BamHI_sites |
|------|-------:|-----------:|-------------|-------------|
| seq1 |     22 |      45.45 | [4]         | [12]        |
```

```
Provenance
1) fasta_parse(input=">seq1 …") → sequence = "GCGCGAATTCATGGATCCATAT", length = 22 bp
2) gc_content(sequence) → gc_percent = 45.45 %, gcCount = 10, totalCount = 22
3) suffix_tree_find_all(pattern="GAATTC") → positions = [4]   (0-based)
4) suffix_tree_find_all(pattern="GGATCC") → positions = [12]  (0-based)
Cross-check: gc_percent 45.45 % consistent with gcCount/totalCount = 10/22; site at 4 confirmed
             by the "GAATTC" substring starting index 4.
Envelope: none guarded (Sequence + Core tools, all within contract).
Caveat: cloning QC only; alpha software — independently verify sites before ordering constructs.
```

## Worked example — PCR primer QC (C# API mode)

Same task via `Seqeron.Genomics`:

```
Provenance
1) FastaParser.Parse(text) → FWD (20 bp), REV (20 bp)
2) DnaSequence.Validate(FWD) → valid; DnaSequence.Validate(REV) → valid
3) SequenceStatistics.GcContent(FWD) → 50.00 %
4) SequenceStatistics.GcContent(REV) → 60.00 %
5) MeltingTemperature.Calculate(FWD) → 51.8 °C
6) MeltingTemperature.Calculate(REV) → 55.9 °C
Cross-check: tm_diff = |51.8 − 55.9| = 4.1 °C derived from the two tool outputs (not re-measured).
Envelope: none guarded.
Caveat: alpha software; validate primers empirically before running a real PCR.
```

The `Method ID`s to cite come from each unit's doc under `docs/algorithms/<Area>/<Unit>.md` and the
per-tool docs under `docs/mcp/tools/<server>/<tool>.md`.
