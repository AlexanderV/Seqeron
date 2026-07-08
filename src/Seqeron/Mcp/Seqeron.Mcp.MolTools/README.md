# Seqeron.Mcp.MolTools

MCP server — **Primer/probe/CRISPR design, codon optimization, restriction analysis, thermodynamics.**

Exposes **47 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/moltools/`](../../../../docs/mcp/tools/moltools). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.MolTools
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.MolTools"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (47)

| Tool | Description |
|------|-------------|
| `analyze_oligo` | Returns Tm, GC fraction, molecular weight (Da), and 260 nm extinction coefficient (M⁻¹·cm⁻¹) for a short oligonucleotide. |
| `blunt_cutters` | Lists all built-in restriction enzymes that produce blunt ends (both strands cut at the same position). |
| `build_codon_table` | Derives a per-organism CodonUsageTable from a reference coding sequence by computing per-amino-acid relative codon frequencies (RNA alpha… |
| `cai_from_organism_table` | Computes the Codon Adaptation Index (Sharp & Li 1987) for a coding sequence against an organism codon-usage FREQUENCY table (distinct fro… |
| `codon_adaptation_index` | Codon Adaptation Index (Sharp & Li 1987) using a caller-supplied reference RSCU table (typically derived from highly expressed genes), co… |
| `codon_usage_statistics` | Aggregate codon-usage report for a coding sequence: per-codon counts, RSCU, Effective Number of Codons (ENC), total codons, GC% at codon… |
| `compare_codon_usage` | Codon-frequency similarity between two coding sequences: 1 − ½·Σ|f1−f2| ∈ [0,1] (1 = identical codon distribution, 0 = disjoint). |
| `compatible_enzymes` | Enumerates all pairs of built-in restriction enzymes whose ends can be ligated to each other — either both produce blunt ends, or both pr… |
| `count_codons` | Counts occurrences of each ACGT codon in a coding DNA sequence (frame 0, non-overlapping triplets). |
| `crispr_specificity_score` | Aggregates off-target hits (≤4 mismatches) for a guide RNA against a genome into a single specificity score in 0..100 (100 = no off-targets; |
| `crispr_system_info` | Returns the metadata record (name, PAM sequence, guide length, PAM placement relative to target, description) for a known CRISPR nuclease… |
| `design_antisense_probes` | Reverse-complements the supplied mRNA-sense sequence and runs the probe designer on it; |
| `design_guide_rnas` | Generates and scores guide-RNA candidates whose Cas9/Cas12a cut site falls inside the requested region. |
| `design_molecular_beacon` | Designs a hairpin molecular-beacon probe: GC-rich complementary stems (stem5 = ⌊stem_length/2⌋ Gs + remaining Cs, stem3 = its reverse com… |
| `design_primers` | Designs forward and reverse PCR primers flanking a target region in a DNA template; |
| `design_probes` | Designs hybridization probes by scanning the target for length-window candidates and ranking by GC%, Tm, homopolymers, self-complementari… |
| `design_tiling_probes` | Generates fixed-length probes covering the entire target with a configurable overlap (step = probe_length − overlap). |
| `digest_summary` | Aggregate statistics over a simulated linear restriction digest: total fragment count, fragment sizes (descending), largest/smallest frag… |
| `effective_number_of_codons` | Effective Number of Codons (Wright's Nc), measuring how far a gene departs from uniform synonymous-codon usage. |
| `enzymes_by_cut_length` | Lists all built-in restriction enzymes whose recognition sequence has exactly the specified length in base pairs (e.g. |
| `enzymes_compatible` | Returns whether two named built-in restriction enzymes produce ligatable (compatible) ends — both blunt, or the same overhang type and se… |
| `evaluate_guide_rna` | Scores a single guide RNA against on-target quality heuristics: overall GC%, seed-region GC%, polyT (Pol III terminator) presence, self-c… |
| `evaluate_primer` | Evaluates a single primer sequence against quality criteria and returns a scored candidate: length, GC%, Tm, longest homopolymer, hairpin… |
| `find_all_restriction_sites` | Finds sites for EVERY built-in restriction enzyme on both strands of a DNA sequence. |
| `find_off_targets` | Naïve genome scan: enumerates all PAM sites in the genome and reports those whose target differs from the guide by 1..max_mismatches (the… |
| `find_pam_sites` | Finds all PAM matches (forward + reverse strand) for the chosen CRISPR system. |
| `find_rare_codons` | Reports every codon in a coding sequence whose frequency in the target organism's codon-usage table is below the threshold (default 0.15)… |
| `find_restriction_sites` | Finds restriction sites for one or more named built-in enzymes on both strands of a DNA sequence. |
| `generate_primer_candidates` | Enumerates all primer candidates of admissible lengths (parameters.MinLength..MaxLength) at every start position within a region of the t… |
| `get_enzyme` | Looks up a built-in restriction enzyme by name (case-insensitive) and returns its recognition sequence, cut positions and organism. |
| `hairpin_potential` | Detects whether a sequence can fold into a hairpin: a self-complementary stem of at least min_stem_length separated by a loop of at least… |
| `longest_dinucleotide_repeat` | Returns the number of repeat units in the longest dinucleotide tandem repeat (e.g. |
| `longest_homopolymer` | Returns the length of the longest run of identical consecutive nucleotides (e.g. |
| `oligo_concentration_from_absorbance` | Computes oligonucleotide concentration in µM from the Beer–Lambert law: c = A₂₆₀ / (ε · path) · 1e6. |
| `oligo_extinction_coefficient` | Sums per-base 260 nm molar extinction contributions (A=15400, C=7400, G=11500, T=8700, U=9900 M⁻¹·cm⁻¹; |
| `optimize_codons` | Optimizes a coding sequence for expression in a target organism using one of five strategies (MaximizeCAI, BalancedOptimization (default)… |
| `primer_dimer` | Heuristic 3'-end primer-dimer check between two primers: reverse-complements primer2 and counts complementary positions in an up-to-8-bp… |
| `primer_melting_temperature` | Computes a primer's melting temperature (Tm, °C): Wallace rule Tm = 2·(A+T) + 4·(G+C) for < 14 valid bases, or Marmur–Doty Tm = 64.9 + 41… |
| `primer_melting_temperature_salt` | Primer Tm with a Schildkraut–Lifson salt correction: adds 16.6·log10([Na+]/1000) to the Wallace/Marmur–Doty Tm, rounded to one decimal. |
| `reduce_secondary_structure` | Greedy synonymous-codon swap that lowers a heuristic local self-complementarity score within a sliding window, reducing mRNA secondary st… |
| `remove_restriction_sites` | Synonymously rewrites codons to eliminate the listed restriction recognition sequences from a coding sequence while preserving the encode… |
| `restriction_digest` | Simulates a restriction digest of a linear DNA molecule with one or more named enzymes and yields the resulting fragments in 5'→3' order… |
| `restriction_map` | Builds a restriction map of a DNA sequence: all forward+reverse sites, sites grouped by enzyme, the total forward-strand site count, uniq… |
| `rscu` | Relative Synonymous Codon Usage (RSCU) per codon: observed count / count-expected-if-uniform among its synonymous codons. |
| `sticky_cutters` | Lists all built-in restriction enzymes that produce sticky (cohesive) ends — a staggered cut leaving a 5' or 3' single-stranded overhang. |
| `three_prime_stability` | SantaLucia (1998) nearest-neighbor ΔG°37 (kcal/mol) of a primer's last 5 bases, including initiation terms (1 M NaCl), matching Primer3 P… |
| `validate_probe` | Validates a probe against a set of reference sequences using ungapped k-mismatch (Hamming) approximate matching. |
