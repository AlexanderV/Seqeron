# Evidence: RNA-ACCESS-001 — Region accessibility (McCaskill unpaired / accessibility probabilities)

**Status:** Stage A / Stage B re-review 2026-11-25 (campaign B12). Engine: RNA-PARTITION-001 ([Evidence](./RNA-PARTITION-001-Evidence.md)).

## 1. Sources (what was actually opened)

| # | Source | Used for |
|---|--------|----------|
| 1 | McCaskill JS (1990) *Biopolymers* 29:1105-1119 | Z = Σ exp(−E/RT); P(i,j); per-base unpaired probability |
| 2 | Bernhart SH, Hofacker IL, Stadler PF (2006) *Bioinformatics* 22:614-615, "Local RNA base pairing probabilities in large sequences"; ViennaRNA RNAplfold man page | accessibility = probability that the length-u stretch ending at i is entirely unpaired (`-u`), computed in windows W with maximum pair span L |
| 3 | Kertesz M, Iovino N, Unnerstall U, Gaul U, Segal E (2007) *Nat Genet* 39:1278-1284 | accessibility as the free-energy cost ΔΔG_open of unpairing the target site: −RT ln P |
| 4 | Agarwal V et al. (2015) *eLife* 4:e05005 | TargetScan SA = log10 unpaired probability of a 14-nt window (RNAplfold `-L 40 -W 80 -u 20`) |
| 5 | ViennaRNA 2.7.2 sdist `partfunc/pf_window.c` (opened lines 682-960: `vrna_probs_window`; `compute_pU`), Python `RNA` 2.7.2 | oracle and the semantics of the local variant |

## 2. Definition implemented

`CalculateRegionUnpairedProbability(seq, windowEnd, windowLength, minLoopSize, temperature)` = `Z_open / Z` of the WHOLE sequence
(pseudoknot-free, Turner 2004, d2): `Z_open` is the same inside recursion with every pair incident to a region position removed.
`CalculateRegionOpeningEnergy` = `−RT ln P` (RT = 1.98717e-3·T). Per-base `p_unpaired(i) = 1 − Σ_j P(i,j)` is its u = 1 case.

Equivalences (all numerically confirmed below): `Z_open/Z = exp(−(G_constrained − G_ensemble)/RT)` with the region hard-constrained
unpaired in ViennaRNA (`hc_add_up(p, CONSTRAINT_CONTEXT_ALL_LOOPS)`); `RNA.pfl_fold_up(seq, u, W = n, L = n)[end+1][u]` (RNAplfold with window and span ≥ n).

## 3. Stage B results (ViennaRNA 2.7.2, NNDB-precision parameter file, `md` explicit, dangles = 2)

| Check | Queries | Result |
|-------|---------|--------|
| exact constrained ratio, 37 °C, 8-100 nt, u ∈ {1, 2, 4, 7, 14, random ≤ n}, regions anywhere | 880 (220 seq) | max |Δ| 0.0056, mean 1.9e-5 |
| same, 20 °C / 60 °C, 20-100 nt | 240 + 240 | max 0.00045 / 0.0036 |
| same, 120-200 nt | 80 | max 6.3e-6 |
| `RNA.pfl_fold_up(seq, 14, n, n)` vs ours | 600 (15-90 nt, u ∈ {2, 4, 7, 14}) | max |Δ| 1.0e-3 |
| constrained ratio vs `pfl_fold_up(W = L = n)` | 42 locked cases | max 2e-6 (same quantity) |
| exhaustive enumeration with `CalculateFreeEnergy` (u 1…6, 8-15 nt, minLoop 3 and 5, 20/37/70 °C) | tests | 1e-9 |

The residual (≤ 0.006) is ViennaRNA's `SMOOTH()` of exterior dangle weights and rounded parameter tables (RNA-PARTITION-001 §9); sequences with
repeated neighbours carry NNDB-only terms ViennaRNA lacks and are validated against enumeration only. Test-tolerance: max(0.0015, 2 %).

Edge cases checked: u = 1 (= per-base unpaired, 1e-9), u = n (= 1/Z), region at the 5′ / 3′ end, region longer than n, `windowEnd` negative or ≥ n,
`windowLength` ≤ 0, null / empty sequence (all `ArgumentOutOfRangeException`), DNA `T` and lower case (identical result), `N` (never pairs ⇒ 1.0),
`minLoopSize` < 3 clamped, larger `minLoopSize` ⇒ more accessible, temperature ≤ 0 throws.

## 4. The joint probability is not a product

`GGGGAAAACCCCUUUUGGGG`, 5′ arm (positions 0-3): joint P > 1.2 × Π p_unpaired(i) (positions of one helix are strongly correlated); for every window
`1 − Σ(1 − p_i) ≤ P ≤ min p_i`.

## 5. Local RNAplfold semantics (W < n, L < n): not implemented

`pf_window.c` computes, per window, exterior contributions from the full sliding matrices (`compute_probs`, `compute_pU` with `MAXLOOP` margins) and averages them. Test of
the obvious reading (mean over windows `[a, a+W−1]` containing the region of `Z_open/Z` of an independent fold with `max_bp_span = L`): 60 nt, W = 30, L = 15, u = 4, end 57:
plfold 0.1364 vs reading 0.2016; end 30: 0.3215 vs 0.3287; end 10: 0.2862 vs 0.2934. The mean of truncated edge windows is worse (end 3, u = 2: 0.187 vs 0.373). The semantics therefore
require a port of the sliding-window recursion (~1500 lines), not a loop over the existing engine. Status: documented, DOABLE-LARGE, registered with the lead.
The TargetScan SA path keeps its explicit `PlfoldUnpairedProbability` override (exact RNAplfold number) and otherwise estimates with the single-window value.

## 6. Reproduction

Oracle scripts (scratchpad, not in the repository): constrained ratio `fc2.hc_add_up(p, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)` for `p` in the region on a fresh `fold_compound(seq, md)`; `G = fc.pf()[1]`;
`ratio = exp(−(G_c − G)/(0.00198717·T))`; `RNA.pfl_fold_up(seq, u, n, n)[end+1][u]`.
