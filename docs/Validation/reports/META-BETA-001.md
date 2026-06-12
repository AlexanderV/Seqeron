# Validation Report: META-BETA-001 — Beta Diversity (Bray–Curtis, Jaccard)

- **Validated:** 2026-06-12   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.CalculateBetaDiversity(string, IReadOnlyDictionary<string,double>, string, IReadOnlyDictionary<string,double>)`
  (helpers `CalculateBrayCurtis`, `CalculateJaccardDistance`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Scope

The implemented unit reports two beta-diversity metrics plus species-overlap counts:
- **Bray–Curtis dissimilarity** (abundance-based)
- **Jaccard distance** (presence/absence; reported as distance = 1 − J, correctly labelled)

`UniFracDistance` is a stub returning 0 (documented: requires a phylogenetic tree, out of scope).
Whittaker beta and Sørensen are **not** implemented; the checklist cites Whittaker (1960) only for
the beta-diversity *concept*, not as an implemented metric. No gap versus the canonical method.

## Stage A — Description

### Sources opened
- Wikipedia, **Bray–Curtis dissimilarity** — confirms both formula forms, range, interpretation, and worked example.
- Wikipedia, **Jaccard index** — confirms similarity/distance formulas, binary form, range, domain (non-empty sets).
- Cited primary literature: Bray & Curtis (1957), Jaccard (1901), Whittaker (1960) — concept/formula origins (via Wikipedia's cited refs).

### Formula check (each verified exactly)
- **Bray–Curtis:** `BC = 1 − 2·Σ min(xᵢ,yᵢ) / Σ(xᵢ+yᵢ) = Σ|xᵢ−yᵢ| / Σ(xᵢ+yᵢ)`.
  Wikipedia states `BC_jk = 1 − 2C_jk/(S_j+S_k)` with `C_jk = Σ min(N_ij,N_ik)`. Range [0,1];
  0 = same composition, 1 = no shared species. **Matches** spec/Evidence exactly.
- **Jaccard:** similarity `J = |A∩B|/|A∪B| = |A∩B|/(|A|+|B|−|A∩B|)`; binary form `J = M11/(M01+M10+M11)`;
  **Jaccard distance** `d_J = 1 − J`. Range [0,1]; J=0 no shared elements, J=1 identical.
  Spec reports the **distance** (1 − J) and labels it "Jaccard Distance" — correct, no sim/distance mislabel.
- **Sørensen / Whittaker:** not implemented — not applicable.

### Edge-case semantics (sourced)
- Identical → BC 0, Jaccard distance 0; disjoint → BC 1, Jaccard distance 1.
- Both empty → returns 0 (documented design decision: identical inputs ⇒ zero distance; denominator-0 guard,
  consistent with scipy convention). One empty → 1.
- Zero-abundance species treated as absent for presence/absence (Jaccard domain) per spec.

### Independent cross-check (hand computation)
- **Wikipedia aquarium example:** Tank1{Goldfish6,Guppy7,Rainbow4}, Tank2{Goldfish10,Guppy0,Rainbow6}:
  C=6+0+4=10, S_j=17, S_k=16 → BC = 1 − 20/33 = **13/33 ≈ 0.3939** (Wikipedia rounds to 0.39). Jaccard: shared={Goldfish,Rainbow}=2,
  union=3 → d_J = 1 − 2/3 = **1/3 ≈ 0.3333**. ✓
- **Independent vector example** x=[1,2,3,0], y=[1,0,3,4]: Σmin = 1+0+3+0 = 4, Σx=6, Σy=8 →
  BC = 1 − 2·4/14 = 6/14 ≈ **0.4286**; Σ|x−y| = 0+2+0+4 = 6, /14 ≈ 0.4286 ✓ (both forms agree).
  Presence A={1,2,3}, B={1,3,4}, ∩=2, ∪=4 → J=0.5, **distance 0.5**. ✓

### Findings
None. Description and Evidence are accurate and faithfully sourced.

## Stage B — Implementation

### Code path
`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs`
- `CalculateBetaDiversity` (lines 389–423): builds `allSpecies = union(keys)`, counts shared/unique using a `>0` guard.
- `CalculateBrayCurtis` (425–442): `sumMin += min(a1,a2)`, `sumTotal += a1+a2`; returns `sumTotal>0 ? 1 − 2·sumMin/sumTotal : 0`.
- `CalculateJaccardDistance` (444–448): `total = shared+unique1+unique2`; returns `total>0 ? 1 − shared/total : 0` (= 1 − J).

### Formula realised correctly?
Yes. BC matches `1 − 2·Σmin/Σtotal` exactly. Jaccard returns `1 − |A∩B|/|A∪B|` (distance), correctly labelled
`JaccardDistance` on the record. Zero-abundance keys contribute 0 to BC sums (harmless) and are excluded from
presence counts via the `>0` guard — consistent with spec. Symmetry holds (min and sum are symmetric in the two args).
All-zero / empty inputs guarded by `sumTotal>0` and `total>0`, returning 0 (no div-by-zero).

### Cross-verification table (recomputed vs code, all pass)
| Case | BC (code) | Jaccard dist (code) | Source/hand value |
|------|-----------|---------------------|-------------------|
| Wikipedia aquarium | 13/33 ≈ 0.3939 | 1/3 ≈ 0.3333 | ✓ Wikipedia |
| Identical {A:.4,B:.6} | 0 | 0 | ✓ |
| Disjoint {A} vs {B} | 1 | 1 | ✓ |
| Symmetry {A:.3,B:.7}/{B:.4,C:.6} | 3/5 | 2/3 | ✓ |
| Both empty | 0 | 0 | ✓ (convention) |
| One empty | 1 | 1 | ✓ |
| Single same diff abund .8/.3 | 5/11 | 0 | ✓ |
| Zero-abundance handling | 0.7 | 2/3 | ✓ |
| Species-count multi | 0.7 | 3/5 | ✓ |
| Abundance skew | 0.4 | (= balanced) | ✓ |

### Variant/delegate consistency
Single canonical method; no `*Fast` variant or instance-property duplicate to reconcile.

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_BetaDiversity_Tests.cs` — 22 tests, all asserting
**exact** sourced/hand-computed values within 1e-10 (not "no-throw" tautologies). Covers every Stage-A edge case:
published example, identity, disjoint, symmetry, empty (both/one), single-species, zero-abundance, counts,
abundance-vs-presence, name preservation.

### Findings / defects
None.

## Verdict & follow-ups

- **Stage A: PASS** — both formulas match Wikipedia and the cited primary literature exactly; Jaccard reported as
  distance and correctly labelled.
- **Stage B: PASS** — code realises the validated formulas; all hand-computed cross-checks match; edge cases and
  symmetry hold.
- **State: CLEAN** — no defects; build succeeds; BetaDiversity filter 22/22 pass; full suite 4484/4484 pass.
- No code changed.
