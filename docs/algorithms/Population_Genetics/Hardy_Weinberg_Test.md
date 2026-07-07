# Hardy-Weinberg Equilibrium Test

> **Baseline / reference method.** The chi-square Hardy-Weinberg test is a fast, deterministic QC screen; exact tests can be preferable for sparse counts. See [Legacy / Baseline Methods](../CANONICAL_MAP.md).

| Field | Value |
|-------|-------|
| Algorithm Group | Population Genetics |
| Test Unit ID | POP-HW-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

The Hardy-Weinberg equilibrium (HWE) test asks whether observed biallelic genotype counts are consistent with the expected `p^2 : 2pq : q^2` proportions under the classical random-mating model.[1][2][5] The repository implements a chi-square goodness-of-fit procedure with one degree of freedom and returns expected counts, chi-square, p-value, and a boolean decision at a caller-supplied significance level.[3][7][8][9] This makes the implementation deterministic and easy to interpret, but narrower than exact HWE procedures used for small-sample or sparse-count settings.[3][6]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The Hardy-Weinberg principle states that allele and genotype frequencies remain constant across generations under a specific set of modeling assumptions.[1][2][5] In practice, departures from the expected genotype proportions are used as a screening signal for non-random mating, population structure, selection, or technical artifacts such as genotyping error.[3][5][7]

### 2.2 Core Model

For a biallelic locus with allele frequencies `p` and `q`, where `p + q = 1`, the expected genotype frequencies are:[1][2][5]

$$
f(AA) = p^2, \quad f(Aa) = 2pq, \quad f(aa) = q^2
$$

From observed genotype counts `n_AA`, `n_Aa`, and `n_aa`, the allele frequencies are estimated as:[5][7]

$$
p = \frac{2n_{AA} + n_{Aa}}{2n}
$$

$$
q = 1 - p = \frac{2n_{aa} + n_{Aa}}{2n}
$$

where `n = n_AA + n_Aa + n_aa`.[5]

The chi-square goodness-of-fit statistic is:[3][5]

$$
\chi^2 = \sum_i \frac{(O_i - E_i)^2}{E_i}
$$

which expands in the biallelic HWE case to

$$
\chi^2 = \frac{(n_{AA} - p^2 n)^2}{p^2 n} + \frac{(n_{Aa} - 2pqn)^2}{2pqn} + \frac{(n_{aa} - q^2 n)^2}{q^2 n}
$$

For the biallelic test documented here, the degrees of freedom are `1 = 3 genotypes - 2 alleles`.[5][7]

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Mating is random with respect to the locus under study | Heterozygote or homozygote proportions can shift away from `p^2 : 2pq : q^2`.[5] |
| ASM-02 | Genotypes have equal fitness | Selection can drive systematic deviations from HWE.[5] |
| ASM-03 | Allele frequencies are not altered by mutation during the modeled interval | Mutation can change `p` and `q` across generations.[5] |
| ASM-04 | The population is effectively closed to migration | Admixture or structure can induce departures from the expected genotype mix.[5] |
| ASM-05 | The population is large enough that drift is negligible in the model | Finite-population drift can move observed counts away from HWE expectations.[5] |
| ASM-06 | Allele frequencies are the same in both sexes | Sex-specific frequency differences can violate the simple pooled expectation.[5] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `p + q = 1` | The two allele frequencies partition the full allele count at the locus.[5] |
| INV-02 | `p^2 + 2pq + q^2 = 1` | The expected genotype frequencies are the expansion of `(p + q)^2`.[1][2][5] |
| INV-03 | Expected genotype counts sum to `n` | Each expected count is one genotype frequency multiplied by the same sample size `n`.[5][7] |
| INV-04 | `chi-square >= 0` | Each term in the chi-square sum is a squared deviation divided by a positive expected count.[3][7] |
| INV-05 | The chi-square tail probability lies in `[0, 1]` | It is a probability derived from the chi-square distribution.[3][7][8] |

### 2.5 Comparison with Related Methods

| Aspect | Chi-square HWE test | Exact HWE test |
|--------|----------------------|----------------|
| Decision basis | Asymptotic chi-square approximation with `df = 1` | Exact tail probability under the null model |
| Typical use | Larger samples and non-sparse expected counts | Small samples or sparse genotype counts |
| Reference in current materials | Primary method documented here.[3][5] | Referenced as an alternative in the source list.[6] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `variantId` | `string` | required | Identifier carried through to the result record | Preserved in the output.[7][8][9] |
| `observedAA` | `int` | required | Observed count of `AA` genotypes | Documented contract assumes non-negative genotype counts.[5][7] |
| `observedAa` | `int` | required | Observed count of `Aa` genotypes | Documented contract assumes non-negative genotype counts.[5][7] |
| `observedaa` | `int` | required | Observed count of `aa` genotypes | Documented contract assumes non-negative genotype counts.[5][7] |
| `significanceLevel` | `double` | `0.05` | Decision threshold for `InEquilibrium` | Default is `0.05` in the repository API.[7][8][9] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `VariantId` | `string` | Input identifier echoed in the result.[7][8][9] |
| `ObservedAA` | `int` | Observed `AA` count passed into the test.[9] |
| `ObservedAa` | `int` | Observed `Aa` count passed into the test.[9] |
| `Observedaa` | `int` | Observed `aa` count passed into the test.[9] |
| `ExpectedAA` | `double` | Expected `AA` count `p^2 n` under HWE.[5][9] |
| `ExpectedAa` | `double` | Expected `Aa` count `2pqn` under HWE.[5][9] |
| `Expectedaa` | `double` | Expected `aa` count `q^2 n` under HWE.[5][9] |
| `ChiSquare` | `double` | Chi-square goodness-of-fit statistic.[3][9] |
| `PValue` | `double` | Upper-tail probability from the chi-square distribution with `df = 1`.[3][8][9] |
| `InEquilibrium` | `bool` | `true` when `PValue >= significanceLevel` in the current implementation.[7][8][9] |

### 3.3 Preconditions and Validation

The documented contract assumes non-negative genotype counts for one biallelic locus.[5][7] When `n = 0`, the repository returns `ChiSquare = 0`, `PValue = 1`, and `InEquilibrium = true`.[7][8][9] If any expected count is zero, that term is omitted from the chi-square sum to avoid division by zero.[7][9] The default significance level is `0.05`.[7][8][9]

## 4. Algorithm

### 4.1 High-Level Steps

1. Sum the three observed genotype counts to obtain `n`.
2. Return an equilibrium result with `ChiSquare = 0` and `PValue = 1` when `n = 0`.
3. Compute allele frequencies `p` and `q` from the observed genotype counts.
4. Compute expected counts `p^2 n`, `2pqn`, and `q^2 n`.
5. Accumulate the chi-square statistic over the genotype categories whose expected counts are positive.
6. Convert the chi-square statistic to a `df = 1` p-value.
7. Set `InEquilibrium` according to whether `PValue >= significanceLevel`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

For the `df = 1` chi-square test documented in the original file, the commonly cited critical values are:[5][7]

| Significance level | Critical chi-square |
|--------------------|---------------------|
| `0.10` | `2.706` |
| `0.05` | `3.841` |
| `0.01` | `6.635` |

The repository exposes the more general p-value form of the same decision rule by returning `PValue` and comparing it with the caller's `significanceLevel`.[7][8][9]

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `TestHardyWeinberg` | `O(1)` | `O(1)` | Uses constant-time arithmetic on one observed genotype triple.[7][9] |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)

- `PopulationGeneticsAnalyzer.TestHardyWeinberg(string, int, int, int, double)`: Computes expected counts, chi-square, p-value, and the equilibrium decision.
- `PopulationGeneticsAnalyzer.ChiSquareCDF(double, int)`: Private helper used for the chi-square distribution CDF.
- `PopulationGeneticsAnalyzer.RegularizedGammaP(double, double)`: Private helper used inside the chi-square CDF implementation.

### 5.2 Current Behavior

The repository computes `PValue` as `1 - ChiSquareCDF(chiSquare, 1)`, where the CDF path is implemented through a regularized lower incomplete gamma helper.[9] Terms with expected count `0` are skipped in the chi-square accumulation.[9] `VariantId` is copied directly into the output record, and `InEquilibrium` is defined strictly as `PValue >= significanceLevel`.[7][8][9] The current method does not perform explicit validation of negative genotype counts or of whether `significanceLevel` lies inside `[0, 1]`.[9]

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Allele-frequency estimation from genotype counts using the standard biallelic formulas.[5][9]
- Expected counts `p^2 n`, `2pqn`, and `q^2 n` under Hardy-Weinberg equilibrium.[1][2][5][9]
- Chi-square goodness-of-fit testing and p-value-based decision at the chosen significance level.[3][7][8][9]

**Intentionally simplified:**

- The repository implements the asymptotic chi-square test for one biallelic locus only; **consequence:** small-sample exact HWE procedures are not available in this API.

**Not implemented:**

- Exact HWE tests such as Wigginton, Cutler, and Abecasis (2005); **users should rely on:** no current alternative.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Non-negative genotype counts are assumed but not validated explicitly | Assumption | Negative counts can produce uninterpretable frequencies and expectations | accepted | The public method accepts `int` values directly and does not guard against negatives.[9] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `n = 0` | Returns `ChiSquare = 0`, `PValue = 1`, `InEquilibrium = true` | The implementation treats no data as no evidence against the null model.[7][8][9] |
| Single monomorphic sample `(1, 0, 0)` | Returns expected counts `(1, 0, 0)` with `ChiSquare = 0` | Observed and expected counts are identical.[7][8] |
| Fixed allele `(100, 0, 0)` or `(0, 0, 100)` | Returns `ChiSquare = 0` and `InEquilibrium = true` | The monomorphic expected counts match the observations exactly.[7][8] |
| Expected genotype count equals `0` | That genotype category is skipped in the chi-square sum | Avoids division by zero in the statistic.[7][9] |
| Borderline p-value case | `InEquilibrium` changes when `significanceLevel` changes | The repository compares `PValue` directly with the supplied threshold.[7][8] |

### 6.2 Limitations

The HWE test reports statistical deviation but not the cause of that deviation.[3][5] Biological causes such as inbreeding, population structure, selection, and non-random mating, as well as technical causes such as genotyping error, null alleles, or contamination, can all lead to rejection of equilibrium.[5] The repository does not provide exact HWE tests or diagnostics that separate biological from technical explanations.[6][9]

## 7. Examples and Related Material

### 7.1 Worked Example

The repository tests preserve Ford's scarlet tiger moth example from the HWE references: `1469` `AA`, `138` `Aa`, and `5` `aa` individuals, for `n = 1612`.[4][5][7][8] The resulting expected counts are approximately `1467.40`, `141.21`, and `3.40`, and the chi-square statistic is approximately `0.8309`.[7][8] Because that value is below the `0.05` critical value `3.841`, the sample is not rejected as out of equilibrium at `alpha = 0.05`.[5][7][8]

## 8. References

1. Hardy, G. H. 1908. Mendelian Proportions in a Mixed Population. Science 28(706):49-50.
2. Weinberg, W. 1908. Uber den Nachweis der Vererbung beim Menschen. Jahreshefte des Vereins fur vaterlandische Naturkunde in Wurttemberg 64:368-382.
3. Emigh, T. H. 1980. A Comparison of Tests for Hardy-Weinberg Equilibrium. Biometrics 36(4):627-642.
4. Ford, E. B. 1971. Ecological Genetics. 3rd ed. Chapman and Hall.
5. Wikipedia contributors. Hardy-Weinberg principle. https://en.wikipedia.org/wiki/Hardy-Weinberg_principle
6. Wigginton, J. E., D. J. Cutler, and G. R. Abecasis. 2005. A Note on Exact Tests of Hardy-Weinberg Equilibrium. American Journal of Human Genetics 76(5):887-893.
7. [POP-HW-001.md](../../../tests/TestSpecs/POP-HW-001.md)
8. [PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs)
9. [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)
