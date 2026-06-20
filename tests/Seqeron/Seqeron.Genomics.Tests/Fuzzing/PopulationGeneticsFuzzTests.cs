using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Population;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the PopGen area — allele frequency (POP-FREQ-001),
/// nucleotide diversity π (POP-DIV-001), Hardy-Weinberg equilibrium
/// (POP-HW-001) and Wright's fixation index F_ST (POP-FST-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (DivideByZeroException, OverflowException,
/// NullReferenceException, …) and no invalid *numeric* result (NaN/±Infinity, a
/// frequency outside its mathematical range, a pair that fails to sum to 1).
/// Every input must result in EITHER a well-defined, theory-correct value, OR a
/// *documented, intentional* validation exception
/// (ArgumentException / ArgumentOutOfRangeException). A raw runtime exception, a
/// silent NaN/garbage frequency, or a hang on garbage input is a bug, not a
/// passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-FREQ-001 — allele frequency (PopGen)
/// Checklist: docs/checklists/03_FUZZING.md, row 43.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the degenerate population boundaries that
///           stress the frequency arithmetic: 0 samples (the DivideByZero
///           boundary — total allele count == 0), 1 individual / 1 allele
///           (the monomorphic, fully-fixed locus), "all the same allele"
///           (every sampled chromosome the same → frequency 1.0 for it, 0 for
///           the other), and negative counts (-1, the biologically impossible
///           input that must be rejected or carried as a *defined* — not a
///           silently-garbage — result).
/// — docs/checklists/03_FUZZING.md §Description ("BE = граничні значення: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The allele-frequency contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Allele frequency is the relative abundance of an allele at one diploid
/// biallelic locus: each individual contributes two chromosome copies, so for
/// genotype counts n_AA (homozygous major), n_Aa (heterozygous), n_aa
/// (homozygous minor):
///     p = (2·n_AA + n_Aa) / (2·(n_AA + n_Aa + n_aa))
///     q = (2·n_aa + n_Aa) / (2·(n_AA + n_Aa + n_aa))
/// and the minor allele frequency MAF = min(f, 1 − f).
///   — docs/algorithms/Population_Genetics/Allele_Frequency.md §2.2, §2.4
///     (INV-01: p + q = 1 for any non-empty sample — THE key invariant;
///      INV-02: 0 ≤ p,q ≤ 1; INV-03: 0 ≤ MAF ≤ 0.5; INV-04: major + minor copies
///      equal the total allele count). Sources: Gillespie (2004) [1];
///      Wikipedia "Allele frequency" [2]; "Genotype frequency" [4].
///
/// POP-FREQ-001 spans TWO documented entry points with DIFFERENT, intentional
/// validation contracts. Fuzzing pins both, and the boundary between them, so
/// neither can silently drift:
///
/// (1) The STRICT genotype-count surface —
///     PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(int, int, int)
///     (src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs
///      lines 138–159). Documented validation (Allele_Frequency.md §3.1, §3.3,
///      §6.1):
///       • ANY negative genotype count (−1, int.MinValue, …) →
///         a *documented, intentional* ArgumentOutOfRangeException — the
///         biologically-impossible-input gate, NOT a silent negative frequency
///         and NOT a crash. This is the KEY negative-counts fuzz target.
///       • 0 samples (all three counts 0) → total allele count == 0 → the
///         explicit `totalAlleles == 0` guard returns the DEFINED (0, 0); there
///         is NO division, so no DivideByZeroException and no NaN. This is the
///         KEY 0-samples fuzz concern.
///       • 1 individual, all the same allele (homozygous major) → (1.0, 0.0);
///         all homozygous minor → (0.0, 1.0); all heterozygous → (0.5, 0.5).
///       • for ANY non-empty, non-negative count triple the result obeys
///         p + q = 1 exactly (INV-01), both in [0, 1] (INV-02), even at scale
///         (large counts must not overflow the int allele-copy arithmetic into a
///         wrong/negative total). The total allele copies = 2·(sum); large
///         inputs are kept under the int boundary so 2·sum cannot overflow — we
///         pin the in-range arithmetic, not an overflowing one.
///
/// (2) The LENIENT genotype-code surface —
///     PopulationGeneticsAnalyzer.CalculateMAF(IEnumerable&lt;int&gt;)
///     (PopulationGeneticsAnalyzer.cs lines 164–177). LENIENT by documented
///     design (Allele_Frequency.md §3.3, §5.2): it assumes the 0/1/2 diploid
///     genotype encoding (0 = hom ref, 1 = het, 2 = hom alt) and does NOT
///     validate each element. Documented behavior:
///       • empty genotype list → 0 (explicit short-circuit; no division by zero).
///       • monomorphic list (all 0 → no alt copies; or all 2 → all alt copies) →
///         MAF 0: a fixed locus has no minor allele (§6.1).
///       • well-formed 0/1/2 input → MAF = min(f, 1 − f) ∈ [0, 0.5] (INV-03).
///       • because the surface does NOT validate, an OUT-OF-ENCODING code (a
///         negative genotype, or a code &gt; 2) is summed verbatim. Fuzzing pins
///         the *actual defined* consequence — the method stays total: it does
///         NOT throw, does NOT divide by zero (the list is non-empty so the
///         denominator 2·count &gt; 0), and does NOT produce NaN/Infinity. It may
///         return a value OUTSIDE [0, 0.5] for such malformed input — a
///         documented consequence of the "caller must normalize" contract
///         (§5.3 "Intentionally simplified"), NOT a crash. We pin that this
///         surface never crashes and never NaNs on the malformed-code boundary,
///         so the strict/lenient split is explicit and cannot drift into a
///         silent DivideByZero.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-DIV-001 — nucleotide diversity π (PopGen)
/// Checklist: docs/checklists/03_FUZZING.md, row 44.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the degenerate sample boundaries that stress
///           the pairwise-diversity arithmetic: 0 sequences and 1 sequence (the
///           n &lt; 2 boundary — no unordered pair exists, so the number of pairs
///           C(n,2) = n(n−1)/2 is 0 → an unguarded denominator would
///           DivideByZero/NaN), "all sequences identical" (no variation → every
///           pairwise difference d_ij is 0 → π must be exactly 0 — THE key
///           identity), and very long sequences (the O(n²·L) cost boundary — must
///           still complete and never hang).
/// — docs/checklists/03_FUZZING.md §Description ("BE = граничні значення: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The nucleotide-diversity contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Nucleotide diversity π is the average number of nucleotide differences PER
/// SITE between all unordered pairs of aligned sequences:
///     π = ( Σ_{i&lt;j} d_ij ) / ( C(n,2) · L )
/// where d_ij is the number of differing positions between sequences i and j,
/// n is the number of sequences, and L is the (common, aligned) sequence length.
///   — docs/algorithms/Population_Genetics/Diversity_Statistics.md §2.2, §2.4
///     (INV-01: π ≥ 0 — it is a normalized count of differences; INV-02: π = 0
///      for identical sequences — every d_ij is 0). Sources: Nei &amp; Li (1979)
///     [1]; Tajima (1989) [3]; Wikipedia "Tajima's D" [5].
///
/// Entry point under test —
///   PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(
///       IEnumerable&lt;IReadOnlyList&lt;char&gt;&gt;)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs
///    lines 205–234). Documented validation/edge behavior
///   (Diversity_Statistics.md §3.3, §6.1):
///     • n &lt; 2 sequences (0 or 1) → the explicit `seqList.Count &lt; 2` guard
///       returns the DEFINED 0 BEFORE any pair loop; there is NO C(n,2) division
///       (comparisons would be 0), so NO DivideByZeroException and NO NaN. This is
///       the KEY 0-/1-sequence fuzz concern — the C(n,2) = 0 boundary.
///     • all sequences identical → every d_ij = 0 → totalDiff = 0 → π = 0 exactly,
///       independent of how many sequences or how long (INV-02).
///     • two sequences differing at every site → π = 1.0 (each of the L sites is
///       one mismatch in the single pair).
///     • very long sequences (n small, L large) → the O(n²·L) scan completes
///       within a bounded time; pinned under `[CancelAfter]` so a hang is a
///       failure, not an infinite wait.
///     • π is ALWAYS ≥ 0 and finite for valid aligned input (INV-01) — a count of
///       mismatches over a positive denominator can never be negative, NaN, or
///       ±Infinity.
///   The implementation does NOT validate equal length (§5.4): callers must supply
///   aligned, equal-length sequences. For n ≥ 2 the denominator is C(n,2)·L; with
///   a positive common L this is &gt; 0, so the division is well-defined.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-HW-001 — Hardy-Weinberg equilibrium (PopGen)
/// Checklist: docs/checklists/03_FUZZING.md, row 45.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the degenerate genotype-count boundaries
///           that stress the allele-frequency arithmetic and the chi-square
///           goodness-of-fit: 0 genotypes (n = 0 → the DivideByZero boundary in
///           BOTH the p/q estimate 1/(2n) AND every expected-count term, whose
///           expected category is 0); a single genotype (the minimal non-empty
///           sample); an all-heterozygous sample (the classic HWE-violation case
///           — observed het excess vs the p=q=0.5 expectation); and an
///           all-homozygous sample (a het DEFICIT, the inbreeding-pattern
///           departure, and the monomorphic fully-fixed sub-case where an
///           expected category is exactly 0).
/// — docs/checklists/03_FUZZING.md §Description ("BE = граничні значення: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The Hardy-Weinberg contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// For a biallelic locus the test estimates allele frequencies from observed
/// genotype counts (n_AA, n_Aa, n_aa), n = n_AA + n_Aa + n_aa:
///     p = (2·n_AA + n_Aa) / (2n),   q = 1 − p
/// computes the HWE-expected counts E_AA = p²n, E_Aa = 2pqn, E_aa = q²n, and a
/// chi-square goodness-of-fit statistic with df = 1:
///     χ² = Σ_i (O_i − E_i)² / E_i   (only over categories with E_i &gt; 0)
/// returning a HardyWeinbergResult(p/q-derived expected counts, χ², p-value,
/// InEquilibrium = PValue ≥ significanceLevel).
///   — docs/algorithms/Population_Genetics/Hardy_Weinberg_Test.md §2.2, §2.4
///     (INV-01: p + q = 1; INV-03: expected counts sum to n; INV-04: χ² ≥ 0;
///      INV-05: p-value ∈ [0, 1]). Sources: Hardy (1908) [1]; Weinberg (1908) [2];
///      Emigh (1980) [3]; Wikipedia "Hardy-Weinberg principle" [5].
///
/// Entry point under test —
///   PopulationGeneticsAnalyzer.TestHardyWeinberg(string, int, int, int, double)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs
///    lines 424–471). Documented validation/edge behavior
///   (Hardy_Weinberg_Test.md §3.3, §5.2, §6.1):
///     • n = 0 (0 genotypes) → the explicit `n == 0` guard returns the DEFINED
///       result (χ² = 0, PValue = 1, InEquilibrium = true) BEFORE any 1/(2n)
///       division; there is NO DivideByZeroException and NO NaN. This is the KEY
///       0-genotypes fuzz concern (Hardy_Weinberg_Test.md §3.3, §6.1).
///     • an expected category equal to 0 (e.g. the monomorphic fully-fixed
///       sample, where E_Aa = E_aa = 0) → that χ² term is SKIPPED by the explicit
///       `expected > 0` guard, so the 0-expected division is never performed —
///       NO DivideByZero, NO NaN (Hardy_Weinberg_Test.md §5.2, §6.1). This is the
///       KEY 0-expected fuzz concern.
///     • a single genotype (e.g. (1,0,0)) → a DEFINED result: monomorphic, so
///       observed equals expected → χ² = 0, InEquilibrium = true (§6.1).
///     • all heterozygous (0, N, 0) → p = q = 0.5, expected (0.25N, 0.5N, 0.25N);
///       the observed het EXCESS gives χ² = N exactly (a meaningful, scaling,
///       nonzero departure — the classic HWE violation, §6.1 / unit test
///       HW-M9 pins N=100 → χ² = 100).
///     • all homozygous mixed (N, 0, N) → p = q = 0.5 again, expected
///       (0.25·2N, 0.5·2N, 0.25·2N); the observed het DEFICIT (inbreeding pattern)
///       gives χ² = 2N exactly (a meaningful nonzero departure).
///     • all homozygous fixed (N, 0, 0) or (0, 0, N) → p = 1 or 0, expected equals
///       observed → χ² = 0, InEquilibrium = true (the monomorphic case where two
///       expected categories are 0 and their χ² terms are skipped, §6.1).
///     • χ² is ALWAYS ≥ 0 and finite for any non-negative count triple (INV-04),
///       and the p-value stays in [0, 1] (INV-05).
///   The method does NOT validate negative genotype counts or significanceLevel
///   range (§5.4 accepted deviation #1); the fuzz targets here are all
///   non-negative, so we pin the defined behavior of the documented boundary.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-FST-001 — Wright's fixation index F_ST (PopGen)
/// Checklist: docs/checklists/03_FUZZING.md, row 46.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the degenerate population boundaries that
///           stress the F_ST variance/heterozygosity ratio: a SINGLE population
///           (one sequence supplied, the other empty → no pair to differentiate);
///           two IDENTICAL populations (no among-population variance → the KEY
///           F_ST = 0 identity); and 0 SAMPLES, in both senses — an EMPTY per-locus
///           sequence (the early guard) AND per-locus sample sizes of 0 (the
///           n1 + n2 == 0 → 0/0 NaN-poisoning boundary inside the weighted-mean
///           arithmetic). The monomorphic-locus sub-case (every population fixed for
///           the same allele → total heterozygosity H_т = 0 → the
///           DivideByZero-on-H_т=0 boundary the prompt flags) is also pinned.
/// — docs/checklists/03_FUZZING.md §Description ("BE = граничні значення: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The F_ST contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// F_ST (the fixation index) measures population differentiation. Wright's
/// variance-based two-population definition, summed across compared loci, is:
///     F_ST = σ²_S / [ p̄·(1 − p̄) ]   (ratio of summed terms)
/// where for each locus p̄ = (n1·p1 + n2·p2)/(n1 + n2) is the sample-size-weighted
/// mean allele frequency, σ²_S = [n1·(p1−p̄)² + n2·(p2−p̄)²]/(n1 + n2) is the
/// weighted among-population variance, and p̄·(1−p̄) is the (expected total)
/// heterozygosity H_т.  F_ST = 0 means no differentiation, F_ST = 1 complete
/// differentiation.
///   — docs/algorithms/Population_Genetics/F_Statistics.md §2.A, §4.A
///     (INV-FST-01: 0 ≤ F_ST ≤ 1; INV-FST-02: identical populations → F_ST = 0
///      because σ²_S = 0 at every locus — THE key identity; INV-FST-03: equal-size
///      fixed-opposite alleles p1=1, p2=0 at every locus → F_ST = 1). Sources:
///      Wright (1965) [2]; Wikipedia "Fixation index" [5].
///
/// Entry point under test —
///   PopulationGeneticsAnalyzer.CalculateFst(
///       IEnumerable&lt;(double AlleleFreq, int SampleSize)&gt;,
///       IEnumerable&lt;(double AlleleFreq, int SampleSize)&gt;)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs
///    lines 609–646). Documented validation/edge behavior (F_Statistics.md §3.3,
///    §5.2, §6.1):
///     • EITHER population sequence empty (the "single population / 0 loci"
///       boundary) → the explicit `pop1.Count == 0 || pop2.Count == 0` guard
///       returns the DEFINED 0 BEFORE any per-locus arithmetic; no DivideByZero,
///       no NaN (§3.3, §6.1; lines 616–617). This is the KEY "single population"
///       fuzz concern — with only one group there is nothing to differentiate.
///     • mismatched per-locus counts → a *documented, intentional* ArgumentException
///       (lines 619–623): the two locus lists must align 1-to-1 (the implementation
///       compares like-for-like loci, not a truncated prefix at this surface).
///     • IDENTICAL populations (same per-locus AlleleFreq and SampleSize) → every
///       locus has σ²_S = 0 → numerator 0 over a positive H_т denominator →
///       F_ST = 0 EXACTLY (INV-FST-02; §6.1). This is THE key identity.
///     • 0 SAMPLES per locus (n1 = n2 = 0): p̄ = 0.0/0 is a 0/0 → NaN; the NaN
///       poisons `het` and hence the accumulated denominator, so the final
///       `denominator > 0` check is FALSE (NaN > 0 is false) and the method returns
///       the DEFINED 0 (line 645). The arithmetic is `double` division, so this is a
///       NaN path, NOT a DivideByZeroException — and the guard ensures NO NaN
///       ESCAPES to the caller. This is the KEY 0-samples fuzz concern.
///     • monomorphic locus across BOTH populations (every locus fixed for the same
///       allele, p1 = p2 ∈ {0, 1}) → σ²_S = 0 AND H_т = p̄(1−p̄) = 0 → the summed
///       denominator is 0 → the `denominator > 0` guard returns the DEFINED 0 with
///       NO DivideByZero (§6.1 "fixed for the same allele … returns 0"). This is the
///       H_т = 0 division boundary the prompt flags.
///     • for any valid two-population input F_ST stays in [0, 1] (INV-FST-01) and is
///       finite — a non-negative variance over a non-negative heterozygosity, both
///       guarded against the 0/0 case.
///   The method does NOT range-check allele frequencies or sample sizes (§3.3,
///   §5.2); the fuzz targets here exercise the documented degenerate boundaries and
///   pin that none of them crash, hang, NaN-escape, or leave [0, 1].
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-LD-001 — linkage disequilibrium (PopGen)
/// Checklist: docs/checklists/03_FUZZING.md, row 47 (the LAST PopGen unit).
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the degenerate locus/sample boundaries that
///           stress the r²/D' arithmetic and its two division denominators:
///           a SINGLE (monomorphic) locus — one of the two loci is invariant, so its
///           genotype variance is 0 → the r² denominator p_A·q_A·p_B·q_B contains a
///           0 factor → an unguarded r² = D²/0 is a DivideByZero / NaN (the KEY
///           monomorphic concern the prompt flags); MONOMORPHIC loci — BOTH loci
///           invariant, so BOTH variances are 0 and the genotype covariance is 0;
///           and 0 SAMPLES — an empty genotype-pair sequence, where there are no
///           haplotypes/individuals to estimate frequencies from, so every mean,
///           covariance and variance would be a 0/0 → the empty-input guard must
///           return the DEFINED (D'=0, r²=0) BEFORE any division.
/// — docs/checklists/03_FUZZING.md §Description ("BE = граничні значення: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The linkage-disequilibrium contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Linkage disequilibrium (LD) is the non-random association of alleles at two
/// loci. The classical coefficient and its two normalized summaries are:
///     D   = p_AB − p_A·p_B                       (disequilibrium coefficient)
///     D'  = D / D_max                            (Lewontin normalization, |D'| ∈ [0,1])
///     r²  = D² / (p_A·q_A·p_B·q_B)               (squared correlation, ∈ [0,1])
/// with q_A = 1 − p_A, q_B = 1 − p_B and
///     D_max = min(p_A·q_B, q_A·p_B)  if D ≥ 0,   min(p_A·p_B, q_A·q_B)  if D &lt; 0.
///   — docs/algorithms/Population_Genetics/Linkage_Disequilibrium.md §2.2, §2.4
///     (INV-02: 0 ≤ r² ≤ 1 — the KEY range; INV-03: 0 ≤ |D'| ≤ 1 — the KEY range;
///      INV-01: D = 0 under random association). Sources: Lewontin (1964) [1];
///      Hill &amp; Robertson (1968) [2]; Wikipedia "Linkage disequilibrium" [4].
///
/// Entry point under test —
///   PopulationGeneticsAnalyzer.CalculateLD(
///       string variant1Id, string variant2Id,
///       IEnumerable&lt;(int Geno1, int Geno2)&gt; genotypes, int distance)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs
///    lines 729–780). The repository works on UNPHASED diploid genotype pairs in the
///    documented 0/1/2 encoding (0 = hom major, 1 = het, 2 = hom minor), NOT on
///    explicit haplotype counts (§5.2, §5.4 deviation #1): it estimates r² as the
///    squared Pearson correlation of the two 0/1/2 genotype vectors and D from the
///    diploid covariance D = Cov(X₁,X₂)/2, then D' = |D|/D_max clamped to [0,1].
///   Documented validation/edge behavior (Linkage_Disequilibrium.md §3.3, §6.1):
///     • EMPTY genotype-pair sequence (0 samples) → the explicit `genoList.Count == 0`
///       guard returns the DEFINED LinkageDisequilibrium with DPrime = 0, RSquared = 0
///       (IDs and distance preserved) BEFORE any mean/covariance/variance — there is
///       NO 0/0, so NO DivideByZeroException and NO NaN (lines 737–740). This is the
///       KEY 0-samples fuzz concern.
///     • MONOMORPHIC locus (either genotype vector invariant → that locus's variance
///       is 0) → the r² denominator p_A·q_A·p_B·q_B has a 0 factor; the explicit
///       `var1 > 0 && var2 > 0` guard returns RSquared = 0 instead of dividing by 0
///       (line 760) — NO DivideByZero, NO NaN, NO ±Infinity. This is the KEY
///       monomorphic / "single locus" fuzz concern (the r² zero-denominator). The
///       same construction also covers BOTH loci monomorphic (both variances 0).
///     • a single genotype PAIR (n = 1) → both variances are 0 (a one-point sample
///       has no spread) → r² = 0, D' guarded to a finite value — a DEFINED, finite
///       (if statistically meaningless) result, never NaN (§6.1; M8/S1).
///     • PERFECT LD — identical 0/1/2 genotype vectors → Pearson correlation 1 →
///       r² = 1 EXACTLY and D' = 1 EXACTLY (§6.1 "Perfect LD … RSquared = 1, DPrime = 1").
///     • for ANY genotype-pair sequence the public DPrime ∈ [0, 1] (INV-03, the
///       implementation returns |D'| clamped to 1) and RSquared ∈ [0, 1] (INV-02) and
///       both are finite — a squared correlation over a guarded denominator can never
///       be negative, exceed 1, NaN, or ±Infinity.
///   The method does NOT validate the 0/1/2 encoding beyond the documented surface
///   (§3.3, §5.2); the fuzz targets here exercise the documented degenerate
///   boundaries and pin that none crash, hang, NaN-escape, or leave [0, 1].
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-ANCESTRY-001 — supervised / projection ADMIXTURE ancestry estimation (PopGen)
/// Checklist: docs/checklists/03_FUZZING.md, row 207.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the degenerate ancestry boundaries that
///           stress the FRAPPE EM ancestry update and its Σ_m q_m f_mj
///           denominators: a SINGLE source population (an individual that is a
///           pure member of one reference panel → its ancestry fraction for that
///           panel must be exactly 1 and 0 for every other — the q_k = 1 corner
///           of the simplex); an ADMIXED 50/50 individual (the symmetric fixed
///           point where two panels contribute equally → equal fractions of 0.5;
///           the K = 2 uniform-q fixed point, INV-04); and EMPTY input — empty
///           individuals OR empty reference panels → the explicit
///           `indList.Count == 0 || refList.Count == 0` guard yields an EMPTY
///           result with no EM, no division and no exception (the
///           "nothing to estimate" boundary, §3.3 / §6.1).
/// — docs/checklists/03_FUZZING.md §Description ("BE = граничні значення: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The ancestry-estimation contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Supervised / projection ADMIXTURE estimates, for each individual i, the
/// vector q_i of ancestry FRACTIONS — what fraction q_ik of i's genome derives
/// from each of the K reference populations — from biallelic SNP genotypes
/// g_ij ∈ {0,1,2} (copies of allele 1) given FIXED reference allele-1
/// frequencies f_kj. With F fixed, the FRAPPE EM update (Alexander, Novembre &amp;
/// Lange 2009, Eq. 4) is:
///     q_ik^{n+1} = (1/2J) Σ_j [ g_ij·a_ijk + (2−g_ij)·b_ijk ],
///     a_ijk = q_ik f_kj / Σ_m q_im f_mj,   b_ijk = q_ik (1−f_kj) / Σ_m q_im (1−f_mj),
/// where J is the count of informative (in-range) SNPs. Iteration stops once the
/// Eq. 2 binomial log-likelihood gain falls below ε = 10⁻⁴ (Eq. 5).
///   — docs/algorithms/Population_Genetics/Ancestry_Estimation.md §2.2, §2.4
///     (INV-01: Σ_k q_ik = 1 for every individual — THE key invariant, because
///      Eq. 4 divides by 2J and Σ_k(g·a + (2−g)·b) = g + (2−g) = 2 per SNP;
///      INV-02: 0 ≤ q_ik ≤ 1 — each fraction is a proportion; INV-03: EM is a
///      monotone log-likelihood ascent; INV-04: identical reference panels keep a
///      uniform q uniform — the 50/50 fixed point at K = 2). Source: Alexander,
///      Novembre &amp; Lange (2009), Genome Research 19(9):1655–1664 [1], Eqs. 2, 4, 5.
///
/// Entry point under test —
///   PopulationGeneticsAnalyzer.EstimateAncestry(
///       IEnumerable&lt;(string IndividualId, IReadOnlyList&lt;int&gt; Genotypes)&gt; individuals,
///       IEnumerable&lt;(string PopulationId, IReadOnlyList&lt;double&gt; AlleleFrequencies)&gt; referencePops,
///       int maxIterations = 100)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs
///    lines 1264–1291). Documented validation/edge behavior
///   (Ancestry_Estimation.md §3.3, §6.1):
///     • EMPTY individuals OR empty reference panels → the explicit
///       `indList.Count == 0 || refList.Count == 0` guard `yield break`s an EMPTY
///       result BEFORE any EM iteration (lines 1272–1273); there is NO 1/2J
///       division, NO exception — "nothing to estimate" (§6.1). This is the KEY
///       empty-input fuzz concern.
///     • SINGLE source population — an individual that is a pure member of one
///       panel (homozygous for the allele that panel is fixed for) → that panel's
///       fraction converges to 1 and every other to 0 (§7.1 diagnostic-individual
///       semantics; INV-02 corner of the simplex).
///     • ADMIXED 50/50 — with identical reference panels the uniform initial
///       q = (1/K,…,1/K) is a fixed point (INV-04): every fraction stays 1/K, i.e.
///       0.5 for K = 2, so the two sources contribute equally. The Eq. 4
///       numerators/denominators are proportional across panels, so the update
///       reproduces the uniform vector exactly.
///     • a genotype whose length ≠ the panel SNP count → that individual is
///       SKIPPED (lines 1280–1281); all-missing genotype (every code &lt; 0 or
///       &gt; 2) → the uniform prior is returned (no informative SNP, lines
///       1344–1345). No exception for these documented input classes (§3.3).
///     • for EVERY emitted individual the ancestry fractions sum to 1 EXACTLY
///       (INV-01), each lies in [0, 1] (INV-02), and none is NaN or ±Infinity —
///       the EM divides accumulated responsibilities by 2J, which keeps Σ_k q_ik
///       = 1 without a separate renormalization (§5.2). The Σ_m q_m f_mj
///       denominators are guarded (`mix &gt; 0 ? … : 0`, lines 1338–1339) so a
///       degenerate panel cannot inject a 0/0 NaN.
///   The fuzz targets here exercise the documented degenerate boundaries and pin
///   that none of them crash, hang, NaN-escape, break Σ_k q_ik = 1, or leave
///   [0, 1].
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-ROH-001 — runs of homozygosity (ROH) detection (PopGen)
/// Checklist: docs/checklists/03_FUZZING.md, row 208.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the degenerate genotype-track boundaries that
///           stress the consecutive-runs ROH scan and its two retention thresholds
///           (minSnps, minLength): an ALL-HETEROZYGOUS track (every genotype is the
///           opposite call 1 → no homozygous SNP can EVER seed a run → ZERO runs,
///           regardless of length — a heterozygous SNP never starts a homozygous
///           run, §6.1 "Leading heterozygotes"); an ALL-HOMOZYGOUS track (no
///           opposite genotype, no gap break → ONE maximal run spanning the whole
///           sorted track, retained iff it clears BOTH minSnps and minLength); and
///           the minLength EDGE (a run whose span End − Start is EXACTLY minLength
///           is retained, one bp below it is discarded — the `End − Start ≥
///           minLength` boundary), plus the twin minSnps edge (SnpCount EXACTLY
///           minSnps retained, one below discarded). The empty-track boundary
///           (0 SNPs → no runs) and the negative/zero argument gates (the BE "0"
///           and "−1" classes that must raise the DOCUMENTED
///           ArgumentOutOfRangeException, never a silent garbage run) are pinned too.
/// — docs/checklists/03_FUZZING.md §Description ("BE = граничні значення: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The ROH-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A run of homozygosity (ROH) is a contiguous stretch over which an individual is
/// homozygous at (nearly) every marker (McQuillan et al. 2008 [1]). FindROH detects
/// runs with the window-free *consecutive-runs* method (Marras et al. 2015 [3]):
/// SNPs are processed in ascending position order; a candidate run grows while it
/// holds at most `maxHeterozygotes` opposite (genotype == 1) calls and no inter-SNP
/// gap exceeds `maxGap`; crossing either limit closes the run at its LAST homozygous
/// SNP and restarts at the breaking SNP. A closed run is reported on the inclusive
/// [first homozygous SNP .. last homozygous SNP] interval ONLY IF it holds
/// SnpCount ≥ minSnps AND End − Start ≥ minLength (PLINK 1.9 `--homozyg-snp` /
/// `--homozyg-kb` thresholds [2]).
///   — docs/algorithms/PopGen/Runs_Of_Homozygosity.md §2.2 (Core Model), §2.4
///     (INV-01: every reported run has SnpCount ≥ minSnps and End − Start ≥
///      minLength — BOTH thresholds checked before emission; INV-02: a reported run
///      has ≤ maxHeterozygotes opposite genotypes and no gap > maxGap; INV-03: runs
///      emitted in ascending Start order with Start ≤ End), §3 (Contract), §6.1
///      (Edge Cases), §7.1 (Worked Example). Genotype 1 is the opposite call; any
///      other integer is treated as homozygous (no missing sentinel, §3.3).
///
/// Entry point under test —
///   PopulationGeneticsAnalyzer.FindROH(
///       IEnumerable&lt;(int Position, int Genotype)&gt; genotypes,
///       int minSnps = 100, int minLength = 1_000_000,
///       int maxHeterozygotes = 1, int maxGap = 1_000_000)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs
///    lines 1485–1581). Documented validation/edge behavior
///   (Runs_Of_Homozygosity.md §3.3, §6.1):
///     • null genotypes → a *documented, intentional* ArgumentNullException; and
///       minSnps &lt; 1, or NEGATIVE minLength / maxHeterozygotes / maxGap →
///       ArgumentOutOfRangeException — the BE "0"/"−1" validation gates (§3.3;
///       lines 1492–1501). These are eager, thrown BEFORE iteration (the deferred
///       iterator is wrapped, §5.2). This is the KEY negative/zero-argument concern.
///     • EMPTY genotype input → NO runs (nothing to scan, §6.1; line 1514).
///     • ALL-HETEROZYGOUS track (every genotype == 1) → NO runs: a het SNP cannot
///       seed a homozygous run, so snpCount never leaves 0 (§6.1 "Leading
///       heterozygotes"; lines 1560–1563). This is the KEY all-het boundary.
///     • ALL-HOMOZYGOUS track (no genotype == 1, consecutive gaps ≤ maxGap) → ONE
///       maximal run [firstPos .. lastPos] with SnpCount = number of SNPs, retained
///       iff SnpCount ≥ minSnps AND lastPos − firstPos ≥ minLength (INV-01). This is
///       the KEY all-homozygous boundary.
///     • minLength EDGE: a run whose span is EXACTLY minLength is RETAINED
///       (End − Start ≥ minLength is a ≥, §2.4 INV-01 / line 1532); one bp below the
///       threshold is DISCARDED. The twin minSnps edge: SnpCount EXACTLY minSnps is
///       retained, one below discarded ("passing count but not length, or
///       vice-versa → discarded", §6.1). This is the KEY minLength-edge boundary.
///     • UNSORTED input is sorted internally and yields the SAME runs (§6.1; line
///       1513). Output runs are non-overlapping, ascending in Start, Start ≤ End
///       (INV-03).
///   FindROH returns positions (int bp); there is no floating-point output, so the
///   fuzz NaN/Infinity concern lives entirely on the F_ROH companion (covered by
///   POP-ROH's CalculateInbreedingFromROH at the same surface; here we pin that
///   FindROH never crashes/hangs and that every emitted run obeys INV-01..INV-03 on
///   random valid tracks, under `[CancelAfter]`). Expected runs are derived BY HAND
///   from the documented contract (§2.2/§2.4), never read back off the code.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PopulationGeneticsFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>Deterministic RNG — seed fixed so generated fuzz inputs are reproducible.</summary>
    private static readonly Random Rng = new(20260620);

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  POP-FREQ-001 — allele frequency : fuzz targets
    //  Strict surface:  CalculateAlleleFrequencies(int, int, int)
    //  Lenient surface: CalculateMAF(IEnumerable<int>)
    // ═══════════════════════════════════════════════════════════════════

    #region POP-FREQ-001 — allele frequency

    #region Positive sanity — known genotype set yields the expected frequencies summing to 1

    /// <summary>
    /// Positive control: the worked example from the algorithm doc
    /// (Allele_Frequency.md §7.1, the four-o'clock plant) — 49 AA, 42 Aa, 9 aa
    /// → 200 allele copies, 140 of A, 60 of a → p(A) = 0.70, q(a) = 0.30. The
    /// frequencies must be exactly those values, must each lie in [0, 1]
    /// (INV-02), and must sum to 1 (INV-01). The same population, expressed as a
    /// 0/1/2 genotype list, must give MAF = min(0.70, 0.30) = 0.30 (INV-03).
    /// This anchors the fuzz battery: before pinning that garbage does no harm,
    /// confirm a known-good input gives the textbook-correct answer.
    /// </summary>
    [Test]
    public void AlleleFrequencies_KnownGenotypeSet_MatchesWorkedExampleAndSumsToOne()
    {
        var (p, q) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor: 49, heterozygous: 42, homozygousMinor: 9);

        p.Should().BeApproximately(0.70, Tolerance,
            because: "p(A) = (2·49 + 42)/200 = 140/200 = 0.70 (Allele_Frequency.md §7.1)");
        q.Should().BeApproximately(0.30, Tolerance,
            because: "q(a) = (2·9 + 42)/200 = 60/200 = 0.30 (Allele_Frequency.md §7.1)");

        p.Should().BeInRange(0.0, 1.0, "each allele frequency is a proportion (INV-02)");
        q.Should().BeInRange(0.0, 1.0, "each allele frequency is a proportion (INV-02)");
        (p + q).Should().BeApproximately(1.0, Tolerance,
            because: "every sampled chromosome is one allele or the other → p + q = 1 (INV-01)");

        // The same population as a 0/1/2 genotype list: 49×0, 42×1, 9×2.
        var genotypes = Enumerable.Repeat(0, 49)
            .Concat(Enumerable.Repeat(1, 42))
            .Concat(Enumerable.Repeat(2, 9));
        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        maf.Should().BeApproximately(0.30, Tolerance,
            because: "MAF = min(p, q) = min(0.70, 0.30) = 0.30 (INV-03)");
        maf.Should().BeInRange(0.0, 0.5, "MAF is the smaller of two complementary frequencies (INV-03)");
    }

    #endregion

    #region BE — 0 samples (the DivideByZero boundary: total allele count == 0)

    /// <summary>
    /// BE / DivideByZero boundary: with 0 samples every genotype count is 0, so
    /// the total allele count 2·(0+0+0) is 0. The contract must NOT divide by
    /// zero and must NOT return NaN — the explicit `totalAlleles == 0` guard
    /// returns the DEFINED (0, 0) (Allele_Frequency.md §3.3, §6.1;
    /// PopulationGeneticsAnalyzer.cs lines 150–153). This is the single most
    /// important fuzz concern for this unit: an unguarded denominator here is a
    /// DivideByZeroException crash in production.
    /// </summary>
    [Test]
    public void AlleleFrequencies_ZeroSamples_ReturnsDefinedZeroPairWithNoDivideByZero()
    {
        (double p, double q) result = default;
        var act = () => result = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(0, 0, 0);

        act.Should().NotThrow<DivideByZeroException>(
            "the zero-total-allele boundary is guarded; the denominator is never divided by zero");
        act.Should().NotThrow("0 samples is a defined boundary input (returns (0,0)), not an error");

        result.p.Should().Be(0.0, "no alleles are present → major frequency is the defined 0");
        result.q.Should().Be(0.0, "no alleles are present → minor frequency is the defined 0");
        double.IsNaN(result.p).Should().BeFalse("the guard avoids a 0/0 NaN");
        double.IsNaN(result.q).Should().BeFalse("the guard avoids a 0/0 NaN");
    }

    /// <summary>
    /// BE / DivideByZero boundary, lenient surface: an EMPTY genotype list is the
    /// 0-samples boundary for CalculateMAF. The denominator would be 2·0 = 0, so
    /// the contract short-circuits to the DEFINED 0 BEFORE any division
    /// (Allele_Frequency.md §3.3, §6.1; PopulationGeneticsAnalyzer.cs lines
    /// 168–169) — no DivideByZeroException, no NaN.
    /// </summary>
    [Test]
    public void Maf_EmptyGenotypeList_ReturnsDefinedZeroWithNoDivideByZero()
    {
        double maf = double.NaN;
        var act = () => maf = PopulationGeneticsAnalyzer.CalculateMAF(Array.Empty<int>());

        act.Should().NotThrow<DivideByZeroException>(
            "the empty-list boundary short-circuits before the 2·count denominator is used");
        act.Should().NotThrow("an empty genotype list is a defined boundary input, not an error");

        maf.Should().Be(0.0, "an empty sample has no minor allele frequency (defined 0)");
        double.IsNaN(maf).Should().BeFalse("the empty short-circuit avoids a 0/0 NaN");
    }

    #endregion

    #region BE — 1 allele / 1 individual (monomorphic, fully-fixed locus)

    /// <summary>
    /// BE: the minimal non-empty population — a single individual. The locus is
    /// monomorphic (fully fixed): one homozygous-major individual fixes the major
    /// allele at frequency 1.0 and the minor at 0.0; one homozygous-minor fixes
    /// the reverse; one heterozygote is the balanced 0.5/0.5. In every case the
    /// pair must sum to 1 (INV-01) and each frequency stays in [0, 1] (INV-02).
    /// This pins that the length-1 boundary neither off-by-ones the 2·count
    /// allele accounting nor escapes the [0,1] range.
    /// </summary>
    [TestCase(1, 0, 0, 1.0, 0.0, TestName = "AlleleFrequencies_SingleHomozygousMajor_IsFixedMajor")]
    [TestCase(0, 0, 1, 0.0, 1.0, TestName = "AlleleFrequencies_SingleHomozygousMinor_IsFixedMinor")]
    [TestCase(0, 1, 0, 0.5, 0.5, TestName = "AlleleFrequencies_SingleHeterozygote_IsBalanced")]
    public void AlleleFrequencies_SingleIndividual_IsMonomorphicOrBalanced(
        int hMajor, int het, int hMinor, double expectedP, double expectedQ)
    {
        var (p, q) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(hMajor, het, hMinor);

        p.Should().BeApproximately(expectedP, Tolerance,
            because: "the major-allele frequency of a single defined genotype is exact");
        q.Should().BeApproximately(expectedQ, Tolerance,
            because: "the minor-allele frequency of a single defined genotype is exact");

        p.Should().BeInRange(0.0, 1.0, "a frequency is a proportion (INV-02)");
        q.Should().BeInRange(0.0, 1.0, "a frequency is a proportion (INV-02)");
        (p + q).Should().BeApproximately(1.0, Tolerance, "p + q = 1 for any non-empty sample (INV-01)");
    }

    /// <summary>
    /// BE: monomorphic locus on the lenient MAF surface. A genotype list that is
    /// all homozygous reference (all 0) carries no alternate copies; a list that
    /// is all homozygous alternate (all 2) carries only alternate copies. Either
    /// way the locus is FIXED, so MAF = 0 — a fixed locus has no minor allele
    /// (Allele_Frequency.md §6.1 "Monomorphic genotype list → 0"). Verified for a
    /// single individual and a longer fixed run; MAF must stay in [0, 0.5].
    /// </summary>
    [TestCase(0, 1, TestName = "Maf_AllHomozygousReference_Single_IsZero")]
    [TestCase(0, 25, TestName = "Maf_AllHomozygousReference_Run_IsZero")]
    [TestCase(2, 1, TestName = "Maf_AllHomozygousAlternate_Single_IsZero")]
    [TestCase(2, 25, TestName = "Maf_AllHomozygousAlternate_Run_IsZero")]
    public void Maf_MonomorphicGenotypeList_IsZero(int code, int count)
    {
        double maf = PopulationGeneticsAnalyzer.CalculateMAF(Enumerable.Repeat(code, count));

        maf.Should().BeApproximately(0.0, Tolerance,
            because: "a fixed (monomorphic) locus has no minor allele → MAF = 0 (INV-03, §6.1)");
        maf.Should().BeInRange(0.0, 0.5, "MAF can never escape [0, 0.5] for valid 0/1/2 input (INV-03)");
    }

    #endregion

    #region BE — all the same allele (every chromosome copy identical → frequency 1.0)

    /// <summary>
    /// BE "all same allele": when every sampled chromosome carries the same
    /// allele, that allele's frequency is exactly 1.0 and the other is exactly
    /// 0.0. A population of N homozygous-major individuals fixes the major allele;
    /// N homozygous-minor fixes the minor. Verified across a fixed-seed sweep of
    /// population sizes so the result is independent of N, stays in [0, 1]
    /// (INV-02), and the pair sums to 1 (INV-01). This is the population-scale
    /// analogue of the homopolymer boundary.
    /// </summary>
    [Test]
    public void AlleleFrequencies_AllSameAllele_IsFrequencyOneForItZeroForOther()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int n = Rng.Next(1, 100_000);

            var (pMaj, qMaj) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(n, 0, 0);
            pMaj.Should().BeApproximately(1.0, Tolerance,
                because: $"all {n} individuals homozygous-major → major frequency 1.0 regardless of N");
            qMaj.Should().BeApproximately(0.0, Tolerance, "the absent allele has frequency 0");
            (pMaj + qMaj).Should().BeApproximately(1.0, Tolerance, "p + q = 1 (INV-01)");
            pMaj.Should().BeInRange(0.0, 1.0, "INV-02"); qMaj.Should().BeInRange(0.0, 1.0, "INV-02");

            var (pMin, qMin) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(0, 0, n);
            pMin.Should().BeApproximately(0.0, Tolerance, "the absent allele has frequency 0");
            qMin.Should().BeApproximately(1.0, Tolerance,
                because: $"all {n} individuals homozygous-minor → minor frequency 1.0 regardless of N");
            (pMin + qMin).Should().BeApproximately(1.0, Tolerance, "p + q = 1 (INV-01)");
            pMin.Should().BeInRange(0.0, 1.0, "INV-02"); qMin.Should().BeInRange(0.0, 1.0, "INV-02");
        }
    }

    /// <summary>
    /// BE "all same allele", lenient surface: a population that is entirely
    /// homozygous alternate (all genotype code 2) carries the alternate allele at
    /// frequency 1.0; the alternate-allele frequency f = 1 makes
    /// MAF = min(1, 0) = 0 — the locus is fixed, so there is no minor allele. The
    /// "all same" boundary must therefore yield MAF = 0, never a value outside
    /// [0, 0.5].
    /// </summary>
    [Test]
    public void Maf_AllSameAlternateAllele_IsZero()
    {
        for (int trial = 0; trial < 16; trial++)
        {
            int n = Rng.Next(1, 50_000);

            double maf = PopulationGeneticsAnalyzer.CalculateMAF(Enumerable.Repeat(2, n));

            maf.Should().BeApproximately(0.0, Tolerance,
                because: $"all {n} genotypes hom-alt → alt frequency 1.0 → MAF = min(1, 0) = 0");
            maf.Should().BeInRange(0.0, 0.5, "MAF ∈ [0, 0.5] (INV-03)");
        }
    }

    #endregion

    #region BE — negative counts (biologically impossible input)

    /// <summary>
    /// BE "−1": a negative genotype count is biologically impossible. The STRICT
    /// surface must reject it with the *documented, intentional*
    /// ArgumentOutOfRangeException — NOT a silent negative/garbage frequency and
    /// NOT a DivideByZero (Allele_Frequency.md §3.1, §3.3;
    /// PopulationGeneticsAnalyzer.cs lines 143–148). Each of the three count
    /// arguments is guarded independently, and the extreme int.MinValue is
    /// covered to pin that no overflow slips a negative past the guard.
    /// </summary>
    [TestCase(-1, 0, 0, TestName = "AlleleFrequencies_NegativeHomozygousMajor_Throws")]
    [TestCase(0, -1, 0, TestName = "AlleleFrequencies_NegativeHeterozygous_Throws")]
    [TestCase(0, 0, -1, TestName = "AlleleFrequencies_NegativeHomozygousMinor_Throws")]
    [TestCase(-5, 3, 2, TestName = "AlleleFrequencies_NegativeMajorAmongPositives_Throws")]
    [TestCase(int.MinValue, 0, 0, TestName = "AlleleFrequencies_IntMinValueMajor_Throws")]
    [TestCase(10, -10, 10, TestName = "AlleleFrequencies_NegativeHetAmongPositives_Throws")]
    public void AlleleFrequencies_NegativeCount_ThrowsDocumentedArgumentOutOfRange(
        int hMajor, int het, int hMinor)
    {
        var act = () => PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(hMajor, het, hMinor);

        act.Should().Throw<ArgumentOutOfRangeException>(
            "a negative genotype count is biologically impossible and is rejected at the validation " +
            "gate, never miscounted into a negative frequency and never crashed");
    }

    /// <summary>
    /// BE "−1", lenient surface: CalculateMAF does NOT validate its genotype codes
    /// (Allele_Frequency.md §3.3, §5.2/§5.3 "caller must normalize"). An
    /// out-of-encoding code (a negative genotype, or a code &gt; 2) is summed
    /// verbatim. The contract guarantee we pin here is the fuzz-safety one: the
    /// method must stay TOTAL — it must NOT throw, must NOT divide by zero (the
    /// list is non-empty, so the denominator 2·count &gt; 0), and must NEVER
    /// produce NaN or ±Infinity. The returned value MAY fall outside [0, 0.5] for
    /// such malformed input — that is the documented consequence of feeding
    /// out-of-encoding codes, not a crash. This makes the strict/lenient split
    /// explicit: the malformed input the strict surface would reject is here
    /// carried as a defined-but-out-of-range number, never a DivideByZero.
    /// </summary>
    [Test]
    public void Maf_OutOfEncodingGenotypeCodes_StaysFiniteAndNeverThrows()
    {
        // A deterministic sweep including negatives and codes > 2, mixed with the
        // valid 0/1/2 encoding. Every list is non-empty (denominator > 0).
        for (int trial = 0; trial < 64; trial++)
        {
            int length = Rng.Next(1, 50);
            var genotypes = new int[length];
            for (int i = 0; i < length; i++)
                genotypes[i] = Rng.Next(-3, 6); // spans -3..-1 (impossible), 0/1/2 (valid), 3..5 (>2)

            double maf = double.NaN;
            var act = () => maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

            act.Should().NotThrow(
                "the lenient MAF surface is total over int genotype codes — it never throws on garbage");
            double.IsNaN(maf).Should().BeFalse(
                "a non-empty list has a non-zero denominator, so the result is never a 0/0 NaN");
            double.IsInfinity(maf).Should().BeFalse(
                "the arithmetic on bounded int sums and a non-zero denominator never yields ±Infinity");
        }
    }

    /// <summary>
    /// BE "−1", explicit defined-result pin for the lenient surface: a single
    /// genotype code of −1 sums to altAlleles = −1 over totalAlleles = 2, giving
    /// altFreq = −0.5 and MAF = min(−0.5, 1.5) = −0.5. We pin this EXACT value to
    /// document that the surface does NOT silently coerce the impossible input
    /// into [0, 0.5] — it carries the malformed code through deterministically
    /// (the "caller must normalize" contract). The point is that the behavior is
    /// *defined and stable*, not that −0.5 is a valid frequency.
    /// </summary>
    [Test]
    public void Maf_SingleNegativeOneCode_ReturnsDefinedOutOfRangeValue_NotACrash()
    {
        double maf = double.NaN;
        var act = () => maf = PopulationGeneticsAnalyzer.CalculateMAF(new[] { -1 });

        act.Should().NotThrow("an out-of-encoding code is carried, not validated, by the lenient surface");
        maf.Should().BeApproximately(-0.5, Tolerance,
            because: "alt copies = -1 over 2 total → altFreq -0.5 → MAF min(-0.5, 1.5) = -0.5; " +
                     "a defined (if biologically meaningless) value, documenting the no-validation contract");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  POP-DIV-001 — nucleotide diversity π : fuzz targets
    //  Surface: CalculateNucleotideDiversity(IEnumerable<IReadOnlyList<char>>)
    // ═══════════════════════════════════════════════════════════════════

    #region POP-DIV-001 — nucleotide diversity

    #region Helpers — sequence builders

    /// <summary>Materializes a string into the IReadOnlyList&lt;char&gt; the diversity API expects.</summary>
    private static IReadOnlyList<char> Seq(string s) => s.ToCharArray();

    /// <summary>A deterministic random aligned ACGT sequence of the given length.</summary>
    private static IReadOnlyList<char> RandomSeq(int length)
    {
        const string alphabet = "ACGT";
        var buffer = new char[length];
        for (int i = 0; i < length; i++)
            buffer[i] = alphabet[Rng.Next(alphabet.Length)];
        return buffer;
    }

    #endregion

    #region Positive sanity — known sample yields the documented π and π ≥ 0

    /// <summary>
    /// Positive control: the worked example pinned in the algorithm doc and the
    /// unit tests (Diversity_Statistics.md §7.1; the Wikipedia Tajima's D dataset)
    /// — five aligned length-20 sequences with 20 total pairwise differences over
    /// C(5,2) = 10 unordered pairs → k̂ = 20/10 = 2.0 → π = k̂ / L = 2.0/20 = 0.1.
    /// The result must be EXACTLY 0.1 (per the documented formula
    /// π = Σ_{i&lt;j} d_ij / (C(n,2)·L)), must be ≥ 0 (INV-01), and must be finite.
    /// This anchors the fuzz battery: before pinning that boundary input does no
    /// harm, confirm a known-good sample gives the textbook-correct π.
    /// </summary>
    [Test]
    public void NucleotideDiversity_KnownSample_MatchesDocumentedPiAndIsNonNegative()
    {
        var sequences = new List<IReadOnlyList<char>>
        {
            Seq("00000000000000000000"), // Y
            Seq("00100000000010000010"), // A
            Seq("00000000000010000010"), // B
            Seq("00000010000000000010"), // C
            Seq("00000010000010000010"), // D
        };

        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        pi.Should().BeApproximately(0.1, Tolerance,
            because: "20 pairwise differences over C(5,2)=10 pairs → k̂=2.0 → π=2.0/20=0.1 " +
                     "(Diversity_Statistics.md §7.1, the documented formula)");
        pi.Should().BeGreaterThanOrEqualTo(0.0, "π is a normalized count of differences (INV-01)");
        double.IsNaN(pi).Should().BeFalse("a known-good sample never yields NaN");
        double.IsInfinity(pi).Should().BeFalse("a known-good sample never yields ±Infinity");

        // A second independent anchor: two sequences differing at EVERY site → π = 1.0
        // (each of the L sites contributes one mismatch in the single unordered pair).
        double piAllDiff = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(
            new List<IReadOnlyList<char>> { Seq("AAAA"), Seq("CCCC") });
        piAllDiff.Should().BeApproximately(1.0, Tolerance,
            because: "every site differs in the only pair → π = 4/(1·4) = 1.0 (§6.1)");
    }

    #endregion

    #region BE — 0 sequences (the C(n,2) = 0 division boundary)

    /// <summary>
    /// BE / DivideByZero boundary: with 0 sequences the number of unordered pairs
    /// C(0,2) = 0 and there is no length to normalize by. The contract must NOT
    /// divide by zero and must NOT return NaN — the explicit `seqList.Count &lt; 2`
    /// guard returns the DEFINED 0 BEFORE the pair loop ever runs
    /// (Diversity_Statistics.md §3.3, §6.1; PopulationGeneticsAnalyzer.cs lines
    /// 210–211). This is the single most important fuzz concern for this unit: an
    /// unguarded C(n,2)·L denominator here is a 0/0 NaN (or DivideByZero) in
    /// production.
    /// </summary>
    [Test]
    public void NucleotideDiversity_ZeroSequences_ReturnsDefinedZeroWithNoDivideByZero()
    {
        double pi = double.NaN;
        var act = () => pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(
            Array.Empty<IReadOnlyList<char>>());

        act.Should().NotThrow<DivideByZeroException>(
            "the empty-sample boundary short-circuits before the C(n,2)·L denominator is used");
        act.Should().NotThrow("0 sequences is a defined boundary input (returns 0), not an error");

        pi.Should().Be(0.0, "no sequences → no pairwise differences → π is the defined 0 (§6.1)");
        double.IsNaN(pi).Should().BeFalse("the n<2 guard avoids a 0/0 NaN");
    }

    #endregion

    #region BE — 1 sequence (no unordered pair exists → π = 0)

    /// <summary>
    /// BE: a single sequence — the minimal-but-still-degenerate sample. There is no
    /// unordered pair (C(1,2) = 0), so pairwise diversity is undefined for n &lt; 2;
    /// the contract returns the DEFINED 0 (Diversity_Statistics.md §3.3, §6.1).
    /// Verified across a fixed-seed sweep of lengths so the single-sequence boundary
    /// returns 0 regardless of the sequence content or length, with no
    /// DivideByZero and no NaN.
    /// </summary>
    [Test]
    public void NucleotideDiversity_SingleSequence_ReturnsDefinedZeroForAnyContent()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int length = Rng.Next(1, 256);
            var single = new List<IReadOnlyList<char>> { RandomSeq(length) };

            double pi = double.NaN;
            var act = () => pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(single);

            act.Should().NotThrow<DivideByZeroException>(
                "a single sequence has no pair, so the C(n,2) denominator is never reached");
            act.Should().NotThrow("one sequence is a defined boundary input, not an error");
            pi.Should().Be(0.0, $"no unordered pair exists for n=1 (length {length}) → defined π = 0");
            double.IsNaN(pi).Should().BeFalse("the n<2 guard avoids a 0/0 NaN");
        }
    }

    #endregion

    #region BE — all sequences identical (no variation → π = 0, the key identity)

    /// <summary>
    /// BE / KEY identity (INV-02): when every sequence in the sample is identical,
    /// every pairwise difference count d_ij is 0, so π = 0 EXACTLY — there is no
    /// variation to measure (Diversity_Statistics.md §2.4 INV-02, §6.1). Verified
    /// across a fixed-seed sweep of sample sizes and lengths so the identity holds
    /// independent of n and L. This is the diversity analogue of the "all same
    /// allele" boundary: a monomorphic alignment must read as zero diversity, never
    /// a spurious non-zero value, never NaN.
    /// </summary>
    [Test]
    public void NucleotideDiversity_AllIdenticalSequences_IsExactlyZero()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int n = Rng.Next(2, 40);          // at least a pair, so a real division happens
            int length = Rng.Next(1, 200);
            var template = RandomSeq(length);

            // n copies of the SAME content (a fresh array per copy, same characters).
            var sequences = Enumerable.Range(0, n)
                .Select(_ => (IReadOnlyList<char>)template.ToArray())
                .ToList();

            double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

            pi.Should().Be(0.0,
                because: $"every one of the C({n},2) pairs is identical (length {length}) → all d_ij = 0 → π = 0 (INV-02)");
            pi.Should().BeGreaterThanOrEqualTo(0.0, "π ≥ 0 always (INV-01)");
            double.IsNaN(pi).Should().BeFalse("an identical, non-empty sample divides 0 by a positive denominator, not 0/0");
        }
    }

    #endregion

    #region BE — random aligned samples are always non-negative and finite

    /// <summary>
    /// BE / INV-01 sweep: across a fixed-seed battery of random, equal-length
    /// aligned samples (varying n ≥ 2 and L ≥ 1) π must ALWAYS be a well-defined,
    /// finite, non-negative per-site rate — a count of mismatches over a positive
    /// C(n,2)·L denominator can never be negative, NaN, or ±Infinity, and since π
    /// is a per-site average difference rate it can never exceed 1.0. This pins the
    /// total-function contract over the random-but-aligned interior, complementing
    /// the degenerate-boundary tests.
    /// </summary>
    [Test]
    public void NucleotideDiversity_RandomAlignedSamples_AlwaysInUnitIntervalAndFinite()
    {
        for (int trial = 0; trial < 64; trial++)
        {
            int n = Rng.Next(2, 25);
            int length = Rng.Next(1, 120);
            var sequences = Enumerable.Range(0, n)
                .Select(_ => RandomSeq(length))
                .ToList();

            double pi = double.NaN;
            var act = () => pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

            act.Should().NotThrow("a random equal-length aligned sample is always a defined input");
            double.IsNaN(pi).Should().BeFalse("a positive C(n,2)·L denominator never yields a 0/0 NaN");
            double.IsInfinity(pi).Should().BeFalse("bounded mismatch counts over a positive denominator never give ±Infinity");
            pi.Should().BeGreaterThanOrEqualTo(0.0, "π is a normalized count of differences (INV-01)");
            pi.Should().BeLessThanOrEqualTo(1.0, "π is a per-site average mismatch rate, capped at one mismatch per site");
        }
    }

    #endregion

    #region BE — very long sequences (the O(n²·L) cost boundary must complete, not hang)

    /// <summary>
    /// BE / hang boundary: nucleotide diversity is O(n²·L). With a small sample of
    /// VERY LONG aligned sequences the pairwise scan must still COMPLETE within a
    /// bounded time and return a finite, in-range π — it must never hang
    /// (Diversity_Statistics.md §4.3). We keep n small (a handful of sequences) and
    /// L large (hundreds of thousands of sites) so the L dimension is stressed
    /// without an n² blow-up, and pin the run under `[CancelAfter]` so a hang is a
    /// hard test failure rather than an infinite wait. The constructed sample
    /// differs at exactly one known site, giving a tiny but exactly-pinned π that
    /// also proves the long scan is arithmetically correct, not just fast.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void NucleotideDiversity_VeryLongSequences_CompletesWithFinitePi()
    {
        const int length = 500_000;
        var a = new char[length];
        var b = new char[length];
        Array.Fill(a, 'A');
        Array.Fill(b, 'A');
        b[length / 2] = 'C'; // exactly one differing site across the single pair

        var sequences = new List<IReadOnlyList<char>> { a, b };

        double pi = double.NaN;
        var act = () => pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        act.Should().NotThrow("the long-sequence O(n²·L) scan must complete without crashing");
        double.IsNaN(pi).Should().BeFalse("a positive denominator never yields NaN even at scale");
        double.IsInfinity(pi).Should().BeFalse("the bounded difference count never overflows into ±Infinity");
        pi.Should().BeApproximately(1.0 / length, Tolerance,
            because: "exactly one of the L=500000 sites differs in the single pair → π = 1/L");
        pi.Should().BeGreaterThanOrEqualTo(0.0, "π ≥ 0 always (INV-01)");
    }

    /// <summary>
    /// BE / hang boundary, identical-at-scale: a small sample of very long but
    /// IDENTICAL sequences must complete and read as exactly π = 0 — the INV-02
    /// identity must survive the long O(n²·L) scan, proving the zero-diversity
    /// short-circuit is by content (all d_ij = 0), not an accident of small inputs.
    /// Pinned under `[CancelAfter]` so the long scan cannot silently hang.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void NucleotideDiversity_VeryLongIdenticalSequences_CompletesAtExactlyZero()
    {
        const int length = 400_000;
        var template = new char[length];
        Array.Fill(template, 'G');

        var sequences = new List<IReadOnlyList<char>>
        {
            (char[])template.Clone(),
            (char[])template.Clone(),
            (char[])template.Clone(),
        };

        double pi = double.NaN;
        var act = () => pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        act.Should().NotThrow("the long identical-sequence scan must complete without crashing");
        pi.Should().Be(0.0, "every one of the C(3,2)=3 pairs is identical at scale → π = 0 (INV-02)");
        double.IsNaN(pi).Should().BeFalse("an identical sample divides 0 by a positive denominator, not 0/0");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  POP-HW-001 — Hardy-Weinberg equilibrium : fuzz targets
    //  Surface: TestHardyWeinberg(string, int, int, int, double)
    // ═══════════════════════════════════════════════════════════════════

    #region POP-HW-001 — Hardy-Weinberg equilibrium

    #region Positive sanity — known genotype set yields the documented chi-square and HWE verdict

    /// <summary>
    /// Positive control: Ford's scarlet tiger moth dataset pinned in the algorithm
    /// doc and the unit tests (Hardy_Weinberg_Test.md §7.1) — 1469 AA, 138 Aa,
    /// 5 aa, n = 1612 → expected counts ≈ (1467.40, 141.21, 3.40) and χ² ≈ 0.8309,
    /// which is below the df=1 α=0.05 critical value 3.841, so the sample is NOT
    /// rejected (InEquilibrium = true). The result must match the documented χ²,
    /// the expected counts must sum to n (INV-03), p + q must equal 1 (INV-01),
    /// χ² must be ≥ 0 (INV-04) and the p-value must lie in [0, 1] (INV-05). This
    /// anchors the fuzz battery: before pinning that boundary input does no harm,
    /// confirm a known-good genotype set gives the textbook-correct verdict.
    /// </summary>
    [Test]
    public void HardyWeinberg_KnownGenotypeSet_MatchesDocumentedChiSquareAndVerdict()
    {
        var r = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "rs-moth", observedAA: 1469, observedAa: 138, observedaa: 5);

        r.ChiSquare.Should().BeApproximately(0.8309, 0.01,
            because: "Ford's moth data → χ² ≈ 0.8309 (Hardy_Weinberg_Test.md §7.1)");
        r.InEquilibrium.Should().BeTrue("χ² ≈ 0.83 < 3.841 → not rejected at α=0.05 (§7.1)");

        r.ExpectedAA.Should().BeApproximately(1467.40, 0.05, "E_AA = p²n (§7.1)");
        r.ExpectedAa.Should().BeApproximately(141.21, 0.05, "E_Aa = 2pqn (§7.1)");
        r.Expectedaa.Should().BeApproximately(3.40, 0.05, "E_aa = q²n (§7.1)");

        (r.ExpectedAA + r.ExpectedAa + r.Expectedaa).Should().BeApproximately(1612.0, Tolerance,
            because: "expected counts sum to n (INV-03)");

        double p = (2.0 * 1469 + 138) / (2.0 * 1612);
        double q = 1 - p;
        (p + q).Should().BeApproximately(1.0, Tolerance, "p + q = 1 (INV-01)");

        r.ChiSquare.Should().BeGreaterThanOrEqualTo(0.0, "χ² ≥ 0 (INV-04)");
        r.PValue.Should().BeInRange(0.0, 1.0, "the chi-square tail probability is a probability (INV-05)");
        double.IsNaN(r.ChiSquare).Should().BeFalse("a known-good sample never yields NaN");
    }

    #endregion

    #region BE — 0 genotypes (the DivideByZero boundary: n = 0)

    /// <summary>
    /// BE / DivideByZero boundary: with 0 genotypes n = 0, so BOTH the allele-
    /// frequency estimate p = (…)/(2n) AND every expected-count term would divide
    /// by zero. The contract must NOT divide by zero and must NOT return NaN — the
    /// explicit `n == 0` guard returns the DEFINED result (χ² = 0, PValue = 1,
    /// InEquilibrium = true) BEFORE any division (Hardy_Weinberg_Test.md §3.3,
    /// §6.1; PopulationGeneticsAnalyzer.cs lines 433–436). This is the single most
    /// important fuzz concern for this unit: an unguarded 1/(2n) here is a
    /// DivideByZero/NaN in production.
    /// </summary>
    [Test]
    public void HardyWeinberg_ZeroGenotypes_ReturnsDefinedEquilibriumWithNoDivideByZero()
    {
        PopulationGeneticsAnalyzer.HardyWeinbergResult r = default;
        var act = () => r = PopulationGeneticsAnalyzer.TestHardyWeinberg("rs-empty", 0, 0, 0);

        act.Should().NotThrow<DivideByZeroException>(
            "the n=0 boundary is guarded; neither 1/(2n) nor an expected-count term is divided by zero");
        act.Should().NotThrow("0 genotypes is a defined boundary input (no evidence against the null), not an error");

        r.ChiSquare.Should().Be(0.0, "no data → no deviation → χ² is the defined 0 (§6.1)");
        r.PValue.Should().Be(1.0, "no data → no evidence against equilibrium → PValue = 1 (§6.1)");
        r.InEquilibrium.Should().BeTrue("no data is treated as no evidence against the null model (§6.1)");
        double.IsNaN(r.ChiSquare).Should().BeFalse("the n=0 guard avoids a 0/0 NaN");
        double.IsNaN(r.PValue).Should().BeFalse("the n=0 guard avoids a 0/0 NaN");
    }

    #endregion

    #region BE — single genotype (the minimal non-empty sample)

    /// <summary>
    /// BE: the minimal non-empty sample — a single individual. Each is monomorphic
    /// for one genotype, so observed equals expected and χ² = 0, InEquilibrium =
    /// true (Hardy_Weinberg_Test.md §6.1 "Single monomorphic sample (1,0,0)").
    /// For the single homozygotes (1,0,0)/(0,0,1) the locus is fixed → two expected
    /// categories are exactly 0 and their χ² terms are skipped (the 0-expected
    /// boundary); for the single heterozygote (0,1,0) p = q = 0.5 but with a single
    /// observation the het-only sample is the all-het case at N=1 → χ² = 1, but
    /// χ²=1 on df=1 has a p-value ≈ 0.317 ≥ 0.05, so a single heterozygote is too
    /// small a sample to be a significant departure → InEquilibrium = true. We pin
    /// each defined result and confirm no DivideByZero and no NaN.
    /// </summary>
    [TestCase(1, 0, 0, 0.0, true, TestName = "HardyWeinberg_SingleHomozygousMajor_IsDefinedEquilibrium")]
    [TestCase(0, 0, 1, 0.0, true, TestName = "HardyWeinberg_SingleHomozygousMinor_IsDefinedEquilibrium")]
    [TestCase(0, 1, 0, 1.0, true, TestName = "HardyWeinberg_SingleHeterozygote_IsDefinedHetExcess")]
    public void HardyWeinberg_SingleGenotype_IsDefinedResult(
        int aa, int ab, int bb, double expectedChi, bool expectedInEq)
    {
        PopulationGeneticsAnalyzer.HardyWeinbergResult r = default;
        var act = () => r = PopulationGeneticsAnalyzer.TestHardyWeinberg("rs-single", aa, ab, bb);

        act.Should().NotThrow("a single genotype is a defined minimal sample, not an error");
        r.ChiSquare.Should().BeApproximately(expectedChi, Tolerance,
            because: "a single genotype gives an exact, documented χ² (§6.1)");
        r.InEquilibrium.Should().Be(expectedInEq, "the single-genotype verdict is defined and stable (§6.1)");

        r.ChiSquare.Should().BeGreaterThanOrEqualTo(0.0, "χ² ≥ 0 (INV-04)");
        double.IsNaN(r.ChiSquare).Should().BeFalse("a defined single-genotype sample never yields NaN");
        r.PValue.Should().BeInRange(0.0, 1.0, "the p-value is a probability (INV-05)");

        (r.ExpectedAA + r.ExpectedAa + r.Expectedaa).Should().BeApproximately(1.0, Tolerance,
            because: "the single expected counts sum to n = 1 (INV-03)");
    }

    #endregion

    #region BE — all heterozygous (the classic HWE-violation case: het excess → χ² = N)

    /// <summary>
    /// BE / KEY HWE-violation case: an all-heterozygous sample (0, N, 0) gives
    /// p = q = 0.5, so the expected counts are (0.25N, 0.5N, 0.25N) — but the
    /// observed are (0, N, 0), a pure het EXCESS. The chi-square then expands to
    /// 0.25N + 0.5N + 0.25N = N EXACTLY, a meaningful, nonzero, scaling departure
    /// from equilibrium (Hardy_Weinberg_Test.md §6.1; unit test HW-M9 pins N=100 →
    /// χ² = 100). Verified across a fixed-seed sweep of N so χ² = N holds at scale,
    /// the departure is always significant (InEquilibrium = false for N large
    /// enough to clear the 3.841 critical value), expected counts sum to n
    /// (INV-03), and χ² stays finite and ≥ 0 (INV-04) — never NaN.
    /// </summary>
    [Test]
    public void HardyWeinberg_AllHeterozygous_HasHetExcessChiSquareEqualToN()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int n = Rng.Next(4, 100_000); // ≥ 4 so χ² = N > 3.841 → a real, significant departure

            var r = PopulationGeneticsAnalyzer.TestHardyWeinberg("rs-allhet", 0, n, 0);

            r.ChiSquare.Should().BeApproximately(n, n * 1e-9 + Tolerance,
                because: $"all {n} heterozygous → p=q=0.5, het excess → χ² = N = {n} exactly (§6.1)");
            r.ChiSquare.Should().BeGreaterThan(3.841,
                "the het excess is a significant departure from HWE at α=0.05 (df=1 critical value 3.841)");
            r.InEquilibrium.Should().BeFalse(
                "a pure het-excess sample is rejected as out of equilibrium — the classic HWE violation (§6.1)");

            r.ExpectedAA.Should().BeApproximately(0.25 * n, Tolerance, "E_AA = p²n = 0.25N");
            r.ExpectedAa.Should().BeApproximately(0.5 * n, Tolerance, "E_Aa = 2pqn = 0.5N");
            r.Expectedaa.Should().BeApproximately(0.25 * n, Tolerance, "E_aa = q²n = 0.25N");
            (r.ExpectedAA + r.ExpectedAa + r.Expectedaa).Should().BeApproximately(n, Tolerance,
                because: "expected counts sum to n (INV-03)");

            r.ChiSquare.Should().BeGreaterThanOrEqualTo(0.0, "χ² ≥ 0 (INV-04)");
            double.IsNaN(r.ChiSquare).Should().BeFalse("a non-empty sample never yields a 0/0 NaN");
            double.IsInfinity(r.ChiSquare).Should().BeFalse("the bounded arithmetic never yields ±Infinity");
            r.PValue.Should().BeInRange(0.0, 1.0, "the p-value is a probability (INV-05)");
        }
    }

    #endregion

    #region BE — all homozygous (het deficit, and the monomorphic 0-expected boundary)

    /// <summary>
    /// BE "all homozygous", het-DEFICIT case: a sample split between the two
    /// homozygotes with NO heterozygotes (N, 0, N) gives p = q = 0.5 over a sample
    /// of 2N, so the expected counts are (0.5N, N, 0.5N) — but the observed are
    /// (N, 0, N), a pure het DEFICIT (the inbreeding pattern). The chi-square
    /// expands to 0.5N + N + 0.5N = 2N EXACTLY, a meaningful, nonzero departure
    /// from equilibrium. Verified across a fixed-seed sweep of N so χ² = 2N holds
    /// at scale, the departure is significant, expected counts sum to 2N (INV-03),
    /// and χ² stays finite and ≥ 0 — never NaN.
    /// </summary>
    [Test]
    public void HardyWeinberg_AllHomozygousMixed_HasHetDeficitChiSquareEqualToTwoN()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int n = Rng.Next(2, 50_000); // χ² = 2N ≥ 4 > 3.841 → a real, significant departure

            var r = PopulationGeneticsAnalyzer.TestHardyWeinberg("rs-homdef", n, 0, n);

            r.ChiSquare.Should().BeApproximately(2.0 * n, n * 1e-9 + Tolerance,
                because: $"(N,0,N) → p=q=0.5 over 2N, het deficit → χ² = 2N = {2 * n} exactly");
            r.ChiSquare.Should().BeGreaterThan(3.841,
                "the het deficit is a significant departure from HWE at α=0.05");
            r.InEquilibrium.Should().BeFalse(
                "a pure het-deficit (inbreeding) sample is rejected as out of equilibrium");

            (r.ExpectedAA + r.ExpectedAa + r.Expectedaa).Should().BeApproximately(2.0 * n, Tolerance,
                because: "expected counts sum to n = 2N (INV-03)");
            r.ChiSquare.Should().BeGreaterThanOrEqualTo(0.0, "χ² ≥ 0 (INV-04)");
            double.IsNaN(r.ChiSquare).Should().BeFalse("a non-empty sample never yields a 0/0 NaN");
            double.IsInfinity(r.ChiSquare).Should().BeFalse("the bounded arithmetic never yields ±Infinity");
            r.PValue.Should().BeInRange(0.0, 1.0, "the p-value is a probability (INV-05)");
        }
    }

    /// <summary>
    /// BE "all homozygous", MONOMORPHIC fully-fixed case — the KEY 0-expected
    /// boundary: when the whole sample is one homozygote (N, 0, 0) or (0, 0, N) the
    /// locus is FIXED (p = 1 or p = 0), so two of the three expected categories are
    /// EXACTLY 0. Their χ² terms (which would divide by the 0 expected count) are
    /// SKIPPED by the explicit `expected > 0` guard, so there is NO DivideByZero
    /// and NO NaN, and because observed equals expected the statistic is χ² = 0 →
    /// InEquilibrium = true (Hardy_Weinberg_Test.md §5.2, §6.1;
    /// PopulationGeneticsAnalyzer.cs lines 450–455). This pins the 0-expected
    /// division boundary the prompt flags as the chi-square danger spot.
    /// </summary>
    [Test]
    public void HardyWeinberg_MonomorphicFixed_SkipsZeroExpectedTermWithNoDivideByZero()
    {
        for (int trial = 0; trial < 16; trial++)
        {
            int n = Rng.Next(1, 100_000);

            foreach (var (aa, ab, bb, label) in new[]
                     {
                         (n, 0, 0, "fixed major"),
                         (0, 0, n, "fixed minor"),
                     })
            {
                PopulationGeneticsAnalyzer.HardyWeinbergResult r = default;
                var act = () => r = PopulationGeneticsAnalyzer.TestHardyWeinberg("rs-fixed", aa, ab, bb);

                act.Should().NotThrow<DivideByZeroException>(
                    $"the 0-expected categories of a {label} locus are skipped, not divided by");
                act.Should().NotThrow($"a {label} monomorphic sample of {n} is a defined input");

                r.ChiSquare.Should().BeApproximately(0.0, Tolerance,
                    because: $"a {label} locus has observed == expected → χ² = 0 (§6.1)");
                r.InEquilibrium.Should().BeTrue($"the {label} monomorphic counts match expectation exactly (§6.1)");

                double.IsNaN(r.ChiSquare).Should().BeFalse(
                    "the expected>0 guard skips the 0-expected term, so no 0/0 NaN arises");
                double.IsInfinity(r.ChiSquare).Should().BeFalse(
                    "skipping the 0-expected term keeps the statistic finite");
                r.PValue.Should().BeInRange(0.0, 1.0, "the p-value is a probability (INV-05)");
            }
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  POP-FST-001 — Wright's fixation index F_ST : fuzz targets
    //  Surface: CalculateFst(
    //      IEnumerable<(double AlleleFreq, int SampleSize)>,
    //      IEnumerable<(double AlleleFreq, int SampleSize)>)
    // ═══════════════════════════════════════════════════════════════════

    #region POP-FST-001 — Fst (fixation index)

    #region Helpers — population builders

    /// <summary>A per-locus (allele frequency, sample size) population vector.</summary>
    private static List<(double AlleleFreq, int SampleSize)> Pop(
        params (double AlleleFreq, int SampleSize)[] loci) => loci.ToList();

    #endregion

    #region Positive sanity — differentiated populations give the documented Fst in (0,1]; identical give 0

    /// <summary>
    /// Positive control: the worked examples pinned in the algorithm doc and the
    /// unit tests (F_Statistics.md §6.1; POP-FST-001 tests). A single differentiated
    /// locus pop1=(0.8, 100) vs pop2=(0.2, 100) gives p̄ = 0.5,
    /// σ²_S = (100·0.09 + 100·0.09)/200 = 0.09, H_т = 0.25, so
    /// F_ST = 0.09/0.25 = 0.36 — a value strictly inside (0, 1]. Equal-size
    /// fixed-opposite alleles (p1 = 1, p2 = 0) give F_ST = 1 (INV-FST-03), and two
    /// IDENTICAL populations give F_ST = 0 (INV-FST-02). The result must match the
    /// documented value, lie in [0, 1] (INV-FST-01), and be finite. This anchors the
    /// fuzz battery: before pinning that boundary input does no harm, confirm a
    /// known-good comparison gives the textbook-correct F_ST.
    /// </summary>
    [Test]
    public void Fst_DifferentiatedAndIdenticalPopulations_MatchDocumentedValues()
    {
        // Differentiated single locus → F_ST = 0.36, strictly inside (0, 1].
        double fstDiff = PopulationGeneticsAnalyzer.CalculateFst(
            Pop((0.8, 100)), Pop((0.2, 100)));

        fstDiff.Should().BeApproximately(0.36, Tolerance,
            because: "p̄=0.5, σ²_S=0.09, H_т=0.25 → F_ST = 0.09/0.25 = 0.36 (F_Statistics.md §6.1)");
        fstDiff.Should().BeGreaterThan(0.0, "two differentiated populations have F_ST > 0");
        fstDiff.Should().BeInRange(0.0, 1.0, "F_ST is a fixation index in [0, 1] (INV-FST-01)");
        double.IsNaN(fstDiff).Should().BeFalse("a known-good comparison never yields NaN");

        // Equal-size fixed-opposite alleles → complete differentiation F_ST = 1.
        double fstFixed = PopulationGeneticsAnalyzer.CalculateFst(
            Pop((1.0, 100)), Pop((0.0, 100)));
        fstFixed.Should().BeApproximately(1.0, Tolerance,
            because: "p̄=0.5, σ²_S=0.25, H_т=0.25 → F_ST = 1 (INV-FST-03, §6.1)");

        // Identical populations → no among-population variance → F_ST = 0.
        double fstSame = PopulationGeneticsAnalyzer.CalculateFst(
            Pop((0.5, 100), (0.3, 100)), Pop((0.5, 100), (0.3, 100)));
        fstSame.Should().Be(0.0,
            because: "identical allele frequencies at every locus → σ²_S = 0 → F_ST = 0 (INV-FST-02)");
    }

    #endregion

    #region BE — single population (one sequence supplied / 0 loci → nothing to differentiate)

    /// <summary>
    /// BE "single population": F_ST is a *between*-population statistic — with only
    /// one group there is nothing to differentiate. The surface models the single-
    /// population boundary as an EMPTY counterpart sequence: `CalculateFst(pop, [])`
    /// (or `[]` vs `pop`). The explicit `pop1.Count == 0 || pop2.Count == 0` guard
    /// returns the DEFINED 0 BEFORE any per-locus arithmetic — no DivideByZero, no
    /// NaN (F_Statistics.md §3.3, §6.1; PopulationGeneticsAnalyzer.cs lines
    /// 616–617). Verified across a fixed-seed sweep of single-population contents and
    /// both argument orders so the boundary is content-independent and symmetric.
    /// </summary>
    [Test]
    public void Fst_SinglePopulation_EmptyCounterpart_ReturnsDefinedZero()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int loci = Rng.Next(1, 12);
            var single = new List<(double, int)>();
            for (int i = 0; i < loci; i++)
                single.Add((Rng.NextDouble(), Rng.Next(1, 500)));

            var empty = new List<(double, int)>();

            double fstA = double.NaN, fstB = double.NaN;
            var actA = () => fstA = PopulationGeneticsAnalyzer.CalculateFst(single, empty);
            var actB = () => fstB = PopulationGeneticsAnalyzer.CalculateFst(empty, single);

            actA.Should().NotThrow<DivideByZeroException>(
                "the empty-counterpart boundary short-circuits before any per-locus division");
            actA.Should().NotThrow("a single population (empty counterpart) is a defined boundary, not an error");
            actB.Should().NotThrow("the empty-first-argument order is equally a defined boundary");

            fstA.Should().Be(0.0, $"only one population ({loci} loci) → nothing to differentiate → defined F_ST = 0 (§6.1)");
            fstB.Should().Be(0.0, "the guard is order-independent: empty-vs-population is also the defined 0");
            double.IsNaN(fstA).Should().BeFalse("the empty-sequence guard avoids any 0/0 NaN");
            double.IsNaN(fstB).Should().BeFalse("the empty-sequence guard avoids any 0/0 NaN");
        }
    }

    #endregion

    #region BE — identical populations (no among-population variance → F_ST = 0, the key identity)

    /// <summary>
    /// BE / KEY identity (INV-FST-02): when the two compared populations are
    /// IDENTICAL (same per-locus allele frequency and sample size) the weighted
    /// among-population variance σ²_S is 0 at every locus, so the numerator is 0 over
    /// a positive total-heterozygosity denominator → F_ST = 0 EXACTLY (no
    /// differentiation; F_Statistics.md §2.A INV-FST-02, §6.1). Verified across a
    /// fixed-seed sweep of random multi-locus populations so the identity holds
    /// independent of the (polymorphic) allele frequencies and locus counts. A
    /// polymorphic locus (0 &lt; p &lt; 1) is used so H_т &gt; 0 and a REAL division
    /// happens — proving F_ST = 0 by the zero variance, not by the H_т = 0 guard.
    /// </summary>
    [Test]
    public void Fst_IdenticalPopulations_IsExactlyZero()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int loci = Rng.Next(1, 15);
            var pop = new List<(double, int)>();
            for (int i = 0; i < loci; i++)
            {
                double p = 0.05 + 0.9 * Rng.NextDouble(); // strictly in (0,1) → H_т > 0
                pop.Add((p, Rng.Next(2, 1000)));
            }

            // A distinct list object with the SAME content, so the comparison is real.
            var copy = pop.ToList();

            double fst = PopulationGeneticsAnalyzer.CalculateFst(pop, copy);

            // σ²_S is the among-population variance of two *equal* frequencies, so it
            // is mathematically 0; with arbitrary (non-dyadic) random p the weighted
            // mean p̄ carries sub-ulp rounding, so (p−p̄)² is a floating-point
            // ~1e-30, not a hard 0. The identity F_ST = 0 is pinned to within the
            // shared 1e-9 tolerance — still an exact-identity assertion, not a
            // weakened one.
            fst.Should().BeApproximately(0.0, Tolerance,
                because: $"identical allele frequencies at every one of {loci} loci → σ²_S = 0 → F_ST = 0 (INV-FST-02)");
            fst.Should().BeInRange(0.0, 1.0, "F_ST stays in [0, 1] (INV-FST-01)");
            double.IsNaN(fst).Should().BeFalse("identical populations divide 0 by a positive H_т, never 0/0");
        }
    }

    /// <summary>
    /// BE: comparing a population against ITSELF (the same list object) is the
    /// canonical single-population/identical case — it must read as exactly 0
    /// differentiation, mirroring the diagonal of the pairwise matrix
    /// (`Fst(i, i) = 0`, F_Statistics.md §5.3.A).
    /// </summary>
    [Test]
    public void Fst_PopulationAgainstItself_IsExactlyZero()
    {
        var pop = Pop((0.7, 250), (0.4, 80), (0.9, 600));

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop, pop);

        fst.Should().Be(0.0, "a population is not differentiated from itself → F_ST = 0 (INV-FST-02)");
        double.IsNaN(fst).Should().BeFalse("self-comparison divides 0 by a positive denominator, not 0/0");
    }

    #endregion

    #region BE — 0 samples (n1 + n2 == 0 → the 0/0 NaN-poisoning boundary, guarded to a defined 0)

    /// <summary>
    /// BE / DivideByZero (here: 0/0 → NaN) boundary: a locus with sample sizes
    /// n1 = n2 = 0 makes the weighted mean p̄ = (0·p1 + 0·p2)/(0 + 0) a 0.0/0 → NaN.
    /// Because the arithmetic is `double`, this is a NaN path, NOT a
    /// DivideByZeroException. The NaN poisons the per-locus H_т and hence the
    /// accumulated denominator; the final `denominator > 0` test is then FALSE
    /// (NaN &gt; 0 is false), so the method returns the DEFINED 0 and NO NaN escapes
    /// to the caller (PopulationGeneticsAnalyzer.cs line 645). This is the single
    /// most important fuzz concern for this unit — an unguarded weighted mean here
    /// would leak a NaN F_ST into production. Both a single all-zero-sample locus and
    /// the pop-vs-empty 0-loci form are pinned.
    /// </summary>
    [Test]
    public void Fst_ZeroSampleSizes_ReturnsDefinedZeroWithNoCrashAndNoNaN()
    {
        // A single locus with zero sample sizes in both populations.
        double fst = double.NaN;
        var act = () => fst = PopulationGeneticsAnalyzer.CalculateFst(
            Pop((0.5, 0)), Pop((0.5, 0)));

        act.Should().NotThrow<DivideByZeroException>(
            "the n1+n2==0 weighted mean is double arithmetic (0/0 → NaN), never a DivideByZeroException");
        act.Should().NotThrow("a 0-sample locus is a defined degenerate boundary, not an error");
        fst.Should().Be(0.0,
            because: "the NaN-poisoned denominator fails the `> 0` guard → defined F_ST = 0 (line 645)");
        double.IsNaN(fst).Should().BeFalse("the denominator guard ensures NO NaN escapes to the caller");
        double.IsInfinity(fst).Should().BeFalse("no ±Infinity escapes the guarded ratio");

        // Different allele frequencies but still zero samples → same guarded 0.
        double fstDiff = double.NaN;
        var actDiff = () => fstDiff = PopulationGeneticsAnalyzer.CalculateFst(
            Pop((1.0, 0)), Pop((0.0, 0)));
        actDiff.Should().NotThrow("zero samples with opposite fixed alleles is still a defined boundary");
        fstDiff.Should().Be(0.0, "the zero-sample 0/0 path returns the guarded 0 regardless of allele frequencies");
        double.IsNaN(fstDiff).Should().BeFalse("no NaN escapes even when p1 ≠ p2 at the zero-sample locus");
    }

    /// <summary>
    /// BE / H_т = 0 (DivideByZero) boundary — the monomorphic-locus case the prompt
    /// flags: when EVERY compared locus is fixed for the SAME allele in both
    /// populations (p1 = p2 ∈ {0, 1}) the total heterozygosity H_т = p̄(1−p̄) = 0 at
    /// every locus, so the accumulated denominator is exactly 0. The explicit
    /// `denominator > 0 ? … : 0` guard returns the DEFINED 0 instead of dividing by
    /// zero (F_Statistics.md §6.1 "both fixed for the same allele … returns 0";
    /// PopulationGeneticsAnalyzer.cs line 645). Verified across a fixed-seed sweep of
    /// monomorphic populations (positive sample sizes, so this is the H_т=0 boundary,
    /// distinct from the n=0 boundary above) so no DivideByZero and no NaN ever
    /// escape, regardless of which allele is fixed or how many loci.
    /// </summary>
    [Test]
    public void Fst_MonomorphicSameAllele_ReturnsDefinedZeroWithNoDivideByZero()
    {
        for (int trial = 0; trial < 16; trial++)
        {
            int loci = Rng.Next(1, 10);
            double fixedAllele = Rng.Next(2) == 0 ? 0.0 : 1.0; // both fixed for 0 or for 1
            var pop1 = new List<(double, int)>();
            var pop2 = new List<(double, int)>();
            for (int i = 0; i < loci; i++)
            {
                pop1.Add((fixedAllele, Rng.Next(1, 1000)));
                pop2.Add((fixedAllele, Rng.Next(1, 1000)));
            }

            double fst = double.NaN;
            var act = () => fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

            act.Should().NotThrow<DivideByZeroException>(
                "the H_т=0 denominator is guarded by `denominator > 0 ? … : 0`, never divided by");
            act.Should().NotThrow("a monomorphic same-allele comparison is a defined boundary, not an error");
            fst.Should().Be(0.0,
                because: $"all {loci} loci fixed for allele {fixedAllele} → H_т = 0 → guarded defined F_ST = 0 (§6.1)");
            double.IsNaN(fst).Should().BeFalse("the H_т=0 guard avoids a 0/0 NaN");
            double.IsInfinity(fst).Should().BeFalse("the H_т=0 guard avoids a ±Infinity");
        }
    }

    #endregion

    #region BE — mismatched locus counts (documented ArgumentException, never a silent miscompare)

    /// <summary>
    /// BE: the two per-locus sequences must align 1-to-1. Two NON-EMPTY sequences of
    /// DIFFERENT length are rejected with the *documented, intentional*
    /// ArgumentException (PopulationGeneticsAnalyzer.cs lines 619–623) — NOT a silent
    /// prefix-truncated miscompare and NOT a crash. (The empty-vs-non-empty case is
    /// the single-population boundary handled by the earlier guard and returns 0; the
    /// length-mismatch gate applies only when BOTH sequences are non-empty.) This
    /// pins the validation boundary so the strict count contract cannot drift.
    /// </summary>
    [Test]
    public void Fst_MismatchedNonEmptyLocusCounts_ThrowsDocumentedArgumentException()
    {
        var pop1 = Pop((0.5, 100), (0.3, 100));
        var pop2 = Pop((0.5, 100)); // one locus short, but non-empty

        var act = () => PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        act.Should().Throw<ArgumentException>(
            "two non-empty populations with mismatched per-locus counts are rejected, " +
            "never silently compared on a truncated prefix (F_Statistics.md §3.3; lines 619–623)");
    }

    #endregion

    #region BE — random differentiated populations stay finite and in [0, 1]

    /// <summary>
    /// BE / INV-FST-01 sweep: across a fixed-seed battery of random, equal-length
    /// two-population comparisons (polymorphic allele frequencies, positive sample
    /// sizes) F_ST must ALWAYS be a well-defined, finite value in [0, 1] — a
    /// non-negative weighted variance over a positive total heterozygosity can never
    /// be negative, exceed 1, NaN, or ±Infinity. This pins the total-function
    /// contract over the random-but-valid interior, complementing the degenerate-
    /// boundary tests.
    /// </summary>
    [Test]
    public void Fst_RandomDifferentiatedPopulations_AlwaysInUnitIntervalAndFinite()
    {
        for (int trial = 0; trial < 64; trial++)
        {
            int loci = Rng.Next(1, 20);
            var pop1 = new List<(double, int)>();
            var pop2 = new List<(double, int)>();
            for (int i = 0; i < loci; i++)
            {
                pop1.Add((0.01 + 0.98 * Rng.NextDouble(), Rng.Next(1, 2000)));
                pop2.Add((0.01 + 0.98 * Rng.NextDouble(), Rng.Next(1, 2000)));
            }

            double fst = double.NaN;
            var act = () => fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

            act.Should().NotThrow("a random equal-length valid two-population input is always defined");
            double.IsNaN(fst).Should().BeFalse("a positive H_т denominator never yields a 0/0 NaN");
            double.IsInfinity(fst).Should().BeFalse("a bounded variance over a positive denominator never gives ±Infinity");
            fst.Should().BeGreaterThanOrEqualTo(0.0, "F_ST is a ratio of non-negative quantities (INV-FST-01)");
            fst.Should().BeLessThanOrEqualTo(1.0, "F_ST never exceeds complete differentiation (INV-FST-01)");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  POP-LD-001 — linkage disequilibrium : fuzz targets
    //  Surface: CalculateLD(string, string, IEnumerable<(int Geno1, int Geno2)>, int)
    // ═══════════════════════════════════════════════════════════════════

    #region POP-LD-001 — linkage disequilibrium

    #region Helpers — genotype-pair builder

    /// <summary>Zips two equal-length 0/1/2 genotype vectors into the paired sequence CalculateLD expects.</summary>
    private static List<(int Geno1, int Geno2)> Pairs(int[] geno1, int[] geno2) =>
        geno1.Zip(geno2, (g1, g2) => (g1, g2)).ToList();

    /// <summary>A deterministic random 0/1/2 diploid genotype vector of the given length.</summary>
    private static int[] RandomGenotypes(int length)
    {
        var buffer = new int[length];
        for (int i = 0; i < length; i++)
            buffer[i] = Rng.Next(0, 3); // 0 = hom major, 1 = het, 2 = hom minor
        return buffer;
    }

    #endregion

    #region Positive sanity — known two-locus genotype set yields the documented r² and D' in range

    /// <summary>
    /// Positive control: a known two-locus genotype distribution with a clean,
    /// exactly-pinnable LD. For the six paired 0/1/2 genotypes
    ///   (0,0) (0,1) (1,1) (1,2) (2,2) (2,0)
    /// the squared Pearson correlation of the two genotype vectors is r² = 1/16 =
    /// 0.0625 EXACTLY and Lewontin's normalized D' = 1/3 (computed from the diploid
    /// covariance D = Cov/2 over D_max per Linkage_Disequilibrium.md §2.2/§5.2).
    /// Both summaries must hit those documented values, r² must lie in [0, 1]
    /// (INV-02) and D' in [0, 1] (INV-03), and neither may be NaN/±Infinity. A
    /// second anchor pins PERFECT LD (identical genotype vectors → r² = 1, D' = 1,
    /// §6.1). This anchors the fuzz battery: before pinning that boundary input does
    /// no harm, confirm a known-good two-locus sample gives the textbook-correct LD.
    /// </summary>
    [Test]
    public void LinkageDisequilibrium_KnownTwoLocusSet_MatchesDocumentedRSquaredAndDPrimeInRange()
    {
        var genotypes = Pairs(
            new[] { 0, 0, 1, 1, 2, 2 },
            new[] { 0, 1, 1, 2, 2, 0 });

        var ld = PopulationGeneticsAnalyzer.CalculateLD("rsA", "rsB", genotypes, distance: 1500);

        ld.RSquared.Should().BeApproximately(0.0625, Tolerance,
            because: "the squared Pearson correlation of the two 0/1/2 vectors is r² = 1/16 = 0.0625 " +
                     "(Linkage_Disequilibrium.md §2.2 r² = D²/(p_A q_A p_B q_B); §5.2 genotype-correlation surface)");
        ld.DPrime.Should().BeApproximately(1.0 / 3.0, Tolerance,
            because: "D = Cov/2 normalized by D_max gives Lewontin D' = 1/3 (Linkage_Disequilibrium.md §2.2, §5.2)");

        ld.RSquared.Should().BeInRange(0.0, 1.0, "r² is a squared correlation measure (INV-02)");
        ld.DPrime.Should().BeInRange(0.0, 1.0, "the public DPrime is the non-negative normalized |D'| (INV-03)");
        double.IsNaN(ld.RSquared).Should().BeFalse("a known-good sample never yields NaN");
        double.IsNaN(ld.DPrime).Should().BeFalse("a known-good sample never yields NaN");

        ld.Variant1.Should().Be("rsA", "the first variant id is preserved verbatim");
        ld.Variant2.Should().Be("rsB", "the second variant id is preserved verbatim");
        ld.Distance.Should().Be(1500, "the input distance is copied into the record");

        // Second independent anchor: identical genotype vectors → perfect LD.
        var perfect = Pairs(
            new[] { 0, 1, 2, 0, 1, 2 },
            new[] { 0, 1, 2, 0, 1, 2 });
        var perfectLd = PopulationGeneticsAnalyzer.CalculateLD("rsC", "rsD", perfect, distance: 50);
        perfectLd.RSquared.Should().BeApproximately(1.0, Tolerance,
            because: "identical 0/1/2 genotype vectors have Pearson correlation 1 → r² = 1 (§6.1 perfect LD)");
        perfectLd.DPrime.Should().BeApproximately(1.0, Tolerance,
            because: "perfect association → D' = 1 (§6.1 perfect LD)");
    }

    #endregion

    #region BE — 0 samples (the empty genotype-pair sequence: every mean/var would be 0/0)

    /// <summary>
    /// BE / DivideByZero boundary: with 0 samples there are no genotype pairs, so
    /// every mean, covariance and variance would be a 0/0. The contract must NOT
    /// divide by zero and must NOT return NaN — the explicit `genoList.Count == 0`
    /// guard returns the DEFINED LinkageDisequilibrium with DPrime = 0, RSquared = 0
    /// (IDs and distance preserved) BEFORE any arithmetic
    /// (Linkage_Disequilibrium.md §3.3, §6.1; PopulationGeneticsAnalyzer.cs lines
    /// 737–740). This is the single most important fuzz concern for this unit: an
    /// unguarded mean/covariance over 0 samples is a 0/0 NaN in production.
    /// </summary>
    [Test]
    public void LinkageDisequilibrium_ZeroSamples_ReturnsDefinedZeroWithNoDivideByZero()
    {
        PopulationGeneticsAnalyzer.LinkageDisequilibrium ld = default;
        var act = () => ld = PopulationGeneticsAnalyzer.CalculateLD(
            "rsEmpty1", "rsEmpty2", Array.Empty<(int, int)>(), distance: 9000);

        act.Should().NotThrow<DivideByZeroException>(
            "the empty-sample boundary short-circuits before any mean/covariance/variance is computed");
        act.Should().NotThrow("0 samples is a defined boundary input (returns zero LD), not an error");

        ld.RSquared.Should().Be(0.0, "no samples → no association measurable → r² is the defined 0 (§6.1)");
        ld.DPrime.Should().Be(0.0, "no samples → no association measurable → D' is the defined 0 (§6.1)");
        double.IsNaN(ld.RSquared).Should().BeFalse("the empty-input guard avoids a 0/0 NaN");
        double.IsNaN(ld.DPrime).Should().BeFalse("the empty-input guard avoids a 0/0 NaN");

        ld.Variant1.Should().Be("rsEmpty1", "the variant ids are preserved even on the empty boundary");
        ld.Variant2.Should().Be("rsEmpty2", "the variant ids are preserved even on the empty boundary");
        ld.Distance.Should().Be(9000, "the distance is preserved even on the empty boundary");
    }

    #endregion

    #region BE — single (monomorphic) locus (the r² zero-denominator: variance 0 → p_A q_A = 0)

    /// <summary>
    /// BE / DivideByZero boundary (KEY): when ONE of the two loci is monomorphic
    /// (its genotype vector is invariant) that locus's variance is 0, so the r²
    /// denominator p_A·q_A·p_B·q_B contains a 0 factor. An unguarded r² = D²/0 is a
    /// DivideByZero / NaN. The contract returns RSquared = 0 via the explicit
    /// `var1 > 0 && var2 > 0` guard (Linkage_Disequilibrium.md §3.3, §6.1;
    /// PopulationGeneticsAnalyzer.cs line 760) — NO DivideByZero, NO NaN, NO
    /// ±Infinity. Verified for the FIRST locus monomorphic and the SECOND locus
    /// monomorphic (the second factor of the denominator), each across a fixed-seed
    /// sweep of the OTHER, polymorphic locus so the guard holds regardless of the
    /// well-behaved partner.
    /// </summary>
    [Test]
    public void LinkageDisequilibrium_OneMonomorphicLocus_ReturnsZeroRSquaredWithNoDivideByZero()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int n = Rng.Next(2, 30);
            int[] polymorphic = RandomGenotypes(n);
            // Force the partner to be genuinely polymorphic so only the fixed locus is degenerate.
            polymorphic[0] = 0;
            polymorphic[1] = 2;
            int fixedAllele = Rng.Next(0, 3);
            var monomorphic = Enumerable.Repeat(fixedAllele, n).ToArray();

            foreach (var genotypes in new[]
                     {
                         Pairs(monomorphic, polymorphic), // first locus fixed → var1 = 0
                         Pairs(polymorphic, monomorphic), // second locus fixed → var2 = 0
                     })
            {
                PopulationGeneticsAnalyzer.LinkageDisequilibrium ld = default;
                var act = () => ld = PopulationGeneticsAnalyzer.CalculateLD("rs1", "rs2", genotypes, 100);

                act.Should().NotThrow<DivideByZeroException>(
                    "a monomorphic locus zeroes a denominator factor; the var>0 guard prevents the division");
                act.Should().NotThrow("a monomorphic locus is a defined boundary input, not an error");

                ld.RSquared.Should().Be(0.0,
                    because: $"a fixed locus (allele {fixedAllele}) makes the r² denominator 0 → guarded RSquared = 0 (§6.1)");
                double.IsNaN(ld.RSquared).Should().BeFalse("the var>0 guard avoids a D²/0 NaN");
                double.IsInfinity(ld.RSquared).Should().BeFalse("the var>0 guard avoids a D²/0 ±Infinity");
                ld.DPrime.Should().BeInRange(0.0, 1.0, "D' stays in [0,1] even when the locus is fixed (INV-03)");
                double.IsNaN(ld.DPrime).Should().BeFalse("a monomorphic locus must not poison D' into NaN");
            }
        }
    }

    #endregion

    #region BE — monomorphic loci / single sample (BOTH variances 0 → covariance 0 too)

    /// <summary>
    /// BE: BOTH loci monomorphic (each genotype vector fixed for one allele) — both
    /// variances are 0 AND the covariance is 0, so BOTH factors of the r²
    /// denominator vanish. The `var1 > 0 && var2 > 0` guard still returns RSquared =
    /// 0 with no division, and D' stays a finite value in [0, 1]. Verified across a
    /// fixed-seed sweep of the two fixed alleles and the sample size. Also pins the
    /// n = 1 single-genotype-pair boundary, where a one-point sample has zero spread
    /// at both loci → the same guarded, finite, defined result (§6.1; M8/S1).
    /// </summary>
    [Test]
    public void LinkageDisequilibrium_BothMonomorphicOrSinglePair_ReturnsDefinedFiniteZeroRSquared()
    {
        // Single genotype pair (n = 1): zero spread at both loci.
        {
            PopulationGeneticsAnalyzer.LinkageDisequilibrium ld = default;
            var act = () => ld = PopulationGeneticsAnalyzer.CalculateLD(
                "rsSingle1", "rsSingle2", new List<(int, int)> { (1, 1) }, 100);

            act.Should().NotThrow("a single genotype pair is a defined boundary input, not an error");
            ld.RSquared.Should().Be(0.0, "a one-point sample has zero variance at both loci → guarded r² = 0");
            double.IsNaN(ld.RSquared).Should().BeFalse("the single-pair boundary must not yield NaN");
            double.IsNaN(ld.DPrime).Should().BeFalse("the single-pair boundary must not yield NaN");
            ld.DPrime.Should().BeInRange(0.0, 1.0, "D' stays in [0,1] for the single-pair boundary (INV-03)");
        }

        // Both loci monomorphic across a sweep of fixed-allele choices and sizes.
        for (int trial = 0; trial < 24; trial++)
        {
            int n = Rng.Next(2, 40);
            int a = Rng.Next(0, 3);
            int b = Rng.Next(0, 3);
            var genotypes = Pairs(
                Enumerable.Repeat(a, n).ToArray(),
                Enumerable.Repeat(b, n).ToArray());

            PopulationGeneticsAnalyzer.LinkageDisequilibrium ld = default;
            var act = () => ld = PopulationGeneticsAnalyzer.CalculateLD("rs1", "rs2", genotypes, 100);

            act.Should().NotThrow<DivideByZeroException>(
                "both loci fixed → both denominator factors 0; the var>0 guard prevents the division");
            act.Should().NotThrow("two monomorphic loci is a defined boundary input, not an error");

            ld.RSquared.Should().Be(0.0,
                because: $"both loci fixed (alleles {a}, {b}) → covariance and both variances 0 → guarded r² = 0 (§6.1)");
            double.IsNaN(ld.RSquared).Should().BeFalse("the var>0 guard avoids a 0/0 NaN");
            double.IsInfinity(ld.RSquared).Should().BeFalse("the var>0 guard avoids a 0/0 ±Infinity");
            ld.DPrime.Should().BeInRange(0.0, 1.0, "D' stays in [0,1] for two fixed loci (INV-03)");
            double.IsNaN(ld.DPrime).Should().BeFalse("two fixed loci must not poison D' into NaN");
        }
    }

    #endregion

    #region BE — random two-locus samples always have r² ∈ [0,1] and D' ∈ [0,1] and are finite

    /// <summary>
    /// BE / INV-02+INV-03 sweep: across a fixed-seed battery of random equal-length
    /// 0/1/2 genotype-pair samples the public RSquared must ALWAYS lie in [0, 1]
    /// (INV-02 — it is a squared correlation) and DPrime in [0, 1] (INV-03 — the
    /// implementation returns |D'| clamped to 1), and both must be finite: a squared
    /// correlation over the guarded p_A·q_A·p_B·q_B denominator can never be
    /// negative, exceed 1, NaN, or ±Infinity. This pins the total-function contract
    /// over the random interior, complementing the degenerate-boundary tests.
    /// </summary>
    [Test]
    public void LinkageDisequilibrium_RandomTwoLocusSamples_AlwaysInUnitIntervalAndFinite()
    {
        for (int trial = 0; trial < 64; trial++)
        {
            int n = Rng.Next(1, 60);
            var genotypes = Pairs(RandomGenotypes(n), RandomGenotypes(n));

            PopulationGeneticsAnalyzer.LinkageDisequilibrium ld = default;
            var act = () => ld = PopulationGeneticsAnalyzer.CalculateLD("rsX", "rsY", genotypes, n);

            act.Should().NotThrow("a random equal-length 0/1/2 genotype-pair sample is always a defined input");
            double.IsNaN(ld.RSquared).Should().BeFalse("the guarded denominator never yields a 0/0 NaN");
            double.IsInfinity(ld.RSquared).Should().BeFalse("a bounded squared correlation never gives ±Infinity");
            double.IsNaN(ld.DPrime).Should().BeFalse("the guarded D_max normalization never yields a 0/0 NaN");
            double.IsInfinity(ld.DPrime).Should().BeFalse("the clamped |D'| never gives ±Infinity");
            ld.RSquared.Should().BeInRange(0.0, 1.0, "r² is a squared correlation measure (INV-02)");
            ld.DPrime.Should().BeInRange(0.0, 1.0, "the public DPrime is |D'| clamped to [0,1] (INV-03)");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  POP-ANCESTRY-001 — supervised / projection ADMIXTURE : fuzz targets
    //  Surface: EstimateAncestry(individuals, referencePops, maxIterations)
    // ═══════════════════════════════════════════════════════════════════

    #region POP-ANCESTRY-001 — ancestry estimation

    #region Helpers — individual / reference-panel builders

    /// <summary>Materializes one individual (id + 0/1/2 genotype vector) for the ancestry API.</summary>
    private static (string, IReadOnlyList<int>) Individual(string id, params int[] genotypes) =>
        (id, (IReadOnlyList<int>)genotypes);

    /// <summary>Materializes one reference panel (id + per-SNP allele-1 frequencies) for the API.</summary>
    private static (string, IReadOnlyList<double>) Panel(string id, params double[] frequencies) =>
        (id, (IReadOnlyList<double>)frequencies);

    /// <summary>Sum of a single individual's ancestry fractions — the INV-01 quantity.</summary>
    private static double FractionSum(PopulationGeneticsAnalyzer.AncestryProportion ap) =>
        ap.Proportions.Values.Sum();

    #endregion

    #region Positive sanity — the documented worked example yields the exact Eq. 4 fractions

    /// <summary>
    /// Positive control: the worked example pinned in the algorithm doc
    /// (Ancestry_Estimation.md §7.1) — start q = (0.5, 0.5), reference panels
    /// f_A = (0.8, 0.2) and f_B = (0.2, 0.8), individual genotype g = (2, 0), one
    /// EM iteration (Eq. 4). The hand-checked walk-through gives
    /// q_A = (1.6 + 1.6)/(2·2) = 0.8 and q_B = (0.4 + 0.4)/4 = 0.2. The result must
    /// be EXACTLY those fractions, each in [0, 1] (INV-02), summing to 1 (INV-01).
    /// This anchors the fuzz battery: before pinning that boundary input does no
    /// harm, confirm a known-good input gives the textbook-derived ancestry vector.
    /// Expected values come from the DOC's arithmetic, not the code.
    /// </summary>
    [Test]
    public void Ancestry_DocumentedWorkedExample_MatchesEq4FractionsAndSumsToOne()
    {
        var individuals = new[] { Individual("ind1", 2, 0) };
        var refs = new[]
        {
            Panel("A", 0.8, 0.2),
            Panel("B", 0.2, 0.8),
        };

        var result = PopulationGeneticsAnalyzer.EstimateAncestry(individuals, refs, maxIterations: 1).Single();

        result.IndividualId.Should().Be("ind1", "the input id is echoed unchanged (§3.2)");
        result.Proportions["A"].Should().BeApproximately(0.8, Tolerance,
            because: "q_A = (1.6 + 1.6)/(2·2) = 0.8 after one EM iteration (Ancestry_Estimation.md §7.1, Eq. 4)");
        result.Proportions["B"].Should().BeApproximately(0.2, Tolerance,
            because: "q_B = (0.4 + 0.4)/(2·2) = 0.2 after one EM iteration (Ancestry_Estimation.md §7.1, Eq. 4)");

        result.Proportions["A"].Should().BeInRange(0.0, 1.0, "an ancestry fraction is a proportion (INV-02)");
        result.Proportions["B"].Should().BeInRange(0.0, 1.0, "an ancestry fraction is a proportion (INV-02)");
        FractionSum(result).Should().BeApproximately(1.0, Tolerance,
            because: "the EM divides by 2J so the ancestry fractions sum to 1 (INV-01)");
    }

    #endregion

    #region BE — single source population (pure member → fraction 1 for one panel, 0 for the rest)

    /// <summary>
    /// BE "single source population": an individual that is a PURE member of one
    /// reference population — homozygous (g = 2 at every SNP) for the allele that
    /// panel A is fixed for (f_A = 1) while panel B is fixed for the other allele
    /// (f_B = 0). The single-source corner of the simplex: A's ancestry fraction
    /// must converge to exactly 1 and B's to exactly 0 (Ancestry_Estimation.md §7.1
    /// diagnostic-individual semantics; INV-02 corner). Hand-checked: from
    /// q = (0.5, 0.5), each SNP (g = 2) contributes A = 2·(0.5·1/0.5) = 2, B = 0,
    /// so q_A = (2+2)/(2·2) = 1.0, q_B = 0.0 after one iteration. The fractions
    /// must stay in [0, 1] (INV-02) and sum to 1 (INV-01).
    /// </summary>
    [Test]
    public void Ancestry_PureSourcePopulation_FractionIsOneForThatPanelZeroForOthers()
    {
        var individuals = new[] { Individual("pureA", 2, 2, 2) };
        var refs = new[]
        {
            Panel("A", 1.0, 1.0, 1.0),  // A fixed for allele 1
            Panel("B", 0.0, 0.0, 0.0),  // B fixed for allele 2
        };

        var result = PopulationGeneticsAnalyzer.EstimateAncestry(individuals, refs, maxIterations: 100).Single();

        result.Proportions["A"].Should().BeApproximately(1.0, Tolerance,
            because: "an individual fixed for the allele panel A is fixed for is pure A → fraction 1 (§7.1, single source)");
        result.Proportions["B"].Should().BeApproximately(0.0, Tolerance,
            because: "no genome derives from panel B → fraction 0 for a pure-A individual");

        result.Proportions["A"].Should().BeInRange(0.0, 1.0, "INV-02");
        result.Proportions["B"].Should().BeInRange(0.0, 1.0, "INV-02");
        FractionSum(result).Should().BeApproximately(1.0, Tolerance, "ancestry fractions sum to 1 (INV-01)");

        // Symmetric pure-B individual: homozygous allele 2 (g = 0) → fraction 1 for B, 0 for A.
        var pureB = PopulationGeneticsAnalyzer.EstimateAncestry(
            new[] { Individual("pureB", 0, 0, 0) }, refs, maxIterations: 100).Single();
        pureB.Proportions["B"].Should().BeApproximately(1.0, Tolerance,
            because: "an individual fixed for the allele panel B is fixed for is pure B → fraction 1");
        pureB.Proportions["A"].Should().BeApproximately(0.0, Tolerance, "no genome derives from panel A");
        FractionSum(pureB).Should().BeApproximately(1.0, Tolerance, "INV-01");
    }

    #endregion

    #region BE — admixed 50/50 (equal contribution → equal fractions, the K=2 uniform fixed point)

    /// <summary>
    /// BE "admixed 50/50": with two IDENTICAL reference panels the uniform initial
    /// vector q = (0.5, 0.5) is a FIXED POINT of the EM (Ancestry_Estimation.md
    /// §2.4 INV-04 — "identical reference panels keep a uniform q uniform"): the
    /// Eq. 4 numerators and denominators are proportional across the two panels, so
    /// the update reproduces (0.5, 0.5) exactly. With nothing to distinguish the
    /// sources the two contribute equally → an admixed-50/50 ancestry vector. The
    /// fractions must be 0.5/0.5 regardless of the genotype, each in [0, 1]
    /// (INV-02), summing to 1 (INV-01). Verified across a fixed-seed sweep of
    /// random 0/1/2 genotypes so the 50/50 fixed point is genotype-independent.
    /// </summary>
    [Test]
    public void Ancestry_IdenticalPanels_StayAdmixedFiftyFifty()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int snps = Rng.Next(1, 24);
            var freqs = new double[snps];
            for (int j = 0; j < snps; j++)
                freqs[j] = Rng.NextDouble();           // any valid in-[0,1] panel
            var genotype = new int[snps];
            for (int j = 0; j < snps; j++)
                genotype[j] = Rng.Next(0, 3);          // valid 0/1/2 genotype

            var individuals = new[] { ("ind", (IReadOnlyList<int>)genotype) };
            var refs = new[]
            {
                ("A", (IReadOnlyList<double>)freqs.ToArray()),
                ("B", (IReadOnlyList<double>)freqs.ToArray()),  // identical panel
            };

            var result = PopulationGeneticsAnalyzer.EstimateAncestry(individuals, refs, maxIterations: 100).Single();

            result.Proportions["A"].Should().BeApproximately(0.5, Tolerance,
                because: "identical panels keep the uniform q uniform → 50/50 admixture (INV-04)");
            result.Proportions["B"].Should().BeApproximately(0.5, Tolerance,
                because: "identical panels keep the uniform q uniform → 50/50 admixture (INV-04)");
            FractionSum(result).Should().BeApproximately(1.0, Tolerance, "INV-01");
        }
    }

    /// <summary>
    /// BE "admixed 50/50", second anchor: a symmetric individual under symmetric,
    /// MIRROR-IMAGE panels also reads as 50/50. With f_A = (0.8, 0.2),
    /// f_B = (0.2, 0.8) and a heterozygote at both SNPs (g = (1, 1)), the Eq. 4
    /// responsibilities are symmetric between A and B at each SNP, so the uniform
    /// q = (0.5, 0.5) is reproduced. Hand-checked: SNP1 g=1, mix1 = 0.5, mix2 = 0.5,
    /// A += 1·(0.5·0.8/0.5) + 1·(0.5·0.2/0.5) = 0.8 + 0.2 = 1.0, B += 0.2 + 0.8 = 1.0;
    /// likewise SNP2; q_A = (1.0+1.0)/(2·2) = 0.5, q_B = 0.5. The mirror symmetry
    /// gives a genuine equal split even with DISTINCT panels.
    /// </summary>
    [Test]
    public void Ancestry_SymmetricHeterozygoteUnderMirrorPanels_IsFiftyFifty()
    {
        var individuals = new[] { Individual("symHet", 1, 1) };
        var refs = new[]
        {
            Panel("A", 0.8, 0.2),
            Panel("B", 0.2, 0.8),
        };

        var result = PopulationGeneticsAnalyzer.EstimateAncestry(individuals, refs, maxIterations: 100).Single();

        result.Proportions["A"].Should().BeApproximately(0.5, Tolerance,
            because: "a heterozygote under mirror-image panels is symmetric between A and B → 50/50 (Eq. 4)");
        result.Proportions["B"].Should().BeApproximately(0.5, Tolerance,
            because: "a heterozygote under mirror-image panels is symmetric between A and B → 50/50 (Eq. 4)");
        FractionSum(result).Should().BeApproximately(1.0, Tolerance, "INV-01");
    }

    #endregion

    #region BE — empty input (no individuals or no panels → empty result, the "nothing to estimate" boundary)

    /// <summary>
    /// BE "empty": with EMPTY individuals OR empty reference panels there is nothing
    /// to estimate. The explicit `indList.Count == 0 || refList.Count == 0` guard
    /// `yield break`s an EMPTY result BEFORE any EM iteration
    /// (Ancestry_Estimation.md §3.3, §6.1; PopulationGeneticsAnalyzer.cs lines
    /// 1272–1273) — there is NO 1/2J division and NO exception. This is the KEY
    /// empty-input fuzz concern: an unguarded `refList[0]` or `1/k` here would crash
    /// (IndexOutOfRange / DivideByZero) in production. Covered: empty individuals
    /// with non-empty panels, non-empty individuals with empty panels, and both
    /// empty. Enumeration is forced with `.ToList()` so the deferred iterator runs.
    /// </summary>
    [Test]
    public void Ancestry_EmptyInput_ReturnsEmptyResultWithNoExceptionOrDivision()
    {
        var emptyInds = Array.Empty<(string, IReadOnlyList<int>)>();
        var emptyPanels = Array.Empty<(string, IReadOnlyList<double>)>();
        var someInds = new[] { Individual("ind", 2, 0) };
        var somePanels = new[] { Panel("A", 0.8, 0.2), Panel("B", 0.2, 0.8) };

        List<PopulationGeneticsAnalyzer.AncestryProportion> r1 = null!, r2 = null!, r3 = null!;

        var actEmptyInds = () => r1 = PopulationGeneticsAnalyzer.EstimateAncestry(emptyInds, somePanels).ToList();
        var actEmptyPanels = () => r2 = PopulationGeneticsAnalyzer.EstimateAncestry(someInds, emptyPanels).ToList();
        var actBothEmpty = () => r3 = PopulationGeneticsAnalyzer.EstimateAncestry(emptyInds, emptyPanels).ToList();

        actEmptyInds.Should().NotThrow("empty individuals is a defined boundary (nothing to estimate), not an error");
        actEmptyPanels.Should().NotThrow<DivideByZeroException>(
            "the empty-panels guard short-circuits before the 1/K uniform-init division");
        actEmptyPanels.Should().NotThrow("empty reference panels is a defined boundary, not an error");
        actBothEmpty.Should().NotThrow("both empty is a defined boundary, not an error");

        r1.Should().BeEmpty("no individuals → nothing to estimate → empty result (§6.1)");
        r2.Should().BeEmpty("no reference panels → nothing to estimate → empty result (§6.1)");
        r3.Should().BeEmpty("both empty → empty result (§6.1)");
    }

    #endregion

    #region BE — random valid inputs always obey INV-01 / INV-02 and stay finite

    /// <summary>
    /// BE / INV sweep: across a fixed-seed battery of random VALID inputs (random
    /// number of individuals, panels K ≥ 1 and SNPs J ≥ 1, random in-[0,1]
    /// reference frequencies and random 0/1/2 genotypes) EVERY emitted ancestry
    /// vector must obey the documented invariants — fractions sum to 1 EXACTLY
    /// (INV-01), each in [0, 1] (INV-02), none NaN or ±Infinity, one entry per
    /// panel. The Σ_m q_m f_mj denominators are guarded so no degenerate panel can
    /// inject a 0/0 NaN. This pins the total-function contract over the random
    /// interior, complementing the degenerate-boundary tests, and runs under
    /// `[CancelAfter]` so a hang in the EM loop is a hard failure, not an
    /// infinite wait. Expected invariants are the doc's (§2.4), not the code's.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void Ancestry_RandomValidInputs_AlwaysObeyInvariantsAndStayFinite()
    {
        for (int trial = 0; trial < 64; trial++)
        {
            int k = Rng.Next(1, 5);       // 1..4 reference panels
            int j = Rng.Next(1, 30);      // 1..29 SNPs
            int n = Rng.Next(1, 6);       // 1..5 individuals

            var refs = new (string, IReadOnlyList<double>)[k];
            for (int pop = 0; pop < k; pop++)
            {
                var freqs = new double[j];
                for (int snp = 0; snp < j; snp++)
                    freqs[snp] = Rng.NextDouble();
                refs[pop] = ($"P{pop}", freqs);
            }

            var inds = new (string, IReadOnlyList<int>)[n];
            for (int i = 0; i < n; i++)
            {
                var g = new int[j];
                for (int snp = 0; snp < j; snp++)
                    g[snp] = Rng.Next(0, 3);
                inds[i] = ($"i{i}", g);
            }

            List<PopulationGeneticsAnalyzer.AncestryProportion> results = null!;
            var act = () => results = PopulationGeneticsAnalyzer.EstimateAncestry(inds, refs, maxIterations: 100).ToList();

            act.Should().NotThrow("random valid genotypes and panels are always a defined input");

            results.Should().HaveCount(n, "every length-matching individual yields exactly one ancestry vector");
            foreach (var ap in results)
            {
                ap.Proportions.Should().HaveCount(k, "one ancestry fraction per reference panel (§3.2)");
                foreach (var (_, fraction) in ap.Proportions)
                {
                    double.IsNaN(fraction).Should().BeFalse("the guarded EM denominators never yield a 0/0 NaN");
                    double.IsInfinity(fraction).Should().BeFalse("bounded responsibilities over 2J never give ±Infinity");
                    fraction.Should().BeInRange(0.0, 1.0, "each ancestry fraction is a proportion (INV-02)");
                }
                FractionSum(ap).Should().BeApproximately(1.0, Tolerance,
                    because: "the EM divides by 2J so the ancestry fractions sum to 1 exactly (INV-01)");
            }
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  POP-ROH-001 — runs of homozygosity detection : fuzz targets
    //  Surface: FindROH(IEnumerable<(int Position, int Genotype)>,
    //                    minSnps, minLength, maxHeterozygotes, maxGap)
    // ═══════════════════════════════════════════════════════════════════

    #region POP-ROH-001 — runs of homozygosity

    #region Helpers — ROH track builders & a LOCALLY-seeded RNG

    /// <summary>
    /// Locally-seeded RNG for the ROH randomized sweeps — kept INDEPENDENT of the
    /// fixture-wide <see cref="Rng"/> so the ROH boundary battery is reproducible on
    /// its own and adding/removing other tests cannot shift these inputs.
    /// </summary>
    private static Random RohRng() => new(20260621);

    /// <summary>
    /// A homozygous track: <paramref name="count"/> SNPs at positions
    /// <c>0, spacing, 2·spacing, …</c>, all homozygous (genotype 0). The last SNP is
    /// at <c>(count − 1)·spacing</c>, so the track spans exactly that many bp.
    /// </summary>
    private static List<(int Position, int Genotype)> HomTrack(int count, int spacing) =>
        Enumerable.Range(0, count).Select(i => (i * spacing, 0)).ToList();

    #endregion

    #region Positive sanity — the documented worked example (§7.1)

    /// <summary>
    /// Positive control: the worked example pinned in the algorithm doc
    /// (Runs_Of_Homozygosity.md §7.1) — 100 homozygous SNPs at positions
    /// 0, 20000, …, 1_980_000 (spacing 20 kb &lt; the 1 Mb default maxGap, so no gap
    /// break; no heterozygotes). With the PLINK 1.9 defaults (minSnps = 100,
    /// minLength = 1_000_000) the single maximal run is
    /// (Start = 0, End = 1_980_000, SnpCount = 100): SnpCount 100 ≥ 100 and span
    /// 1_980_000 ≥ 1_000_000 (INV-01). Expected values are the doc's, hand-derived
    /// from the consecutive-runs model, NOT read off the implementation. This anchors
    /// the fuzz battery: before pinning that boundary input does no harm, confirm a
    /// known-good track reproduces the documented run exactly.
    /// </summary>
    [Test]
    public void FindRoh_DocumentedWorkedExample_YieldsSingleSpanningRun()
    {
        var snps = HomTrack(count: 100, spacing: 20_000); // 0 .. 1_980_000

        var runs = PopulationGeneticsAnalyzer.FindROH(snps).ToList(); // defaults: minSnps 100, minLength 1 Mb

        runs.Should().ContainSingle("an unbroken homozygous track is one maximal ROH (§7.1)");
        runs[0].Start.Should().Be(0, "the run starts at the first homozygous SNP (§7.1)");
        runs[0].End.Should().Be(1_980_000, "the run ends at the last homozygous SNP, 99·20000 (§7.1)");
        runs[0].SnpCount.Should().Be(100, "all 100 SNPs are homozygous and in one run (§7.1)");
    }

    #endregion

    #region BE — all-heterozygous track (no homozygous SNP can seed a run → ZERO runs)

    /// <summary>
    /// BE "all-heterozygous": when every genotype is the opposite call (1), no
    /// homozygous SNP exists to seed a run, so FindROH must return ZERO runs no
    /// matter how long the track or how permissive the thresholds — a heterozygous
    /// SNP never starts a homozygous run (Runs_Of_Homozygosity.md §6.1 "Leading
    /// heterozygotes"; lines 1560–1563). Verified across a locally-seeded sweep of
    /// track lengths and (deliberately trivial) minSnps = 1, minLength = 0 so the
    /// emptiness comes purely from the all-het structure, not from a threshold cut.
    /// No crash, no hang, and crucially never a phantom run over heterozygous SNPs.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void FindRoh_AllHeterozygous_YieldsNoRunsRegardlessOfLength()
    {
        var rng = RohRng();
        for (int trial = 0; trial < 64; trial++)
        {
            int count = rng.Next(1, 500);
            int spacing = rng.Next(1, 5_000);
            var allHet = Enumerable.Range(0, count).Select(i => (i * spacing, 1)).ToList();

            List<(int Start, int End, int SnpCount)> runs = null!;
            var act = () => runs = PopulationGeneticsAnalyzer
                .FindROH(allHet, minSnps: 1, minLength: 0, maxHeterozygotes: 0, maxGap: int.MaxValue)
                .ToList();

            act.Should().NotThrow("an all-heterozygous track is a defined input, not an error");
            runs.Should().BeEmpty(
                $"every one of {count} genotypes is the opposite call 1 → no homozygous SNP can seed a " +
                "run → zero ROH (§6.1); a het SNP never starts a homozygous run");
        }
    }

    #endregion

    #region BE — all-homozygous track (ONE maximal run, retained iff it clears both thresholds)

    /// <summary>
    /// BE "all-homozygous", retention PASSES: a homozygous track with no gap break is
    /// ONE maximal run [firstPos .. lastPos] of SnpCount = count. With minSnps and
    /// minLength chosen so BOTH thresholds are cleared, the single run is emitted with
    /// the hand-derived boundaries (INV-01). Here 50 SNPs spaced 100 bp span 0..4900;
    /// with minSnps = 50 and minLength = 4900 the SnpCount 50 ≥ 50 and span 4900 ≥
    /// 4900 both hold (the span is exactly the threshold), so the run is retained.
    /// Mixing homozygous codes (0 = hom-ref and 2 = hom-alt) must NOT break the run —
    /// only genotype 1 is the opposite call (§3.3).
    /// </summary>
    [Test]
    public void FindRoh_AllHomozygous_PassingBothThresholds_YieldsOneSpanningRun()
    {
        // 50 homozygous SNPs (alternating hom-ref/hom-alt codes), positions 0..4900.
        var track = Enumerable.Range(0, 50).Select(i => (i * 100, i % 2 == 0 ? 0 : 2)).ToList();

        var runs = PopulationGeneticsAnalyzer
            .FindROH(track, minSnps: 50, minLength: 4_900, maxHeterozygotes: 0, maxGap: 1_000)
            .ToList();

        runs.Should().ContainSingle("an unbroken all-homozygous track is exactly one ROH (§2.2)");
        runs[0].Start.Should().Be(0, "the run starts at the first homozygous SNP");
        runs[0].End.Should().Be(4_900, "the run ends at the last homozygous SNP, 49·100 (§2.2)");
        runs[0].SnpCount.Should().Be(50, "both hom-ref (0) and hom-alt (2) are homozygous; only 1 is opposite (§3.3)");
    }

    /// <summary>
    /// BE "all-homozygous", retention FAILS on count: an all-homozygous track that is
    /// maximal but holds FEWER than minSnps SNPs is DISCARDED, even though its physical
    /// span clears minLength ("a run passing length but not count → discarded", §6.1;
    /// both thresholds required, INV-01). 10 SNPs spaced 1000 bp span 0..9000 (≥
    /// minLength 0) but SnpCount 10 &lt; minSnps 11 → no run.
    /// </summary>
    [Test]
    public void FindRoh_AllHomozygous_BelowMinSnps_YieldsNoRun()
    {
        var track = HomTrack(count: 10, spacing: 1_000); // 0 .. 9000, 10 SNPs

        var runs = PopulationGeneticsAnalyzer
            .FindROH(track, minSnps: 11, minLength: 0, maxHeterozygotes: 0, maxGap: 1_000_000)
            .ToList();

        runs.Should().BeEmpty(
            "SnpCount 10 < minSnps 11 → the run fails the count threshold and is discarded, " +
            "even though its span 9000 ≥ minLength 0 (both thresholds required, §6.1 / INV-01)");
    }

    #endregion

    #region BE — minLength edge (span EXACTLY at / one below the threshold)

    /// <summary>
    /// BE "minLength edge": the retention rule is End − Start ≥ minLength, a
    /// NON-strict ≥ (Runs_Of_Homozygosity.md §2.4 INV-01; line 1532). A homozygous
    /// track of 5 SNPs spaced 100 bp spans 0..400, i.e. End − Start = 400. So:
    ///   • minLength = 400 (EXACTLY the span) → RETAINED (400 ≥ 400);
    ///   • minLength = 399 (one below)        → RETAINED (400 ≥ 399);
    ///   • minLength = 401 (one above the span) → DISCARDED (400 ≥ 401 is false).
    /// minSnps is held at 5 (= the SnpCount) so the COUNT threshold never confounds the
    /// LENGTH boundary. This pins the exact ≥ semantics of the length cut — an
    /// off-by-one to a strict &gt; would wrongly drop the exactly-at-threshold run.
    /// </summary>
    [TestCase(399, true, TestName = "FindRoh_MinLengthEdge_OneBelowSpan_Retained")]
    [TestCase(400, true, TestName = "FindRoh_MinLengthEdge_ExactlySpan_Retained")]
    [TestCase(401, false, TestName = "FindRoh_MinLengthEdge_OneAboveSpan_Discarded")]
    public void FindRoh_MinLengthEdge_RetainsIffSpanAtLeastThreshold(int minLength, bool retained)
    {
        var track = HomTrack(count: 5, spacing: 100); // 0,100,200,300,400 → span 400, 5 SNPs

        var runs = PopulationGeneticsAnalyzer
            .FindROH(track, minSnps: 5, minLength: minLength, maxHeterozygotes: 0, maxGap: 1_000)
            .ToList();

        if (retained)
        {
            runs.Should().ContainSingle(
                $"span 400 ≥ minLength {minLength} → the run is retained (End − Start ≥ minLength, INV-01)");
            runs[0].Start.Should().Be(0);
            runs[0].End.Should().Be(400, "the run ends at the last homozygous SNP, 4·100");
            runs[0].SnpCount.Should().Be(5);
        }
        else
        {
            runs.Should().BeEmpty(
                $"span 400 < minLength {minLength} → the run fails the length threshold (§6.1 / INV-01)");
        }
    }

    /// <summary>
    /// BE twin edge — minSnps boundary: the count rule is SnpCount ≥ minSnps, also a
    /// non-strict ≥ (§2.4 INV-01; line 1532). A 6-SNP homozygous track (span large
    /// enough to clear minLength) is RETAINED at minSnps = 6 (6 ≥ 6) and at
    /// minSnps = 5, but DISCARDED at minSnps = 7 (6 ≥ 7 is false). This pins the exact
    /// ≥ semantics of the count cut, the companion of the length edge above.
    /// </summary>
    [TestCase(5, true, TestName = "FindRoh_MinSnpsEdge_OneBelowCount_Retained")]
    [TestCase(6, true, TestName = "FindRoh_MinSnpsEdge_ExactlyCount_Retained")]
    [TestCase(7, false, TestName = "FindRoh_MinSnpsEdge_OneAboveCount_Discarded")]
    public void FindRoh_MinSnpsEdge_RetainsIffCountAtLeastThreshold(int minSnps, bool retained)
    {
        var track = HomTrack(count: 6, spacing: 100); // 6 SNPs, span 500

        var runs = PopulationGeneticsAnalyzer
            .FindROH(track, minSnps: minSnps, minLength: 0, maxHeterozygotes: 0, maxGap: 1_000)
            .ToList();

        if (retained)
            runs.Should().ContainSingle(
                $"SnpCount 6 ≥ minSnps {minSnps} → retained (SnpCount ≥ minSnps, INV-01)")
                .Which.SnpCount.Should().Be(6);
        else
            runs.Should().BeEmpty(
                $"SnpCount 6 < minSnps {minSnps} → fails the count threshold (§6.1 / INV-01)");
    }

    #endregion

    #region BE — empty track & the 0/−1 argument-validation gates

    /// <summary>
    /// BE "empty": an empty genotype track has nothing to scan, so FindROH yields NO
    /// runs (Runs_Of_Homozygosity.md §6.1 "Empty genotype input → no runs"; line 1514).
    /// No crash, no hang, no phantom run.
    /// </summary>
    [Test]
    public void FindRoh_EmptyTrack_YieldsNoRuns()
    {
        var runs = PopulationGeneticsAnalyzer
            .FindROH(Array.Empty<(int, int)>(), minSnps: 1, minLength: 0)
            .ToList();

        runs.Should().BeEmpty("an empty track has nothing to scan → no runs (§6.1)");
    }

    /// <summary>
    /// BE "−1"/"0" validation gates: null genotypes raise ArgumentNullException; and
    /// minSnps &lt; 1 (0 or negative) or a NEGATIVE minLength / maxHeterozygotes /
    /// maxGap raise the *documented, intentional* ArgumentOutOfRangeException
    /// (Runs_Of_Homozygosity.md §3.3; lines 1492–1501). These are eager — thrown when
    /// FindROH is CALLED, before the deferred iterator is enumerated (§5.2) — so the
    /// act under test is the bare call, not a later .ToList(). This pins the BE
    /// "0"/"−1" boundary as a clean validation error, never a silent garbage run.
    /// </summary>
    [Test]
    public void FindRoh_InvalidArguments_ThrowDocumentedValidationExceptions()
    {
        // null genotypes — ArgumentNullException, eagerly.
        var nullCall = () => PopulationGeneticsAnalyzer.FindROH(null!);
        nullCall.Should().Throw<ArgumentNullException>("null genotypes are rejected eagerly (§3.3)");

        var track = HomTrack(count: 3, spacing: 100);

        // minSnps must be ≥ 1: 0 and −1 are out of range (the BE "0"/"−1" class).
        var minSnpsZero = () => PopulationGeneticsAnalyzer.FindROH(track, minSnps: 0);
        minSnpsZero.Should().Throw<ArgumentOutOfRangeException>("minSnps < 1 is invalid (§3.1, §3.3)");
        var minSnpsNeg = () => PopulationGeneticsAnalyzer.FindROH(track, minSnps: -1);
        minSnpsNeg.Should().Throw<ArgumentOutOfRangeException>("minSnps < 1 is invalid (§3.1, §3.3)");

        // Negative minLength / maxHeterozygotes / maxGap are each rejected independently.
        var negLength = () => PopulationGeneticsAnalyzer.FindROH(track, minLength: -1);
        negLength.Should().Throw<ArgumentOutOfRangeException>("minLength cannot be negative (§3.3)");
        var negHet = () => PopulationGeneticsAnalyzer.FindROH(track, maxHeterozygotes: -1);
        negHet.Should().Throw<ArgumentOutOfRangeException>("maxHeterozygotes cannot be negative (§3.3)");
        var negGap = () => PopulationGeneticsAnalyzer.FindROH(track, maxGap: -1);
        negGap.Should().Throw<ArgumentOutOfRangeException>("maxGap cannot be negative (§3.3)");
    }

    #endregion

    #region BE — randomized boundary sweep (every emitted run obeys INV-01 / INV-02 / INV-03)

    /// <summary>
    /// BE / INV sweep: across a locally-seeded battery of random VALID tracks (random
    /// SNP count, random positions, random genotypes drawn from {0,1,2} with the
    /// het-rich and hom-rich extremes both represented, random thresholds) EVERY
    /// emitted run must obey the documented invariants, derived independently from the
    /// model (§2.4), NOT read off the code:
    ///   • INV-01 — SnpCount ≥ minSnps AND End − Start ≥ minLength;
    ///   • INV-02 — the emitted [Start, End] interval contains ≤ maxHeterozygotes
    ///     opposite (genotype == 1) SNPs and no inter-SNP gap &gt; maxGap;
    ///   • INV-03 — runs are emitted in ascending, NON-overlapping Start order with
    ///     Start ≤ End.
    /// Positions are deduplicated and re-sorted so the "same runs on unsorted input"
    /// contract (§6.1) is implicitly exercised by feeding the analyzer a shuffled copy.
    /// Run under `[CancelAfter]` so any hang in the linear scan is a hard failure.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void FindRoh_RandomValidTracks_EveryEmittedRunObeysInvariants()
    {
        var rng = RohRng();
        for (int trial = 0; trial < 200; trial++)
        {
            int count = rng.Next(0, 60);
            // Distinct ascending positions, then shuffle so FindROH must sort internally.
            var positions = new SortedSet<int>();
            while (positions.Count < count)
                positions.Add(rng.Next(0, 50_000));
            var sortedPos = positions.ToList();

            // Bias genotypes toward homozygous so non-trivial runs actually form,
            // but include heterozygous (1) calls to exercise the opposite-genotype break.
            var byPos = sortedPos.ToDictionary(p => p, _ => rng.Next(10) < 7 ? rng.Next(2) * 2 : 1);
            var shuffled = sortedPos.OrderBy(_ => rng.Next()).Select(p => (p, byPos[p])).ToList();

            int minSnps = rng.Next(1, 6);
            int minLength = rng.Next(0, 5_000);
            int maxHet = rng.Next(0, 3);
            int maxGap = rng.Next(1, 20_000);

            List<(int Start, int End, int SnpCount)> runs = null!;
            var act = () => runs = PopulationGeneticsAnalyzer
                .FindROH(shuffled, minSnps, minLength, maxHet, maxGap)
                .ToList();

            act.Should().NotThrow("a random valid track with valid thresholds is always a defined input");

            int previousEnd = int.MinValue;
            foreach (var run in runs)
            {
                run.Start.Should().BeLessThanOrEqualTo(run.End, "a run spans first → last SNP (INV-03)");
                run.SnpCount.Should().BeGreaterThanOrEqualTo(minSnps,
                    "a retained run holds at least minSnps homozygous SNPs (INV-01)");
                (run.End - run.Start).Should().BeGreaterThanOrEqualTo(minLength,
                    "a retained run spans at least minLength bp (INV-01)");
                run.Start.Should().BeGreaterThan(previousEnd,
                    "runs are emitted in ascending, non-overlapping Start order (INV-03)");
                previousEnd = run.End;

                // INV-02: independently recount the opposite genotypes and the max gap
                // inside the emitted [Start, End] interval from the ORIGINAL sorted data.
                var inside = sortedPos.Where(p => p >= run.Start && p <= run.End).ToList();
                int hetInside = inside.Count(p => byPos[p] == 1);
                hetInside.Should().BeLessThanOrEqualTo(maxHet,
                    "the emitted interval holds at most maxHeterozygotes opposite genotypes (INV-02)");

                int maxGapInside = 0;
                for (int k = 1; k < inside.Count; k++)
                    maxGapInside = Math.Max(maxGapInside, inside[k] - inside[k - 1]);
                maxGapInside.Should().BeLessThanOrEqualTo(maxGap,
                    "no inter-SNP gap inside the emitted run exceeds maxGap (INV-02)");
            }
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  POP-SELECT-001 — iHS selection-signature detection : fuzz targets
    //  Doc: docs/algorithms/Population_Genetics/Integrated_Haplotype_Score.md
    //  Surface:
    //    CalculateEhh(IReadOnlyList<string>)                       → EHH ∈ [0,1]
    //    CalculateIHS(IReadOnlyList<string>, IReadOnlyList<int>,   → IhsResult
    //                 int coreIndex)
    //  Strategy: BE (Boundary Exploitation) — checklist row 209
    //    "neutral, fixed locus, single locus".
    // ═══════════════════════════════════════════════════════════════════

    #region POP-SELECT-001 — integrated haplotype score (iHS)

    // ───────────────────────────────────────────────────────────────────
    //  Unit: POP-SELECT-001 — integrated haplotype score (iHS) selection
    //        signature detection (PopGen).
    //  Checklist: docs/checklists/03_FUZZING.md, row 209 (BE).
    //
    //  What the algorithm doc fixes (Integrated_Haplotype_Score.md):
    //    • §2.2  EHH_c = Σ_h C(n_h,2)/C(n_c,2); iHH = trapezoidal area under
    //            the EHH-vs-position curve, both directions, truncated where
    //            EHH first drops below 0.05; unstandardized iHS = ln(iHH_A/iHH_D)
    //            (Voight et al. 2006 sign — ancestral numerator, derived denom).
    //    • §2.4  INV-01 EHH ∈ [0,1] (1 for a single chromosome, 0 for all-distinct);
    //            INV-02 iHH ≥ 0; INV-03 balanced EHH decay ⇒ unstandardized iHS = 0;
    //            INV-04 ln(iHH_A/iHH_D) = −ln(iHH_D/iHH_A).
    //    • §3.3  core ∈ {'0','1'} else ArgumentException; null → ArgumentNullException;
    //            coreIndex out of range → ArgumentOutOfRangeException; inconsistent
    //            lengths or a MONOMORPHIC focal SNP → ArgumentException; empty EHH
    //            sample → 0; single-chromosome sample → EHH = 1.
    //    • §6.1  edge table: single chromosome ⇒ EHH 1; empty sample ⇒ EHH 0;
    //            balanced decay ⇒ iHS 0; monomorphic focal SNP ⇒ ArgumentException;
    //            core ∉ {'0','1'} ⇒ ArgumentException.
    //    • §7.1  worked example: haplotypes {AA1GG×3, TC0TC, GA0AG, CT0CA},
    //            positions {0,10,20,30,40}, coreIndex 2 ⇒ iHH_A = 10, iHH_D = 40,
    //            unstandardized iHS = ln(10/40) = −1.386294361.
    //
    //  BE boundaries exercised for THIS row (row hint "neutral, fixed locus,
    //  single locus"):
    //    • neutral / balanced decay → unstandardized iHS = 0 (INV-03);
    //    • fixed (monomorphic) focal SNP → documented ArgumentException, NEVER a
    //      NaN/Infinity or a silent divide-by-zero;
    //    • single locus (one marker, no flanks) → iHH_A = iHH_D = 0 ⇒ a DEFINED
    //      0 score, not NaN from ln(0/0);
    //    • the EHH degenerate boundaries (empty sample → 0, single chromosome → 1);
    //    • a randomized boundary sweep: every random valid panel yields finite,
    //      in-range EHH/iHS obeying INV-01/02/03/04 — no crash, hang, NaN or ∞.
    //  Fuzz bar — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing": every input is
    //  EITHER a theory-correct value OR a documented validation exception; a raw
    //  runtime exception, a NaN/±Infinity score, or a hang is a bug.
    //
    //  Expected values below are derived INDEPENDENTLY from the doc's EHH/iHH model
    //  and hand-computed, NOT read off the implementation.
    // ───────────────────────────────────────────────────────────────────

    #region Helpers — a LOCALLY-seeded RNG for the iHS sweeps

    /// <summary>
    /// Locally-seeded RNG for the iHS randomized sweep — kept INDEPENDENT of the
    /// fixture-wide <see cref="Rng"/> (and of <see cref="RohRng"/>) so this boundary
    /// battery is reproducible on its own and adding/removing other tests cannot
    /// shift its inputs.
    /// </summary>
    private static Random IhsRng() => new(20260620);

    #endregion

    #region Positive sanity — the documented §7.1 worked example (iHH_A=10, iHH_D=40, iHS=ln 0.25)

    /// <summary>
    /// Positive control: the worked example pinned in the algorithm doc
    /// (Integrated_Haplotype_Score.md §7.1). Three identical derived chromosomes
    /// "AA1GG" and three all-distinct ancestral chromosomes ("TC0TC","GA0AG",
    /// "CT0CA"), markers spaced 10 at positions {0,10,20,30,40}, focal SNP at index 2.
    ///
    /// Hand-derivation (from §2.2, NOT from the code):
    ///   • Derived set is identical ⇒ EHH_D = 1 at every flank marker. Each
    ///     direction contributes ½(1+1)·10 + ½(1+1)·10 = 20 ⇒ iHH_D = 40.
    ///   • Ancestral set is all-distinct ⇒ EHH_A = 0 at the first flank marker;
    ///     each direction contributes ½(1+0)·10 = 5 then truncates (EHH < 0.05)
    ///     ⇒ iHH_A = 10.
    ///   • unstandardized iHS = ln(iHH_A/iHH_D) = ln(10/40) = ln(0.25)
    ///                        = −1.386294361 (Voight sign).
    /// The derived allele frequency is 3/6 = 0.5. This anchors the fuzz battery:
    /// before pinning the degenerate boundaries, confirm the canonical case is exact.
    /// </summary>
    [Test]
    public void CalculateIhs_DocumentedWorkedExample_MatchesHandDerivedScore()
    {
        var haplotypes = new[] { "AA1GG", "AA1GG", "AA1GG", "TC0TC", "GA0AG", "CT0CA" };
        var positions = new[] { 0, 10, 20, 30, 40 };

        var ihs = PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 2);

        ihs.IhhAncestral.Should().BeApproximately(10.0, 1e-9,
            "3 all-distinct ancestral chromosomes ⇒ EHH_A = 0 at the first flank, ½(1+0)·10 each side, truncated (§7.1)");
        ihs.IhhDerived.Should().BeApproximately(40.0, 1e-9,
            "3 identical derived chromosomes ⇒ EHH_D = 1 throughout, 20 per side (§7.1)");
        ihs.UnstandardizedIHS.Should().BeApproximately(Math.Log(0.25), 1e-9,
            "unstandardized iHS = ln(iHH_A/iHH_D) = ln(10/40) = ln(0.25) = −1.386294361 (Voight sign, §7.1)");
        ihs.DerivedAlleleFrequency.Should().BeApproximately(0.5, 1e-12,
            "3 of 6 core chromosomes carry the derived allele (§7.1)");
    }

    #endregion

    #region BE — neutral / balanced decay ⇒ unstandardized iHS = 0 (INV-03)

    /// <summary>
    /// BE "neutral": when the two core-allele subsets decay IDENTICALLY, iHH_A = iHH_D
    /// so the unstandardized iHS = ln(iHH_A/iHH_D) = ln(1) = 0 (INV-03, §2.4 / §6.1
    /// "Balanced decay ⇒ unstandardized iHS = 0"). Independent hand-derivation:
    /// two identical derived chromosomes "A1A" and two identical ancestral "C0C",
    /// positions {0,10,20}, core at index 1. BOTH subsets are within-group identical,
    /// so EHH = 1 at every flank marker for each ⇒ each side contributes ½(1+1)·10 = 10
    /// ⇒ iHH_A = iHH_D = 20. Hence iHS = ln(20/20) = 0 EXACTLY — the neutral, no-signal
    /// boundary — with derived frequency 2/4 = 0.5 and no NaN/Infinity.
    /// </summary>
    [Test]
    public void CalculateIhs_BalancedDecay_YieldsZeroScore()
    {
        var haplotypes = new[] { "A1A", "A1A", "C0C", "C0C" };
        var positions = new[] { 0, 10, 20 };

        var ihs = PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 1);

        ihs.IhhAncestral.Should().BeApproximately(20.0, 1e-9,
            "two identical ancestral chromosomes ⇒ EHH_A = 1 throughout ⇒ iHH_A = 10 + 10 (§2.2)");
        ihs.IhhDerived.Should().BeApproximately(20.0, 1e-9,
            "two identical derived chromosomes ⇒ EHH_D = 1 throughout ⇒ iHH_D = 10 + 10 (§2.2)");
        ihs.UnstandardizedIHS.Should().Be(0.0,
            "iHH_A = iHH_D ⇒ ln(iHH_A/iHH_D) = ln(1) = 0 — balanced decay, the neutral no-signal case (INV-03)");
        ihs.UnstandardizedIHS.Should().NotBe(double.NaN);
        ihs.DerivedAlleleFrequency.Should().BeApproximately(0.5, 1e-12);
    }

    /// <summary>
    /// BE "neutral", reciprocal sign (INV-04): swapping which group is identical must
    /// flip the score sign but not its magnitude, because ln(iHH_A/iHH_D) =
    /// −ln(iHH_D/iHH_A) (§2.4 INV-04, the Voight-vs-selscan sign note §2.2). Using the
    /// §7.1 structure with the core allele labels swapped — three identical "0" carriers
    /// and three all-distinct "1" carriers — yields iHH_A = 40, iHH_D = 10, so the score
    /// is ln(40/10) = +ln(4) = +1.386294361: the exact negation of the §7.1 score, and a
    /// POSITIVE value indicating a long ancestral haplotype. This pins the contract's sign
    /// semantics, not just |iHS|.
    /// </summary>
    [Test]
    public void CalculateIhs_AncestralLongHaplotype_YieldsPositiveReciprocalScore()
    {
        // §7.1 mirror: the identical trio now carries the ANCESTRAL allele '0'.
        var haplotypes = new[] { "AA0GG", "AA0GG", "AA0GG", "TC1TC", "GA1AG", "CT1CA" };
        var positions = new[] { 0, 10, 20, 30, 40 };

        var ihs = PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 2);

        ihs.IhhAncestral.Should().BeApproximately(40.0, 1e-9, "the identical trio now carries the ancestral allele");
        ihs.IhhDerived.Should().BeApproximately(10.0, 1e-9, "the all-distinct trio now carries the derived allele");
        ihs.UnstandardizedIHS.Should().BeApproximately(Math.Log(4.0), 1e-9,
            "ln(iHH_A/iHH_D) = ln(40/10) = +ln(4) = −ln(10/40): the exact negation of §7.1 (INV-04)");
        ihs.UnstandardizedIHS.Should().BeApproximately(-Math.Log(0.25), 1e-9,
            "the reciprocal-sign identity INV-04 pins this as the negation of the §7.1 score");
    }

    #endregion

    #region BE — fixed (monomorphic) focal SNP ⇒ documented ArgumentException (no NaN)

    /// <summary>
    /// BE "fixed locus": a MONOMORPHIC focal SNP — every core chromosome carrying the
    /// SAME allele — is rejected with the documented <see cref="ArgumentException"/>
    /// (§3.3, §6.1 "Monomorphic focal SNP ⇒ ArgumentException"; "iHS defined only for
    /// polymorphic SNPs"). iHS contrasts the two alleles, so with only one allele present
    /// there is nothing to contrast and the score is undefined — the analyzer must throw
    /// cleanly, NEVER return a 0/0 NaN or a divide-by-zero. Verified for an all-derived
    /// ('1') core AND an all-ancestral ('0') core, so neither fixation direction slips
    /// through.
    /// </summary>
    [TestCase('1', TestName = "CalculateIhs_FixedLocus_AllDerived_Throws")]
    [TestCase('0', TestName = "CalculateIhs_FixedLocus_AllAncestral_Throws")]
    public void CalculateIhs_MonomorphicFocalSnp_ThrowsArgumentException(char core)
    {
        // All four chromosomes carry the SAME core allele at index 1 → monomorphic.
        var haplotypes = new[] { $"A{core}A", $"C{core}T", $"G{core}G", $"A{core}C" };
        var positions = new[] { 0, 10, 20 };

        var act = () => PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 1);

        act.Should().Throw<ArgumentException>(
            "a fixed/monomorphic focal SNP carries only one allele → iHS is undefined and must be rejected, " +
            "never reported as a 0/0 NaN (§3.3 / §6.1)");
    }

    /// <summary>
    /// BE "core ∉ {'0','1'}": the core allele must be the polarized encoding '0'
    /// (ancestral) or '1' (derived); ANY other core character is rejected with the
    /// documented <see cref="ArgumentException"/> (§3.3, §6.1 "Core allele ∉ {'0','1'}
    /// ⇒ ArgumentException"). Flanking markers may use any alphabet (§5.2), so the guard
    /// must fire on the CORE column specifically, not on the flanks. Here the flanks are
    /// ordinary DNA but the core column holds 'A' — an invalid polarization symbol.
    /// </summary>
    [Test]
    public void CalculateIhs_NonBinaryCoreAllele_ThrowsArgumentException()
    {
        var haplotypes = new[] { "GAT", "CAT", "GAC" }; // core column (index 1) is all 'A'
        var positions = new[] { 0, 10, 20 };

        var act = () => PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 1);

        act.Should().Throw<ArgumentException>(
            "the core allele must be '0' or '1'; any other symbol is an undefined polarization and is rejected (§3.3)");
    }

    #endregion

    #region BE — single locus (one marker, no flanks) ⇒ a DEFINED zero score, not NaN

    /// <summary>
    /// BE "single locus": when there is exactly ONE marker — the focal SNP itself, with no
    /// flanking markers to decay over — neither direction has a marker to integrate, so
    /// iHH_A = iHH_D = 0. The doc's score is ln(iHH_A/iHH_D); with both areas 0 that is the
    /// degenerate ln(0/0). The contract (§3.3 "An empty EHH sample returns 0"; INV-02
    /// "iHH ≥ 0") requires a DEFINED result, so the implementation collapses the 0/0 to a
    /// score of 0 rather than emitting NaN. The focal SNP is still polymorphic ('1' and
    /// '0' both present, derived frequency 1/2) so the monomorphic guard does NOT fire.
    /// This pins the single-locus boundary: no flank ⇒ zero integrated area ⇒ a clean 0,
    /// never NaN or ±Infinity.
    /// </summary>
    [Test]
    public void CalculateIhs_SingleLocus_NoFlanks_YieldsDefinedZeroScore()
    {
        var haplotypes = new[] { "1", "0" }; // one marker = the core; polymorphic
        var positions = new[] { 5 };

        var ihs = PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 0);

        ihs.IhhAncestral.Should().Be(0.0, "no flanking marker exists to integrate over → iHH_A = 0 (INV-02)");
        ihs.IhhDerived.Should().Be(0.0, "no flanking marker exists to integrate over → iHH_D = 0 (INV-02)");
        ihs.UnstandardizedIHS.Should().Be(0.0,
            "with both integrated areas 0 the score collapses to a DEFINED 0, not the ln(0/0) NaN (§3.3)");
        ihs.UnstandardizedIHS.Should().NotBe(double.NaN);
        double.IsInfinity(ihs.UnstandardizedIHS).Should().BeFalse("a single-locus score must be finite");
        ihs.DerivedAlleleFrequency.Should().BeApproximately(0.5, 1e-12, "1 of 2 core chromosomes is derived");
    }

    #endregion

    #region BE — EHH degenerate boundaries (empty sample → 0, single chromosome → 1)

    /// <summary>
    /// BE "EHH boundaries" (§3.3 / §6.1 / INV-01): the EHH primitive itself must be
    /// defined at its degenerate endpoints —
    ///   • an EMPTY sample has no pairs ⇒ EHH = 0 (boundary of the C(n,2) ratio);
    ///   • a SINGLE chromosome is trivially homozygous ⇒ EHH = 1;
    ///   • all-identical chromosomes ⇒ every pair matches ⇒ EHH = 1;
    ///   • all-DISTINCT chromosomes ⇒ no pair matches ⇒ EHH = 0;
    ///   • a hand-checkable mixed case: {AB, AB, CD} ⇒ one matching pair of three
    ///     ⇒ EHH = C(2,2)/C(3,2) = 1/3 (§2.2 EHH = Σ C(n_h,2)/C(n_c,2)).
    /// A null sample throws the documented ArgumentNullException. Every EHH value lies in
    /// [0,1] (INV-01) with no NaN/Infinity.
    /// </summary>
    [Test]
    public void CalculateEhh_DegenerateSamples_AreDefinedAndInUnitInterval()
    {
        PopulationGeneticsAnalyzer.CalculateEhh(Array.Empty<string>())
            .Should().Be(0.0, "an empty sample has no pairs ⇒ EHH = 0 (§3.3 / §6.1)");

        PopulationGeneticsAnalyzer.CalculateEhh(new[] { "ACGT" })
            .Should().Be(1.0, "a single chromosome is trivially homozygous ⇒ EHH = 1 (§6.1 / INV-01)");

        PopulationGeneticsAnalyzer.CalculateEhh(new[] { "ACGT", "ACGT", "ACGT" })
            .Should().Be(1.0, "all-identical chromosomes ⇒ every pair matches ⇒ EHH = 1 (INV-01)");

        PopulationGeneticsAnalyzer.CalculateEhh(new[] { "AAAA", "CCCC", "GGGG" })
            .Should().Be(0.0, "all-distinct chromosomes ⇒ no pair matches ⇒ EHH = 0 (INV-01)");

        PopulationGeneticsAnalyzer.CalculateEhh(new[] { "AB", "AB", "CD" })
            .Should().BeApproximately(1.0 / 3.0, 1e-12,
                "one identical pair of three ⇒ EHH = C(2,2)/C(3,2) = 1/3 (§2.2)");

        var nullCall = () => PopulationGeneticsAnalyzer.CalculateEhh(null!);
        nullCall.Should().Throw<ArgumentNullException>("a null sample is rejected (§3.3)");
    }

    #endregion

    #region BE — null / out-of-range argument-validation gates

    /// <summary>
    /// BE "0"/"−1"/null validation gates for CalculateIHS (§3.3): null haplotypes or null
    /// positions raise <see cref="ArgumentNullException"/>; a coreIndex of −1 or one at/past
    /// positions.Count raises <see cref="ArgumentOutOfRangeException"/>; an empty haplotype
    /// list or rows of inconsistent length raise <see cref="ArgumentException"/>. Each is the
    /// documented, intentional validation error — never a silent miscompute or an
    /// IndexOutOfRange crash.
    /// </summary>
    [Test]
    public void CalculateIhs_InvalidArguments_ThrowDocumentedValidationExceptions()
    {
        var haps = new[] { "A1A", "C0C" };
        var positions = new[] { 0, 10, 20 }.Take(haps[0].Length).ToArray(); // {0,10,20}

        var nullHaps = () => PopulationGeneticsAnalyzer.CalculateIHS(null!, positions, 1);
        nullHaps.Should().Throw<ArgumentNullException>("null haplotypes are rejected (§3.3)");

        var nullPos = () => PopulationGeneticsAnalyzer.CalculateIHS(haps, null!, 1);
        nullPos.Should().Throw<ArgumentNullException>("null positions are rejected (§3.3)");

        var negCore = () => PopulationGeneticsAnalyzer.CalculateIHS(haps, positions, -1);
        negCore.Should().Throw<ArgumentOutOfRangeException>("coreIndex −1 is out of range (§3.3)");

        var bigCore = () => PopulationGeneticsAnalyzer.CalculateIHS(haps, positions, positions.Length);
        bigCore.Should().Throw<ArgumentOutOfRangeException>("coreIndex == positions.Count is out of range (§3.3)");

        var emptyHaps = () => PopulationGeneticsAnalyzer.CalculateIHS(Array.Empty<string>(), positions, 1);
        emptyHaps.Should().Throw<ArgumentException>("an empty haplotype list has nothing to score (§3.3)");

        // One row is shorter than positions.Count → inconsistent length.
        var ragged = () => PopulationGeneticsAnalyzer.CalculateIHS(new[] { "A1A", "C0" }, positions, 1);
        ragged.Should().Throw<ArgumentException>("each haplotype must have one allele per position (§3.3)");
    }

    #endregion

    #region BE — randomized boundary sweep (every random panel: finite, in-range, INV-01/02/03/04)

    /// <summary>
    /// BE / INV sweep: across a locally-seeded battery of random VALID phased panels
    /// (random chromosome count, random marker count, random binary core column kept
    /// polymorphic, random flank alleles, random ascending positions, random coreIndex)
    /// EVERY result must obey the invariants derived independently from §2.4, NOT read off
    /// the code:
    ///   • INV-01 — every intermediate EHH ∈ [0,1] and is finite;
    ///   • INV-02 — iHH_A ≥ 0 and iHH_D ≥ 0;
    ///   • the unstandardized score is finite (never NaN/±Infinity) and equals
    ///     ln(iHH_A/iHH_D) whenever both areas are positive, else the documented 0;
    ///   • INV-04 — recomputing with the core column's allele labels FLIPPED negates the
    ///     score (ln(iHH_A/iHH_D) = −ln(iHH_D/iHH_A)), within fp tolerance.
    /// Run under [CancelAfter] so any hang in the O(n·h) outward integration is a hard
    /// failure. This is the core fuzz bar: garbage-shaped but valid input never crashes,
    /// hangs, or produces a non-finite or out-of-range result.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void CalculateIhs_RandomValidPanels_AlwaysFiniteInRangeAndObeyInvariants()
    {
        const string flankAlphabet = "ACGT";
        var rng = IhsRng();

        for (int trial = 0; trial < 200; trial++)
        {
            int markers = rng.Next(1, 12);
            int chromosomes = rng.Next(2, 16);
            int coreIndex = rng.Next(0, markers);

            // Ascending, strictly-increasing positions.
            var positions = new int[markers];
            int p = rng.Next(0, 50);
            for (int m = 0; m < markers; m++)
            {
                positions[m] = p;
                p += rng.Next(1, 40);
            }

            // Build a polymorphic core column: at least one '0' and one '1'.
            var coreColumn = new char[chromosomes];
            for (int c = 0; c < chromosomes; c++)
                coreColumn[c] = rng.Next(2) == 0 ? '0' : '1';
            coreColumn[0] = '0';
            coreColumn[1] = '1';

            // Assemble each haplotype: random DNA flanks, the binary core at coreIndex.
            var haplotypes = new string[chromosomes];
            for (int c = 0; c < chromosomes; c++)
            {
                var chars = new char[markers];
                for (int m = 0; m < markers; m++)
                    chars[m] = m == coreIndex ? coreColumn[c] : flankAlphabet[rng.Next(flankAlphabet.Length)];
                haplotypes[c] = new string(chars);
            }

            PopulationGeneticsAnalyzer.IhsResult result = default;
            var act = () => result = PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex);
            act.Should().NotThrow("a random polymorphic, length-consistent panel is always a defined input");

            result.IhhAncestral.Should().BeGreaterThanOrEqualTo(0.0, "iHH is a sum of non-negative trapezoid areas (INV-02)");
            result.IhhDerived.Should().BeGreaterThanOrEqualTo(0.0, "iHH is a sum of non-negative trapezoid areas (INV-02)");
            double.IsNaN(result.IhhAncestral).Should().BeFalse("iHH_A must be finite");
            double.IsNaN(result.IhhDerived).Should().BeFalse("iHH_D must be finite");

            double s = result.UnstandardizedIHS;
            double.IsNaN(s).Should().BeFalse("a score over a valid panel is never NaN (the fuzz bar)");
            double.IsInfinity(s).Should().BeFalse("a score over a valid panel is never ±Infinity (the fuzz bar)");

            if (result.IhhAncestral > 0 && result.IhhDerived > 0)
                s.Should().BeApproximately(
                    Math.Log(result.IhhAncestral / result.IhhDerived), 1e-9,
                    "with both areas positive the score is exactly ln(iHH_A/iHH_D) (§2.2)");
            else
                s.Should().Be(0.0, "with a zero integrated area the documented score is a defined 0 (§3.3)");

            result.DerivedAlleleFrequency.Should().BeInRange(0.0, 1.0, "a derived allele frequency is a proportion");

            // INV-04: flip the core labels → the score must negate.
            var flipped = new string[chromosomes];
            for (int c = 0; c < chromosomes; c++)
            {
                var chars = haplotypes[c].ToCharArray();
                chars[coreIndex] = chars[coreIndex] == '0' ? '1' : '0';
                flipped[c] = new string(chars);
            }

            var flippedResult = PopulationGeneticsAnalyzer.CalculateIHS(flipped, positions, coreIndex);

            // ln(iHH_A/iHH_D) = −ln(iHH_D/iHH_A): flipping swaps the two integrated areas,
            // so the score must be the exact negation (both 0 when either area vanishes).
            if ((result.IhhAncestral > 0 && result.IhhDerived > 0))
                flippedResult.UnstandardizedIHS.Should().BeApproximately(-s, 1e-9,
                    "flipping the ancestral/derived labels negates ln(iHH_A/iHH_D) (INV-04)");
        }
    }

    #endregion

    #endregion
}
