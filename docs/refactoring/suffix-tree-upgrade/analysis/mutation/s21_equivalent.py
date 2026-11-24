# Mutants of the audit of S2.1 that no test kills, argued (id -> argument). Quoted in PROGRESS-S2.md (S2 2026-11-24 work package S2.1).
# EQUIVALENT: the mutant behaves like the original on every input that a caller or a test can reach (the argument says why). All three survived stage 1 (the owning class) and stage 2 (the four classes).
EQUIVALENT = {
    's21_K03': 'PeriodicAkB.RankOf, `position >= q p` -> `position > q p`: the two differ only at position = q p. There the original answers n - position = t (the run a^t of the last partial period), the mutant falls to the block formula '
               'with phase 0 and c = q, Rank(0, q - 1 - q) = 1 + t + 0 - 1 = t: the same number. (For t = 0, q p = n and the first line of the method answers 0 before either condition.)',
    's21_D03': 'DeBruijnRepeat.WindowLcp, `common < d` -> `common < d - 1`: the loop counts the equal leading digits of two windows. It is called only for suffixes of different rotations (PairLcp tests the same rotation first), whose '
               'windows are different words of d digits, so they differ somewhere and the count is at most d - 1: the loop ends at the same value under either bound.',
    's21_S15': 'SparseContainer.PsiOf, `position + 1 >= n` -> `position + 1 > n`: the two differ only at position + 1 = n, the suffix of the last symbol, where the original answers 0 and the mutant answers RankOf(n), which is 0 '
               '(the empty suffix has rank 0).',
}
