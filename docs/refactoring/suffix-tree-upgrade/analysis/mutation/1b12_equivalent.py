# Mutants of the audit of 1b.12 that no test kills, argued (id -> argument). Quoted in PROGRESS-S1.md (S1b 2026-11-23 work package 1b.12).
# EQUIVALENT: the mutant behaves like the original on every input that a caller or a test can reach (the argument says why). Both survived stage 1 (the three new classes) and stage 2 (the existing contract and facade classes).
EQUIVALENT = {
    '1b12_r2_Q01': 'the facade hands the index a copy of the locus (IndexLocus copy = locus; TryParent(in copy, out parent)): the index reads its locus argument into a handle before it writes parent, so with the copy it reads a variable that is never the result; the answer is the one of the original in the aliased and in the separate-variable call, and nothing else observes the copy. It is the fix applied one layer higher, a pure refactoring.',
    '1b12_r2_Q02': 'the symbol overload of TryExtendRight hands Extend a copy of the locus: Extend reads the locus before its first write to result (A01 restores the bug there and is killed), so the copy changes nothing that is observable in either call form; the same refactoring, one layer lower.',
}
