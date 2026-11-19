# Mutants of the audit of 1b.8 that the harness does not kill, in two groups (id -> argument). Read by summarize_1b8.py and quoted in PROGRESS-S1.md.
# EQUIVALENT: the mutant behaves like the original in every state that a caller or a test can reach (the argument says why); DEFENCE: a second line of defence behind a check that holds the same rule one call later, or an
# effect that no assertion can see: not equivalent in every input, not killed, reported as such.
EQUIVALENT = {
    '1b8_A01': 'a run of the maximum LCP value can only begin at a rank of at least 2: LCP[1] is 0 (rank 1 follows the empty suffix) and the maximum is at least 1 whenever there is a longest repeat, so runStart = rank - 1 is at least 1 and `runStart >= 0` and `runStart > 0` select the same runs',
    '1b8_M01': 'SuffixTreeAlgorithms.FindMaximalExactMatches begins with ArgumentNullException.ThrowIfNull(query) and ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 1): the same exception type and parameter name come from the generic algorithm one call deeper',
    '1b8_M02': 'SuffixTreeAlgorithms.FindMaximalUniqueMatches validates minLength and the value of uniqueness at its first lines (SuffixTreeAlgorithms.cs lines 285-288): the same exception, the same parameter name',
    '1b8_M04': 'SuffixTreeAlgorithms.FindExactMatchAnchors begins with ArgumentNullException.ThrowIfNull(query) (SuffixTreeAlgorithms.cs line 105): the same exception and parameter name',
    '1b8_S07': 'for the empty text the alphabet is empty, so the first symbol of any non-empty pattern is not in the symbol map and both descents (in memory and mapped) return false before the tables are touched: the guard is an early exit',
    '1b8_S09': 'EdgeEquals begins with `expected.Length > text.Length - start` and returns false: with expected = pat.Slice(c) and start the same, it is the test `m - c > n - start` of the line before; both give the same false',
    '1b8_N01': 'CollectLeaves is called by the generic algorithms for the node of a match of length at least 1 (the longest common substrings, the maximal exact matches, the child of such a node): never for the root or for the terminator leaf of the root, the only nodes whose interval holds rank 0',
    '1b8_L01': 'PrevLess(0, bound) is reached only with bound 0 (PreviousSmaller(1): LCP[1] = 0), where it returns -1 either way, or from Lca with lo = 0, which returns the root first because the minimum over [1, hi] is LCP[1] = 0',
    '1b8_L02': 'NextLessCore begins with `if (pos >= _length) return _length;`: the test of NextLess at the same value is a copy of it',
    '1b8_O02': 'a negative LCP value cannot reach LeftmostOccurrence in the Checked mode: the 32-bit reader of CompactLcp refuses `(uint)v > (uint)Length` at every read (a negative value is a large unsigned one), and the 8- and 16-bit readers have no negative values; the check is a second line for the Trusted mode',
}
DEFENCE = {
}
