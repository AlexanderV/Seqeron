# Mutants of the audit of 1b.9b that neither the unit tests of the new code nor the budgeted harness kill, argued (id -> argument). Read by summarize_1b9b.py and quoted in PROGRESS-S1.md.
# EQUIVALENT: the mutant behaves like the original on every input that a caller or a test can reach (the argument says why); DEFENCE: a second line of defence behind a check that holds the same rule at the same read: not
# equivalent on every input, not killed, reported as such; KILLED_BY_UNIT_TESTS: killed by an existing test class that the harness cannot replace (none in this audit: the survivors of round 1 that could be killed got a test, see the round 2 rows).
EQUIVALENT = {
    '1b9b_B06': 'BUILD_INFO is a constant text of about 300 bytes with a few numbers (the format, the writer, n, sigma, the widths, the codec names): no input of the builder makes it reach the reserved 4,096 bytes, so the guard that the mutant leaves out is never true on a reachable input; the comparison itself is tested through its own function (B03, B04 and B05 are killed); the removal of the call is the one mutation of this guard that has no observable effect',
    '1b9b_F04': 'a call of an in-memory tree that leaves stripe 0 without having entered it lowers the count of that stripe below zero at its first call, and it never returns to zero; the count is read only by ReleaseWhenIdle, which an in-memory tree reaches only through Leave after its disposal, and which returns while any count is not zero, so nothing is ever released and the observable behaviour (a call in flight completes, a later call throws, no release) is that of the original; the difference is a counter that nothing reads',
}
DEFENCE = {}
KILLED_BY_UNIT_TESTS = {}
