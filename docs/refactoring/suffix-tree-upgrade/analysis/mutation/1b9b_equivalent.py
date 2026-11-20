# Mutants of the audit of 1b.9b that neither the unit tests of the new code nor the budgeted harness kill, argued (id -> argument). Read by summarize_1b9b.py and quoted in PROGRESS-S1.md.
# EQUIVALENT: the mutant behaves like the original on every input that a caller or a test can reach (the argument says why); DEFENCE: a second line of defence behind a check that holds the same rule at the same read: not
# equivalent on every input, not killed, reported as such; KILLED_BY_UNIT_TESTS: killed by an existing test class that the harness cannot replace (none in this audit: the survivors of round 1 that could be killed got a test, see the round 2 rows).
EQUIVALENT = {
}
DEFENCE = {}
KILLED_BY_UNIT_TESTS = {}
