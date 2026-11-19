# Mutants of the audit of 1b.7f that no test kills, in two groups (id -> argument). Read by summarize_1b7f.py and quoted in the ADR subsection.
# EQUIVALENT: the mutant behaves like the original in every state that a caller or a test can reach (the argument says why); DEFENCE: a second line of defence behind a check that is tested at another level, or an effect that no assertion
# of the language runtime can see: not equivalent, not killed, reported as such.
EQUIVALENT = {}
DEFENCE = {}
