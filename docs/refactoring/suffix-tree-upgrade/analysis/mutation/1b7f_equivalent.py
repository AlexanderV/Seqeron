# Mutants of the audit of 1b.7f that no test kills, in two groups (id -> argument). Read by summarize_1b7f.py and quoted in the ADR subsection.
# EQUIVALENT: the mutant behaves like the original in every state that a caller or a test can reach (the argument says why); DEFENCE: a second line of defence behind a check that is tested at another level, or an effect that no assertion
# of the language runtime can see: not equivalent, not killed, reported as such. Ids V10, B102 and B118 are the survivors of the audit of 1b.7e that its orchestrator left for this package.
EQUIVALENT = {
    'D01': 'a key equal to the symbol leaves the loop at the `break` two lines above, so `key > symbol` and `key >= symbol` select the same inputs',
    'D08': 'EdgeEquals begins with the same bound (`expected.Length > text.Length - start` returns false): the leaf test and EdgeEquals give the same false',
    'O08': 'the same refusal with the same message is the first line of SesaHeaderReader.Parse (the file is mapped with its 4,095 bytes and refused there); the line of SesaFile.Open only keeps CreateFromFile away from a file of length 0',
    'W09': 'the hardware CRC and the table CRC give the same value (tested against each other and against the bitwise definition in 1b.7e on every section size, and the table path is forced in TheTablePath_OfTheChecksum_GivesTheSameVerdict): only the speed differs, which no assertion can see',
    'W32': 'the bit n of the bitmap is never asked: every entry at a rank above 0 is below n (the rule on the line before), so no position can reach the bit, and marking it changes nothing',
    'W37': 'LCP[1] = 0 is checked by the rule on LCP[0] and LCP[1] three lines above, and r = 1 is the only rank at which the previous entry is the empty suffix',
    'W39': 'two different suffixes are never equal: after the duplicate check and the rule l <= n - max(previous, b), the comparison from the stored LCP ends at a mismatch or at the end of the shorter suffix, never at 0',
    'W46': 'rank 0 holds n, which lies outside every chunk of text positions (p - c0 >= span), so starting the scan at rank 0 adds no entry',
    'W49': 'the common prefix with the empty suffix is 0, the value that the special case assigns; and the branch with a carry above 64 cannot meet f = n in a file that passed the permutation and order check (a predecessor of the empty suffix shares nothing, so no carry above 0 reaches it)',
    'W50': 'reached only for a suffix array that passed the permutation and the order check of the same pass; a permutation whose adjacent suffixes are in strictly ascending order is the suffix array, for which the carry of Kasai\'s scan always holds: the check cannot fail',
    'W66': 'the child-table scan runs on an LCP array that the order check and the Phi pass have shown to be the true one, and the scan of a true LCP array defines every cell exactly once',
    'W73': 'Spill moves exactly _half entries, so the number of spilled entries is always a multiple of _half and Math.Min(_half, _spilled) is _half',
    'W79': 'the LCP value of every rank was checked at 0 .. the length of the shorter suffix by the order check before the child-table scan reads it',
    'R02': 'the checks of the frame (b above the right bound, a next boundary that does not increase or leaves the interval) repeat the checks that CountChildren made on the same chain of links, and CheckFirst on its first element, when the frame was created: with those in place the frame cannot hold such a boundary',
    'R03': 'the same chain of l-indices was walked and bounded by CountChildren before the frame was created (it throws the same exception with the same message for a next boundary that does not increase)',
    'Y02': 'SesaFile.Open, which the line after it calls, begins with the same ArgumentException.ThrowIfNullOrEmpty(path): the same exception type, the same parameter name',
    'B118': 'BUILD_INFO is a constant text of about 300 bytes with a few numbers: it cannot reach the reserved 4,096 bytes, so the comparison at equality is never evaluated (survivor of the 1b.7e audit)',
}
DEFENCE = {
    'R07': 'a first l-index at the left bound of its interval is stopped one line later by the rule on the depth of the child (the lcp at the left bound is not deeper than its parent) and, if the lcp is damaged as well, by the next iteration (b < a): the same exception with the same message; only what the visitor was told before the exception differs',
    'R08': 'a first l-index beyond the right bound is stopped by the rule on the frame (b above the right bound) one iteration later: the same exception with the same message; only the events delivered before it differ',
    'N02': 'the first test of the scatter (position >= count) is what keeps the read of the inverse inside its buffer: with `>` a position equal to the count reads one entry past the end and is refused only because that entry is not -1, which is whatever lies there (in the checked build, -p:EsaBounds=true, the Buf bound refuses the read with another exception, so the equality is observable there only); no test can choose the byte behind a pinned array',
    'N05': 'a link of a node outside the text is refused one call later by Isa, which has its own bound: the same damaged read, another message',
    'F02': 'the failure path of FromFile is taken only for a text beyond the 32-bit surface (EnsureTextLength, 2^31 - 2 symbols or more): a container of that size is 20 GB and no test makes one; the disposal that the mutant leaves out is a release of the mapping, which only a mapping count can see',
    'N07': 'the strict side of the bound (a target whose parent is deeper than itself) is what the hostile rounds exercise; the equality side needs one LCP value set to exactly the depth of the target at a boundary of the target\'s interval and the child table left alone, and no such construction was made for the pyramid\'s Borders',
    'N08': 'the same for the leaf link: the equality side (the target leaf as deep as the longer of its neighbours) needs a damaged LCP array at one rank; not constructed',
    'B102': 'a guard only poisons the arena of a pass and changes no output of a correct pass; the build tests run with ArenaOptions.Guarded for the whole build, so the poison of every other pass is exercised and the BWT pass\'s is not observable (survivor of the 1b.7e audit)',
}
