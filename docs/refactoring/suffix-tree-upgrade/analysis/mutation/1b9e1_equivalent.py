#!/usr/bin/env python3
"""Work package 1b.9e.1: the mutants that survived every test of the default set (and the tests written for the survivors) and are argued equivalent, one by one. Each argument was also checked by the whole default set of
SuffixTree.Esa.Tests (stage 2 of run_1b9e1_audit.py): an equivalent mutant must survive it. summarize_1b9e1.py reads EQUIVALENT."""
EQUIVALENT = {
    'Z04': ('ZeroFill.ClearOutside: a gap of no bytes is also visited (>= instead of >)',
            'The condition guards the call of ClearSpan(file, at, gapEnd) and the lazy allocation of the buffer. With gapEnd == at the span is empty: the loop of ClearSpan, for (chunk = start; chunk < end; ...), does not run, so nothing is read, '
            'written or counted. The only difference is that the 1 MiB buffer is allocated although nothing is read, which no result, counter or test can see.'),
    'Z10': ('ZeroFill.ClearOutside: a null list of ranges is not refused by an explicit check',
            'The first use of the list is keep.OrderBy(...), and LINQ refuses a null source with an ArgumentNullException (parameter "source"): the same exception type for the same input. Only the name of the parameter in the message differs. '
            'The explicit check is kept for the clearer message (the test asserts the type).'),
    'Z18': ('ZeroFill.ClearSpan: after a short read the next read goes into the start of the buffer instead of behind the bytes already read',
            'The first read of a chunk starts at got = 0, where the two forms are the same, and a read of a regular file returns every byte that was asked for unless it meets the end of the file. At the end of the file the first short read is followed by a read '
            'that returns 0 and ends the call with an IOException in both forms; the test of the file that is shorter than its length passes with either. Only a file system that returns a short read in the middle of a file (a pipe, some network volumes) can tell them apart, '
            'and nothing in a test can make a regular file do that. The form with the offset got is the one that is right when it happens.'),
    'Z22': ('ZeroFill.ClearSpan: the loop over the units has one step more (unit <= take instead of unit < take)',
            'When take is a multiple of the unit the extra step is at unit == take: the slice is empty, IndexOfAnyExcept of an empty span is -1, the unit counts as clean, and a run that is still open is written there, with the length unit - dirtyFrom == take - dirtyFrom, '
            'which is what the line after the loop writes in the original. When take is not a multiple of the unit the steps are multiples of 4,096 and the last one that is below take is the same in both loops; the next one is above take and ends the loop in both. '
            'The bytes written, the bytes counted and the order of the writes are the same.'),
}
