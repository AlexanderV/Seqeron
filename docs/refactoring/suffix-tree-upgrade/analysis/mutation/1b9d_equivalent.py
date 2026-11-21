#!/usr/bin/env python3
"""Work package 1b.9d: the mutants that survived every test of the default set (and the tests written for the survivors) and are argued equivalent, one by one. Each argument was also checked by the whole default set of
SuffixTree.Esa.Tests (stage 2 of run_1b9d_audit.py): an equivalent mutant must survive it. summarize_1b9d.py reads EQUIVALENT."""
EQUIVALENT = {
    'Z03': ('ZeroFill: the loop condition at < end becomes at <= end',
            'The only difference is one more iteration with at == end, which happens when the length is zero or a multiple of the chunk: take = min(chunk, end - at) = 0, so the call is a write of an empty span '
            'and the counter gets 0. A write of no bytes does not change the file (pwrite of zero bytes does not extend it: the test ZeroFill_WritesExactlyTheRange... writes zero bytes at an offset 1 MB behind the end '
            'and asserts that the length stays 3,000,000; the mutant passes it too). The mutant survived the whole default set.'),
    'R04': ('MappedOutputFile.Grow: the lock of the mapping is not taken',
            'The lock keeps OpenView (another thread) from using the mapping while Grow replaces it. Grow requires that no writer is open (it throws InvalidOperationException otherwise, tested), and a writer is the only '
            'thing that calls OpenView, so under the contract of the type there is no other thread in the mapping while Grow runs; the lock is a defence against a misuse (CreateWriter on another thread during Grow) that '
            'no test of the project does and that the type does not promise to support. Not observable by a deterministic test; the mutant survived the whole default set.'),
}
