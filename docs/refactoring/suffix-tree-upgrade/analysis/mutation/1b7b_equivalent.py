# Mutants of the audit of 1b.7b that are argued equivalent (id -> argument). Read by make_1b7b_round2.py (they are reported, not re-run), summarize_1b7b.py and quoted in the ADR subsection.
EQUIVALENT = {
    'O02': 'Auto is not resolved: MappedOutputFile.AutoMode is Preallocate, and an unresolved Auto falls to the same branch as Preallocate in Allocate (anything but Sparse); the two differ only if AutoMode is changed to Sparse, which O01 does and the tests kill',
    'O15': 'fsync always skipped (replaced by the compiler-rejected return): the only difference is durability after a crash, which no assertion on a correct file system can observe (same class as O16 and P35, and as E06 of 1b.7a)',
    'O16': 'fsync always done: the same bytes, only the time differs',
    'O34': 'CanPreallocate is true everywhere: on Linux, where the tests run, OperatingSystem.IsLinux() is true as well; the difference exists only on macOS and Windows, which were not run (the fake source tests the false branch)',
    'R04': 'the setter of UseVectors without the Avx2.IsSupported guard: on this machine (AVX2) the guarded and the unguarded value are equal; the difference exists only on hardware without AVX2, which cannot be tested here (the test pins UseVectors == Avx2.IsSupported and would fail there)',
    'R25': 'Lce: a difference in the partial symbol at the end of the 256-bit window returns k + vs instead of continuing: matched == vs means the first differing bit lies in the symbol with index vs, which therefore differs, so the common prefix is exactly k + vs either way (the continuing path finds it in the next step); only the step count differs',
    'R26': 'Lce: after a difference in the partial symbol the vector path restarts one symbol early (k += vs - 1): the symbol is compared again, the result is the same, one step more is possible',
    'R27': 'Lce: the same argument for the scalar path (a difference in the last, partial bit of the 64-bit window): the symbol at index ks differs, the result is k + ks either way',
    'R31': 'Lce: the vector result without Math.Min(., max): the loop runs only while k + vs <= max and matched < vs, so k + matched < max already',
    'F12': 'PackInto emits a byte when more than 8 bits are pending instead of at 8: the pending bits are written by the final partial-byte statement (acc << (8 - accBits) with accBits = 8 is acc << 0), the bytes and their order are the same',
    'P26': 'Pack without the data-length check: if every symbol count equals the scanned one, the number of symbols (and so of data bytes) is equal; the check is defence in depth after the count comparison',
    'P35': 'WriteFile always without fsync: durability only (the same class as O16)',
}
