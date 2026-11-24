# Mutants of the audit of S2.3 that no test kills, argued (id -> argument). Quoted in PROGRESS-S2.md (S2 2026-11-24 work package S2.3).
# EQUIVALENT: the mutant behaves like the original on every input that a caller or a test can reach today (the argument says why). Each survived stage 1 (the owning class), stage 2 (the default-set classes of the package)
# and stage 3 (the Large generator tests), in round 1 and in round 2 (after the tests that close the other survivors).
EQUIVALENT = {
    's23_S14': 'SaView.CopyTo(long start, Span<int>), the pointer copy `new ReadOnlySpan<int>(_pointer + start, ...)` with `(int)start`: the copy through the pointer is taken only by an array-free section that is not marked FourBytesWide, i.e. one of at most int.MaxValue entries, '
               'and `start <= Length - destination.Length` (the range test above it) makes start at most int.MaxValue, so the cast loses nothing. A section of more entries copies entry by entry through the checked bridge (the test '
               'TheSuffixArrayView_OfMoreThan2To31FourByteEntries_* copies at a start above 2^31 and kills the narrowing of that path: S31).',
    's23_P04': 'LcpPyramid.NextLess(long pos, long bound), `pos >= _longLength` -> `pos > _longLength`: the two differ only at pos = Length. There the original answers Length; the mutant calls the core, and the first line of both cores is the same '
               'test: NextLessWide starts with `if (pos >= _longLength) return _longLength;`, the 32-bit core with `if (pos >= _length) return _length;` (_length = _longLength for a pyramid that the 32-bit core serves). The same number.',
    's23_P17': 'LcpPyramid.Borders(long a, long b, ...) of a pyramid that the 32-bit cores serve, `b >= _longLength` -> `b > _longLength`: the two differ only at b = Length. The original throws ArgumentOutOfRangeException(nameof(a)); the mutant passes '
               'b to the 32-bit Borders, whose first test is `(uint)a >= (uint)_length || b < a || b >= _length` and throws the identical exception (same type, same parameter name, same message). The long test is needed for what the cast would '
               'wrap (b above int.MaxValue), and that case is covered (Borders(5, Length) and above Length are refused, tested).',
    's23_L07': 'LcpStore.Find, `if (v < rank) lo = mid + 1;` -> `if (v <= rank) lo = mid + 1;`: the line above is `if (v == rank) return mid;`, so the comparison is never reached with v == rank; for v != rank `<` and `<=` are the same.',
    's23_F02': 'SesaFile, `bool narrow = n < int.MaxValue && posWidth == 4 && ...` without `posWidth == 4`: the header reader refuses POS_WIDTH 5 until S2.7 (its acceptance rule is unchanged by S2.3: width 5 stays refused), so every container that SesaFile opens '
               'has posWidth 4 here and the conjunct is always true. It is a mutant that S2.7 must kill: the first container of width 5 that is opened with n < int.MaxValue takes the legacy 4-byte readers with this mutant (the tests that need it, '
               'LargeRedTests.Width5_*, are the red tests of S2.7).',
    's23_F05': 'SesaFile, the wide LCP array opened with the width 4 instead of POS_WIDTH: the same argument as F02. A container that is opened has POS_WIDTH 4 until S2.7, so `posWidth` and the constant 4 are the same number; S2.7 opens width 5 and '
               'must kill it (the overflow pairs of width 5 are 10 bytes).',
}
