#!/usr/bin/env python3
"""Work package 1b.9e.1: single-site mutants of the lines that the package adds or changes: the normalisation of the reopened container (ZeroFill.ClearOutside and its helpers, the keep ranges and the result of
MappedOutputFile.OpenExisting), the ranges that the kept stages vouch for and the report (SesaResumeRun), and the call in SesaBuilder.Build. Written after the tests. Each is built and run by run_1b9c_mutants.py (Release,
watchdog, analyzers off, no NUnit.StopOnError) with stage 1 = the tests of the package (ZeroFillClearingTests for the file level; ZeroFillClearingTests, SesaResume1b9e1Tests, SesaResumeFormatTests for the rest), stage 1b = the
survivors of the quick set against the whole set, closed = the survivors run again at the head that holds the tests written for them, stage 2 = the whole default set. Writes 1b9e1_<file>_listed.json. A pattern that does not
match exactly once is reported and must be fixed before the list is run."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '../../../../..'))
E = 'src/SuffixTree/Algorithms/SuffixTree.Esa/External/'
NEW = 'FullyQualifiedName~ZeroFillClearingTests|FullyQualifiedName~SesaResume1b9e1Tests|FullyQualifiedName~SesaResumeFormatTests'
RS, MO, SB = 'Sesa/SesaResume.cs', 'MappedOutputFile.cs', 'Sesa/SesaBuilder.cs'
M = [
    # file, id, description, old, new
    # ---- ZeroFill.ClearOutside
    (MO, 'Z01', 'the ranges are not sorted', 'keep.OrderBy(k => k.Offset).Append((length, 0))', 'keep.Append((length, 0))'),
    (MO, 'Z02', 'the last gap (up to the end of the file) is not cleared', 'keep.OrderBy(k => k.Offset).Append((length, 0))', 'keep.OrderBy(k => k.Offset)'),
    (MO, 'Z03', 'a range that starts beyond the file does not end the gap at the file', 'long gapEnd = Math.Min(offset, length);', 'long gapEnd = offset;'),
    (MO, 'Z04', 'a gap of no bytes is also visited', 'if (gapEnd > at)', 'if (gapEnd >= at)'),
    (MO, 'Z05', 'the end of a range does not move the cursor over a nested range (the last range wins)', 'at = Math.Max(at, offset + len);', 'at = offset + len;'),
    (MO, 'Z06', 'the cursor is not moved past a kept range', 'at = Math.Max(at, offset + len);', 'at = Math.Max(at, offset);'),
    (MO, 'Z07', 'the cursor is moved one byte past a kept range', 'at = Math.Max(at, offset + len);', 'at = Math.Max(at, offset + len + 1);'),
    (MO, 'Z08', 'the result swaps its two numbers', 'return new ClearResult(scanned, cleared);', 'return new ClearResult(cleared, scanned);'),
    (MO, 'Z09', 'the lengths of the file and of the ranges are not checked for sign (a negative length)', 'ArgumentOutOfRangeException.ThrowIfNegative(length);\n        long scanned = 0, cleared = 0;', 'long scanned = 0, cleared = 0;'),
    (MO, 'Z10', 'a null list of ranges is not refused', 'ArgumentNullException.ThrowIfNull(keep);\n        ArgumentOutOfRangeException.ThrowIfNegative(length);\n        long scanned', 'ArgumentOutOfRangeException.ThrowIfNegative(length);\n        long scanned'),
    (MO, 'Z11', 'a null handle is not refused', 'ArgumentNullException.ThrowIfNull(file);\n        ArgumentNullException.ThrowIfNull(keep);', 'ArgumentNullException.ThrowIfNull(keep);'),
    (MO, 'Z12', 'the chunk is one byte shorter than the step', 'for (long chunk = start; chunk < end; chunk += ChunkBytes)\n        {\n            int take', 'for (long chunk = start; chunk < end; chunk += ChunkBytes - 1)\n        {\n            int take'),
    (MO, 'Z13', 'a chunk is not cut at the end of the gap', 'int take = (int)Math.Min(ChunkBytes, end - chunk);\n            for (int got', 'int take = ChunkBytes;\n            for (int got'),
    (MO, 'Z14', 'a chunk is one byte shorter', 'int take = (int)Math.Min(ChunkBytes, end - chunk);\n            for (int got', 'int take = (int)Math.Min(ChunkBytes, end - chunk - 1);\n            for (int got'),
    (MO, 'Z15', 'a read of no bytes is not an error (the loop does not end)', 'if (n <= 0) throw new IOException', 'if (n < 0) throw new IOException'),
    (MO, 'Z16', 'a short read is taken as a full one', 'got += n;', 'got = take;'),
    (MO, 'Z17', 'a short read is read again at the same place', 'int n = RandomAccess.Read(file, buffer.AsSpan(got, take - got), chunk + got);', 'int n = RandomAccess.Read(file, buffer.AsSpan(got, take - got), chunk);'),
    (MO, 'Z18', 'a short read is read into the start of the buffer', 'int n = RandomAccess.Read(file, buffer.AsSpan(got, take - got), chunk + got);', 'int n = RandomAccess.Read(file, buffer.AsSpan(0, take - got), chunk + got);'),
    (MO, 'Z19', 'the scanned bytes are one short per chunk', 'scanned += take;', 'scanned += take - 1;'),
    (MO, 'Z20', 'the scanned bytes are not counted', 'scanned += take;', ''),
    (MO, 'Z21', 'the last unit of a chunk is not looked at', 'for (int unit = 0; unit < take; unit += UnitBytes)', 'for (int unit = 0; unit < take - 1; unit += UnitBytes)'),
    (MO, 'Z22', 'the loop over the units has one more step', 'for (int unit = 0; unit < take; unit += UnitBytes)', 'for (int unit = 0; unit <= take; unit += UnitBytes)'),
    (MO, 'Z23', 'a unit is looked at in full also at the end of the chunk (the stale tail of the buffer)', 'if (buffer.AsSpan(unit, Math.Min(UnitBytes, take - unit)).IndexOfAnyExcept((byte)0) >= 0)', 'if (buffer.AsSpan(unit, UnitBytes).IndexOfAnyExcept((byte)0) >= 0)'),
    (MO, 'Z24', 'a unit of zeros is taken for data', 'if (buffer.AsSpan(unit, Math.Min(UnitBytes, take - unit)).IndexOfAnyExcept((byte)0) >= 0)', 'if (buffer.AsSpan(unit, Math.Min(UnitBytes, take - unit)).IndexOfAnyExcept((byte)0) >= -1)'),
    (MO, 'Z25', 'only the first non-zero unit of a run opens it: a run starts at the last unit', 'if (dirtyFrom < 0) dirtyFrom = unit;', 'dirtyFrom = unit;'),
    (MO, 'Z26', 'a run of non-zero units is never started', 'if (dirtyFrom < 0) dirtyFrom = unit;', ''),
    (MO, 'Z27', 'a clean unit writes even when no run is open', 'else if (dirtyFrom >= 0)', 'else'),
    (MO, 'Z28', 'a run is written one unit late', 'cleared += WriteZeros(file, chunk + dirtyFrom, unit - dirtyFrom);', 'cleared += WriteZeros(file, chunk + dirtyFrom, unit - dirtyFrom + UnitBytes);'),
    (MO, 'Z29', 'a run is written from the wrong place', 'cleared += WriteZeros(file, chunk + dirtyFrom, unit - dirtyFrom);', 'cleared += WriteZeros(file, chunk + dirtyFrom + 1, unit - dirtyFrom - 1);'),
    (MO, 'Z30', 'a run is written from the start of the chunk', 'cleared += WriteZeros(file, chunk + dirtyFrom, unit - dirtyFrom);', 'cleared += WriteZeros(file, chunk, unit - dirtyFrom);'),
    (MO, 'Z31', 'a run is not closed after it is written', 'dirtyFrom = -1;\n                }\n            }', '}\n            }'),
    (MO, 'Z32', 'the run that is open at the end of the chunk is not written', 'if (dirtyFrom >= 0) cleared += WriteZeros(file, chunk + dirtyFrom, take - dirtyFrom);', ''),
    (MO, 'Z33', 'the run at the end of the chunk is written one byte short', 'if (dirtyFrom >= 0) cleared += WriteZeros(file, chunk + dirtyFrom, take - dirtyFrom);', 'if (dirtyFrom >= 0) cleared += WriteZeros(file, chunk + dirtyFrom, take - dirtyFrom - 1);'),
    (MO, 'Z34', 'the run at the end of the chunk is written a whole unit', 'if (dirtyFrom >= 0) cleared += WriteZeros(file, chunk + dirtyFrom, take - dirtyFrom);', 'if (dirtyFrom >= 0) cleared += WriteZeros(file, chunk + dirtyFrom, UnitBytes);'),
    (MO, 'Z35', 'the zeros are written one byte short', 'RandomAccess.Write(file, s_zeros.AsSpan(0, length), offset);\n        return length;', 'RandomAccess.Write(file, s_zeros.AsSpan(0, length - 1), offset);\n        return length;'),
    (MO, 'Z36', 'the zeros are written one byte late', 'RandomAccess.Write(file, s_zeros.AsSpan(0, length), offset);\n        return length;', 'RandomAccess.Write(file, s_zeros.AsSpan(0, length), offset + 1);\n        return length;'),
    (MO, 'Z37', 'the cleared bytes are not counted', 'RandomAccess.Write(file, s_zeros.AsSpan(0, length), offset);\n        return length;', 'RandomAccess.Write(file, s_zeros.AsSpan(0, length), offset);\n        return 0;'),
    # ---- MappedOutputFile.OpenExisting
    (MO, 'O01', 'the clearing is done only together with a cut', 'if (truncateTo is not null || keep is not null)', 'if (truncateTo is not null)'),
    (MO, 'O02', 'the keep ranges are not applied', 'if (keep is not null) cleared = ZeroFill.ClearOutside(handle, RandomAccess.GetLength(handle), keep);', ''),
    (MO, 'O03', 'a file that has no keep ranges is cleared whole', 'if (keep is not null) cleared = ZeroFill.ClearOutside(handle, RandomAccess.GetLength(handle), keep);', 'cleared = ZeroFill.ClearOutside(handle, RandomAccess.GetLength(handle), keep ?? []);'),
    (MO, 'O04', 'the file is cleared before it is cut', 'if (truncateTo is { } cut && RandomAccess.GetLength(handle) > cut) RandomAccess.SetLength(handle, cut);\n            if (keep is not null) cleared = ZeroFill.ClearOutside(handle, RandomAccess.GetLength(handle), keep);', 'if (keep is not null) cleared = ZeroFill.ClearOutside(handle, RandomAccess.GetLength(handle), keep);\n            if (truncateTo is { } cut && RandomAccess.GetLength(handle) > cut) RandomAccess.SetLength(handle, cut);'),
    (MO, 'O05', 'the length that is cleared is the cut, not the length of the file', 'cleared = ZeroFill.ClearOutside(handle, RandomAccess.GetLength(handle), keep);', 'cleared = ZeroFill.ClearOutside(handle, truncateTo ?? RandomAccess.GetLength(handle) + 1, keep);'),
    (MO, 'O06', 'the result of the clearing is not kept', '{ Cleared = cleared }', '{ }'),
    (MO, 'O07', 'a file is cut when it is not longer (it grows)', 'if (truncateTo is { } cut && RandomAccess.GetLength(handle) > cut) RandomAccess.SetLength(handle, cut);', 'if (truncateTo is { } cut) RandomAccess.SetLength(handle, cut);'),
    (MO, 'O08', 'a file is cut when it is exactly as long as the cut', 'if (truncateTo is { } cut && RandomAccess.GetLength(handle) > cut) RandomAccess.SetLength(handle, cut);', 'if (truncateTo is { } cut && RandomAccess.GetLength(handle) >= cut) RandomAccess.SetLength(handle, cut - 1);'),
    # ---- the ranges that the kept stages vouch for, and the report
    (RS, 'R01', 'the ranges are those of every stage that the state names', 'internal IReadOnlyList<(long Offset, long Length)> KeptContainerRanges() => [.. _skipped.SelectMany(', 'internal IReadOnlyList<(long Offset, long Length)> KeptContainerRanges() => [.. _state.Stages.SelectMany('),
    (RS, 'R02', 'the scratch files of the kept stages are taken for ranges of the container', '.Where(c => c.InContainer).Select(c => (c.Offset, c.Length))];', '.Select(c => (c.Offset, c.Length))];'),
    (RS, 'R03', 'only the ranges of scratch files are kept', '.Where(c => c.InContainer).Select(c => (c.Offset, c.Length))];', '.Where(c => !c.InContainer).Select(c => (c.Offset, c.Length))];'),
    (RS, 'R04', 'a kept range is one byte shorter', '.Where(c => c.InContainer).Select(c => (c.Offset, c.Length))];', '.Where(c => c.InContainer).Select(c => (c.Offset, c.Length - 1))];'),
    (RS, 'R05', 'a kept range starts one byte late', '.Where(c => c.InContainer).Select(c => (c.Offset, c.Length))];', '.Where(c => c.InContainer).Select(c => (c.Offset + 1, c.Length))];'),
    (RS, 'R06', 'the result of the normalisation is not noted', '        _cleared = cleared;\n        _clearSeconds = seconds;', '        _clearSeconds = seconds;'),
    (RS, 'R07', 'the seconds of the normalisation are not noted', '        _cleared = cleared;\n        _clearSeconds = seconds;', '        _cleared = cleared;'),
    (RS, 'R08', 'the report swaps the two numbers', 'ClearScannedBytes = _cleared.ScannedBytes,\n        ClearedBytes = _cleared.ClearedBytes,', 'ClearScannedBytes = _cleared.ClearedBytes,\n        ClearedBytes = _cleared.ScannedBytes,'),
    (RS, 'R09', 'the report has no seconds of the normalisation', 'ClearSeconds = _clearSeconds,', 'ClearSeconds = 0,'),
    (RS, 'R10', 'the report has no scanned bytes', 'ClearScannedBytes = _cleared.ScannedBytes,', 'ClearScannedBytes = 0,'),
    (RS, 'R11', 'the report has no cleared bytes', 'ClearedBytes = _cleared.ClearedBytes,', 'ClearedBytes = 0,'),
    # ---- SesaBuilder.Build
    (SB, 'B01', 'the keep ranges are not handed to the opening of the container', 'container = MappedOutputFile.OpenExisting(tmp, cut, run.KeptContainerRanges());', 'container = MappedOutputFile.OpenExisting(tmp, cut);'),
    (SB, 'B02', 'the keep ranges are empty (every byte is cleared)', 'container = MappedOutputFile.OpenExisting(tmp, cut, run.KeptContainerRanges());', 'container = MappedOutputFile.OpenExisting(tmp, cut, []);'),
    (SB, 'B03', 'the result of the normalisation is not handed to the run', '                run.NoteCleared(container.Cleared, Elapsed(tClear));\n', ''),
    (SB, 'B04', 'the seconds of the normalisation are zero', 'run.NoteCleared(container.Cleared, Elapsed(tClear));', 'run.NoteCleared(container.Cleared, 0);'),
    (SB, 'B05', 'the seconds of the normalisation include the whole build so far', 'run.NoteCleared(container.Cleared, Elapsed(tClear));', 'run.NoteCleared(container.Cleared, Elapsed(tStart));'),
]
by_file = {}
bad = 0
for f, i, desc, old, new in M:
    src = open(os.path.join(ROOT, E + f)).read()
    if src.count(old) != 1:
        print('PATTERN', i, f, src.count(old))
        bad += 1
    by_file.setdefault(f, []).append({'id': i, 'desc': desc, 'file': E + f, 'old': old, 'new': new, 'filter': '(%s)' % NEW})
for f, ms in by_file.items():
    json.dump(ms, open(os.path.join(HERE, '1b9e1_%s_listed.json' % os.path.splitext(os.path.basename(f))[0].lower()), 'w'), indent=1)
print(len(M), 'mutants in', len(by_file), 'files;', bad, 'patterns to fix')
