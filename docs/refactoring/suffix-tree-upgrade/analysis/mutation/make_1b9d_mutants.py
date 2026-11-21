#!/usr/bin/env python3
"""Work package 1b.9d: single-site mutants of the code that the package adds or changes (the reservation of the growth of a file, the zero fill, the default of the interface, the counter of the reserved bytes, the target of
Dist, the call of the builder), written after the tests. Each is built and run by run_1b9c_mutants.py (Release, watchdog, analyzers off, no NUnit.StopOnError) with stage 1 = the tests of the package (the classes whose name
contains 1b9dTests plus the classes that name the mutated type, found by the runner) and, for a survivor, stage 2 = the whole default set. Writes 1b9d_<file>_listed.json. A pattern that does not match exactly once is reported."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '../../../../..'))
E = 'src/SuffixTree/Algorithms/SuffixTree.Esa/External/'
NEW = 'FullyQualifiedName~1b9dTests'
M = [
    # file, id, description, old, new
    ('MappedOutputFile.cs', 'G01', 'the growth is never reserved', 'bool reserve = (options is null || options.Reservation == ReservationMode.Auto ? AutoMode : options.Reservation) != ReservationMode.Sparse;', 'bool reserve = false;'),
    ('MappedOutputFile.cs', 'G02', 'the growth is reserved only in the sparse mode', 'options.Reservation) != ReservationMode.Sparse;', 'options.Reservation) == ReservationMode.Sparse;'),
    ('MappedOutputFile.cs', 'G03', 'the automatic mode grows sparse', 'options.Reservation == ReservationMode.Auto ? AutoMode : options.Reservation)', 'options.Reservation == ReservationMode.Auto ? ReservationMode.Sparse : options.Reservation)'),
    ('MappedOutputFile.cs', 'G04', 'no options mean the sparse mode', '(options is null || options.Reservation == ReservationMode.Auto ?', '(options is null ? ReservationMode.Sparse : options.Reservation == ReservationMode.Auto ?'),
    ('MappedOutputFile.cs', 'G05', 'the reserved range starts at the start of the file', 'source.ReserveRange(handle, previous, growth, counters);', 'source.ReserveRange(handle, 0, newLength, counters);'),
    ('MappedOutputFile.cs', 'G06', 'the reserved range is as long as the new file', 'source.ReserveRange(handle, previous, growth, counters);', 'source.ReserveRange(handle, previous, newLength, counters);'),
    ('MappedOutputFile.cs', 'G07', 'the reserved range starts one byte late', 'source.ReserveRange(handle, previous, growth, counters);', 'source.ReserveRange(handle, previous + 1, growth - 1, counters);'),
    ('MappedOutputFile.cs', 'G08', 'the reserved range ends one byte early', 'source.ReserveRange(handle, previous, growth, counters);', 'source.ReserveRange(handle, previous, growth - 1, counters);'),
    ('MappedOutputFile.cs', 'G09', 'the counters are not given to the source', 'source.ReserveRange(handle, previous, growth, counters);', 'source.ReserveRange(handle, previous, growth, null);'),
    ('MappedOutputFile.cs', 'G10', 'a failed reservation does not give the length back', '                        RandomAccess.SetLength(handle, previous);\n                        throw;', '                        throw;'),
    ('MappedOutputFile.cs', 'G11', 'a failed reservation gives back the new length', 'RandomAccess.SetLength(handle, previous);\n                        throw;', 'RandomAccess.SetLength(handle, newLength);\n                        throw;'),
    ('MappedOutputFile.cs', 'G12', 'a failed reservation is a different exception', '                        RandomAccess.SetLength(handle, previous);\n                        throw;', '                        RandomAccess.SetLength(handle, previous);\n                        throw new InvalidOperationException("the reservation failed");'),
    ('MappedOutputFile.cs', 'G13', 'the check of the free space refuses exactly enough', 'if (free >= 0 && free < growth) throw new IOException($"not enough free space to grow', 'if (free >= 0 && free <= growth) throw new IOException($"not enough free space to grow'),
    ('MappedOutputFile.cs', 'G14', 'the check of the free space of the growth is not made', 'if (free >= 0 && free < growth) throw new IOException($"not enough free space to grow', 'if (free >= 0 && free < 0) throw new IOException($"not enough free space to grow'),
    ('MappedOutputFile.cs', 'G15', 'the file is one byte short after the growth', 'RandomAccess.SetLength(handle, newLength);', 'RandomAccess.SetLength(handle, newLength - 1);'),
    ('MappedOutputFile.cs', 'G16', 'a reserved growth says Sparse', 'if (reserve) return ReservationOutcome.Filled;', 'if (reserve) return ReservationOutcome.Sparse;'),
    ('MappedOutputFile.cs', 'G17', 'a reserved growth says what the check did', 'if (reserve) return ReservationOutcome.Filled;', 'if (reserve) return free >= 0 ? ReservationOutcome.FreeSpaceChecked : ReservationOutcome.Unchecked;'),
    ('MappedOutputFile.cs', 'G18', 'the sparse growth with a known free space says Unchecked', 'return free >= 0 ? ReservationOutcome.FreeSpaceChecked : ReservationOutcome.Unchecked;\n    }\n\n    /// <summary>A writer for the bytes', 'return free >= 0 ? ReservationOutcome.Unchecked : ReservationOutcome.FreeSpaceChecked;\n    }\n\n    /// <summary>A writer for the bytes'),
    ('MappedOutputFile.cs', 'G19', 'the file is opened without sharing', 'File.OpenHandle(Path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite))\n            {\n                RandomAccess.SetLength(handle, newLength);', 'File.OpenHandle(Path, FileMode.Open, FileAccess.ReadWrite, FileShare.None))\n            {\n                RandomAccess.SetLength(handle, newLength);'),
    ('MappedOutputFile.cs', 'G20', 'the file is not extended before the reservation', 'RandomAccess.SetLength(handle, newLength);\n                if (reserve)', 'if (reserve)'),
    ('MappedOutputFile.cs', 'G21', 'the length of the object is not updated', '            _map = map;\n            Length = newLength;', '            _map = map;'),
    ('MappedOutputFile.cs', 'G22', 'a growth of nothing reserves', 'if (newLength == Length) return ReservationOutcome.Sparse;', 'if (newLength == Length) return ReservationOutcome.Filled;'),
    ('MappedOutputFile.cs', 'Z01', 'the chunk of the zero fill is 512 KiB', 'internal const int ChunkBytes = 1 << 20;', 'internal const int ChunkBytes = 1 << 19;'),
    ('MappedOutputFile.cs', 'Z02', 'the zero fill writes a chunk that is a byte short of the next', 'for (long at = offset; at < end; at += ChunkBytes)', 'for (long at = offset; at < end; at += ChunkBytes - 1)'),
    ('MappedOutputFile.cs', 'Z03', 'the zero fill also writes at the end of the range', 'for (long at = offset; at < end; at += ChunkBytes)', 'for (long at = offset; at <= end; at += ChunkBytes)'),
    ('MappedOutputFile.cs', 'Z04', 'the last chunk of the zero fill is a full chunk', 'int take = (int)Math.Min(ChunkBytes, end - at);', 'int take = ChunkBytes;'),
    ('MappedOutputFile.cs', 'Z05', 'the zero fill counts one byte per write', 'counters?.AddReserved(take);', 'counters?.AddReserved(1);'),
    ('MappedOutputFile.cs', 'Z06', 'the zero fill does not count', '            counters?.AddReserved(take);\n', ''),
    ('MappedOutputFile.cs', 'Z07', 'the range ends one byte late', 'long end = offset + length;\n        for (long at = offset;', 'long end = offset + length + 1;\n        for (long at = offset;'),
    ('MappedOutputFile.cs', 'Z08', 'the zeros are ones', 'private static readonly byte[] s_zeros = new byte[ChunkBytes];', 'private static readonly byte[] s_zeros = [.. Enumerable.Repeat((byte)1, ChunkBytes)];'),
    ('MappedOutputFile.cs', 'Z09', 'a negative offset is accepted', '        ArgumentOutOfRangeException.ThrowIfNegative(offset);\n        ArgumentOutOfRangeException.ThrowIfNegative(length);\n        long end', '        ArgumentOutOfRangeException.ThrowIfNegative(length);\n        long end'),
    ('MappedOutputFile.cs', 'Z10', 'a negative length is accepted', '        ArgumentOutOfRangeException.ThrowIfNegative(length);\n        long end = offset + length;', '        long end = offset + length;'),
    ('MappedOutputFile.cs', 'Z11', 'a null handle is accepted', '        ArgumentNullException.ThrowIfNull(file);\n        ArgumentOutOfRangeException.ThrowIfNegative(offset);', '        ArgumentOutOfRangeException.ThrowIfNegative(offset);'),
    ('MappedOutputFile.cs', 'Z12', 'the zeros are written from the start of the file', 'RandomAccess.Write(file, s_zeros.AsSpan(0, take), at);', 'RandomAccess.Write(file, s_zeros.AsSpan(0, take), at - offset);'),
    ('MappedOutputFile.cs', 'D01', 'the default of the interface reserves nothing', 'void ReserveRange(SafeFileHandle file, long offset, long length, IoCounters? counters) => ZeroFill.Write(file, offset, length, counters);', 'void ReserveRange(SafeFileHandle file, long offset, long length, IoCounters? counters) { }'),
    ('MappedOutputFile.cs', 'D02', 'the default of the interface reserves half', 'IoCounters? counters) => ZeroFill.Write(file, offset, length, counters);', 'IoCounters? counters) => ZeroFill.Write(file, offset, length / 2, counters);'),
    ('MappedOutputFile.cs', 'D03', 'the default of the interface does not count', 'IoCounters? counters) => ZeroFill.Write(file, offset, length, counters);', 'IoCounters? counters) => ZeroFill.Write(file, offset, length, null);'),
    ('IoCounters.cs', 'I01', 'the reserved bytes are counted one more', 'internal void AddReserved(long bytes) => Interlocked.Add(ref _reserved, bytes);', 'internal void AddReserved(long bytes) => Interlocked.Add(ref _reserved, bytes + 1);'),
    ('IoCounters.cs', 'I02', 'Reset keeps the reserved bytes', '        Interlocked.Exchange(ref _reserved, 0);\n', ''),
    ('IoCounters.cs', 'I03', 'the reserved bytes are the explicit writes', 'internal long ReservedWriteBytes => Interlocked.Read(ref _reserved);', 'internal long ReservedWriteBytes => Interlocked.Read(ref _writeBytes);'),
    ('IoCounters.cs', 'I04', 'a reserved byte is also an explicit write', 'internal void AddReserved(long bytes) => Interlocked.Add(ref _reserved, bytes);', 'internal void AddReserved(long bytes) { Interlocked.Add(ref _reserved, bytes); Interlocked.Add(ref _writeBytes, bytes); }'),
    ('Dist.cs', 'S01', 'the target of Dist is created sparse', 'MappedOutputFile.Create(targetPath, slots * sizeof(TValue), options.Files);', 'MappedOutputFile.Create(targetPath, slots * sizeof(TValue), new MappedFileOptions { Reservation = ReservationMode.Sparse });'),
    ('Dist.cs', 'S02', 'the target of Dist ignores the options', 'MappedOutputFile.Create(targetPath, slots * sizeof(TValue), options.Files);', 'MappedOutputFile.Create(targetPath, slots * sizeof(TValue));'),
    ('Dist.cs', 'S03', 'Dist reports a sparse target', 'return new DistStats(chunks, bytesWritten, bytesRead, created.Outcome);', 'return new DistStats(chunks, bytesWritten, bytesRead, ReservationOutcome.Sparse);'),
    ('Sesa/SesaBuilder.cs', 'B01', 'the builder gives no counters to the growth', 'growth = container.Grow(tail.FinalLength, options.Files, counters);', 'growth = container.Grow(tail.FinalLength, options.Files);'),
    ('Sesa/SesaBuilder.cs', 'B02', 'the builder gives no file options to the growth', 'growth = container.Grow(tail.FinalLength, options.Files, counters);', 'growth = container.Grow(tail.FinalLength, null, counters);'),
    ('Sesa/SesaBuilder.cs', 'B03', 'the builder forgets what the growth did', 'growth = container.Grow(tail.FinalLength, options.Files, counters);', 'container.Grow(tail.FinalLength, options.Files, counters);'),
    # the second batch: written after the tests that the survivors of the first (Z03, Z09, Z11) and the review of the tests asked for
    ('MappedOutputFile.cs', 'R01', 'the explicit preallocate mode grows sparse', 'options.Reservation == ReservationMode.Auto ? AutoMode : options.Reservation)', 'options.Reservation == ReservationMode.Auto ? AutoMode : ReservationMode.Sparse)'),
    ('MappedOutputFile.cs', 'R02', 'only an IOException of the reservation gives the length back', '                    catch\n                    {\n                        RandomAccess.SetLength(handle, previous);', '                    catch (IOException)\n                    {\n                        RandomAccess.SetLength(handle, previous);'),
    ('MappedOutputFile.cs', 'R03', 'the zero fill counts a chunk before it is written', '            RandomAccess.Write(file, s_zeros.AsSpan(0, take), at);\n            counters?.AddReserved(take);', '            counters?.AddReserved(take);\n            RandomAccess.Write(file, s_zeros.AsSpan(0, take), at);'),
    ('MappedOutputFile.cs', 'R04', 'the growth does not take the lock of the mapping', '        long previous = Length;\n        lock (_gate)\n        {\n            using (SafeFileHandle handle', '        long previous = Length;\n        {\n            using (SafeFileHandle handle'),
    ('MappedOutputFile.cs', 'R05', 'a growth that is not longer than the file is also refused', 'ArgumentOutOfRangeException.ThrowIfLessThan(newLength, Length);\n        if (newLength == Length) return ReservationOutcome.Sparse;', 'ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(newLength, Length);\n        if (newLength == Length) return ReservationOutcome.Sparse;'),
    ('Dist.cs', 'R06', 'the target writer of Dist is not given the counters', 'created.CreateWriter(options.WindowBytes, counters);', 'created.CreateWriter(options.WindowBytes, null);'),
    ('Sesa/SesaBuilder.cs', 'R07', 'the container is grown also when it does not grow', 'if (growthBytes > 0)\n                {\n                    growth = container.Grow', 'if (growthBytes >= 0)\n                {\n                    growth = container.Grow'),
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
    json.dump(ms, open(os.path.join(HERE, '1b9d_%s_listed.json' % os.path.splitext(os.path.basename(f))[0].lower()), 'w'), indent=1)
print(len(M), 'mutants in', len(by_file), 'files;', bad, 'patterns to fix')
