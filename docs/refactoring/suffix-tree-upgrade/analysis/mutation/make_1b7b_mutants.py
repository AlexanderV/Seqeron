#!/usr/bin/env python3
"""Writes the hand-made single-site mutants of the new production code of work package 1b.7b (src/SuffixTree/Algorithms/SuffixTree.Esa/External: PackedTextFormat, MappedTextReader, PackedTextWriter,
   MappedOutputFile, and the changes to MappedWindowWriter) as one JSON list per source file (1b7b_<file>_listed.json). Each pattern must occur exactly once in its file (checked here). Round 1 is written against
   the sources at the commit of the first audit run; later rounds (1b7b_r2_*, 1b7b_r3_*) are written by make_1b7b_round2.py for the mutants that did not compile or survived, after the tests were strengthened.
   usage: make_1b7b_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1]
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/External/'


def m(i, desc, old, new):
    return {'id': i, 'desc': desc, 'old': old, 'new': new}


FILES = {
    'PackedTextFormat.cs': [
        m('F01', 'BitsFor: sigma 3 gets one bit', 'if (sigma <= 2) return 1;', 'if (sigma <= 3) return 1;'),
        m('F02', 'BitsFor: counts the bits of sigma, not of sigma - 1', '(ulong)(sigma - 1)', '(ulong)sigma'),
        m('F03', 'BitsFor: sigma 2 is not special-cased below (always the formula)', 'if (sigma <= 2) return 1;', 'if (sigma <= 0) return 1;'),
        m('F04', 'DataBytes rounds up one byte too much', '((n * bits) + 7) / 8', '((n * bits) + 8) / 8'),
        m('F05', 'DataBytes rounds down', '((n * bits) + 7) / 8', '((n * bits) + 6) / 8'),
        m('F06', 'TotalBytes forgets one padding byte', 'DataBytes(n, bits) + PadBytes;', 'DataBytes(n, bits) + PadBytes - 1;'),
        m('F07', 'PackInto: a destination that is too small is not refused', 'if (destination.Length < need) throw', 'if (destination.Length < 0) throw'),
        m('F08', 'PackInto: a symbol that does not fit is not refused', 'if ((s >> bits) != 0) throw', 'if ((s >> bits) > 255) throw'),
        m('F09', 'PackInto: the fit test misses symbols just above the range', 'if ((s >> bits) != 0) throw', 'if ((s >> bits) > 1) throw'),
        m('F10', 'PackInto: the last partial byte is not written', '        if (accBits > 0) destination[o] = (byte)(acc << (8 - accBits));', '        _ = accBits;'),
        m('F11', 'PackInto: the last partial byte is shifted the wrong way', '(byte)(acc << (8 - accBits))', '(byte)(acc >> (8 - accBits))'),
        m('F12', 'PackInto: bytes are emitted one symbol late', 'while (accBits >= 8)', 'while (accBits > 8)'),
        m('F13', 'PackInto: symbols are added instead of or-ed', 'acc = (acc << bits) | s;', 'acc = (acc << bits) + s;'),
        m('F14', 'PackInto: the upper width is not checked', '        ArgumentOutOfRangeException.ThrowIfGreaterThan(bits, MaxBits);\n        long need', '        long need'),
        m('F15', 'PackInto: the lower width is not checked', '        ArgumentOutOfRangeException.ThrowIfLessThan(bits, MinBits);\n', ''),
        m('F16', 'the padding is 63 bytes', 'internal const int PadBytes = 64;', 'internal const int PadBytes = 63;'),
        m('F17', 'the design limit is one more', 'uint.MaxValue - 3L;', 'uint.MaxValue - 2L;'),
        m('F18', 'the first byte of a block is not cleared by the emission: bits are emitted from the wrong end', '(byte)(acc >> (accBits - 8))', '(byte)(acc >> (accBits - 7))'),
    ],
    'MappedTextReader.cs': [
        m('R01', 'the key holds one symbol less (63 / b)', '_keySymbols = 64 / bits;', '_keySymbols = 63 / bits;'),
        m('R02', 'the key shift counts the unused bits the wrong way', '_keyShift = 64 - (_keySymbols * bits);', '_keyShift = 64 + (_keySymbols * bits);'),
        m('R03', 'the vector step covers 128 bits', '_vectorSymbols = 256 / bits;', '_vectorSymbols = 128 / bits;'),
        m('R04', 'the vector flag is not guarded by the hardware', 'set => _useVectors = value && Avx2.IsSupported;', 'set => _useVectors = value;'),
        m('R05', 'a negative position is not recognised by the window', '(ulong)i < (ulong)_limit ? WindowAt(i)', 'i < _limit ? WindowAt(i)'),
        m('R06', 'position n is read through the window', '(ulong)i < (ulong)_limit ? WindowAt(i)', '(ulong)i <= (ulong)_limit ? WindowAt(i)'),
        m('R07', 'the key drops one bit too many', 'Window(i) >> _keyShift;', 'Window(i) >> (_keyShift + 1);'),
        m('R08', 'the symbol is taken from the wrong bits', '(uint)(Window(i) >> (64 - _bits))', '(uint)(Window(i) >> (63 - _bits))'),
        m('R09', 'ReadSymbols does not check the disposal', '        ObjectDisposedException.ThrowIf(_disposed, this);\n        if (start < 0', '        if (start < 0'),
        m('R10', 'ReadSymbols refuses start = n', 'start > _n ||', 'start >= _n ||'),
        m('R11', 'ReadSymbols refuses a range that ends at n', 'destination.Length > _n - start', 'destination.Length >= _n - start'),
        m('R12', 'the byte fast path of ReadSymbols is taken for 7 bits', 'if (_bits == 8)', 'if (_bits == 7)'),
        m('R13', 'the byte fast path of ReadSymbols is never taken', 'if (_bits == 8)', 'if (_bits == 99)'),
        m('R14', 'ReadSymbols: shift off by one', 'int shift = 16 - (int)(bit & 7) - bits;', 'int shift = 15 - (int)(bit & 7) - bits;'),
        m('R15', 'ReadSymbols: no mask', '>> shift) & mask)', '>> shift))'),
        m('R16', 'ReadSymbols: the bit position advances by one per symbol', 'for (int k = 0; k < destination.Length; k++, bit += bits)', 'for (int k = 0; k < destination.Length; k++, bit += 1)'),
        m('R17', 'Lce: an empty room is not an early return', '        if (room <= 0)\n        {', '        if (room < 0)\n        {'),
        m('R18', 'Lce: the limit is not clamped to the room', '        if (max > room) max = room;\n', ''),
        m('R19', 'Lce: the room is measured from the smaller position', 'long room = _limit - Math.Max(i, j);', 'long room = _limit - Math.Min(i, j);'),
        m('R20', 'Lce: the vector path starts one step late', 'if (_useVectors && max >= _vectorSymbols)', 'if (_useVectors && max > _vectorSymbols)'),
        m('R21', 'Lce: the vector loop leaves the last full step to the scalar path', 'while (k + vs <= max)', 'while (k + vs < max)'),
        m('R22', 'Lce: the vector steps are not counted', '                steps++;\n                Vector256<ulong> a', '                Vector256<ulong> a'),
        m('R23', 'Lce: the lane of the first difference is found with a 4-bit unit', 'TrailingZeroCount(~mask) >> 3', 'TrailingZeroCount(~mask) >> 2'),
        m('R24', 'Lce: the matched symbols of a lane are counted from its end', '((lane * 64) + BitOperations.LeadingZeroCount(x)) / b', '((lane * 64) - BitOperations.LeadingZeroCount(x)) / b'),
        m('R25', 'Lce: a difference in the partial symbol at the end of the window ends the common prefix', 'if (matched < vs) return', 'if (matched <= vs) return'),
        m('R26', 'Lce: the vector path continues one symbol short', 'k += vs; // the difference lies', 'k += vs - 1; // the difference lies'),
        m('R27', 'Lce: scalar: a difference in the padding bit ends the common prefix', 'if (matched < ks) return', 'if (matched <= ks) return'),
        m('R28', 'Lce: scalar steps are not counted', '            steps++;\n            ulong x = WindowAt(i + k)', '            ulong x = WindowAt(i + k)'),
        m('R29', 'Lce: scalar advances one symbol short', '            k += ks;\n        }\n        return max;', '            k += ks - 1;\n        }\n        return max;'),
        m('R30', 'Lce: the scalar result is not clamped to max', 'return Math.Min(k + matched, max);\n            }\n            k += ks;', 'return k + matched;\n            }\n            k += ks;'),
        m('R31', 'Lce: the vector result is not clamped to max', 'if (matched < vs) return Math.Min(k + matched, max);', 'if (matched < vs) return k + matched;'),
        m('R32', 'Lce: a negative i is not refused', 'if ((i | j | max) < 0) ThrowNegative(i, j, max);', 'if ((j | max) < 0) ThrowNegative(i, j, max);'),
        m('R33', 'Lce: a negative max is not refused', 'if ((i | j | max) < 0) ThrowNegative(i, j, max);', 'if ((i | j) < 0) ThrowNegative(i, j, max);'),
        m('R71', 'Lce: a read after disposal is not refused', '            ObjectDisposedException.ThrowIf(_disposed, this);\n            return 0;', '            return 0;'),
        m('R72', 'Lce: a negative j is not refused', 'if ((i | j | max) < 0) ThrowNegative(i, j, max);', 'if ((i | max) < 0) ThrowNegative(i, j, max);'),
        m('R34', 'CompareSuffixes: positions beyond n are not clamped', '        i = Math.Min(i, _n);\n        j = Math.Min(j, _n);\n', ''),
        m('R35', 'CompareSuffixes: the same position is not special', '        if (i == j)\n        {', '        if (i == j && false)\n        {'),
        m('R36', 'CompareSuffixes: the skip is not clamped', '        skip = Math.Min(skip, maxLen);\n', ''),
        m('R37', 'CompareSuffixes: the shorter suffix is the larger', 'return i > j ? -1 : 1; // one suffix ended', 'return i > j ? 1 : -1; // one suffix ended'),
        m('R38', 'CompareSuffixes: the symbol order is reversed', 'return Symbol(i + m) < Symbol(j + m) ? -1 : 1;', 'return Symbol(i + m) > Symbol(j + m) ? -1 : 1;'),
        m('R39', 'CompareSuffixes: the skipped symbols are not counted', 'long m = skip + Lce(', 'long m = Lce('),
        m('R40', 'Verify: the length is not compared', '        if (info.N != _n) throw', '        if (info.N < 0) throw'),
        m('R41', 'Verify: the width is not compared', '        if (info.Bits != _bits) throw', '        if (info.Bits < 0) throw'),
        m('R42', 'Verify: an alphabet that is not as long as the counts is accepted', '        if (info.Counts.Length != sigma) throw', '        if (info.Counts.Length < 0) throw'),
        m('R43', 'Verify: sigma = 2^b is refused', 'if (sigma > (1 << _bits)) throw', 'if (sigma >= (1 << _bits)) throw'),
        m('R44', 'Verify: a symbol equal to sigma is accepted', 'if (s >= sigma) throw', 'if (s > sigma) throw'),
        m('R45', 'Verify: the symbols are not counted', '                seen[s]++;\n', ''),
        m('R46', 'Verify: the symbols are not mapped back for the hash', '                block[k] = info.Letters[s];\n', ''),
        m('R47', 'Verify: the counts are not compared', 'if (seen[s] != info.Counts[s]) throw', 'if (seen[s] < 0) throw'),
        m('R48', 'Verify: the hash is not compared', 'if (hash is not null && !string.Equals(', 'if (hash is null && !string.Equals('),
        m('R49', 'Verify: the hash comparison is case sensitive', 'StringComparison.OrdinalIgnoreCase', 'StringComparison.Ordinal'),
        m('R50', 'Verify: the blocks are read from the start every time', 'ReadSymbols(at, block.AsSpan(0, take));', 'ReadSymbols(0, block.AsSpan(0, take));'),
        m('R51', 'Dispose: the limit is not cleared (a read after disposal touches the unmapped view)', '        _limit = 0;\n        if (_view is not null)', '        if (_view is not null)'),
        m('R52', 'OpenCore: a file of exactly the length is refused', 'if (offset > fileLength - total) throw', 'if (offset >= fileLength - total) throw'),
        m('R53', 'OpenCore: a longer file is accepted as a whole file', 'if (sectionLength is null && fileLength != total) throw', 'if (sectionLength is null && fileLength < total) throw'),
        m('R54', 'OpenCore: a section length that differs is accepted', 'if (sectionLength is long declared && declared != total) throw', 'if (sectionLength is long declared && declared < 0) throw'),
        m('R55', 'OpenCore: a negative offset is accepted', 'if (offset < 0) throw', 'if (offset < -1) throw'),
        m('R56', 'OpenCore: the view starts at 0 instead of the offset', 'map.CreateViewAccessor(offset, total,', 'map.CreateViewAccessor(0, total,'),
        m('R57', 'CheckLayout: a short region is accepted', '        if (length < total) throw', '        if (length < total - 1) throw'),
        m('R58', 'CheckLayout: a long region is accepted', '        if (length > total) throw', '        if (length > total + 1) throw'),
        m('R59', 'CheckLayout: the lowest unused bit is not checked', '((1 << unused) - 1)) != 0)', '((1 << unused) - 2)) != 0)'),
        m('R60', 'CheckLayout: the unused bits are not checked', 'if (unused > 0 && (data[dataBytes - 1]', 'if (unused > 99 && (data[dataBytes - 1]'),
        m('R61', 'CheckLayout: the padding is not checked', 'if (bad >= 0) throw', 'if (bad >= 1) throw'),
        m('R62', 'CheckShape: the lower width is not checked', 'bits < PackedTextFormat.MinBits ||', 'bits < 0 ||'),
        m('R63', 'CheckShape: the upper width is not checked', '|| bits > PackedTextFormat.MaxBits)', '|| bits > 32)'),
        m('R64', 'CheckShape: a negative length is accepted', 'if (n < 0 || n >', 'if (n < -1 || n >'),
        m('R65', 'CheckShape: the length is not bounded', '|| n > long.MaxValue / PackedTextFormat.MaxBits - PackedTextFormat.PadBytes)', '|| n > long.MaxValue)'),
        m('R66', 'the window carries the wrong bits from the ninth byte', '(ulong)q[8] >> (8 - s)', '(ulong)q[8] >> (7 - s)'),
        m('R67', 'the window is not byte-reversed', 'BinaryPrimitives.ReverseEndianness(Unsafe.ReadUnaligned<ulong>(q))', 'Unsafe.ReadUnaligned<ulong>(q)'),
        m('R68', 'the window drops the shift when the bit offset is zero', 'return s == 0 ? w : (w << s)', 'return s == 1 ? w : (w << s)'),
        m('R69', 'the outside window does not refuse a negative position', '        ArgumentOutOfRangeException.ThrowIfNegative(i);\n        return 0;', '        return 0;'),
        m('R70', 'the outside window does not check the disposal', '        ObjectDisposedException.ThrowIf(_disposed, this);\n        ArgumentOutOfRangeException.ThrowIfNegative(i);\n        return 0;', '        ArgumentOutOfRangeException.ThrowIfNegative(i);\n        return 0;'),
    ],
}

FILES['MappedOutputFile.cs'] = [
    m('O01', 'Auto means sparse', 'internal const ReservationMode AutoMode = ReservationMode.Preallocate;', 'internal const ReservationMode AutoMode = ReservationMode.Sparse;'),
    m('O02', 'Auto is not resolved', 'options.Reservation == ReservationMode.Auto ? AutoMode : options.Reservation', 'options.Reservation'),
    m('O03', 'a negative length is not refused up front', '        ArgumentOutOfRangeException.ThrowIfNegative(length);\n        options ??=', '        options ??='),
    m('O04', 'a null path is not refused up front', '        ArgumentNullException.ThrowIfNull(path);\n        ArgumentOutOfRangeException.ThrowIfNegative(length);\n        options ??=', '        ArgumentOutOfRangeException.ThrowIfNegative(length);\n        options ??='),
    m('O05', 'an empty file is mapped', 'length > 0 ? MemoryMappedFile.CreateFromFile', 'length >= 0 ? MemoryMappedFile.CreateFromFile'),
    m('O06', 'a failed creation leaves the file', '        catch\n        {\n            DeleteQuietly(path);\n            throw;', '        catch\n        {\n            throw;'),
    m('O07', 'CreateWriter: a range that ends at the end of the file is refused', 'if (start > Length - length) throw', 'if (start >= Length - length) throw'),
    m('O08', 'CreateWriter: a negative start is not refused', '        ArgumentOutOfRangeException.ThrowIfNegative(start);\n', ''),
    m('O09', 'CreateWriter: a negative length is not refused', '        ArgumentOutOfRangeException.ThrowIfNegative(length);\n        ArgumentOutOfRangeException.ThrowIfLessThan(windowBytes, 1);', '        ArgumentOutOfRangeException.ThrowIfLessThan(windowBytes, 1);'),
    m('O10', 'CreateWriter: a window of 0 bytes is accepted', '        ArgumentOutOfRangeException.ThrowIfLessThan(windowBytes, 1);\n        if (start >', '        if (start >'),
    m('O11', 'CreateWriter: the writer is not counted', '        Interlocked.Increment(ref _writers);\n        return new MappedWindowWriter(this, start', '        return new MappedWindowWriter(this, start'),
    m('O12', 'CreateWriter does not check the disposal', '        ObjectDisposedException.ThrowIf(_disposed, this);\n        ArgumentOutOfRangeException.ThrowIfNegative(start);', '        ArgumentOutOfRangeException.ThrowIfNegative(start);'),
    m('O13', 'Complete tolerates one open writer', 'if (open != 0) throw', 'if (open > 1) throw'),
    m('O14', 'Complete does not check the disposal', '        ObjectDisposedException.ThrowIf(_disposed, this);\n        int open', '        int open'),
    m('O15', 'Complete: fsync is always skipped', '        if (!fsync) return;\n', '        return;\n'),
    m('O16', 'Complete: fsync is always done', '        if (!fsync) return;\n', ''),
    m('O17', 'Dispose is not idempotent-guarded', '        if (_disposed) return;\n        _disposed = true;\n        _map?.Dispose();', '        _disposed = true;\n        _map?.Dispose();'),
    m('O18', 'Dispose does not mark the file disposed', '        _disposed = true;\n        _map?.Dispose();', '        _map?.Dispose();'),
    m('O19', 'Dispose does not release the mapping', '        _disposed = true;\n        _map?.Dispose();', '        _disposed = true;'),
    m('O20', 'OpenView does not check the disposal', '            ObjectDisposedException.ThrowIf(_disposed, this);\n            return _map!', '            return _map!'),
    m('O21', 'OpenView opens a read-only view', 'CreateViewAccessor(fileOffset, size, MemoryMappedFileAccess.Write)', 'CreateViewAccessor(fileOffset, size, MemoryMappedFileAccess.Read)'),
    m('O22', 'WriterClosed does not decrement', 'internal void WriterClosed() => Interlocked.Decrement(ref _writers);', 'internal void WriterClosed() { }'),
    m('O23', 'Sparse: the length is not set', '            sparse.SetLength(length);\n', ''),
    m('O24', 'Sparse: the file is not truncated first', 'using var sparse = new FileStream(path, FileMode.Create,', 'using var sparse = new FileStream(path, FileMode.OpenOrCreate,'),
    m('O25', 'Preallocate: the real source is used when it cannot preallocate', 'if (source.CanPreallocate)\n        {\n            using FileStream reserved', 'if (!source.CanPreallocate)\n        {\n            using FileStream reserved'),
    m('O26', 'Preallocate: the length is not set after the reservation', '            reserved.SetLength(length);\n', ''),
    m('O27', 'Fallback: an unknown free space refuses the file', 'if (free >= 0 && free < length) throw', 'if (free < length) throw'),
    m('O28', 'Fallback: free space equal to the length is refused', 'if (free >= 0 && free < length) throw', 'if (free >= 0 && free <= length) throw'),
    m('O29', 'Fallback: the free space is never compared', 'if (free >= 0 && free < length) throw', 'if (free >= 0 && free < 0) throw'),
    m('O30', 'Fallback: the outcome labels are swapped', 'return free >= 0 ? ReservationOutcome.FreeSpaceChecked : ReservationOutcome.Unchecked;', 'return free >= 0 ? ReservationOutcome.Unchecked : ReservationOutcome.FreeSpaceChecked;'),
    m('O31', 'Fallback: the free space is asked for the file, not its directory', 'long free = source.AvailableFreeBytes(directory);', 'long free = source.AvailableFreeBytes(path);'),
    m('O32', 'Fallback: the file is not truncated first', 'using var checkedFile = new FileStream(path, FileMode.Create,', 'using var checkedFile = new FileStream(path, FileMode.OpenOrCreate,'),
    m('O33', 'Fallback: the length is not set', '        checkedFile.SetLength(length);\n', ''),
    m('O34', 'the system source preallocates everywhere', 'public bool CanPreallocate => OperatingSystem.IsLinux();', 'public bool CanPreallocate => true;'),
    m('O35', 'the system source never preallocates', 'public bool CanPreallocate => OperatingSystem.IsLinux();', 'public bool CanPreallocate => false;'),
    m('O36', 'the system source reserves nothing', 'PreallocationSize = length }', 'PreallocationSize = 0 }'),
    m('O37', 'the system source does not truncate', 'Mode = FileMode.Create, Access = FileAccess.ReadWrite, Share = FileShare.ReadWrite, PreallocationSize', 'Mode = FileMode.OpenOrCreate, Access = FileAccess.ReadWrite, Share = FileShare.ReadWrite, PreallocationSize'),
    m('O38', 'the free space is wrong by a factor', 'return new DriveInfo(directory).AvailableFreeSpace;', 'return new DriveInfo(directory).AvailableFreeSpace / 2;'),
    m('O39', 'an unknown free space is 0', '            return -1;\n        }\n    }\n}', '            return 0;\n        }\n    }\n}'),
    m('O40', 'DeleteQuietly deletes nothing', '            File.Delete(path);\n        }\n        catch (Exception e) when (e is IOException or UnauthorizedAccessException)\n        {\n            // the path is not ours', '            _ = path;\n        }\n        catch (Exception e) when (e is IOException or UnauthorizedAccessException)\n        {\n            // the path is not ours'),
    m('O41', 'the writer counter starts at one', 'private int _writers;', 'private int _writers = 1;'),
]

FILES['MappedWindowWriter.cs'] = [
    m('W01', 'a range writer maps its window at the range-relative offset', '_view = _file.OpenView(_start + offset, size);', '_view = _file.OpenView(offset, size);'),
    m('W02', 'Flush does nothing', 'internal void Flush() => CloseWindow();', 'internal void Flush() { }'),
    m('W03', 'Dispose is not idempotent-guarded', '        if (_closed) return;\n        _closed = true;', '        _closed = true;'),
    m('W04', 'Dispose does not close the window', '        _closed = true;\n        CloseWindow();', '        _closed = true;'),
    m('W05', 'Dispose of a range writer disposes the shared file', '        if (_ownsFile) _file.Dispose();\n        else _file.WriterClosed();', '        _file.Dispose();'),
    m('W06', 'Dispose of a range writer does not report the closing', '        if (_ownsFile) _file.Dispose();\n        else _file.WriterClosed();', '        if (_ownsFile) _file.Dispose();'),
    m('W07', 'the owning writer does not dispose its file', '        if (_ownsFile) _file.Dispose();\n        else _file.WriterClosed();', '        if (!_ownsFile) _file.WriterClosed();'),
    m('W08', 'FileOffset is 0', 'internal long FileOffset => _start;', 'internal long FileOffset => 0;'),
    m('W09', 'a write that ends at the end of the range is refused', 'if (offset < 0 || offset > Length - data.Length) throw', 'if (offset < 0 || offset >= Length - data.Length) throw'),
    m('W10', 'the owning writer creates its file with a reservation', 'new MappedFileOptions { Reservation = ReservationMode.Sparse }', 'new MappedFileOptions { Reservation = ReservationMode.Auto }'),
    m('W11', 'the owning writer does not count itself', '        _file.WriterOpened();\n', ''),
    m('W12', 'the owning writer forgets that it owns the file', '        _ownsFile = true;\n', '        _ownsFile = false;\n'),
    m('W13', 'the window of a range writer is as large as the rest of the file', 'long size = Math.Min(_window, Length - offset);', 'long size = Math.Min(_window, _file.Length - offset);'),
    m('W14', 'Path is the empty string', 'internal string Path => _file.Path;', 'internal string Path => string.Empty;'),
]

FILES['PackedTextWriter.cs'] = [
    m('P01', 'a block that is not a multiple of eight is accepted', 'o.BlockSymbols < BlockMultiple || o.BlockSymbols % BlockMultiple != 0', 'o.BlockSymbols < BlockMultiple'),
    m('P02', 'a block of 8 is refused', 'o.BlockSymbols < BlockMultiple ||', 'o.BlockSymbols <= BlockMultiple ||'),
    m('P03', 'the lower MinBits is not checked', '        ArgumentOutOfRangeException.ThrowIfLessThan(o.MinBits, PackedTextFormat.MinBits);\n', ''),
    m('P04', 'the upper MinBits is not checked', '        ArgumentOutOfRangeException.ThrowIfGreaterThan(o.MinBits, PackedTextFormat.MaxBits);\n', ''),
    m('P05', 'the window is not checked', '        ArgumentOutOfRangeException.ThrowIfLessThan(o.WindowBytes, 1);\n', ''),
    m('P06', 'Scan: the limit is exclusive', 'if (n > options.MaxSymbols) throw', 'if (n >= options.MaxSymbols) throw'),
    m('P07', 'Scan: the limit is not applied', 'if (n > options.MaxSymbols) throw', 'if (n > long.MaxValue) throw'),
    m('P08', 'Scan: the symbols are not counted', '                foreach (byte b in block.AsSpan(0, take)) counts[b]++;\n', ''),
    m('P09', 'Scan: the hash is not fed', '                hash.AppendData(block.AsSpan(0, take));\n', ''),
    m('P10', 'Scan: the hash is not reset before reading (it is read once)', 'sha = Convert.ToHexStringLower(hash.GetHashAndReset());', 'sha = Convert.ToHexString(hash.GetHashAndReset());'),
    m('P11', 'Scan: the letters are the ones that occur more than once', 'if (counts[v] != 0) letters.Add((byte)v);', 'if (counts[v] > 1) letters.Add((byte)v);'),
    m('P12', 'Scan: the width is the larger of nothing', 'int bits = Math.Max(options.MinBits, PackedTextFormat.BitsFor(letters.Count));', 'int bits = Math.Min(options.MinBits, PackedTextFormat.BitsFor(letters.Count));'),
    m('P13', 'Scan: the width ignores the minimum', 'int bits = Math.Max(options.MinBits, PackedTextFormat.BitsFor(letters.Count));', 'int bits = PackedTextFormat.BitsFor(letters.Count);'),
    m('P14', 'Scan: the counters are not given to the input', 'using (var input = new CountedFile(inputPath, FileMode.Open, FileAccess.Read, options.Counters))\n        using (var hash', 'using (var input = new CountedFile(inputPath, FileMode.Open, FileAccess.Read))\n        using (var hash'),
    m('P15', 'Scan: the blocks are read from the start', 'input.ReadExactAt(block.AsSpan(0, take), at);\n                foreach', 'input.ReadExactAt(block.AsSpan(0, take), 0);\n                foreach'),
    m('P16', 'Pack: a target with more room is accepted', 'if (room != info.TotalBytes) throw', 'if (room < info.TotalBytes) throw'),
    m('P17', 'Pack: a target with less room is accepted', 'if (room != info.TotalBytes) throw', 'if (room > info.TotalBytes) throw'),
    m('P18', 'Pack: a changed length is not noticed', 'if (n != info.N) throw', 'if (n < 0) throw'),
    m('P19', 'Pack: the alphabet map is not built', 'for (int s = 0; s < info.Letters.Length; s++) map[info.Letters[s]] = s;', 'for (int s = 0; s < info.Letters.Length; s++) map[info.Letters[s]] = 0;'),
    m('P20', 'Pack: a byte that is not in the alphabet is not noticed', 'if (symbol < 0) throw', 'if (symbol < -1) throw'),
    m('P21', 'Pack: the symbols are not counted', '                    seen[symbol]++;\n', ''),
    m('P22', 'Pack: the counts are not compared', 'if (seen[s] != info.Counts[s]) throw', 'if (seen[s] < 0) throw'),
    m('P23', 'Pack: a block is appended one byte short', 'target.Append(packed.AsSpan(0, bytes));', 'target.Append(packed.AsSpan(0, bytes - 1));'),
    m('P24', 'Pack: the padding is not written', '        target.Append(new byte[PackedTextFormat.PadBytes]);\n', ''),
    m('P25', 'Pack: the padding is one byte short', 'target.Append(new byte[PackedTextFormat.PadBytes]);', 'target.Append(new byte[PackedTextFormat.PadBytes - 1]);'),
    m('P26', 'Pack: the data length is not checked', 'if (target.Position - start != info.DataBytes) throw', 'if (target.Position - start < 0) throw'),
    m('P27', 'Pack: the packed buffer holds one byte less than a block needs', 'var packed = new byte[(options.BlockSymbols * (long)bits / 8) + 1];', 'var packed = new byte[(options.BlockSymbols * (long)bits / 8)];'),
    m('P28', 'Pack: the counters are not given to the input', 'using (var input = new CountedFile(inputPath, FileMode.Open, FileAccess.Read, options.Counters))\n        {\n            for', 'using (var input = new CountedFile(inputPath, FileMode.Open, FileAccess.Read))\n        {\n            for'),
    m('P29', 'WriteFile: the directory is not checked first', 'if (directory is not null && !Directory.Exists(directory)) throw', 'if (directory is not null && Directory.Exists(directory) && false) throw'),
    m('P30', 'WriteFile: the hook between the passes is not called', '        options.BetweenPasses?.Invoke();\n', ''),
    m('P31', 'WriteFile: the file is written in place', 'string temporary = packedPath + ".tmp";', 'string temporary = packedPath;'),
    m('P32', 'WriteFile: the window option is not used', 'file.CreateWriter(options.WindowBytes, options.Counters)', 'file.CreateWriter(MappedWindowWriter.DefaultWindowBytes, options.Counters)'),
    m('P33', 'WriteFile: the counters are not given to the writer', 'file.CreateWriter(options.WindowBytes, options.Counters)', 'file.CreateWriter(options.WindowBytes)'),
    m('P34', 'WriteFile: the file options are not used', 'MappedOutputFile.Create(temporary, info.TotalBytes, options.File)', 'MappedOutputFile.Create(temporary, info.TotalBytes)'),
    m('P35', 'WriteFile: fsync is always skipped', 'file.Complete(options.Fsync);', 'file.Complete(false);'),
    m('P36', 'WriteFile: an existing target is not replaced', 'File.Move(temporary, packedPath, overwrite: true);', 'File.Move(temporary, packedPath, overwrite: false);'),
    m('P37', 'WriteFile: a failure leaves the temporary file', '            file?.Dispose();\n            DeleteQuietly(temporary);\n            throw;', '            file?.Dispose();\n            throw;'),
    m('P38', 'WriteFile: the file is not disposed before the rename', '            file.Dispose();\n            file = null;\n', '            file = null;\n'),
    m('P39', 'WriteFile: the file is disposed again by the handler (null guard lost)', '            file.Dispose();\n            file = null;\n', '            file.Dispose();\n'),
    m('P40', 'InputLength: the length of the link, not of the text', '        using var input = new CountedFile(inputPath, FileMode.Open, FileAccess.Read);\n        return input.Length;', '        return new FileInfo(inputPath).Length;'),
    m('P41', 'WriteFile: the info is of an empty text', '        return info;\n    }\n\n    private static void CheckOptions', '        return info with { Sha256 = string.Empty };\n    }\n\n    private static void CheckOptions'),
]

# Mutants that were written and dropped before the first run to bound the time of the audit (economy cap of the brief): near duplicates of a kept mutant, or equivalent by inspection of the code that the
# kept ones already stress (the 5-hour cap leaves about 1.5 hours of machine time for the audit). They were never run and are not counted.
DROPPED = {'F03', 'F13', 'F18', 'R06', 'R13', 'R17', 'R22', 'R28', 'R55', 'R60', 'R63', 'R64', 'R65', 'R67', 'R68', 'O03', 'O04', 'O12', 'O14', 'O17', 'O18', 'O19', 'O20', 'O24', 'O32', 'O37', 'O38', 'O39',
           'O40', 'O41', 'W11', 'W12', 'W13', 'W14', 'P02', 'P05', 'P07', 'P10', 'P13', 'P15', 'P27', 'P29', 'P34', 'P38', 'P39', 'P41'}

# Fresh mutants of round 3 (written after the review of the round 1 survivors; the vector loads, the verification loop, the layout check, the scan).
FRESH = {
    'MappedTextReader.cs': [
        m('N01', 'Lce: the vector step is the key length', 'int vs = _vectorSymbols;', 'int vs = _keySymbols;'),
        m('N02', 'Lce: the lane shuffle swaps two bytes', '7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8,\n        7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8);', '7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8,\n        7, 6, 5, 4, 3, 2, 0, 1, 15, 14, 13, 12, 11, 10, 9, 8);'),
        m('N03', 'LoadShifted: the carry shift is off by one', 'Vector128.CreateScalar(64 - s)', 'Vector128.CreateScalar(63 - s)'),
        m('N04', 'LoadShifted: the carry is loaded from one byte earlier', 'Vector256.Load(q + 8)', 'Vector256.Load(q + 7)'),
        m('N05', 'Lce: the first difference is always taken from lane 0', 'ulong x = a.GetElement(lane) ^ c.GetElement(lane);', 'ulong x = a.GetElement(0) ^ c.GetElement(0);'),
        m('N06', 'Lce: equal vectors advance one symbol short', '                    k += vs;\n                    continue;', '                    k += vs - 1;\n                    continue;'),
        m('N08', 'Verify: counts every symbol twice', '                seen[s]++;\n', '                seen[s] += 2;\n'),
        m('N09', 'Verify: one letter too many is accepted', 'if (sigma > (1 << _bits)) throw', 'if (sigma > (1 << _bits) + 1) throw'),
        m('N10', 'ReadSymbols: the mask is one bit short', 'int mask = (1 << bits) - 1;', 'int mask = (1 << bits) - 2;'),
        m('N11', 'CheckLayout: the number of unused bits is one too many', 'int unused = (int)((dataBytes * 8) - (n * bits));', 'int unused = (int)((dataBytes * 8) - (n * bits)) + 1;'),
        m('N12', 'Verify: the blocks overlap by one symbol', 'for (long at = 0; at < _n; at += block.Length)', 'for (long at = 0; at < _n; at += block.Length - 1)'),
    ],
    'PackedTextWriter.cs': [
        m('N13', 'Scan: the blocks overlap by one symbol', 'for (long at = 0; at < n; at += block.Length)\n            {\n                int take = (int)Math.Min(block.Length, n - at);\n                input.ReadExactAt(block.AsSpan(0, take), at);\n                foreach', 'for (long at = 0; at < n; at += block.Length - 1)\n            {\n                int take = (int)Math.Min(block.Length, n - at);\n                input.ReadExactAt(block.AsSpan(0, take), at);\n                foreach'),
        m('N14', 'WriteFile: a directory that exists is refused', 'if (directory is not null && !Directory.Exists(directory)) throw', 'if (directory is not null && Directory.Exists(directory)) throw'),
    ],
    'MappedWindowWriter.cs': [
        m('N17', 'the window is as large as the larger of the window and the rest', 'long size = Math.Min(_window, Length - offset);', 'long size = Math.Max(_window, Length - offset);'),
    ],
    'PackedTextFormat.cs': [
        m('N19', 'PackInto starts writing at byte 1', '        int o = 0;\n        foreach (byte s in symbols)', '        int o = 1;\n        foreach (byte s in symbols)'),
    ],
}

HERE = os.path.dirname(os.path.abspath(__file__))
for name in sorted(set(FILES) | set(FRESH)):
    text = open(os.path.join(ROOT, SRC + name)).read()
    all_ = FILES.get(name, [])
    kept = [x for x in all_ if x['id'] not in DROPPED]
    r3 = [x for x in all_ if x['id'] in DROPPED] + FRESH.get(name, [])
    for tag, lst in (('', kept), ('r3_', r3)):
        for x in lst:
            if text.count(x['old']) != 1:
                raise SystemExit('%s %s: pattern occurs %d times' % (name, x['id'], text.count(x['old'])))
            x['file'] = SRC + name
        if not lst:
            continue
        out = os.path.join(HERE, '1b7b_%s%s_listed.json' % (tag, name[:-3].lower()))
        if tag == '' and os.path.exists(out) and '--force' not in sys.argv:
            continue  # round 1 was run from the files that exist; they are not rewritten
        json.dump(lst, open(out, 'w'), indent=1)
        print(tag or 'r1_', name, len(lst))
