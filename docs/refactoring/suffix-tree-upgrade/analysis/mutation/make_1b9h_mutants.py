#!/usr/bin/env python3
"""Mutants of the audit of work package 1b.9h (prefix 1b9h_): the code that the package adds or changes in src/ and in the test support. Four files, all single-site, every `old` occurs exactly once at the head:
   RankTable.cs (R): the bounds checks of the checked build (run with CHECKED=1, run_1b9h_mutants.py), MappedTextReader.cs (T): the same for the mapped text (CHECKED=1), DifferenceCover.cs (C): the table of the default period
   and its check when it is used (the production build), HangLimitAttribute.cs (H): the bookkeeping of the time limit (the production build; test support, so the mutant is a change of a test file and its tests are HangBreakerTests).
   usage: make_1b9h_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/External/'
REG = '(FullyQualifiedName~MappedFileRegistryTests|FullyQualifiedName~RankTableBounds1b9hTests|FullyQualifiedName~MappedTextBounds1b9hTests)'
FILES = {
    'ranktable': (SRC + 'Dcs/RankTable.cs', '(FullyQualifiedName~RankTableBounds1b9hTests|FullyQualifiedName~RankTableTests|FullyQualifiedName~SuffixOrderTests)'),
    'mappedtextreader': (SRC + 'MappedTextReader.cs', '(FullyQualifiedName~MappedTextBounds1b9hTests|FullyQualifiedName~MappedTextReaderTests)'),
    'differencecover': (SRC + 'Dcs/DifferenceCover.cs', '(FullyQualifiedName~DifferenceCoverTable1b9hTests|FullyQualifiedName~DifferenceCoverTests|FullyQualifiedName~SaCoverGap1b9cTests)'),
    'hanglimit': ('tests/SuffixTree/SuffixTree.Esa.Tests/TestSupport/HangLimitAttribute.cs', '(FullyQualifiedName~HangBreakerTests)'),
    'mappedfileregistry': (SRC + 'MappedFileRegistry.cs', REG),
    'countedfile': (SRC + 'CountedFile.cs', REG),
    'mappedoutputfile': (SRC + 'MappedOutputFile.cs', REG),
    'esaanalysis': ('src/SuffixTree/Algorithms/SuffixTree.Esa/EsaIndex.Analysis.cs', '(FullyQualifiedName~EsaDistinctSubstringTests)'),
    'discoverylimit': ('tests/SuffixTree/SuffixTree.Esa.Tests/TestSupport/DiscoveryLimit.cs', '(FullyQualifiedName~DiscoveryLimitTests)'),
    'allocationmeter': ('tests/SuffixTree/SuffixTree.Esa.Tests/TestSupport/AllocationMeter.cs', '(FullyQualifiedName~AllocationMeterTests|FullyQualifiedName~EsaDistinctSubstringTests)'),
}
G = {k: [] for k in FILES}


def m(key, i, desc, old, new):
    G[key].append({'id': '1b9h_' + i, 'desc': desc, 'old': old, 'new': new})


R = lambda *a: m('ranktable', *a)  # noqa: E731
R('R01', 'ReadChecked: the last position of the text is refused', '(ulong)p > (ulong)N ||', '(ulong)p >= (ulong)N ||')
R('R02', 'ReadChecked: one position beyond the text is accepted', '(ulong)p > (ulong)N ||', '(ulong)p > (ulong)(N + 1) ||')
R('R03', 'ReadChecked: no check of the position (a negative or a far position is folded into the table)', '(ulong)p > (ulong)N || index < 0 ||', 'index < 0 ||')
R('R04', 'ReadChecked: a residue outside the sample reads an arbitrary slot', ' || index < 0 || (ulong)slot >= (ulong)Slots', ' || (ulong)slot >= (ulong)Slots')
R('R05', 'ReadChecked: the residue check refuses only below minus one', 'index < 0 || (ulong)slot', 'index < -1 || (ulong)slot')
R('R06', 'ReadChecked: the last slot of the table is refused', '(ulong)slot >= (ulong)Slots || _view is null) ThrowOutOfBounds(p, r, slot);', '(ulong)slot >= (ulong)Slots - 1 || _view is null) ThrowOutOfBounds(p, r, slot);')
R('R07', 'ReadChecked: no check of the slot', ' || (ulong)slot >= (ulong)Slots || _view is null) ThrowOutOfBounds(p, r, slot);', ' || _view is null) ThrowOutOfBounds(p, r, slot);')
R('R08', 'ReadChecked: a read after Dispose is accepted', '(ulong)slot >= (ulong)Slots || _view is null) ThrowOutOfBounds(p, r, slot);', '(ulong)slot >= (ulong)Slots) ThrowOutOfBounds(p, r, slot);')
R('R09', 'ReadChecked: the slot ignores the class of the residue', 'long slot = ((long)block * _classes) + index;\n        if ((ulong)p', 'long slot = ((long)block * _classes);\n        if ((ulong)p')
R('R10', 'ReadChecked: the slot is of the next block', 'long slot = ((long)block * _classes) + index;\n        if ((ulong)p', 'long slot = (((long)block + 1) * _classes) + index;\n        if ((ulong)p')
R('R11', 'Rank: the checked build reads through the raw pointer', '        return ReadChecked(p, block, r);\n#else', '        return _ranks[((long)block * _classes) + _indexOf[r]];\n#else')
R('R12', 'ThrowOutOfBounds: the message does not say that the table is disposed', '{(_view is null ? "; the table is disposed" : string.Empty)}', '{(_view is null ? string.Empty : string.Empty)}')
R('R13', 'ThrowOutOfBounds: the violation is of the other kind', 'MemoryViolationKind.OutOfBounds,\n            $"out-of-bounds read of the rank table', 'MemoryViolationKind.CanaryOverwritten,\n            $"out-of-bounds read of the rank table')
R('R14', 'ReadSlot: the last slot is refused', 'if ((ulong)slot >= (ulong)Slots || _view is null) ThrowOutOfBounds(-1, 0, slot);', 'if ((ulong)slot >= (ulong)Slots - 1 || _view is null) ThrowOutOfBounds(-1, 0, slot);')
R('R15', 'ReadSlot: no check of the slot', 'if ((ulong)slot >= (ulong)Slots || _view is null) ThrowOutOfBounds(-1, 0, slot);', 'if (_view is null) ThrowOutOfBounds(-1, 0, slot);')
R('R16', 'ReadSlot: a read after Dispose is accepted', 'if ((ulong)slot >= (ulong)Slots || _view is null) ThrowOutOfBounds(-1, 0, slot);', 'if ((ulong)slot >= (ulong)Slots) ThrowOutOfBounds(-1, 0, slot);')
R('R17', 'Verify: the checked build reads the slots through the raw pointer', '            uint v = ReadSlot(slot);\n#else', '            uint v = _ranks[slot];\n#else')
R('R18', 'Verify after Dispose does not throw the disposed exception', 'ObjectDisposedException.ThrowIf(_view is null, this);\n        long samples', 'long samples')

T = lambda *a: m('mappedtextreader', *a)  # noqa: E731
T('T01', 'the checked region is the data without the padding', '        _bytes = PackedTextFormat.TotalBytes(n, bits);', '        _bytes = PackedTextFormat.DataBytes(n, bits);')
T('T02', 'the checked region is one byte short', '        _bytes = PackedTextFormat.TotalBytes(n, bits);', '        _bytes = PackedTextFormat.TotalBytes(n, bits) - 1;')
T('T03', 'the checked region is one byte long', '        _bytes = PackedTextFormat.TotalBytes(n, bits);', '        _bytes = PackedTextFormat.TotalBytes(n, bits) + 1;')
T('T04', 'CheckRead: a read that ends exactly at the end of the region is refused', 'if (offset < 0 || count > _bytes - offset) ThrowOutOfBounds(offset, count);', 'if (offset < 0 || count >= _bytes - offset) ThrowOutOfBounds(offset, count);')
T('T05', 'CheckRead: a read of one byte beyond the region is accepted', 'if (offset < 0 || count > _bytes - offset) ThrowOutOfBounds(offset, count);', 'if (offset < 0 || count > _bytes - offset + 1) ThrowOutOfBounds(offset, count);')
T('T06', 'CheckRead: a negative offset is accepted', 'if (offset < 0 || count > _bytes - offset) ThrowOutOfBounds(offset, count);', 'if (count > _bytes - offset) ThrowOutOfBounds(offset, count);')
T('T07', 'CheckRead: the count is not looked at', 'if (offset < 0 || count > _bytes - offset) ThrowOutOfBounds(offset, count);', 'if (offset < 0 || 0 > _bytes - offset) ThrowOutOfBounds(offset, count);')
T('T08', 'WindowAt: the window is checked as eight bytes (the ninth carries the shift)', '        CheckRead(bit >> 3, 9);', '        CheckRead(bit >> 3, 8);')
T('T09', 'WindowAt: the window is checked as ten bytes', '        CheckRead(bit >> 3, 9);', '        CheckRead(bit >> 3, 10);')
T('T10', 'WindowAt: the check is of the symbol and not of the byte', '        CheckRead(bit >> 3, 9);', '        CheckRead(i, 9);')
T('T11', 'WindowAt: no check', '        CheckRead(bit >> 3, 9);\n', '')
T('T12', 'LoadShifted: the load of 32 + 8 bytes is checked as 32', '        CheckRead(bit >> 3, 40);', '        CheckRead(bit >> 3, 32);')
T('T13', 'LoadShifted: the load is checked as 41 bytes', '        CheckRead(bit >> 3, 40);', '        CheckRead(bit >> 3, 41);')
T('T14', 'LoadShifted: no check', '        CheckRead(bit >> 3, 40);\n', '')
T('T15', 'ReadSymbols: the two bytes of the window are not counted', '(((long)destination.Length * _bits) + 7 >> 3) + 2', '(((long)destination.Length * _bits) + 7 >> 3)')
T('T16', 'ReadSymbols: the bytes are rounded down', '(((long)destination.Length * _bits) + 7 >> 3) + 2', '(((long)destination.Length * _bits) >> 3) + 2')
T('T17', 'ReadSymbols: the first byte is taken from the symbol number', 'CheckRead(start * _bits >> 3,', 'CheckRead(start,')
T('T18', 'ReadSymbols: no check', '        CheckRead(start * _bits >> 3, (((long)destination.Length * _bits) + 7 >> 3) + 2);\n', '')
T('T19', 'ThrowOutOfBounds: the violation is of the other kind', 'MemoryViolationKind.OutOfBounds,\n            $"out-of-bounds read of the mapped text', 'MemoryViolationKind.CanaryOverwritten,\n            $"out-of-bounds read of the mapped text')

C = lambda *a: m('differencecover', *a)  # noqa: E731
C('C01', 'Create: the table is used for every period (and checked for none)', 'if (period == DefaultSingerPeriod && FromTable(period, tableOfTheDefaultPeriod) is { } known) return known;', 'if (FromTable(period, tableOfTheDefaultPeriod) is { } known) return known;')
C('C02', 'Create: the table is never used', 'if (period == DefaultSingerPeriod && FromTable(period, tableOfTheDefaultPeriod) is { } known) return known;', 'if (period != period && FromTable(period, tableOfTheDefaultPeriod) is { } known) return known;')
C('C03', 'Create: the table is used without the check', 'if (period == DefaultSingerPeriod && FromTable(period, tableOfTheDefaultPeriod) is { } known) return known;', 'if (period == DefaultSingerPeriod) return new DifferenceCover(period, tableOfTheDefaultPeriod.ToArray(), CoverKind.Singer);')
C('C04', 'Create: the default period is a different one', 'internal const int DefaultSingerPeriod = 993;', 'internal const int DefaultSingerPeriod = 553;')
C('C05', 'Create: the cover of the table is of the kind Pruned', 'var cover = new DifferenceCover(period, residues.ToArray(), CoverKind.Singer);', 'var cover = new DifferenceCover(period, residues.ToArray(), CoverKind.Pruned);')
C('C07', 'FromTable: the table of q residues is accepted', 'residues.Length != q + 1) return null;', 'residues.Length != q) return null;')
C('C08', 'FromTable: a table of any length is accepted', 'if (SingerPrime(period) is not int q || residues.Length != q + 1) return null;', 'if (SingerPrime(period) is not int q) return null;')
C('C09', 'FromTable: a period that is no Singer period is accepted (the length is not looked at)', 'if (SingerPrime(period) is not int q || residues.Length != q + 1) return null;', 'int q = 31; if (residues.Length != q + 1) return null;')
C('C10', 'FromTable: a residue of zero is refused', 'if (residues[i] < 0 || residues[i] >= period ||', 'if (residues[i] < 1 || residues[i] >= period ||')
C('C11', 'FromTable: a residue below zero is accepted', 'if (residues[i] < 0 || residues[i] >= period ||', 'if (residues[i] < -1 || residues[i] >= period ||')
C('C12', 'FromTable: a residue equal to the period is accepted', 'if (residues[i] < 0 || residues[i] >= period ||', 'if (residues[i] < 0 || residues[i] > period ||')
C('C13', 'FromTable: a repeated residue is accepted (sorted, not strictly)', '(i > 0 && residues[i] <= residues[i - 1])', '(i > 0 && residues[i] < residues[i - 1])')
C('C14', 'FromTable: the order is not looked at', '|| (i > 0 && residues[i] <= residues[i - 1])) return null;', ') return null;')
C('C15', 'FromTable: the second residue is not compared with the first', '(i > 0 && residues[i] <= residues[i - 1])', '(i > 1 && residues[i] <= residues[i - 1])')
C('C16', 'FromTable: either check is enough', 'return cover.Verify() && cover.IsPerfect() ? cover : null;', 'return cover.Verify() || cover.IsPerfect() ? cover : null;')
C('C17', 'FromTable: the perfect-set check is dropped', 'return cover.Verify() && cover.IsPerfect() ? cover : null;', 'return cover.Verify() ? cover : null;')
C('C18', 'FromTable: the cover check is dropped', 'return cover.Verify() && cover.IsPerfect() ? cover : null;', 'return cover.IsPerfect() ? cover : null;')
C('C19', 'FromTable: the table is not checked at all', 'return cover.Verify() && cover.IsPerfect() ? cover : null;', 'return cover;')
C('C20', 'SingerSet: the search is not counted', '        Interlocked.Increment(ref s_searches);\n', '')
C('C21', 'SingerSet: the search is counted twice', '        Interlocked.Increment(ref s_searches);\n', '        Interlocked.Add(ref s_searches, 2);\n')
C('C22', 'CreateBySearch: the table is used', 'if (SingerPrime(period) is int q) return new DifferenceCover(period, SingerSet(q), CoverKind.Singer);\n        return new DifferenceCover(period, Pruned(period), CoverKind.Pruned);\n    }\n\n    /// <summary>\n    /// The Singer cover', 'if (SingerPrime(period) is int q) return FromTable(period, DefaultSingerResidues) ?? new DifferenceCover(period, SingerSet(q), CoverKind.Singer);\n        return new DifferenceCover(period, Pruned(period), CoverKind.Pruned);\n    }\n\n    /// <summary>\n    /// The Singer cover')
C('C23', 'the table: one residue is another perfect set (shifted by one)', '[1, 2, 21, 28, 67, 108, 209, 222, 233, 291, 300, 306, 314, 350, 458, 468, 486, 543, 588, 631, 686, 734, 805, 810, 845, 857, 861, 882, 899, 931, 961, 992]', '[0, 2, 3, 22, 29, 68, 109, 210, 223, 234, 292, 301, 307, 315, 351, 459, 469, 487, 544, 589, 632, 687, 735, 806, 811, 846, 858, 862, 883, 900, 932, 962]')
C('C24', 'the table: one residue is changed', '[1, 2, 21, 28, 67, 108,', '[1, 3, 21, 28, 67, 108,')
C('C25', 'the table: the last residue is dropped', '931, 961, 992];', '931, 961];')

H = lambda *a: m('hanglimit', *a)  # noqa: E731
H('H01', 'Begin: a test begun again under its own name is taken for a runaway', 'if (before is not null && before != name)', 'if (before is not null)')
H('H02', 'Begin: a runaway is never recorded', 'if (before is not null && before != name) Interlocked.CompareExchange(ref _runaway, before, null);', 'if (before is not null && before != name) { }')
H('H03', 'Begin: the first runaway is replaced by a later one', 'Interlocked.CompareExchange(ref _runaway, before, null);', 'Interlocked.Exchange(ref _runaway, before);')
H('H04', 'Begin: the test that begins is marked as inside also when it is refused', 'if (runaway is null) Volatile.Write(ref _inside, name);', 'Volatile.Write(ref _inside, name);')
H('H05', 'Begin: the test that begins is never marked as inside', 'if (runaway is null) Volatile.Write(ref _inside, name);', 'if (runaway is null) { }')
H('H06', 'Begin: no test is refused', 'string? runaway = Volatile.Read(ref _runaway);\n        if (runaway is null)', 'string? runaway = null;\n        if (runaway is null)')
H('H07', 'End: the end of any test clears the mark', 'Interlocked.CompareExchange(ref _inside, null, name);', 'Interlocked.Exchange(ref _inside, null);')
H('H08', 'End: the end of a test does not clear the mark', 'internal void End(string name) => Interlocked.CompareExchange(ref _inside, null, name);', 'internal void End(string name) { }')
H('H09', 'End: the end of the runaway opens the breaker', 'internal void End(string name) => Interlocked.CompareExchange(ref _inside, null, name);', 'internal void End(string name) { Interlocked.CompareExchange(ref _inside, null, name); if (_runaway == name) _runaway = null; }')
H('H10', 'ApplyToContext: the limit is in tenths of seconds', 'checked(Seconds * 1000);', 'checked(Seconds * 100);')
H('H11', 'ApplyToContext: the product is not checked', 'checked(Seconds * 1000);', 'unchecked(Seconds * 1000);')
H('H13', 'the attribute: the limit of a method is in tenths of seconds', 'PropertyAttribute(PropertyNames.Timeout, checked(seconds * 1000))', 'PropertyAttribute(PropertyNames.Timeout, checked(seconds * 100))')
H('H14', 'the attribute: the limit of a method is not the limit of the property Timeout', 'PropertyAttribute(PropertyNames.Timeout, checked(seconds * 1000))', 'PropertyAttribute(PropertyNames.Category, checked(seconds * 1000))')
H('H12', 'BeforeTest: a refused test is not failed', '            Assert.Fail($"not started:', '            Assert.Pass($"not started:')

# the registry of the mapped files (the checked build; run with CHECKED=1) and where it is called from
Gm = lambda *a: m('mappedfileregistry', *a)  # noqa: E731
Gm('G01', 'Map: the path is not made absolute', '        string full = Path.GetFullPath(path);\n        lock (s_mapped) s_mapped[full]', '        string full = path;\n        lock (s_mapped) s_mapped[full]')
Gm('G02', 'Map: a second mapping of the file does not count', 's_mapped[full] = s_mapped.GetValueOrDefault(full) + 1;', 's_mapped[full] = 1;')
Gm('G03', 'Map: a mapping counts twice', 's_mapped[full] = s_mapped.GetValueOrDefault(full) + 1;', 's_mapped[full] = s_mapped.GetValueOrDefault(full) + 2;')
Gm('G04', 'Unmap: the file is removed while a mapping is left', 'if (count <= 1) s_mapped.Remove(full);', 'if (count <= 2) s_mapped.Remove(full);')
Gm('G05', 'Unmap: the last mapping is not removed', 'if (count <= 1) s_mapped.Remove(full);', 'if (count < 1) s_mapped.Remove(full);')
Gm('G06', 'Unmap: the count is not reduced', 'else s_mapped[full] = count - 1;', 'else s_mapped[full] = count;')
Gm('G07', 'Unmap: an unknown file is an error', 'if (!s_mapped.TryGetValue(full, out int count)) return;', 'int count = s_mapped[full];')
Gm('G08', 'RequireNotMapped: the path is not made absolute', '        string full = Path.GetFullPath(path);\n        bool mapped;', '        string full = path;\n        bool mapped;')
Gm('G09', 'RequireNotMapped: nothing is refused', '        if (mapped)\n        {\n            throw', '        if (false && mapped)\n        {\n            throw')
Gm('G10', 'RequireNotMapped: everything is refused', '        if (mapped)\n        {\n            throw', '        if (true)\n        {\n            throw')
Gm('G11', 'RequireNotMapped: the violation is of the other kind', 'throw new EsaMemoryViolationException(MemoryViolationKind.OutOfBounds, $"the file', 'throw new EsaMemoryViolationException(MemoryViolationKind.CanaryOverwritten, $"the file')
Gm('G12', 'Count: always zero', '            lock (s_mapped) return s_mapped.Count;', '            lock (s_mapped) return 0;')
Gm('G13', 'Count: the mappings, not the files', '            lock (s_mapped) return s_mapped.Count;', '            lock (s_mapped) return s_mapped.Values.Sum();')
Gm('G14', 'the names are compared without regard to case', 'new(StringComparer.Ordinal);', 'new(StringComparer.OrdinalIgnoreCase);')

Km = lambda *a: m('countedfile', *a)  # noqa: E731
Km('K01', 'CountedFile: a file that is mapped may be created', 'if (mode is FileMode.Create or FileMode.Truncate) MappedFileRegistry.RequireNotMapped(path);', 'if (mode is FileMode.Create or FileMode.Truncate) { }')
Km('K02', 'CountedFile: only Create is looked at (not Truncate)', 'if (mode is FileMode.Create or FileMode.Truncate) MappedFileRegistry.RequireNotMapped(path);', 'if (mode is FileMode.Create) MappedFileRegistry.RequireNotMapped(path);')
Km('K03', 'CountedFile: only Truncate is looked at (not Create)', 'if (mode is FileMode.Create or FileMode.Truncate) MappedFileRegistry.RequireNotMapped(path);', 'if (mode is FileMode.Truncate) MappedFileRegistry.RequireNotMapped(path);')
Km('K04', 'CountedFile: a file that is mapped may not be opened for reading either', 'if (mode is FileMode.Create or FileMode.Truncate) MappedFileRegistry.RequireNotMapped(path);', 'MappedFileRegistry.RequireNotMapped(path);')
Km('K05', 'CountedFile: OpenOrCreate is looked at too', 'if (mode is FileMode.Create or FileMode.Truncate) MappedFileRegistry.RequireNotMapped(path);', 'if (mode is FileMode.Create or FileMode.Truncate or FileMode.OpenOrCreate) MappedFileRegistry.RequireNotMapped(path);')

Om = lambda *a: m('mappedoutputfile', *a)  # noqa: E731
Om('K06', 'MappedOutputFile.Create: a file that is mapped may be created', '        MappedFileRegistry.RequireNotMapped(path);\n', '')

Rm = lambda *a: m('ranktable', *a)  # noqa: E731
Rm('R19', 'Open: the table is not put in the list of the mapped files', 'table._registered = MappedFileRegistry.Map(path);', 'table._registered = null;')
Rm('R20', 'Dispose: the table is not taken out of the list', 'if (view is not null && _registered is not null) MappedFileRegistry.Unmap(_registered);', 'if (view is not null && _registered is not null) { }')
Rm('R21', 'Dispose: a second Dispose takes it out again', 'if (view is not null && _registered is not null) MappedFileRegistry.Unmap(_registered);', 'if (_registered is not null) MappedFileRegistry.Unmap(_registered);')

Tm = lambda *a: m('mappedtextreader', *a)  # noqa: E731
Tm('T20', 'OpenCore: the reader is not put in the list of the mapped files', 'reader._registered = MappedFileRegistry.Map(path);', 'reader._registered = null;')
Tm('T21', 'Dispose: the reader is not taken out of the list', '        if (_registered is not null) MappedFileRegistry.Unmap(_registered);\n', '')
Tm('T22', 'Dispose: a second Dispose takes it out again (the field is not cleared first)', '        if (_registered is not null) MappedFileRegistry.Unmap(_registered);\n', '        if (_registered is not null) { MappedFileRegistry.Unmap(_registered); MappedFileRegistry.Unmap(_registered); }\n')

# the discovery limit (test support)
Dm = lambda *a: m('discoverylimit', *a)  # noqa: E731
Dm('D01', 'the exception of the computation is lost', 'error = ExceptionDispatchInfo.Capture(e);', 'error = null;')
Dm('D02', 'the exception of the computation is not thrown again', '        error?.Throw();\n', '')
Dm('D03', 'the computation runs on a foreground thread', 'IsBackground = true, Name = "discovery-limit"', 'IsBackground = false, Name = "discovery-limit"')
Dm('D04', 'the limit is not used', 'if (!done.Wait(limit))', 'if (!done.Wait(Timeout.InfiniteTimeSpan))')
Dm('D05', 'a computation that does not return gives true', '            value = default;\n            return false;', '            value = default;\n            return true;')
Dm('D06', 'a computation that returns gives false', '        value = result;\n        return true;', '        value = result;\n        return false;')
Dm('D07', 'a computation that does not return leaves its (partial) value', '            value = default;\n            return false;', '            value = result;\n            return false;')
Dm('D08', 'the done signal is not set when the computation fails', '            finally\n            {\n                done.Set();\n            }', '            finally\n            {\n            }\n            if (error is null) done.Set();')
Dm('D09', 'the computation is run on the calling thread', '        thread.Start();\n        if (!done.Wait(limit))', '        thread.Start(); thread.Join();\n        if (!done.Wait(limit))')

# the allocation meter (test support)
Am = lambda *a: m('allocationmeter', *a)  # noqa: E731
Am('A01', 'AssertLessThan: the limit is inclusive', 'if (Bytes < limit) return;', 'if (Bytes <= limit) return;')
Am('A02', 'AssertLessThan: the limit is one byte lower', 'if (Bytes < limit) return;', 'if (Bytes < limit - 1) return;')
Am('A03', 'AssertLessThan: a reading in a window without collections is inconclusive', 'if (GcFree) Assert.Fail(', 'if (!GcFree) Assert.Fail(')
Am('A04', 'AssertLessThan: a reading in a window with a collection is a pass', '        Assert.Inconclusive(', '        Assert.Pass(')
Am('A05', 'AssertLessThan: a reading in a window without collections is a pass', 'if (GcFree) Assert.Fail(', 'if (GcFree) Assert.Pass(')
Am('A06', 'Measure: no region is opened', 'bool opened = !alreadyOpen && TryStartRegion();', 'bool opened = false;')
Am('A07', 'Measure: the region is never left', '            if (opened) EndRegion();', '            if (opened) { }')
Am('A08', 'Measure: the reading is the counter after the code', 'return new AllocationReading(after - before, gcFree);', 'return new AllocationReading(after, gcFree);')
Am('A09', 'Measure: the counter is read before the region is opened', '        bool alreadyOpen = GCSettings.LatencyMode == GCLatencyMode.NoGCRegion;\n        bool opened = !alreadyOpen && TryStartRegion();\n        try\n        {\n            int collections = GC.CollectionCount(0);                // every collection counts in generation 0, whatever its generation\n            long before = GC.GetAllocatedBytesForCurrentThread();\n            code();', '        long before = GC.GetAllocatedBytesForCurrentThread();\n        bool alreadyOpen = GCSettings.LatencyMode == GCLatencyMode.NoGCRegion;\n        bool opened = !alreadyOpen && TryStartRegion();\n        try\n        {\n            int collections = GC.CollectionCount(0);                // every collection counts in generation 0, whatever its generation\n            code();')
Am('A10', 'Measure: the window is free of collections whatever happened', 'bool gcFree = (alreadyOpen || opened) && GCSettings.LatencyMode == GCLatencyMode.NoGCRegion && GC.CollectionCount(0) == collections;', 'bool gcFree = true;')
Am('A11', 'Measure: the window is free of collections when a region was opened', 'bool gcFree = (alreadyOpen || opened) && GCSettings.LatencyMode == GCLatencyMode.NoGCRegion && GC.CollectionCount(0) == collections;', 'bool gcFree = alreadyOpen || opened;')
Am('A12', 'Measure: the region need not have been opened', 'bool gcFree = (alreadyOpen || opened) && GCSettings.LatencyMode == GCLatencyMode.NoGCRegion && GC.CollectionCount(0) == collections;', 'bool gcFree = GCSettings.LatencyMode == GCLatencyMode.NoGCRegion && GC.CollectionCount(0) == collections;')
Am('A13', 'Measure: the latency mode is not looked at', 'bool gcFree = (alreadyOpen || opened) && GCSettings.LatencyMode == GCLatencyMode.NoGCRegion && GC.CollectionCount(0) == collections;', 'bool gcFree = (alreadyOpen || opened) && GC.CollectionCount(0) == collections;')
Am('A14', 'Measure: the number of collections is not looked at', 'bool gcFree = (alreadyOpen || opened) && GCSettings.LatencyMode == GCLatencyMode.NoGCRegion && GC.CollectionCount(0) == collections;', 'bool gcFree = (alreadyOpen || opened) && GCSettings.LatencyMode == GCLatencyMode.NoGCRegion;')
Am('A15', 'Measure: the code is not checked for null', '        ArgumentNullException.ThrowIfNull(code);\n        bool alreadyOpen', '        bool alreadyOpen')
Am('A16', 'the budget of the region is too small for the large array', 'RegionBytes = 64L << 20', 'RegionBytes = 1L << 20')
Am('A17', 'the region has no room for large objects', 'RegionLargeObjectBytes = 32L << 20', 'RegionLargeObjectBytes = 1L << 10')
Am('A18', 'TryStartRegion: a region that another thread opened is an error', '        catch (InvalidOperationException)           // another thread opened a region between the check and the start: the window is measured without the promise (and that region is gone: a failed start ends it)\n        {\n            return false;\n        }\n', '')
Am('A21', 'Measure: a region that is already open is not looked for (a start inside it ends it)', 'bool opened = !alreadyOpen && TryStartRegion();', 'bool opened = TryStartRegion();')
Am('A22', 'Measure: a region that was open is ended by the meter', '            if (opened) EndRegion();', '            if (opened || alreadyOpen) EndRegion();')
Am('A19', 'EndRegion: the region is ended whatever the mode', '            if (GCSettings.LatencyMode == GCLatencyMode.NoGCRegion) GC.EndNoGCRegion();', '            GC.EndNoGCRegion();')
Am('A20', 'EndRegion: the exception of a region that was left is not caught', '        catch (InvalidOperationException)           // the region was left by a collection (the budget was exceeded or a collection was induced) between the check and the call: nothing is left to end\n        {\n        }\n', '        finally { }\n')

# the code that the 1b.3 test measures: mutants that allocate (the test must stay red against each)
Zm = lambda *a: m('esaanalysis', *a)  # noqa: E731
Zm('Z01', 'CountDistinctSubstrings allocates a long per symbol', '        long n = Length;\n        long total = n * (n + 1) / 2;\n        CompactLcp.Cursor lcp', '        long n = Length;\n        GC.KeepAlive(new long[Length]);\n        long total = n * (n + 1) / 2;\n        CompactLcp.Cursor lcp')
Zm('Z02', 'CountDistinctSubstrings allocates a byte per symbol', '        long n = Length;\n        long total = n * (n + 1) / 2;\n        CompactLcp.Cursor lcp', '        long n = Length;\n        GC.KeepAlive(new byte[Length]);\n        long total = n * (n + 1) / 2;\n        CompactLcp.Cursor lcp')
Zm('Z03', 'CountDistinctSubstrings allocates 2 KB', '        long n = Length;\n        long total = n * (n + 1) / 2;\n        CompactLcp.Cursor lcp', '        long n = Length;\n        GC.KeepAlive(new byte[2000]);\n        long total = n * (n + 1) / 2;\n        CompactLcp.Cursor lcp')
Zm('Z04', 'ByLength allocates a second array of the size of the result', '        var counts = new long[m + 1];\n', '        var counts = new long[m + 1];\n        GC.KeepAlive(new long[m + 1]);\n')
Zm('Z05', 'ByLength(10) allocates a byte per symbol', '        var counts = new long[m + 1];\n', '        var counts = new long[m + 1];\n        GC.KeepAlive(new byte[n]);\n')
Zm('Z06', 'ByLength allocates a copy of the result on return', '        counts[0] = 1;\n        return counts;', '        counts[0] = 1;\n        return (long[])counts.Clone();')
Zm('Z07', 'CountDistinctSubstrings allocates a boxed value per rank', '        for (int rank = 1; rank <= Length; rank++) total -= lcp.Get(rank);', '        for (int rank = 1; rank <= Length; rank++) { object box = lcp.Get(rank); total -= (int)box; }')

out = []
for key, (rel, filt) in FILES.items():
    src = open(os.path.join(ROOT, rel), encoding='utf-8').read()
    lst = []
    for x in G[key]:
        n = src.count(x['old'])
        if n != 1:
            raise SystemExit('%s: pattern occurs %d times: %r' % (x['id'], n, x['old'][:90]))
        if x['old'] == x['new']:
            raise SystemExit('%s: not a change' % x['id'])
        d = dict(x)
        d['file'] = rel
        d['filter'] = filt
        lst.append(d)
    json.dump(lst, open(os.path.join(HERE, '1b9h_%s_listed.json' % key), 'w'), indent=1)
    out.append((key, len(lst)))
print('listed', out, sum(n for _, n in out), 'mutants')
