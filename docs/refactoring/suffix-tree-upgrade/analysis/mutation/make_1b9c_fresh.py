#!/usr/bin/env python3
"""Work package 1b.9c: fresh single-site mutants, written after the new tests (none of them is in the lists of 1b.7d or 1b.7e or of earlier packages), over the code that the new tests of this package touch: the test seams that
were added, and the checks, the pass-throughs and the cleanups that the new tests pin. They test the new tests themselves. Each is built and run by run_1b9c_mutants.py with stage 1 (the listed filter = the new test classes) and stage 2
(the whole default set) when it survives. Writes 1b9c_f_<file>_listed.json. A pattern that does not match exactly once is reported."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '../../../../..'))
E = 'src/SuffixTree/Algorithms/SuffixTree.Esa/External/'
NEW = 'FullyQualifiedName~Gap1b9cTests'
M = [
    # file, id, description, old, new
    ('Ph6/PsiStream.cs', 'F01', 'the live-block counter is not raised at an allocation', '            Interlocked.Increment(ref s_liveBlocks);\n', ''),
    ('Ph6/PsiStream.cs', 'F02', 'the live-block counter is raised at a release', 'Interlocked.Decrement(ref s_liveBlocks);', 'Interlocked.Increment(ref s_liveBlocks);'),
    ('Ph6/PsiStream.cs', 'F03', 'the live-block counter is not lowered at a release', '            Interlocked.Decrement(ref s_liveBlocks);\n', ''),
    ('IoCounters.cs', 'F04', 'a flush is counted twice', 'internal void AddFlush() => Interlocked.Increment(ref _flushes);', 'internal void AddFlush() => Interlocked.Add(ref _flushes, 2);'),
    ('IoCounters.cs', 'F05', 'Reset keeps the flushes', '        Interlocked.Exchange(ref _flushes, 0);\n', ''),
    ('IoCounters.cs', 'F06', 'Flushes reads the counter of the write calls', 'internal long Flushes => Interlocked.Read(ref _flushes);', 'internal long Flushes => Interlocked.Read(ref _writeCalls);'),
    ('CountedFile.cs', 'F07', 'a flush is not counted', '        _counters?.AddFlush();\n', ''),
    ('Sort/BucketSorter.cs', 'F08', 'no scan point in the partition', '                    _o.Trace?.Invoke("scan"); // a point at which a test fails the partition after some of the writers hold files\n', ''),
    ('Sort/BucketSorter.cs', 'F09', 'no hand-over point in the partition', '                        _o.Trace?.Invoke("handover"); // a point at which a test fails the hand-over after some children exist\n', ''),
    ('Sort/BucketSorter.cs', 'F10', 'the children are disposed after a hand-over that succeeded', 'if (!handedOver) foreach (BuildTask child in children) child.Bucket.Dispose();', 'if (handedOver) foreach (BuildTask child in children) child.Bucket.Dispose();'),
    ('Sort/BucketSorter.cs', 'F11', 'the entry limit of the sort is at the limit', 'if (text.Length > PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{text.Length} symbols: positions are 32-bit, the limit is {PackedTextFormat.MaxSymbols}");\n        string dir', 'if (text.Length >= PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{text.Length} symbols: positions are 32-bit, the limit is {PackedTextFormat.MaxSymbols}");\n        string dir'),
    ('Dcs/ExternalDcs.cs', 'F12', 'no rank-written point after the merge of unique names', '            ls.MergeSeconds = Seconds(t0);\n            o.Trace?.Invoke("rank-written");\n            return;', '            ls.MergeSeconds = Seconds(t0);\n            return;'),
    ('Dcs/ExternalDcs.cs', 'F13', 'no rank-written point after the base case', '            ls.BelowSeconds = Seconds(t0);\n            o.Trace?.Invoke("rank-written");\n            return;', '            ls.BelowSeconds = Seconds(t0);\n            return;'),
    ('Dcs/ExternalDcs.cs', 'F14', 'no rank-written point after the inversion', '        ls.BelowSeconds = Seconds(t0);\n        o.Trace?.Invoke("rank-written");\n    }', '        ls.BelowSeconds = Seconds(t0);\n    }'),
    ('Dcs/ExternalDcs.cs', 'F15', 'the entry limit of PH1 is at the limit', 'if (text.Length > PackedTextFormat.MaxSymbols) throw new NotSupportedException', 'if (text.Length >= PackedTextFormat.MaxSymbols) throw new NotSupportedException'),
    ('Sa/SaBuilder.cs', 'F16', 'the rank file is checked when it is not asked for', '            if (options.VerifyRanks)', '            if (!options.VerifyRanks)'),
    ('Sa/SaBuilder.cs', 'F17', 'the length of the suffix array file is checked with a smaller-than', 'if (got != 4 * (n + 1)) throw', 'if (got < 4 * (n + 1)) throw'),
    ('Sa/SaBuilder.cs', 'F18', 'PH1 gets no trace', '                Files = options.Files,\n                Trace = options.Trace,\n            });\n            long rankBytes', '                Files = options.Files,\n            });\n            long rankBytes'),
    ('Sesa/SesaBuilder.cs', 'F19', 'the first sync only when it is not asked for', '                if (options.Fsync) file.Flush();\n                Hook(SesaHooks.HeaderSynced);', '                if (!options.Fsync) file.Flush();\n                Hook(SesaHooks.HeaderSynced);'),
    ('Sesa/SesaBuilder.cs', 'F20', 'the second sync only when it is not asked for', '                if (options.Fsync) file.Flush();\n                Hook(SesaHooks.MagicWritten);', '                if (!options.Fsync) file.Flush();\n                Hook(SesaHooks.MagicWritten);'),
    ('Sesa/SesaBuilder.cs', 'F21', 'the counts of the BWT are compared with themselves', 'if (!bwt.Counts.AsSpan().SequenceEqual(info.Counts)) throw', 'if (!bwt.Counts.AsSpan().SequenceEqual(bwt.Counts)) throw'),
    ('Sesa/SesaBuilder.cs', 'F22', 'a container that is shorter than planned is accepted', 'if (actual != planned) throw new InvalidOperationException($"the container has', 'if (actual > planned) throw new InvalidOperationException($"the container has'),
    ('Sesa/SesaBuilder.cs', 'F23', 'the message of the length check names the planned length twice', 'the container has {actual} bytes, expected {planned}', 'the container has {planned} bytes, expected {planned}'),
    ('Sesa/SesaBuilder.cs', 'F24', 'BUILD_INFO is refused only above twice its room', 'if (length > capacity) throw new InvalidOperationException($"BUILD_INFO has', 'if (length > 2 * capacity) throw new InvalidOperationException($"BUILD_INFO has'),
    ('Sesa/SesaChecksums.cs', 'F25', 'the seconds of a section are negative', '(Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency, Crc32C.HardwareAccelerated', '(t0 - Stopwatch.GetTimestamp()) / (double)Stopwatch.Frequency, Crc32C.HardwareAccelerated'),
    ('Sesa/SesaChecksums.cs', 'F26', 'the seconds of a section are in ticks', '(Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency, Crc32C.HardwareAccelerated', '(Stopwatch.GetTimestamp() - t0), Crc32C.HardwareAccelerated'),
    ('Sesa/SesaChecksums.cs', 'F27', 'the row says the hardware path was not used', 'Crc32C.HardwareAccelerated));', '!Crc32C.HardwareAccelerated));'),
    ('Sesa/SesaChecksums.cs', 'F28', 'the bytes of a row are one short', 'new SesaCrcRow(kind, length,', 'new SesaCrcRow(kind, length - 1,'),
    ('Lcp/LcpBuild.cs', 'F29', 'the overflow writer is not disposed after a failed output writer', '                baseWriter.Dispose();\n                overflowWriter.Dispose();\n                throw;', '                baseWriter.Dispose();\n                throw;'),
    ('Lcp/LcpBuild.cs', 'F30', 'the base writer is not disposed after a failed output writer', '                baseWriter.Dispose();\n                overflowWriter.Dispose();\n                throw;', '                overflowWriter.Dispose();\n                throw;'),
    ('MappedTextReader.cs', 'F31', 'a length of -1 is a length of a text', 'if (n < 0 || n > long.MaxValue / PackedTextFormat.MaxBits - PackedTextFormat.PadBytes) throw new InvalidDataException($"{n} symbols: not a length of a text");', 'if (n < -1 || n > long.MaxValue / PackedTextFormat.MaxBits - PackedTextFormat.PadBytes) throw new InvalidDataException($"{n} symbols: not a length of a text");'),
    ('Dcs/SampleNamer.cs', 'F32', 'a sample one record above the capacity is a leaf', '_g.Samples > _o.BucketCapacity, 1)', '_g.Samples > _o.BucketCapacity + 1, 1)'),
    ('Dcs/NameSinks.cs', 'F33', 'the file of the rank sink is not disposed', '        _writer.Dispose();\n        _file.Dispose();', '        _writer.Dispose();'),
    ('Sort/BucketSorter.cs', 'F34', 'the bound of the leaf capacity is half of it', 'ArgumentOutOfRangeException.ThrowIfGreaterThan(options.LeafCapacity, int.MaxValue / 2);', 'ArgumentOutOfRangeException.ThrowIfGreaterThan(options.LeafCapacity, int.MaxValue / 4);'),
    ('Sesa/SesaBuilder.cs', 'F35', 'PH6 is run for a build with threads above 1000 only', '            if (options.IncludePsi)\n            {\n                using var ph6', '            if (options.IncludePsi && options.Threads <= 1000)\n            {\n                using var ph6'),
    ('Sesa/SesaBuilder.cs', 'F36', 'the BWT scratch file is not given the file options of the build', 'Counters = counters,\n            Files = options.Files,\n        });\n        double secBwt', 'Counters = counters,\n        });\n        double secBwt'),
    ('Sesa/SectionLcpChildBuilder.cs', 'F37', 'PH5 flushes only the last window', 'FlushWindows = true,', 'FlushWindows = options.Threads > 100,'),
    ('Sesa/SectionLcpChildBuilder.cs', 'F38', 'PH4 gets the guards of PH5 only', 'ArenaGuards = options.ArenaGuards,\n                Counters = counters,\n                Files = options.Files,\n                Trace', 'ArenaGuards = options.ArenaGuards is null ? null : null,\n                Counters = counters,\n                Files = options.Files,\n                Trace'),
    # the second batch: the code of the tests written after the first batch (the position store, the report of the suffix array build, the fsync and the growth of a file, the pruned covers, the BWT pass)
    ('Sort/PositionStore.cs', 'F39', 'the floor of the capacity of a bucket buffer is 17 positions', '_capacity = Math.Max(spillBytes / 4, 16);', '_capacity = Math.Max(spillBytes / 4, 17);'),
    ('Sort/PositionStore.cs', 'F40', 'a bucket in a file is handed over with its file open', '        Flush();\n        _file.Dispose();\n        _file = null;', '        Flush();\n        _file = null;'),
    ('Sort/PositionStore.cs', 'F41', 'a position file gives one entry less back to the ledger', '        files.Ledger.Remove(count * 4);', '        files.Ledger.Remove((count - 1) * 4);'),
    ('Sort/PositionStore.cs', 'F42', 'a writer that is disposed uncompleted keeps its buffer', '        if (_buffer.Length > 0) ArrayPool<uint>.Shared.Return(_buffer);\n', ''),
    ('Sort/PositionStore.cs', 'F43', 'a memory source can give its array back twice', 'uint[]? d = Interlocked.Exchange(ref _data, null);', 'uint[]? d = _data;'),
    ('Sort/PositionStore.cs', 'F44', 'the second half of a flushed buffer is written at the wrong offset', '(_total - _filled) * 4);', '(_total - _filled) * 2);'),
    ('Sa/SaBuilder.cs', 'F45', 'PH1 time is the time of PH3', 'Ph1Seconds = ph1Seconds,', 'Ph1Seconds = ph23Seconds,'),
    ('Sa/SaBuilder.cs', 'F46', 'the entry limit of the build is at the limit', 'if (n > PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{n} symbols: positions are 32-bit', 'if (n >= PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{n} symbols: positions are 32-bit'),
    ('Sa/SaBuilder.cs', 'F47', 'the whole-run RssAnon figure is in 1000 byte units', 'double? anonMb = all.Available ? all.MaxRssAnonKb * 1024.0 / 1e6 : null;', 'double? anonMb = all.Available ? all.MaxRssAnonKb * 1000.0 / 1e6 : null;'),
    ('Sa/SaBuilder.cs', 'F48', 'the fsync of the output is not counted', 'saFile.Complete(options.SyncOutput, counters);', 'saFile.Complete(options.SyncOutput, null);'),
    ('Sa/SaBuilder.cs', 'F49', 'the rank file is not booked in the ledger', '            ledger.Add(rankBytes); // the rank file is scratch of the pipeline from here on\n', ''),
    ('MappedOutputFile.cs', 'F50', 'the fsync of a file is not counted', '        counters?.AddFlush();\n', ''),
    ('MappedOutputFile.cs', 'F51', 'a failed mapping sets the length again to the new length', 'undo.SetLength(previous);', 'undo.SetLength(newLength);'),
    ('MappedOutputFile.cs', 'F52', 'the mapping hook of Grow is ignored', 'map = options?.Mapper is { } mapper ? mapper(Path) : MemoryMappedFile.CreateFromFile(Path, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite);', 'map = MemoryMappedFile.CreateFromFile(Path, FileMode.Open, null, 0, MemoryMappedFileAccess.ReadWrite);'),
    ('Dcs/DifferenceCover.cs', 'F53', 'the 40 rounds of the small periods end one period earlier', 'if (x <= 256) rounds = 40;', 'if (x <= 255) rounds = 40;'),
    ('Dcs/DifferenceCover.cs', 'F54', 'the large periods get 7 rounds', 'int rounds = 6;', 'int rounds = 7;'),
    ('Dcs/DifferenceCover.cs', 'F55', 'the run of the start set is one short', 'for (int j = 0; j < m; j++) start.Add(j % x);', 'for (int j = 0; j < m - 1; j++) start.Add(j % x);'),
    ('Dcs/DifferenceCover.cs', 'F56', 'a later round of the same size replaces the best', 'if (d.Count < best.Count) best = d;', 'if (d.Count <= best.Count) best = d;'),
    ('Ph6/BwtPass.cs', 'F57', 'a window of one byte is refused by the BWT pass', 'ArgumentOutOfRangeException.ThrowIfLessThan(options.WindowBytes, 1);', 'ArgumentOutOfRangeException.ThrowIfLessThan(options.WindowBytes, 2);'),
    ('Ph6/BwtPass.cs', 'F58', 'a read of exactly four bytes is refused by the BWT pass', 'ArgumentOutOfRangeException.ThrowIfLessThan(options.ReadBytes, sizeof(uint));', 'ArgumentOutOfRangeException.ThrowIfLessThan(options.ReadBytes, sizeof(uint) + 1);'),
    ('Sesa/SesaBuilder.cs', 'F59', 'the BWT default read is two mebibytes', 'ReadBytes = tune.ReadBytes ?? 1 << 20,\n            WindowBytes = window,', 'ReadBytes = tune.ReadBytes ?? 1 << 21,\n            WindowBytes = window,'),
    ('Sort/BucketSorter.cs', 'F60', 'the positions of the scans of the partitions are not counted', '                work.ScanPositions += m;\n', ''),
    ('Sort/BucketSorter.cs', 'F61', 'the chunk of the partition is rented one bit short', 'chunk = ArrayPool<uint>.Shared.Rent(1 << 14);', 'chunk = ArrayPool<uint>.Shared.Rent(1 << 13);'),
    ('Sa/SaBuilder.cs', 'F62', 'the PH3 figure of the scratch peak is that of PH1', 'Ph23ScratchPeakBytes = io2.ScratchPeakBytes,', 'Ph23ScratchPeakBytes = ph1ScratchPeak,'),
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
    json.dump(ms, open(os.path.join(HERE, '1b9c_f_%s_listed.json' % os.path.splitext(os.path.basename(f))[0].lower()), 'w'), indent=1)
print(len(M), 'fresh mutants in', len(by_file), 'files;', bad, 'patterns to fix')
