#!/usr/bin/env python3
"""Writes the hand-made mutant lists of the new code of E6 (the driver: Pipeline/EndToEnd.cs, Pipeline/Prepare.cs, the scale rules in Pipeline/EndToEndOptions.cs, Experiments/E6.cs) as e6_<name>_mutants.json
   (id, desc, old, new[, also]): old must occur exactly once in the file; also = further (old, new) pairs applied together (a mutant that needs a helper class). The stage code is not mutated here (audited with E3b-E5)."""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '../../../../..'))
APP = os.path.join(ROOT, 'apps/SuffixTree.EmProto')

def m(i, desc, old, new, also=None):
    d = {'id': i, 'desc': desc, 'old': old, 'new': new}
    if also: d['also'] = also
    return d

END_TO_END = [
    m('P01', 'the statistics file is written after the sync phase instead of before it',
      "            File.WriteAllText(files.Stats, StatsText(result));\n            if (o.SyncOutputs) Phase(o, phases, dirty, measure, PhaseSync, () => SyncAll(files));",
      "            if (o.SyncOutputs) Phase(o, phases, dirty, measure, PhaseSync, () => SyncAll(files));\n            File.WriteAllText(files.Stats, StatsText(result));"),
    m('P02', 'PH1 is given the bucket capacity of PH3', "BucketCapacity = p.DcsBucketCapacity,", "BucketCapacity = p.BucketCapacity,"),
    m('P03', 'PH1 is given the sort capacity as its base case', "BaseCase = p.DcsBaseCase,", "BaseCase = p.DcsSortCapacity,"),
    m('P04', 'PH1 is given the chunk of PH4 as its chunk slots', "ChunkSlots = p.DcsChunkSlots,", "ChunkSlots = p.LcpChunk,"),
    m('P05', 'PH3 is given the bucket capacity of PH1', "LeafCapacity = p.BucketCapacity,", "LeafCapacity = p.DcsBucketCapacity,"),
    m('P06', 'PH4 is given the chunk slots of PH1 as its chunk', "ChunkPositions = p.LcpChunk,", "ChunkPositions = p.DcsChunkSlots,"),
    m('P07', 'the thread count is ignored by PH3', "Threads = o.Threads,\n                    Seed = 1,\n                    Oversample = p.Oversample,", "Threads = 1,\n                    Seed = 1,\n                    Oversample = p.Oversample,"),
    m('P08', 'the thread count is ignored by PH1', "Threads = o.Threads,\n                    Seed = 1,\n                    FanOut = p.DcsFanOut,", "Threads = 1,\n                    Seed = 1,\n                    FanOut = p.DcsFanOut,"),
    m('P09', 'the thread count is ignored by PH4', "ChunkPositions = p.LcpChunk,\n                    Threads = o.Threads,", "ChunkPositions = p.LcpChunk,\n                    Threads = 1,"),
    m('P10', 'PH4 runs in the bounded mode', "Mode = LcpMode.Distributed,", "Mode = LcpMode.Bounded,"),
    m('P11', 'the codec option is ignored (always tagged)', "Codec = o.Codec,", "Codec = ChildCodec.Tagged,"),
    m('P12', 'the scratch directory is not deleted after a successful run', "        DeleteScratch(full);\n        r.ScratchDirectoryExists", "        r.ScratchDirectoryExists"),
    m('P13', 'a failing run leaves its scratch directory', "        catch\n        {\n            DeleteScratch(full);\n            throw;\n        }", "        catch\n        {\n            throw;\n        }"),
    m('P14', 'the rank file is not deleted after PH3', "                File.Delete(ranksPath); // the rank file has served its only consumer\n", ""),
    m('P15', 'the tracker is not told that the rank file is gone', "                ScratchTracker.Remove(rankBytes);\n", ""),
    m('P16', 'the tracker is not told about the rank file', "                ScratchTracker.Add(new FileInfo(ranksPath).Length); // the rank file is scratch of the pipeline from here on\n", ""),
    m('P17', 'the suffix array file has n entries instead of n + 1', "new MappedSaSink(files.Sa, n + 1)", "new MappedSaSink(files.Sa, n)"),
    m('P18', 'the final-length check of the packed text forgets the padding', "Expect(r.Files.Packed, PackedText.TotalBytes(n, r.Prepare.Bits), \"the packed text\");", "Expect(r.Files.Packed, PackedText.DataBytes(n, r.Prepare.Bits), \"the packed text\");"),
    m('P19', 'the metadata file is not written', "                File.WriteAllText(files.PackedMeta, prep.MetaText(Path.GetFileName(o.InputPath)));\n", ""),
    m('P20', 'the sync is skipped although it was asked for', "if (o.SyncOutputs) Phase(", "if (false && o.SyncOutputs) Phase("),
    m('P21', 'PH1 gets the arena cap of PH5', "ArenaBudget dcsBudget = Budget(o.ArenaCapBytes);", "ArenaBudget dcsBudget = Budget(o.ChildArenaCapBytes);"),
    m('P22', 'the CHILD windows are not flushed when released', "FlushWindows = true,", "FlushWindows = false,"),
    m('P23', 'PH5 is not attached to the hook of PH4', "files.Overflow, lcp, scan);", "files.Overflow, lcp, null);"),
    m('P24', 'the hook consumer is not called for the last chunk',
      "files.Overflow, lcp, scan);", "files.Overflow, lcp, new LastChunkDropper(scan, n));",
      [["    internal static string StatsText(", "    private sealed class LastChunkDropper(ILcpChunkConsumer inner, long n) : ILcpChunkConsumer\n    {\n        public void Consume(long firstRank, ReadOnlySpan<uint> lcp)\n        {\n            if (firstRank + lcp.Length <= n) inner.Consume(firstRank, lcp);\n        }\n    }\n\n    internal static string StatsText("]]),
    m('P25', 'off by one at the chunk boundary handed to PH5: every chunk loses its last value',
      "files.Overflow, lcp, scan);", "files.Overflow, lcp, new ShortChunk(scan));",
      [["    internal static string StatsText(", "    private sealed class ShortChunk(ILcpChunkConsumer inner) : ILcpChunkConsumer\n    {\n        public void Consume(long firstRank, ReadOnlySpan<uint> lcp) => inner.Consume(firstRank, lcp.Length > 1 ? lcp[..^1] : lcp);\n    }\n\n    internal static string StatsText("]]),
    m('P26', 'NODE_COUNT is written as the number of internal nodes', "$\"nodecount={r.Scan.NodeCount}\\ninternal", "$\"nodecount={r.Scan.InternalCount}\\ninternal"),
    m('P27', 'LRS_DEPTH and LRS_RANK are swapped in the statistics file', "lrsdepth={r.Scan.LrsDepth}\\nlrsrank={r.Scan.LrsRank}", "lrsdepth={r.Scan.LrsRank}\\nlrsrank={r.Scan.LrsDepth}"),
    m('P28', 'the scratch tracker is not reset at the start of a run', "        ScratchTracker.Reset();\n        var budgets", "        var budgets"),
    m('P29', 'the exact period and the Singer period are swapped', "o.ExactPeriod ? DifferenceCover.Create(o.Period) : DifferenceCover.Create(DifferenceCover.NearestSingerPeriod(o.Period))", "o.ExactPeriod ? DifferenceCover.Create(DifferenceCover.NearestSingerPeriod(o.Period)) : DifferenceCover.Create(o.Period)"),
    m('P30', 'a non-empty scratch directory is accepted', "if (Directory.Exists(full) && Directory.EnumerateFileSystemEntries(full).Any()) throw", "if (false && Directory.Exists(full) && Directory.EnumerateFileSystemEntries(full).Any()) throw"),
    m('P31', 'an empty input is accepted', "if (n == 0) throw new ArgumentException(\"the input text is empty\", nameof(o));", ""),
    m('P32', 'the stage time includes the sync', "Phases.Where(p => p.Name != EndToEnd.PhaseSync).Sum(p => p.Seconds)", "Phases.Sum(p => p.Seconds)"),
    m('P33', 'the SA stage includes PH4', "p.Name is EndToEnd.PhasePrepare or EndToEnd.PhaseDcs or EndToEnd.PhaseSa).Sum", "p.Name is EndToEnd.PhasePrepare or EndToEnd.PhaseDcs or EndToEnd.PhaseSa or EndToEnd.PhaseLcpChild).Sum"),
    m('P34', 'the RssAnon peak of the run is not reported', "r.RssAnonPeakMb = all.MaxRssAnonMb;", "r.RssAnonPeakMb = 0;"),
    m('P35', 'PH4 workers are not passed on', "Workers = p.LcpWorkers,", ""),
    m('P36', 'the CHILD windows use the PH4 window size', "WindowBytes = p.ChildWindowBytes,", "WindowBytes = p.WindowBytes,"),
    m('P38', 'the input file is not checked for existence', "if (!File.Exists(o.InputPath)) throw new FileNotFoundException(\"the input text does not exist\", o.InputPath);", ""),
    m('P39', 'the scratch of PH5 is the scratch of PH4', "ScratchDir = dirPh5,", "ScratchDir = dirPh4,"),
]

PREPARE = [
    m('R01', 'the symbol counts are not computed', "foreach (byte b in block.AsSpan(0, take)) counts[b]++;", "foreach (byte b in block.AsSpan(0, take)) { }"),
    m('R02', 'the hash misses the last byte of every block', "hash.AppendData(block.AsSpan(0, take));", "hash.AppendData(block.AsSpan(0, take - 1));"),
    m('R03', 'the width is computed for sigma - 1', "int bits = Math.Max(1, PackedText.BitsFor(letters.Count));", "int bits = Math.Max(1, PackedText.BitsFor(letters.Count - 1));"),
    m('R04', 'the dense symbols follow decreasing byte order', "for (int v = 0; v < 256; v++)\n        {\n            if (counts[v] == 0) continue;", "for (int v = 255; v >= 0; v--)\n        {\n            if (counts[v] == 0) continue;"),
    m('R05', 'a block that is not a multiple of 8 is accepted', "blockSymbols % BlockMultiple != 0", "false"),
    m('R06', 'a block of zero symbols is accepted', "blockSymbols < BlockMultiple ||", "blockSymbols < 0 ||"),
    m('R07', 'the second pass does not compare its counts with the first', "if (seen[s] != counts[letters[s]]) throw", "if (false) throw"),
    m('R08', 'the second pass compares the wrong counts', "seen[s] != counts[letters[s]]", "seen[s] != counts[s]"),
    m('R09', 'the hook between the passes is not called', "        betweenPasses?.Invoke();\n", ""),
    m('R10', 'the position check always fails', "if (writer.Position != dataBytes) throw", "if (writer.Position == dataBytes) throw"),
    m('R11', 'the last partial byte of a block is truncated', "int bytes = (int)PackedText.DataBytes(take, bits);", "int bytes = (int)((long)take * bits / 8);"),
    m('R12', 'the packed file has no padding', "PackedText.TotalBytes(n, bits), windowBytes)", "PackedText.DataBytes(n, bits), windowBytes)"),
    m('R13', 'the letters are not in the metadata', "letters_hex={Convert.ToHexString(Letters)}\\n", ""),
    m('R14', 'the metadata gives the width as sigma', "sigma={Sigma}\\n", "sigma={Bits}\\n"),
    m('R15', 'the second pass reads one byte less per block', "input.ReadExactAt(block.AsSpan(0, take), at);\n                for (int i = 0; i < take; i++)", "input.ReadExactAt(block.AsSpan(0, take - 1), at);\n                for (int i = 0; i < take; i++)"),
    m('R16', 'the sigma-th letter is lost: the counts list drops the last symbol', "[.. letters.Select(l => counts[l])]", "[.. letters.Select(l => counts[l]).SkipLast(1)]"),
    m('R17', 'the number of symbols is the length of the file minus one', "long n = new FileInfo(inputPath).Length;", "long n = new FileInfo(inputPath).Length - 1;"),
]

SCALE = [
    m('S01', 'the references are truncated instead of rounded', "Math.Max(floor, (long)Math.Round(reference * f))", "Math.Max(floor, (long)(reference * f))"),
    m('S02', 'the factor is n / (ScaleN + 1)', "double f = scaled ? n / (double)o.ScaleN : 1.0;", "double f = scaled ? n / (double)(o.ScaleN + 1) : 1.0;"),
    m('S03', 'R of PH4 is rounded down', "(long)Math.Ceiling(E4.ChunkRef * f)", "(long)Math.Floor(E4.ChunkRef * f)"),
    m('S04', 'the reference K of PH3 is 2^20', "(long)((decimal)n * BucketReference / o.ScaleN)", "(long)((decimal)n * BucketUnscaled / o.ScaleN)"),
    m('S05', 'K has a floor of one', "Math.Max(2, t.BucketCapacity", "Math.Max(1, t.BucketCapacity"),
    m('S06', 'the B2 workers are not limited by the threads', "Math.Clamp(o.ArenaCapBytes / Math.Max(1, perTask), 1, o.Threads)", "Math.Clamp(o.ArenaCapBytes / Math.Max(1, perTask), 1, 64)"),
    m('S07', 'the B2 workers are computed from the arena cap of PH5', "Math.Clamp(o.ArenaCapBytes / Math.Max(1, perTask), 1, o.Threads)", "Math.Clamp(o.ChildArenaCapBytes / Math.Max(1, perTask), 1, o.Threads)"),
    m('S08', 'an explicit R of PH4 is ignored', "long chunk = t.LcpChunk ?? (scaled", "long chunk = (scaled"),
    m('S09', 'the run spill has a floor of 3', "Scale(RunSpillReference, 0)", "Scale(RunSpillReference, 3)"),
    m('S10', 'N_base has no floor of 4', "Math.Max(4, t.DcsBaseCase ?? Scale(E3b.BaseRef, 4))", "Math.Max(1, t.DcsBaseCase ?? Scale(E3b.BaseRef, 1))"),
    m('S11', 'unscaled runs use the scale rule with a factor of one on ScaleN 0 (division by zero)', "bool scaled = o.ScaleN > 0;", "bool scaled = o.ScaleN >= 0;"),
    m('S12', 'the PH1 fan-out default is 64', "t.DcsFanOut ?? 1024", "t.DcsFanOut ?? 64"),
    m('S13', 'the PH3 fan-out default is 16', "FanOut = t.FanOut ?? 64,", "FanOut = t.FanOut ?? 16,"),
    m('S14', 'the stack block default is 64 K entries', "t.BlockEntries ?? 256 * 1024", "t.BlockEntries ?? 64 * 1024"),
    m('S15', 'the unscaled K of PH3 is the reference of the scaled one', "(scaled ? (long)((decimal)n * BucketReference / o.ScaleN) : BucketUnscaled)", "(scaled ? (long)((decimal)n * BucketReference / o.ScaleN) : BucketReference)"),
    m('S16', 'the windows of CHILD default to 8 slots of 64 KiB', "t.ChildWindowBytes ?? 1 << 20", "t.ChildWindowBytes ?? 1 << 16"),
    m('S17', 'the spill flush default is 1 MiB', "t.SpillFlushBytes ?? 16L << 20", "t.SpillFlushBytes ?? 1L << 20"),
    m('S18', 'the PH1 sort capacity is scaled from the PH1 chunk reference', "Scale(E3b.SortKRef, 3)", "Scale(E3b.ChunkRef, 3)"),
    m('S19', 'the PH1 bucket capacity is scaled from the base-case reference', "Scale(E3b.KsRef, 3)", "Scale(E3b.BaseRef, 3)"),
    m('S20', 'the reservoir default is 1024', "t.DcsReservoir ?? 16384", "t.DcsReservoir ?? 1024"),
]

E6CLI = [
    m('C01', 'the default thread count is one', 'Threads = o.Int("threads", 4),', 'Threads = o.Int("threads", 1),'),
    m('C02', 'the sync is off by default and on with --no-sync', 'SyncOutputs = !o.Flag("no-sync"),', 'SyncOutputs = o.Flag("no-sync"),'),
    m('C03', 'the scale-down is on by default', 'ScaleN = o.Long("scale-n", 0),', 'ScaleN = o.Long("scale-n", 3_000_000_000),'),
    m('C04', 'the default codec is tagged', 'o.Str("codec", "packed")', 'o.Str("codec", "tagged")'),
    m('C05', '--spill-kb is read in thousands of bytes', 'SpillBytes = Kb("spill-kb"),', 'SpillBytes = I("spill-kb") is { } sk ? sk * 1000 : null,'),
    m('C06', '--chunk sets the chunk slots of PH1', 'LcpChunk = L("chunk"),', 'LcpChunk = L("dcs-chunk"),'),
    m('C07', '--k-s sets the bucket capacity of PH3', 'DcsBucketCapacity = L("k-s"),', 'DcsBucketCapacity = L("k"),'),
    m('C08', 'the PH5 arena cap is read in KiB', 'o.Long("ph5-budget-mb", 64) << 20', 'o.Long("ph5-budget-mb", 64) << 10'),
    m('C09', 'the packed text is not compared with its reference', '        Compare("pk", res.Files.Packed, o.Get("verify-pk"));\n', ''),
    m('C10', 'the suffix array is not compared with its reference', '        Compare("sa", res.Files.Sa, o.Get("verify-sa"));\n', ''),
    m('C11', 'a difference does not fail the run', 'ok &= equal;\n        }\n        Compare("pk"', 'if (false) ok &= equal;\n        }\n        Compare("pk"'),
    m('C12', 'LRS_RANK is not compared', 'new[] { "nodecount", "internal", "lrsdepth", "lrsrank" }.All', 'new[] { "nodecount", "internal", "lrsdepth" }.All'),
    m('C13', 'the overflow list is not compared', '            Compare("ovf", res.Files.Overflow, lcp + ".ovf");\n', ''),
    m('C14', 'the exit code needs only one of the RssAnon cap and a clean run', 'bool ok = res.RssAnonPeakMb <= o.Double("rss-cap-mb", RssCapMb) && res.Clean;', 'bool ok = res.RssAnonPeakMb <= o.Double("rss-cap-mb", RssCapMb) || res.Clean;'),
    m('C15', 'the default RssAnon cap is 25,600 MB', 'internal const double RssCapMb = 256;', 'internal const double RssCapMb = 25600;'),
    m('C16', 'the CHILD bytes are 4 n instead of 4 (n + 1)', 'long child = 4 * (n + 1);', 'long child = 4 * n;'),
    m('C17', 'total_s is the wall time of the whole run', 'r.Row("E6", c, "total_s", total, "s"', 'r.Row("E6", c, "total_s", wall, "s"'),
    m('C18', 'the text file holds the dense codes instead of the letters', 'File.WriteAllBytes(outPath, Gen.Families.ToAscii(text));', 'File.WriteAllBytes(outPath, text.Symbols);'),
    m('C19', 'the hash comparison accepts any hash', 'bool equal = string.Equals(hash, expectedHash, StringComparison.OrdinalIgnoreCase);', 'bool equal = true;'),
    m('C20', 'the SA stage row is the total', 'r.Row("E6", c, "sa_stage_s", res.SaStageSeconds, "s"', 'r.Row("E6", c, "sa_stage_s", res.StageSeconds, "s"'),
    m('C21', 'ns per symbol uses the wall time', 'r.Row("E6", c, "ns_per_symbol", total * 1e9 / n, "ns"', 'r.Row("E6", c, "ns_per_symbol", wall * 1e9 / n, "ns"'),
    m('C22', 'the output footprint leaves out the metadata files', 'foreach (string f in res.Files.All)\n        {\n            long len', 'foreach (string f in res.Files.Big)\n        {\n            long len'),
    m('C23', 'the packed codec is selected by --codec tagged', '"tagged" => ChildCodec.Tagged,\n            "packed" => ChildCodec.Packed,', '"tagged" => ChildCodec.Packed,\n            "packed" => ChildCodec.Packed,'),
    m('C24', 'the exact flag is not read', 'ExactPeriod = o.Flag("exact"),', 'ExactPeriod = false,'),
]

FILES = {
    'endtoend': ('Pipeline/EndToEnd.cs', END_TO_END),
    'prepare': ('Pipeline/Prepare.cs', PREPARE),
    'scalerules': ('Pipeline/EndToEndOptions.cs', SCALE),
    'e6': ('Experiments/E6.cs', E6CLI),
}

if __name__ == '__main__':
    bad = 0
    for name, (rel, ms) in FILES.items():
        src = open(os.path.join(APP, rel)).read()
        for x in ms:
            pairs = [(x['old'], x['new'])] + [tuple(a) for a in x.get('also', [])]
            for o, _ in pairs:
                if src.count(o) != 1:
                    print('NOT UNIQUE', name, x['id'], src.count(o), repr(o[:60]))
                    bad += 1
        json.dump(ms, open(os.path.join(HERE, 'e6_%s_mutants.json' % name), 'w'), indent=1)
        print(name, len(ms))
    sys.exit(1 if bad else 0)
