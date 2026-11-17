using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Verify;

/// <summary>
/// E8 of ADR-004 section 5: the checker of the outputs of <c>e6</c> (packed text, suffix array, compact LCP, child table, statistics) that needs no in-RAM reference and so works above 2^31 symbols (positions are 32-bit: n below 2^32 - 3).
/// Memory: a few buffers and the stack of the child scan (16 bytes per open interval, bounded by <see cref="VerifyOptions.StackCapEntries"/>); the permutation bitmap (n + 1 bits) is a temporary mapped file, not anonymous memory.
/// Cost (n symbols, p threads): the rank pass reads SA and LCP once and the text at random twice per rank (about 2n random reads and, per rank, the common-prefix comparison up to the cap), the scan reads LCP once and
/// the child table at the cells it expects, the hashes read every file once. See <see cref="Names"/> for the checks.
/// </summary>
internal static unsafe class OutputVerifier
{
    /// <summary>The checks in the order they are listed.</summary>
    internal static readonly string[] Names = ["files", "meta", "text", "input", "sa_permutation", "sa_order", "lcp_values", "lcp_width", "stats", "child", "sha256"];

    private static readonly string[] Extensions = ["pk", "pk.meta", "sa", "lcp", "ovf", "cld", "stats"];

    /// <summary>Collects the messages of one check (empty = the check holds).</summary>
    private sealed class Problems
    {
        private readonly System.Text.StringBuilder _text = new();

        internal bool Any => _text.Length > 0;

        internal void Add(string message)
        {
            if (_text.Length > 0) _text.Append("; ");
            _text.Append(message);
        }

        public override string ToString() => _text.ToString();
    }

    internal static VerifyResult Run(VerifyOptions o, Action<string>? progress = null)
    {
        ArgumentNullException.ThrowIfNull(o);
        var res = new VerifyResult();
        var done = new Dictionary<string, CheckResult>();
        void Log(string s) => progress?.Invoke(s);
        void Add(string name, CheckStatus status, string detail, double seconds)
        {
            done[name] = new CheckResult(name, status, detail, seconds);
            Log($"check {name}: {status} {detail} ({seconds:F2} s)");
        }
        void Skip(string name, string why)
        {
            if (!done.ContainsKey(name)) Add(name, CheckStatus.Skipped, "not run: " + why, 0);
        }
        string P(string ext) => o.Prefix + "." + ext;
        var wanted = new HashSet<string>(Names);
        if (o.InputPath is null) wanted.Remove("input");
        if (!o.Child)
        {
            wanted.Remove("stats");
            wanted.Remove("child");
        }
        if (!o.Hash) wanted.Remove("sha256");
        string scratchDir = o.ScratchDir ?? Path.GetDirectoryName(Path.GetFullPath(o.Prefix))!;

        // ---- meta and files
        long t0 = Stopwatch.GetTimestamp();
        MetaInfo? meta = null;
        string metaError = "the metadata file does not exist";
        bool metaOk = File.Exists(P("pk.meta")) && MetaFiles.TryParseMeta(P("pk.meta"), out meta, out metaError);
        if (metaOk) metaError = CheckMeta(meta!);
        Add("meta", metaOk && metaError.Length == 0 ? CheckStatus.Pass : CheckStatus.Fail, metaOk && metaError.Length == 0 ? $"n {meta!.N}, sigma {meta.Sigma}, {meta.Bits} bits" : metaError, Elapsed(t0));

        t0 = Stopwatch.GetTimestamp();
        var filesError = new Problems();
        StatsInfo? stats = null;
        long n = 0;
        int bits = 0;
        int width = 0;
        long overflow = 0;
        foreach (string e in Extensions)
        {
            if (!File.Exists(P(e))) filesError.Add($"{e} is missing");
        }
        if (!filesError.Any)
        {
            if (!metaOk || metaError.Length != 0) filesError.Add("the metadata is unusable, the expected file lengths cannot be derived");
            else if (!MetaFiles.TryParseStats(P("stats"), out stats, out string statsError)) filesError.Add("the statistics file is malformed: " + statsError);
            else
            {
                n = meta!.N;
                bits = meta.Bits;
                width = stats!.LcpWidth;
                overflow = stats.Overflow;
                if (stats.N != n) filesError.Add($"stats n {stats.N} differs from meta n {n}");
                if (stats.Sigma != meta.Sigma) filesError.Add($"stats sigma {stats.Sigma} differs from meta sigma {meta.Sigma}");
                if (width is not (1 or 2 or 4)) filesError.Add($"lcpwidth {width} is not 1, 2 or 4");
                else
                {
                    void Length(string ext, long expected)
                    {
                        long got = new FileInfo(P(ext)).Length;
                        if (got != expected) filesError.Add($"{ext} has {got} bytes, expected {expected}");
                    }
                    Length("pk", PackedText.TotalBytes(n, bits));
                    Length("sa", 4 * (n + 1));
                    Length("cld", 4 * (n + 1));
                    Length("lcp", width * (n + 1));
                    Length("ovf", 8 * overflow);
                }
            }
        }
        bool filesOk = !filesError.Any && metaOk && metaError.Length == 0;
        Add("files", filesOk ? CheckStatus.Pass : CheckStatus.Fail, filesOk ? "all seven files exist with the lengths of n, bits, lcpwidth and overflow" : filesError.ToString(), Elapsed(t0));
        if (filesOk)
        {
            res.Facts["n"] = n.ToString(CultureInfo.InvariantCulture);
            res.Facts["bits"] = bits.ToString(CultureInfo.InvariantCulture);
            res.Facts["sigma"] = meta!.Sigma.ToString(CultureInfo.InvariantCulture);
            res.Facts["lcpwidth"] = width.ToString(CultureInfo.InvariantCulture);
            res.Facts["overflow"] = overflow.ToString(CultureInfo.InvariantCulture);
            res.Facts["lce_cap"] = o.LceCap.ToString(CultureInfo.InvariantCulture);
        }
        if (filesOk && n >= uint.MaxValue - 2)
        {
            filesOk = false;
            Add("text", CheckStatus.Fail, $"n = {n} is beyond the 32-bit positions of the prototype (below 2^32 - 3)", 0);
        }

        if (filesOk)
        {
            try
            {
                RunChecks(o, res, meta!, stats!, n, bits, width, overflow, scratchDir, wanted, Add, Skip, Log);
            }
            catch (Exception ex) when (ex is IOException or InvalidDataException or UnauthorizedAccessException or ArgumentException or OutOfMemoryException)
            {
                Add("files", CheckStatus.Fail, "a file could not be read: " + ex.Message, 0);
            }
        }
        foreach (string name in Names.Where(wanted.Contains))
        {
            if (name == "sha256" && !done.ContainsKey(name) && o.Hash) HashAll(o, res, Add, Log);
            Skip(name, "the files or the metadata check failed");
        }
        foreach (string name in Names)
        {
            if (done.TryGetValue(name, out CheckResult? c) && (wanted.Contains(name) || c.Status == CheckStatus.Fail)) res.Checks.Add(c);
        }
        return res;
    }

    private static double Elapsed(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;

    /// <summary>Consistency of the metadata with itself; empty when it holds.</summary>
    private static string CheckMeta(MetaInfo m)
    {
        var e = new Problems();
        if (m.N < 1) e.Add("n is below 1");
        if (m.Bits != Math.Max(1, BitsFor(m.Sigma))) e.Add($"bits {m.Bits} is not the width of {m.Sigma} symbols");
        for (int i = 1; i < m.Letters.Length; i++)
        {
            if (m.Letters[i] <= m.Letters[i - 1]) e.Add("letters are not strictly increasing");
        }
        if (m.Counts.Sum() != m.N) e.Add($"counts add up to {m.Counts.Sum()}, not n");
        if (m.Counts.Any(c => c < 1)) e.Add("a letter of the alphabet does not occur");
        return e.ToString();
    }

    private static int BitsFor(int sigma)
    {
        int b = 1;
        while ((1 << b) < sigma) b++;
        return b;
    }

    private static void RunChecks(VerifyOptions o, VerifyResult res, MetaInfo meta, StatsInfo stats, long n, int bits, int width, long overflow, string scratchDir,
        HashSet<string> wanted, Action<string, CheckStatus, string, double> add, Action<string, string> skip, Action<string> log)
    {
        string P(string ext) => o.Prefix + "." + ext;
        using var mappedText = new MappedFile(P("pk"), writable: false);
        using var text = new PackedText(mappedText.Pointer, n, bits, mappedText, ownsOwner: false);
        // ---- text
        long t0 = Stopwatch.GetTimestamp();
        string textError = CheckText(mappedText, text, meta, o.Threads, o.ScanBlock);
        add("text", textError.Length == 0 ? CheckStatus.Pass : CheckStatus.Fail, textError.Length == 0 ? "symbol counts equal the metadata, symbols below sigma, padding and the unused bits of the last byte are zero" : textError, Elapsed(t0));

        // ---- input
        if (wanted.Contains("input"))
        {
            t0 = Stopwatch.GetTimestamp();
            string inputError = CheckInput(o.InputPath!, text, meta, o.IoBufferBytes);
            add("input", inputError.Length == 0 ? CheckStatus.Pass : CheckStatus.Fail, inputError.Length == 0 ? "the packed text equals the input symbol by symbol and its SHA-256 equals the metadata's" : inputError, Elapsed(t0));
        }

        // ---- the rank pass: permutation, order, LCP values
        using var mappedSa = new MappedFile(P("sa"), writable: false);
        using var mappedLcp = new MappedFile(P("lcp"), writable: false);
        MappedFile? mappedOvf = overflow > 0 ? new MappedFile(P("ovf"), writable: false) : null;
        try
        {
            var lcp = new LcpView(mappedLcp.Pointer, width, mappedOvf is null ? null : (uint*)mappedOvf.Pointer, overflow);
            string bitmapPath = Path.Combine(scratchDir, Path.GetFileName(o.Prefix) + ".verify.bits");
            PassTotals totals;
            t0 = Stopwatch.GetTimestamp();
            try
            {
                totals = RankPass(o, text, (uint*)mappedSa.Pointer, lcp, n, bitmapPath, log);
            }
            finally
            {
                if (File.Exists(bitmapPath)) File.Delete(bitmapPath);
            }
            double passSeconds = Elapsed(t0);
            res.Facts["long_pairs"] = totals.LongPairs.ToString(CultureInfo.InvariantCulture);
            res.Facts["lce_steps"] = totals.Steps.ToString(CultureInfo.InvariantCulture);
            add("sa_permutation", totals.Permutation.Count == 0 ? CheckStatus.Pass : CheckStatus.Fail,
                totals.Permutation.Count == 0 ? $"SA[0] = n and the {n + 1} entries are distinct and at most n" : $"{totals.Permutation.Count} violations, first {totals.Permutation.First}", passSeconds);
            add("sa_order", totals.Order.Count == 0 ? CheckStatus.Pass : CheckStatus.Fail,
                totals.Order.Count == 0 ? $"all {n} adjacent pairs of suffixes are in strict order ({totals.LongPairs} long pairs checked at their ends and by probes)" : $"{totals.Order.Count} violations, first {totals.Order.First}", passSeconds);
            add("lcp_values", totals.Lcp.Count == 0 ? CheckStatus.Pass : CheckStatus.Fail,
                totals.Lcp.Count == 0 ? $"every stored LCP equals the common prefix of its adjacent pair ({totals.LongPairs} long pairs above the cap {o.LceCap}), the overflow list matches the escapes" : $"{totals.Lcp.Count} violations, first {totals.Lcp.First}", passSeconds);
            if (totals.Lcp.Count == 0)
            {
                LcpKind rule = Rule(n, totals.Escapes8, totals.Escapes16);
                int ruleWidth = rule switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 4 };
                add("lcp_width", ruleWidth == width ? CheckStatus.Pass : CheckStatus.Fail,
                    ruleWidth == width ? $"width {width} is the rule's choice (values >= 255: {totals.Escapes8}, >= 65535: {totals.Escapes16})" : $"width {width}, the rule gives {ruleWidth} (values >= 255: {totals.Escapes8}, >= 65535: {totals.Escapes16})", 0);
            }
            else
            {
                skip("lcp_width", "lcp_values failed");
            }
            bool structure = totals.Structure.Count == 0;

            // ---- the scan: statistics and child table
            if (wanted.Contains("child"))
            {
                if (!structure)
                {
                    skip("stats", "the overflow list of the LCP is inconsistent");
                    skip("child", "the overflow list of the LCP is inconsistent");
                }
                else
                {
                    using var mappedCld = new MappedFile(P("cld"), writable: false);
                    t0 = Stopwatch.GetTimestamp();
                    ScanResult scan = Scan(lcp, (uint*)mappedCld.Pointer, n, o.StackCapEntries, o.Threads, o.ScanBlock, log);
                    double scanSeconds = Elapsed(t0);
                    res.Facts["stack_peak"] = scan.StackPeak.ToString(CultureInfo.InvariantCulture);
                    res.Facts["internal"] = scan.Internal.ToString(CultureInfo.InvariantCulture);
                    res.Facts["maxlcp"] = scan.MaxLcp.ToString(CultureInfo.InvariantCulture);
                    if (scan.Aborted.Length > 0)
                    {
                        skip("stats", scan.Aborted);
                        skip("child", scan.Aborted);
                    }
                    else
                    {
                        var expected = new StatsInfo(1 + scan.Internal + n + 1, scan.Internal, scan.LrsDepth, scan.LrsRank, n, meta.Sigma, width, overflow, scan.MaxLcp);
                        string expectedText = MetaFiles.StatsText(expected);
                        string found = MetaFiles.ReadText(P("stats"));
                        string statsDetail = "";
                        if (expectedText != found) statsDetail = MetaFiles.FirstDifference(expectedText, found);
                        else if (scan.LrsDepth != scan.MaxLcp || scan.LrsRank != (scan.MaxLcp == 0 ? 0 : scan.FirstMaxRank - 1)) statsDetail = $"the two formulations of LRS disagree (depth {scan.LrsDepth} / max LCP {scan.MaxLcp}, rank {scan.LrsRank} / first max at {scan.FirstMaxRank})";
                        res.Facts["lrsdepth"] = scan.LrsDepth.ToString(CultureInfo.InvariantCulture);
                        res.Facts["lrsrank"] = scan.LrsRank.ToString(CultureInfo.InvariantCulture);
                        add("stats", statsDetail.Length == 0 ? CheckStatus.Pass : CheckStatus.Fail,
                            statsDetail.Length == 0 ? $"NODE_COUNT {expected.NodeCount}, internal {scan.Internal}, LRS {scan.LrsDepth} at rank {scan.LrsRank}, max LCP {scan.MaxLcp}: recomputed from the LCP, the file is byte-identical to it" : statsDetail, scanSeconds);
                        string childDetail = "";
                        if (scan.Mismatches != 0 || scan.NonZeroCells != scan.ExpectedWrites)
                        {
                            childDetail = $"{scan.Mismatches} cells differ from the definition, {scan.NonZeroCells} non-zero cells against {scan.ExpectedWrites} expected";
                            if (scan.FirstMismatch.Length > 0) childDetail += ", first: " + scan.FirstMismatch;
                        }
                        add("child", childDetail.Length == 0 ? CheckStatus.Pass : CheckStatus.Fail,
                            childDetail.Length == 0 ? $"all {scan.ExpectedWrites} written cells equal the lcp-interval definition (next links, up and down fields) and no other cell is set; stack peak {scan.StackPeak}" : childDetail, scanSeconds);
                    }
                }
            }
        }
        finally
        {
            mappedOvf?.Dispose();
        }
        if (wanted.Contains("sha256")) HashAll(o, res, add, log);
    }

    private static LcpKind Rule(long n, long escapes8, long escapes16)
    {
        if (16 * escapes8 <= n) return LcpKind.U8;
        if (8 * escapes16 <= n) return LcpKind.U16;
        return LcpKind.U32;
    }

    private enum LcpKind
    {
        U8,
        U16,
        U32,
    }

    // ---------------------------------------------------------------- text

    private static string CheckText(MappedFile file, PackedText text, MetaInfo meta, int threads, long block)
    {
        long n = text.Length;
        int bits = text.Bits;
        var total = new long[256];
        Parallel.For(0, (int)((n + block - 1) / block), new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, threads) }, () => new long[256], (b, _, local) =>
        {
            long end = Math.Min(n, (b + 1) * block);
            for (long i = b * block; i < end; i++) local[text.Symbol(i)]++;
            return local;
        }, local =>
        {
            lock (total) for (int s = 0; s < 256; s++) total[s] += local[s];
        });
        var error = new Problems();
        for (int s = 0; s < 256; s++)
        {
            long expected = s < meta.Sigma ? meta.Counts[s] : 0;
            if (total[s] != expected) error.Add($"symbol {s} occurs {total[s]} times, the metadata says {expected}");
        }
        long dataBytes = PackedText.DataBytes(n, bits);
        for (long i = dataBytes; i < file.Bytes; i++)
        {
            if (file.Pointer[i] != 0)
            {
                error.Add($"padding byte {i - dataBytes} is not zero");
                break;
            }
        }
        int usedBits = (int)((n * bits) & 7);
        if (usedBits != 0 && (file.Pointer[dataBytes - 1] & ((1 << (8 - usedBits)) - 1)) != 0) error.Add("the unused bits of the last data byte are not zero");
        return error.ToString();
    }

    private static string CheckInput(string path, PackedText text, MetaInfo meta, int bufferBytes)
    {
        if (!File.Exists(path)) return "the input file does not exist";
        using var handle = File.OpenHandle(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite, FileOptions.SequentialScan);
        long length = RandomAccess.GetLength(handle);
        if (length != text.Length) return $"the input has {length} bytes, n is {text.Length}";
        var map = new int[256];
        Array.Fill(map, -1);
        for (int s = 0; s < meta.Sigma; s++) map[meta.Letters[s]] = s;
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        var buffer = new byte[bufferBytes];
        long bad = 0, first = -1;
        for (long at = 0; at < length; at += buffer.Length)
        {
            int take = (int)Math.Min(buffer.Length, length - at);
            int got = 0;
            while (got < take)
            {
                int r = RandomAccess.Read(handle, buffer.AsSpan(got, take - got), at + got);
                if (r <= 0) return "the input could not be read to its end";
                got += r;
            }
            hash.AppendData(buffer.AsSpan(0, take));
            for (int i = 0; i < take; i++)
            {
                if (map[buffer[i]] != (int)text.Symbol(at + i))
                {
                    bad++;
                    if (first < 0) first = at + i;
                }
            }
        }
        var error = new Problems();
        if (bad > 0) error.Add($"{bad} symbols of the packed text differ from the input, first at {first}");
        string sha = Convert.ToHexStringLower(hash.GetHashAndReset());
        if (sha != meta.Sha256) error.Add($"the SHA-256 of the input is {sha}, the metadata says {meta.Sha256}");
        return error.ToString();
    }

    // ---------------------------------------------------------------- the rank pass

    /// <summary>Count and first occurrence of one kind of violation (the lowest rank over all threads).</summary>
    private sealed class Violations
    {
        private readonly object _gate = new();
        private long _count;
        private long _rank = long.MaxValue;

        internal long Count => Interlocked.Read(ref _count);

        internal string First { get; private set; } = "";

        internal void Add(long rank, string detail)
        {
            Interlocked.Increment(ref _count);
            lock (_gate)
            {
                if (rank < _rank)
                {
                    _rank = rank;
                    First = $"at rank {rank}: {detail}";
                }
            }
        }
    }

    private sealed class PassTotals
    {
        internal Violations Permutation { get; } = new();

        internal Violations Order { get; } = new();

        internal Violations Lcp { get; } = new();

        internal Violations Structure { get; } = new();

        internal long LongPairs;
        internal long Steps;
        internal long Escapes8;
        internal long Escapes16;
    }

    private sealed class PassLocal
    {
        internal long LongPairs;
        internal long Steps;
        internal long Escapes8;
        internal long Escapes16;
    }

    private static PassTotals RankPass(VerifyOptions o, PackedText text, uint* sa, LcpView lcp, long n, string bitmapPath, Action<string> log)
    {
        long ranks = n + 1;
        long blockRanks = Math.Max(1, o.BlockRanks);
        long blocks = (ranks + blockRanks - 1) / blockRanks;
        long bitmapBytes = ((ranks + 63) / 64) * 8;
        using (var fs = new FileStream(bitmapPath, FileMode.Create, FileAccess.Write, FileShare.None)) fs.SetLength(bitmapBytes);
        using var bitmap = new MappedFile(bitmapPath, writable: true);
        ulong* bits = (ulong*)bitmap.Pointer;
        var totals = new PassTotals();
        long finished = 0;
        long step = Math.Max(1, blocks / 20);
        Parallel.For(0L, blocks, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, o.Threads) }, () => new PassLocal(), (b, _, local) =>
        {
            RankBlock(o, text, sa, lcp, n, bits, b * blockRanks, Math.Min(ranks, (b + 1) * blockRanks), totals, local);
            long f = Interlocked.Increment(ref finished);
            if (f % step == 0) log($"rank pass: {100.0 * f / blocks:F0} % of the blocks");
            return local;
        }, local =>
        {
            Interlocked.Add(ref totals.LongPairs, local.LongPairs);
            Interlocked.Add(ref totals.Steps, local.Steps);
            Interlocked.Add(ref totals.Escapes8, local.Escapes8);
            Interlocked.Add(ref totals.Escapes16, local.Escapes16);
        });
        return totals;
    }

    private static void RankBlock(VerifyOptions o, PackedText text, uint* sa, LcpView lcp, long n, ulong* bits, long first, long end, PassTotals totals, PassLocal local)
    {
        long cap = o.LceCap;
        long pair = lcp.HasEscape && lcp.PairCount > 0 ? lcp.LowerBound(first) : 0;
        for (long i = first; i < end; i++)
        {
            uint q = sa[i];
            // permutation: the entry is in range and new
            if (q > n) totals.Permutation.Add(i, $"SA[{i}] = {q} exceeds n = {n}");
            else
            {
                ulong mask = 1UL << (int)(q & 63);
                ulong old = Interlocked.Or(ref bits[q >> 6], mask);
                if ((old & mask) != 0) totals.Permutation.Add(i, $"SA[{i}] = {q} occurs more than once");
            }
            // the stored LCP, with its overflow entry
            uint stored = lcp.BaseAt(i);
            if (lcp.HasEscape && stored == lcp.Escape)
            {
                if (pair < lcp.PairCount && lcp.RankOfPair(pair) == i && lcp.ValueOfPair(pair) >= lcp.Escape)
                {
                    stored = lcp.ValueOfPair(pair);
                    pair++;
                }
                else
                {
                    totals.Structure.Add(i, "an escape entry without its overflow pair, or the pair is out of order or below the escape value");
                    totals.Lcp.Add(i, $"the escape entry at rank {i} has no matching overflow pair");
                    continue;
                }
            }
            if (stored >= 255) local.Escapes8++;
            if (stored >= 65535) local.Escapes16++;
            if (i == 0)
            {
                if (q != n) totals.Permutation.Add(0, $"SA[0] = {q}, expected n = {n} (the empty suffix)");
                if (stored != 0) totals.Lcp.Add(0, $"LCP[0] = {stored}, expected 0");
                continue;
            }
            uint p = sa[i - 1];
            if (p > n || q > n || p == q) continue; // reported by the permutation check
            long maxLen = n - Math.Max(p, q);
            long lim = cap > 0 ? Math.Min(maxLen, cap) : maxLen;
            long l = text.Lce(p, q, lim, ref local.Steps);
            if (l < lim || lim == maxLen)
            {
                // decided: the common prefix is exactly l (a mismatch, or one suffix ended)
                bool ordered = l == maxLen ? p > q : text.Symbol(p + l) < text.Symbol(q + l);
                if (!ordered) totals.Order.Add(i, $"suffix {p} is not smaller than suffix {q} (common prefix {l})");
                if (stored != l) totals.Lcp.Add(i, $"LCP[{i}] = {stored}, the suffixes {p} and {q} share {l} symbols");
                continue;
            }
            // a long pair: the first `cap` symbols agree; check the stored length at its ends and at a few places in between
            local.LongPairs++;
            if (stored < cap || stored > maxLen)
            {
                totals.Lcp.Add(i, $"LCP[{i}] = {stored}, but the suffixes {p} and {q} share at least {cap} and at most {maxLen} symbols");
                continue;
            }
            bool ok;
            if (stored == maxLen) ok = p > q;
            else
            {
                uint a = text.Symbol(p + stored), b = text.Symbol(q + stored);
                if (a == b)
                {
                    totals.Lcp.Add(i, $"LCP[{i}] = {stored}, but the suffixes {p} and {q} agree at offset {stored}");
                    continue;
                }
                ok = a < b;
            }
            if (!ok) totals.Order.Add(i, $"suffix {p} is not smaller than suffix {q} (stored common prefix {stored})");
            long span = stored - cap;
            for (int k = 0; k < o.Probes && span > 0; k++)
            {
                long off = k == 0 ? stored - 1 : cap + (long)(Mix(((ulong)i * 0x9E3779B97F4A7C15UL) + (ulong)k) % (ulong)span);
                if (text.Symbol(p + off) != text.Symbol(q + off))
                {
                    totals.Lcp.Add(i, $"LCP[{i}] = {stored}, but the suffixes {p} and {q} differ at offset {off}");
                    break;
                }
            }
        }
        if (lcp.HasEscape && lcp.PairCount > 0 && pair != lcp.LowerBound(end)) totals.Structure.Add(first, "an overflow pair of this block has no escape entry");
        if (lcp.HasEscape && lcp.PairCount > 0 && pair != lcp.LowerBound(end)) totals.Lcp.Add(first, $"overflow pairs in the rank range {first}..{end - 1} without an escape entry");
    }

    private static ulong Mix(ulong x)
    {
        unchecked
        {
            x += 0x9E3779B97F4A7C15UL;
            x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9UL;
            x = (x ^ (x >> 27)) * 0x94D049BB133111EBUL;
            return x ^ (x >> 31);
        }
    }

    // ---------------------------------------------------------------- the scan

    private sealed class ScanResult
    {
        internal string Aborted = "";
        internal long Internal;
        internal long MaxLcp;
        internal long FirstMaxRank;
        internal long LrsDepth;
        internal long LrsRank;
        internal long ExpectedWrites;
        internal long Mismatches;
        internal long NonZeroCells;
        internal string FirstMismatch = "";
        internal long StackPeak;
    }

    /// <summary>
    /// One pass over the LCP in rank order, the textbook bottom-up traversal of the lcp-intervals with a stack: counts the internal intervals, finds the deepest, and for every field that the layout A / packed child table must hold
    /// (the next l-index link of an l-index whose successor is not the next rank, the first l-index at the right boundary of a child that is not the last, at the left boundary of the last child) compares the cell of the file with it.
    /// Afterwards the non-zero cells of the file are counted: they must be exactly the expected ones.
    /// </summary>
    private static ScanResult Scan(LcpView lcp, uint* cld, long n, long stackCap, int threads, long block, Action<string> log)
    {
        var r = new ScanResult();
        var stack = new IntervalStack(stackCap);
        const uint None = IntervalStack.None;
        stack.TryPush(new OpenInterval { Ell = 0, Lb = 0, First = None, Last = None });
        long pair = 0;
        long step = Math.Max(1, (n + 1) / 20);
        void Expect(long cell, uint value)
        {
            r.ExpectedWrites++;
            uint got = cld[cell];
            if (got != value)
            {
                r.Mismatches++;
                if (r.FirstMismatch.Length == 0) r.FirstMismatch = $"cell {cell} holds {got}, expected {value}";
            }
        }
        void Close(uint ell, uint lb, long rb, uint first, bool lastChild)
        {
            if (ell > 0)
            {
                r.Internal++;
                if (ell > r.LrsDepth)
                {
                    r.LrsDepth = ell;
                    r.LrsRank = lb;
                }
            }
            Expect(lastChild ? lb : rb, first);
        }
        for (long i = 0; i <= n; i++)
        {
            uint cur = lcp.BaseAt(i);
            if (lcp.HasEscape && cur == lcp.Escape)
            {
                if (pair >= lcp.PairCount || lcp.RankOfPair(pair) != i)
                {
                    r.Aborted = $"the overflow list does not match the escape entry at rank {i}";
                    return r;
                }
                cur = lcp.ValueOfPair(pair);
                pair++;
            }
            if (cur > r.MaxLcp)
            {
                r.MaxLcp = cur;
                r.FirstMaxRank = i;
            }
            if (i == 0) continue;
            if (i == 1 && cur != 0)
            {
                r.Aborted = "LCP[1] is not 0";
                return r;
            }
            if (cur > n)
            {
                r.Aborted = $"LCP[{i}] = {cur} exceeds n";
                return r;
            }
            uint left = (uint)(i - 1);
            while (stack.Top.Ell > cur)
            {
                OpenInterval e = stack.Top;
                stack.Pop();
                Close(e.Ell, e.Lb, i - 1, e.First, cur < stack.Top.Ell);
                left = e.Lb;
            }
            ref OpenInterval top = ref stack.Top;
            if (top.Ell < cur)
            {
                if (!stack.TryPush(new OpenInterval { Ell = cur, Lb = left, First = (uint)i, Last = (uint)i }))
                {
                    r.Aborted = $"the stack of open intervals would exceed {stackCap} entries (LCP values as large as {r.MaxLcp}); raise --stack-cap-m";
                    r.StackPeak = stack.Peak;
                    return r;
                }
            }
            else if (top.First == None)
            {
                top.First = (uint)i;
                top.Last = (uint)i;
            }
            else
            {
                if (i != top.Last + 1L) Expect(top.Last, (uint)i);
                top.Last = (uint)i;
            }
            if (i % step == 0) log($"child scan: {100.0 * i / (n + 1):F0} % of the ranks");
        }
        while (stack.Count > 0)
        {
            OpenInterval e = stack.Top;
            stack.Pop();
            Close(e.Ell, e.Lb, n, e.First, true);
        }
        if (pair != lcp.PairCount)
        {
            r.Aborted = $"{lcp.PairCount - pair} overflow pairs were not used";
            return r;
        }
        r.StackPeak = stack.Peak;
        // the cells that are set must be exactly the expected ones
        long cells = n + 1;
        long nonZero = 0;
        Parallel.For(0, (int)((cells + block - 1) / block), new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, threads) }, b =>
        {
            long count = 0;
            long endCell = Math.Min(cells, (b + 1) * block);
            for (long k = b * block; k < endCell; k++)
            {
                if (cld[k] != 0) count++;
            }
            Interlocked.Add(ref nonZero, count);
        });
        r.NonZeroCells = nonZero;
        return r;
    }

    // ---------------------------------------------------------------- hashes

    private static void HashAll(VerifyOptions o, VerifyResult res, Action<string, CheckStatus, string, double> add, Action<string> log)
    {
        long t0 = Stopwatch.GetTimestamp();
        var found = new Dictionary<string, string>();
        var error = new Problems();
        Parallel.ForEach(Extensions, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, o.Threads) }, e =>
        {
            try
            {
                string sha = HashFile(o.Prefix + "." + e, o.IoBufferBytes);
                lock (found) found[e] = sha;
                log($"sha256 {e}: {sha[..16]}");
            }
            catch (Exception ex) when (ex is IOException or UnauthorizedAccessException)
            {
                lock (found) error.Add($"{e}: {ex.Message}");
            }
        });
        foreach (string e in Extensions)
        {
            if (found.TryGetValue(e, out string? sha)) res.Sha256[e] = sha;
        }
        add("sha256", error.Any ? CheckStatus.Fail : CheckStatus.Pass, error.Any ? error.ToString() : "SHA-256 of all seven files recorded", Elapsed(t0));
    }

    internal static string HashFile(string path, int bufferBytes = 1 << 22)
    {
        using var handle = File.OpenHandle(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite, FileOptions.SequentialScan);
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        long length = RandomAccess.GetLength(handle);
        var buffer = new byte[bufferBytes];
        for (long at = 0; at < length;)
        {
            int got = RandomAccess.Read(handle, buffer, at);
            if (got <= 0) break;
            hash.AppendData(buffer.AsSpan(0, got));
            at += got;
        }
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }
}
