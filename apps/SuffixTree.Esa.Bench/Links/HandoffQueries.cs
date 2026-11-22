using System.Buffers.Binary;
using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using Microsoft.Win32.SafeHandles;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Bench.Links;

// Work package 1b.9f: the query check of the hand-off package (`sesa handoff queries`, analysis/s1b-handoff). It lives in Links/ for the reason given in HandoffOracle.cs.

/// <summary>What the query check takes.</summary>
internal sealed record HandoffQueryOptions
{
    internal required string Container { get; init; }

    internal required string Input { get; init; }

    internal ulong Seed { get; init; } = 20261122;

    internal IReadOnlyList<int> Lengths { get; init; } = [12, 20, 32];

    /// <summary>Queries cut from the text, per length.</summary>
    internal int Present { get; init; } = 300;

    /// <summary>Random queries over the alphabet of the text, per length.</summary>
    internal int Absent { get; init; } = 200;

    /// <summary>Locate sets are compared when the count is at most this; above it the counts are compared (and the positions that were reported are not).</summary>
    internal int Cap { get; init; } = 1000;

    /// <summary><c>auto</c>: the engine answers when the container is within its 32-bit surface; <c>true</c>: it must; <c>false</c>: not at all.</summary>
    internal string Engine { get; init; } = "auto";

    /// <summary>A seam for the tests (a tree that answers wrongly on purpose); <c>null</c> in the product: <see cref="SuffixIndexFactory.Open"/>.</summary>
    internal Func<string, ISuffixTree>? OpenEngine { get; init; }
}

/// <summary>What an answerer said about one query: the count, whether it exists, and the positions (sorted, complete) when the count is within the cap.</summary>
internal sealed record QueryAnswer(long Count, bool Exists, IReadOnlyList<long>? Positions);

internal sealed record HandoffMismatch(int Index, string Answerer, string Kind, string Detail);

/// <summary>Where an answerer and the oracle disagree.</summary>
internal static class HandoffCompare
{
    internal static List<HandoffMismatch> Compare(string answerer, IReadOnlyList<HandoffQuery> queries, ScanResult oracle, IReadOnlyList<QueryAnswer> answers, int cap)
    {
        var found = new List<HandoffMismatch>();
        foreach (HandoffQuery q in queries)
        {
            QueryAnswer a = answers[q.Index];
            long truth = oracle.Counts[q.Index];
            string who = $"len {q.Bytes.Length} {(q.Present ? "cut at " + q.Origin.ToString(CultureInfo.InvariantCulture) : "absent")}";
            if (a.Count != truth) found.Add(new HandoffMismatch(q.Index, answerer, "count", $"{who}: answered {a.Count}, the text has {truth}"));
            if (a.Exists != (truth > 0)) found.Add(new HandoffMismatch(q.Index, answerer, "exists", $"{who}: answered {a.Exists}, the text has {truth} occurrences"));
            if (truth > cap) continue;
            if (a.Positions is null)
            {
                if (a.Count == truth) found.Add(new HandoffMismatch(q.Index, answerer, "locate", $"{who}: no positions answered for {truth} occurrences"));
            }
            else if (!a.Positions.SequenceEqual(oracle.Positions[q.Index]))
            {
                found.Add(new HandoffMismatch(q.Index, answerer, "locate", $"{who}: {a.Positions.Count} positions answered, {oracle.Positions[q.Index].Count} in the text, first difference at {FirstDifference(a.Positions, oracle.Positions[q.Index])}"));
            }
        }
        return found;
    }

    private static string FirstDifference(IReadOnlyList<long> a, IReadOnlyList<long> b)
    {
        for (int i = 0; i < Math.Min(a.Count, b.Count); i++)
        {
            if (a[i] != b[i]) return $"index {i} ({a[i]} against {b[i]})";
        }
        return $"index {Math.Min(a.Count, b.Count)} (one list ends)";
    }
}

/// <summary>
/// One query check: the query set from the text, the oracle's single pass over the text, the search over the container's suffix array (reads the SA and the input text; works for any container the builder makes, also
/// above the 32-bit surface of the engine), and the engine (<see cref="SuffixIndexFactory.Open"/>) where it can open the container. Rows go to <c>row(metric, value, notes)</c>; the result is 0 for PASS and 3 for FAIL.
/// </summary>
internal static class HandoffQueryRun
{
    private static string Num(long v) => v.ToString(CultureInfo.InvariantCulture);

    private static string Sec(long start) => ((Stopwatch.GetTimestamp() - start) / (double)Stopwatch.Frequency).ToString("R", CultureInfo.InvariantCulture);

    internal static int Execute(HandoffQueryOptions options, Action<string, string, string> row)
    {
        ArgumentNullException.ThrowIfNull(options);
        ArgumentNullException.ThrowIfNull(row);
        using SafeFileHandle input = File.OpenHandle(options.Input, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
        long n = RandomAccess.GetLength(input);
        long saOffset;
        long textLength;
        using (SesaFile header = SesaFile.Open(options.Container, new EsaOpenOptions()))
        {
            textLength = header.Header.TextLength;
            saOffset = header.Parsed.Sa.Offset;
        }
        row("q_input_length", $"{n} bytes in the input, {textLength} symbols in the container", n == textLength ? "equal" : "NOT EQUAL");
        if (n != textLength)
        {
            row("q_result", "FAIL", "the input file is not the text of the container");
            return 3;
        }

        // the set
        long t = Stopwatch.GetTimestamp();
        byte[] alphabet = HandoffQuerySet.AlphabetOf(ReadBytes(input, 0, (int)Math.Min(n, 1 << 20)));
        IReadOnlyList<HandoffQuery> queries = HandoffQuerySet.Build(n, (pos, len) => ReadBytes(input, pos, len), alphabet, options.Lengths, options.Present, options.Absent, options.Seed);
        row("q_total", Num(queries.Count), $"lengths {string.Join(',', options.Lengths)}, {options.Present} cut from the text and {options.Absent} random per length, seed {options.Seed}, alphabet of {alphabet.Length} letters");
        row("q_present", Num(queries.Count(q => q.Present)), "cut from the text");
        row("q_absent", Num(queries.Count(q => !q.Present)), "random strings over the alphabet");
        row("q_set_sha256", SetHash(queries), "SHA-256 of the lengths, origins and bytes of the set: equal on every machine for the same text and parameters");
        row("q_set_s", Sec(t), "seconds to cut the set");

        // the oracle: one pass over the text
        t = Stopwatch.GetTimestamp();
        ScanResult oracle;
        using (var stream = new FileStream(options.Input, FileMode.Open, FileAccess.Read, FileShare.ReadWrite, 1 << 20, FileOptions.SequentialScan)) oracle = RollingScan.Run(stream, queries, options.Cap);
        row("q_oracle_scan_s", Sec(t), $"one pass over {n} symbols, {queries.Select(q => q.Bytes.Length).Distinct().Count()} lengths");
        row("q_oracle_hits_total", Num(oracle.Counts.Sum()), "occurrences of all queries together");
        row("q_oracle_present_without_hit", Num(queries.Count(q => q.Present && oracle.Counts[q.Index] == 0)), "a query cut from the text that the scan does not find (0 unless the oracle is wrong)");
        row("q_locate_skipped_over_cap", Num(oracle.Counts.Count(c => c > options.Cap)), $"queries with more than {options.Cap} occurrences: counts compared, positions not");
        var all = new List<HandoffMismatch>();
        all.AddRange(queries.Where(q => q.Present && oracle.Counts[q.Index] == 0).Select(q => new HandoffMismatch(q.Index, "oracle", "self", $"a query cut at {q.Origin} was not found by the scan")));

        // the search over the container's suffix array
        t = Stopwatch.GetTimestamp();
        using SafeFileHandle container = File.OpenHandle(options.Container, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
        var four = new byte[4];
        long Sa(long rank)
        {
            RandomAccess.Read(container, four, saOffset + (4 * rank));
            return BinaryPrimitives.ReadUInt32LittleEndian(four);
        }

        int ReadText(long pos, byte[] buffer)
        {
            int count = (int)Math.Max(0, Math.Min(buffer.Length, n - pos));
            return count == 0 ? 0 : RandomAccess.Read(input, buffer.AsSpan(0, count), pos);
        }

        var viaSearch = new QueryAnswer[queries.Count];
        foreach (HandoffQuery q in queries)
        {
            (long lo, long hi) = SaSearch.Range(q.Bytes, n, Sa, ReadText);
            long count = hi - lo;
            List<long>? positions = null;
            if (count <= options.Cap)
            {
                positions = [];
                for (long r = lo; r < hi; r++) positions.Add(Sa(r));
                positions.Sort();
            }
            viaSearch[q.Index] = new QueryAnswer(count, count > 0, positions);
        }
        List<HandoffMismatch> saMismatches = HandoffCompare.Compare("sa-search", queries, oracle, viaSearch, options.Cap);
        row("q_sa_s", Sec(t), "binary searches over the suffix array of the container and the input text");
        row("q_sa_mismatches", Num(saMismatches.Count), "count, exists and locate answers of the search that differ from the oracle");
        all.AddRange(saMismatches);

        // the engine
        string status = "skipped";
        if (options.Engine != "false")
        {
            t = Stopwatch.GetTimestamp();
            ISuffixTree? tree = null;
            try
            {
                tree = options.OpenEngine is { } open ? open(options.Container) : SuffixIndexFactory.Open(options.Container, new EsaOpenOptions());
                status = "ok";
            }
            catch (NotSupportedException e)
            {
                status = "not_applicable";
                row("q_engine_reason", e.Message.Replace('\t', ' ').Replace('\n', ' '), "the engine refuses a container beyond its 32-bit surface");
                if (options.Engine == "true") all.Add(new HandoffMismatch(-1, "engine", "open", "the engine cannot open the container: " + e.Message));
            }
            if (tree is not null)
            {
                using (tree as IDisposable)
                {
                    var viaEngine = new QueryAnswer[queries.Count];
                    string? refusal = null;
                    int refused = 0;
                    foreach (HandoffQuery q in queries)
                    {
                        try
                        {
                            string s = System.Text.Encoding.Latin1.GetString(q.Bytes);
                            long count = tree.CountOccurrences(s);
                            bool exists = tree.Contains(s);
                            List<long>? positions = null;
                            if (count <= options.Cap)
                            {
                                positions = [.. tree.FindAllOccurrences(s).Select(p => (long)p)];
                                positions.Sort();
                            }
                            viaEngine[q.Index] = new QueryAnswer(count, exists, positions);
                        }
                        catch (InvalidDataException e)
                        {
                            // the Checked mode refuses to read a damaged file: that is an answer that differs from the oracle's (the sentinel count is -1)
                            refused++;
                            refusal ??= e.Message.Replace('\t', ' ').Replace('\n', ' ');
                            viaEngine[q.Index] = new QueryAnswer(-1, false, null);
                        }
                    }
                    row("q_engine_refused", Num(refused), refusal ?? "queries that the Checked mode refused as a damaged file");
                    List<HandoffMismatch> engineMismatches = HandoffCompare.Compare("engine", queries, oracle, viaEngine, options.Cap);
                    row("q_engine_s", Sec(t), "Contains, CountOccurrences and FindAllOccurrences of every query through the file-backed tree (Header, Checked)");
                    row("q_engine_mismatches", Num(engineMismatches.Count), "answers of the engine that differ from the oracle");
                    all.AddRange(engineMismatches);
                }
            }
        }
        row("q_engine_status", status, "ok: the engine answered; not_applicable: the container is above the 32-bit surface (a text of 2^31 - 1 symbols or more); skipped: switched off");
        if (status != "ok")
        {
            row("q_engine_mismatches", "0", "the engine did not answer");
        }

        foreach (HandoffMismatch m in all.Take(10)) row("q_mismatch", $"{m.Answerer} {m.Kind} query {m.Index}", m.Detail);
        bool pass = all.Count == 0;
        row("q_result", pass ? "PASS" : "FAIL", pass ? "every answer equals the oracle's" : $"{all.Count} answers differ from the oracle");
        return pass ? 0 : 3;
    }

    private static byte[] ReadBytes(SafeFileHandle handle, long position, int length)
    {
        var buffer = new byte[length];
        int done = 0;
        while (done < length)
        {
            int got = RandomAccess.Read(handle, buffer.AsSpan(done), position + done);
            if (got <= 0) throw new EndOfStreamException($"the input ended at {position + done}");
            done += got;
        }
        return buffer;
    }

    private static string SetHash(IReadOnlyList<HandoffQuery> queries)
    {
        using var sha = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        Span<byte> head = stackalloc byte[13];
        foreach (HandoffQuery q in queries)
        {
            BinaryPrimitives.WriteInt32LittleEndian(head, q.Bytes.Length);
            BinaryPrimitives.WriteInt64LittleEndian(head[4..], q.Origin);
            head[12] = q.Present ? (byte)1 : (byte)0;
            sha.AppendData(head);
            sha.AppendData(q.Bytes);
        }
        return Convert.ToHexStringLower(sha.GetHashAndReset());
    }
}
