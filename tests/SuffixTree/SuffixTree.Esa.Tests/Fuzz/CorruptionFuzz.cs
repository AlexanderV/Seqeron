using System.Buffers.Binary;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>The kinds of damage of the container corruption fuzz (the kinds of the corruption, hostile and isolation tests of 1b.7f, drawn at random).</summary>
internal enum DamageKind
{
    /// <summary>One byte of the header page flipped.</summary>
    HeaderFlip,

    /// <summary>One byte of the body flipped, no checksum recomputed.</summary>
    BodyFlip,

    /// <summary>One byte of the suffix array, the LCP array, the child table or the PSI section flipped and every checksum recomputed.</summary>
    SealedByteFlip,

    /// <summary>The file cut short.</summary>
    Truncation,

    /// <summary>A share of the words of one section replaced by random or extreme values, sealed.</summary>
    SealedWords,

    /// <summary>Two adjacent suffix-array entries exchanged, sealed.</summary>
    SealedSwap,

    /// <summary>Bytes appended behind the last section.</summary>
    TrailingBytes,

    /// <summary>One byte of a section-table entry flipped, the header checksum recomputed.</summary>
    TableField,

    /// <summary>One byte of the header fields (the magic, the sizes, the statistics, the hash of the text) flipped and the header checksum recomputed: what the checksum no longer guards.</summary>
    HeaderField,

    /// <summary>One byte of the packed text, the alphabet, the symbol counts or the build information flipped and every checksum recomputed.</summary>
    SealedPayloadFlip,

    /// <summary>One entry of the suffix array, the LCP array or the child table replaced by the entry next to it, by its own index, by the index next to it or by itself plus or minus one, sealed: the values that are plausible and wrong (equal depths, a cell that points at itself, an off-by-one).</summary>
    SealedNeighbour,
}

/// <summary>Counts of the corruption lane, per kind of damage and outcome, thread safe.</summary>
internal sealed class CorruptionStats
{
    private static readonly string[] Outcomes = ["refused@Header", "refused@Sections", "refused@Structure", "opened-equal", "skipped"];
    private readonly long[,] _outcomes = new long[Enum.GetValues<DamageKind>().Length, Outcomes.Length];
    private long _batteryOk, _batteryInvalid, _batteryRuns, _guarded;

    internal void Count(DamageKind kind, string outcome) => Interlocked.Increment(ref _outcomes[(int)kind, Array.IndexOf(Outcomes, outcome)]);

    internal void Battery(int ok, int invalid, bool guarded)
    {
        Interlocked.Increment(ref _batteryRuns);
        Interlocked.Add(ref _batteryOk, ok);
        Interlocked.Add(ref _batteryInvalid, invalid);
        if (guarded) Interlocked.Increment(ref _guarded);
    }

    internal long Total(DamageKind kind)
    {
        long sum = 0;
        for (int o = 0; o < Outcomes.Length - 1; o++) sum += Interlocked.Read(ref _outcomes[(int)kind, o]);
        return sum;
    }

    internal long Cases => Enum.GetValues<DamageKind>().Sum(Total);

    internal long Skipped => Enum.GetValues<DamageKind>().Sum(k => Interlocked.Read(ref _outcomes[(int)k, Outcomes.Length - 1]));

    internal long BatteryRuns => Interlocked.Read(ref _batteryRuns);

    internal long BatteryOk => Interlocked.Read(ref _batteryOk);

    internal long BatteryInvalid => Interlocked.Read(ref _batteryInvalid);

    internal long GuardedRuns => Interlocked.Read(ref _guarded);

    public override string ToString() => string.Join("; ", Enum.GetValues<DamageKind>().Select(k =>
        $"{k}: {Total(k)} (" + string.Join(", ", Outcomes.Take(Outcomes.Length - 1).Select((o, i) => $"{o} {Interlocked.Read(ref _outcomes[(int)k, i])}")) + ")"));
}

/// <summary>What a corruption case needs besides its text.</summary>
internal sealed class CorruptionContext
{
    internal required MemberCounters Counters { get; init; }

    internal required string WorkDir { get; init; }

    internal CorruptionStats Stats { get; init; } = new();
}

/// <summary>
/// The container corruption fuzz of the harness (work package 1b.8): a text from the generator families, its container, and one damage drawn from the kinds of 1b.7f. The outcomes that are allowed: every open at
/// every level either succeeds or is refused with <see cref="InvalidDataException"/> or <see cref="NotSupportedException"/> (never another exception), a refusal at a level is a refusal at every higher level, a file
/// that every level accepts answers like the undamaged text (the digest equals that of the in-memory tree), and the kinds whose level was fixed by the tests of 1b.7f are refused at that level (a header flip, a cut,
/// bytes behind the last section and a table flip at <c>Header</c>, a body flip at <c>Sections</c> at the latest, a sealed flip or swap at <c>Structure</c> at the latest). A file that the <c>Header</c> level accepts is then
/// read in the <c>Checked</c> mode from memory between two inaccessible pages with every kind of operation, each of which must return or throw <see cref="InvalidDataException"/> (a locus that a damaged table produced may fail the
/// argument check of the next call): a read outside the mapping ends the process, a loop without end is stopped by the watchdog of the runner.
/// </summary>
internal static class CorruptionFuzz
{
    private static readonly EsaValidationLevel[] Levels = [EsaValidationLevel.Header, EsaValidationLevel.Sections, EsaValidationLevel.Structure];

    /// <summary>Case <paramref name="index"/> of a run: the text (as a generator case), the kind of damage and the seed of the damage.</summary>
    internal static (FuzzCase Text, DamageKind Kind, ulong DamageSeed) Draw(ulong masterSeed, long index)
    {
        ulong h = masterSeed ^ unchecked(0x1B8C0AAFEUL * 0xD1B54A32D192ED03UL) ^ unchecked((ulong)index * 0xD6E8FEB86659FD93UL);
        var rng = new SplitMix64(new SplitMix64(h).Next());
        int band = rng.NextInt(100);
        int n = band < 70 ? 1 + rng.NextInt(64) : band < 95 ? 65 + rng.NextInt(536) : 601 + rng.NextInt(3400);
        FuzzCase c = FuzzTexts.Pick(rng, n);
        var kind = (DamageKind)rng.NextInt(Enum.GetValues<DamageKind>().Length);
        return (c, kind, rng.Next());
    }

    /// <summary>The operation name of a corruption failure, as written in the regression corpus: the kind and the seed of the damage.</summary>
    internal static string Op(DamageKind kind, ulong damageSeed) => $"corruption:{kind}:{damageSeed:x16}";

    /// <summary>Parses <see cref="Op"/>; null for any other operation.</summary>
    internal static (DamageKind Kind, ulong DamageSeed)? ParseOp(string op)
    {
        string[] f = op.Split(':');
        if (f.Length != 3 || f[0] != "corruption" || !Enum.TryParse(f[1], out DamageKind kind) || !ulong.TryParse(f[2], System.Globalization.NumberStyles.AllowHexSpecifier, System.Globalization.CultureInfo.InvariantCulture, out ulong seed)) return null;
        return (kind, seed);
    }

    /// <summary>Runs one corruption case; null when every outcome was allowed, otherwise the description of the first that was not.</summary>
    internal static string? Run(string label, string text, DamageKind kind, ulong damageSeed, CorruptionContext ctx)
    {
        if (text.Length == 0 || text.Any(ch => ch > 255))
        {
            ctx.Stats.Count(kind, "skipped");
            return null;
        }
        var rng = new SplitMix64(damageSeed);
        bool psi = rng.NextInt(4) != 0;
        byte[] clean = SesaTestFile.Image(LcpChildCase.Of(text), psi);
        (byte[] image, string what, EsaValidationLevel? must) = Damage(clean, kind, rng);
        bool noop = image.AsSpan().SequenceEqual(clean);
        string where = $"{label} n={text.Length} {kind} (damage seed {damageSeed:x16}, psi {psi}): {what}";
        string path = Path.Combine(ctx.WorkDir, "d.sesa");
        File.WriteAllBytes(path, image);

        EsaValidationLevel? refusedAt = null;
        bool dispatcher = (damageSeed & 2) != 0; // half of the cases are opened through SuffixIndexFactory.Open, the others through EsaSuffixTreeFactory.Open
        foreach (EsaValidationLevel level in Levels)
        {
            try
            {
                using EsaSuffixTree tree = OpenAt(path, level, EsaAccessMode.Checked, level == EsaValidationLevel.Structure, dispatcher);
                if (refusedAt is { } earlier) return $"{where}: refused at {earlier} but opened at {level}";
            }
            catch (Exception e) when (e is InvalidDataException or NotSupportedException)
            {
                refusedAt ??= level;
            }
            catch (Exception e)
            {
                return $"{where}: the open at {level} threw {e.GetType().Name}: {e.Message}";
            }
        }

        if (!noop && must is { } mustLevel && (refusedAt is null || refusedAt > mustLevel))
            return $"{where}: this damage must be refused at {mustLevel} at the latest, the file was {(refusedAt is null ? "accepted" : "refused only at " + refusedAt)}";

        if (refusedAt is null)
        {
            string? different = CompareWithTheTrueText(path, text, damageSeed, dispatcher);
            if (different is not null) return $"{where}: every level accepted the file, but {different}";
            ctx.Stats.Count(kind, "opened-equal");
            return null;
        }
        ctx.Stats.Count(kind, refusedAt switch { EsaValidationLevel.Header => "refused@Header", EsaValidationLevel.Sections => "refused@Sections", _ => "refused@Structure" });
        if (refusedAt == EsaValidationLevel.Header) return null;
        return Battery(image, text, damageSeed, where, ctx);
    }

    private static EsaSuffixTree OpenAt(string path, EsaValidationLevel level, EsaAccessMode access, bool sha, bool dispatcher) => dispatcher
        ? (EsaSuffixTree)SuffixIndexFactory.Open(path, new EsaOpenOptions { Validation = level, Access = access, VerifyTextSha256 = sha })
        : SesaTestFile.Open(path, level, access, sha);

    private static string? CompareWithTheTrueText(string path, string text, ulong seed, bool dispatcher)
    {
        FuzzScript script = FuzzScript.Create(text, seed);
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        List<string> expected = FuzzDigest.Compute(memory, text, script, DigestLimits.Default);
        foreach (EsaAccessMode access in new[] { EsaAccessMode.Checked, EsaAccessMode.Trusted })
        {
            using EsaSuffixTree tree = OpenAt(path, EsaValidationLevel.Structure, access, sha: true, dispatcher);
            if (FuzzDigest.Diff(FuzzDigest.Compute(tree, text, script, DigestLimits.Default), expected, $"the file ({access})", "the in-memory tree") is { } diff) return diff;
        }
        return null;
    }

    private static (byte[] Image, string What, EsaValidationLevel? Must) Damage(byte[] clean, DamageKind kind, SplitMix64 rng)
    {
        byte[] image = (byte[])clean.Clone();
        SesaOracle o = SesaOracle.FromBytes(image);
        switch (kind)
        {
            case DamageKind.HeaderFlip:
            {
                int at = rng.NextInt(SesaFormat.PageBytes);
                image[at] ^= (byte)(1 + rng.NextInt(255));
                return (image, $"flip of page byte {at}", EsaValidationLevel.Header);
            }
            case DamageKind.BodyFlip:
            {
                int at = SesaFormat.PageBytes + rng.NextInt(image.Length - SesaFormat.PageBytes);
                image[at] ^= (byte)(1 + rng.NextInt(255));
                return (image, $"flip of body byte {at}", EsaValidationLevel.Sections);
            }
            case DamageKind.SealedByteFlip:
            {
                OracleSection s = Pick(o, rng, SesaFormat.KindSa, SesaFormat.KindLcp, SesaFormat.KindChild, SesaFormat.KindPsi);
                int at = (int)s.Offset + rng.NextInt((int)s.Length);
                image[at] ^= (byte)(1 + rng.NextInt(255));
                SesaTestFile.Reseal(image);
                return (image, $"sealed flip of byte {at} of section {s.Kind}", EsaValidationLevel.Structure);
            }
            case DamageKind.Truncation:
            {
                int length = rng.NextInt(image.Length);
                return (image[..length], $"truncation to {length}", EsaValidationLevel.Header);
            }
            case DamageKind.SealedWords:
            {
                OracleSection s = Pick(o, rng, SesaFormat.KindSa, SesaFormat.KindChild, SesaFormat.KindLcp, SesaFormat.KindPsi, SesaFormat.KindLcpOverflow);
                int words = (int)(s.Length / 4);
                if (words == 0) return (image, "no words", null);
                long n = o.TextLen;
                uint[] extreme = [0, 1, (uint)n, (uint)n + 1, (uint)n + 2, 0x7FFFFFFF, 0x80000000, 0xFFFFFFFF, 0xFFFFFFFE, 0x40000000, 0x3FFFFFFF, 0xC0000001];
                int hits = 1 + rng.NextInt(Math.Max(1, words / 5));
                for (int h = 0; h < hits; h++)
                {
                    int w = rng.NextInt(words);
                    uint v = rng.NextInt(3) switch { 0 => extreme[rng.NextInt(extreme.Length)], 1 => (uint)rng.Next(), _ => (uint)rng.NextInt((int)Math.Max(2, n + 2)) };
                    BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan((int)s.Offset + (4 * w)), v);
                }
                SesaTestFile.Reseal(image);
                return (image, $"{hits} word(s) of section {s.Kind} replaced, sealed", null);
            }
            case DamageKind.SealedSwap:
            {
                OracleSection s = o.Get(SesaFormat.KindSa);
                long n = o.TextLen;
                if (n < 2) return (image, "text too short for a swap", null);
                int r = 1 + rng.NextInt((int)n - 1);
                int at = (int)s.Offset + (4 * r);
                int a = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan(at));
                int b = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan(at + 4));
                BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(at), b);
                BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(at + 4), a);
                SesaTestFile.Reseal(image);
                return (image, $"entries {r} and {r + 1} of the suffix array exchanged, sealed", EsaValidationLevel.Structure);
            }
            case DamageKind.TrailingBytes:
            {
                int extra = 1 + rng.NextInt(64);
                var longer = new byte[image.Length + extra];
                image.CopyTo(longer, 0);
                for (int i = image.Length; i < longer.Length; i++) longer[i] = (byte)rng.NextInt(256);
                return (longer, $"{extra} bytes appended", EsaValidationLevel.Header);
            }
            case DamageKind.HeaderField:
            {
                int at = rng.NextInt(124);
                image[at] ^= (byte)(1 + rng.NextInt(255));
                int count = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan(48));
                if (count is >= 1 and <= SesaFormat.MaxSections) SesaTestFile.ResealHeader(image); // a damaged SECTION_COUNT has no checksum to compute: the stale one refuses the file
                // the magic, the endian probe, the file size (it must be the length of the file) and the reserved word are fixed by the format (ADR-003): a change is refused at Header
                return (image, $"flip of header byte {at}, header checksum recomputed", at < 8 || at is >= 12 and < 16 || at is >= 32 and < 40 || at is >= 120 and < 124 ? EsaValidationLevel.Header : null);
            }
            case DamageKind.SealedPayloadFlip:
            {
                OracleSection s = Pick(o, rng, SesaFormat.KindText, SesaFormat.KindAlphabet, SesaFormat.KindCounts, SesaFormat.KindBuildInfo);
                int at = (int)s.Offset + rng.NextInt((int)s.Length);
                image[at] ^= (byte)(1 + rng.NextInt(255));
                SesaTestFile.Reseal(image);
                // the symbol counts must sum to the length of the text: any single changed byte breaks the sum
                return (image, $"sealed flip of byte {at} of section {s.Kind}", s.Kind == SesaFormat.KindCounts ? EsaValidationLevel.Header : null);
            }
            case DamageKind.SealedNeighbour:
            {
                OracleSection s = Pick(o, rng, SesaFormat.KindSa, SesaFormat.KindLcp, SesaFormat.KindChild);
                int width = s.Kind == SesaFormat.KindLcp ? Math.Max(1, s.ElemWidth) : 4;
                int entries = (int)(s.Length / width);
                if (entries < 3) return (image, "section too short", null);
                int r = 1 + rng.NextInt(entries - 2);
                int mode = rng.NextInt(7);
                if (mode == 6 && s.Kind == SesaFormat.KindLcp && width <= 2 && PeakOf(image, s, width, entries, rng) is { } peak) r = peak; // an interior peak of the LCP array: one less keeps the shape of the interval tree
                int at = (int)s.Offset + (r * width);
                if (mode >= 4)
                {
                    // the entry itself, one more or one less (an off-by-one: the damage that keeps the structure of the tree and changes a depth or a position)
                    long current = 0;
                    for (int b = 0; b < width; b++) current |= (long)image[at + b] << (8 * b);
                    long changed = mode == 4 ? current + 1 : current - 1; // 4: plus one; 5 and 6: minus one
                    for (int b = 0; b < width; b++) image[at + b] = (byte)(changed >> (8 * b));
                }
                else if (mode < 2)
                {
                    int from = (int)s.Offset + ((r + (mode == 0 ? -1 : 1)) * width);
                    image.AsSpan(from, width).CopyTo(image.AsSpan(at, width));
                }
                else
                {
                    long value = mode == 2 ? r : r + 1;
                    for (int b = 0; b < width; b++) image[at + b] = (byte)(value >> (8 * b));
                }
                SesaTestFile.Reseal(image);
                return (image, $"entry {r} of section {s.Kind} replaced by {(mode switch { 0 => "the previous entry", 1 => "the next entry", 2 => "its own index", 3 => "its index + 1", 4 => "itself + 1", _ => "itself - 1" })}, sealed", null);
            }
            default:
            {
                int index = rng.NextInt(o.SectionCount);
                int b = rng.NextInt(64);
                int at = 128 + (64 * index) + b;
                image[at] ^= (byte)(1 + rng.NextInt(255));
                SesaTestFile.ResealHeader(image);
                return (image, $"flip of byte {b} of section-table entry {index}, header checksum recomputed", null);
            }
        }
    }

    /// <summary>A rank whose LCP value stands at least two above both neighbours and below the maximum: decreasing it by one changes no comparison with another value, so the child table, the node count and the longest repeat stay what they are and only the comparison with the text can see it.</summary>
    private static int? PeakOf(byte[] image, OracleSection s, int width, int entries, SplitMix64 rng)
    {
        var value = new int[entries];
        for (int i = 0; i < entries; i++)
        {
            int v = 0;
            for (int b = 0; b < width; b++) v |= image[(int)s.Offset + (i * width) + b] << (8 * b);
            value[i] = v;
        }
        int max = value.Max();
        var peaks = new List<int>();
        for (int i = 2; i < entries - 1; i++)
            if (value[i] >= Math.Max(value[i - 1], value[i + 1]) + 2 && value[i] < max && value[i] < (1 << (8 * width)) - 1) peaks.Add(i);
        return peaks.Count == 0 ? null : peaks[rng.NextInt(peaks.Count)];
    }

    private static OracleSection Pick(SesaOracle o, SplitMix64 rng, params int[] kinds)
    {
        OracleSection[] present = [.. kinds.Select(o.Find).Where(s => s is not null && s.Length > 0).Select(s => s!)];
        return present[rng.NextInt(present.Length)];
    }

    // ------------------------------------------------------------------ the operations on a damaged file that passed the Header level

    private sealed class Counting : ISuffixTreeVisitor
    {
        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
        }

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }

    private static string? Battery(byte[] image, string text, ulong damageSeed, string where, CorruptionContext ctx)
    {
        using var guarded = new GuardedImage(image, flush: (damageSeed & 1) == 0);
        int ok = 0, invalid = 0;
        string? failure = null;
        void Step(string what, bool navigation, Action action)
        {
            if (failure is not null) return;
            try
            {
                action();
                ok++;
            }
            catch (InvalidDataException)
            {
                invalid++;
            }
            catch (ArgumentException e) when (navigation && e.GetType() == typeof(ArgumentException))
            {
                invalid++; // a locus that a damaged table produced fails the checks of the next call
            }
            catch (Exception e)
            {
                failure = $"{where}: {what} threw {e.GetType().Name}: {e.Message} at {string.Join(" <- ", (e.StackTrace ?? string.Empty).Split('\n').Take(6).Select(l => l.Trim()))}";
            }
        }

        FuzzScript script = FuzzScript.Create(text, damageSeed);
        FuzzSteps.Run(FuzzSteps.Corruption, () =>
        {
            unsafe
            {
                EsaSuffixTree tree;
                try
                {
                    SesaFile file = SesaFile.FromMemory(guarded.Data, guarded.Length, new EsaOpenOptions { Validation = EsaValidationLevel.Header, Access = EsaAccessMode.Checked });
                    tree = EsaSuffixTree.FromFile(file);
                }
                catch (Exception e) when (e is InvalidDataException or NotSupportedException)
                {
                    return; // refused from memory at the Header level: the same refusal as through the file
                }
                using (tree)
                {
                    var counted = new CountingTree(tree, ctx.Counters, ctx.Counters.EngineIndex("file-header-checked"));
                    RunSteps(counted, text, script, Step);
                }
            }
        });
        ctx.Stats.Battery(ok, invalid, guarded.Guarded);
        return failure;
    }

    private static void RunSteps(CountingTree tree, string text, FuzzScript script, Action<string, bool, Action> step)
    {
        foreach (string p in script.Patterns.Take(10))
        {
            step("Contains", false, () => tree.Contains(p));
            step("CountOccurrences", false, () => tree.CountOccurrences(p.AsSpan()));
            step("FindAllOccurrences", false, () => tree.FindAllOccurrences(p));
        }
        step("LongestRepeatedSubstring", false, () => tree.LongestRepeatedSubstring());
        step("LongestRepeatedSubstringMemory", false, () => tree.LongestRepeatedSubstringMemory());
        step("FindAllLongestRepeatedSubstrings", false, () => tree.FindAllLongestRepeatedSubstrings());
        step("CountDistinctSubstrings", false, () => tree.CountDistinctSubstrings());
        step("CountDistinctSubstringsByLength", false, () => tree.CountDistinctSubstringsByLength(9));
        if (text.Length <= 600) step("GetAllSuffixes", false, () => tree.GetAllSuffixes());
        step("EnumerateSuffixes", false, () => tree.EnumerateSuffixes().Take(5).ToList());
        step("Traverse", false, () => tree.Traverse(new Counting()));
        if (text.Length <= 300) step("PrintTree", false, () => tree.PrintTree());
        foreach (string query in script.Queries.Take(4))
        {
            step("LongestCommonSubstring", false, () => tree.LongestCommonSubstring(query));
            step("LongestCommonSubstringInfo", false, () => tree.LongestCommonSubstringInfo(query));
            step("FindAllLongestCommonSubstrings", false, () => tree.FindAllLongestCommonSubstrings(query));
            step("FindAllDistinctLongestCommonSubstrings", false, () => tree.FindAllDistinctLongestCommonSubstrings(query));
            step("FindMaximalExactMatches", false, () => tree.FindMaximalExactMatches(query, 3));
            step("FindMaximalUniqueMatches", false, () => tree.FindMaximalUniqueMatches(query, 3));
            step("FindExactMatchAnchors", false, () => tree.FindExactMatchAnchors(query, 3));
        }
        if (text.Length <= 400) step("FindMaximalRepeatedPairs", false, () => tree.FindMaximalRepeatedPairs(3));
        step("Text", false, () =>
        {
            ITextSource source = tree.Text;
            _ = source.ToString();
            if (text.Length > 4) _ = source.Slice(1, Math.Min(3, text.Length - 1));
        });
        step("IIndexText", false, () =>
        {
            IIndexText indexText = ((ISuffixIndex)tree).Text;
            _ = indexText.SymbolAt(text.Length / 2);
            _ = indexText.CopyTo(0, new char[8]);
        });
        long rank = Math.Min(3, tree.Text.Length);
        step("SuffixAt", true, () => tree.SuffixAt(rank));
        step("LcpAt", true, () => tree.LcpAt(rank));
        step("RankOf", true, () => tree.RankOf(rank));
        step("Lcp", true, () => tree.Lcp(1, Math.Max(1, tree.Text.Length)));
        step("Lca", true, () => tree.Lca(1, Math.Max(1, tree.Text.Length)));
        step("walk", true, () => Walk(tree, tree.Root, 0, 6));
        step("TryExtendLeft", true, () => tree.TryExtendLeft(tree.Root, 'a', out _));
        step("Support", true, () =>
        {
            foreach (IndexOperation op in Enum.GetValues<IndexOperation>())
            {
                _ = tree.Support(op);
                _ = tree.SupportNote(op);
            }
            _ = tree.Engine;
            _ = tree.TextLength;
        });
    }

    private static void Walk(ISuffixTreeNavigation nav, IndexLocus node, int depth, int maxDepth)
    {
        if (depth >= maxDepth) return;
        var children = new IndexLocus[64];
        int count = nav.GetChildren(node, children);
        _ = nav.Count(node);
        _ = nav.IsExplicitNode(node);
        var positions = new long[4];
        _ = nav.Locate(node, 0, positions);
        for (int i = 0; i < Math.Min(count, children.Length); i++)
        {
            _ = nav.TryParent(children[i], out _);
            _ = nav.TrySuffixLink(children[i], out _);
            _ = nav.EdgeSymbol(children[i], 0);
            _ = nav.TryChild(node, nav.EdgeSymbol(children[i], 0), out _);
            _ = nav.TryExtendRight(node, 'a', out _);
            Walk(nav, children[i], depth + 1, maxDepth);
        }
    }
}
