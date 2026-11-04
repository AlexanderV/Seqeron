using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>
/// 1o.3b: the suffix-array-order LCP kernel (<see cref="LcpKernel.BuildDirect{TSym}"/>) beyond the generic registry run: the vector-width boundaries, each reason
/// to give up (and the fallback result), the success path with escapes, the argument checks, and memory safety on corrupt input.
/// </summary>
[TestFixture]
public class LcpDirectTests
{
    private static LcpEngine Engine => KernelRegistry.DirectLcp;

    private static CompactLcp Reference(PreparedText p, int[] sa)
        => p.Width == SymbolWidth.Byte ? LcpReference.Build(p.ByteReader(), sa) : LcpReference.Build(p.UShortReader(), sa);

    /// <summary>Runs the direct attempt alone on a text with its true suffix array; null when it gave up.</summary>
    private static CompactLcp? Attempt(string text, out int[] sa, out PreparedText prepared, DirectReport? report = null)
    {
        prepared = TextPreparation.Prepare(text);
        sa = SaisReference.BuildSuffixArray(text);
        return Attempt(prepared, sa, report);
    }

    /// <summary>Runs the direct attempt alone on the given (possibly corrupt) suffix array in a guarded arena and verifies the canaries.</summary>
    private static CompactLcp? Attempt(PreparedText prepared, int[] sa, DirectReport? report = null)
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        int n = prepared.Length;
        Buf<int> saBuf = arena.Allocate<int>(n + 1, zero: false, "sa");
        sa.CopyTo(saBuf.AsSpan(n + 1));
        CompactLcp? result;
        int blocks;
        if (prepared.Width == SymbolWidth.Byte)
        {
            Buf<byte> t = arena.Allocate<byte>(n, zero: false, "text");
            prepared.Bytes.CopyTo(t.AsSpan(n));
            blocks = arena.BlockCount;
            result = LcpKernel.TryDirect(t, n, saBuf, arena, report: report);
        }
        else
        {
            Buf<ushort> t = arena.Allocate<ushort>(n, zero: false, "text");
            prepared.UShorts.CopyTo(t.AsSpan(n));
            blocks = arena.BlockCount;
            result = LcpKernel.TryDirect(t, n, saBuf, arena, report: report);
        }
        arena.VerifyCanaries();
        Assert.That(arena.BlockCount, Is.EqualTo(blocks), "the attempt, successful or not, returns all its workspace to the arena");
        return result;
    }

    private static void AssertAttempt(string text, bool succeeds)
    {
        CompactLcp? direct = Attempt(text, out int[] sa, out PreparedText p);
        Assert.That(direct is not null, Is.EqualTo(succeeds), $"direct attempt on {text.Length} symbols");
        if (direct is not null)
        {
            CompactLcp reference = Reference(p, sa);
            Assert.That(direct.Kind, Is.EqualTo(reference.Kind));
            Assert.That(direct.OverflowCount, Is.EqualTo(reference.OverflowCount));
            Assert.That(direct.ToArray(), Is.EqualTo(reference.ToArray()));
        }
        // the resolution in text order alone (forced) gives the same result whatever the budget would have decided, and gives up exactly where the scan gives up
        var forced = new DirectReport { ForceTextOrder = true };
        CompactLcp? viaTextOrder = Attempt(p, sa, forced);
        Assert.That(viaTextOrder is not null, Is.EqualTo(succeeds), $"forced text order on {text.Length} symbols");
        if (viaTextOrder is not null)
        {
            CompactLcp reference = Reference(p, sa);
            Assert.That(forced.Outcome, Is.EqualTo(DirectOutcome.ResolvedInTextOrder));
            Assert.That(viaTextOrder.Kind, Is.EqualTo(reference.Kind));
            Assert.That(viaTextOrder.OverflowCount, Is.EqualTo(reference.OverflowCount));
            Assert.That(viaTextOrder.ToArray(), Is.EqualTo(reference.ToArray()));
        }
        KernelHarness.AssertLcp(Engine, text, allPoisons: false); // with or without the attempt, the engine gives the reference result
        KernelHarness.AssertLcp(KernelRegistry.DirectTextOrderLcp, text, allPoisons: false);
    }

    [Test]
    public void TinyTexts_n0_To_n3_AndTheTwoEntryPoints()
    {
        foreach (string text in new[] { "", "a", "aa", "ab", "ba", "aaa", "aba", "abc", "cab" })
        {
            KernelHarness.AssertLcp(Engine, text);
            PreparedText p = TextPreparation.Prepare(text);
            int[] sa = SaisReference.BuildSuffixArray(text);
            Assert.That(LcpKernel.Build(LcpKernelKind.Direct, p, sa).ToArray(), Is.EqualTo(NaiveLcp.Build(text, sa)), $"selector on '{text}'");
        }
    }

    [Test]
    public void VectorWidthBoundaries_PeriodicTexts_AtEveryLengthAndPeriod()
    {
        // 16 symbols per vector (bytes) / 8 (16-bit): suffix pairs with 15/16/17 symbols left, matches that end in the first vector, exactly at its end, and in the next one
        for (int period = 1; period <= 20; period++)
            for (int n = 0; n <= 80; n++)
            {
                string s = StringFamilies.PeriodicAaB(n, period + 1);
                KernelHarness.AssertLcp(Engine, s, allPoisons: false);
                KernelHarness.AssertLcp(Engine, new string([.. s.Select(c => (char)(c + 0x200))]), allPoisons: false);
            }
        string distinct = new([.. Enumerable.Range(0, 300).Select(i => (char)(0x100 + i))]);
        Assert.That(TextPreparation.Prepare(distinct + "ab").Width, Is.EqualTo(SymbolWidth.UShort));
        for (int period = 1; period <= 12; period++)
            for (int n = 0; n <= 50; n++)
                KernelHarness.AssertLcp(Engine, distinct + StringFamilies.PeriodicAaB(n, period + 1), allPoisons: false);
    }

    [Test]
    public void RandomDna_SucceedsDirectly_AndLongRepeatsAreResolved()
    {
        AssertAttempt(StringFamilies.RandomDna(100_000, 5), succeeds: true);
        // a block of 1,000 twice inside random filler: the pair gives one LCP entry per offset, 746 escapes (16 * 746 <= 40,000: U8), 984 long ranks (<= n / 4)
        string block = StringFamilies.Random(1000, 4, 7);
        string text = StringFamilies.Random(19_000, 4, 8) + block + StringFamilies.Random(19_000, 4, 9) + block;
        AssertAttempt(text, succeeds: true);
        Assert.That(Attempt(text, out _, out _)!.OverflowCount, Is.InRange(700, 800));
    }

    [Test]
    public void TooManyLongRanks_GiveUp_AndTheFallbackIsTheReference()
    {
        AssertAttempt(StringFamilies.AllSame(1000), succeeds: false);
        AssertAttempt(StringFamilies.AbRepeat(5000), succeeds: false);
        AssertAttempt(StringFamilies.Fibonacci(5000), succeeds: false);
    }

    [Test]
    public void TheLongRankCap_IsExactlyAQuarterOfTheRanks()
    {
        // a^m followed by a random tail: the ranks with an LCP of at least 16 are counted from the reference and compared with n / 4 on both sides of the limit
        int below = 0, above = 0;
        for (int m = 20; m <= 90; m += 7)
            for (int tail = 20; tail <= 400; tail += 13)
            {
                string text = StringFamilies.AllSame(m) + StringFamilies.Random(tail, 3, (ulong)(m * 1000 + tail)).Replace('a', 'b');
                PreparedText p = TextPreparation.Prepare(text);
                int[] sa = SaisReference.BuildSuffixArray(text);
                int[] lcp = Reference(p, sa).ToArray();
                int longs = lcp.Count(v => v >= 16);
                int n = text.Length;
                bool fitsCap = longs <= n / 4;
                CompactLcp? direct = Attempt(text, out _, out _);
                if (fitsCap) { Assert.That(direct, Is.Not.Null, $"m {m} tail {tail}: {longs} long of {n}"); below++; }
                else { Assert.That(direct, Is.Null, $"m {m} tail {tail}: {longs} long of {n}"); above++; }
            }
        Assert.That(below, Is.GreaterThan(5));
        Assert.That(above, Is.GreaterThan(5));
    }

    [Test]
    public void TheKindRule_IsTheReferenceRule_OnBothSidesOfTheU8Limit()
    {
        // Before 1o.3c the attempt gave up when the escapes did not fit the U8 rule and the three-pass kernel chose the wider kind; now the long ranks are always resolved and the kind
        // comes from the exact counts, so the direct result itself must have the reference kind.
        // 1,500 symbols copied from the start of a random text of 10,000: 1,245 escapes against a limit of n / 16 = 718, 1,500 long ranks (<= n / 4)
        string x = StringFamilies.Random(10_000, 4, 21);
        string text = x + x[..1500];
        AssertAttempt(text, succeeds: true);
        Assert.That(Attempt(text, out _, out _)!.Kind, Is.EqualTo(LcpKind.U16));
        // the limit itself passes (16 * escapes <= n), one escape more does not; a random text of f symbols followed by a copy of its first c symbols has c - 254 escapes:
        // n = 4,096 = 16 * 256 with (f, c) = (3,586, 510) is exactly on the limit, (3,585, 511) one escape over it
        foreach (var (filler, copy, u8) in new[] { (3586, 510, true), (3585, 511, false) })
        {
            string t = StringFamilies.Random(filler, 4, 77);
            t += t[..copy];
            Assert.That(t.Length, Is.EqualTo(4096));
            int escapes = Reference(TextPreparation.Prepare(t), SaisReference.BuildSuffixArray(t)).ToArray().Count(v => v >= 255);
            Assert.That(escapes, Is.EqualTo(copy - 254).Or.EqualTo(copy - 253), $"copy {copy}");
            AssertAttempt(t, succeeds: true);
            Assert.That(Attempt(t, out _, out _)!.Kind, Is.EqualTo(u8 ? LcpKind.U8 : LcpKind.U16), $"copy {copy}: escapes {escapes}");
        }
    }

    [Test]
    public void TheEscapeThresholds_At255_And65535_AreExact()
    {
        // two copies of a block of L symbols have the LCP L at one rank (the suffix of the later copy is shorter than the earlier one): that value is an escape from L = 255 on
        // (the escape value itself is an overflow entry), so there are L - 254 escapes of the byte array; one level up, 65,535 is the escape of the U16 array
        foreach (int length in new[] { 253, 254, 255, 256, 257 })
        {
            string block = StringFamilies.Random(length, 4, 6);
            string text = StringFamilies.Random(30_000, 4, 5) + block + StringFamilies.Random(5000, 4, 15) + block;
            AssertAttempt(text, succeeds: true);
            CompactLcp direct = Attempt(text, out int[] sa, out PreparedText p)!;
            CompactLcp reference = Reference(p, sa);
            Assert.That(direct.Kind, Is.EqualTo(LcpKind.U8));
            Assert.That(direct.OverflowCount, Is.EqualTo(reference.OverflowCount), $"block {length}");
            Assert.That(direct.OverflowCount, Is.EqualTo(Math.Max(0, length - 254)), $"block {length}");
        }

        foreach (int length in new[] { 65_533, 65_534, 65_535, 65_536, 65_537 })
        {
            string block = StringFamilies.Random(length, 4, 8);
            string text = StringFamilies.Random(200_000, 4, 9) + block + StringFamilies.Random(60_000, 4, 10) + block; // n about 391,000, 17 % of the ranks long
            AssertAttempt(text, succeeds: true);
            CompactLcp direct = Attempt(text, out int[] sa, out PreparedText p)!;
            CompactLcp reference = Reference(p, sa);
            Assert.That(direct.Kind, Is.EqualTo(LcpKind.U16), $"block {length}"); // 16 * escapes8 > n, 8 * escapes16 <= n
            Assert.That(direct.OverflowCount, Is.EqualTo(reference.OverflowCount), $"block {length}");
            Assert.That(direct.OverflowCount, Is.EqualTo(Math.Max(0, length - 65_534)), $"block {length}");
        }
    }

    /// <summary>
    /// Replays the per-rank resolution on the reference data without the kernel: the minimum over all its vector steps of (allowance - steps taken so far), the number of steps at the
    /// end (or at the first step over the allowance) and whether the allowance was exceeded. The allowance is <see cref="LcpKernel.ResolveBaseSteps"/> plus
    /// <see cref="LcpKernel.ResolveStepsPerLong"/> for every long rank processed so far, checked after each full-vector step.
    /// </summary>
    private static (long MinSlack, long Steps, bool Exhausted, int Processed) PerRankReplay(PreparedText p, int[] sa)
    {
        int n = p.Length, width = p.Width == SymbolWidth.Byte ? 16 : 8;
        long work = 0, budget = LcpKernel.ResolveBaseSteps, minSlack = long.MaxValue;
        int previous = n, processed = 0;
        for (int r = 1; r <= n; r++)
        {
            int s = sa[r];
            int limit = n - Math.Max(previous, s);
            bool isLong = limit >= width;
            for (int i = 0; isLong && i < width; i++) isLong = p.SymbolAt(previous + i) == p.SymbolAt(s + i);
            if (isLong)
            {
                processed++;
                budget += LcpKernel.ResolveStepsPerLong;
                int offset = width;
                while (offset + width <= limit)
                {
                    int c = 0;
                    while (c < width && p.SymbolAt(previous + offset + c) == p.SymbolAt(s + offset + c)) c++;
                    offset += c;
                    if (c < width) break;
                    work++;
                    minSlack = Math.Min(minSlack, budget - work);
                    if (work > budget) return (minSlack, work, true, processed);
                }
            }

            previous = s;
        }

        return (minSlack, work, false, processed);
    }

    [Test]
    public void ALongRepeat_ExhaustsThePerRankAllowance_AndIsResolvedInTextOrder()
    {
        // 2,000-symbol block twice in 40,000: 1,984 long ranks (<= n / 4), 1,746 escapes, about 125,000 vector steps for the per-rank resolution against an allowance of 4,096 + 4 per
        // processed rank; the attempt of 1o.3b (budget 2 n = 80,000) gave up here, now the text-order resolution takes over after a few thousand steps
        string block = StringFamilies.Random(2000, 4, 31);
        string text = StringFamilies.Random(18_000, 4, 32) + block + StringFamilies.Random(18_000, 4, 33) + block;
        Assert.That(text.Length, Is.EqualTo(40_000));
        var report = new DirectReport();
        CompactLcp? direct = Attempt(text, out int[] sa, out PreparedText p, report);
        Assert.That(direct, Is.Not.Null);
        Assert.That(report.Outcome, Is.EqualTo(DirectOutcome.ResolvedInTextOrder));
        var replay = PerRankReplay(p, sa);
        Assert.That(replay.Exhausted, Is.True);
        Assert.That(report.PerRankSteps, Is.EqualTo(replay.Steps), "the per-rank resolution stopped at the first step over its cumulative allowance");
        Assert.That(replay.Steps, Is.EqualTo(LcpKernel.ResolveBaseSteps + (LcpKernel.ResolveStepsPerLong * (long)replay.Processed) + 1));
        Assert.That(replay.Processed, Is.LessThan(report.LongRanks / 2), "the allowance is cumulative: the attempt ends long before the end of the list");
        Assert.That(report.LongPositions, Is.EqualTo(report.LongRanks));
        int[] expected = Reference(p, sa).ToArray();
        Assert.That(report.LongRanks, Is.EqualTo(expected.Count(v => v >= 16)));
        Assert.That(direct!.ToArray(), Is.EqualTo(expected));
        AssertAttempt(text, succeeds: true);
    }

    [Test]
    public void ThePerRankAllowance_IsCumulativeAndExactToTheStep()
    {
        // repeats around the length where the steps (about (m - 16)^2 / 32) cross the allowance (4,096 + 4 per processed rank, i.e. m of about 440): the outcome and the number of steps must
        // equal the replay for every text, and some texts end exactly on the allowance (kept, not exceeded) or one step over it
        int over = 0, kept = 0;
        for (int length = 380; length <= 500; length += 1)
            for (ulong seed = 0; seed < 3; seed++)
            {
                string block = StringFamilies.Random(length, 4, 1000 + seed);
                string text = StringFamilies.Random(5000, 4, 2000 + seed) + block + StringFamilies.Random(1000, 4, 3000 + seed) + block;
                var report = new DirectReport();
                Assert.That(Attempt(text, out int[] sa, out PreparedText p, report), Is.Not.Null);
                var replay = PerRankReplay(p, sa);
                Assert.That(report.Outcome, Is.EqualTo(replay.Exhausted ? DirectOutcome.ResolvedInTextOrder : DirectOutcome.ResolvedPerRank), $"length {length} seed {seed}");
                Assert.That(report.PerRankSteps, Is.EqualTo(replay.Steps), $"length {length} seed {seed}");
                if (replay.Exhausted) over++;
                else kept++;
            }

        Assert.That(kept, Is.GreaterThan(20));
        Assert.That(over, Is.GreaterThan(20));

        // a second, short repeat of q symbols moves the allowance by 4 per rank and the steps by (q - 16)^2 / 32: with it a text can be tuned until the per-rank resolution ends exactly on its
        // allowance (kept: slack 0), where a check that is one step too strict would give it up
        // (the slack of the main repeat falls by about 22 per symbol of its length, from 344 at 441 symbols: it is near zero at 457)
        bool foundExact = false;
        for (int length = 452; length <= 461 && !foundExact; length++)
            for (int q = 16; q <= 90 && !foundExact; q++)
            {
                string block = StringFamilies.Random(length, 4, 77), small = StringFamilies.Random(q, 4, 78);
                string text = StringFamilies.Random(5000, 4, 79) + block + StringFamilies.Random(1000, 4, 80) + block + StringFamilies.Random(2000, 4, 81) + small + StringFamilies.Random(1000, 4, 82) + small;
                PreparedText p = TextPreparation.Prepare(text);
                int[] sa = SaisReference.BuildSuffixArray(text);
                var replay = PerRankReplay(p, sa);
                if (replay.Exhausted || replay.MinSlack != 0) continue;
                foundExact = true;
                var report = new DirectReport();
                Assert.That(Attempt(p, sa, report), Is.Not.Null);
                Assert.That(report.Outcome, Is.EqualTo(DirectOutcome.ResolvedPerRank), $"length {length} q {q}: steps {replay.Steps} on an allowance of exactly that many");
                Assert.That(report.PerRankSteps, Is.EqualTo(replay.Steps));
            }

        Assert.That(foundExact, Is.True, "no text ended exactly on the allowance: widen the family");
    }

    [Test]
    public void ManyLongRanks_GoStraightToTextOrder_EvenWhenTheirMatchesAreShort()
    {
        // random text over three letters: about 7 % of the ranks have an LCP of 16 or more (more than TextOrderMinLongs of them at 3 million), none of them a long one, so the per-rank loop would
        // need almost no step but two random cache lines per rank; the text-order resolution is chosen up front
        const int N = 3_000_000;
        string text = StringFamilies.Random(N, 3, 91);
        PreparedText prepared = TextPreparation.Prepare(text);
        int[] sa = SaisKernel.Build(SaisKernelKind.Pointer, prepared);
        var report = new DirectReport();
        CompactLcp? direct = Attempt(prepared, sa, report);
        Assert.That(direct, Is.Not.Null);
        Assert.That(report.LongRanks, Is.GreaterThanOrEqualTo(LcpKernel.TextOrderMinLongs));
        Assert.That(report.Outcome, Is.EqualTo(DirectOutcome.ResolvedInTextOrder));
        Assert.That(report.PerRankSteps, Is.EqualTo(0), "the per-rank loop was not tried");
        CompactLcp reference = LcpKernel.Build(LcpKernelKind.Reference, prepared, sa);
        Assert.That(direct!.Kind, Is.EqualTo(reference.Kind));
        Assert.That(direct.OverflowCount, Is.EqualTo(reference.OverflowCount));
        Assert.That(direct.ToArray(), Is.EqualTo(reference.ToArray()));
    }

    [Test]
    public void ShortRepeats_StayOnThePerRankPath_WithinTheBudget()
    {
        // a repeat of 100 symbols: about 85 long ranks and (100 - 16)^2 / 32 = 220 steps, far inside 4 * 85 + 4,096; random DNA needs almost no step
        string block = StringFamilies.Random(100, 4, 41);
        string text = StringFamilies.Random(20_000, 4, 42) + block + StringFamilies.Random(20_000, 4, 43) + block;
        var report = new DirectReport();
        Assert.That(Attempt(text, out _, out _, report), Is.Not.Null);
        Assert.That(report.Outcome, Is.EqualTo(DirectOutcome.ResolvedPerRank));
        Assert.That(report.PerRankSteps, Is.InRange(100, 400));
        report = new DirectReport();
        Assert.That(Attempt(StringFamilies.RandomDna(50_000, 44), out _, out _, report), Is.Not.Null);
        Assert.That(report.Outcome, Is.EqualTo(DirectOutcome.ResolvedPerRank));
        Assert.That(report.PerRankSteps, Is.LessThanOrEqualTo(report.LongRanks), "random DNA: at most one step per long rank");
    }

    [Test]
    public void RepeatArrangements_AreResolvedExactly_WhenTheyFitTheCap()
    {
        // repeats of lengths around the vector widths in front, behind and inside the filler, at the very start and end of the text (the suffixes there have fewer symbols than a vector),
        // adjacent (tandem), mutated copies (the run of long positions breaks and restarts), three copies (the predecessor of a position changes), in byte and 16-bit symbols
        int succeeded = 0, gaveUp = 0;
        foreach (int length in new[] { 15, 16, 17, 23, 31, 32, 33, 64, 100, 255, 256, 257, 700 })
            foreach (int gap in new[] { 0, 1, 40 })
                foreach (int tail in new[] { 0, 3, 40 })
                    foreach (int head in new[] { 0, 9, 200 })
                    {
                        string block = StringFamilies.Random(length, 4, (ulong)(length * 7 + gap));
                        string text = StringFamilies.Random(head + 2000, 4, (ulong)(head + 1)) + block + StringFamilies.Random(gap, 4, (ulong)(gap + 2)) + block + StringFamilies.Random(tail, 4, (ulong)(tail + 3));
                        CheckArrangement(text, ref succeeded, ref gaveUp);
                        string substituted = block.Remove(length / 2, 1).Insert(length / 2, block[length / 2] == 'a' ? "b" : "a");
                        CheckArrangement(StringFamilies.Random(head + 2000, 4, (ulong)(head + 4)) + block + StringFamilies.Random(gap + 20, 4, (ulong)(gap + 5)) + substituted + StringFamilies.Random(tail, 4, (ulong)(tail + 6)), ref succeeded, ref gaveUp);
                        CheckArrangement(StringFamilies.Random(head + 2000, 4, (ulong)(head + 7)) + block + StringFamilies.Random(gap + 3, 4, 8) + block + StringFamilies.Random(gap + 5, 4, 9) + block + StringFamilies.Random(tail, 4, (ulong)(tail + 10)), ref succeeded, ref gaveUp);
                    }
        Assert.That(succeeded, Is.GreaterThan(1000));
        Assert.That(gaveUp, Is.LessThan(succeeded / 4));
    }

    private static void CheckArrangement(string text, ref int succeeded, ref int gaveUp)
    {
        foreach (string t in new[] { text, new string([.. text.Select(c => (char)(c + 0x300))]) })
        {
            PreparedText p = TextPreparation.Prepare(t);
            int[] sa = SaisReference.BuildSuffixArray(t);
            int width = p.Width == SymbolWidth.Byte ? 16 : 8;
            int longs = Reference(p, sa).ToArray().Count(v => v >= width);
            bool fits = longs <= t.Length / 4;
            AssertAttempt(t, fits);
            if (fits) succeeded++;
            else gaveUp++;
        }
    }

    [Test]
    public void TheCap_AndTheBudget_DecideWhereTheyShould()
    {
        var report = new DirectReport();
        // few long ranks and no long repeat: the per-rank path, the report counts the long ranks
        string random = StringFamilies.RandomDna(100_000, 12);
        Assert.That(Attempt(random, out int[] sa, out PreparedText p, report), Is.Not.Null);
        Assert.That(report.Outcome, Is.EqualTo(DirectOutcome.ResolvedPerRank));
        Assert.That(report.LongRanks, Is.EqualTo(Reference(p, sa).ToArray().Count(v => v >= 16)));

        // one block of f * n occurring twice gives about f * n long ranks: 20 % fits the cap and its repeat is long (text order), 26 % and 32 % exceed the cap of 25 %
        foreach (var (fraction, outcome) in new[] { (0.20, DirectOutcome.ResolvedInTextOrder), (0.26, DirectOutcome.GaveUp), (0.32, DirectOutcome.GaveUp) })
        {
            string text = Review.ReviewS1aOptWorstCaseTests.EmbeddedRepeat(100_000, fraction, 17);
            report = new DirectReport();
            CompactLcp? result = Attempt(text, out sa, out p, report);
            Assert.That(report.Outcome, Is.EqualTo(outcome), $"fraction {fraction}");
            Assert.That(result is not null, Is.EqualTo(outcome != DirectOutcome.GaveUp));
            KernelHarness.AssertLcp(Engine, text, allPoisons: false); // whatever the decision, the engine returns the reference result
        }

        // more than half of the ranks long: the early exit after the first n / 32 ranks (the result is the same as for the cap: gave up)
        report = new DirectReport();
        Assert.That(Attempt(StringFamilies.AllSame(100_000), out _, out _, report), Is.Null);
        Assert.That(report.Outcome, Is.EqualTo(DirectOutcome.GaveUp));
    }

    [Test]
    public void CorruptSuffixArrays_WithLongRanks_StayInsideTheBuffers()
    {
        // The resolution walks records that came from the suffix array: duplicated and displaced entries must not lead outside the text, the bitmap, the sorted records or the value array.
        // The texts have a few percent of long ranks so that the attempt reaches the resolution; the guarded arena and, in the checked build, the accessors detect any stray access.
        var rng = new SplitMix64(0x51CEUL);
        int reached = 0;
        for (int i = 0; i < 400; i++)
        {
            int n = 400 + rng.NextInt(2500);
            int block = 20 + rng.NextInt(n / 8);
            string text = StringFamilies.Random(n - block, 4, rng.Next());
            text += text.Substring(rng.NextInt(text.Length - block), block);
            PreparedText p = TextPreparation.Prepare(text);
            int[] sa = SaisReference.BuildSuffixArray(text);
            switch (i % 5)
            {
                case 0:
                    for (int k = 0, swaps = 1 + rng.NextInt(6); k < swaps; k++) { int a = 1 + rng.NextInt(n), b = 1 + rng.NextInt(n); (sa[a], sa[b]) = (sa[b], sa[a]); }
                    break;
                case 1:
                    for (int k = 0, copies = 1 + rng.NextInt(6); k < copies; k++) { int to = 1 + rng.NextInt(n), from = 1 + rng.NextInt(n); sa[to] = sa[from]; }
                    break;
                case 2:
                    { int a = 1 + rng.NextInt(n - 1); sa[a + 1] = sa[a]; }
                    break;
                case 3:
                    Array.Reverse(sa, 1 + rng.NextInt(n / 2), 2 + rng.NextInt(n / 4));
                    break;
                default:
                    for (int k = 1; k <= n; k++) sa[k] = (sa[k] + (n - 1 - sa[k] > 0 ? 1 : 0)) % n; // every position shifted by one (still in range, not a permutation of the suffixes)
                    break;
            }

            try
            {
                var report = new DirectReport { ForceTextOrder = (i & 1) == 0 };
                CompactLcp? result = Attempt(p, sa, report);
                if (result is not null)
                {
                    reached++;
                    Assert.That(result.Length, Is.EqualTo(n + 1));
                    if (report.Outcome == DirectOutcome.ResolvedInTextOrder) Assert.That(report.LongPositions, Is.LessThanOrEqualTo(report.LongRanks));
                }
            }
            catch (ArgumentException)
            {
                // an entry outside [0, n) is rejected like the reference does
            }
        }

        Assert.That(reached, Is.GreaterThan(200));
    }

    [Test]
    public void ARepeatedTextPosition_KeepsTheMinimum_AndIsCountedOnce()
    {
        string block = StringFamilies.Random(300, 4, 3);
        string text = StringFamilies.Random(3000, 4, 4) + block + StringFamilies.Random(3000, 4, 5) + block;
        PreparedText p = TextPreparation.Prepare(text);
        int[] sa = SaisReference.BuildSuffixArray(text);
        int[] lcp = Reference(p, sa).ToArray();
        int r = Array.FindIndex(lcp, 2, v => v >= 16 && v < 100); // a long rank somewhere in the middle of the repeat
        Assert.That(r, Is.GreaterThan(1));
        sa[r + 1] = sa[r]; // the same text position at two consecutive ranks: the second pair compares a suffix with itself (long, as far as the text goes)
        var report = new DirectReport { ForceTextOrder = true };
        CompactLcp? result = Attempt(p, sa, report);
        Assert.That(result, Is.Not.Null);
        Assert.That(report.LongRanks, Is.GreaterThan(report.LongPositions));
        // both ranks are long: whichever record the counting sort kept, the other one still holds the minimum of a long rank (the width), not stale memory
        int[] values = result!.ToArray();
        Assert.That(values[r], Is.GreaterThanOrEqualTo(16));
        Assert.That(values[r + 1], Is.GreaterThanOrEqualTo(16));
    }

    [Test]
    public void ArgumentChecks_AreTheOnesOfTheReference()
    {
        PreparedText p = TextPreparation.Prepare("aba");
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 3, 2, 0 }));        // wrong length
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 0, 1, 2, 3 }));     // sa[0] != n
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 2, 0, 1, 2 }));     // sa[0] != n, every other entry in range
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 3, 0, 1, 3 }));     // entry n at rank >= 1
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 3, 0, -1, 2 }));    // negative entry
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 3, 0, 1, 99 }));    // entry far outside the text
        // a bad entry far from the start, after many valid ones (the scan, not a pre-check, finds it) and in a text longer than the vector width
        string text = StringFamilies.RandomDna(500, 4);
        PreparedText q = TextPreparation.Prepare(text);
        int[] sa = SaisReference.BuildSuffixArray(text);
        sa[400] = 500;
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, q, sa));
        sa[400] = -1;
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, q, sa));
    }

    [Test]
    public void CorruptSuffixArrays_NeverLeaveTheBuffers_AndEitherFailCleanlyOrReturnAnArray()
    {
        // arbitrary permutations and arrays with repeated entries: the result is unspecified, the accesses are not (guarded arena here, bounds-checked accessors in the checked build)
        var rng = new SplitMix64(0x9753UL);
        for (int i = 0; i < 4000; i++)
        {
            int n = 1 + rng.NextInt(i < 2000 ? 40 : 400);
            string text = StringFamilies.Random(n, 1 + rng.NextInt(3), rng.Next());
            PreparedText p = TextPreparation.Prepare(text);
            int[] sa = SaisReference.BuildSuffixArray(text);
            if (i % 2 == 0)
                for (int k = n; k > 1; k--) { int r = 1 + rng.NextInt(k); (sa[k], sa[r]) = (sa[r], sa[k]); }
            else
                for (int k = 0, copies = 1 + rng.NextInt(4); k < copies; k++)
                {
                    int from = 1 + rng.NextInt(n), to = 1 + rng.NextInt(n);
                    sa[to] = sa[from];
                }
            using var arena = new NativeArena(ArenaOptions.Guarded((byte)(i & 0xFF)));
            try
            {
                LcpRun run = Engine.Run(p, sa, arena);
                Assert.That(run.Lcp.Length, Is.EqualTo(n + 1));
            }
            catch (InvalidOperationException)
            {
                // the fallback path rejects arrays whose escape counts do not match, like the reference
            }
            arena.VerifyCanaries();
        }
    }
}
