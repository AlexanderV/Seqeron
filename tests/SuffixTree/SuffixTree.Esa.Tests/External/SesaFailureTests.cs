using System.Text;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>A reservation source that behaves like a disk that fills up: the n-th creation of a file with reserved blocks fails with an <see cref="IOException"/>, the others succeed.</summary>
internal sealed class FailAtSource(int failAt) : IReservationSource
{
    private int _calls;

    public bool CanPreallocate => true;

    internal int Calls => _calls;

    public FileStream CreatePreallocated(string path, long length)
    {
        int call = Interlocked.Increment(ref _calls);
        if (call == failAt) throw new IOException($"disk full (injected at creation {call}: {Path.GetFileName(path)})");
        return new FileStream(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite);
    }

    public long AvailableFreeBytes(string directory) => -1;
}

/// <summary>The thrown exception of a hook that makes the build fail at a boundary.</summary>
public sealed class InjectedFailureException(string at) : Exception($"injected at {at}");

/// <summary>Failure modes, atomic publication and cleanup of the container writer (gate G1 of work package 1b.7e): nothing is left behind, no descriptor and no mapping stays open, a target appears only complete.</summary>
[TestFixture]
public sealed class SesaFailureTests
{
    private static string Text => StringFamilies.RandomDna(1500, 11) + new string('A', 300) + StringFamilies.RandomDna(900, 12);

    private static byte[] ReadShared(string path)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete);
        var bytes = new byte[fs.Length];
        fs.ReadExactly(bytes);
        return bytes;
    }

    private static void AssertNothingLeft(string label, string input, string target, string outDir, IoCounters? counters, bool targetMayExist = false)
    {
        string tmp = target + ".tmp";
        if (!targetMayExist) Assert.That(File.Exists(target), Is.False, $"{label}: a target was published");
        Assert.That(File.Exists(tmp), Is.False, $"{label}: the temporary file is left");
        Assert.That(Directory.Exists(tmp + ".d"), Is.False, $"{label}: the scratch directory is left");
        Assert.That(Directory.GetFileSystemEntries(outDir).Where(e => !targetMayExist || e != target), Is.Empty, $"{label}: something is left in the output directory");
        if (counters is not null) Assert.That(counters.ScratchLiveBytes, Is.Zero, $"{label}: scratch bytes are still on the account");
        ProcFds.AssertClosed(label, [input, target, tmp], [tmp + ".d"]);
        Assert.That(ProcMaps.Contains(tmp) ?? false, Is.False, $"{label}: the temporary file is still mapped");
        Assert.That(ProcMaps.Contains(target) ?? false, Is.False, $"{label}: the target is mapped");
    }

    private static (string Input, string Target, string OutDir) Setup(TempDir dir, string text, string tag = "")
    {
        string input = SesaRun.WriteInput(dir, text, $"in{tag}.txt");
        string outDir = dir.Sub($"out{tag}");
        Directory.CreateDirectory(outDir);
        return (input, Path.Combine(outDir, "x.sesa"), outDir);
    }

    private static Exception? Run(string input, string target, SesaBuildOptions options)
    {
        try
        {
            SesaBuilder.Build(input, target, options);
            return null;
        }
        catch (Exception e)
        {
            return e;
        }
    }

    [Test]
    public void ADiskThatFillsAtEveryCreationOfAFile_IsAnIoException_AndLeavesNothing()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, Text);
        var counting = new FailAtSource(int.MaxValue);
        Assert.That(Run(input, target, SesaRun.Tiny(2) with { Files = new MappedFileOptions { Source = counting } }), Is.Null);
        int creations = counting.Calls;
        Assert.That(creations, Is.GreaterThanOrEqualTo(3), "the container, the rank file and the BWT scratch are created with a reservation");
        File.Delete(target);
        for (int k = 1; k <= creations; k++)
        {
            var counters = new IoCounters();
            Exception? e = Run(input, target, SesaRun.Tiny(2) with { Counters = counters, Files = new MappedFileOptions { Source = new FailAtSource(k) } });
            Assert.That(e, Is.InstanceOf<IOException>(), $"creation {k}: {e}");
            Assert.That(e!.Message, Does.Contain("disk full"), $"creation {k}");
            AssertNothingLeft($"creation {k} of {creations}", input, target, outDir, counters);
        }
    }

    [Test]
    public void ADiskThatFillsAtEveryCreation_AlsoWithoutPsi_AndWithThreads()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, "mississippi mississippi mississippi");
        foreach ((bool psi, int threads) in new[] { (false, 3), (true, 4) })
        {
            var counting = new FailAtSource(int.MaxValue);
            Assert.That(Run(input, target, SesaRun.Tiny(5, threads, psi) with { Files = new MappedFileOptions { Source = counting } }), Is.Null);
            File.Delete(target);
            for (int k = 1; k <= counting.Calls; k++)
            {
                Exception? e = Run(input, target, SesaRun.Tiny(5, threads, psi) with { Files = new MappedFileOptions { Source = new FailAtSource(k) } });
                Assert.That(e, Is.InstanceOf<IOException>(), $"psi {psi}, creation {k}: {e}");
                AssertNothingLeft($"psi {psi}, creation {k}", input, target, outDir, null);
            }
        }
    }

    [Test]
    public void ARealFullDisk_AtTheContainer_IsRefusedAtCreation()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, Text);
        var small = new FakeReservationSource { CanPreallocate = false, Free = 1000 };
        Exception? e = Run(input, target, SesaRun.Tiny(1) with { Files = new MappedFileOptions { Source = small } });
        Assert.That(e, Is.InstanceOf<IOException>(), e?.ToString());
        Assert.That(e!.Message, Does.Contain("free space"));
        AssertNothingLeft("free space 1000", input, target, outDir, null);
    }

    [Test]
    public void ThePreflight_RefusesATextWhoseScratchDoesNotFit_BeforeAnythingIsCreated()
    {
        using var dir = new TempDir();
        string text = "mississippi";
        (string input, string target, string outDir) = Setup(dir, text);
        PackedTextInfo info = PackedTextWriter.Scan(input);
        SesaPlan plan = SesaLayout.Plan(info, true);
        long need = plan.ReservedLength + (12 * info.N);
        var tight = new FakeReservationSource { CanPreallocate = true, Free = need - 1 };
        Exception? e = Run(input, target, SesaRun.Tiny(1) with { Files = new MappedFileOptions { Source = tight } });
        Assert.That(e, Is.InstanceOf<IOException>(), e?.ToString());
        Assert.That(e!.Message, Does.Contain("free space").And.Contain(need.ToString(System.Globalization.CultureInfo.InvariantCulture)));
        Assert.That(tight.PreallocateCalls, Is.Zero, "nothing was created");
        AssertNothingLeft("preflight", input, target, outDir, null);
        var enough = new FakeReservationSource { CanPreallocate = true, Free = need };
        Assert.That(Run(input, target, SesaRun.Tiny(1) with { Files = new MappedFileOptions { Source = enough } }), Is.Null, "exactly enough is enough");
    }

    [Test]
    public void TheLatePreflight_AfterB2_CountsTheGrowthAndTheSpillBoundOfTheStack_ExactlyToTheByte()
    {
        using var dir = new TempDir();
        string text = new('a', 2000);
        (string input, string target, string outDir) = Setup(dir, text);
        SesaBuildResult baseline = SesaBuilder.Build(input, target, SesaRun.Tiny(1));
        File.Delete(target);
        long need = baseline.GrownBytes + (16L * (baseline.MaxLcp + 1));
        Assert.That(baseline.GrownBytes, Is.GreaterThan(0), "a^n has an LCP wider than one byte, so the file grows");
        Assert.That(baseline.MaxLcp, Is.EqualTo(1999));
        foreach ((long free, bool fits) in new[] { (need - 1, false), (need, true), (need + 1, true), (0L, false) })
        {
            var source = new FakeReservationSource { CanPreallocate = true, Free = 1L << 40 };
            SesaBuildOptions options = SesaRun.Tiny(1) with
            {
                Files = new MappedFileOptions { Source = source },
                LcpTrace = name =>
                {
                    if (name == "b3") source.Free = free; // after B2: the width, the growth and the depth of the stack are known
                },
            };
            Exception? e = Run(input, target, options);
            Assert.That(source.Asked, Is.Not.Empty.And.All.EqualTo(Path.GetDirectoryName(Path.GetFullPath(target))), $"free {free}: the free space of the directory of the target is what is asked, early (before the creation), late (before the growth) and by the growth itself");
            if (fits)
            {
                Assert.That(e, Is.Null, $"free {free} = need {need} + {free - need}: {e}");
                File.Delete(target);
            }
            else
            {
                Assert.That(e, Is.InstanceOf<IOException>(), $"free {free}: {e}");
                Assert.That(e!.Message, Does.Contain("free space"));
                AssertNothingLeft("late preflight", input, target, outDir, null);
            }
        }
    }

    [Test]
    public void AFailureAtEveryHookPoint_LeavesNothing_AndKeepsTheTypeOfTheException()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, Text);
        foreach (string point in SesaHooks.All)
        {
            foreach (int threads in new[] { 1, 3 })
            {
                var counters = new IoCounters();
                SesaBuildOptions options = SesaRun.Tiny(threads, threads) with { Counters = counters, Hook = name => { if (name == point) throw new InjectedFailureException(point); } };
                Exception? e = Run(input, target, options);
                Assert.That(e, Is.InstanceOf<InjectedFailureException>(), $"{point}: {e}");
                AssertNothingLeft($"hook {point}, {threads} threads", input, target, outDir, counters);
            }
        }
    }

    [Test]
    public void AFailureInsideAStage_LeavesNothing_AtEveryTracePoint()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, Text);
        var names = new List<string>();
        Assert.That(Run(input, target, SesaRun.Tiny(3) with { SaTrace = names.Add, LcpTrace = names.Add }), Is.Null);
        File.Delete(target);
        Assert.That(names, Is.Not.Empty);
        foreach (string point in names.Distinct())
        {
            foreach (bool inLcp in new[] { false, true })
            {
                var counters = new IoCounters();
                Action<string> hook = name => { if (name == point) throw new InjectedFailureException(point); };
                SesaBuildOptions options = SesaRun.Tiny(3) with { Counters = counters, SaTrace = inLcp ? null : hook, LcpTrace = inLcp ? hook : null };
                Exception? e = Run(input, target, options);
                if (e is null)
                {
                    File.Delete(target);
                    continue; // the point belongs to the other stage
                }
                Assert.That(e, Is.InstanceOf<InjectedFailureException>(), $"{point}: {e}");
                AssertNothingLeft($"trace {point}", input, target, outDir, counters);
            }
        }
    }

    [Test]
    public void AnArenaCapThatIsTooSmall_FailsWithTheBudgetException_AtEveryStage_AndLeavesNothing()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, Text);
        (string Label, Func<SesaBuildOptions, SesaBuildOptions> Make)[] caps =
        [
            ("PH1 and PH3 cap", o => o with { ArenaCapBytes = 1000 }),
            ("PH5 cap", o => o with { ChildArenaCapBytes = 1 }),
            ("PH6 cap", o => o with { Ph6CapBytes = 1 }),
        ];
        foreach ((string label, Func<SesaBuildOptions, SesaBuildOptions> make) in caps)
        {
            var counters = new IoCounters();
            Exception? e = Run(input, target, make(SesaRun.Tiny(2) with { Counters = counters }));
            Assert.That(e, Is.InstanceOf<ArenaBudgetExceededException>(), $"{label}: {e}");
            AssertNothingLeft(label, input, target, outDir, counters);
        }
    }

    [Test]
    public void AnInputThatShrinksBetweenThePasses_IsRefused_AndLeavesNothing()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, Text);
        foreach (string change in new[] { "truncate", "grow", "other-symbol" })
        {
            File.WriteAllBytes(input, Encoding.Latin1.GetBytes(Text));
            SesaBuildOptions options = SesaRun.Tiny(2) with
            {
                BetweenPasses = () =>
                {
                    byte[] bytes = File.ReadAllBytes(input);
                    File.WriteAllBytes(input, change switch { "truncate" => bytes[..^50], "grow" => [.. bytes, .. bytes[..10]], _ => [.. bytes[..^1], (byte)'Z'] });
                },
            };
            Exception? e = Run(input, target, options);
            Assert.That(e, Is.InstanceOf<InvalidOperationException>(), $"{change}: {e}");
            Assert.That(e!.Message, Does.Contain("changed"));
            AssertNothingLeft(change, input, target, outDir, null);
        }
    }

    [Test]
    public void AMissingInput_AMissingTargetDirectory_AndAnInputAboveTheLimit_AreRefusedBeforeAnythingIsCreated()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, "abc");
        Assert.Throws<FileNotFoundException>(() => SesaBuilder.Build(dir.File("nothing.txt"), target, SesaRun.Tiny(1)));
        AssertNothingLeft("missing input", input, target, outDir, null);
        string missing = Path.Combine(dir.Path, "no-such-dir", "x.sesa");
        Assert.Throws<DirectoryNotFoundException>(() => SesaBuilder.Build(input, missing, SesaRun.Tiny(1)));
        Assert.That(Directory.Exists(Path.Combine(dir.Path, "no-such-dir")), Is.False);
        string sparse = dir.File("huge.txt");
        using (var fs = new FileStream(sparse, FileMode.Create)) fs.SetLength(PackedTextFormat.MaxSymbols + 1);
        Assert.Throws<NotSupportedException>(() => SesaBuilder.Build(sparse, target, SesaRun.Tiny(1)));
        AssertNothingLeft("above the limit", input, target, outDir, null);
        Assert.Throws<ArgumentNullException>(() => SesaBuilder.Build(null!, target));
        Assert.Throws<ArgumentNullException>(() => SesaBuilder.Build(input, null!));
    }

    [Test]
    public void BetweenTheHeaderAndTheMagic_TheTemporaryFileHasNoMagic_AndThereIsNoTarget()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, Text);
        string tmp = target + ".tmp";
        var snapshots = new Dictionary<string, byte[]>();
        var states = new Dictionary<string, (bool TargetExists, bool ScratchExists)>();
        SesaBuildOptions options = SesaRun.Tiny(4) with
        {
            Fsync = true,
            Hook = name =>
            {
                if (name is SesaHooks.Created or SesaHooks.Packed or SesaHooks.Checksummed or SesaHooks.HeaderWritten or SesaHooks.HeaderSynced or SesaHooks.MagicWritten)
                {
                    snapshots[name] = ReadShared(tmp);
                    states[name] = (File.Exists(target), Directory.Exists(tmp + ".d"));
                }
            },
        };
        Assert.That(Run(input, target, options), Is.Null);
        byte[] published = File.ReadAllBytes(target);
        byte[] magic = [0x89, 0x53, 0x45, 0x53, 0x41, 0x0D, 0x0A, 0x1A];
        foreach (string name in new[] { SesaHooks.Created, SesaHooks.Packed, SesaHooks.Checksummed, SesaHooks.HeaderWritten, SesaHooks.HeaderSynced })
        {
            Assert.That(snapshots[name].AsSpan(0, 8).ToArray(), Is.EqualTo(new byte[8]), $"{name}: the magic is written last");
            Assert.That(states[name].TargetExists, Is.False, $"{name}: no target before the rename");
            Assert.That(states[name].ScratchExists, Is.True, $"{name}: the scratch directory exists while the build runs");
        }
        Assert.That(snapshots[SesaHooks.MagicWritten].AsSpan(0, 8).ToArray(), Is.EqualTo(magic), "the magic is there after the second sync");
        Assert.That(states[SesaHooks.MagicWritten].TargetExists, Is.False, "the rename is after the magic");
        foreach (string name in new[] { SesaHooks.HeaderWritten, SesaHooks.HeaderSynced })
        {
            SesaOracle o = SesaOracle.FromBytes(snapshots[name]);
            byte[] withMagic = (byte[])snapshots[name].Clone();
            magic.CopyTo(withMagic, 0);
            Assert.That(snapshots[name].AsSpan(8).SequenceEqual(published.AsSpan(8)), Is.True, $"{name}: everything but the magic is already the final file");
            SesaOracle.FromBytes(withMagic).AssertConforms($"the file at {name} plus the magic", true);
            Assert.That(o.HeaderCrc, Is.EqualTo(SesaOracle.FromBytes(published).HeaderCrc));
        }
        Assert.That(snapshots[SesaHooks.Checksummed].AsSpan(0, 128).ToArray(), Is.EqualTo(new byte[128]), "the header page is written after the checksums");
        Assert.That(snapshots[SesaHooks.Packed].AsSpan(4096).ToArray().Any(b => b != 0), Is.True);
        Assert.That(File.Exists(tmp), Is.False);
        Assert.That(Directory.Exists(tmp + ".d"), Is.False);
        SesaOracle.FromBytes(published).AssertConforms("published", true);
    }

    [Test]
    public void ABuildThatIsHeldAtTheHeaderPoint_CanBeAbandoned_AndALaterBuildReplacesTheResidue()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, Text);
        string tmp = target + ".tmp";
        using var reached = new ManualResetEventSlim();
        using var release = new ManualResetEventSlim();
        byte[]? snapshot = null;
        SesaBuildOptions held = SesaRun.Tiny(2) with
        {
            Hook = name =>
            {
                if (name != SesaHooks.HeaderSynced) return;
                snapshot = ReadShared(tmp);
                reached.Set();
                release.Wait(TimeSpan.FromSeconds(60));
                throw new InjectedFailureException("abandoned");
            },
        };
        Task<Exception?> first = Task.Run(() => Run(input, target, held));
        SpinWait.SpinUntil(() => reached.IsSet || first.IsCompleted, TimeSpan.FromSeconds(60));
        Assert.That(reached.IsSet, Is.True, $"the build never reached the point: {(first.IsCompleted ? first.Result : null)}");
        // the state of a hard kill at this point: the temporary file stays, without magic, and there is no target
        Assert.That(File.Exists(target), Is.False);
        Assert.That(snapshot!.AsSpan(0, 8).ToArray(), Is.EqualTo(new byte[8]));
        release.Set();
        Assert.That(first.Wait(TimeSpan.FromSeconds(60)), Is.True);
        Assert.That(first.Result, Is.InstanceOf<InjectedFailureException>());
        AssertNothingLeft("abandoned", input, target, outDir, null);
        // a residue of a hard kill is replaced by the next build
        File.WriteAllBytes(tmp, snapshot!);
        Directory.CreateDirectory(tmp + ".d/old/deeper");
        File.WriteAllText(tmp + ".d/old/deeper/x.bin", "stale");
        File.WriteAllText(tmp + ".d/y.bin", "stale");
        Assert.That(Run(input, target, SesaRun.Tiny(2)), Is.Null);
        AssertNothingLeft("after the residue", input, target, outDir, null, targetMayExist: true);
        SesaOracle.Read(target).AssertConforms("rebuilt", true);
    }

    [Test]
    public void AnExistingTarget_StaysUntouchedByAFailedBuild_AndIsReplacedByASuccessfulOne()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, Text);
        File.WriteAllText(target, "OLD CONTENT");
        Exception? e = Run(input, target, SesaRun.Tiny(2) with { Hook = name => { if (name == SesaHooks.PsiBuilt) throw new InjectedFailureException(name); } });
        Assert.That(e, Is.InstanceOf<InjectedFailureException>());
        Assert.That(File.ReadAllText(target), Is.EqualTo("OLD CONTENT"), "a failed build leaves a good old file alone");
        AssertNothingLeft("old target", input, target, outDir, null, targetMayExist: true);
        Assert.That(Run(input, target, SesaRun.Tiny(2)), Is.Null);
        SesaOracle.Read(target).AssertConforms("replaced", true);
        AssertNothingLeft("replaced", input, target, outDir, null, targetMayExist: true);
    }

    [Test]
    public void ATargetThatIsADirectory_FailsAtTheRename_AndLeavesNoTemporaryFile()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, "mississippi");
        Directory.CreateDirectory(target);
        File.WriteAllText(Path.Combine(target, "keep.txt"), "x");
        Exception? e = Run(input, target, SesaRun.Tiny(2));
        Assert.That(e, Is.InstanceOf<IOException>(), e?.ToString());
        Assert.That(File.Exists(Path.Combine(target, "keep.txt")), Is.True);
        Assert.That(File.Exists(target + ".tmp"), Is.False);
        Assert.That(Directory.Exists(target + ".tmp.d"), Is.False);
        ProcFds.AssertClosed("directory target", [input, target + ".tmp"], [target + ".tmp.d"]);
    }

    [Test]
    public void AScratchDirectoryThatCannotBeCreated_FailsCleanly()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, "mississippi");
        File.WriteAllText(target + ".tmp.d", "a file where the scratch directory goes");
        Exception? e = Run(input, target, SesaRun.Tiny(2));
        // a stale entry of that name is ours to replace: the build succeeds and the file is gone
        Assert.That(e, Is.Null, e?.ToString());
        Assert.That(File.Exists(target + ".tmp.d"), Is.False);
        SesaOracle.Read(target).AssertConforms("scratch was a file", true);
    }

    [Test]
    public void TheDisposalOfEveryResource_HoldsAfterASuccess_Too()
    {
        using var dir = new TempDir();
        (string input, string target, string outDir) = Setup(dir, Text);
        var counters = new IoCounters();
        Assert.That(Run(input, target, SesaRun.Tiny(2, 3) with { Counters = counters }), Is.Null);
        AssertNothingLeft("success", input, target, outDir, counters, targetMayExist: true);
        Assert.That(ProcMaps.Contains(target + ".tmp") ?? false, Is.False);
    }
}
