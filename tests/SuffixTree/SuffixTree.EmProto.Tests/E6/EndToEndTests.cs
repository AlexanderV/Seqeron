using SplitMix64 = SuffixTree.Esa.Tests.TestSupport.SplitMix64;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph5;
using SuffixTree.EmProto.Pipeline;
using SuffixTree.EmProto.Tests.Ph5;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.WholePipeline;

/// <summary>
/// The whole pipeline PH0 to PH5 of the prototype E6 (ADR-004 section 5, row E6, gate G1) against the in-RAM structures of the product: small texts, the families, exhaustive strings, identity over threads and over
/// K, the order of the phases, the deletion of the scratch, the final length of every output, the counters of the phases.
/// </summary>
[TestFixture]
public class EndToEndTests : TempDirFixture
{
    [Test]
    public void Banana_And_Friends_AreEqualToTheInRamStructures_AndTheirEnumerationToTheUkkonenHash()
    {
        int v = 0;
        foreach (string s in new[] { "a", "aa", "ab", "aba", "banana", "mississippi", "abracadabra", "aaaaaaaaab", "abababababab", "\0ÿ\0ÿ" })
        {
            ChildReference r = ChildReference.Of(s);
            PipelineChecks.RunAndCompare(r, PipelineChecks.Options(Run($"c{v}"), s, v), $"'{s}' variant {v}", childTraversal: true);
            v++;
        }
    }

    [Test]
    public void Exhaustive_Binary1To12_Ternary1To6_WholePipeline_TinyParameters()
    {
        int count = 0;
        foreach ((int sigma, int max) in new[] { (2, 12), (3, 6) })
        {
            foreach (string s in ExhaustiveStrings.All(sigma, max))
            {
                if (s.Length == 0) continue;
                ChildReference r = ChildReference.Of(s);
                string d = Run("x");
                PipelineChecks.RunAndCompare(r, PipelineChecks.Options(d, s, count, 1 + (count % 3)), $"'{s}' variant {count}");
                Directory.Delete(d, true);
                count++;
            }
        }
        Assert.That(count, Is.EqualTo(8190 + 1092), "binary 1-12 and ternary 1-6");
    }

    [Test]
    public void TinyTexts_EnumerationThroughChild_EqualsTheInRamIndex_AndTheUkkonenLogicalHash()
    {
        int count = 0;
        var texts = new List<string> { "banana", "mississippi", "abracadabra", "aaaa", "abababab", "aabaabaab" };
        foreach ((int sigma, int max) in new[] { (2, 7), (3, 4) })
            texts.AddRange(ExhaustiveStrings.All(sigma, max).Where(s => s.Length > 0));
        foreach (string s in texts)
        {
            string d = Run("t");
            PipelineChecks.RunAndCompare(ChildReference.Of(s), PipelineChecks.Options(d, s, count), $"'{s}'", childTraversal: true);
            Directory.Delete(d, true);
            count++;
        }
        Assert.That(count, Is.GreaterThan(250));
    }

    [Test]
    public void AdversarialFamilies_1000_And_5000_TinyParameters()
    {
        int cases = 0;
        foreach (int n in new[] { 1000, 5000 })
        {
            foreach ((string name, string text) in StringFamilies.Adversarial(n, 3))
            {
                if (text.Any(c => c > 255) || text.Length == 0) continue;
                string d = Run("f");
                PipelineChecks.RunAndCompare(ChildReference.Of(text), PipelineChecks.Options(d, text, cases, 1 + (cases % 4)), $"{name}-n{n} variant {cases}");
                Directory.Delete(d, true);
                cases++;
            }
        }
        Assert.That(cases, Is.GreaterThan(40));
    }

    [Test]
    public void AdversarialFamilies_At100000_ModerateParameters()
    {
        int cases = 0;
        foreach ((string name, string text) in StringFamilies.Adversarial(100_000, 7))
        {
            if (text.Any(c => c > 255) || text.Length == 0) continue;
            string d = Run("f");
            EndToEndOptions o = PipelineChecks.Options(d, text, 0, 4) with
            {
                Period = 64,
                Tuning = new EndToEndTuning
                {
                    BucketCapacity = 2000,
                    DcsBucketCapacity = 1500,
                    DcsBaseCase = 3000,
                    DcsChunkSlots = 700,
                    DcsSortCapacity = 900,
                    LcpChunk = 9000,
                    BlockEntries = 4096,
                    ChildWindowBytes = 1 << 14,
                    ChildSlots = 4,
                },
            };
            PipelineChecks.RunAndCompare(ChildReference.Of(text), o, $"{name}-n100000");
            Directory.Delete(d, true);
            cases++;
        }
        Assert.That(cases, Is.GreaterThan(25));
    }

    [Test]
    public void TheDeepFamilies_At1000000_AreEqual()
    {
        foreach (string family in new[] { "a", "fib", "tm" })
        {
            string text = family switch
            {
                "a" => StringFamilies.AllSame(1_000_000),
                "fib" => StringFamilies.Fibonacci(1_000_000),
                _ => StringFamilies.ThueMorse(1_000_000),
            };
            string d = Run(family);
            EndToEndOptions o = PipelineChecks.Options(d, text, 0, 4) with
            {
                Period = 64,
                Tuning = new EndToEndTuning
                {
                    BucketCapacity = 20_000,
                    DcsBucketCapacity = 15_000,
                    DcsBaseCase = 30_000,
                    DcsChunkSlots = 7000,
                    DcsSortCapacity = 9000,
                    LcpChunk = 90_000,
                    BlockEntries = 8192,
                    ChildWindowBytes = 1 << 16,
                    ChildSlots = 4,
                    WindowBytes = 1 << 16,
                    DcsIoBytes = 4096,
                    LcpReadBytes = 1 << 14,
                    LcpTupleBytes = 1 << 12,
                    SpillFlushBytes = 1 << 20,
                },
            };
            PipelineChecks.RunAndCompare(ChildReference.Of(text), o, $"{family}-n1000000");
            Directory.Delete(d, true);
        }
    }

    [Test]
    public void OutputIsIdentical_ForOneTwoAndFourThreads_AndForTwoValuesOfK_AndEqualsTheInRamStructures()
    {
        foreach ((string name, string text) in new[]
        {
            ("random-dna", StringFamilies.RandomDna(100_000, 5)),
            ("periodic-7", StringFamilies.Periodic(100_000, 7, 3)),
            ("fibonacci", StringFamilies.Fibonacci(100_000)),
            ("runs-L256", StringFamilies.Runs(100_000, 256, 4)),
            ("all-same", StringFamilies.AllSame(60_000)),
        })
        {
            ChildReference r = ChildReference.Of(text);
            Dictionary<string, string>? first = null;
            string firstLabel = "";
            foreach ((int threads, int k, string label) in new[] { (1, 1000, "1 thread K=1000"), (2, 1000, "2 threads K=1000"), (4, 1000, "4 threads K=1000"), (4, 500, "4 threads K=500"), (3, 1900, "3 threads K=1900") })
            {
                string d = Run("i");
                EndToEndOptions o = PipelineChecks.Options(d, text, 2, threads);
                o = o with
                {
                    Period = 64,
                    Tuning = o.Tuning with
                    {
                        BucketCapacity = k,
                        DcsBucketCapacity = k / 2,
                        DcsBaseCase = 3000,
                        DcsSortCapacity = k / 2,
                        DcsChunkSlots = 500 + k,
                        LcpChunk = 3000 + (2 * k),
                        LcpWorkers = null,
                        BlockEntries = 64,
                    },
                };
                PipelineChecks.RunAndCompare(r, o, $"{name} {label}");
                Dictionary<string, string> hashes = PipelineChecks.OutputHashes(o.OutPrefix);
                if (first is null)
                {
                    first = hashes;
                    firstLabel = label;
                }
                else
                {
                    foreach (string ext in first.Keys) Assert.That(hashes[ext], Is.EqualTo(first[ext]), $"{name}: .{ext} of '{label}' differs from '{firstLabel}'");
                }
                Directory.Delete(d, true);
            }
        }
    }

    [Test]
    public void ThePhases_RunInOrder_AndTheScratchIsReleasedStageByStage()
    {
        string text = StringFamilies.RandomDna(20_000, 9);
        string d = Run("p");
        EndToEndOptions o = PipelineChecks.Options(d, text, 3, 2) with { SyncOutputs = true };
        string scratch = o.ScratchDir;
        string ranks = Path.Combine(scratch, "ranks.bin");
        var events = new List<string>();
        EndToEndOptions traced = o with
        {
            Trace = (what, phase) =>
            {
                events.Add($"{what} {phase}");
                bool isBegin = what == "begin";
                switch (phase)
                {
                    case EndToEnd.PhasePrepare when isBegin:
                        Assert.That(Directory.GetFileSystemEntries(scratch), Is.Empty, "the scratch directory exists and is empty when the run starts");
                        Assert.That(File.Exists(o.OutPrefix + ".pk"), Is.False);
                        break;
                    case EndToEnd.PhaseDcs when isBegin:
                        Assert.That(new FileInfo(o.OutPrefix + ".pk").Length, Is.EqualTo(((20_000 * 2) + 7) / 8 + 64), "the packed text is final before PH1");
                        Assert.That(File.Exists(o.OutPrefix + ".pk.meta"), Is.True);
                        Assert.That(File.Exists(ranks), Is.False);
                        break;
                    case EndToEnd.PhaseSa when isBegin:
                        Assert.That(File.Exists(ranks), Is.True, "PH1 leaves the rank file for PH3");
                        Assert.That(ScratchTracker.Live, Is.EqualTo(new FileInfo(ranks).Length), "the rank file is the only tracked scratch between PH1 and PH3");
                        Assert.That(Directory.GetFileSystemEntries(Path.Combine(scratch, "ph1")), Is.Empty);
                        Assert.That(File.Exists(o.OutPrefix + ".sa"), Is.False);
                        break;
                    case EndToEnd.PhaseLcpChild when isBegin:
                        Assert.That(File.Exists(ranks), Is.False, "the rank file is deleted as soon as PH3 is done");
                        Assert.That(ScratchTracker.Live, Is.Zero);
                        Assert.That(new FileInfo(o.OutPrefix + ".sa").Length, Is.EqualTo(4L * 20_001), "the suffix array is final before PH4");
                        Assert.That(Directory.GetFileSystemEntries(Path.Combine(scratch, "ph3")), Is.Empty);
                        Assert.That(File.Exists(o.OutPrefix + ".lcp"), Is.False, "the LCP files appear only when PH4 knows the width");
                        Assert.That(File.Exists(o.OutPrefix + ".cld"), Is.False);
                        break;
                    case EndToEnd.PhaseSync when isBegin:
                        Assert.That(File.Exists(o.OutPrefix + ".stats"), Is.True);
                        Assert.That(ScratchTracker.Live, Is.Zero);
                        Assert.That(Directory.Exists(Path.Combine(scratch, "ph4")) ? Directory.GetFileSystemEntries(Path.Combine(scratch, "ph4")) : [], Is.Empty, "PH4 scratch is gone when PH4 ends");
                        break;
                    default:
                        break;
                }
            },
        };
        EndToEndResult res = EndToEnd.Run(traced);
        Assert.That(events, Is.EqualTo(new[]
        {
            "begin ph0_prepare", "end ph0_prepare", "begin ph1_dcs", "end ph1_dcs", "begin ph23_sa", "end ph23_sa", "begin ph45_lcp_child", "end ph45_lcp_child", "begin sync", "end sync",
        }));
        Assert.That(res.Phases.Select(x => x.Name), Is.EqualTo(new[] { EndToEnd.PhasePrepare, EndToEnd.PhaseDcs, EndToEnd.PhaseSa, EndToEnd.PhaseLcpChild, EndToEnd.PhaseSync }));
        Assert.That(res.Clean, Is.True);
        Assert.That(res.StageSeconds, Is.EqualTo(res.Phases.Take(4).Sum(x => x.Seconds)).Within(1e-9), "the stage time excludes the sync");
        Assert.That(res.SaStageSeconds, Is.EqualTo(res.Phases.Take(3).Sum(x => x.Seconds)).Within(1e-9), "the SA stage is PH0 to PH3");
    }

    [Test]
    public void NoSync_MeansNoSyncPhase()
    {
        string d = Run("n");
        EndToEndResult res = EndToEnd.Run(PipelineChecks.Options(d, "mississippi", 1));
        Assert.That(res.Phases.Select(x => x.Name), Is.EqualTo(EndToEnd.StagePhases));
    }

    [TestCase(EndToEnd.PhasePrepare, "begin")]
    [TestCase(EndToEnd.PhasePrepare, "end")]
    [TestCase(EndToEnd.PhaseDcs, "begin")]
    [TestCase(EndToEnd.PhaseDcs, "end")]
    [TestCase(EndToEnd.PhaseSa, "begin")]
    [TestCase(EndToEnd.PhaseSa, "end")]
    [TestCase(EndToEnd.PhaseLcpChild, "begin")]
    [TestCase(EndToEnd.PhaseLcpChild, "end")]
    [TestCase(EndToEnd.PhaseSync, "begin")]
    public void AFailingStage_PropagatesTheException_AndLeavesNoScratch(string phase, string when)
    {
        string d = Run("e");
        EndToEndOptions o = PipelineChecks.Options(d, StringFamilies.RandomDna(3000, 2), 4, 2) with
        {
            SyncOutputs = true,
            Trace = (w, p) =>
            {
                if (w == when && p == phase) throw new InvalidOperationException("injected failure");
            },
        };
        Assert.Throws<InvalidOperationException>(() => EndToEnd.Run(o));
        Assert.That(Directory.Exists(o.ScratchDir), Is.False, "a failed run deletes its scratch directory");
    }

    [Test]
    public void ANonEmptyScratchDirectory_IsRefused_AndLeftAlone()
    {
        string d = Run("s");
        EndToEndOptions o = PipelineChecks.Options(d, "banana", 0);
        Directory.CreateDirectory(o.ScratchDir);
        File.WriteAllText(Path.Combine(o.ScratchDir, "keep.me"), "x");
        Assert.Throws<ArgumentException>(() => EndToEnd.Run(o));
        Assert.That(File.Exists(Path.Combine(o.ScratchDir, "keep.me")), Is.True);
        Assert.That(File.Exists(o.OutPrefix + ".pk"), Is.False, "nothing was written");
    }

    [Test]
    public void AnEmptyInput_IsRefused_BeforeAnythingIsWritten()
    {
        string d = Run("z");
        EndToEndOptions o = PipelineChecks.Options(d, "", 0);
        Assert.Throws<ArgumentException>(() => EndToEnd.Run(o));
        Assert.That(Directory.Exists(o.ScratchDir), Is.False);
        Assert.That(File.Exists(o.OutPrefix + ".pk"), Is.False);
    }

    [Test]
    public void AnInputThatIsASymbolicLink_IsTheWholeText()
    {
        string text = StringFamilies.RandomDna(5000, 17);
        string d = Run("l");
        EndToEndOptions o = PipelineChecks.Options(d, text, 2, 2);
        string link = Path.Combine(d, "link.txt");
        File.CreateSymbolicLink(link, o.InputPath);
        PipelineChecks.RunAndCompare(ChildReference.Of(text), o with { InputPath = link }, "symbolic link");
    }

    [Test]
    public void AMissingInput_IsRefused()
    {
        string d = Run("m");
        EndToEndOptions o = PipelineChecks.Options(d, "banana", 0) with { InputPath = Path.Combine(d, "nope.txt") };
        Assert.Throws<FileNotFoundException>(() => EndToEnd.Run(o));
    }

    [Test]
    public void TheCountersOfThePhases_AccountForTheInputTheOutputsAndTheMappedWrites()
    {
        const int n = 30_000;
        string text = StringFamilies.RandomDna(n, 13);
        string d = Run("c");
        EndToEndResult res = EndToEnd.Run(PipelineChecks.Options(d, text, 5, 2));
        PhaseRecord ph0 = res.Phase(EndToEnd.PhasePrepare);
        Assert.That(ph0.ExplicitReadBytes, Is.EqualTo(2L * n), "PH0 reads the input twice (alphabet scan, then the pack)");
        Assert.That(ph0.MappedWriteBytes, Is.EqualTo((n * 2L) / 8 + (((n * 2L) % 8) != 0 ? 1 : 0)), "the packed text goes through the mapped windows");
        Assert.That(ph0.ExplicitWriteBytes, Is.Zero);
        PhaseRecord ph1 = res.Phase(EndToEnd.PhaseDcs);
        Assert.That(ph1.MappedWriteBytes, Is.GreaterThanOrEqualTo(res.Dcs.Levels[0].SampleSize * 4), "the rank file");
        PhaseRecord ph3 = res.Phase(EndToEnd.PhaseSa);
        Assert.That(ph3.MappedWriteBytes, Is.EqualTo(4L * (n + 1)), "the SA goes through the mapped windows exactly once");
        Assert.That(ph3.ExplicitReadBytes, Is.GreaterThan(0));
        PhaseRecord ph4 = res.Phase(EndToEnd.PhaseLcpChild);
        Assert.That(ph4.MappedWriteBytes, Is.EqualTo(new FileInfo(res.Files.Lcp).Length + new FileInfo(res.Files.Overflow).Length), "the LCP files go through the mapped windows exactly once");
        Assert.That(ph4.ExplicitReadBytes, Is.GreaterThanOrEqualTo(4L * (n + 1)), "PH4 reads the suffix array");
        Assert.That(res.Phases.All(x => x.Seconds > 0), Is.True);
        Assert.That(res.ScratchPeakBytes, Is.GreaterThanOrEqualTo(12L * n), "PH4 base mode peaks at 12 B/char of tuples");
        Assert.That(res.ScratchPeakBytes, Is.LessThan(13L * n), "and not more");
        Assert.That(res.ScratchLiveBytes, Is.Zero);
        Assert.That(res.ArenaPeakBytes, Is.GreaterThan(0));
        Assert.That(res.Lcp.Rounds, Is.EqualTo(1), "PH4 runs in base mode");
        Assert.That(res.Parameters.N, Is.EqualTo(n));
        Assert.That(res.CoverPeriod, Is.EqualTo(PipelineChecks.Period(5)));
    }

    [Test]
    public void TheFusedPhase_FeedsPh5ThroughTheHook_EveryRankOnce_InOrder()
    {
        // the child table and the statistics are only right if every LCP chunk, the last one included, reaches PH5 exactly once and in order
        string text = StringFamilies.Periodic(5000, 3, 5) + "x";
        string d = Run("h");
        EndToEndOptions o = PipelineChecks.Options(d, text, 7, 3);
        ChildReference r = ChildReference.Of(text);
        EndToEndResult res = PipelineChecks.RunAndCompare(r, o, "hook");
        Assert.That(res.Lcp.RankChunks, Is.GreaterThan(5), "several chunks, so the chunk boundaries are crossed");
        Assert.That(res.Scan.NodeCount, Is.EqualTo(r.Statistics.NodeCount));
    }

    [Test]
    public void LargeAlphabet_AllBytes_Bits8()
    {
        var chars = new char[1200];
        var rng = new SplitMix64(77);
        for (int i = 0; i < chars.Length; i++) chars[i] = (char)rng.NextInt(256);
        for (int i = 0; i < 256; i++) chars[i * 3] = (char)i; // every byte value occurs
        string text = new(chars);
        string d = Run("w");
        PipelineChecks.RunAndCompare(ChildReference.Of(text), PipelineChecks.Options(d, text, 3, 2), "all 256 bytes");
    }

    [Test]
    public void ProteinLikeAlphabet_Bits5()
    {
        string text = StringFamilies.RandomProtein(4000, 3);
        string d = Run("w");
        EndToEndResult res = PipelineChecks.RunAndCompare(ChildReference.Of(text), PipelineChecks.Options(d, text, 6, 2), "protein");
        Assert.That(res.Prepare.Bits, Is.EqualTo(5));
        Assert.That(res.Prepare.Sigma, Is.EqualTo(20));
    }

    [Test]
    public void TaggedCodec_IsSupportedToo()
    {
        string text = StringFamilies.RandomDna(3000, 4);
        string d = Run("w");
        EndToEndOptions o = PipelineChecks.Options(d, text, 1, 2) with { Codec = ChildCodec.Tagged };
        PipelineChecks.RunAndCompare(ChildReference.Of(text), o, "tagged");
    }
}
