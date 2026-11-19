using System.Diagnostics;
using System.Text;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>What a case builds from one text that is not a function of the text alone: the seed-derived knobs of the container builder and of the layout of the in-memory trees.</summary>
/// <param name="Seed">The seed of the case (the probes of the script are drawn from it).</param>
/// <param name="Variant">Rotates the tiny parameters of the builder, the period of the difference cover and the way a tree is created.</param>
/// <param name="Psi">Whether the container has the optional PSI section (the builder's <c>IncludePsi</c>; the harness builds with and without it, the default is not changed).</param>
internal sealed record FullFuzzConfig(ulong Seed, int Variant, bool Psi)
{
    /// <summary>The configuration of a case: a pure function of its seed.</summary>
    internal static FullFuzzConfig Of(FuzzCase c) => Of(c.Seed);

    /// <summary>The configuration for a seed.</summary>
    internal static FullFuzzConfig Of(ulong seed) => new(seed, (int)(seed % 4093), ((seed >> 29) & 1) == 0);
}

/// <summary>An engine under test: its name (as in <see cref="FuzzSubjects.EngineNames"/>), the counted tree, whether it is file backed.</summary>
internal sealed record FuzzSubject(string Engine, CountingTree Tree, bool IsFile, long BuildTicks);

/// <summary>
/// The engines of one text (work package 1b.8): the in-memory trees (the child layouts), the file-backed trees (a container built by the product from the text, opened at the Header and the Structure level in the
/// Checked and the Trusted mode), and the references (the Ukkonen tree and the persistent v6 tree). A text with a code unit above 255 cannot be held by the container (one byte per symbol): it has no file-backed
/// engines and the case is counted as skipped for them.
/// </summary>
internal sealed class FuzzSubjects : IDisposable
{
    /// <summary>The names of the engines under test, in the order of the counters.</summary>
    internal static readonly string[] EngineNames =
    [
        "memory-default", "memory-packed", "memory-reference3",
        "file-header-checked", "file-header-trusted", "file-structure-checked", "file-structure-trusted",
    ];

    /// <summary>The reference engines (not counted per member).</summary>
    internal static readonly string[] ReferenceNames = ["ukkonen", "v6"];

    /// <summary>The longest text of the Tiny and the Small bands for which the container is built with the tiny parameters of the builder (they make the buckets, runs and windows of the build many); above it the default tuning with a small difference cover.</summary>
    internal const int TinyBuilderUpTo = 300;

    private static readonly int[] Periods = [31, 57, 133];

    private readonly List<FuzzSubject> _subjects = [];

    private FuzzSubjects(UkkonenTree ukkonen, PersistentSuffixTree? v6, bool fileSkipped)
    {
        Ukkonen = ukkonen;
        V6 = v6;
        FileSkipped = fileSkipped;
    }

    /// <summary>The engines under test.</summary>
    internal IReadOnlyList<FuzzSubject> Subjects => _subjects;

    /// <summary>The Ukkonen tree of the text (oracle L2).</summary>
    internal UkkonenTree Ukkonen { get; }

    /// <summary>The v6 tree of the text (oracle L3), null when the text is longer than the harness builds it for.</summary>
    internal PersistentSuffixTree? V6 { get; }

    /// <summary>True when the text has a code unit above 255 and the file-backed engines were left out.</summary>
    internal bool FileSkipped { get; }

    /// <summary>Ticks spent building the Ukkonen tree.</summary>
    internal long ReferenceBuildTicks { get; private set; }

    /// <summary>Ticks spent building the v6 tree.</summary>
    internal long V6BuildTicks { get; private set; }

    /// <summary>The options of the container builder for a text of <paramref name="n"/> symbols.</summary>
    internal static SesaBuildOptions ContainerOptions(int n, FullFuzzConfig cfg) => n <= TinyBuilderUpTo
        ? SesaRun.Tiny(cfg.Variant, 1, cfg.Psi) with { MeasureRss = false }
        : new SesaBuildOptions
        {
            Threads = 1,
            Fsync = false,
            MeasureRss = false,
            IncludePsi = cfg.Psi,
            Period = Periods[cfg.Variant % Periods.Length],
            ExactPeriod = true,
            Seed = (ulong)(cfg.Variant + 1),
        };

    /// <summary>Builds the engines of <paramref name="text"/>; the container lives in <paramref name="workDir"/> only while the trees are opened.</summary>
    internal static FuzzSubjects Build(string text, FullFuzzConfig cfg, MemberCounters counters, string workDir, int v6UpTo)
    {
        int n = text.Length;
        long t0 = Stopwatch.GetTimestamp();
        UkkonenTree ukkonen = UkkonenTree.Build(text);
        long t1 = Stopwatch.GetTimestamp();
        PersistentSuffixTree? v6 = n <= v6UpTo ? PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text)) : null;
        var set = new FuzzSubjects(ukkonen, v6, text.Any(c => c > 255)) { ReferenceBuildTicks = t1 - t0, V6BuildTicks = Stopwatch.GetTimestamp() - t1 };
        try
        {
            FuzzSteps.Run(FuzzSteps.Open, () => set.BuildMemory(text, cfg, counters));
            if (!set.FileSkipped) set.BuildFiles(text, cfg, counters, workDir);
            return set;
        }
        catch
        {
            set.Dispose();
            throw;
        }
    }

    private void BuildMemory(string text, FullFuzzConfig cfg, MemberCounters counters)
    {
        int n = text.Length;
        ChildLayoutKind?[] layouts = n <= 64
            ? [null, ChildLayoutKind.Packed, ChildLayoutKind.Reference3]
            : n <= 2000 ? [null, cfg.Variant % 2 == 0 ? ChildLayoutKind.Packed : ChildLayoutKind.Reference3] : [null];
        foreach (ChildLayoutKind? layout in layouts)
        {
            string engine = layout switch { null => "memory-default", ChildLayoutKind.Packed => "memory-packed", _ => "memory-reference3" };
            int e = counters.EngineIndex(engine);
            long t = Stopwatch.GetTimestamp();
            EsaSuffixTree tree;
            if (layout is null)
            {
                switch (cfg.Variant % 4)
                {
                    case 0:
                        counters.Hit(e, "EsaSuffixTreeFactory.Create(String)");
                        tree = (EsaSuffixTree)EsaSuffixTreeFactory.Create(text);
                        break;
                    case 1:
                        counters.Hit(e, "EsaSuffixTreeFactory.CreateEsa(String)");
                        tree = EsaSuffixTreeFactory.CreateEsa(text);
                        break;
                    case 2:
                        counters.Hit(e, "EsaSuffixTreeFactory.Create(ITextSource)");
                        tree = (EsaSuffixTree)EsaSuffixTreeFactory.Create(new StringTextSource(text));
                        break;
                    default:
                        counters.Hit(e, "EsaSuffixTreeFactory.CreateEsa(ITextSource)");
                        tree = EsaSuffixTreeFactory.CreateEsa(new StringTextSource(text));
                        break;
                }
            }
            else
            {
                tree = EsaSuffixTree.Build(text, layout);
            }
            _subjects.Add(new FuzzSubject(engine, new CountingTree(tree, counters, e), false, Stopwatch.GetTimestamp() - t));
        }
    }

    private void BuildFiles(string text, FullFuzzConfig cfg, MemberCounters counters, string workDir)
    {
        string input = Path.Combine(workDir, "in.txt");
        string target = Path.Combine(workDir, "c.sesa");
        File.WriteAllBytes(input, Encoding.Latin1.GetBytes(text));
        long t = Stopwatch.GetTimestamp();
        try
        {
            SesaBuilder.Build(input, target, ContainerOptions(text.Length, cfg));
            long built = Stopwatch.GetTimestamp() - t;
            (string Engine, EsaValidationLevel Level, EsaAccessMode Access)[] combos =
            [
                ("file-header-checked", EsaValidationLevel.Header, EsaAccessMode.Checked),
                ("file-header-trusted", EsaValidationLevel.Header, EsaAccessMode.Trusted),
                ("file-structure-checked", EsaValidationLevel.Structure, EsaAccessMode.Checked),
                ("file-structure-trusted", EsaValidationLevel.Structure, EsaAccessMode.Trusted),
            ];
            foreach (var (engine, level, access) in combos)
            {
                int e = counters.EngineIndex(engine);
                long o = Stopwatch.GetTimestamp();
                EsaSuffixTree tree = FuzzSteps.Run(FuzzSteps.Open, () => Open(target, level, access, cfg, counters, e));
                _subjects.Add(new FuzzSubject(engine, new CountingTree(tree, counters, e), true, built + Stopwatch.GetTimestamp() - o));
            }
        }
        finally
        {
            File.Delete(input);
            File.Delete(target);
        }
    }

    private static EsaSuffixTree Open(string path, EsaValidationLevel level, EsaAccessMode access, FullFuzzConfig cfg, MemberCounters counters, int engine)
    {
        counters.Hit(engine, "EsaOpenOptions.Validation");
        counters.Hit(engine, "EsaOpenOptions.Access");
        bool sha = level == EsaValidationLevel.Structure && cfg.Variant % 2 == 0;
        if (sha) counters.Hit(engine, "EsaOpenOptions.VerifyTextSha256");
        var options = new EsaOpenOptions { Validation = level, Access = access, VerifyTextSha256 = sha };
        if (level == EsaValidationLevel.Header && access == EsaAccessMode.Checked && cfg.Variant % 3 == 0)
        {
            counters.Hit(engine, "EsaSuffixTreeFactory.Open(String,EsaOpenOptions)");
            return EsaSuffixTreeFactory.Open(path, options);
        }
        counters.Hit(engine, "SuffixIndexFactory.Open(String,EsaOpenOptions)");
        return (EsaSuffixTree)SuffixIndexFactory.Open(path, options);
    }

    /// <summary>Disposes every tree.</summary>
    public void Dispose()
    {
        FuzzSteps.Run(FuzzSteps.Open, () =>
        {
            foreach (FuzzSubject subject in _subjects) subject.Tree.Dispose();
        });
        _subjects.Clear();
        V6?.Dispose();
    }
}
