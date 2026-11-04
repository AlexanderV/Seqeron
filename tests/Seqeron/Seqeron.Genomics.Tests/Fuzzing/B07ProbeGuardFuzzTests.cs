using System.Reflection;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Heavy-tier guard fuzz for review batch B07 F67–F71 on the probe side: every public static method of
/// <see cref="ProbeDesigner"/> is invoked with random degenerate arguments (null / empty / junk / over-long strings,
/// negative / zero / int.MinValue / int.MaxValue / long extremes, NaN / ±∞ / negative / zero / huge doubles, undefined
/// enum values, default / illegal parameter records, null and null-element collections, degenerate scoring schemes,
/// base-frequency vectors). A call may succeed or throw an <see cref="ArgumentException"/> (incl.
/// <see cref="ArgumentNullException"/> / <see cref="ArgumentOutOfRangeException"/>) raised by the library itself whose
/// <see cref="ArgumentException.ParamName"/>, when set, names one of the called method's parameters (F71: each public
/// method reports its own parameter name). Any other exception (NullReference, IndexOutOfRange, Overflow,
/// InvalidOperation, KeyNotFound, …), an <see cref="ArgumentException"/> thrown from inside the framework (the "raw
/// framework AOORE" of F68) or a foreign ParamName is a defect. The per-run seed comes from NUnit's randomizer and is
/// printed on failure. Lazy results are enumerated (up to 50 items).
///
/// Test Units: PROBE-DESIGN-001, PROBE-VALID-001, PROBE-EVALUE-001.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("MolTools")]
public class B07ProbeGuardFuzzTests
{
    private const int CallsPerMethod = 60;

    private static IEnumerable<TestCaseData> Methods() =>
        typeof(ProbeDesigner)
            .GetMethods(BindingFlags.Public | BindingFlags.Static | BindingFlags.DeclaredOnly)
            .Where(m => !m.IsGenericMethodDefinition && !m.IsSpecialName)
            .Select(m => new TestCaseData(m).SetArgDisplayNames(
                $"ProbeDesigner.{m.Name}({string.Join(",", m.GetParameters().Select(p => p.ParameterType.Name))})"));

    private static readonly string Long70 = string.Concat(Enumerable.Repeat("ACGTTGCAAGGCCTTA", 5))[..70];

    private static readonly string?[] Strings =
    [
        null, "", " ", "A", "AC", "ACGT", "acgtn", "NNNN", "ACGU", "XYZ!", "GGGGCCCC", "ACGTACGTACGTACGTACGTAC",
        "ATGCATGCATGCATGCATGCAT", Long70, "ACGT.ACGT", "RYKMSWBDHV", "AUGCAUGCAUGCAUGCAUGC",
    ];

    private static readonly int[] Ints = [int.MinValue, -5, -1, 0, 1, 2, 3, 4, 5, 18, 25, 36, 60, 100, int.MaxValue];

    private static readonly long[] Longs = [long.MinValue, -1, 0, 1, 100, 1_000_000, long.MaxValue];

    private static readonly double[] Doubles =
    [
        double.NaN, double.NegativeInfinity, double.PositiveInfinity, -1, 0, 1e-12, 0.05, 0.25, 0.5, 0.75, 1, 37, 50,
        100, 1e9, double.MaxValue,
    ];

    [TestCaseSource(nameof(Methods))]
    public void PublicStaticMethod_DegenerateArguments_OnlyOwnArgumentExceptionsFromTheLibrary(MethodInfo method)
    {
        int seed = TestContext.CurrentContext.Random.Next();
        var rng = new Random(seed);
        var parameters = method.GetParameters();
        var failures = new List<string>();
        for (int call = 0; call < CallsPerMethod && failures.Count < 3; call++)
        {
            object?[] args;
            try
            {
                args = parameters.Select(p => Sample(p.ParameterType, rng)).ToArray();
            }
            catch (ArgumentException)
            {
                continue; // an argument object itself rejected the sample
            }

            var task = Task.Run(() => Invoke(method, args));
            if (!task.Wait(TimeSpan.FromSeconds(60)))
            {
                failures.Add($"timeout with ({Describe(args)})");
                break;
            }
            if (task.Result is { } error && !IsOwnLibraryArgumentException(error, parameters))
                failures.Add($"{error.GetType().Name} '{error.Message}' (param {(error as ArgumentException)?.ParamName}, " +
                             $"site {error.TargetSite?.DeclaringType?.FullName}.{error.TargetSite?.Name}) with ({Describe(args)})");
        }
        Assert.That(failures, Is.Empty, $"seed {seed}: ProbeDesigner.{method.Name}\n" + string.Join("\n", failures));
    }

    private static Exception? Invoke(MethodInfo method, object?[] args)
    {
        try
        {
            object? result = method.Invoke(null, args);
            if (result is System.Collections.IEnumerable seq and not string)
                foreach (var _ in seq.Cast<object?>().Take(50)) { }
            return null;
        }
        catch (TargetInvocationException e) when (e.InnerException is not null)
        {
            return e.InnerException;
        }
        catch (Exception e)
        {
            return e;
        }
    }

    private static bool IsOwnLibraryArgumentException(Exception e, ParameterInfo[] parameters)
    {
        if (e is not ArgumentException ae)
            return false;
        string site = e.TargetSite?.DeclaringType?.FullName ?? "";
        bool frameworkSite = site.StartsWith("System.", StringComparison.Ordinal)
                             && !site.StartsWith("System.Argument", StringComparison.Ordinal);
        bool ownParam = ae.ParamName is null || parameters.Any(p => p.Name == ae.ParamName);
        return !frameworkSite && ownParam;
    }

    private static object? Sample(Type t, Random rng)
    {
        var underlying = Nullable.GetUnderlyingType(t);
        if (underlying is not null)
            return rng.Next(3) == 0 ? null : Sample(underlying, rng);
        if (t == typeof(string)) return Strings[rng.Next(Strings.Length)];
        if (t == typeof(int)) return Ints[rng.Next(Ints.Length)];
        if (t == typeof(long)) return Longs[rng.Next(Longs.Length)];
        if (t == typeof(double)) return Doubles[rng.Next(Doubles.Length)];
        if (t == typeof(bool)) return rng.Next(2) == 0;
        if (t.IsEnum)
        {
            var values = Enum.GetValues(t);
            return rng.Next(5) == 0 ? Enum.ToObject(t, 999) : values.GetValue(rng.Next(values.Length));
        }
        if (t == typeof(ProbeDesigner.ProbeParameters))
            return SampleParameters(rng);
        if (t == typeof(ProbeDesigner.Primer3ProbeSettings))
            return rng.Next(5) switch
            {
                0 => null,
                1 => new ProbeDesigner.Primer3ProbeSettings(),
                2 => new ProbeDesigner.Primer3ProbeSettings(MinSize: 0, OptSize: 0, MaxSize: 0),
                3 => new ProbeDesigner.Primer3ProbeSettings(MinSize: 30, OptSize: 20, MaxSize: 10, MinTm: double.NaN),
                _ => new ProbeDesigner.Primer3ProbeSettings(MonovalentMillimolar: -1, DnaConcentrationNanomolar: 0),
            };
        if (t == typeof(global::SuffixTree.ISuffixTree))
            return rng.Next(3) switch
            {
                0 => null,
                1 => global::SuffixTree.SuffixTree.Build(Long70),
                _ => global::SuffixTree.SuffixTree.Build("A"),
            };
        if (t == typeof(IEnumerable<string>))
            return rng.Next(5) switch
            {
                0 => null,
                1 => Array.Empty<string>(),
                2 => new[] { Long70 },
                3 => new[] { "ACGT", null!, "" },
                _ => new[] { Long70 + Long70, "acgtnacgtn", "XYZ" },
            };
        if (t == typeof(IReadOnlyList<int>))
            return rng.Next(3) switch { 0 => null, 1 => Array.Empty<int>(), _ => Enumerable.Repeat(rng.Next(-5, 120), 70).ToArray() };
        if (t == typeof(IReadOnlyList<double>))
            return rng.Next(6) switch
            {
                0 => null,
                1 => Array.Empty<double>(),
                2 => new[] { 0.25, 0.25, 0.25, 0.25 },
                3 => new[] { 1.0, 0, 0, 0 },
                4 => new[] { double.NaN, 1, 1, 1 },
                _ => new[] { -1.0, 1, 1, 1 },
            };
        if (t == typeof(ScoringMatrix))
            return rng.Next(7) switch
            {
                0 => null,
                1 => SequenceAligner.BlastDna,
                2 => new ScoringMatrix(1, -3, -5, -2),
                3 => new ScoringMatrix(0, 0, 0, 0),
                4 => new ScoringMatrix(5, -1, -5, -2),
                5 => new ScoringMatrix(2, -3, 5, 2),
                _ => new ScoringMatrix(int.MaxValue, int.MinValue, int.MinValue, int.MinValue),
            };
        if (t.IsValueType)
            return Activator.CreateInstance(t);
        return null;
    }

    private static ProbeDesigner.ProbeParameters SampleParameters(Random rng) =>
        rng.Next(9) switch
        {
            0 => default,
            1 => ProbeDesigner.Defaults.qPCR,
            2 => ProbeDesigner.Defaults.qPCR with { Ranking = ProbeDesigner.ProbeRanking.Primer3Penalty, OptLength = 99, OptTm = double.NaN },
            3 => ProbeDesigner.Defaults.qPCR with { Ranking = (ProbeDesigner.ProbeRanking)999 },
            4 => ProbeDesigner.Defaults.qPCR with { ThermodynamicScreenMaxLength = 59 },
            5 => ProbeDesigner.Defaults.qPCR with { MonovalentMillimolar = -1, DnaConcentrationNanomolar = 0, MaxNearestNeighborLength = -1 },
            6 => ProbeDesigner.Defaults.qPCR with { StructureScreen = ProbeDesigner.ProbeStructureScreen.Heuristic, MinLength = -5, MaxLength = 3 },
            7 => ProbeDesigner.Defaults.qPCR with { MinLength = 4, MaxLength = 8, MinGc = double.NaN, MaxTm = double.NegativeInfinity },
            _ => ProbeDesigner.Defaults.Microarray with { MinLength = 50, MaxLength = 51 },
        };

    private static string Describe(object?[] args) =>
        string.Join(", ", args.Select(a => a switch
        {
            null => "null",
            string s => s.Length > 24 ? $"\"{s[..20]}…\"({s.Length})" : $"\"{s}\"",
            double d => d.ToString("R", System.Globalization.CultureInfo.InvariantCulture),
            IEnumerable<string> e => $"[{string.Join(",", e.Select(x => x is null ? "null" : x.Length > 12 ? x[..10] + "…" : x))}]",
            _ => a.ToString() ?? "?",
        }));
}
