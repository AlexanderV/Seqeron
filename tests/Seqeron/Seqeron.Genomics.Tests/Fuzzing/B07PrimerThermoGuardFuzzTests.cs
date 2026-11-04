using System.Reflection;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Heavy-tier guard fuzz for review batch B07 F68 / F68a / F38: every public static method of
/// <see cref="PrimerDesigner"/> and <see cref="ThermoConstants"/> is invoked with random degenerate arguments
/// (null / empty / junk / over-long strings, negative / zero / int.MinValue / int.MaxValue integers, NaN / ±∞ /
/// negative / zero / huge doubles, undefined enum values, default structs). A call may succeed or throw an
/// <see cref="ArgumentException"/> (incl. <see cref="ArgumentNullException"/> / <see cref="ArgumentOutOfRangeException"/>)
/// raised by the library itself; any other exception (NullReference, IndexOutOfRange, Overflow, InvalidOperation,
/// KeyNotFound, OutOfMemory, …) or an <see cref="ArgumentException"/> thrown from inside the framework
/// (String.Substring, collections, Math.Clamp, … — the "raw framework AOORE" of F68) is a defect.
/// The per-run seed comes from NUnit's randomizer and is printed on failure.
///
/// Test Units: PRIMER-STRUCT-001, PRIMER-DESIGN-001, PRIMER-TM-001, PRIMER-NNTM-001, PRIMER-HAIRPIN-001, PRIMER-DIMER-001.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("MolTools")]
public class B07PrimerThermoGuardFuzzTests
{
    private const int CallsPerMethod = 60;

    // File-system loaders are not argument guards (a junk path is an I/O error by design).
    private static readonly HashSet<string> Excluded = ["FromGenomeTester4Files"];

    private static IEnumerable<TestCaseData> Methods() =>
        new[] { typeof(PrimerDesigner), typeof(ThermoConstants) }
            .SelectMany(t => t.GetMethods(BindingFlags.Public | BindingFlags.Static | BindingFlags.DeclaredOnly))
            .Where(m => !m.IsGenericMethodDefinition && !m.IsSpecialName && !Excluded.Contains(m.Name))
            .Select(m => new TestCaseData(m).SetArgDisplayNames(
                $"{m.DeclaringType!.Name}.{m.Name}({string.Join(",", m.GetParameters().Select(p => p.ParameterType.Name))})"));

    private static readonly string Long70 = string.Concat(Enumerable.Repeat("ACGTTGCAAGGCCTTA", 5))[..70];

    private static readonly string?[] Strings =
    [
        null, "", " ", "A", "AC", "ACGT", "acgtn", "NNNN", "ACGU", "XYZ!", "GGGGCCCC", "ACGTACGTACGTACGTACGTAC",
        "ATGCATGCATGCATGCATGCAT", Long70, "ACGT.ACGT", "RYKMSWBDHV",
    ];

    private static readonly int[] Ints = [int.MinValue, -5, -1, 0, 1, 2, 3, 4, 5, 18, 25, 36, 60, 100, int.MaxValue];

    private static readonly double[] Doubles =
    [
        double.NaN, double.NegativeInfinity, double.PositiveInfinity, -1, 0, 1e-12, 0.05, 1, 37, 50, 100, 1e9,
        double.MaxValue,
    ];

    [TestCaseSource(nameof(Methods))]
    public void PublicStaticMethod_DegenerateArguments_OnlyArgumentExceptionsFromTheLibrary(MethodInfo method)
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
                args = parameters.Select(p => Sample(p.ParameterType.IsByRef ? p.ParameterType.GetElementType()! : p.ParameterType, rng)).ToArray();
            }
            catch (ArgumentException)
            {
                continue; // an argument object itself rejected the sample (e.g. a DnaSequence of junk)
            }

            var task = Task.Run(() => Invoke(method, args));
            if (!task.Wait(TimeSpan.FromSeconds(30)))
            {
                failures.Add($"timeout with ({Describe(args)})");
                break;
            }
            if (task.Result is { } error && !IsLibraryArgumentException(error, parameters))
                failures.Add($"{error.GetType().Name} '{error.Message}' (param {(error as ArgumentException)?.ParamName}, " +
                             $"site {error.TargetSite?.DeclaringType?.FullName}.{error.TargetSite?.Name}) with ({Describe(args)})");
        }
        Assert.That(failures, Is.Empty, $"seed {seed}: {method.DeclaringType!.Name}.{method.Name}\n" + string.Join("\n", failures));
    }

    private static Exception? Invoke(MethodInfo method, object?[] args)
    {
        try
        {
            object? result = method.Invoke(null, args);
            if (result is System.Collections.IEnumerable seq and not string)
                foreach (var _ in seq.Cast<object?>().Take(5000)) { }
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

    private static bool IsLibraryArgumentException(Exception e, ParameterInfo[] parameters)
    {
        if (e is not ArgumentException ae)
            return false;
        // Thrown by a framework member other than the Argument*Exception throw helpers ⇒ an unguarded raw failure.
        string site = e.TargetSite?.DeclaringType?.FullName ?? "";
        bool frameworkSite = site.StartsWith("System.", StringComparison.Ordinal)
                             && !site.StartsWith("System.Argument", StringComparison.Ordinal);
        bool foreignParam = ae.ParamName is "startIndex" or "index" or "count" or "length" or "capacity" or "destination"
                            && parameters.All(p => p.Name != ae.ParamName);
        return !frameworkSite && !foreignParam;
    }

    private static object? Sample(Type t, Random rng)
    {
        var underlying = Nullable.GetUnderlyingType(t);
        if (underlying is not null)
            return rng.Next(3) == 0 ? null : Sample(underlying, rng);
        if (t == typeof(string)) return Strings[rng.Next(Strings.Length)];
        if (t == typeof(int)) return Ints[rng.Next(Ints.Length)];
        if (t == typeof(double)) return Doubles[rng.Next(Doubles.Length)];
        if (t == typeof(bool)) return rng.Next(2) == 0;
        if (t == typeof(char)) return "ACGTN.x"[rng.Next(7)];
        if (t.IsEnum)
        {
            var values = Enum.GetValues(t);
            return rng.Next(5) == 0 ? Enum.ToObject(t, 999) : values.GetValue(rng.Next(values.Length));
        }
        if (t == typeof(DnaSequence))
            return rng.Next(4) switch
            {
                0 => null,
                1 => new DnaSequence("ACGTACGTAC"),
                2 => new DnaSequence(Long70 + Long70 + Long70),
                _ => new DnaSequence(Long70),
            };
        if (t == typeof(PrimerParameters))
            return rng.Next(3) switch
            {
                0 => default(PrimerParameters),
                1 => PrimerDesigner.DefaultParameters,
                _ => PrimerDesigner.Primer3DefaultParameters,
            };
        if (t == typeof(PrimerPairOptions))
            return rng.Next(2) == 0 ? null : PrimerPairOptions.Primer3Defaults with { ProductSizeRanges = [new ProductSizeRange(40, 200)] };
        if (t == typeof(IReadOnlyList<int>))
            return rng.Next(3) switch { 0 => null, 1 => Array.Empty<int>(), _ => Enumerable.Repeat(rng.Next(-5, 120), 70).ToArray() };
        if (t == typeof(PrimerMisprimingLibrary))
            return rng.Next(2) == 0 ? null : new PrimerMisprimingLibrary([new KeyValuePair<string, string>("rep*2", Long70)]);
        if (t == typeof(PrimerMaskingKmerLists))
            return rng.Next(2) == 0
                ? null
                : new PrimerMaskingKmerLists(new Dictionary<string, int> { [Long70[..11]] = 500 },
                    new Dictionary<string, int> { [Long70[..16]] = 500 });
        if (t.IsValueType)
            return Activator.CreateInstance(t);
        return null;
    }

    private static string Describe(object?[] args) =>
        string.Join(", ", args.Select(a => a switch
        {
            null => "null",
            string s => s.Length > 24 ? $"\"{s[..20]}…\"({s.Length})" : $"\"{s}\"",
            DnaSequence d => $"Dna({d.Length})",
            double d => d.ToString("R", System.Globalization.CultureInfo.InvariantCulture),
            _ => a.ToString() ?? "?",
        }));
}
