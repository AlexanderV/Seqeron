using System.Reflection;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// The operation coverage list of the full differential fuzz (work package 1b.8), derived from the types and not typed in: every public instance member of the audited interfaces, and the public entry points of the factories and of
/// <see cref="EsaOpenOptions"/> that open and describe a file-backed tree. A key is <c>Type.Member(ParameterTypes)</c> (a property is <c>Type.Member</c>), so an overload that is added later is a new key. The committed table
/// (<c>docs/refactoring/suffix-tree-upgrade/analysis/s1b/harness_member_coverage.tsv</c>) must have one row per key (a reflection test fails otherwise), and the counting decorator <see cref="CountingTree"/> reports every call under its key.
/// </summary>
internal static class FuzzMemberKeys
{
    /// <summary>The interfaces whose members the harness calls on the tree under test: the four search and analysis interfaces, the contract interfaces of ADR-001, what a tree returns (<see cref="ITextSource"/>, <see cref="IIndexText"/>), the visitor that <c>Traverse</c> calls, and <see cref="IDisposable"/>.</summary>
    internal static readonly Type[] Interfaces =
    [
        typeof(ISuffixTreeSearch), typeof(ISuffixTreeAnalysis), typeof(ISuffixTreeDiagnostics), typeof(ISuffixTree),
        typeof(ISuffixIndex), typeof(ISuffixTreeNavigation), typeof(ISuffixArrayAccess), typeof(IIndexText), typeof(ITextSource),
        typeof(ISuffixTreeVisitor), typeof(IDisposable),
    ];

    /// <summary>The classes whose public static methods or properties are the way to a tree or its options (the file-backed entry points).</summary>
    internal static readonly Type[] ApiTypes = [typeof(SuffixIndexFactory), typeof(EsaSuffixTreeFactory), typeof(EsaOpenOptions)];

    private static readonly Lazy<(string[] Keys, Dictionary<string, int> Index)> Table = new(Build);

    /// <summary>Every key, in a stable order.</summary>
    internal static IReadOnlyList<string> All => Table.Value.Keys;

    /// <summary>The position of <paramref name="key"/> in <see cref="All"/>; an unknown key is a defect of the caller (a typo in the decorator).</summary>
    internal static int IndexOf(string key) =>
        Table.Value.Index.TryGetValue(key, out int i) ? i : throw new ArgumentException($"'{key}' is not a member key of the audited interfaces", nameof(key));

    /// <summary>The key of a member.</summary>
    internal static string KeyOf(MemberInfo member)
    {
        string owner = TypeName(member.DeclaringType!);
        if (member is PropertyInfo p)
        {
            ParameterInfo[] indexParameters = p.GetIndexParameters();
            return indexParameters.Length == 0 ? $"{owner}.{p.Name}" : $"{owner}.{p.Name}[{string.Join(",", indexParameters.Select(x => TypeName(x.ParameterType)))}]";
        }
        var m = (MethodInfo)member;
        return $"{owner}.{m.Name}({string.Join(",", m.GetParameters().Select(x => TypeName(x.ParameterType)))})";
    }

    private static string TypeName(Type t)
    {
        if (t.IsByRef) return TypeName(t.GetElementType()!);
        if (!t.IsGenericType) return t.Name;
        string name = t.Name[..t.Name.IndexOf('`', StringComparison.Ordinal)];
        return $"{name}<{string.Join(",", t.GetGenericArguments().Select(TypeName))}>";
    }

    private static (string[], Dictionary<string, int>) Build()
    {
        const BindingFlags instance = BindingFlags.Public | BindingFlags.Instance | BindingFlags.DeclaredOnly;
        const BindingFlags statics = BindingFlags.Public | BindingFlags.Static | BindingFlags.DeclaredOnly;
        var keys = new SortedSet<string>(StringComparer.Ordinal);
        foreach (Type i in Interfaces)
            foreach (MemberInfo member in i.GetMembers(instance))
                if (member is PropertyInfo || (member is MethodInfo m && !m.IsSpecialName)) keys.Add(KeyOf(member));
        foreach (Type t in ApiTypes)
        {
            foreach (MemberInfo member in t.GetMembers(statics))
                if (member is MethodInfo m && !m.IsSpecialName) keys.Add(KeyOf(member));
            foreach (PropertyInfo p in t.GetProperties(instance)) keys.Add(KeyOf(p));
        }
        string[] all = [.. keys];
        return (all, all.Select((k, i) => (k, i)).ToDictionary(x => x.k, x => x.i, StringComparer.Ordinal));
    }
}

/// <summary>
/// The harness steps that call members (the classes of the coverage table): a worker names the step it is in, and every counted call records it, so the committed table says for each member which steps call it.
/// </summary>
internal static class FuzzSteps
{
    internal const string Definitions = "FacadeChecks.CheckEverything";
    internal const string Search = "FacadeChecks.CheckSearch";
    internal const string Digest = "FuzzDigest.Compute";
    internal const string Contract = "EsaContractChecker.CheckAll";
    internal const string ContractDigest = "FuzzContractDigest.Compute";
    internal const string Arrays = "FullFuzzCheck.SuffixArrayLevel";
    internal const string Concurrency = "FullFuzzConcurrency.Script";
    internal const string Corruption = "CorruptionFuzz.Battery";
    internal const string Open = "FuzzSubjects.Build";
    internal const string Other = "other";

    /// <summary>Every step name, in the order of the table columns.</summary>
    internal static readonly string[] All = [Definitions, Search, Digest, Contract, ContractDigest, Arrays, Concurrency, Corruption, Open, Other];

    [ThreadStatic]
    private static string? _current;

    [ThreadStatic]
    private static string? _failed;

    /// <summary>The step of the calling thread.</summary>
    internal static string Current => _current ?? Other;

    /// <summary>Runs <paramref name="body"/> as <paramref name="step"/>.</summary>
    internal static void Run(string step, Action body)
    {
        string? previous = _current;
        _current = step;
        try
        {
            body();
        }
        catch (Exception)
        {
            _failed ??= step;
            throw;
        }
        finally
        {
            _current = previous;
        }
    }

    /// <summary>Runs <paramref name="body"/> as <paramref name="step"/> and returns its value.</summary>
    internal static T Run<T>(string step, Func<T> body)
    {
        string? previous = _current;
        _current = step;
        try
        {
            return body();
        }
        catch (Exception)
        {
            _failed ??= step;
            throw;
        }
        finally
        {
            _current = previous;
        }
    }

    /// <summary>The innermost step that an exception has passed through since the last call of this method, or null (the step is part of the operation name of a failure).</summary>
    internal static string? TakeFailedStep()
    {
        string? step = _failed;
        _failed = null;
        return step;
    }

    internal static int IndexOf(string step) => Array.IndexOf(All, step);
}

/// <summary>The calls of a run, per engine and per member key, and which steps made them. Thread safe; a count is an <see cref="Interlocked"/> increment.</summary>
internal sealed class MemberCounters(IReadOnlyList<string> engines)
{
    private readonly long[][] _calls = [.. engines.Select(_ => new long[FuzzMemberKeys.All.Count])];
    private readonly bool[,] _steps = new bool[FuzzMemberKeys.All.Count, FuzzSteps.All.Length];

    /// <summary>The engine names, in the order of <see cref="Calls"/>.</summary>
    internal IReadOnlyList<string> Engines { get; } = engines;

    /// <summary>Counts one call of the member <paramref name="key"/> on the engine <paramref name="engine"/>.</summary>
    internal void Hit(int engine, string key)
    {
        int k = FuzzMemberKeys.IndexOf(key);
        Interlocked.Increment(ref _calls[engine][k]);
        int s = FuzzSteps.IndexOf(FuzzSteps.Current);
        _steps[k, s < 0 ? FuzzSteps.All.Length - 1 : s] = true;
    }

    /// <summary>The number of calls of a member on an engine.</summary>
    internal long Calls(int engine, string key) => Interlocked.Read(ref _calls[engine][FuzzMemberKeys.IndexOf(key)]);

    /// <summary>The number of calls of a member on every engine.</summary>
    internal long Calls(string key) => Enumerable.Range(0, Engines.Count).Sum(e => Calls(e, key));

    /// <summary>The number of calls of every member on an engine.</summary>
    internal long Total(int engine) => _calls[engine].Sum(c => Interlocked.Read(ref c));

    /// <summary>The names of the steps that called a member.</summary>
    internal IReadOnlyList<string> StepsOf(string key)
    {
        int k = FuzzMemberKeys.IndexOf(key);
        return [.. FuzzSteps.All.Where((_, s) => _steps[k, s])];
    }

    /// <summary>The engine index of a name; -1 when it is not an engine of this run.</summary>
    internal int EngineIndex(string name)
    {
        for (int i = 0; i < Engines.Count; i++)
            if (Engines[i] == name) return i;
        return -1;
    }
}
