using System.Globalization;

namespace SuffixTree.EmProto.Harness;

/// <summary>Minimal <c>--key value</c> parser of the command line.</summary>
internal sealed class Options
{
    private readonly Dictionary<string, string> _values = new(StringComparer.Ordinal);

    internal static Options Parse(string[] args)
    {
        var o = new Options();
        int i = 0;
        while (i < args.Length)
        {
            if (!args[i].StartsWith("--", StringComparison.Ordinal)) throw new ArgumentException($"unexpected argument '{args[i]}'");
            string key = args[i][2..];
            bool hasValue = i + 1 < args.Length && !args[i + 1].StartsWith("--", StringComparison.Ordinal);
            o._values[key] = hasValue ? args[i + 1] : "true";
            i += hasValue ? 2 : 1;
        }
        return o;
    }

    internal string? Get(string key) => _values.TryGetValue(key, out string? v) ? v : null;

    internal string Require(string key) => Get(key) ?? throw new ArgumentException($"missing --{key}");

    internal string Str(string key, string fallback) => Get(key) ?? fallback;

    internal long Long(string key, long fallback) => Get(key) is { } v ? long.Parse(v, CultureInfo.InvariantCulture) : fallback;

    internal int Int(string key, int fallback) => (int)Long(key, fallback);

    internal double Double(string key, double fallback) => Get(key) is { } v ? double.Parse(v, CultureInfo.InvariantCulture) : fallback;

    internal bool Flag(string key) => Get(key) is { } v && v != "false";
}
