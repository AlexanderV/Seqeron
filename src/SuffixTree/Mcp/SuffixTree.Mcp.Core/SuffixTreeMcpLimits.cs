using System.Globalization;
using Microsoft.Extensions.Configuration;

namespace SuffixTree.Mcp.Core;

/// <summary>
/// Size limit of the suffix-tree MCP tools (ST-S0 WP 0.6).
/// </summary>
/// <remarks>
/// <para>
/// Every <c>suffix_tree_*</c> tool (and the two genomics tools that use a tree) rebuilds an in-memory suffix tree
/// from its JSON string arguments on every call. Measured on 10 Mbp of random DNA, that tree takes
/// <b>169 bytes of managed heap per character</b> (190 for a 20-letter alphabet, 193 for a highly repetitive text;
/// <c>docs/refactoring/suffix-tree-upgrade/BENCHMARK.md</c>), plus about 15 s of build time per 10 Mbp.
/// </para>
/// <para>
/// <b>How the default was chosen.</b> A single call should stay near 1 GiB of heap, so that a server on a small
/// machine, or several concurrent calls, do not exhaust memory. Taking the worst measured 193 bytes per character
/// and rounding up to a design constant of 200 bytes per character, 1 GiB allows
/// 1 073 741 824 / 200 = 5.37 million characters; the default is rounded down to
/// <see cref="DefaultMaxTextChars"/> = 5 000 000 characters (about 0.85 GiB for random DNA, 0.97 GiB for the
/// worst measured text, 7 s of build time). The limit counts the total number of characters of all text
/// arguments of one call (texts, queries and patterns; for the k-strings tool every text plus one separator).
/// </para>
/// <para>
/// <b>Configuration.</b> The limit can be raised or lowered when the server starts with the configuration key
/// <c>SuffixTree:MaxTextChars</c> (command line <c>--SuffixTree:MaxTextChars=20000000</c>, environment variable
/// <c>SuffixTree__MaxTextChars=20000000</c>, or <c>appsettings.json</c>). A value that is not a positive integer stops
/// the server at startup. Sessions that keep one index across calls are planned for S8; until then the limit is the
/// protection.
/// </para>
/// </remarks>
public static class SuffixTreeMcpLimits
{
    /// <summary>Design constant: bytes of heap per input character, the worst measured 193 rounded up.</summary>
    public const int DesignBytesPerChar = 200;

    /// <summary>Default per-call limit: 5 million characters (about 1 GiB at <see cref="DesignBytesPerChar"/>).</summary>
    public const int DefaultMaxTextChars = 5_000_000;

    /// <summary>Configuration key that overrides the limit.</summary>
    public const string ConfigurationKey = "SuffixTree:MaxTextChars";

    private static int s_maxTextChars = DefaultMaxTextChars;

    /// <summary>The current per-call limit on the total number of characters of all text arguments.</summary>
    public static int MaxTextChars => Volatile.Read(ref s_maxTextChars);

    /// <summary>
    /// Reads <see cref="ConfigurationKey"/> from <paramref name="configuration"/>; keeps the default when it is absent.
    /// </summary>
    /// <exception cref="InvalidOperationException">The value is not a positive 32-bit integer.</exception>
    public static void Configure(IConfiguration configuration)
    {
        ArgumentNullException.ThrowIfNull(configuration);
        string? raw = configuration[ConfigurationKey];
        if (string.IsNullOrWhiteSpace(raw))
        {
            SetMaxTextChars(DefaultMaxTextChars);
            return;
        }
        if (!int.TryParse(raw, NumberStyles.None, CultureInfo.InvariantCulture, out int value) || value < 1)
            throw new InvalidOperationException(
                $"Configuration value '{ConfigurationKey}' must be a positive integer (characters per call), but was '{raw}'.");
        SetMaxTextChars(value);
    }

    internal static void SetMaxTextChars(int value) => Volatile.Write(ref s_maxTextChars, value);

    /// <summary>
    /// Throws <see cref="ArgumentException"/> when the call's text arguments exceed the limit; called before any tree is built.
    /// </summary>
    internal static void EnsureWithinLimit(long totalChars, string paramName)
    {
        int max = MaxTextChars;
        if (totalChars <= max)
            return;

        throw new ArgumentException(
            $"Input too large: {totalChars.ToString("N0", CultureInfo.InvariantCulture)} characters in total exceed the limit of " +
            $"{max.ToString("N0", CultureInfo.InvariantCulture)} characters per call. Each call builds a suffix tree in memory " +
            $"(about {DesignBytesPerChar} bytes per character). Split the input, or raise the limit when starting the server with " +
            $"'--{ConfigurationKey}=<n>' or the environment variable '{ConfigurationKey.Replace(":", "__", StringComparison.Ordinal)}=<n>'.",
            paramName);
    }

    /// <summary>Convenience overload for one or more text arguments.</summary>
    internal static void EnsureWithinLimit(string paramName, params string[] texts)
    {
        long total = 0;
        foreach (string text in texts)
            total += text.Length;
        EnsureWithinLimit(total, paramName);
    }
}
