using System;

namespace SuffixTree;

/// <summary>
/// Hard limits that follow from the 32-bit counters used by every suffix-tree backend.
/// </summary>
public static class SuffixTreeLimits
{
    /// <summary>
    /// The largest text length, in characters, that a suffix tree accepts.
    /// </summary>
    /// <remarks>
    /// A suffix tree over <c>n</c> characters plus the terminator has <c>n + 1</c> leaves and at most
    /// <c>n</c> internal nodes including the root, i.e. at most <c>2n + 1</c> nodes, and node counts and
    /// node indexes are <see cref="int"/> values. <c>2n + 1 &lt;= int.MaxValue</c> holds exactly for
    /// <c>n &lt;= int.MaxValue / 2 = 1,073,741,823</c>. A <see cref="string"/> cannot get anywhere near this
    /// length (about 1,073,741,791 characters), so only a custom <see cref="ITextSource"/> can exceed it.
    /// The bound is the worst case; a text below it may still need gigabytes (about 170 bytes per character
    /// in memory, see <c>SuffixTree.Build</c>).
    /// </remarks>
    public const int MaxTextLength = int.MaxValue / 2;

    /// <summary>
    /// Throws <see cref="ArgumentException"/> when <paramref name="length"/> is larger than <see cref="MaxTextLength"/>.
    /// </summary>
    /// <param name="length">The reported length of the text.</param>
    /// <param name="paramName">The name of the argument that carries the text.</param>
    /// <exception cref="ArgumentException">The text is longer than <see cref="MaxTextLength"/> characters.</exception>
    public static void ThrowIfTextTooLong(int length, string paramName)
    {
        if (length > MaxTextLength)
            throw new ArgumentException(TooLongMessage(length), paramName);
    }

    /// <summary>
    /// Returns <see langword="true"/> when a text of <paramref name="length"/> characters is within <see cref="MaxTextLength"/>.
    /// </summary>
    /// <param name="length">The reported length of the text.</param>
    /// <returns><see langword="true"/> if the text can be indexed; otherwise <see langword="false"/>.</returns>
    public static bool IsTextLengthSupported(int length) => length <= MaxTextLength;

    private static string TooLongMessage(int length) =>
        "The text has " + length + " characters; a suffix tree supports at most " + MaxTextLength +
        " characters because its node count (up to 2n + 1) must fit in a 32-bit integer.";
}
