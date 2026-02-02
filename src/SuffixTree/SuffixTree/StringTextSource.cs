using System;
using System.Collections;
using System.Collections.Generic;

namespace SuffixTree;

/// <summary>
/// A string-based implementation of <see cref="ITextSource"/>.
/// </summary>
public sealed class StringTextSource : ITextSource
{
    private readonly string _value;

    /// <summary>
    /// Initializes a new instance of the <see cref="StringTextSource"/> class.
    /// </summary>
    /// <param name="value">The string value to wrap. Cannot be null.</param>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="value"/> is null.</exception>
    public StringTextSource(string value)
    {
        _value = value ?? throw new ArgumentNullException(nameof(value));
    }

    /// <inheritdoc/>
    public int Length => _value.Length;

    /// <inheritdoc/>
    public char this[int index] => _value[index];

    /// <inheritdoc/>
    public string Substring(int start, int length) => _value.Substring(start, length);

    /// <inheritdoc/>
    public ReadOnlySpan<char> Slice(int start, int length) => _value.AsSpan(start, length);

    /// <inheritdoc/>
    public IEnumerator<char> GetEnumerator() => _value.GetEnumerator();

    /// <inheritdoc/>
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    /// <summary>
    /// Returns the underlying string.
    /// </summary>
    public override string ToString() => _value;

    /// <summary>
    /// Implicitly converts a string to a <see cref="StringTextSource"/>.
    /// </summary>
    public static implicit operator StringTextSource(string value) => new(value);
}
