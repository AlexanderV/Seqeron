using System;

namespace SuffixTree.Indexing;

/// <summary>
/// A place in the virtual suffix tree of the indexed text: an explicit node or a point inside an edge. A locus is a small value; engines hand
/// it out and take it back and are the only parties that interpret <see cref="Lo"/>, <see cref="Hi"/> and <see cref="Aux"/>. Callers use
/// <see cref="Length"/> (the number of matched symbols) and compare loci only for equality.
/// </summary>
public readonly struct IndexLocus : IEquatable<IndexLocus>
{
    /// <summary>Creates a locus from its four engine-defined parts.</summary>
    /// <param name="lo">First rank of the locus interval (suffix-array engines), otherwise engine-defined.</param>
    /// <param name="hi">Last rank of the locus interval (inclusive), otherwise engine-defined.</param>
    /// <param name="length">String depth: the number of matched symbols.</param>
    /// <param name="aux">Engine-defined (enhanced suffix array: the depth of the lowest explicit node at or below the locus).</param>
    public IndexLocus(long lo, long hi, long length, long aux)
    {
        Lo = lo;
        Hi = hi;
        Length = length;
        Aux = aux;
    }

    /// <summary>Engine-defined; for suffix-array engines the first rank of the locus interval.</summary>
    public long Lo { get; }

    /// <summary>Engine-defined; for suffix-array engines the last rank of the locus interval (inclusive).</summary>
    public long Hi { get; }

    /// <summary>String depth: the number of matched symbols (0 for the root).</summary>
    public long Length { get; }

    /// <summary>Engine-defined; for the enhanced suffix array the depth of the lowest explicit node at or below the locus.</summary>
    public long Aux { get; }

    /// <summary>Component-wise equality.</summary>
    public bool Equals(IndexLocus other) => Lo == other.Lo && Hi == other.Hi && Length == other.Length && Aux == other.Aux;

    /// <inheritdoc />
    public override bool Equals(object? obj) => obj is IndexLocus other && Equals(other);

    /// <inheritdoc />
    public override int GetHashCode() => HashCode.Combine(Lo, Hi, Length, Aux);

    /// <summary>Component-wise equality.</summary>
    public static bool operator ==(IndexLocus left, IndexLocus right) => left.Equals(right);

    /// <summary>Component-wise inequality.</summary>
    public static bool operator !=(IndexLocus left, IndexLocus right) => !left.Equals(right);

    /// <inheritdoc />
    public override string ToString() => $"[{Lo}..{Hi}] depth {Length} aux {Aux}";
}
