namespace Seqeron.Genomics.Alignment;

/// <summary>Options of <see cref="MummerDeltaFilter.Run"/> (<c>delta-filter</c> of MUMmer 3.23).</summary>
public sealed record MummerDeltaFilterOptions
{
    /// <summary><c>-q</c>: map each position of each query to its best hit in the reference (allowing reference overlaps).</summary>
    public bool QueryLis { get; init; }

    /// <summary><c>-r</c>: map each position of each reference to its best hit in the query.</summary>
    public bool ReferenceLis { get; init; }

    /// <summary><c>-g</c>: 1-to-1 global alignment not allowing rearrangements.</summary>
    public bool GlobalLis { get; init; }

    /// <summary><c>-1</c>: 1-to-1 alignment allowing rearrangements (intersection of -r and -q).</summary>
    public bool OneToOne { get; init; }

    /// <summary><c>-m</c>: many-to-many alignment allowing rearrangements (union of -r and -q).</summary>
    public bool ManyToMany { get; init; }

    /// <summary><c>-l</c>: minimum alignment length (default 0).</summary>
    public long MinLength { get; init; }

    /// <summary><c>-i</c>: minimum alignment identity, 0 - 100 (default 0).</summary>
    public float MinIdentity { get; init; }

    /// <summary><c>-u</c>: minimum alignment uniqueness, 0 - 100 (default 0).</summary>
    public float MinUnique { get; init; }

    /// <summary><c>-o</c>: maximum alignment overlap for -r and -q as a percent of the alignment length, 0 - 100 (default 100).</summary>
    public float MaxOverlap { get; init; } = 100f;

    /// <summary><c>-e</c>: negligible alignment score (a negative value, the default -1, picks randomly among equally good chains with the program's <c>rand()</c>).</summary>
    public float Epsilon { get; init; } = -1f;
}

/// <summary>
/// <c>delta-filter</c> of MUMmer 3.23 (<c>src/tigr/delta-filter.cc</c>): reads a delta file and keeps the alignments that survive, in the order <c>-i -l -u -q -r -g -m -1</c>; the output is a delta file again
/// (grouped by query record). The choice among equally good chains uses glibc <c>rand()</c> seeded with 1, as the program does (<c>srand (1)</c>), so the output is reproducible.
/// </summary>
public static class MummerDeltaFilter
{
    /// <summary>Filters the delta file text.</summary>
    /// <param name="deltaText">The delta file.</param>
    /// <param name="options">The switches (null = none: the file is rewritten grouped by query).</param>
    /// <exception cref="ArgumentOutOfRangeException">A percentage outside [0, 100] or a negative length (the program's range errors).</exception>
    /// <exception cref="FormatException">The delta file cannot be parsed.</exception>
    public static string Run(string deltaText, MummerDeltaFilterOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(deltaText);
        var o = options ?? new MummerDeltaFilterOptions();
        if (o.MinIdentity < 0.0 || o.MinIdentity > 100.0)
            throw new ArgumentOutOfRangeException(nameof(options), "ERROR: Minimum identity must be within the range [0, 100]");
        if (o.MinLength < 0)
            throw new ArgumentOutOfRangeException(nameof(options), "ERROR: Minimum length must be greater than or equal to zero");
        if (o.MinUnique < 0.0 || o.MinUnique > 100.0)
            throw new ArgumentOutOfRangeException(nameof(options), "ERROR: Minimum uniqueness must be within the range [0, 100]");
        if (o.MaxOverlap < 0.0 || o.MaxOverlap > 100.0)
            throw new ArgumentOutOfRangeException(nameof(options), "ERROR: Maximum overlap must be within the range [0, 100]");

        var graph = MummerDeltaGraph.Build(deltaText, true);
        if (o.MinIdentity > 0 || o.MinLength > 0)
            graph.FlagScore(o.MinLength, o.MinIdentity);
        if (o.MinUnique > 0)
            graph.FlagUniq(o.MinUnique);
        if (o.QueryLis)
            graph.FlagQLis(o.Epsilon, o.MaxOverlap);
        if (o.ReferenceLis)
            graph.FlagRLis(o.Epsilon, o.MaxOverlap);
        if (o.GlobalLis)
            graph.FlagGlis(o.Epsilon);
        if (o.ManyToMany)
            graph.FlagMtoM(o.Epsilon, o.MaxOverlap);
        if (o.OneToOne)
            graph.Flag1To1(o.Epsilon, o.MaxOverlap);
        return graph.OutputDelta();
    }
}
