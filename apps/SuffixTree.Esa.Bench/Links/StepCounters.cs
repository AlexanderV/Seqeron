using System.Runtime.CompilerServices;

namespace SuffixTree.Esa.Bench.Links;

/// <summary>
/// Step accounting of the link kernels (prototype P1, work package 1b.4a). The kernels are generic over the counter, so the timed instantiation
/// (<see cref="NoCount"/>) carries no accounting code at all, the counting one (<see cref="CountSteps"/>) is used for the step statistics and the budgeted one
/// (<see cref="BudgetSteps"/>) stops a kernel that is not bounded (a^n, Fibonacci) instead of letting it run for minutes.
/// </summary>
internal interface IStepCounter
{
    /// <summary>Adds <paramref name="steps"/> elementary steps (one array entry read of the search structure each).</summary>
    void Add(int steps);

    /// <summary>True when the budget is used up; kernels check it once per loop iteration and give up.</summary>
    bool Exceeded { get; }

    /// <summary>Steps added so far (0 for <see cref="NoCount"/>).</summary>
    long Steps { get; }
}

/// <summary>No accounting; every member compiles away.</summary>
internal struct NoCount : IStepCounter
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public readonly void Add(int steps)
    {
    }

    public readonly bool Exceeded
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => false;
    }

    public readonly long Steps => 0;
}

/// <summary>Counts steps, no limit.</summary>
internal struct CountSteps : IStepCounter
{
    private long _steps;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Add(int steps) => _steps += steps;

    public readonly bool Exceeded
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => false;
    }

    public readonly long Steps => _steps;
}

/// <summary>Counts steps and reports <see cref="Exceeded"/> once more than <c>budget</c> steps were added.</summary>
internal struct BudgetSteps(long budget) : IStepCounter
{
    private long _steps;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Add(int steps) => _steps += steps;

    public readonly bool Exceeded
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _steps > budget;
    }

    public readonly long Steps => _steps;
}
