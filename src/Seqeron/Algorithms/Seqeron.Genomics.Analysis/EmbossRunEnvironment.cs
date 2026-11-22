namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The clock a run of an EMBOSS program sees (B09 WP62, F281): the date the writers print (<c>ajTimeToday</c>, <c>time()</c>) and the microseconds of <c>gettimeofday</c> that <c>ajRandomSeed</c> turns into the seed of the random
/// number generator (<c>seed = tv_usec % 9999 + 1</c>; with the EMBOSS variable <c>timetoday</c> defined the seed is 0). A run with a fixed environment is deterministic.
/// </summary>
/// <param name="Now">The date and time of the run (UTC); null = the clock of the machine.</param>
/// <param name="Microseconds">The microseconds of <c>gettimeofday</c> at the seeding of the random number generator; null = the clock of the machine.</param>
/// <param name="TimeToday">The EMBOSS variable <c>timetoday</c> is defined (<c>SET timetoday 1</c> in <c>.embossrc</c>): the seed is 0 whatever the clock says.</param>
public sealed record EmbossRunEnvironment(DateTime? Now = null, int? Microseconds = null, bool TimeToday = false)
{
    private static readonly AsyncLocal<EmbossRunEnvironment?> CurrentValue = new();

    /// <summary>The environment of the run in progress (the default environment outside a run).</summary>
    internal static EmbossRunEnvironment Current => CurrentValue.Value ?? new EmbossRunEnvironment();

    internal static IDisposable Enter(EmbossRunEnvironment environment)
    {
        var previous = CurrentValue.Value;
        CurrentValue.Value = environment;
        return new Scope(previous);
    }

    /// <summary>The seed of <c>ajRandomSeed</c> for this environment.</summary>
    internal int Seed()
    {
        if (TimeToday)
            return 0;                                                       // "always zero microseconds in a defined time"
        long usec = Microseconds ?? (DateTime.UtcNow.Ticks / 10) % 1_000_000;
        return (int)(usec % 9999) + 1;
    }

    private sealed class Scope(EmbossRunEnvironment? previous) : IDisposable
    {
        public void Dispose() => CurrentValue.Value = previous;
    }
}
