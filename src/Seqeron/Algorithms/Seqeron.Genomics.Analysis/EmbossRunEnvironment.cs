namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The clock a run of an EMBOSS program sees (B09 WP62, F281): the date the writers print (<c>ajTimeToday</c>, <c>time()</c>) and the microseconds of <c>gettimeofday</c> that <c>ajRandomSeed</c> turns into the seed of the random
/// number generator (<c>seed = tv_usec % 9999 + 1</c>; with the EMBOSS variable <c>timetoday</c> defined the seed is 0). A run with a fixed environment is deterministic.
/// </summary>
/// <param name="Now">The date and time of the run (UTC); null = the clock of the machine.</param>
/// <param name="Microseconds">The microseconds of <c>gettimeofday</c> at the seeding of the random number generator; null = the clock of the machine.</param>
/// <param name="TimeToday">The value of the EMBOSS variable <c>timetoday</c> (<c>SET timetoday "2012-03-04 05:06:07"</c> in <c>.embossrc</c>), null when it is not defined: the seed is 0 whatever the clock says and the date
/// of every time object is the one of the variable (<c>ajTimeSetS</c>, read as <c>yyyy-mm-dd hh:mm:ss</c> by <c>ajFmtScanS</c>; a text that gives no number leaves the clock).</param>
public sealed record EmbossRunEnvironment(DateTime? Now = null, int? Microseconds = null, string? TimeToday = null)
{
    private static readonly AsyncLocal<EmbossRunEnvironment?> CurrentValue = new();

    /// <summary>The data directory of the run (<c>EMBOSS_DATA</c>: where <c>ajDatafileNewInNameS</c> looks after the work directory and where <c>ajDatafileNewOutNameS</c> writes; a file by its path in the files of the run); null = the installation of EMBOSS 6.6.0 (<see cref="DefaultDataDirectory"/>).</summary>
    public string? DataDirectory { get; init; }

    /// <summary>The data directory of the Debian / Ubuntu package <c>emboss-data</c> 6.6.0.</summary>
    public const string DefaultDataDirectory = "/usr/share/EMBOSS/data/";

    /// <summary>The data directory with the trailing slash (<c>ajDirnameFix</c>).</summary>
    internal string DataDir
    {
        get
        {
            if (string.IsNullOrEmpty(DataDirectory))
                return DefaultDataDirectory;
            return DataDirectory.EndsWith('/') ? DataDirectory : DataDirectory + "/";
        }
    }

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
        if (TimeToday is not null)
            return 0;                                                       // "always zero microseconds in a defined time"
        long usec = Microseconds ?? (DateTime.UtcNow.Ticks / 10) % 1_000_000;
        return (int)(usec % 9999) + 1;
    }

    /// <summary>The date of the run: <c>ajTimeSetLocal</c> (the variable <c>timetoday</c> when it reads as a time, else the clock); null = the clock of the machine.</summary>
    internal DateTime? RunDate => TimeToday is not null && ParseTimeToday(TimeToday) is { } t ? t : Now;

    // ajTimeSetS: "%4d-%2d-%2d %2d:%2d:%2d" (the fields that the text does not give are 0), the year is years since 1900 (a year above 1899 is reduced), mktime normalises the fields; no field at all: not a time
    private static DateTime? ParseTimeToday(string text)
    {
        var field = new int[6];
        int[] width = [4, 2, 2, 2, 2, 2];
        char[] literal = ['-', '-', ' ', ':', ':', '\0'];
        int pos = 0, n = 0;
        for (int i = 0; i < 6; i++)
        {
            while (pos < text.Length && char.IsWhiteSpace(text[pos]))
                pos++;
            int start = pos, digits = 0;
            if (pos < text.Length && text[pos] is '+' or '-' && width[i] > 1)
                pos++;
            while (pos < text.Length && pos - start < width[i] && char.IsAsciiDigit(text[pos]))
            {
                pos++;
                digits++;
            }

            if (digits == 0)
                break;
            field[i] = int.Parse(text.AsSpan(start, pos - start), System.Globalization.NumberStyles.AllowLeadingSign, System.Globalization.CultureInfo.InvariantCulture);
            n++;
            if (i == 5)
                break;
            if (literal[i] == ' ')
            {
                while (pos < text.Length && char.IsWhiteSpace(text[pos]))
                    pos++;
            }
            else if (pos < text.Length && text[pos] == literal[i])
            {
                pos++;
            }
            else
            {
                break;
            }
        }

        if (n == 0)
            return null;
        int year = field[0] > 1899 ? field[0] - 1900 : field[0];
        try
        {
            return new DateTime(1900, 1, 1, 0, 0, 0, DateTimeKind.Unspecified).AddYears(year).AddMonths(field[1] - 1).AddDays(field[2] - 1).AddHours(field[3]).AddMinutes(field[4]).AddSeconds(field[5]);
        }
        catch (ArgumentOutOfRangeException)
        {
            return null;
        }
    }

    private sealed class Scope(EmbossRunEnvironment? previous) : IDisposable
    {
        public void Dispose() => CurrentValue.Value = previous;
    }
}
