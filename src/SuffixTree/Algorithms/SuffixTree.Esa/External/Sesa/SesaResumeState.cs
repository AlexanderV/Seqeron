using System.Buffers;
using System.Globalization;
using System.Text;
using System.Text.Json;

namespace SuffixTree.Esa.External;

/// <summary>The names of the stages of a resumable build (work package 1b.9e, ADR-004 section 4.3 "Crash safety and atomicity"), in the order in which they complete.</summary>
internal static class SesaStages
{
    /// <summary>PH0: the container is created at its reserved length and TEXT, ALPHABET and COUNTS are written.</summary>
    internal const string Packed = "packed";

    /// <summary>PH1: the rank file of the difference-cover sample (scratch). Superseded by <see cref="Sa"/> (the rank file is deleted when PH3 is done).</summary>
    internal const string Ranks = "ranks";

    /// <summary>PH2 / PH3: the SA section.</summary>
    internal const string Sa = "sa";

    /// <summary>PH4 / PH5: the container grown behind its reserved length, the LCP, the overflow list and CHILD sections.</summary>
    internal const string Lcp = "lcp";

    /// <summary>PH6a: the BWT file (scratch).</summary>
    internal const string Bwt = "bwt";

    /// <summary>PH6b: the PSI section; the BWT file is deleted.</summary>
    internal const string Psi = "psi";

    /// <summary>All stage names in the order of completion.</summary>
    internal static readonly string[] Order = [Packed, Ranks, Sa, Lcp, Bwt, Psi];

    /// <summary>The position of <paramref name="name"/> in <see cref="Order"/>, or -1.</summary>
    internal static int IndexOf(string name) => Array.IndexOf(Order, name);
}

/// <summary>One file range whose bytes a completed stage vouches for: a range of the container (<see cref="InContainer"/>) or a scratch file, its length and its CRC-32C.</summary>
/// <param name="Role">What it is: a section name (<c>text</c>, <c>alphabet</c>, <c>counts</c>, <c>sa</c>, <c>child</c>, <c>psi</c>, <c>lcp</c>, <c>overflow</c>) or the path of a scratch file below the scratch directory (<c>ranks.bin</c>, <c>ph6/bwt.bin</c>).</param>
/// <param name="InContainer">True for a range of the container, false for a scratch file.</param>
/// <param name="Offset">Offset in the container (0 for a scratch file).</param>
/// <param name="Length">Bytes.</param>
/// <param name="Crc">CRC-32C of those bytes.</param>
internal sealed record SesaFileCheck(string Role, bool InContainer, long Offset, long Length, uint Crc);

/// <summary>A completed stage: what it wrote (to be verified on resume) and the few numbers the later stages and the header need that only this stage knows.</summary>
/// <param name="Name">One of <see cref="SesaStages"/>.</param>
/// <param name="Checks">The outputs.</param>
/// <param name="Values">Named numbers (for <see cref="SesaStages.Lcp"/>: the width, the escapes, the final length, the statistics of the scan; for <see cref="SesaStages.Bwt"/>: the primary row; for <see cref="SesaStages.Ranks"/>: the period and a digest of the cover).</param>
internal sealed record SesaStageRecord(string Name, IReadOnlyList<SesaFileCheck> Checks, IReadOnlyDictionary<string, long> Values)
{
    internal long Value(string key) => Values.TryGetValue(key, out long v) ? v : throw new InvalidDataException($"the record of the stage {Name} has no value {key}");
}

/// <summary>
/// The contents of <c>state.json</c> (ADR-004 section 4.3): what a resume has to know to trust the files of an earlier run. The <em>key</em> is what changes a byte of the output or of a skipped stage's files: the writer, the
/// input (length, alphabet, counts and SHA-256 from PH0), whether PSI is written, the codec of CHILD and, through them, the layout (<see cref="ReservedLength"/>). Thread count, bucket sizes, windows, the seed and the like
/// are not in it: the output does not depend on them (ADR-004 A8, owner answer Q5), so a build with other values of them resumes from the same state. The one parameter that does decide the bytes of an intermediate file,
/// the difference-cover period, is a value of the stage that depends on it (<see cref="SesaStages.Ranks"/>).
/// </summary>
internal sealed record SesaResumeState
{
    /// <summary>The format of this file.</summary>
    internal const int CurrentVersion = 1;

    internal int Version { get; init; } = CurrentVersion;

    /// <summary>The writer whose files these are (<see cref="SesaResumeFormat.WriterId"/>); another writer is not trusted.</summary>
    internal string Writer { get; init; } = SesaResumeFormat.WriterId;

    internal long N { get; init; }

    internal int Bits { get; init; }

    internal int Sigma { get; init; }

    /// <summary>Lower-case hex SHA-256 of the input.</summary>
    internal string Sha256 { get; init; } = string.Empty;

    /// <summary>The symbols of the input in increasing byte order and their counts (part of the fingerprint: they decide the packed text).</summary>
    internal IReadOnlyList<int> Letters { get; init; } = [];

    internal IReadOnlyList<long> Counts { get; init; } = [];

    internal bool Psi { get; init; }

    /// <summary><c>tagged</c> or <c>packed</c>.</summary>
    internal string ChildCodec { get; init; } = "tagged";

    /// <summary>The length at which the container was created (everything up to the end of the narrowest LCP).</summary>
    internal long ReservedLength { get; init; }

    /// <summary>The completed stages in order.</summary>
    internal IReadOnlyList<SesaStageRecord> Stages { get; init; } = [];

    internal SesaStageRecord? Find(string name) => Stages.FirstOrDefault(s => s.Name == name);

    /// <summary>The state with one more completed stage; a stage that this one supersedes (<see cref="SesaStages.Sa"/> replaces <see cref="SesaStages.Ranks"/>, <see cref="SesaStages.Psi"/> replaces <see cref="SesaStages.Bwt"/>) is dropped.</summary>
    internal SesaResumeState With(SesaStageRecord stage)
    {
        var next = Stages.Where(s => s.Name != stage.Name && !(s.Name == SesaStages.Ranks && stage.Name == SesaStages.Sa) && !(s.Name == SesaStages.Bwt && stage.Name == SesaStages.Psi)).ToList();
        next.Add(stage);
        next.Sort((a, b) => SesaStages.IndexOf(a.Name).CompareTo(SesaStages.IndexOf(b.Name)));
        return this with { Stages = next };
    }

    /// <summary>The state with only the first <paramref name="count"/> stages.</summary>
    internal SesaResumeState Truncate(int count) => this with { Stages = [.. Stages.Take(count)] };
}

/// <summary>
/// The file format of <c>state.json</c>: one JSON object (UTF-8, no byte order mark) whose last member is <c>"crc32c"</c>, the CRC-32C of every byte before the comma that precedes it, eight lower-case hex digits; a newline
/// follows. The file stays valid JSON (a person can read it) and a damaged, cut or edited file fails the CRC. Parsing never throws: anything wrong is <c>false</c>, and the caller starts clean.
/// </summary>
internal static class SesaResumeFormat
{
    /// <summary>The writer of the files whose state this is; it must change whenever a byte that the builder writes changes (the same string as in BUILD_INFO).</summary>
    internal const string WriterId = "SuffixTree.Esa.External.SesaBuilder 1";

    /// <summary>The name of the state file in the scratch directory.</summary>
    internal const string FileName = "state.json";

    /// <summary>The name that the next state is written under before it replaces <see cref="FileName"/>.</summary>
    internal const string TempName = "state.json.tmp";

    /// <summary>The longest state file that is read (a state with six stages is about 2 KiB).</summary>
    internal const int MaxBytes = 1 << 20;

    private static readonly byte[] s_marker = Encoding.ASCII.GetBytes(",\"crc32c\":\"");

    /// <summary>The bytes of the state file.</summary>
    internal static byte[] Serialize(SesaResumeState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        var buffer = new ArrayBufferWriter<byte>(1024);
        using (var w = new Utf8JsonWriter(buffer))
        {
            w.WriteStartObject();
            w.WriteString("format", "sesa-resume");
            w.WriteNumber("version", state.Version);
            w.WriteString("writer", state.Writer);
            w.WriteNumber("n", state.N);
            w.WriteNumber("bits", state.Bits);
            w.WriteNumber("sigma", state.Sigma);
            w.WriteString("sha256", state.Sha256);
            w.WriteStartArray("letters");
            foreach (int l in state.Letters) w.WriteNumberValue(l);
            w.WriteEndArray();
            w.WriteStartArray("counts");
            foreach (long c in state.Counts) w.WriteNumberValue(c);
            w.WriteEndArray();
            w.WriteBoolean("psi", state.Psi);
            w.WriteString("child_codec", state.ChildCodec);
            w.WriteNumber("reserved_length", state.ReservedLength);
            w.WriteStartArray("stages");
            foreach (SesaStageRecord s in state.Stages)
            {
                w.WriteStartObject();
                w.WriteString("name", s.Name);
                w.WriteStartArray("files");
                foreach (SesaFileCheck f in s.Checks)
                {
                    w.WriteStartObject();
                    w.WriteString("role", f.Role);
                    w.WriteBoolean("container", f.InContainer);
                    w.WriteNumber("offset", f.Offset);
                    w.WriteNumber("length", f.Length);
                    w.WriteString("crc32c", Hex(f.Crc));
                    w.WriteEndObject();
                }
                w.WriteEndArray();
                w.WriteStartObject("values");
                foreach (KeyValuePair<string, long> v in s.Values.OrderBy(v => v.Key, StringComparer.Ordinal)) w.WriteNumber(v.Key, v.Value);
                w.WriteEndObject();
                w.WriteEndObject();
            }
            w.WriteEndArray();
            w.WriteEndObject();
        }
        return Seal(buffer.WrittenSpan[..^1]);
    }

    /// <summary>The file for a JSON object whose closing brace has been left off: the CRC member, the brace and a newline are appended.</summary>
    internal static byte[] Seal(ReadOnlySpan<byte> bodyWithoutClosingBrace)
    {
        uint crc = Crc32C.Compute(bodyWithoutClosingBrace);
        var file = new byte[bodyWithoutClosingBrace.Length + s_marker.Length + 8 + 3];
        bodyWithoutClosingBrace.CopyTo(file);
        s_marker.CopyTo(file, bodyWithoutClosingBrace.Length);
        Encoding.ASCII.GetBytes(Hex(crc), file.AsSpan(bodyWithoutClosingBrace.Length + s_marker.Length));
        int at = bodyWithoutClosingBrace.Length + s_marker.Length + 8;
        file[at] = (byte)'"';
        file[at + 1] = (byte)'}';
        file[at + 2] = (byte)'\n';
        return file;
    }

    /// <summary>True when <paramref name="file"/> ends with a CRC member that matches the bytes before it.</summary>
    internal static bool CrcMatches(ReadOnlySpan<byte> file)
    {
        int tail = s_marker.Length + 8 + 2;
        if (file.Length > 0 && file[^1] == (byte)'\n') file = file[..^1];
        if (file.Length < tail + 1) return false;
        int at = file.Length - tail;
        if (!file.Slice(at, s_marker.Length).SequenceEqual(s_marker)) return false;
        if (file[^2] != (byte)'"' || file[^1] != (byte)'}') return false;
        ReadOnlySpan<byte> digits = file.Slice(at + s_marker.Length, 8);
        foreach (byte d in digits)
        {
            if (!(d is >= (byte)'0' and <= (byte)'9' or >= (byte)'a' and <= (byte)'f')) return false; // exactly what Hex writes: lower case
        }
        uint stored = uint.Parse(Encoding.ASCII.GetString(digits), NumberStyles.AllowHexSpecifier, CultureInfo.InvariantCulture);
        return stored == Crc32C.Compute(file[..at]);
    }

    /// <summary>Parses a state file. False (and a null state) for anything that is not exactly what <see cref="Serialize"/> writes: a cut or edited file, another version, wrong types, stages out of order or repeated, negative lengths.</summary>
    internal static bool TryParse(ReadOnlySpan<byte> file, out SesaResumeState? state)
    {
        state = null;
        if (file.Length == 0 || file.Length > MaxBytes || !CrcMatches(file)) return false;
        try
        {
            using JsonDocument doc = JsonDocument.Parse(file.ToArray());
            JsonElement root = doc.RootElement;
            if (root.ValueKind != JsonValueKind.Object) return false;
            if (root.GetProperty("format").GetString() != "sesa-resume") return false;
            if (root.GetProperty("version").GetInt32() != SesaResumeState.CurrentVersion) return false;
            var letters = new List<int>();
            foreach (JsonElement e in root.GetProperty("letters").EnumerateArray()) letters.Add(e.GetInt32());
            var counts = new List<long>();
            foreach (JsonElement e in root.GetProperty("counts").EnumerateArray()) counts.Add(e.GetInt64());
            var stages = new List<SesaStageRecord>();
            foreach (JsonElement s in root.GetProperty("stages").EnumerateArray())
            {
                var checks = new List<SesaFileCheck>();
                foreach (JsonElement f in s.GetProperty("files").EnumerateArray())
                {
                    if (!uint.TryParse(f.GetProperty("crc32c").GetString(), NumberStyles.AllowHexSpecifier, CultureInfo.InvariantCulture, out uint crc)) return false;
                    checks.Add(new SesaFileCheck(f.GetProperty("role").GetString() ?? string.Empty, f.GetProperty("container").GetBoolean(), f.GetProperty("offset").GetInt64(), f.GetProperty("length").GetInt64(), crc));
                }
                var values = new Dictionary<string, long>(StringComparer.Ordinal);
                foreach (JsonProperty v in s.GetProperty("values").EnumerateObject()) values[v.Name] = v.Value.GetInt64();
                stages.Add(new SesaStageRecord(s.GetProperty("name").GetString() ?? string.Empty, checks, values));
            }
            var parsed = new SesaResumeState
            {
                Version = SesaResumeState.CurrentVersion,
                Writer = root.GetProperty("writer").GetString() ?? string.Empty,
                N = root.GetProperty("n").GetInt64(),
                Bits = root.GetProperty("bits").GetInt32(),
                Sigma = root.GetProperty("sigma").GetInt32(),
                Sha256 = root.GetProperty("sha256").GetString() ?? string.Empty,
                Letters = letters,
                Counts = counts,
                Psi = root.GetProperty("psi").GetBoolean(),
                ChildCodec = root.GetProperty("child_codec").GetString() ?? string.Empty,
                ReservedLength = root.GetProperty("reserved_length").GetInt64(),
                Stages = stages,
            };
            if (!IsWellFormed(parsed)) return false;
            state = parsed;
            return true;
        }
        catch (Exception e) when (e is JsonException or InvalidOperationException or FormatException or KeyNotFoundException or ArgumentException or OverflowException)
        {
            return false;
        }
    }

    /// <summary>What no honest writer produces: negative numbers, a fingerprint of the wrong shape, unknown, repeated or unordered stages, a stage that cannot exist next to another.</summary>
    private static bool IsWellFormed(SesaResumeState s)
    {
        if (s.N < 0 || s.Bits < 0 || s.Sigma < 0 || s.ReservedLength < 0) return false;
        if (s.Sha256.Length != 64 || s.Sha256.Any(c => !Uri.IsHexDigit(c) || char.IsUpper(c))) return false;
        if (s.Letters.Count != s.Sigma || s.Counts.Count != s.Sigma || s.Counts.Any(c => c < 0) || s.Letters.Any(l => l is < 0 or > 255)) return false;
        if (s.ChildCodec is not ("tagged" or "packed")) return false;
        int previous = -1;
        foreach (SesaStageRecord st in s.Stages)
        {
            int at = SesaStages.IndexOf(st.Name);
            if (at <= previous) return false;
            previous = at;
            foreach (SesaFileCheck f in st.Checks)
            {
                if (f.Offset < 0 || f.Length < 0 || f.Role.Length == 0 || (!f.InContainer && (f.Offset != 0 || f.Role.Contains("..", StringComparison.Ordinal) || Path.IsPathRooted(f.Role)))) return false;
            }
        }
        bool ranks = s.Find(SesaStages.Ranks) is not null, sa = s.Find(SesaStages.Sa) is not null;
        if (ranks && sa) return false;
        if (s.Find(SesaStages.Packed) is null && s.Stages.Count > 0) return false;
        if (!sa && (s.Find(SesaStages.Lcp) is not null || s.Find(SesaStages.Bwt) is not null || s.Find(SesaStages.Psi) is not null)) return false;
        if (s.Find(SesaStages.Lcp) is null && (s.Find(SesaStages.Bwt) is not null || s.Find(SesaStages.Psi) is not null)) return false;
        if (!s.Psi && (s.Find(SesaStages.Bwt) is not null || s.Find(SesaStages.Psi) is not null)) return false;
        return true;
    }

    private static string Hex(uint crc) => crc.ToString("x8", CultureInfo.InvariantCulture);
}
