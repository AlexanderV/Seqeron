using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit A2-29b, B09 F310): the ksw2 extension functions <see cref="Ksw2Aligner.ExtendZ"/>, <see cref="Ksw2Aligner.ExtendD"/>, <see cref="Ksw2Aligner.ExtendZ2"/>, <see cref="Ksw2Aligner.ExtendD2"/> against the REAL
/// programs: <c>ksw_extz</c> / <c>ksw_extd</c> of lh3/ksw2 (scalar) and <c>ksw_extz2_sse</c> / <c>ksw_extd2_sse</c> of minimap2 v2.26 (SSE4.1; the SSE2-only build printed the same lines), run on the user's machine
/// (Docker package docs/Validation/review-2026-09/handoff/A2-29b-ksw2, run 2: 2 000 cases, 500 per variant, band, Z-drop, end bonus, generic matrices, approximate maximum / drop, score only, left / right gaps, reversed CIGAR;
/// TestData/Alignment/ksw2_real_oracle.tsv.gz, scripts/alignment/convert_ksw2_returned.py). Every field of <c>ksw_extz_t</c> (max, zdropped, max_q, max_t, mqe, mqe_t, mte, mte_q, score, n_cigar, reach_end) and the CIGAR text are compared.
/// The cases on which the C program itself has undefined behaviour (<see cref="Ksw2ReferenceUndefined"/>: ASan / valgrind: heap buffer overflow at ksw2_extz.c:116 and the same line of ksw2_extd.c, 33 aborted runs) are refused, and
/// the refused set is derived from the condition, not listed.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class Ksw2Aligner_RealProgram_Tests
{
    private const string ResourceName = "Seqeron.Genomics.Tests.TestData.Alignment.ksw2_real_oracle.tsv";

    private sealed record Row(int Id, string Variant, Ksw2Option Flag, int W, int ZDrop, int EndBonus, int Q, int E, int Q2, int E2, int M, sbyte[] Mat, byte[] Query, byte[] Target,
        string Status, int Max, bool ZDropped, int MaxQ, int MaxT, int Mqe, int MqeT, int Mte, int MteQ, int Score, int NCigar, bool ReachEnd, string Cigar, string Sse2)
    {
        public bool Simd => Variant.EndsWith('2');
        public bool Crashed => Status != "ok";
    }

    private static int I(string s) => int.Parse(s, CultureInfo.InvariantCulture);

    private static byte[] Digits(string s) => s == "-" ? [] : [.. s.Select(c => (byte)(c - '0'))];

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        using var stream = TestDataResource.Open(ResourceName) ?? throw new InvalidOperationException($"Embedded oracle '{ResourceName}' not found.");
        using var reader = new StreamReader(stream);
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#')
                continue;
            var f = line.Split('\t');
            bool ok = f[14] == "ok";
            list.Add(new Row(I(f[0]), f[1], (Ksw2Option)I(f[2]), I(f[3]), I(f[4]), I(f[5]), I(f[6]), I(f[7]), I(f[8]), I(f[9]), I(f[10]), [.. f[11].Split(',').Select(x => sbyte.Parse(x, CultureInfo.InvariantCulture))],
                Digits(f[12]), Digits(f[13]), f[14],
                ok ? I(f[15]) : 0, ok && I(f[16]) != 0, ok ? I(f[17]) : 0, ok ? I(f[18]) : 0, ok ? I(f[19]) : 0, ok ? I(f[20]) : 0, ok ? I(f[21]) : 0, ok ? I(f[22]) : 0, ok ? I(f[23]) : 0, ok ? I(f[24]) : 0,
                ok && I(f[25]) != 0, f[26], f[27]));
        }

        return list;
    });

    private static Ksw2Result Run(Row r) => r.Variant switch
    {
        "extz" => Ksw2Aligner.ExtendZ(r.Query, r.Target, r.M, r.Mat, r.Q, r.E, r.W, r.ZDrop, r.Flag),
        "extz2" => Ksw2Aligner.ExtendZ2(r.Query, r.Target, r.M, r.Mat, r.Q, r.E, r.W, r.ZDrop, r.EndBonus, r.Flag),
        "extd" => Ksw2Aligner.ExtendD(r.Query, r.Target, r.M, r.Mat, r.Q, r.E, r.Q2, r.E2, r.W, r.ZDrop, r.Flag),
        _ => Ksw2Aligner.ExtendD2(r.Query, r.Target, r.M, r.Mat, r.Q, r.E, r.Q2, r.E2, r.W, r.ZDrop, r.EndBonus, r.Flag),
    };

    private static string Head(Row r) => $"#{r.Id} {r.Variant} flag={(int)r.Flag} w={r.W} zdrop={r.ZDrop} eb={r.EndBonus} q={r.Q} e={r.E} q2={r.Q2} e2={r.E2} len={r.Query.Length}x{r.Target.Length}";

    /// <summary>Every field of the reference line against the port (empty when equal).</summary>
    private static List<string> Differences(Row r, Ksw2Result g)
    {
        var d = new List<string>();
        void Check<T>(string name, T got, T expected)
        {
            if (!EqualityComparer<T>.Default.Equals(got, expected))
                d.Add($"{name} got {got} expected {expected}");
        }

        Check("max", g.Max, r.Max);
        Check("zdropped", g.ZDropped, r.ZDropped);
        Check("max_q", g.MaxQ, r.MaxQ);
        Check("max_t", g.MaxT, r.MaxT);
        Check("mqe", g.Mqe, r.Mqe);
        Check("mqe_t", g.MqeT, r.MqeT);
        Check("mte", g.Mte, r.Mte);
        Check("mte_q", g.MteQ, r.MteQ);
        Check("score", g.Score, r.Score);
        Check("n_cigar", g.Cigar.Count, r.NCigar);
        Check("reach_end", g.ReachEnd, r.ReachEnd);
        Check("cigar", g.CigarString, r.Cigar);
        return d;
    }

    /// <summary>
    /// The refusal of a scalar case, derived from the condition and the run order of the harness (one process; after an abort it resumed at the next id). A case is refused when
    /// <list type="bullet">
    /// <item><b>culprit</b>: the port reports <see cref="Ksw2ReferenceUndefined"/> (the C code writes beyond <c>eh</c>; its printed result is not a reference value: the writes land on the header and the first bytes of the traceback matrix allocated right after <c>eh</c>);</item>
    /// <item><b>victim</b>: the real program aborted in this case (<c>crash:N</c>, no result line) although the case itself is well defined: the glibc abort (<c>malloc(): invalid size</c>, <c>double free or corruption</c>, <c>free(): invalid next size</c>) comes from the heap that
    /// a culprit of the same process (since the last abort) corrupted; the refusal needs such a culprit at or before the victim in its segment.</item>
    /// </list>
    /// </summary>
    private static Dictionary<int, string> RefusedScalar(List<Row> allRows, Dictionary<int, Ksw2Result> results)
    {
        var refused = new Dictionary<int, string>();
        int segmentStart = 0;
        bool culpritInSegment = false;
        foreach (var r in allRows.Where(r => !r.Simd))       // the SSE ids are skipped by the scalar harness
        {
            bool culprit = results[r.Id].ReferenceUndefined != Ksw2ReferenceUndefined.None;
            if (culprit)
            {
                refused[r.Id] = "culprit";
                culpritInSegment = true;
            }

            if (r.Crashed)
            {
                Assert.That(culpritInSegment, Is.True, $"{Head(r)}: the real program aborted (segment starting at id {segmentStart}) without any earlier undefined-behaviour case in the process: the abort is not explained");
                refused.TryAdd(r.Id, "victim");
                segmentStart = r.Id + 1;
                culpritInSegment = false;          // the next process starts with a clean heap
            }
        }

        return refused;
    }

    private static readonly Lazy<Dictionary<int, Ksw2Result>> Results = new(() => Oracle.Value.ToDictionary(r => r.Id, Run));

    private static readonly Lazy<Dictionary<int, string>> Refused = new(() => RefusedScalar(Oracle.Value, Results.Value));

    private static void CheckVariant(string variant)
    {
        var failures = new List<string>();
        var rows = Oracle.Value.Where(r => r.Variant == variant).ToList();
        int compared = 0;
        foreach (var r in rows)
        {
            var g = Results.Value[r.Id];
            if (r.Simd)
                Assert.That(g.ReferenceUndefined, Is.EqualTo(Ksw2ReferenceUndefined.None), $"{Head(r)}: a SIMD function is never flagged (ASan / UBSan clean on all 1 000 SIMD cases)");
            if (Refused.Value.ContainsKey(r.Id))
                continue;
            Assert.That(r.Crashed, Is.False, Head(r));
            compared++;
            var d = Differences(r, g);
            if (d.Count > 0)
                failures.Add($"{Head(r)}: {string.Join("; ", d)}");
        }

        Assert.That(rows, Has.Count.EqualTo(500));
        Assert.That(compared, Is.GreaterThan(400));
        TestContext.Out.WriteLine($"{variant}: {rows.Count} cases, {compared} compared and equal, {rows.Count - compared} refused");
        Assert.That(failures, Is.Empty, $"{failures.Count} of {compared} compared\n" + string.Join("\n", failures.Take(25)));
    }

    [Test]
    public void ExtZ_ScalarMatchesTheRealProgram() => CheckVariant("extz");

    [Test]
    public void ExtZ2Sse_MatchesTheRealProgram() => CheckVariant("extz2");

    [Test]
    public void ExtD_ScalarMatchesTheRealProgram() => CheckVariant("extd");

    [Test]
    public void ExtD2Sse_MatchesTheRealProgram() => CheckVariant("extd2");

    [Test]
    public void TheOracleHoldsTheTwoThousandCasesAndBothSseBuildsAgree()
    {
        var rows = Oracle.Value;
        Assert.That(rows, Has.Count.EqualTo(2000));
        Assert.That(rows.Select(r => r.Id), Is.EqualTo(Enumerable.Range(0, 2000)));
        Assert.That(rows.GroupBy(r => r.Variant).Select(g => (g.Key, g.Count())), Is.EquivalentTo(new[] { ("extz", 500), ("extz2", 500), ("extd", 500), ("extd2", 500) }));
        // the SSE2-only build (KSW_SSE2_ONLY) printed the same line as the SSE4.1 build on all 1 000 SIMD cases
        Assert.That(rows.Where(r => r.Simd).Select(r => r.Sse2).Distinct(), Is.EqualTo(new[] { "=" }));
        Assert.That(rows.Where(r => !r.Simd).Select(r => r.Sse2).Distinct(), Is.EqualTo(new[] { "-" }));
        Assert.That(rows.Where(r => r.Simd).Any(r => r.Crashed), Is.False);
    }

    // ---- the refused set -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// <summary>Independent of the port: the scalar function reaches the overflowing row <c>i = w + qlen + 1</c> (the first row with <c>st &gt; qlen</c>) only when it exists (<c>tlen &gt; w + qlen + 1</c>).</summary>
    private static bool OverflowRowExists(Row r)
    {
        int w = r.W < 0 || r.W > Math.Max(r.Query.Length, r.Target.Length) ? Math.Max(r.Query.Length, r.Target.Length) : r.W;
        return r.Target.Length - 1 >= w + r.Query.Length + 1;
    }

    [Test]
    public void RefusedSet_IsExactlyTheHeapOverflowCondition()
    {
        var rows = Oracle.Value.Where(r => !r.Simd).ToList();
        var culprits = rows.Where(r => Results.Value[r.Id].ReferenceUndefined != Ksw2ReferenceUndefined.None).ToList();
        var victims = Refused.Value.Where(kv => kv.Value == "victim").Select(kv => kv.Key).ToList();

        // 1. only the heap-overflow flag occurs, only on cases that have the overflowing row (tlen > w + qlen + 1)
        Assert.That(culprits.All(r => Results.Value[r.Id].ReferenceUndefined == Ksw2ReferenceUndefined.WritesBeyondRowBuffer), Is.True);
        Assert.That(culprits.All(OverflowRowExists), Is.True);

        // 2. without a Z-drop every row is run, so the condition decides alone: the cases with the row are exactly the culprits
        var noZdrop = rows.Where(r => r.ZDrop < 0).ToList();
        Assert.That(noZdrop.Where(OverflowRowExists).Select(r => r.Id), Is.EquivalentTo(noZdrop.Where(culprits.Contains).Select(r => r.Id)));

        // 3. with a Z-drop the loop may stop before the row: then the port says the extension was dropped
        foreach (var r in rows.Where(r => r.ZDrop >= 0 && OverflowRowExists(r) && !culprits.Contains(r)))
            Assert.That(Results.Value[r.Id].ZDropped, Is.True, $"{Head(r)}: the overflow row exists but is not reached although the loop did not stop by Z-drop");

        // 4. every abort of the real program (33, exit 134 / 139) is explained: the culprit itself or a victim of an earlier culprit of the same process (RefusedScalar asserts this for each of them)
        var crashed = rows.Where(r => r.Crashed).ToList();
        Assert.That(crashed, Has.Count.EqualTo(33));
        Assert.That(crashed.Select(r => r.Status).Distinct().Order(), Is.EqualTo(new[] { "crash:134", "crash:139" }));
        Assert.That(crashed.Select(r => r.Id).Except(Refused.Value.Keys), Is.Empty);

        // 5. the refused set is exactly culprits + victims, nothing else: every other scalar case has a result line and is compared
        Assert.That(Refused.Value.Keys.Order(), Is.EqualTo(culprits.Select(r => r.Id).Union(victims).Order()));
        Assert.That(rows.Where(r => !Refused.Value.ContainsKey(r.Id)).Any(r => r.Crashed), Is.False);
        // information only: the culprits that printed a result (the 8 bytes landed in the allocator's slack or on the traceback matrix): whether the port happens to equal that print
        var printed = culprits.Where(r => !r.Crashed).ToList();
        TestContext.Out.WriteLine($"culprits that printed a result: {printed.Count}, port equal on {printed.Count(r => Differences(r, Results.Value[r.Id]).Count == 0)}");
        TestContext.Out.WriteLine($"scalar cases {rows.Count}: refused {Refused.Value.Count} (culprits {culprits.Count}, aborted {crashed.Count}, victims {victims.Count}), compared {rows.Count - Refused.Value.Count}");
    }
}
