using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit R14-10, B09 F303): the <c>Bio.Align.Alignment</c> operations (<see cref="AlignmentOperations"/>, <see cref="AlignmentCounts"/>) against Biopython 1.88 on 572 recorded cases
/// (<c>scripts/alignment/generate_alignment_ops_oracle.py</c>): shape / length, rows, <c>__array__</c>, aligned, indices, inverse_indices, frequencies, substitutions, counts with every kind of scoring, reverse_complement, sort,
/// <c>str(alignment)</c>, 12 index expressions per alignment (row / column / slice / list combinations), <c>+</c>, the ordering of alignments, map, mapall, from_alignments_with_same_reference, parse_printed_alignment.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class AlignmentOperations_Oracle_Tests
{
    private static readonly Lazy<List<JsonObject>> Records = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.alignment_ops_oracle.jsonl";
        using var stream = TestDataResource.Open(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        var list = new List<JsonObject>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length > 0)
                list.Add((JsonObject)JsonNode.Parse(line)!);
        }

        return list;
    });

    private static IEnumerable<TestCaseData> Of(params string[] kinds) =>
        Records.Value.Where(r => kinds.Contains(r["kind"]!.GetValue<string>())).Select(r => new TestCaseData(r.ToJsonString()).SetName($"{r["kind"]}_{r["id"]}"));

    private static IEnumerable<TestCaseData> Alignments() => Of("pairwise", "pairwise_rc", "msa", "msa_rc", "protein_pairwise", "protein_msa");

    private static BioCoordAlignment Make(JsonNode spec)
    {
        var sequences = new List<BioCoordSequence>();
        foreach (var s in spec["sequences"]!.AsArray())
        {
            string seq = s!["seq"]!.GetValue<string>();
            var annotations = new JsonObject();
            if (s["weight"] is JsonNode w)
                annotations["weight"] = w.GetValue<double>();
            sequences.Add(new BioCoordSequence
            {
                Id = s["id"]?.GetValue<string>(),
                Length = seq.Length,
                Parts = AlignmentCoordinateFormats.SeqParts([new KeyValuePair<long, string>(0, seq)], seq.Length),
                Annotations = annotations,
            });
        }

        var coordinates = spec["coordinates"]!.AsArray().Select(r => r!.AsArray().Select(x => x!.GetValue<long>()).ToArray()).ToArray();
        return new BioCoordAlignment { Sequences = sequences, Coordinates = coordinates };
    }

    private static string Text(BioCoordSequence s) => s.Slice(0, s.Length)!;

    private static void AssertSpec(JsonNode expected, BioCoordAlignment actual, bool ids, string label)
    {
        var seqs = expected["sequences"]!.AsArray();
        Assert.That(actual.Sequences, Has.Count.EqualTo(seqs.Count), label + " sequence count");
        for (int i = 0; i < seqs.Count; i++)
        {
            Assert.That(Text(actual.Sequences[i]), Is.EqualTo(seqs[i]!["seq"]!.GetValue<string>()), $"{label} sequence {i}");
            if (ids && seqs[i]!["id"] is JsonNode id)
                Assert.That(actual.Sequences[i].Id, Is.EqualTo(id.GetValue<string>()), $"{label} id {i}");
        }

        var coordinates = expected["coordinates"]!.AsArray().Select(r => r!.AsArray().Select(x => x!.GetValue<long>()).ToArray()).ToArray();
        Assert.That(actual.Coordinates, Is.EqualTo(coordinates), label + " coordinates");
    }

    private static bool IsError(JsonNode? n) => n is JsonObject o && o.ContainsKey("error");

    private static void Expect(JsonNode? expected, Func<object?> actual, Action<JsonNode, object?> compare, string label)
    {
        if (IsError(expected))
        {
            Assert.That(() => actual(), NUnit.Framework.Throws.Exception, label + " should fail like Biopython (" + expected!["error"] + ")");
            return;
        }

        compare(expected!, actual());
    }

    [TestCaseSource(nameof(Alignments))]
    public void ShapeRowsArrayAndIndices_EqualBiopython(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var a = Make(rec["spec"]!);
        var ops = rec["ops"]!;
        Expect(ops["shape"], () => AlignmentOperations.Shape(a), (e, v) => Assert.That(v, Is.EqualTo((e[0]!.GetValue<int>(), e[1]!.GetValue<int>()))), "shape");
        Expect(ops["length"], () => AlignmentOperations.Length(a), (e, v) => Assert.That(v, Is.EqualTo(e.GetValue<int>())), "length");
        Expect(ops["rows"], () => Enumerable.Range(0, AlignmentOperations.Count(a)).Select(k => AlignmentOperations.Row(a, k)).ToArray(), (e, v) => Assert.That(v, Is.EqualTo(e.AsArray().Select(x => x!.GetValue<string>()).ToArray())), "rows");
        Expect(ops["array"], () => AlignmentOperations.Rows(a), (e, v) => Assert.That(v, Is.EqualTo(e.AsArray().Select(x => x!.GetValue<string>()).ToArray())), "array");
        Expect(ops["aligned"], () => AlignmentOperations.Aligned(a), (e, v) =>
        {
            var expected = e.AsArray().Select(r => r!.AsArray().Select(s => s!.AsArray().Select(x => x!.GetValue<long>()).ToArray()).ToArray()).ToArray();
            Assert.That(v, Is.EqualTo(expected));
        }, "aligned");
        Expect(ops["indices"], () => AlignmentOperations.Indices(a), (e, v) => Assert.That(v, Is.EqualTo(e.AsArray().Select(r => r!.AsArray().Select(x => x!.GetValue<int>()).ToArray()).ToArray())), "indices");
        Expect(ops["inverse_indices"], () => AlignmentOperations.InverseIndices(a), (e, v) => Assert.That(v, Is.EqualTo(e.AsArray().Select(r => r!.AsArray().Select(x => x!.GetValue<int>()).ToArray()).ToArray())), "inverse_indices");
    }

    [TestCaseSource(nameof(Alignments))]
    public void FrequenciesAndSubstitutions_EqualBiopython(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var a = Make(rec["spec"]!);
        var ops = rec["ops"]!;
        Expect(ops["frequencies"], () => AlignmentOperations.Frequencies(a), (e, v) =>
        {
            var freq = (IReadOnlyDictionary<char, double[]>)v!;
            var expected = e.AsObject();
            Assert.That(freq.Keys.Select(k => k.ToString()).OrderBy(x => x, StringComparer.Ordinal), Is.EqualTo(expected.Select(kv => kv.Key).OrderBy(x => x, StringComparer.Ordinal)));
            foreach (var (key, value) in expected)
                Assert.That(freq[key[0]], Is.EqualTo(value!.AsArray().Select(x => x!.GetValue<double>()).ToArray()).Within(1e-9), "frequency " + key);
        }, "frequencies");
        Expect(ops["substitutions"], () => AlignmentOperations.Substitutions(a), (e, v) =>
        {
            var subs = (SubstitutionCounts)v!;
            Assert.That(subs.Letters, Is.EqualTo(e["letters"]!.GetValue<string>()));
            var values = e["values"]!.AsArray();
            for (int i = 0; i < values.Count; i++)
                for (int j = 0; j < values.Count; j++)
                    Assert.That(subs.Values[i, j], Is.EqualTo(values[i]![j]!.GetValue<double>()), $"substitutions[{i},{j}]");
        }, "substitutions");
    }

    private static double? Number(JsonNode? n) => n is null ? null : n.GetValue<double>();

    private static object? CountValue(AlignmentCounts c, string name) => name switch
    {
        "open_left_insertions" => c.OpenLeftInsertions,
        "extend_left_insertions" => c.ExtendLeftInsertions,
        "open_left_deletions" => c.OpenLeftDeletions,
        "extend_left_deletions" => c.ExtendLeftDeletions,
        "open_internal_insertions" => c.OpenInternalInsertions,
        "extend_internal_insertions" => c.ExtendInternalInsertions,
        "open_internal_deletions" => c.OpenInternalDeletions,
        "extend_internal_deletions" => c.ExtendInternalDeletions,
        "open_right_insertions" => c.OpenRightInsertions,
        "extend_right_insertions" => c.ExtendRightInsertions,
        "open_right_deletions" => c.OpenRightDeletions,
        "extend_right_deletions" => c.ExtendRightDeletions,
        "aligned" => c.Aligned,
        "identities" => c.Identities,
        "mismatches" => c.Mismatches,
        "positives" => c.Positives,
        "gap_score" => c.GapScore,
        "substitution_score" => c.SubstitutionScore,
        "score" => c.Score,
        "left_gaps" => c.LeftGaps,
        "internal_gaps" => c.InternalGaps,
        "right_gaps" => c.RightGaps,
        "gaps" => c.Gaps,
        "left_insertions" => c.LeftInsertions,
        "left_deletions" => c.LeftDeletions,
        "internal_insertions" => c.InternalInsertions,
        "internal_deletions" => c.InternalDeletions,
        "right_insertions" => c.RightInsertions,
        "right_deletions" => c.RightDeletions,
        "insertions" => c.Insertions,
        "deletions" => c.Deletions,
        "open_gaps" => c.OpenGaps,
        "extend_gaps" => c.ExtendGaps,
        "str" => c.ToString(),
        _ => throw new InvalidOperationException(name),
    };

    private static void AssertCounts(JsonNode expected, AlignmentCounts counts, string label)
    {
        foreach (var (name, value) in expected.AsObject())
        {
            object? actual = CountValue(counts, name);
            if (name == "str")
            {
                Assert.That(actual, Is.EqualTo(value!.GetValue<string>()), label + " str");
                continue;
            }

            if (value is null)
            {
                Assert.That(actual, Is.Null, $"{label} {name}");
                continue;
            }

            Assert.That(actual, Is.Not.Null, $"{label} {name}");
            Assert.That(Convert.ToDouble(actual, System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(value.GetValue<double>()).Within(1e-9), $"{label} {name}");
        }
    }

    private static PairwiseAligner MakeAligner(JsonNode cfg)
    {
        var aligner = new PairwiseAligner();
        foreach (var (key, value) in cfg["props"]!.AsObject())
        {
            double v = value!.GetValue<double>();
            if (key == "match_score")
                aligner.MatchScore = v;
            else if (key == "mismatch_score")
                aligner.MismatchScore = v;
            else
                aligner.SetGapScore(key, v);
        }

        if (cfg["matrix"] is JsonNode m)
            aligner.SubstitutionMatrix = SubstitutionMatrix.Get(m.GetValue<string>());
        if (cfg["wildcard"] is JsonNode w)
            aligner.Wildcard = w.GetValue<string>()[0];
        if (cfg["function"] is not null)
        {
            aligner.SetGapFunction("insertion_score", (position, length) => -1.0 - 0.5 * length - 0.01 * position);
            aligner.SetGapFunction("deletion_score", (position, length) => -2.0 - 0.25 * length);
        }

        return aligner;
    }

    [TestCaseSource(nameof(Alignments))]
    public void Counts_EqualBiopython_ForEveryKindOfScoring(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var a = Make(rec["spec"]!);
        var counts = rec["ops"]!["counts"]!;
        bool protein = rec["kind"]!.GetValue<string>().StartsWith("protein", StringComparison.Ordinal);
        Expect(counts["plain"], () => AlignmentOperations.Counts(a), (e, v) => AssertCounts(e, (AlignmentCounts)v!, "plain"), "counts()");
        Expect(counts["wildcard"], () => AlignmentOperations.Counts(a, 'N'), (e, v) => AssertCounts(e, (AlignmentCounts)v!, "wildcard"), "counts('N')");
        var matrix = SubstitutionMatrix.Get(protein ? "BLOSUM62" : "NUC.4.4");
        Expect(counts["matrix"], () => AlignmentOperations.Counts(a, matrix), (e, v) => AssertCounts(e, (AlignmentCounts)v!, "matrix"), "counts(matrix)");
        foreach (var cfg in counts["aligner_configs"]!.AsArray())
        {
            string name = cfg!["name"]!.GetValue<string>();
            Expect(counts["aligners"]![name], () => AlignmentOperations.Counts(a, MakeAligner(cfg)), (e, v) => AssertCounts(e, (AlignmentCounts)v!, "aligner " + name), "counts(aligner " + name + ")");
        }
    }

    [TestCaseSource(nameof(Alignments))]
    public void ReverseComplementSortAndText_EqualBiopython(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var a = Make(rec["spec"]!);
        var ops = rec["ops"]!;
        bool protein = rec["kind"]!.GetValue<string>().StartsWith("protein", StringComparison.Ordinal);
        Expect(ops["reverse_complement"], () => AlignmentOperations.ReverseComplement(a), (e, v) => AssertSpec(e, (BioCoordAlignment)v!, false, "reverse_complement"), "reverse_complement");
        Expect(ops["sort"], () => AlignmentOperations.Sort(a), (e, v) => AssertSpec(e, (BioCoordAlignment)v!, true, "sort"), "sort");
        Expect(ops["sort_reverse"], () => AlignmentOperations.Sort(a, reverse: true), (e, v) => AssertSpec(e, (BioCoordAlignment)v!, true, "sort reverse"), "sort reverse");
        Expect(ops["str"], () => AlignmentOperations.Format(a), (e, v) => Assert.That(v, Is.EqualTo(e.GetValue<string>())), "str");
        var matrix = SubstitutionMatrix.Get(protein ? "BLOSUM62" : "NUC.4.4");
        Expect(ops["str_matrix"], () => AlignmentOperations.Format(a, matrix), (e, v) => Assert.That(v, Is.EqualTo(e.GetValue<string>())), "str with matrix");
    }

    private static PySlice SliceOf(JsonNode n) => new(n[0]?.GetValue<int>(), n[1]?.GetValue<int>(), n[2]?.GetValue<int>());

    [TestCaseSource(nameof(Alignments))]
    public void IndexExpressions_EqualBiopython(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var a = Make(rec["spec"]!);
        int k = 0;
        foreach (var entry in rec["ops"]!["slices"]!.AsArray())
        {
            k++;
            int row = entry!["row"]!.GetValue<int>(), col = entry["col"]!.GetValue<int>();
            var rows = SliceOf(entry["rowsl"]!);
            var cols = SliceOf(entry["colsl"]!);
            var list = entry["collist"]!.AsArray().Select(x => x!.GetValue<int>()).ToList();
            string label = $"entry {k}";
            Expect(entry["row_col"], () => AlignmentOperations.Letter(a, row, col).ToString(), (e, v) => Assert.That(v, Is.EqualTo(e.GetValue<string>())), label + " [row, col]");
            Expect(entry["row_cols"], () => AlignmentOperations.Row(a, row, cols), (e, v) => Assert.That(v, Is.EqualTo(e.GetValue<string>())), label + " [row, slice]");
            Expect(entry["row_collist"], () => AlignmentOperations.Row(a, row, list), (e, v) => Assert.That(v, Is.EqualTo(e.GetValue<string>())), label + " [row, list]");
            Expect(entry["rows_col"], () => AlignmentOperations.Column(a, rows, col), (e, v) => Assert.That(v, Is.EqualTo(e.GetValue<string>())), label + " [slice, col]");
            Expect(entry["row_int"], () => AlignmentOperations.Row(a, row), (e, v) => Assert.That(v, Is.EqualTo(e.GetValue<string>())), label + " [row]");
            Expect(entry["rows"], () => AlignmentOperations.Slice(a, rows), (e, v) => AssertSpec(e, (BioCoordAlignment)v!, true, label + " [slice]"), label + " [slice]");
            Expect(entry["rows_cols"], () => AlignmentOperations.Slice(a, rows, cols), (e, v) => AssertSpec(e, (BioCoordAlignment)v!, true, label + " [slice, slice]"), label + " [slice, slice]");
            Expect(entry["rows_collist"], () => AlignmentOperations.Slice(a, rows, list), (e, v) => AssertSpec(e, (BioCoordAlignment)v!, true, label + " [slice, list]"), label + " [slice, list]");
        }
    }

    [TestCaseSource(nameof(Alignments))]
    public void AddAndOrdering_EqualBiopython(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var ops = rec["ops"]!;
        if (ops["spec2"] is not JsonNode spec2)
            return;
        var a = Make(rec["spec"]!);
        var b = Make(spec2);
        Expect(ops["add"], () => AlignmentOperations.Add(a, b), (e, v) => AssertSpec(e, (BioCoordAlignment)v!, false, "add"), "add");
        var cmp = ops["compare"]!;
        if (!IsError(cmp["lt"]))
        {
            int c = AlignmentOperations.Compare(a, b);
            Assert.That(c < 0, Is.EqualTo(cmp["lt"]!.GetValue<bool>()), "lt");
            Assert.That(c <= 0, Is.EqualTo(cmp["le"]!.GetValue<bool>()), "le");
            Assert.That(c > 0, Is.EqualTo(cmp["gt"]!.GetValue<bool>()), "gt");
            Assert.That(c >= 0, Is.EqualTo(cmp["ge"]!.GetValue<bool>()), "ge");
            Assert.That(AlignmentOperations.AreEqual(a, b), Is.EqualTo(cmp["eq"]!.GetValue<bool>()), "eq");
        }

        Assert.That(AlignmentOperations.AreEqual(a, a), Is.True);
    }

    [TestCaseSource(nameof(MapCases))]
    public void Map_EqualBiopython(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var first = Make(rec["spec"]!);
        var second = Make(rec["spec2"]!);
        Expect(rec["map"], () => AlignmentOperations.Map(first, second), (e, v) =>
        {
            var actual = (BioCoordAlignment)v!;
            var seqs = e["sequences"]!.AsArray();
            // the mapped sequences are the target of the first and the query of the second alignment
            Assert.That(actual.Sequences, Has.Count.EqualTo(2));
            Assert.That(Text(actual.Sequences[0]), Is.EqualTo(seqs[0]!["seq"]!.GetValue<string>()));
            Assert.That(Text(actual.Sequences[1]), Is.EqualTo(seqs[1]!["seq"]!.GetValue<string>()));
            Assert.That(actual.Coordinates, Is.EqualTo(e["coordinates"]!.AsArray().Select(r => r!.AsArray().Select(x => x!.GetValue<long>()).ToArray()).ToArray()));
        }, "map");
    }

    private static IEnumerable<TestCaseData> MapCases() => Of("map");

    private static IEnumerable<TestCaseData> MapAllCases() => Of("mapall");

    [TestCaseSource(nameof(MapAllCases))]
    public void MapAll_EqualBiopython(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var msa = Make(rec["msa"]!);
        var alignments = rec["alignments"]!.AsArray().Select(x => Make(x!)).ToList();
        Expect(rec["mapall"], () => AlignmentOperations.MapAll(msa, alignments), (e, v) => AssertSpec(e, (BioCoordAlignment)v!, false, "mapall"), "mapall");
    }

    private static IEnumerable<TestCaseData> SameReferenceCases() => Of("same_reference");

    [TestCaseSource(nameof(SameReferenceCases))]
    public void FromAlignmentsWithSameReference_EqualBiopython(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var alignments = rec["alignments"]!.AsArray().Select(x => Make(x!)).ToList();
        Expect(rec["result"], () => AlignmentOperations.FromAlignmentsWithSameReference(alignments), (e, v) => AssertSpec(e, (BioCoordAlignment)v!, false, "from_alignments_with_same_reference"), "from_alignments_with_same_reference");
    }

    private static IEnumerable<TestCaseData> ParseCases() => Of("parse_printed");

    [TestCaseSource(nameof(ParseCases))]
    public void ParsePrintedAlignment_EqualBiopython(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var rows = rec["rows"]!.AsArray().Select(x => x!.GetValue<string>()).ToList();
        Expect(rec["result"], () => AlignmentOperations.ParsePrintedAlignment(rows), (e, v) =>
        {
            var (sequences, coordinates) = ((IReadOnlyList<string>, long[][]))v!;
            Assert.That(sequences, Is.EqualTo(e[0]!.AsArray().Select(x => x!.GetValue<string>()).ToList()));
            Assert.That(coordinates, Is.EqualTo(e[1]!.AsArray().Select(r => r!.AsArray().Select(x => x!.GetValue<long>()).ToArray()).ToArray()));
        }, "parse_printed_alignment");
    }

    private static IEnumerable<TestCaseData> SortingCases() => Of("sorting");

    [TestCaseSource(nameof(SortingCases))]
    public void SortedAlignments_HaveBiopythonOrder(string json)
    {
        var rec = (JsonObject)JsonNode.Parse(json)!;
        var pool = rec["alignments"]!.AsArray().Select(x => Make(x!)).ToList();
        var order = rec["order"]!.AsArray().Select(x => x!.GetValue<int>()).ToList();
        // Python's sorted() is stable; indices of equal alignments (duplicates) are those of the first occurrence in the recorded order
        var sorted = Enumerable.Range(0, pool.Count).OrderBy(i => pool[i], Comparer<BioCoordAlignment>.Create(AlignmentOperations.Compare)).ToList();
        for (int i = 0; i < sorted.Count; i++)
            Assert.That(AlignmentOperations.AreEqual(pool[sorted[i]], pool[order[i]]), Is.True, $"position {i}");
    }

    [Test]
    public void Oracle_HasTheAnnouncedInventory()
    {
        var kinds = Records.Value.GroupBy(r => r["kind"]!.GetValue<string>()).ToDictionary(g => g.Key, g => g.Count());
        Assert.That(Records.Value, Has.Count.EqualTo(572));
        Assert.That(kinds["pairwise"], Is.EqualTo(44));
        Assert.That(kinds["map"], Is.EqualTo(120));
        Assert.That(kinds["mapall"], Is.EqualTo(40));
        Assert.That(kinds["same_reference"], Is.EqualTo(62));
        Assert.That(kinds["parse_printed"], Is.EqualTo(60));
        Assert.That(kinds["sorting"], Is.EqualTo(30));
    }
}
