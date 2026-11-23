using System.Collections;
using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The legacy <c>Bio.Align.MultipleSeqAlignment</c> of Biopython 1.88 (B09 WP67, F307, R14-10b): a list of <see cref="SeqRecord"/> rows of one length with alignment annotations and per-column annotations. Every member is a port of the Python code with its
/// errors: construction / <see cref="Append"/> / <see cref="Extend"/> (rows of equal length), <c>get_alignment_length</c>, <c>align[i]</c>, <c>align[a:b:c]</c>, <c>align[i, j]</c>, <c>align[:, j]</c>, <c>align[r, c]</c> (<see cref="Row(int)"/>, <see cref="Rows"/>,
/// <see cref="Row(int, PySlice)"/>, <see cref="Cell"/>, <see cref="Column"/>, <see cref="Slice"/>), <c>del</c>, <c>+</c>, <c>sort</c>, <c>substitutions</c>, <c>alignment</c> (the modern <see cref="BioCoordAlignment"/>) and <c>str()</c>.
/// The per-column annotations are Biopython's restricted dictionary (a string or an array with one entry per column; the length is the alignment length at the time the dictionary was created, as in the original). File output through
/// <c>Bio.AlignIO</c> is not part of the class: <see cref="ToBioAlignment"/> gives the object <see cref="AlignmentFormats.Write"/> takes.
/// </summary>
public sealed class MultipleSeqAlignment : IReadOnlyList<SeqRecord>
{
    private readonly List<SeqRecord> _records = [];
    private JsonObject? _perColumn;
    private int _perColumnLength;

    /// <summary>Creates an alignment from rows of equal length.</summary>
    /// <param name="records">The rows (may be empty).</param>
    /// <param name="annotations">Information about the whole alignment.</param>
    /// <param name="columnAnnotations">Per-column annotation: strings or arrays with one entry per column.</param>
    /// <exception cref="ArgumentException">The rows differ in length ("Sequences must all be the same length"), a column annotation has the wrong length, or column annotations are given without rows.</exception>
    public MultipleSeqAlignment(IEnumerable<SeqRecord>? records = null, JsonObject? annotations = null, JsonObject? columnAnnotations = null)
    {
        if (records is not null)
            Extend(records);
        Annotations = annotations is null ? [] : (JsonObject)annotations.DeepClone();
        SetColumnAnnotations(columnAnnotations ?? []);
    }

    /// <summary><c>alignment.annotations</c>.</summary>
    public JsonObject Annotations { get; }

    /// <summary><c>len(alignment)</c>: the number of rows.</summary>
    public int Count => _records.Count;

    /// <summary><c>alignment[i]</c>: a row (negative indices count from the end).</summary>
    /// <exception cref="ArgumentOutOfRangeException">Python's <c>IndexError: list index out of range</c>.</exception>
    public SeqRecord this[int index] => Row(index);

    /// <inheritdoc/>
    public IEnumerator<SeqRecord> GetEnumerator() => _records.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    /// <summary><c>alignment.get_alignment_length()</c>: the length of the longest row.</summary>
    public int GetAlignmentLength() => _records.Count == 0 ? 0 : _records.Max(r => r.Length);

    /// <summary><c>alignment.column_annotations</c>: a snapshot of the entries (set them with <see cref="SetColumnAnnotation"/> / <see cref="SetColumnAnnotations"/>).</summary>
    public JsonObject ColumnAnnotations => (JsonObject)PerColumn().DeepClone();

    // _get_per_column_annotations: the dictionary is created with the alignment length of that moment and keeps it
    private JsonObject PerColumn()
    {
        if (_perColumn is null)
        {
            _perColumnLength = _records.Count > 0 ? GetAlignmentLength() : 0;
            _perColumn = [];
        }

        return _perColumn;
    }

    /// <summary><c>alignment.column_annotations = value</c>: replaces the dictionary (validated against the current alignment length).</summary>
    /// <exception cref="ArgumentException">An entry is not a string or an array of the alignment length, or entries are given for an alignment without rows ("Can't set per-column-annotations without an alignment").</exception>
    public void SetColumnAnnotations(JsonObject value)
    {
        ArgumentNullException.ThrowIfNull(value);
        if (_records.Count > 0)
        {
            _perColumnLength = GetAlignmentLength();
            _perColumn = [];
            foreach (var kv in value)
                SetColumnAnnotation(kv.Key, kv.Value);
        }
        else
        {
            _perColumn = null;
            if (value.Count > 0)
                throw new ArgumentException("Can't set per-column-annotations without an alignment");
        }
    }

    /// <summary><c>alignment.column_annotations[key] = value</c> of the restricted dictionary.</summary>
    /// <exception cref="ArgumentException">The value is not a string or an array of the length the dictionary was created with.</exception>
    public void SetColumnAnnotation(string key, JsonNode? value)
    {
        var dict = PerColumn();
        if (!SeqRecord.IsSequenceOfLength(value, _perColumnLength))
            throw new ArgumentException($"Any per-letter annotation should be a Python sequence (list, tuple or string) of the same length as the biological sequence, here {_perColumnLength}.");
        dict[key] = value!.DeepClone();
    }

    /// <summary><c>alignment.extend(records)</c>: adds rows of the length of the alignment (of the first row for an empty alignment, which also resets the column annotations).</summary>
    /// <exception cref="ArgumentException">"Sequences must all be the same length".</exception>
    public void Extend(IEnumerable<SeqRecord> records)
    {
        ArgumentNullException.ThrowIfNull(records);
        using var it = records.GetEnumerator();
        int expected;
        if (_records.Count > 0)
        {
            expected = GetAlignmentLength();
        }
        else
        {
            if (!it.MoveNext())
                return;
            expected = it.Current.Length;
            AppendRecord(it.Current, expected);
            SetColumnAnnotations([]);
        }

        while (it.MoveNext())
            AppendRecord(it.Current, expected);
    }

    /// <summary><c>alignment.append(record)</c>.</summary>
    /// <exception cref="ArgumentException">"Sequences must all be the same length".</exception>
    public void Append(SeqRecord record)
    {
        if (_records.Count > 0)
            AppendRecord(record, GetAlignmentLength());
        else
            AppendRecord(record, null);
    }

    private void AppendRecord(SeqRecord record, int? expected)
    {
        ArgumentNullException.ThrowIfNull(record);
        if (expected is { } e && record.Length != e)
            throw new ArgumentException("Sequences must all be the same length");
        _records.Add(record);
    }

    /// <summary><c>alignment + other</c>: the rows are added pairwise (columns appended); common annotations are kept, per-column annotations present in both are concatenated.</summary>
    /// <exception cref="ArgumentException">The alignments have different numbers of rows.</exception>
    public MultipleSeqAlignment Add(MultipleSeqAlignment other)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (Count != other.Count)
            throw new ArgumentException("When adding two alignments they must have the same length (i.e. same number of rows)");
        var merged = _records.Zip(other._records, (l, r) => l.Add(r)).ToList();
        var annotations = new JsonObject();
        foreach (var kv in Annotations)
            if (other.Annotations.TryGetPropertyValue(kv.Key, out var v) && JsonNode.DeepEquals(v, kv.Value))
                annotations[kv.Key] = kv.Value?.DeepClone();
        var columns = new JsonObject();
        var otherColumns = other.PerColumn();
        foreach (var kv in PerColumn())
            if (otherColumns.TryGetPropertyValue(kv.Key, out var v))
                columns[kv.Key] = SeqRecord.Concatenate(kv.Value!, v!);
        return new MultipleSeqAlignment(merged, annotations, columns);
    }

    private int RowIndex(int index, string message = "list index out of range")
    {
        int i = index < 0 ? index + _records.Count : index;
        if (i < 0 || i >= _records.Count)
            throw new ArgumentOutOfRangeException(nameof(index), message);
        return i;
    }

    /// <summary><c>alignment[i]</c>: a row.</summary>
    public SeqRecord Row(int index) => _records[RowIndex(index)];

    /// <summary><c>alignment[a:b:c]</c>: a sub-alignment of those rows (the same records); the column annotations are kept when all rows are.</summary>
    public MultipleSeqAlignment Rows(PySlice rows)
    {
        var picked = rows.Range(_records.Count).Select(i => _records[i]).ToList();
        var result = new MultipleSeqAlignment(picked);
        var mine = PerColumn();
        if (mine.Count > 0 && result.Count == Count)
            foreach (var kv in mine)
                result.SetColumnAnnotation(kv.Key, kv.Value);
        return result;
    }

    /// <summary><c>alignment[r, a:b:c]</c>: the part of a row as a record.</summary>
    public SeqRecord Row(int row, PySlice columns) => Row(row).Slice(columns);

    /// <summary><c>alignment[r, c]</c>: one letter.</summary>
    public string Cell(int row, int column) => Row(row)[column];

    /// <summary><c>alignment[a:b:c, j]</c>: the letters of column <paramref name="column"/> of those rows.</summary>
    public string Column(PySlice rows, int column) => string.Concat(rows.Range(_records.Count).Select(i => _records[i][column]));

    /// <summary><c>alignment[:, j]</c>: the letters of column <paramref name="column"/>.</summary>
    public string Column(int column) => Column(PySlice.All, column);

    /// <summary><c>alignment[a:b:c, d:e:f]</c>: a sub-alignment of the rows and columns; the column annotations are sliced and kept when all rows are.</summary>
    public MultipleSeqAlignment Slice(PySlice rows, PySlice columns)
    {
        var picked = rows.Range(_records.Count).Select(i => _records[i].Slice(columns)).ToList();
        var result = new MultipleSeqAlignment(picked);
        var mine = PerColumn();
        if (mine.Count > 0 && result.Count == Count)
        {
            foreach (var kv in mine)
            {
                var indices = columns.Range(kv.Value is JsonArray a ? a.Count : kv.Value!.GetValue<string>().Length);
                result.SetColumnAnnotation(kv.Key, SeqRecord.Pick(kv.Value!, indices));
            }
        }

        return result;
    }

    /// <summary><c>del alignment[i]</c>.</summary>
    public void RemoveRow(int index) => _records.RemoveAt(RowIndex(index, "list assignment index out of range"));

    /// <summary><c>del alignment[a:b:c]</c>.</summary>
    public void RemoveRows(PySlice rows)
    {
        foreach (int i in rows.Range(_records.Count).OrderByDescending(i => i))
            _records.RemoveAt(i);
    }

    /// <summary><c>alignment.sort(key, reverse)</c>: stable, by the identifier unless <paramref name="key"/> maps a row to another sort value (strings compare by code point, numbers by value).</summary>
    public void Sort(Func<SeqRecord, IComparable>? key = null, bool reverse = false)
    {
        Func<SeqRecord, IComparable> k = key ?? (r => r.Id);
        var comparer = Comparer<IComparable>.Create((a, b) => a is string x && b is string y ? string.CompareOrdinal(x, y) : a.CompareTo(b));
        var sorted = reverse ? _records.OrderByDescending(k, comparer).ToList() : _records.OrderBy(k, comparer).ToList();
        _records.Clear();
        _records.AddRange(sorted);
    }

    /// <summary>
    /// <c>alignment.substitutions</c>: how often each letter is aligned to each other letter, weighted by the product of the <c>weight</c> annotations of the two rows (default 1.0); gaps are skipped, the matrix is symmetric with the counts shared between both sides.
    /// </summary>
    /// <exception cref="InvalidOperationException">No rows (Python raises <c>TypeError</c> from <c>set.union()</c>).</exception>
    public SubstitutionCounts Substitutions()
    {
        if (_records.Count == 0)
            throw new InvalidOperationException("unbound method set.union() needs an argument");
        var letters = new SortedSet<char>(_records.SelectMany(r => r.Seq));
        letters.Remove('-');
        string alphabet = new([.. letters]);
        int n = alphabet.Length;
        var m = new double[n, n];
        for (int r1 = 0; r1 < _records.Count; r1++)
        {
            string s1 = _records[r1].Seq;
            double w1 = Weight(_records[r1]);
            for (int r2 = 0; r2 < r1; r2++)
            {
                string s2 = _records[r2].Seq;
                double w2 = Weight(_records[r2]);
                int len = Math.Min(s1.Length, s2.Length);
                for (int c = 0; c < len; c++)
                {
                    if (s1[c] == '-' || s2[c] == '-')
                        continue;
                    m[alphabet.IndexOf(s1[c]), alphabet.IndexOf(s2[c])] += w1 * w2;
                }
            }
        }

        var result = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result[i, j] = (m[i, j] + m[j, i]) / 2.0;
        return new SubstitutionCounts(alphabet, result);

        static double Weight(SeqRecord r) => r.Annotations.TryGetPropertyValue("weight", out var w) && w is not null ? w.GetValue<double>() : 1.0;
    }

    /// <summary>
    /// <c>alignment.alignment</c>: the modern alignment object (gap-free sequences with block coordinates; letter annotations reduced to the letters kept; the alignment annotations and column annotations as attributes). Biopython 1.88 copies each row
    /// with <c>copy.copy</c>, which shares the letter-annotation dictionary with the source row, and then clears it: the property EMPTIES the letter annotations of the source rows (a defect of the reference, reproduced by default).
    /// </summary>
    /// <param name="keepSourceLetterAnnotations"><c>true</c> leaves the letter annotations of the rows of this alignment untouched (the result is the same).</param>
    public BioCoordAlignment ToCoordinateAlignment(bool keepSourceLetterAnnotations = false)
    {
        var attributes = new JsonObject { ["annotations"] = Annotations.DeepClone(), ["column_annotations"] = PerColumn().DeepClone() };
        if (_records.Count == 0)
            return new BioCoordAlignment { Attributes = attributes };
        var (sequences, coordinates) = AlignmentOperations.ParsePrintedAlignment([.. _records.Select(r => r.Seq)]);
        var rows = new List<BioCoordSequence>();
        for (int k = 0; k < _records.Count; k++)
        {
            var record = _records[k];
            var letters = new JsonObject();
            var kept = Enumerable.Range(0, record.Length).Where(i => record.Seq[i] != '-').ToList();
            foreach (var kv in record.LetterAnnotations)
                letters[kv.Key] = SeqRecord.Pick(kv.Value!, kept);
            if (!keepSourceLetterAnnotations)
                record.LetterAnnotations.Clear();
            rows.Add(new BioCoordSequence
            {
                Id = record.Id,
                Length = sequences[k].Length,
                Parts = AlignmentCoordinateFormats.SeqParts([new KeyValuePair<long, string>(0, sequences[k])], sequences[k].Length),
                Description = record.Description,
                Annotations = (JsonObject)record.Annotations.DeepClone(),
                LetterAnnotations = letters,
            });
        }

        return new BioCoordAlignment { Sequences = rows, Coordinates = coordinates, Attributes = attributes };
    }

    /// <summary>The alignment as the object of the text-format writers (<see cref="AlignmentFormats.Write"/>).</summary>
    public BioAlignment ToBioAlignment() => new()
    {
        Records = [.. _records.Select(r => new BioAlignmentRecord
        {
            Id = r.Id,
            Description = r.Description,
            Row = r.Seq,
            Annotations = (JsonObject)r.Annotations.DeepClone(),
            LetterAnnotations = (JsonObject)r.LetterAnnotations.DeepClone(),
            DbXrefs = [.. r.DbXrefs],
        })],
        Annotations = (JsonObject)Annotations.DeepClone(),
        ColumnAnnotations = (JsonObject)PerColumn().DeepClone(),
    };

    /// <summary><c>str(alignment)</c>: the summary line and up to 20 rows (above 20, the first 18, <c>...</c> and the last) with at most 50 letters each.</summary>
    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.Append($"Alignment with {_records.Count} rows and {GetAlignmentLength()} columns");
        if (_records.Count <= 20)
        {
            foreach (var r in _records)
                sb.Append('\n').Append(r.ToRowText());
        }
        else
        {
            foreach (var r in _records.Take(18))
                sb.Append('\n').Append(r.ToRowText());
            sb.Append("\n...\n").Append(_records[^1].ToRowText());
        }

        return sb.ToString();
    }
}
