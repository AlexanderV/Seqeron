using System.Text.Json.Nodes;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// <c>Bio.SeqRecord.SeqRecord</c> of Biopython 1.88 as the rows of a <see cref="MultipleSeqAlignment"/> use it (B09 WP67, F307, R14-10b): the sequence text with identifier, name, description, annotations, per-letter annotations (a string or an
/// array of one entry per letter; the restricted dictionary of Biopython: any other value or a length that differs from the sequence is refused) and database cross-references. Slicing keeps id / name / description and <c>molecule_type</c> and slices the
/// letter annotations; adding two records merges them as <c>SeqRecord.__add__</c> does. Sequence features (<c>SeqFeature</c>) are not modelled.
/// </summary>
public sealed class SeqRecord
{
    /// <summary>The identifier of a record without one (<c>"&lt;unknown id&gt;"</c>).</summary>
    public const string UnknownId = "<unknown id>";

    /// <summary>The name of a record without one.</summary>
    public const string UnknownName = "<unknown name>";

    /// <summary>The description of a record without one.</summary>
    public const string UnknownDescription = "<unknown description>";

    private readonly JsonObject _letterAnnotations = [];

    /// <summary>Creates a record.</summary>
    /// <param name="seq">The sequence text (<c>record.seq</c>).</param>
    /// <param name="id">The identifier (default <c>&lt;unknown id&gt;</c>).</param>
    /// <param name="name">The name (default <c>&lt;unknown name&gt;</c>).</param>
    /// <param name="description">The description (default <c>&lt;unknown description&gt;</c>).</param>
    /// <param name="annotations"><c>record.annotations</c>.</param>
    /// <param name="letterAnnotations"><c>record.letter_annotations</c>: strings or arrays with one entry per letter.</param>
    /// <param name="dbxrefs"><c>record.dbxrefs</c>.</param>
    /// <exception cref="ArgumentException">A letter annotation is not a string or an array of the length of the sequence (Biopython's <c>TypeError</c>).</exception>
    public SeqRecord(string seq, string? id = null, string? name = null, string? description = null, JsonObject? annotations = null, JsonObject? letterAnnotations = null, IEnumerable<string>? dbxrefs = null)
    {
        Seq = seq ?? throw new ArgumentNullException(nameof(seq));
        Id = id ?? UnknownId;
        Name = name ?? UnknownName;
        Description = description ?? UnknownDescription;
        Annotations = annotations is null ? [] : (JsonObject)annotations.DeepClone();
        DbXrefs = dbxrefs is null ? [] : [.. dbxrefs];
        if (letterAnnotations is not null)
        {
            // the property setter of Biopython: every entry must have the length of the sequence, the message lists the entries that do not
            var wrong = letterAnnotations.Where(kv => !IsSequenceOfLength(kv.Value, Length)).ToList();
            if (wrong.Count > 0)
                throw new ArgumentException("The per-letter-annotations have the same length as the sequence, but found: \n " + string.Join(",", wrong.Select(kv => $"{kv.Key}={PyStr(kv.Value)}")));
            foreach (var kv in letterAnnotations)
                _letterAnnotations[kv.Key] = kv.Value?.DeepClone();
        }
    }

    // Python's str(value) of a string, a list of numbers / strings or a number
    private static string PyStr(JsonNode? v) => v switch
    {
        JsonArray a => "[" + string.Join(", ", a.Select(PyRepr)) + "]",
        JsonValue x when x.TryGetValue<string>(out var s) => s,
        _ => PyRepr(v),
    };

    private static string PyRepr(JsonNode? v)
    {
        if (v is JsonValue x)
        {
            if (x.TryGetValue<string>(out var s))
                return "'" + s.Replace("\\", "\\\\").Replace("'", "\\'") + "'";
            if (x.TryGetValue<bool>(out var b))
                return b ? "True" : "False";
            string text = x.ToJsonString();
            if (!text.Contains('.') && !text.Contains('E') && !text.Contains('e'))
                return text;
            double d = x.GetValue<double>();
            string r = d.ToString("R", System.Globalization.CultureInfo.InvariantCulture);
            return d == Math.Floor(d) && !r.Contains('E') ? r + ".0" : r;
        }

        return v is null ? "None" : v.ToJsonString();
    }

    /// <summary>The sequence (<c>str(record.seq)</c>).</summary>
    public string Seq { get; }

    /// <summary><c>len(record)</c>.</summary>
    public int Length => Seq.Length;

    /// <summary><c>record.id</c>.</summary>
    public string Id { get; set; }

    /// <summary><c>record.name</c>.</summary>
    public string Name { get; set; }

    /// <summary><c>record.description</c>.</summary>
    public string Description { get; set; }

    /// <summary><c>record.annotations</c>.</summary>
    public JsonObject Annotations { get; }

    /// <summary><c>record.dbxrefs</c>.</summary>
    public List<string> DbXrefs { get; }

    /// <summary><c>record.letter_annotations</c> (read; set entries with <see cref="SetLetterAnnotation"/>).</summary>
    public JsonObject LetterAnnotations => _letterAnnotations;

    /// <summary><c>record.letter_annotations[key] = value</c> of the restricted dictionary.</summary>
    /// <exception cref="ArgumentException">The value is not a string or an array of the length of the sequence.</exception>
    public void SetLetterAnnotation(string key, JsonNode? value)
    {
        if (!IsSequenceOfLength(value, Length))
            throw new ArgumentException($"Any per-letter annotation should be a Python sequence (list, tuple or string) of the same length as the biological sequence, here {Length}.");
        _letterAnnotations[key] = value!.DeepClone();
    }

    internal static bool IsSequenceOfLength(JsonNode? value, int length) => value switch
    {
        JsonArray a => a.Count == length,
        JsonValue v when v.TryGetValue<string>(out var s) => s.Length == length,
        _ => false,
    };

    /// <summary><c>record[i]</c>: one letter (negative indices count from the end).</summary>
    /// <exception cref="ArgumentOutOfRangeException">Python's <c>IndexError: index out of range</c>.</exception>
    public string this[int index]
    {
        get
        {
            int i = index < 0 ? index + Length : index;
            if (i < 0 || i >= Length)
                throw new ArgumentOutOfRangeException(nameof(index), "index out of range");
            return Seq[i].ToString();
        }
    }

    /// <summary><c>record[start:stop:step]</c>: the sliced letters with id, name, description, <c>molecule_type</c> and the sliced letter annotations.</summary>
    public SeqRecord Slice(PySlice slice)
    {
        var indices = slice.Range(Length);
        var answer = new SeqRecord(new string([.. indices.Select(i => Seq[i])]), Id, Name, Description);
        if (Annotations.TryGetPropertyValue("molecule_type", out var molecule))
            answer.Annotations["molecule_type"] = molecule?.DeepClone();
        foreach (var kv in _letterAnnotations)
            answer._letterAnnotations[kv.Key] = Pick(kv.Value!, indices);
        return answer;
    }

    internal static JsonNode Pick(JsonNode value, IReadOnlyList<int> indices)
    {
        if (value is JsonArray array)
            return new JsonArray([.. indices.Select(i => array[i]?.DeepClone())]);
        string s = value.GetValue<string>();
        return JsonValue.Create(new string([.. indices.Select(i => s[i])]))!;
    }

    /// <summary><c>record + other</c> (<c>SeqRecord.__add__</c> of two records): the sequences are joined; id, name, description and every annotation the two have in common are kept, the cross-references are united,
    /// letter annotations present in both are concatenated.</summary>
    /// <exception cref="InvalidOperationException">A letter annotation is a string in one record and an array in the other (Python's <c>TypeError</c>).</exception>
    public SeqRecord Add(SeqRecord other)
    {
        ArgumentNullException.ThrowIfNull(other);
        var answer = new SeqRecord(Seq + other.Seq, dbxrefs: DbXrefs);
        foreach (var r in other.DbXrefs)
            if (!answer.DbXrefs.Contains(r))
                answer.DbXrefs.Add(r);
        if (Id == other.Id)
            answer.Id = Id;
        if (Name == other.Name)
            answer.Name = Name;
        if (Description == other.Description)
            answer.Description = Description;
        foreach (var kv in Annotations)
            if (other.Annotations.TryGetPropertyValue(kv.Key, out var v) && JsonNode.DeepEquals(v, kv.Value))
                answer.Annotations[kv.Key] = kv.Value?.DeepClone();
        foreach (var kv in _letterAnnotations)
            if (other._letterAnnotations.TryGetPropertyValue(kv.Key, out var v))
                answer._letterAnnotations[kv.Key] = Concatenate(kv.Value!, v!);
        return answer;
    }

    // Python's v + w of two strings or two lists
    internal static JsonNode Concatenate(JsonNode a, JsonNode b)
    {
        if (a is JsonArray x && b is JsonArray y)
            return new JsonArray([.. x.Select(n => n?.DeepClone()), .. y.Select(n => n?.DeepClone())]);
        if (a is JsonValue && b is JsonValue && a.GetValueKind() == System.Text.Json.JsonValueKind.String && b.GetValueKind() == System.Text.Json.JsonValueKind.String)
            return JsonValue.Create(a.GetValue<string>() + b.GetValue<string>())!;
        throw new InvalidOperationException(a is JsonArray ? "can only concatenate list (not \"str\") to list" : "can only concatenate str (not \"list\") to str");
    }

    /// <summary>The record as text of one line (<c>"&lt;seq&gt; &lt;id&gt;"</c>, the sequence cut to 44 letters, <c>...</c> and the last 3 above 50).</summary>
    internal string ToRowText(int length = 50) => Seq.Length <= length ? $"{Seq} {Id}" : $"{Seq[..(length - 6)]}...{Seq[^3..]} {Id}";
}
