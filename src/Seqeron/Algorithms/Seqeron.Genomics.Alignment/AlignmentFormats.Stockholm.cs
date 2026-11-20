using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// Bio/Align/stockholm.py AlignmentIterator
public static partial class AlignmentFormats
{
    private sealed class StockholmReference
    {
        internal long Number;
        internal string? Comment;
        internal string? Medline;
        internal List<string> Title = [];
        internal List<string> Author = [];
        internal List<string> Location = [];
        internal List<string> Order = [];
    }

    // a number that is not a 64-bit integer is a FormatException (int.Parse threw OverflowException / FormatException of the BCL)
    private static long StockholmInt(string text, string what)
    {
        // Python's int is unbounded; 64 bits cover every real file, a larger number is a FormatException
        if (!long.TryParse(text, NumberStyles.Integer, CultureInfo.InvariantCulture, out long value))
            throw new FormatException($"Stockholm {what}: '{text}' is not a valid integer (at most {long.MaxValue}).");
        return value;
    }

    private static List<BioAlignment> ParseStockholm(IReadOnlyList<string> fileLines)
    {
        var result = new List<BioAlignment>();
        List<string> ids = [];
        List<string> aligned = [];
        List<StockholmReference> references = [];
        List<string> referenceComments = [];
        List<JsonObject> databaseReferences = [];
        List<JsonObject> nestedDomains = [];
        List<KeyValuePair<string, List<string>>> gf = [];
        List<KeyValuePair<string, string>> gc = [];
        List<KeyValuePair<string, Dictionary<string, List<string>>>> gs = [];
        List<KeyValuePair<string, List<KeyValuePair<string, string>>>> gr = [];
        int length = -1;
        string seqname = "";
        bool started = false;
        foreach (string raw in fileLines)
        {
            string line = raw.Trim();
            if (line.Length == 0)
                continue;
            if (line == "# STOCKHOLM 1.0")
            {
                started = true;
                ids = [];
                aligned = [];
                references = [];
                referenceComments = [];
                databaseReferences = [];
                nestedDomains = [];
                gf = [];
                gc = [];
                gs = [];
                gr = [];
                length = -1;
            }
            else if (!started)
            {
                continue;
            }
            else if (line == "//")
            {
                result.Add(FinishStockholm(ids, aligned, references, databaseReferences, nestedDomains, gf, gc, gs, gr));
                started = false;
            }
            else if (!line.StartsWith('#'))
            {
                var parts = SplitWhitespace(line, 1);
                if (parts.Length != 2)
                    throw new FormatException("Could not split line into sequence name and aligned sequence:\n" + line);
                seqname = parts[0];
                string row = parts[1];
                if (length < 0)
                    length = row.Length;
                else if (length != row.Length)
                    throw new FormatException($"Aligned sequence {seqname} consists of {row.Length} letters, expected {length} letters)");
                aligned.Add(row.Replace('.', '-'));
                ids.Add(seqname);
            }
            else if (line.StartsWith("#=GF ", StringComparison.Ordinal))
            {
                var parts = SplitWhitespace(line[5..].Trim(), 1);
                string feature = parts[0];
                string text = parts.Length > 1 ? parts[1] : throw new FormatException("Empty #=GF line.");
                switch (feature)
                {
                    case "RN":
                        references.Add(new StockholmReference { Number = StockholmInt(text.Length >= 2 ? text[1..^1] : text, "#=GF RN reference number"), Comment = referenceComments.Count > 0 ? string.Join(' ', referenceComments) : null });
                        referenceComments = [];
                        break;
                    case "RM":
                        references[^1].Medline = text;
                        references[^1].Order.Add("medline");
                        break;
                    case "RT":
                        references[^1].Title.Add(text);
                        break;
                    case "RA":
                        references[^1].Author.Add(text);
                        break;
                    case "RL":
                        references[^1].Location.Add(text);
                        break;
                    case "RC":
                        referenceComments.Add(text);
                        break;
                    case "DR":
                        databaseReferences.Add(new JsonObject { ["reference"] = text });
                        break;
                    case "DC":
                        databaseReferences[^1]["comment"] = text;
                        break;
                    case "NE":
                        nestedDomains.Add(new JsonObject { ["accession"] = text });
                        break;
                    case "NL":
                        nestedDomains[^1]["location"] = text;
                        break;
                    default:
                        {
                            int k = gf.FindIndex(x => x.Key == feature);
                            if (k < 0)
                                gf.Add(new KeyValuePair<string, List<string>>(feature, [text]));
                            else
                                gf[k].Value.Add(text);
                            break;
                        }
                }
            }
            else if (line.StartsWith("#=GC ", StringComparison.Ordinal))
            {
                var parts = SplitWhitespace(line[5..].Trim(), 2);
                if (parts.Length != 2)
                    throw new FormatException("Could not split the #=GC line: " + line);
                int k = gc.FindIndex(x => x.Key == parts[0]);
                if (k < 0)
                    gc.Add(new KeyValuePair<string, string>(parts[0], parts[1].Trim()));
                else
                    gc[k] = new KeyValuePair<string, string>(parts[0], gc[k].Value + parts[1].Trim());
            }
            else if (line.StartsWith("#=GS ", StringComparison.Ordinal))
            {
                var parts = SplitWhitespace(line[5..].Trim(), 2);
                string name = parts[0];
                string feature = parts[1];
                string text = parts.Length > 2 ? parts[2] : "";
                int k = gs.FindIndex(x => x.Key == name);
                if (k < 0)
                {
                    gs.Add(new KeyValuePair<string, Dictionary<string, List<string>>>(name, new Dictionary<string, List<string>> { ["DR"] = [] }));
                    k = gs.Count - 1;
                }

                var dict = gs[k].Value;
                if (feature == "DR")
                {
                    dict["DR"].Add(text);
                }
                else
                {
                    if (dict.ContainsKey(feature))
                        throw new FormatException($"Duplicate #=GS feature {feature} for {name}.");
                    dict[feature] = [text];
                }
            }
            else if (line.StartsWith("#=GR ", StringComparison.Ordinal))
            {
                var terms = SplitWhitespace(line[5..], 2);
                if (terms[0] != seqname)
                    throw new FormatException("#=GR line does not follow its sequence.");
                int k = gr.FindIndex(x => x.Key == seqname);
                if (k < 0)
                {
                    gr.Add(new KeyValuePair<string, List<KeyValuePair<string, string>>>(seqname, []));
                    k = gr.Count - 1;
                }

                gr[k].Value.Add(new KeyValuePair<string, string>(terms[1], terms[2].Trim()));
            }
        }

        return result;
    }

    private static BioAlignment FinishStockholm(List<string> ids, List<string> aligned, List<StockholmReference> references, List<JsonObject> databaseReferences, List<JsonObject> nestedDomains,
        List<KeyValuePair<string, List<string>>> gf, List<KeyValuePair<string, string>> gc, List<KeyValuePair<string, Dictionary<string, List<string>>>> gs, List<KeyValuePair<string, List<KeyValuePair<string, string>>>> gr)
    {
        int length = aligned.Count > 0 ? aligned[0].Length : 0;
        var skipped = new HashSet<int>();
        for (int c = 0; c < length; c++)
        {
            if (aligned.All(r => r[c] == '-'))
                skipped.Add(c);
        }

        var rows = aligned.Select(r => new string([.. r.Where((_, i) => !skipped.Contains(i))])).ToList();
        int columns = rows.Count > 0 ? rows[0].Length : 0;
        var annotations = new JsonObject();
        if (references.Count > 0)
        {
            var list = new JsonArray();
            foreach (var r in references)
            {
                var o = new JsonObject { ["number"] = r.Number };
                if (r.Comment is not null)
                    o["comment"] = r.Comment;
                if (r.Medline is not null)
                    o["medline"] = r.Medline;
                o["title"] = string.Join(' ', r.Title);
                o["author"] = string.Join(' ', r.Author);
                o["location"] = string.Join(' ', r.Location);
                list.Add(o);
            }

            annotations["references"] = list;
        }

        if (databaseReferences.Count > 0)
            annotations["database references"] = new JsonArray([.. databaseReferences.Select(d => d.DeepClone())]);
        if (nestedDomains.Count > 0)
            annotations["nested domains"] = new JsonArray([.. nestedDomains.Select(d => d.DeepClone())]);
        foreach (var (key, value) in gf)
        {
            JsonNode node;
            switch (key)
            {
                case "WK":
                    {
                        var refs = new JsonArray();
                        var reference = new StringBuilder();
                        foreach (string item in value)
                        {
                            // a line ending in '/' continues in the next one
                            if (item.EndsWith('/'))
                            {
                                reference.Append(item, 0, item.Length - 1);
                            }
                            else
                            {
                                reference.Append(item);
                                refs.Add(reference.ToString());
                                reference.Clear();
                            }
                        }

                        node = refs;
                        break;
                    }

                case "SM" or "CC" or "**":
                    node = string.Join(' ', value);
                    break;
                case "SQ":
                    if (StockholmInt(value[0], "#=GS SQ sequence count") != rows.Count)
                        throw new FormatException("Inconsistent number of sequences in alignment");
                    continue;
                case "AU":
                    node = new JsonArray([.. value.Select(v => (JsonNode)JsonValue.Create(v)!)]);
                    break;
                default:
                    node = value[^1];
                    break;
            }

            if (StockholmGfByCode.TryGetValue(key, out var longName))
                annotations[longName] = node;
        }

        var column = new JsonObject();
        foreach (var (key, v) in gc)
        {
            string value = skipped.Count > 0 ? new string([.. v.Where((_, i) => !skipped.Contains(i))]) : v;
            if (value.Length != columns)
                throw new FormatException($"{key} length is {value.Length}, expected {columns}");
            column[StockholmGcByCode.TryGetValue(key, out var name) ? name : key] = value;
        }

        var recordAnnotations = new Dictionary<int, JsonObject>();
        var descriptions = ids.Select(_ => "").ToList();
        var dbxrefs = new Dictionary<int, List<string>>();
        var letters = new Dictionary<int, JsonObject>();
        foreach (var (name, features) in gs)
        {
            int idx = ids.IndexOf(name);
            if (idx < 0)
                throw new FormatException($"Failed to find seqname {name}");
            var ann = new JsonObject();
            foreach (var (key, values) in features)
            {
                if (key == "DE")
                    descriptions[idx] = values[0];
                else if (key == "DR")
                    dbxrefs[idx] = values;
                else
                    ann[StockholmGsByCode.TryGetValue(key, out var n) ? n : key] = values[0];
            }

            recordAnnotations[idx] = ann;
        }

        foreach (var (name, features) in gr)
        {
            int idx = ids.IndexOf(name);
            if (idx < 0)
                throw new FormatException($"Failed to find seqname {name}");
            var o = new JsonObject();
            foreach (var (key, value) in features)
                o[StockholmGrByCode.TryGetValue(key, out var n) ? n : key] = key == "CSA" ? value.Replace("-", "", StringComparison.Ordinal) : value.Replace(".", "", StringComparison.Ordinal);
            letters[idx] = o;
        }

        var records = new List<BioAlignmentRecord>();
        for (int i = 0; i < ids.Count; i++)
        {
            records.Add(new BioAlignmentRecord
            {
                Id = ids[i],
                Description = descriptions[i],
                Row = rows[i],
                Annotations = recordAnnotations.TryGetValue(i, out var a) ? a : [],
                LetterAnnotations = letters.TryGetValue(i, out var l) ? l : [],
                DbXrefs = dbxrefs.TryGetValue(i, out var d) ? d : [],
            });
        }

        return new BioAlignment { Records = records, Annotations = annotations, ColumnAnnotations = column };
    }
}
