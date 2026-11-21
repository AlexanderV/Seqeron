using System.Globalization;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    /// <summary>
    /// Reads the JSON description of a pairwise coordinate alignment: <c>{"sequences": [{"id", "length", "parts": {"0": "ACGT..."} | null}, ...], "coordinates": [[...], [...]], "attrs": {"score": 12, ...}}</c>
    /// (<c>parts</c> = the known residues by start position, null = the contents are unknown; <c>attrs</c> = the attributes of the alignment object).
    /// </summary>
    public static BioCoordAlignment FromJson(JsonNode spec)
    {
        ArgumentNullException.ThrowIfNull(spec);
        var seqs = new List<BioCoordSequence>();
        foreach (var s in spec["sequences"]?.AsArray() ?? throw new ArgumentException("\"sequences\" is missing.", nameof(spec)))
        {
            if (s is null)
            {
                seqs.Add(new BioCoordSequence { Id = null, Length = -1 });
                continue;
            }

            SortedDictionary<long, string>? parts = null;
            if (s!["parts"] is JsonObject po)
            {
                parts = [];
                foreach (var kv in po)
                    parts[long.Parse(kv.Key, CultureInfo.InvariantCulture)] = kv.Value!.GetValue<string>();
            }

            seqs.Add(new BioCoordSequence
            {
                Id = s["id"]?.GetValue<string>(),
                Length = s["length"] is null ? -1 : s["length"]!.GetValue<long>(),
                Parts = parts,
                Description = s["description"]?.GetValue<string>() ?? "",
                Annotations = s["annotations"] is JsonObject an ? (JsonObject)an.DeepClone() : [],
                LetterAnnotations = s["letter_annotations"] is JsonObject la ? (JsonObject)la.DeepClone() : [],
            });
        }

        if (!spec.AsObject().ContainsKey("coordinates"))
            throw new ArgumentException("\"coordinates\" is missing.", nameof(spec));
        bool none = spec["coordinates"] is null;
        var co = none ? [] : spec["coordinates"]!.AsArray().Select(r => r!.AsArray().Select(x => x!.GetValue<long>()).ToArray()).ToArray();
        return new BioCoordAlignment { Sequences = seqs, Coordinates = co, NoCoordinates = none, Attributes = spec["attrs"] is JsonObject a ? (JsonObject)a.DeepClone() : [] };
    }

    /// <summary>The JSON description of an alignment (<see cref="FromJson"/>).</summary>
    public static JsonNode ToJson(BioCoordAlignment a)
    {
        ArgumentNullException.ThrowIfNull(a);
        var seqs = new JsonArray();
        foreach (var s in a.Sequences)
        {
            if (a.NoCoordinates && ReferenceEquals(s, a.Sequences[0]) && s.Length < 0 && s.Id is null)
            {
                seqs.Add(null);
                continue;
            }

            var o = new JsonObject { ["id"] = s.Id, ["length"] = s.Length < 0 ? null : s.Length };
            if (s.Parts is null)
            {
                o["parts"] = null;
            }
            else
            {
                var po = new JsonObject();
                foreach (var kv in s.Parts)
                    po[kv.Key.ToString(CultureInfo.InvariantCulture)] = kv.Value;
                o["parts"] = po;
            }

            if (s.Description.Length > 0)
                o["description"] = s.Description;
            if (s.Annotations.Count > 0)
                o["annotations"] = s.Annotations.DeepClone();
            if (s.LetterAnnotations.Count > 0)
                o["letter_annotations"] = s.LetterAnnotations.DeepClone();
            seqs.Add(o);
        }

        return new JsonObject
        {
            ["sequences"] = seqs,
            ["coordinates"] = a.NoCoordinates ? null : new JsonArray(a.Coordinates.Select(r => (JsonNode)new JsonArray(r.Select(x => (JsonNode)JsonValue.Create(x)!).ToArray())).ToArray()),
            ["attrs"] = a.Attributes.DeepClone(),
        };
    }
}
