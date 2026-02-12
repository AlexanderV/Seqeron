using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.IO;

/// <summary>
/// Parser for EMBL flat file format (.embl, .dat).
/// EMBL format uses two-letter line type prefixes.
/// </summary>
public static partial class EmblParser
{
    #region Records

    /// <summary>Represents a complete EMBL record</summary>
    public readonly record struct EmblRecord(
        string Accession,
        string SequenceVersion,
        string DataClass,
        string MoleculeType,
        string Topology,
        string TaxonomicDivision,
        int SequenceLength,
        string Description,
        IReadOnlyList<string> Keywords,
        string Organism,
        IReadOnlyList<string> OrganismClassification,
        IReadOnlyList<Reference> References,
        IReadOnlyList<Feature> Features,
        string Sequence,
        IReadOnlyDictionary<string, string> AdditionalFields);

    /// <summary>Literature reference</summary>
    public readonly record struct Reference(
        int Number,
        string Citation,
        string Authors,
        string Title,
        string Journal,
        string CrossReference,
        string Comment);

    /// <summary>Sequence feature with location and qualifiers</summary>
    public readonly record struct Feature(
        string Key,
        Location Location,
        IReadOnlyDictionary<string, string> Qualifiers);

    /// <summary>Feature location</summary>
    public readonly record struct Location(
        int Start,
        int End,
        bool IsComplement,
        bool IsJoin,
        IReadOnlyList<(int Start, int End)> Parts,
        string RawLocation);

    #endregion

    #region Line Type Prefixes

    // Standard EMBL line types
    private const string ID = "ID";   // Identification
    private const string AC = "AC";   // Accession number
    private const string SV = "SV";   // Sequence version
    private const string DT = "DT";   // Date
    private const string DE = "DE";   // Description
    private const string KW = "KW";   // Keywords
    private const string OS = "OS";   // Organism species
    private const string OC = "OC";   // Organism classification
    private const string OG = "OG";   // Organelle
    private const string RN = "RN";   // Reference number
    private const string RC = "RC";   // Reference comment
    private const string RP = "RP";   // Reference positions
    private const string RX = "RX";   // Reference cross-reference
    private const string RG = "RG";   // Reference group
    private const string RA = "RA";   // Reference authors
    private const string RT = "RT";   // Reference title
    private const string RL = "RL";   // Reference location
    private const string DR = "DR";   // Database cross-reference
    private const string CC = "CC";   // Comments
    private const string FH = "FH";   // Feature header
    private const string FT = "FT";   // Feature table
    private const string SQ = "SQ";   // Sequence header
    private const string XX = "XX";   // Spacer line
    private const string END = "//";  // Entry terminator

    #endregion

    #region Main Parsing Methods

    /// <summary>
    /// Parses EMBL records from a file.
    /// </summary>
    public static IEnumerable<EmblRecord> ParseFile(string filePath)
    {
        if (string.IsNullOrEmpty(filePath) || !File.Exists(filePath))
            yield break;

        foreach (var record in Parse(File.ReadAllText(filePath)))
        {
            yield return record;
        }
    }

    /// <summary>
    /// Parses EMBL records from text content.
    /// Multiple records are separated by // delimiter.
    /// </summary>
    public static IEnumerable<EmblRecord> Parse(string content)
    {
        if (string.IsNullOrEmpty(content))
            yield break;

        // Split by record delimiter
        var recordTexts = content.Split(new[] { "\n//" }, StringSplitOptions.RemoveEmptyEntries);

        foreach (var recordText in recordTexts)
        {
            var trimmed = recordText.Trim();
            if (!string.IsNullOrEmpty(trimmed) && trimmed.StartsWith("ID", StringComparison.Ordinal))
            {
                var record = ParseRecord(trimmed);
                if (record.HasValue)
                    yield return record.Value;
            }
        }
    }

    /// <summary>
    /// Parses a single EMBL record.
    /// </summary>
    private static EmblRecord? ParseRecord(string text)
    {
        var lines = text.Split('\n')
            .Select(l => l.TrimEnd('\r'))
            .Where(l => !string.IsNullOrEmpty(l) && l.Length >= 2)
            .ToList();

        // Extract line groups by prefix
        var lineGroups = GroupLinesByPrefix(lines);

        // Parse ID line
        var (accession, version, dataClass, moleculeType, topology, division, seqLength) =
            ParseIdLine(GetFirstLine(lineGroups, ID));

        // Parse version from separate SV line if not in ID line
        if (string.IsNullOrEmpty(version))
        {
            version = ParseAccessionLine(GetFirstLine(lineGroups, SV));
        }
        if (string.IsNullOrEmpty(accession))
        {
            accession = ParseAccessionLine(GetLines(lineGroups, AC).FirstOrDefault() ?? "");
        }

        // Description
        var description = JoinLines(GetLines(lineGroups, DE));

        // Keywords
        var keywords = ParseKeywords(JoinLines(GetLines(lineGroups, KW)));

        // Organism
        var organism = JoinLines(GetLines(lineGroups, OS)).TrimEnd('.');
        var classification = ParseClassification(JoinLines(GetLines(lineGroups, OC)));

        // References
        var references = ParseReferences(lines);

        // Features - parse directly from lines, not via GroupLinesByPrefix
        var features = ParseFeaturesFromLines(lines);

        // Sequence
        var sequence = ParseSequence(lines);

        // Additional fields
        var additionalFields = new Dictionary<string, string>();
        foreach (var (prefix, content) in lineGroups)
        {
            if (!IsStandardPrefix(prefix))
            {
                additionalFields[prefix] = content;
            }
        }

        return new EmblRecord(
            accession,
            version,
            dataClass,
            moleculeType,
            topology,
            division,
            seqLength,
            description,
            keywords,
            organism,
            classification,
            references,
            features,
            sequence,
            additionalFields
        );
    }

    #endregion

    #region Line Parsing Helpers

    private static Dictionary<string, string> GroupLinesByPrefix(List<string> lines)
    {
        var groups = new Dictionary<string, StringBuilder>();

        foreach (var line in lines)
        {
            if (line.Length < 2) continue;

            var prefix = line.Length >= 2 ? line[..2] : line;
            var content = line.Length > 5 ? line[5..] : "";

            if (!groups.ContainsKey(prefix))
                groups[prefix] = new StringBuilder();

            if (groups[prefix].Length > 0)
                groups[prefix].Append(' ');
            groups[prefix].Append(content.Trim());
        }

        return groups.ToDictionary(kv => kv.Key, kv => kv.Value.ToString());
    }

    private static string GetFirstLine(Dictionary<string, string> groups, string prefix)
    {
        return groups.GetValueOrDefault(prefix, "");
    }

    private static IEnumerable<string> GetLines(Dictionary<string, string> groups, string prefix)
    {
        if (groups.TryGetValue(prefix, out var content))
            yield return content;
    }

    private static string JoinLines(IEnumerable<string> lines)
    {
        return string.Join(" ", lines).Trim();
    }

    private static (string Accession, string Version, string DataClass, string MoleculeType, string Topology,
        string Division, int Length) ParseIdLine(string line)
    {
        if (string.IsNullOrEmpty(line))
            return ("", "", "", "", "", "", 0);

        // Format: ACCESSION; SV VERSION; TOPOLOGY; MOLECULE; DATA_CLASS; DIVISION; LENGTH BP.
        var parts = line.Split(';').Select(p => p.Trim()).ToArray();

        string accession = parts.Length > 0 ? parts[0].Split(' ')[0] : "";
        string version = "";
        string dataClass = "";
        string moleculeType = "";
        string topology = "";
        string division = "";
        int length = 0;

        foreach (var part in parts.Skip(1))
        {
            var trimmed = part.Trim();
            if (trimmed.EndsWith("BP", StringComparison.OrdinalIgnoreCase) ||
                trimmed.EndsWith("BP.", StringComparison.OrdinalIgnoreCase))
            {
                var lengthMatch = LengthRegex().Match(trimmed);
                if (lengthMatch.Success)
                    length = int.Parse(lengthMatch.Groups[1].Value, CultureInfo.InvariantCulture);
            }
            else if (trimmed.StartsWith("SV ", StringComparison.OrdinalIgnoreCase))
            {
                version = trimmed[3..].Trim();
            }
            else if (trimmed is "linear" or "circular")
            {
                topology = trimmed;
            }
            else if (trimmed is "DNA" or "RNA" or "mRNA" or "genomic DNA" or "genomic RNA")
            {
                moleculeType = trimmed;
            }
            else if (trimmed is "STD" or "CON" or "PAT" or "EST" or "GSS" or "HTC" or "HTG" or "TSA" or "WGS")
            {
                dataClass = trimmed;
            }
            else if (trimmed.Length == 3 && char.IsUpper(trimmed[0]))
            {
                division = trimmed;
            }
        }

        return (accession, version, dataClass, moleculeType, topology, division, length);
    }

    private static string ParseAccessionLine(string line)
    {
        if (string.IsNullOrEmpty(line))
            return "";

        // Take first accession number (primary)
        return line.Split(new[] { ';', ' ' }, StringSplitOptions.RemoveEmptyEntries)
            .FirstOrDefault()?.Trim() ?? "";
    }

    private static IReadOnlyList<string> ParseKeywords(string text)
    {
        if (string.IsNullOrWhiteSpace(text) || text.Trim() == ".")
            return Array.Empty<string>();

        return text
            .Split(new[] { ';' }, StringSplitOptions.RemoveEmptyEntries)
            .Select(k => k.Trim().TrimEnd('.'))
            .Where(k => !string.IsNullOrEmpty(k))
            .ToList();
    }

    private static IReadOnlyList<string> ParseClassification(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return Array.Empty<string>();

        return text
            .Split(new[] { ';' }, StringSplitOptions.RemoveEmptyEntries)
            .Select(c => c.Trim().TrimEnd('.'))
            .Where(c => !string.IsNullOrEmpty(c))
            .ToList();
    }

    private static IReadOnlyList<Reference> ParseReferences(List<string> lines)
    {
        var references = new List<Reference>();

        int currentRefNum = 0;
        string citation = "", authors = "", title = "", journal = "", xref = "", comment = "";

        foreach (var line in lines)
        {
            if (line.Length < 5) continue;

            var prefix = line[..2];
            var content = line.Length > 5 ? line[5..].Trim() : "";

            switch (prefix)
            {
                case RN:
                    // Save previous reference
                    if (currentRefNum > 0)
                    {
                        references.Add(new Reference(currentRefNum, citation, authors, title, journal, xref, comment));
                    }
                    // Parse new reference number
                    var numMatch = ReferenceNumberRegex().Match(content);
                    currentRefNum = numMatch.Success ? int.Parse(numMatch.Groups[1].Value, CultureInfo.InvariantCulture) : 0;
                    citation = authors = title = journal = xref = comment = "";
                    break;
                case RC:
                    comment = string.IsNullOrEmpty(comment) ? content : comment + " " + content;
                    break;
                case RA:
                    authors = string.IsNullOrEmpty(authors) ? content : authors + " " + content;
                    break;
                case RT:
                    title = string.IsNullOrEmpty(title) ? content : title + " " + content;
                    break;
                case RL:
                    journal = string.IsNullOrEmpty(journal) ? content : journal + " " + content;
                    break;
                case RX:
                    xref = string.IsNullOrEmpty(xref) ? content : xref + "; " + content;
                    break;
            }
        }

        // Save last reference
        if (currentRefNum > 0)
        {
            references.Add(new Reference(currentRefNum, citation, authors.TrimEnd(';', ' '),
                title.Trim('"', ';', ' '), journal.TrimEnd(';', ' '), xref, comment));
        }

        return references;
    }

    /// <summary>
    /// Parses features directly from FT lines, preserving whitespace structure.
    /// EMBL FT format:
    /// FT   key             location
    /// FT                   /qualifier="value"
    /// </summary>
    private static IReadOnlyList<Feature> ParseFeaturesFromLines(List<string> lines)
    {
        var features = new List<Feature>();

        string? currentKey = null;
        string currentLocation = "";
        var currentQualifiers = new Dictionary<string, string>();
        var qualifierBuilder = new StringBuilder();
        string? currentQualifierName = null;

        foreach (var line in lines)
        {
            if (!line.StartsWith("FT", StringComparison.Ordinal) || line.Length < 5)
                continue;

            // FT lines: positions 0-1 = "FT", 2-4 = spaces, 5+ = content
            // Feature key starts at position 5, qualifiers are indented further
            var content = line.Length > 5 ? line[5..] : "";

            if (string.IsNullOrWhiteSpace(content))
                continue;

            // Check if this is a new feature (key starts at column 5, not indented)
            // or a continuation (qualifiers/location continuation are indented)
            bool isNewFeature = content.Length > 0 && !char.IsWhiteSpace(content[0]) && !content.TrimStart().StartsWith('/');
            bool isQualifier = content.TrimStart().StartsWith('/');

            if (isNewFeature)
            {
                // Save previous feature
                if (currentKey != null)
                {
                    // Finish any pending qualifier
                    FinishQualifier(currentQualifierName, qualifierBuilder, currentQualifiers);
                    features.Add(CreateFeature(currentKey, currentLocation, currentQualifiers));
                }

                // Parse new feature: key and location
                var trimmed = content.Trim();
                var spaceIdx = trimmed.IndexOf(' ');
                if (spaceIdx > 0)
                {
                    currentKey = trimmed[..spaceIdx];
                    currentLocation = trimmed[(spaceIdx + 1)..].Trim();
                }
                else
                {
                    currentKey = trimmed;
                    currentLocation = "";
                }
                currentQualifiers = new Dictionary<string, string>();
                qualifierBuilder.Clear();
                currentQualifierName = null;
            }
            else if (isQualifier)
            {
                // Finish previous qualifier if any
                FinishQualifier(currentQualifierName, qualifierBuilder, currentQualifiers);

                // Parse new qualifier
                var qualContent = content.TrimStart();
                var match = QualifierRegex().Match(qualContent);
                if (match.Success)
                {
                    currentQualifierName = match.Groups[1].Value;
                    var value = match.Groups[2].Success ? match.Groups[2].Value : "true";
                    qualifierBuilder.Clear();
                    qualifierBuilder.Append(value);
                }
            }
            else
            {
                // Continuation line - could be location or qualifier value continuation
                var trimmed = content.Trim();
                if (currentQualifierName != null)
                {
                    // Qualifier value continuation
                    qualifierBuilder.Append(trimmed);
                }
                else if (currentKey != null)
                {
                    // Location continuation
                    currentLocation += trimmed;
                }
            }
        }

        // Save last feature
        if (currentKey != null)
        {
            FinishQualifier(currentQualifierName, qualifierBuilder, currentQualifiers);
            features.Add(CreateFeature(currentKey, currentLocation, currentQualifiers));
        }

        return features;
    }

    private static void FinishQualifier(string? name, StringBuilder value, Dictionary<string, string> qualifiers)
    {
        if (name != null)
        {
            var val = value.ToString().Trim().Trim('"');
            qualifiers[name] = val;
        }
    }

    private static IReadOnlyList<Feature> ParseFeatures(IEnumerable<string> featureLines)
    {
        var features = new List<Feature>();
        var allContent = string.Join(" ", featureLines);

        if (string.IsNullOrWhiteSpace(allContent))
            return features;

        // Parse FT lines - format: key location /qualifier=value
        string currentKey = "";
        string currentLocation = "";
        var qualifiers = new Dictionary<string, string>();

        // Re-process raw lines
        var lines = allContent.Split(new[] { "FT " }, StringSplitOptions.RemoveEmptyEntries);

        foreach (var segment in lines)
        {
            var trimmed = segment.Trim();
            if (string.IsNullOrEmpty(trimmed)) continue;

            // Check if this is a new feature (starts with key) or qualifier
            if (!trimmed.StartsWith('/'))
            {
                // Save previous feature
                if (!string.IsNullOrEmpty(currentKey))
                {
                    features.Add(CreateFeature(currentKey, currentLocation, qualifiers));
                }

                // Parse new feature
                var spaceIdx = trimmed.IndexOf(' ');
                if (spaceIdx > 0)
                {
                    currentKey = trimmed[..spaceIdx];
                    var rest = trimmed[(spaceIdx + 1)..];

                    // Find location (up to first / or end)
                    var qualStart = rest.IndexOf('/');
                    if (qualStart > 0)
                    {
                        currentLocation = rest[..qualStart].Trim();
                        ParseQualifierString(rest[qualStart..], qualifiers = new Dictionary<string, string>());
                    }
                    else
                    {
                        currentLocation = rest.Trim();
                        qualifiers = new Dictionary<string, string>();
                    }
                }
            }
            else
            {
                // Qualifier continuation
                ParseQualifierString(trimmed, qualifiers);
            }
        }

        // Save last feature
        if (!string.IsNullOrEmpty(currentKey))
        {
            features.Add(CreateFeature(currentKey, currentLocation, qualifiers));
        }

        return features;
    }

    private static void ParseQualifierString(string text, Dictionary<string, string> qualifiers)
    {
        var matches = QualifierRegex().Matches(text);
        foreach (Match match in matches)
        {
            var name = match.Groups[1].Value;
            var value = match.Groups[2].Success ? match.Groups[2].Value.Trim('"') : "true";
            qualifiers[name] = value;
        }
    }

    private static Feature CreateFeature(string key, string rawLocation, Dictionary<string, string> qualifiers)
    {
        var location = ParseLocation(rawLocation);
        return new Feature(key, location, qualifiers);
    }

    /// <summary>
    /// Parses a feature location string.
    /// </summary>
    public static Location ParseLocation(string locationStr)
    {
        if (string.IsNullOrEmpty(locationStr))
            return new Location(0, 0, false, false, Array.Empty<(int, int)>(), locationStr);

        var (start, end, isComplement, isJoin, parts) =
            SequenceFormatHelper.ParseLocationParts(locationStr, useStartsWithForComplement: false);
        return new Location(start, end, isComplement, isJoin, parts, locationStr);
    }

    private static string ParseSequence(List<string> lines)
    {
        var sb = new StringBuilder();
        bool inSequence = false;

        foreach (var line in lines)
        {
            if (line.StartsWith("SQ", StringComparison.Ordinal))
            {
                inSequence = true;
                continue;
            }

            if (inSequence && !line.StartsWith("//", StringComparison.Ordinal))
            {
                foreach (var c in line)
                {
                    if (char.IsLetter(c))
                    {
                        sb.Append(char.ToUpperInvariant(c));
                    }
                }
            }
        }

        return sb.ToString();
    }

    private static bool IsStandardPrefix(string prefix)
    {
        return prefix is ID or AC or SV or DT or DE or KW or OS or OC or OG or
               RN or RC or RP or RX or RG or RA or RT or RL or DR or CC or FH or FT or SQ or XX;
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Converts an EMBL record to GenBank format.
    /// </summary>
    public static GenBankParser.GenBankRecord ToGenBank(EmblRecord embl)
    {
        var gbFeatures = embl.Features.Select(f =>
            new GenBankParser.Feature(
                f.Key,
                new GenBankParser.Location(f.Location.Start, f.Location.End, f.Location.IsComplement,
                    f.Location.IsJoin, f.Location.Parts, f.Location.RawLocation),
                f.Qualifiers
            )).ToList();

        var gbReferences = embl.References.Select(r =>
            new GenBankParser.Reference(r.Number, r.Authors, r.Title, r.Journal, r.CrossReference, null, null)
        ).ToList();

        return new GenBankParser.GenBankRecord(
            embl.Accession,
            embl.SequenceLength,
            embl.MoleculeType,
            embl.Topology,
            embl.TaxonomicDivision,
            null,
            embl.Description,
            embl.Accession,
            embl.SequenceVersion,
            embl.Keywords,
            embl.Organism,
            string.Join("; ", embl.OrganismClassification),
            gbReferences,
            gbFeatures,
            embl.Sequence,
            embl.AdditionalFields
        );
    }

    /// <summary>
    /// Extracts features by type.
    /// </summary>
    public static IEnumerable<Feature> GetFeatures(EmblRecord record, string featureKey)
    {
        return record.Features.Where(f =>
            f.Key.Equals(featureKey, StringComparison.OrdinalIgnoreCase));
    }

    /// <summary>
    /// Gets CDS features.
    /// </summary>
    public static IEnumerable<Feature> GetCDS(EmblRecord record)
    {
        return GetFeatures(record, "CDS");
    }

    /// <summary>
    /// Gets gene features.
    /// </summary>
    public static IEnumerable<Feature> GetGenes(EmblRecord record)
    {
        return GetFeatures(record, "gene");
    }

    /// <summary>
    /// Extracts subsequence based on location.
    /// </summary>
    public static string ExtractSequence(EmblRecord record, Location location)
        => FeatureLocationHelper.ExtractSequence(record.Sequence, location);

    #endregion

    #region Regex Patterns

    [GeneratedRegex(@"(\d+)\s*BP", RegexOptions.IgnoreCase)]
    private static partial Regex LengthRegex();

    [GeneratedRegex(@"\[(\d+)\]")]
    private static partial Regex ReferenceNumberRegex();

    [GeneratedRegex(@"/(\w+)(?:=([^/]+))?")]
    private static partial Regex QualifierRegex();

    #endregion
}

