using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace SuffixTree.Genomics;

/// <summary>
/// Parser for GenBank flat file format (.gb, .gbk, .genbank).
/// Supports parsing of sequence records with annotations, features, and references.
/// </summary>
public static partial class GenBankParser
{
    #region Records

    /// <summary>Represents a complete GenBank record</summary>
    public readonly record struct GenBankRecord(
        string Locus,
        int SequenceLength,
        string MoleculeType,
        string Topology,
        string Division,
        DateTime? Date,
        string Definition,
        string Accession,
        string Version,
        IReadOnlyList<string> Keywords,
        string Organism,
        string Taxonomy,
        IReadOnlyList<Reference> References,
        IReadOnlyList<Feature> Features,
        string Sequence,
        IReadOnlyDictionary<string, string> AdditionalFields);

    /// <summary>Literature reference</summary>
    public readonly record struct Reference(
        int Number,
        string Authors,
        string Title,
        string Journal,
        string PubMed,
        int? BaseFrom,
        int? BaseTo);

    /// <summary>Sequence feature with location and qualifiers</summary>
    public readonly record struct Feature(
        string Key,
        Location Location,
        IReadOnlyDictionary<string, string> Qualifiers);

    /// <summary>Feature location (supports joins, complements, etc.)</summary>
    public readonly record struct Location(
        int Start,
        int End,
        bool IsComplement,
        bool IsJoin,
        IReadOnlyList<(int Start, int End)> Parts,
        string RawLocation);

    #endregion

    #region Main Parsing Methods

    /// <summary>
    /// Parses GenBank records from a file.
    /// </summary>
    public static IEnumerable<GenBankRecord> ParseFile(string filePath)
    {
        if (string.IsNullOrEmpty(filePath))
            yield break;

        if (!File.Exists(filePath))
            yield break;

        foreach (var record in Parse(File.ReadAllText(filePath)))
        {
            yield return record;
        }
    }

    /// <summary>
    /// Parses GenBank records from text content.
    /// Multiple records are separated by // delimiter.
    /// </summary>
    public static IEnumerable<GenBankRecord> Parse(string content)
    {
        if (string.IsNullOrEmpty(content))
            yield break;

        // Split by record delimiter
        var recordTexts = content.Split(new[] { "\n//" }, StringSplitOptions.RemoveEmptyEntries);

        foreach (var recordText in recordTexts)
        {
            var trimmed = recordText.Trim();
            if (!string.IsNullOrEmpty(trimmed) && trimmed.StartsWith("LOCUS", StringComparison.Ordinal))
            {
                var record = ParseRecord(trimmed);
                if (record.HasValue)
                    yield return record.Value;
            }
        }
    }

    /// <summary>
    /// Parses a single GenBank record.
    /// </summary>
    private static GenBankRecord? ParseRecord(string text)
    {
        var lines = text.Split('\n').Select(l => l.TrimEnd('\r')).ToList();

        // Parse LOCUS line
        var locusInfo = ParseLocusLine(lines.FirstOrDefault(l => l.StartsWith("LOCUS", StringComparison.Ordinal)) ?? "");

        // Collect sections
        var sections = ExtractSections(lines);

        var definition = sections.GetValueOrDefault("DEFINITION", "").Trim();
        var accession = sections.GetValueOrDefault("ACCESSION", "").Trim().Split(' ')[0];
        var version = sections.GetValueOrDefault("VERSION", "").Trim();

        // Keywords
        var keywordsText = sections.GetValueOrDefault("KEYWORDS", "");
        var keywords = ParseKeywords(keywordsText);

        // Organism and taxonomy from SOURCE section
        var sourceSection = sections.GetValueOrDefault("SOURCE", "");
        var (organism, taxonomy) = ParseSource(sourceSection);

        // References
        var references = ParseReferences(sections.GetValueOrDefault("REFERENCE", ""));

        // Features
        var features = ParseFeatures(sections.GetValueOrDefault("FEATURES", ""));

        // Sequence
        var sequence = ParseSequence(sections.GetValueOrDefault("ORIGIN", ""));

        // Additional fields
        var additionalFields = new Dictionary<string, string>();
        foreach (var (key, value) in sections)
        {
            if (!IsStandardField(key))
            {
                additionalFields[key] = value.Trim();
            }
        }

        return new GenBankRecord(
            locusInfo.Locus,
            locusInfo.Length,
            locusInfo.MoleculeType,
            locusInfo.Topology,
            locusInfo.Division,
            locusInfo.Date,
            definition,
            accession,
            version,
            keywords,
            organism,
            taxonomy,
            references,
            features,
            sequence,
            additionalFields
        );
    }

    #endregion

    #region Section Parsing

    private static Dictionary<string, string> ExtractSections(List<string> lines)
    {
        var sections = new Dictionary<string, string>();
        string currentSection = "";
        var currentContent = new StringBuilder();

        foreach (var line in lines)
        {
            if (line.Length >= 12 && !char.IsWhiteSpace(line[0]))
            {
                // New section
                if (!string.IsNullOrEmpty(currentSection))
                {
                    sections[currentSection] = currentContent.ToString();
                }

                var match = SectionHeaderRegex().Match(line);
                if (match.Success)
                {
                    currentSection = match.Groups[1].Value;
                    currentContent = new StringBuilder(line.Length > 12 ? line[12..] : "");
                }
                else
                {
                    currentSection = line.Trim();
                    currentContent = new StringBuilder();
                }
            }
            else if (!string.IsNullOrEmpty(currentSection))
            {
                // Continuation of current section
                currentContent.AppendLine();
                currentContent.Append(line.TrimStart());
            }
        }

        if (!string.IsNullOrEmpty(currentSection))
        {
            sections[currentSection] = currentContent.ToString();
        }

        return sections;
    }

    private static (string Locus, int Length, string MoleculeType, string Topology, string Division, DateTime? Date)
        ParseLocusLine(string line)
    {
        if (string.IsNullOrEmpty(line))
            return ("", 0, "", "", "", null);

        // LOCUS format: LOCUS name length bp type topology division date
        var parts = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

        string locus = parts.Length > 1 ? parts[1] : "";
        int length = 0;
        string moleculeType = "";
        string topology = "";
        string division = "";
        DateTime? date = null;

        if (parts.Length > 2 && int.TryParse(parts[2], out var len))
        {
            length = len;
        }

        // Find molecule type (DNA, RNA, etc.)
        foreach (var part in parts.Skip(3))
        {
            if (part is "DNA" or "RNA" or "mRNA" or "rRNA" or "tRNA" or "ss-DNA" or "ds-DNA")
            {
                moleculeType = part;
            }
            else if (part is "linear" or "circular")
            {
                topology = part;
            }
            else if (part.Length == 3 && char.IsUpper(part[0]))
            {
                division = part;
            }
            else if (DateTime.TryParseExact(part, new[] { "dd-MMM-yyyy", "dd-MMM-yy" },
                     CultureInfo.InvariantCulture, DateTimeStyles.None, out var d))
            {
                date = d;
            }
        }

        return (locus, length, moleculeType, topology, division, date);
    }

    private static IReadOnlyList<string> ParseKeywords(string text)
    {
        if (string.IsNullOrWhiteSpace(text) || text.Trim() == ".")
            return Array.Empty<string>();

        return text
            .Replace("\n", " ")
            .Split(new[] { ';', ',' }, StringSplitOptions.RemoveEmptyEntries)
            .Select(k => k.Trim().TrimEnd('.'))
            .Where(k => !string.IsNullOrEmpty(k))
            .ToList();
    }

    private static (string Organism, string Taxonomy) ParseSource(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return ("", "");

        var lines = text.Split('\n');
        var organism = "";
        var taxonomy = new StringBuilder();
        var inTaxonomy = false;

        foreach (var line in lines)
        {
            var trimmed = line.Trim();
            if (trimmed.StartsWith("ORGANISM", StringComparison.Ordinal))
            {
                organism = trimmed.Length > 9 ? trimmed[9..].Trim() : "";
                inTaxonomy = true;
            }
            else if (inTaxonomy && !string.IsNullOrEmpty(trimmed))
            {
                if (taxonomy.Length > 0) taxonomy.Append(' ');
                taxonomy.Append(trimmed);
            }
        }

        return (organism, taxonomy.ToString().TrimEnd('.'));
    }

    private static IReadOnlyList<Reference> ParseReferences(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return Array.Empty<Reference>();

        var references = new List<Reference>();
        var refBlocks = text.Split(new[] { "\nREFERENCE" }, StringSplitOptions.None);

        // Handle first block specially
        var blocks = new List<string> { refBlocks[0] };
        blocks.AddRange(refBlocks.Skip(1).Select(b => "REFERENCE" + b));

        foreach (var block in blocks.Where(b => !string.IsNullOrWhiteSpace(b)))
        {
            var refNum = 0;
            string authors = "", title = "", journal = "", pubmed = "";
            int? baseFrom = null, baseTo = null;

            var lines = block.Split('\n');
            string currentField = "";
            var currentValue = new StringBuilder();

            foreach (var line in lines)
            {
                var trimmed = line.Trim();
                if (line.Length >= 12 && !char.IsWhiteSpace(line[0]) && line.Length > 2)
                {
                    // Save previous field
                    SaveRefField(ref authors, ref title, ref journal, ref pubmed, currentField, currentValue.ToString());

                    var spaceIdx = trimmed.IndexOf(' ');
                    currentField = spaceIdx > 0 ? trimmed[..spaceIdx] : trimmed;
                    currentValue = new StringBuilder(spaceIdx > 0 ? trimmed[(spaceIdx + 1)..] : "");

                    // Parse reference number and bases
                    if (currentField == "REFERENCE" || currentField.StartsWith("REFERENCE", StringComparison.Ordinal))
                    {
                        var refMatch = ReferenceNumberRegex().Match(currentValue.ToString());
                        if (refMatch.Success)
                        {
                            refNum = int.Parse(refMatch.Groups[1].Value, CultureInfo.InvariantCulture);
                            if (refMatch.Groups[2].Success)
                                baseFrom = int.Parse(refMatch.Groups[2].Value, CultureInfo.InvariantCulture);
                            if (refMatch.Groups[3].Success)
                                baseTo = int.Parse(refMatch.Groups[3].Value, CultureInfo.InvariantCulture);
                        }
                    }
                }
                else if (!string.IsNullOrEmpty(currentField))
                {
                    currentValue.Append(' ');
                    currentValue.Append(trimmed);
                }
            }

            SaveRefField(ref authors, ref title, ref journal, ref pubmed, currentField, currentValue.ToString());

            if (refNum > 0 || !string.IsNullOrEmpty(title))
            {
                references.Add(new Reference(refNum, authors, title, journal, pubmed, baseFrom, baseTo));
            }
        }

        return references;
    }

    private static void SaveRefField(ref string authors, ref string title, ref string journal,
        ref string pubmed, string field, string value)
    {
        switch (field)
        {
            case "AUTHORS": authors = value.Trim(); break;
            case "TITLE": title = value.Trim(); break;
            case "JOURNAL": journal = value.Trim(); break;
            case "PUBMED": pubmed = value.Trim(); break;
        }
    }

    private static IReadOnlyList<Feature> ParseFeatures(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return Array.Empty<Feature>();

        var features = new List<Feature>();
        var lines = text.Split('\n').Where(l => l.Length > 0).ToList();

        string currentKey = "";
        string currentLocation = "";
        var qualifiers = new Dictionary<string, string>();

        foreach (var line in lines)
        {
            // Skip header line
            if (line.TrimStart().StartsWith("Location/Qualifiers", StringComparison.OrdinalIgnoreCase))
                continue;

            // New feature (starts at column 5)
            if (line.Length > 5 && !char.IsWhiteSpace(line[5]) && char.IsWhiteSpace(line[0]))
            {
                // Save previous feature
                if (!string.IsNullOrEmpty(currentKey))
                {
                    features.Add(CreateFeature(currentKey, currentLocation, qualifiers));
                }

                var trimmed = line.Trim();
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
                qualifiers = new Dictionary<string, string>();
            }
            // Qualifier (starts at column 21, with /)
            else if (line.Length > 21 && line.TrimStart().StartsWith('/'))
            {
                var qualLine = line.Trim()[1..]; // Remove leading /
                var eqIdx = qualLine.IndexOf('=');
                if (eqIdx > 0)
                {
                    var qualName = qualLine[..eqIdx];
                    var qualValue = qualLine[(eqIdx + 1)..].Trim('"');
                    qualifiers[qualName] = qualValue;
                }
                else
                {
                    qualifiers[qualLine] = "true";
                }
            }
            // Continuation of qualifier value
            else if (!string.IsNullOrEmpty(currentKey) && line.Length > 21)
            {
                var lastQual = qualifiers.Keys.LastOrDefault();
                if (lastQual != null)
                {
                    qualifiers[lastQual] += " " + line.Trim().Trim('"');
                }
            }
        }

        // Save last feature
        if (!string.IsNullOrEmpty(currentKey))
        {
            features.Add(CreateFeature(currentKey, currentLocation, qualifiers));
        }

        return features;
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

        bool isComplement = locationStr.StartsWith("complement(", StringComparison.OrdinalIgnoreCase);
        bool isJoin = locationStr.Contains("join(", StringComparison.OrdinalIgnoreCase);

        var parts = new List<(int Start, int End)>();

        // Extract ranges using regex
        var rangeMatches = LocationRangeRegex().Matches(locationStr);
        foreach (Match match in rangeMatches)
        {
            int start = int.Parse(match.Groups[1].Value, CultureInfo.InvariantCulture);
            int end = match.Groups[2].Success
                ? int.Parse(match.Groups[2].Value, CultureInfo.InvariantCulture)
                : start;
            parts.Add((start, end));
        }

        int overallStart = parts.Count > 0 ? parts.Min(p => p.Start) : 0;
        int overallEnd = parts.Count > 0 ? parts.Max(p => p.End) : 0;

        return new Location(overallStart, overallEnd, isComplement, isJoin, parts, locationStr);
    }

    private static string ParseSequence(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return "";

        var sb = new StringBuilder();
        foreach (var line in text.Split('\n'))
        {
            foreach (var c in line)
            {
                if (char.IsLetter(c))
                {
                    sb.Append(char.ToUpperInvariant(c));
                }
            }
        }
        return sb.ToString();
    }

    private static bool IsStandardField(string field)
    {
        return field is "LOCUS" or "DEFINITION" or "ACCESSION" or "VERSION" or
               "KEYWORDS" or "SOURCE" or "REFERENCE" or "FEATURES" or "ORIGIN" or
               "DBLINK" or "DBSOURCE" or "COMMENT";
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Extracts a specific feature type from a GenBank record.
    /// </summary>
    public static IEnumerable<Feature> GetFeatures(GenBankRecord record, string featureKey)
    {
        return record.Features.Where(f =>
            f.Key.Equals(featureKey, StringComparison.OrdinalIgnoreCase));
    }

    /// <summary>
    /// Gets all CDS (coding sequence) features.
    /// </summary>
    public static IEnumerable<Feature> GetCDS(GenBankRecord record)
    {
        return GetFeatures(record, "CDS");
    }

    /// <summary>
    /// Gets all gene features.
    /// </summary>
    public static IEnumerable<Feature> GetGenes(GenBankRecord record)
    {
        return GetFeatures(record, "gene");
    }

    /// <summary>
    /// Extracts a subsequence based on a feature location.
    /// </summary>
    public static string ExtractSequence(GenBankRecord record, Location location)
    {
        if (string.IsNullOrEmpty(record.Sequence))
            return "";

        var sb = new StringBuilder();

        if (location.Parts.Count > 0)
        {
            foreach (var (start, end) in location.Parts)
            {
                var realStart = Math.Max(0, start - 1); // GenBank is 1-based
                var realEnd = Math.Min(record.Sequence.Length, end);
                if (realStart < realEnd)
                {
                    sb.Append(record.Sequence[realStart..realEnd]);
                }
            }
        }
        else if (location.Start > 0 && location.End > 0)
        {
            var realStart = Math.Max(0, location.Start - 1);
            var realEnd = Math.Min(record.Sequence.Length, location.End);
            if (realStart < realEnd)
            {
                sb.Append(record.Sequence[realStart..realEnd]);
            }
        }

        var seq = sb.ToString();
        return location.IsComplement ? new DnaSequence(seq).ReverseComplement().Sequence : seq;
    }

    /// <summary>
    /// Gets qualifier value from a feature.
    /// </summary>
    public static string? GetQualifier(Feature feature, string qualifierName)
    {
        return feature.Qualifiers.TryGetValue(qualifierName, out var value) ? value : null;
    }

    /// <summary>
    /// Translates a CDS feature to protein sequence.
    /// </summary>
    public static string? TranslateCDS(GenBankRecord record, Feature cds)
    {
        // Check if translation is already in qualifiers
        if (GetQualifier(cds, "translation") is { } existingTranslation)
            return existingTranslation;

        var dnaSeq = ExtractSequence(record, cds.Location);
        if (string.IsNullOrEmpty(dnaSeq))
            return null;

        // Simple translation
        return Translate(dnaSeq);
    }

    private static string Translate(string dna)
    {
        var sb = new StringBuilder();
        for (int i = 0; i + 2 < dna.Length; i += 3)
        {
            var codon = dna.Substring(i, 3);
            sb.Append(CodonToAminoAcid(codon));
        }
        return sb.ToString();
    }

    private static char CodonToAminoAcid(string codon)
    {
        return codon.ToUpperInvariant() switch
        {
            "TTT" or "TTC" => 'F',
            "TTA" or "TTG" or "CTT" or "CTC" or "CTA" or "CTG" => 'L',
            "ATT" or "ATC" or "ATA" => 'I',
            "ATG" => 'M',
            "GTT" or "GTC" or "GTA" or "GTG" => 'V',
            "TCT" or "TCC" or "TCA" or "TCG" or "AGT" or "AGC" => 'S',
            "CCT" or "CCC" or "CCA" or "CCG" => 'P',
            "ACT" or "ACC" or "ACA" or "ACG" => 'T',
            "GCT" or "GCC" or "GCA" or "GCG" => 'A',
            "TAT" or "TAC" => 'Y',
            "TAA" or "TAG" or "TGA" => '*',
            "CAT" or "CAC" => 'H',
            "CAA" or "CAG" => 'Q',
            "AAT" or "AAC" => 'N',
            "AAA" or "AAG" => 'K',
            "GAT" or "GAC" => 'D',
            "GAA" or "GAG" => 'E',
            "TGT" or "TGC" => 'C',
            "TGG" => 'W',
            "CGT" or "CGC" or "CGA" or "CGG" or "AGA" or "AGG" => 'R',
            "GGT" or "GGC" or "GGA" or "GGG" => 'G',
            _ => 'X'
        };
    }

    #endregion

    #region Regex Patterns

    [GeneratedRegex(@"^([A-Z]+)\s+")]
    private static partial Regex SectionHeaderRegex();

    [GeneratedRegex(@"(\d+)(?:\s+\(bases\s+(\d+)\s+to\s+(\d+)\))?")]
    private static partial Regex ReferenceNumberRegex();

    [GeneratedRegex(@"(\d+)(?:\.\.(\d+))?")]
    private static partial Regex LocationRangeRegex();

    #endregion
}
