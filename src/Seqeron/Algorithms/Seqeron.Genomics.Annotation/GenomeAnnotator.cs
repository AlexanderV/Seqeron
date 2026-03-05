using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Provides genome annotation algorithms including ORF finding and gene prediction.
/// </summary>
public static class GenomeAnnotator
{
    /// <summary>
    /// Represents an open reading frame (ORF).
    /// </summary>
    public readonly record struct OpenReadingFrame(
        int Start,
        int End,
        int Frame,
        bool IsReverseComplement,
        string Sequence,
        string ProteinSequence);

    /// <summary>
    /// Represents a gene annotation.
    /// </summary>
    public readonly record struct GeneAnnotation(
        string GeneId,
        int Start,
        int End,
        char Strand,
        string Type,
        string Product,
        IReadOnlyDictionary<string, string> Attributes);

    /// <summary>
    /// Represents a genomic feature.
    /// </summary>
    public readonly record struct GenomicFeature(
        string FeatureId,
        string Type,
        int Start,
        int End,
        char Strand,
        double? Score,
        int? Phase,
        IReadOnlyDictionary<string, string> Attributes);

    /// <summary>
    /// Start codons to look for.
    /// </summary>
    private static readonly HashSet<string> StartCodons = new(StringComparer.OrdinalIgnoreCase)
    {
        "ATG", "GTG", "TTG"
    };

    /// <summary>
    /// Stop codons.
    /// </summary>
    private static readonly HashSet<string> StopCodons = new(StringComparer.OrdinalIgnoreCase)
    {
        "TAA", "TAG", "TGA"
    };

    /// <summary>
    /// Complement mapping for nucleotides.
    /// </summary>
    private static readonly Dictionary<char, char> ComplementMap = new()
    {
        ['A'] = 'T',
        ['T'] = 'A',
        ['C'] = 'G',
        ['G'] = 'C',
        ['a'] = 't',
        ['t'] = 'a',
        ['c'] = 'g',
        ['g'] = 'c',
        ['N'] = 'N',
        ['n'] = 'n'
    };

    /// <summary>
    /// Finds all open reading frames in a DNA sequence.
    /// </summary>
    /// <param name="dnaSequence">The DNA sequence to search.</param>
    /// <param name="minLength">Minimum ORF length in amino acids.</param>
    /// <param name="searchBothStrands">Whether to search reverse complement.</param>
    /// <param name="requireStartCodon">Whether to require ATG start.</param>
    /// <returns>All found ORFs.</returns>
    public static IEnumerable<OpenReadingFrame> FindOrfs(
        string dnaSequence,
        int minLength = 100,
        bool searchBothStrands = true,
        bool requireStartCodon = true)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            yield break;

        var geneticCode = GeneticCode.Standard;

        // Search forward strand
        foreach (var orf in FindOrfsInStrand(dnaSequence, geneticCode, minLength, requireStartCodon, false))
        {
            yield return orf;
        }

        // Search reverse complement
        if (searchBothStrands)
        {
            string revComp = DnaSequence.GetReverseComplementString(dnaSequence);
            foreach (var orf in FindOrfsInStrand(revComp, geneticCode, minLength, requireStartCodon, true))
            {
                // Adjust coordinates for reverse strand
                int adjStart = dnaSequence.Length - orf.End;
                int adjEnd = dnaSequence.Length - orf.Start;
                yield return orf with { Start = adjStart, End = adjEnd };
            }
        }
    }

    private static IEnumerable<OpenReadingFrame> FindOrfsInStrand(
        string sequence,
        GeneticCode geneticCode,
        int minLength,
        bool requireStartCodon,
        bool isReverseComplement)
    {
        // Search all three reading frames
        for (int frame = 0; frame < 3; frame++)
        {
            foreach (var orf in FindOrfsInFrame(sequence, frame, geneticCode, minLength, requireStartCodon, isReverseComplement))
            {
                yield return orf;
            }
        }
    }

    private static IEnumerable<OpenReadingFrame> FindOrfsInFrame(
        string sequence,
        int frame,
        GeneticCode geneticCode,
        int minLength,
        bool requireStartCodon,
        bool isReverseComplement)
    {
        var currentOrfStarts = new List<int>();

        for (int i = frame; i <= sequence.Length - 3; i += 3)
        {
            string codon = sequence.Substring(i, 3).ToUpperInvariant();

            if (StartCodons.Contains(codon))
            {
                currentOrfStarts.Add(i);
            }

            if (StopCodons.Contains(codon))
            {
                // Complete all current ORFs
                foreach (int start in currentOrfStarts)
                {
                    int orfNucleotideLength = i - start;
                    int orfAaLength = orfNucleotideLength / 3;

                    if (orfAaLength >= minLength)
                    {
                        string orfSeq = sequence.Substring(start, i + 3 - start);
                        string protein = Translator.Translate(orfSeq, geneticCode).Sequence;

                        yield return new OpenReadingFrame(
                            Start: start,
                            End: i + 3,
                            Frame: frame + 1,
                            IsReverseComplement: isReverseComplement,
                            Sequence: orfSeq,
                            ProteinSequence: protein);
                    }
                }

                currentOrfStarts.Clear();

                if (!requireStartCodon)
                {
                    currentOrfStarts.Add(i + 3);
                }
            }
        }

        // Handle ORFs that extend to end of sequence (no stop codon found)
        if (!requireStartCodon && currentOrfStarts.Count > 0)
        {
            foreach (int start in currentOrfStarts)
            {
                int endPos = sequence.Length - ((sequence.Length - start) % 3);
                int orfAaLength = (endPos - start) / 3;

                if (orfAaLength >= minLength)
                {
                    string orfSeq = sequence.Substring(start, endPos - start);
                    string protein = Translator.Translate(orfSeq, geneticCode).Sequence;

                    yield return new OpenReadingFrame(
                        Start: start,
                        End: endPos,
                        Frame: frame + 1,
                        IsReverseComplement: isReverseComplement,
                        Sequence: orfSeq,
                        ProteinSequence: protein);
                }
            }
        }
    }

    /// <summary>
    /// Finds the longest ORF in each reading frame.
    /// </summary>
    public static IReadOnlyDictionary<int, OpenReadingFrame?> FindLongestOrfsPerFrame(
        string dnaSequence,
        bool searchBothStrands = true)
    {
        var orfs = FindOrfs(dnaSequence, minLength: 1, searchBothStrands, requireStartCodon: true).ToList();

        var result = new Dictionary<int, OpenReadingFrame?>();

        for (int frame = 1; frame <= 3; frame++)
        {
            var frameOrfs = orfs.Where(o => o.Frame == frame && !o.IsReverseComplement);
            result[frame] = frameOrfs.OrderByDescending(o => o.ProteinSequence.Length).FirstOrDefault();
        }

        if (searchBothStrands)
        {
            for (int frame = 1; frame <= 3; frame++)
            {
                var frameOrfs = orfs.Where(o => o.Frame == frame && o.IsReverseComplement);
                result[-frame] = frameOrfs.OrderByDescending(o => o.ProteinSequence.Length).FirstOrDefault();
            }
        }

        return result;
    }

    /// <summary>
    /// Finds potential Shine-Dalgarno (ribosome binding site) sequences.
    /// </summary>
    public static IEnumerable<(int position, string sequence, double score)> FindRibosomeBindingSites(
        string dnaSequence,
        int upstreamWindow = 20,
        int minDistance = 4,
        int maxDistance = 15)
    {
        // Consensus Shine-Dalgarno: AGGAGG (binds to 3' end of 16S rRNA)
        string[] sdMotifs = { "AGGAGG", "GGAGG", "AGGAG", "GAGG", "AGGA" };

        var orfs = FindOrfs(dnaSequence, minLength: 30).ToList();

        foreach (var orf in orfs.Where(o => !o.IsReverseComplement))
        {
            int searchStart = Math.Max(0, orf.Start - upstreamWindow);
            int searchEnd = orf.Start - minDistance;

            if (searchEnd <= searchStart) continue;

            string upstream = dnaSequence.Substring(searchStart, searchEnd - searchStart).ToUpperInvariant();

            foreach (string motif in sdMotifs)
            {
                int pos = upstream.IndexOf(motif, StringComparison.Ordinal);
                while (pos >= 0)
                {
                    int genomicPos = searchStart + pos;
                    int distanceToStart = orf.Start - genomicPos - motif.Length;

                    if (distanceToStart >= minDistance && distanceToStart <= maxDistance)
                    {
                        double score = (double)motif.Length / 6.0; // Normalize to consensus length
                        yield return (genomicPos, motif, score);
                    }

                    pos = upstream.IndexOf(motif, pos + 1, StringComparison.Ordinal);
                }
            }
        }
    }

    /// <summary>
    /// Predicts genes using a simple ORF-based approach.
    /// </summary>
    public static IEnumerable<GeneAnnotation> PredictGenes(
        string dnaSequence,
        int minOrfLength = 100,
        string prefix = "gene")
    {
        var orfs = FindOrfs(dnaSequence, minOrfLength, searchBothStrands: true, requireStartCodon: true)
            .OrderBy(o => o.Start)
            .ToList();

        int geneCount = 0;

        foreach (var orf in orfs)
        {
            geneCount++;
            char strand = orf.IsReverseComplement ? '-' : '+';

            var attributes = new Dictionary<string, string>
            {
                ["frame"] = orf.Frame.ToString(),
                ["protein_length"] = orf.ProteinSequence.TrimEnd('*').Length.ToString(),
                ["translation"] = orf.ProteinSequence
            };

            yield return new GeneAnnotation(
                GeneId: $"{prefix}_{geneCount:D4}",
                Start: orf.Start,
                End: orf.End,
                Strand: strand,
                Type: "CDS",
                Product: "hypothetical protein",
                Attributes: attributes);
        }
    }

    /// <summary>
    /// Parses a GFF3 format annotation file.
    /// </summary>
    public static IEnumerable<GenomicFeature> ParseGff3(IEnumerable<string> lines)
    {
        int featureCount = 0;

        foreach (string line in lines)
        {
            if (string.IsNullOrWhiteSpace(line) || line.StartsWith("#"))
                continue;

            string[] parts = line.Split('\t');
            if (parts.Length < 9)
                continue;

            featureCount++;

            int start = int.Parse(parts[3]);
            int end = int.Parse(parts[4]);
            double? score = parts[5] == "." ? null : double.Parse(parts[5]);
            char strand = parts[6][0];
            int? phase = parts[7] == "." ? null : int.Parse(parts[7]);

            var attributes = ParseGff3Attributes(parts[8]);

            string featureId = attributes.GetValueOrDefault("ID", $"feature_{featureCount}");

            yield return new GenomicFeature(
                FeatureId: featureId,
                Type: parts[2],
                Start: start,
                End: end,
                Strand: strand,
                Score: score,
                Phase: phase,
                Attributes: attributes);
        }
    }

    private static Dictionary<string, string> ParseGff3Attributes(string attributeString)
    {
        var attributes = new Dictionary<string, string>();

        foreach (string attr in attributeString.Split(';', StringSplitOptions.RemoveEmptyEntries))
        {
            int eqPos = attr.IndexOf('=');
            if (eqPos > 0)
            {
                string key = attr.Substring(0, eqPos).Trim();
                string value = Uri.UnescapeDataString(attr.Substring(eqPos + 1).Trim());
                attributes[key] = value;
            }
        }

        return attributes;
    }

    /// <summary>
    /// Exports annotations in GFF3 format.
    /// Phase: "0" for CDS features, "." for all others.
    /// Source: GFF3 Spec v1.26 — NOTE 4: "The phase is REQUIRED for all CDS features."
    /// Encoding: GFF3-specific (only tab, newline, CR, %, control chars, ;, =, &amp;, , in column 9).
    /// Source: GFF3 Spec v1.26 — "no other characters may be encoded."
    /// </summary>
    public static IEnumerable<string> ToGff3(
        IEnumerable<GeneAnnotation> annotations,
        string seqId = "seq1")
    {
        yield return "##gff-version 3";

        foreach (var ann in annotations)
        {
            string attributes = FormatGff3Attributes(ann.Attributes, ann.GeneId, ann.Product);
            string phase = string.Equals(ann.Type, "CDS", StringComparison.OrdinalIgnoreCase) ? "0" : ".";

            yield return $"{seqId}\t.\t{ann.Type}\t{ann.Start + 1}\t{ann.End}\t.\t{ann.Strand}\t{phase}\t{attributes}";
        }
    }

    private static string FormatGff3Attributes(
        IReadOnlyDictionary<string, string> attributes,
        string id,
        string product)
    {
        var sb = new StringBuilder();
        sb.Append($"ID={EncodeGff3Value(id)}");
        sb.Append($";product={EncodeGff3Value(product)}");

        foreach (var kvp in attributes)
        {
            if (kvp.Key != "translation") // Skip large translation
            {
                sb.Append($";{kvp.Key}={EncodeGff3Value(kvp.Value)}");
            }
        }

        return sb.ToString();
    }

    /// <summary>
    /// Encodes a value for GFF3 column 9 (attributes) per the GFF3 Specification v1.26.
    /// Only encodes characters required by the spec; "no other characters may be encoded."
    /// Source: GFF3 Spec — tab (%09), newline (%0A), CR (%0D), percent (%25), control chars,
    /// plus column 9 reserved: semicolon (%3B), equals (%3D), ampersand (%26), comma (%2C).
    /// Spaces are explicitly allowed unencoded: "unescaped spaces are allowed within fields."
    /// </summary>
    private static string EncodeGff3Value(string value)
    {
        var sb = new StringBuilder(value.Length);
        foreach (char c in value)
        {
            switch (c)
            {
                case '\t': sb.Append("%09"); break;
                case '\n': sb.Append("%0A"); break;
                case '\r': sb.Append("%0D"); break;
                case '%':  sb.Append("%25"); break;
                case ';':  sb.Append("%3B"); break;
                case '=':  sb.Append("%3D"); break;
                case '&':  sb.Append("%26"); break;
                case ',':  sb.Append("%2C"); break;
                default:
                    if (c < 0x20 || c == 0x7F)
                        sb.Append($"%{(int)c:X2}");
                    else
                        sb.Append(c);
                    break;
            }
        }
        return sb.ToString();
    }

    /// <summary>
    /// Finds promoter motifs (-10 and -35 boxes in bacterial promoters).
    /// Scoring uses E. coli position-specific nucleotide occurrence probabilities.
    /// Source: Wikipedia "Promoter (genetics)" / Harley &amp; Reynolds (1987) NAR 15(5):2343-2361.
    /// </summary>
    public static IEnumerable<(int position, string type, string sequence, double score)> FindPromoterMotifs(
        string dnaSequence)
    {
        // Probability-weighted scoring from E. coli nucleotide occurrence data.
        // Source: Wikipedia "Promoter (genetics)" — Probability of occurrence of each nucleotide
        // Reference: Harley & Reynolds (1987) NAR 15(5):2343-2361
        //
        // -35 box consensus TTGACA: T(69%) T(79%) G(61%) A(56%) C(54%) A(54%) → total 3.73
        // -10 box consensus TATAAT: T(77%) A(76%) T(60%) A(61%) A(56%) T(82%) → total 4.12
        //
        // Score = sum(matched position probabilities) / sum(all 6 consensus probabilities)
        // Variants: full consensus, prefix-5bp, suffix-5bp, prefix-4bp.

        // -35 box variants with probability-weighted scores
        (string motif, double score)[] minus35Motifs =
        {
            ("TTGACA", 1.000),  // Pos 1–6: 3.73 / 3.73
            ("TTGAC",  0.855),  // Pos 1–5: (0.69+0.79+0.61+0.56+0.54) / 3.73
            ("TGACA",  0.815),  // Pos 2–6: (0.79+0.61+0.56+0.54+0.54) / 3.73
            ("TTGA",   0.710),  // Pos 1–4: (0.69+0.79+0.61+0.56) / 3.73
        };

        // -10 box (Pribnow box) variants with probability-weighted scores
        (string motif, double score)[] minus10Motifs =
        {
            ("TATAAT", 1.000),  // Pos 1–6: 4.12 / 4.12
            ("TATAA",  0.801),  // Pos 1–5: (0.77+0.76+0.60+0.61+0.56) / 4.12
            ("ATAAT",  0.813),  // Pos 2–6: (0.76+0.60+0.61+0.56+0.82) / 4.12
            ("TATA",   0.665),  // Pos 1–4: (0.77+0.76+0.60+0.61) / 4.12
        };

        string seq = dnaSequence.ToUpperInvariant();

        foreach (var (motif, score) in minus35Motifs)
        {
            for (int i = 0; i <= seq.Length - motif.Length; i++)
            {
                if (seq.Substring(i, motif.Length) == motif)
                {
                    yield return (i, "-35 box", motif, score);
                }
            }
        }

        foreach (var (motif, score) in minus10Motifs)
        {
            for (int i = 0; i <= seq.Length - motif.Length; i++)
            {
                if (seq.Substring(i, motif.Length) == motif)
                {
                    yield return (i, "-10 box", motif, score);
                }
            }
        }
    }

    /// <summary>
    /// Calculates coding potential using hexamer frequency bias.
    /// </summary>
    public static double CalculateCodingPotential(string sequence)
    {
        if (sequence.Length < 6) return 0;

        // Simplified approach: measure in-frame vs out-of-frame hexamer usage
        // In real implementation, would use pre-computed coding/non-coding hexamer tables

        sequence = sequence.ToUpperInvariant();

        // Count in-frame vs out-of-frame codon usage patterns
        int validCodons = 0;
        int totalCodons = 0;

        for (int i = 0; i <= sequence.Length - 3; i += 3)
        {
            string codon = sequence.Substring(i, 3);
            if (codon.All(c => "ACGT".Contains(c)))
            {
                totalCodons++;
                // Stop codons indicate non-coding in middle
                if (!StopCodons.Contains(codon) || i == sequence.Length - 3)
                {
                    validCodons++;
                }
            }
        }

        if (totalCodons == 0) return 0;

        // Additional factors: GC content in third position (wobble)
        int gc3 = 0;
        for (int i = 2; i < sequence.Length; i += 3)
        {
            if (sequence[i] == 'G' || sequence[i] == 'C')
                gc3++;
        }

        double gc3ratio = totalCodons > 0 ? (double)gc3 / totalCodons : 0;
        double validRatio = (double)validCodons / totalCodons;

        // Combine factors (simplified scoring)
        return validRatio * 0.7 + Math.Abs(gc3ratio - 0.5) * 0.6;
    }

    /// <summary>
    /// Finds repetitive elements in the sequence.
    /// </summary>
    public static IEnumerable<(int start, int end, string type, string sequence)> FindRepetitiveElements(
        string dnaSequence,
        int minRepeatLength = 10,
        int minCopies = 2)
    {
        // Find tandem repeats
        foreach (var repeat in FindTandemRepeats(dnaSequence, minRepeatLength, minCopies))
        {
            yield return repeat;
        }

        // Find inverted repeats
        foreach (var repeat in FindInvertedRepeats(dnaSequence, minRepeatLength))
        {
            yield return repeat;
        }
    }

    private static IEnumerable<(int start, int end, string type, string sequence)> FindTandemRepeats(
        string sequence,
        int minLength,
        int minCopies)
    {
        sequence = sequence.ToUpperInvariant();

        for (int unitLen = 1; unitLen <= sequence.Length / minCopies && unitLen <= 50; unitLen++)
        {
            for (int start = 0; start <= sequence.Length - unitLen * minCopies; start++)
            {
                string unit = sequence.Substring(start, unitLen);
                int copies = 1;
                int pos = start + unitLen;

                while (pos + unitLen <= sequence.Length &&
                       sequence.Substring(pos, unitLen) == unit)
                {
                    copies++;
                    pos += unitLen;
                }

                if (copies >= minCopies && unitLen * copies >= minLength)
                {
                    int end = start + unitLen * copies;
                    string repeatSeq = sequence.Substring(start, end - start);
                    yield return (start, end, "tandem_repeat", repeatSeq);

                    // Skip past this repeat
                    start = end - 1;
                }
            }
        }
    }

    private static IEnumerable<(int start, int end, string type, string sequence)> FindInvertedRepeats(
        string sequence,
        int minLength)
    {
        sequence = sequence.ToUpperInvariant();

        for (int i = 0; i <= sequence.Length - minLength * 2; i++)
        {
            for (int len = minLength; len <= (sequence.Length - i) / 2 && len <= 100; len++)
            {
                string arm1 = sequence.Substring(i, len);
                string arm1RevComp = DnaSequence.GetReverseComplementString(arm1);

                // Look for reverse complement downstream
                for (int j = i + len; j <= sequence.Length - len; j++)
                {
                    string arm2 = sequence.Substring(j, len);

                    if (arm2 == arm1RevComp)
                    {
                        yield return (i, j + len, "inverted_repeat", arm1 + "..." + arm2);
                    }
                }
            }
        }
    }

    /// <summary>
    /// Calculates codon usage statistics for a sequence.
    /// </summary>
    public static IReadOnlyDictionary<string, int> GetCodonUsage(string dnaSequence)
    {
        var usage = new Dictionary<string, int>();
        string seq = dnaSequence.ToUpperInvariant();

        for (int i = 0; i <= seq.Length - 3; i += 3)
        {
            string codon = seq.Substring(i, 3);
            if (codon.All(c => "ACGT".Contains(c)))
            {
                usage[codon] = usage.GetValueOrDefault(codon, 0) + 1;
            }
        }

        return usage;
    }
}
