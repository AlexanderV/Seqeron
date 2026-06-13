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
            double? score = parts[5] == "." ? null : double.Parse(parts[5], System.Globalization.CultureInfo.InvariantCulture);
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
                case '%': sb.Append("%25"); break;
                case ';': sb.Append("%3B"); break;
                case '=': sb.Append("%3D"); break;
                case '&': sb.Append("%26"); break;
                case ',': sb.Append("%2C"); break;
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
    /// Default k-mer (hexamer) word size used by the CPAT hexamer-usage score.
    /// </summary>
    /// <remarks>
    /// Per Wang et al. (2013), CPAT scores hexamers (6-mers); the reference
    /// implementation calls <c>kmer_ratio(seq, word_size, step_size, ...)</c>
    /// with <c>word_size = 6</c>.
    /// CPAT: Wang L et al. (2013). Nucleic Acids Res 41(6):e74. https://doi.org/10.1093/nar/gkt006
    /// FrameKmer.kmer_ratio: https://github.com/WGLab/lncScore/blob/master/tools/cpmodule/FrameKmer.py
    /// </remarks>
    private const int HexamerWordSize = 6;

    /// <summary>
    /// Default step size for the in-frame hexamer sliding window.
    /// </summary>
    /// <remarks>
    /// CPAT extracts hexamers in-frame: the sliding window advances by 3
    /// nucleotides (one codon) so successive 6-mers stay on codon boundaries
    /// (FrameKmer.word_generator with step_size = 3, frame = 0).
    /// </remarks>
    private const int HexamerStepSize = 3;

    /// <summary>
    /// Calculates the coding potential of a DNA sequence using the CPAT hexamer
    /// usage-bias log-likelihood score (Wang et al. 2013).
    /// </summary>
    /// <param name="sequence">DNA sequence to score (case-insensitive).</param>
    /// <param name="codingHexamerFrequencies">
    /// In-frame hexamer frequency table built from a coding (CDS) training set.
    /// Keys are uppercase 6-mers over A/C/G/T; values are non-negative counts or
    /// proportions.
    /// </param>
    /// <param name="noncodingHexamerFrequencies">
    /// In-frame hexamer frequency table built from a non-coding (background)
    /// training set, in the same units as <paramref name="codingHexamerFrequencies"/>.
    /// </param>
    /// <param name="wordSize">Hexamer length (default 6, per CPAT).</param>
    /// <param name="stepSize">In-frame window step (default 3, per CPAT).</param>
    /// <returns>
    /// The mean per-hexamer log-likelihood ratio. Positive values indicate a
    /// coding sequence, negative values a non-coding sequence. Returns 0 when the
    /// sequence is shorter than <paramref name="wordSize"/> or when no scorable
    /// hexamer is found.
    /// </returns>
    /// <remarks>
    /// Implements <c>FrameKmer.kmer_ratio</c> (frame 0) from the CPAT reference
    /// implementation. For each in-frame hexamer present in both tables the score
    /// adds <c>ln(coding[k] / noncoding[k])</c>; a hexamer found only in the
    /// coding table adds +1, one found only in the non-coding table subtracts 1;
    /// hexamers missing from either table (e.g. containing N) are skipped. The sum
    /// is divided by the number of scored hexamers.
    /// CPAT: Wang L et al. (2013). Nucleic Acids Res 41(6):e74. https://doi.org/10.1093/nar/gkt006
    /// </remarks>
    public static double CalculateCodingPotential(
        string sequence,
        IReadOnlyDictionary<string, double> codingHexamerFrequencies,
        IReadOnlyDictionary<string, double> noncodingHexamerFrequencies,
        int wordSize = HexamerWordSize,
        int stepSize = HexamerStepSize)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(codingHexamerFrequencies);
        ArgumentNullException.ThrowIfNull(noncodingHexamerFrequencies);
        if (wordSize <= 0) throw new ArgumentOutOfRangeException(nameof(wordSize));
        if (stepSize <= 0) throw new ArgumentOutOfRangeException(nameof(stepSize));

        // kmer_ratio: sequences shorter than one word yield no hexamer → score 0.
        if (sequence.Length < wordSize) return 0;

        string seq = sequence.ToUpperInvariant();

        double sumOfLogRatio = 0.0;
        int scoredHexamers = 0;

        // word_generator(frame=0): start at 0, advance by stepSize, keep only
        // full-length words.
        for (int i = 0; i + wordSize <= seq.Length; i += stepSize)
        {
            string kmer = seq.Substring(i, wordSize);

            bool inCoding = codingHexamerFrequencies.TryGetValue(kmer, out double coding);
            bool inNoncoding = noncodingHexamerFrequencies.TryGetValue(kmer, out double noncoding);

            // Skip hexamers absent from either table (matches has_key guard).
            if (!inCoding || !inNoncoding) continue;

            if (coding > 0 && noncoding > 0)
                sumOfLogRatio += Math.Log(coding / noncoding);
            else if (coding > 0 && noncoding == 0)
                sumOfLogRatio += 1;
            else if (coding == 0 && noncoding > 0)
                sumOfLogRatio -= 1;
            // coding == 0 && noncoding == 0 contributes 0 but still counts.

            scoredHexamers++;
        }

        if (scoredHexamers == 0) return 0;

        return sumOfLogRatio / scoredHexamers;
    }

    // Tandem-repeat unit upper bound. Microsatellite (STR) motifs are 1-6 bp and
    // minisatellite (VNTR) motifs 10-60 bp per Wikipedia "Tandem repeat" (cites
    // Duitama et al. 2014, Nucleic Acids Res 42(9):5728-5741). 60 bp covers both.
    private const int MaxTandemUnitLength = 60;

    // Inverted-repeat arm upper bound (each arm of the WGW^R stem). Bounded for
    // tractability of the O(n^2) scan; arms beyond this are not reported.
    private const int MaxInvertedArmLength = 100;

    // Default maximum gap (loop/spacer length |G|) between inverted-repeat arms.
    // The IUPACpal definition allows |G| >= 0 (Hampson et al. 2021); the cited
    // hairpin/stem-loop literature uses loops up to ~50 nt, matching the repo's
    // RepeatFinder default loop bound.
    private const int DefaultMaxInvertedGap = 50;

    /// <summary>
    /// Finds repetitive elements in a DNA sequence: tandem repeats (head-to-tail
    /// consecutive copies of a primitive motif) and inverted repeats (a sequence
    /// followed downstream by its reverse complement, form WGW^R).
    /// </summary>
    /// <param name="dnaSequence">DNA sequence to scan (case-insensitive).</param>
    /// <param name="minRepeatLength">Minimum total span (bp) of a reported element.</param>
    /// <param name="minCopies">Minimum number of adjacent copies for a tandem repeat (>= 2).</param>
    /// <returns>
    /// Tuples of (start, end, type, sequence) where start is 0-based inclusive,
    /// end is exclusive, and type is "tandem_repeat" or "inverted_repeat".
    /// </returns>
    /// <exception cref="ArgumentNullException">dnaSequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">minRepeatLength &lt; 1 or minCopies &lt; 2.</exception>
    public static IEnumerable<(int start, int end, string type, string sequence)> FindRepetitiveElements(
        string dnaSequence,
        int minRepeatLength = 10,
        int minCopies = 2)
    {
        ArgumentNullException.ThrowIfNull(dnaSequence);
        // A tandem repeat requires "two or more" adjacent copies (Wikipedia: Tandem repeat).
        if (minCopies < 2) throw new ArgumentOutOfRangeException(nameof(minCopies));
        if (minRepeatLength < 1) throw new ArgumentOutOfRangeException(nameof(minRepeatLength));

        return FindRepetitiveElementsCore(dnaSequence, minRepeatLength, minCopies);
    }

    private static IEnumerable<(int start, int end, string type, string sequence)> FindRepetitiveElementsCore(
        string dnaSequence,
        int minRepeatLength,
        int minCopies)
    {
        if (dnaSequence.Length == 0)
            yield break;

        string sequence = dnaSequence.ToUpperInvariant();

        foreach (var repeat in FindTandemRepeats(sequence, minRepeatLength, minCopies))
            yield return repeat;

        foreach (var repeat in FindInvertedRepeats(sequence, minRepeatLength, DefaultMaxInvertedGap))
            yield return repeat;
    }

    private static IEnumerable<(int start, int end, string type, string sequence)> FindTandemRepeats(
        string sequence,
        int minLength,
        int minCopies)
    {
        // Span de-duplication: a single maximal tandem array can be discovered from
        // several start/unit-length combinations (e.g. AAAAAA via unit "A" or "AA").
        // We report each maximal array once, using its primitive (shortest) period.
        var reported = new HashSet<(int Start, int End)>();

        for (int unitLen = 1; unitLen <= sequence.Length / minCopies && unitLen <= MaxTandemUnitLength; unitLen++)
        {
            int lastStart = sequence.Length - unitLen * minCopies;
            for (int start = 0; start <= lastStart; start++)
            {
                string unit = sequence.Substring(start, unitLen);

                // Only consider primitive units: a unit that is itself a tandem array
                // of a shorter period (e.g. "AA" = "A"x2) double-counts. Skip it so the
                // primitive period is the one reported (Wikipedia: non-primitive corner case).
                if (!IsPrimitive(unit))
                    continue;

                int copies = 1;
                int pos = start + unitLen;
                while (pos + unitLen <= sequence.Length &&
                       string.CompareOrdinal(sequence, pos, sequence, start, unitLen) == 0)
                {
                    copies++;
                    pos += unitLen;
                }

                if (copies < minCopies)
                    continue;

                int end = start + unitLen * copies;
                if (end - start < minLength)
                    continue;

                // Only report a left-maximal array: if the previous unit-length block
                // equals the unit, this array was already counted from an earlier start.
                if (start - unitLen >= 0 &&
                    string.CompareOrdinal(sequence, start - unitLen, sequence, start, unitLen) == 0)
                    continue;

                if (reported.Add((start, end)))
                {
                    string repeatSeq = sequence.Substring(start, end - start);
                    yield return (start, end, "tandem_repeat", repeatSeq);
                }
            }
        }
    }

    /// <summary>
    /// Returns true if <paramref name="unit"/> is primitive (not itself a tandem
    /// repeat of a shorter period). "ATAT" is not primitive ("AT"x2); "AT" is.
    /// </summary>
    private static bool IsPrimitive(string unit)
    {
        int n = unit.Length;
        for (int period = 1; period <= n / 2; period++)
        {
            if (n % period != 0) continue;
            bool isRepeat = true;
            for (int i = period; i < n; i++)
            {
                if (unit[i] != unit[i % period]) { isRepeat = false; break; }
            }
            if (isRepeat) return false;
        }
        return true;
    }

    private static IEnumerable<(int start, int end, string type, string sequence)> FindInvertedRepeats(
        string sequence,
        int minArmLength,
        int maxGap)
    {
        // Inverted repeat = left arm W followed by gap G (|G|>=0) then right arm W^R
        // (reverse complement of W). Per Hampson et al. (2021): form WGW^R.
        var reported = new HashSet<(int Start, int End)>();

        for (int i = 0; i + minArmLength * 2 <= sequence.Length; i++)
        {
            int maxLen = Math.Min(MaxInvertedArmLength, (sequence.Length - i) / 2);
            for (int len = minArmLength; len <= maxLen; len++)
            {
                string arm1 = sequence.Substring(i, len);
                string arm1RevComp = DnaSequence.GetReverseComplementString(arm1);

                // Right arm must start after the left arm plus a gap of 0..maxGap.
                int firstRight = i + len;
                int lastRight = Math.Min(sequence.Length - len, firstRight + maxGap);
                for (int j = firstRight; j <= lastRight; j++)
                {
                    if (string.CompareOrdinal(sequence, j, arm1RevComp, 0, len) != 0)
                        continue;

                    int end = j + len;
                    if (!reported.Add((i, end)))
                        continue;

                    string gap = sequence.Substring(firstRight, j - firstRight);
                    string arm2 = sequence.Substring(j, len);
                    string composite = gap.Length == 0 ? arm1 + arm2 : arm1 + "[" + gap + "]" + arm2;
                    yield return (i, end, "inverted_repeat", composite);
                }
            }
        }
    }

    /// <summary>
    /// Classifies a repeat sequence against a repeat-element library, mirroring the
    /// RepeatMasker convention of assigning the class of the best-matching library
    /// entry. Matching here is exact substring containment (either direction); when
    /// no library entry matches, the query falls back to a simple-repeat class by
    /// primitive motif size (STR 1-6 bp) or "Unknown".
    /// </summary>
    /// <param name="sequence">Repeat sequence to classify (case-insensitive).</param>
    /// <param name="repeatDb">
    /// Library mapping a known repeat element sequence to its class label
    /// (e.g. "SINE/Alu", "LINE/L1", "LTR", "DNA", "Satellite").
    /// </param>
    /// <returns>The class label of the matched library entry, or a fallback class.</returns>
    /// <exception cref="ArgumentNullException">sequence or repeatDb is null.</exception>
    public static string ClassifyRepeat(string sequence, IReadOnlyDictionary<string, string> repeatDb)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(repeatDb);

        if (sequence.Length == 0)
            return SimpleRepeatClass;

        string query = sequence.ToUpperInvariant();

        // Best library match: prefer the longest known element that is contained in
        // the query (or that contains the query) — approximates RepeatMasker's
        // "best match above a minimum score" against the Repbase library.
        string? bestClass = null;
        int bestLen = -1;
        foreach (var entry in repeatDb)
        {
            string element = entry.Key.ToUpperInvariant();
            if (element.Length == 0) continue;
            bool matches = query.Contains(element, StringComparison.Ordinal) ||
                           element.Contains(query, StringComparison.Ordinal);
            if (matches && element.Length > bestLen)
            {
                bestLen = element.Length;
                bestClass = entry.Value;
            }
        }

        if (bestClass is not null)
            return bestClass;

        // No library match: if the query is a tandem repeat of a short primitive
        // motif (1-6 bp), it is a simple repeat (microsatellite); otherwise Unknown.
        return IsSimpleRepeat(query) ? SimpleRepeatClass : UnknownClass;
    }

    // RepeatMasker output class labels (Smit/Hubley/Green; Repbase Update).
    private const string SimpleRepeatClass = "Simple_repeat";
    private const string UnknownClass = "Unknown";

    // STR/microsatellite motifs are 1-6 bp (Wikipedia: Tandem repeat).
    private const int MaxSimpleRepeatMotif = 6;

    /// <summary>
    /// True if <paramref name="query"/> is a perfect tandem repeat of a primitive
    /// motif of length 1-6 bp (a microsatellite / simple repeat).
    /// </summary>
    private static bool IsSimpleRepeat(string query)
    {
        int n = query.Length;
        for (int period = 1; period <= MaxSimpleRepeatMotif && period <= n / 2; period++)
        {
            if (n % period != 0) continue;
            bool isRepeat = true;
            for (int i = period; i < n; i++)
            {
                if (query[i] != query[i % period]) { isRepeat = false; break; }
            }
            if (isRepeat) return true;
        }
        return false;
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
