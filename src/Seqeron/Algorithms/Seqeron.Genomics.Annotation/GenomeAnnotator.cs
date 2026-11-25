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
    /// <param name="Source">
    /// GFF3 column 2 (source): "a free text qualifier intended to describe the algorithm or
    /// operating procedure that generated this feature" (Sequence Ontology GFF3 Spec v1.26).
    /// Defaults to "." (the undefined-field placeholder) so features constructed without a
    /// source still serialise to spec-valid output.
    /// </param>
    /// <param name="Score">
    /// GFF3 column 6 (score): "a floating point number" whose semantics are ill-defined
    /// (Sequence Ontology GFF3 Spec v1.26). <see langword="null"/> means undefined and is
    /// serialised as the "." placeholder.
    /// </param>
    /// <param name="AttributeValues">
    /// Optional multi-valued GFF3 column-9 attributes ("Multiple attributes of the same type are
    /// indicated by separating the values with the comma ',' character, as in
    /// <c>Parent=AF2312,AB2812,abc-3</c>" — GFF3 Spec v1.26). Each value is percent-encoded on its
    /// own and the values are joined with a literal comma, so a value that itself contains a comma
    /// stays one value. A tag present here replaces the same tag in <paramref name="Attributes"/>
    /// (whose values are always written as one value, with commas encoded as <c>%2C</c>).
    /// </param>
    public readonly record struct GeneAnnotation(
        string GeneId,
        int Start,
        int End,
        char Strand,
        string Type,
        string Product,
        IReadOnlyDictionary<string, string> Attributes,
        string Source = UndefinedFieldPlaceholder,
        double? Score = null,
        IReadOnlyDictionary<string, IReadOnlyList<string>>? AttributeValues = null);

    /// <summary>
    /// Represents a genomic feature parsed from one GFF3 feature line. <see cref="Start"/> and
    /// <see cref="End"/> stay in GFF3 space (1-based, inclusive).
    /// </summary>
    /// <param name="Attributes">
    /// Column 9 as tag → percent-decoded value text. For a multi-valued tag the decoded values are
    /// joined with ',' (so <c>Note=a%2Cb</c> and <c>Note=a,b</c> both read "a,b" here); use
    /// <paramref name="AttributeValues"/> to tell them apart.
    /// </param>
    /// <param name="Seqid">Column 1 (landmark ID), percent-decoded.</param>
    /// <param name="Source">Column 2 (source), percent-decoded; "." when undefined.</param>
    /// <param name="AttributeValues">
    /// Column 9 as tag → list of values: the raw value is split on the unescaped ',' separator and
    /// each value is percent-decoded separately (GFF3 Spec v1.26, column 9), e.g.
    /// <c>Parent=mRNA00001,mRNA00002</c> → ["mRNA00001", "mRNA00002"] and <c>Note=a%2Cb</c> → ["a,b"].
    /// Always set by <see cref="ParseGff3"/>.
    /// </param>
    public readonly record struct GenomicFeature(
        string FeatureId,
        string Type,
        int Start,
        int End,
        char Strand,
        double? Score,
        int? Phase,
        IReadOnlyDictionary<string, string> Attributes,
        string Seqid = UndefinedFieldPlaceholder,
        string Source = UndefinedFieldPlaceholder,
        IReadOnlyDictionary<string, IReadOnlyList<string>>? AttributeValues = null);

    /// <summary>
    /// Initiation codons recognised by <see cref="FindOrfs"/>: ATG, GTG and TTG — the
    /// prokaryotic start set used by Prodigal (Hyatt et al. 2010, BMC Bioinformatics 11:119;
    /// <c>is_start()</c> in Prodigal <c>sequence.c</c> accepts exactly ATG/GTG/TTG for
    /// translation table 11). This is a start-selection policy, not a genetic-code table:
    /// stop recognition and codon translation come from the canonical
    /// <see cref="GeneticCode.Standard"/> via <see cref="GeneticCode.IsStopCodon"/> and
    /// <see cref="Translator.Translate(string, GeneticCode?, int, bool)"/>.
    /// </summary>
    private static readonly HashSet<string> StartCodons = new(StringComparer.OrdinalIgnoreCase)
    {
        "ATG", "GTG", "TTG"
    };

    /// <summary>Genetic code used for stop recognition and translation (NCBI table 1).</summary>
    private static GeneticCode OrfGeneticCode => GeneticCode.Standard;

    /// <summary>
    /// Finds all open reading frames in a DNA sequence.
    /// </summary>
    /// <param name="dnaSequence">The DNA sequence to search.</param>
    /// <param name="minLength">Minimum ORF length in amino acids.</param>
    /// <param name="searchBothStrands">Whether to search reverse complement.</param>
    /// <param name="requireStartCodon">
    /// <see langword="true"/> (default): start-to-stop ORFs — every ATG/GTG/TTG start is paired
    /// with its nearest downstream in-frame stop (Rosalind ORF semantics, so nested ORFs sharing
    /// a stop are all reported); an initiating codon is translated as Met, as in EMBOSS getorf
    /// (<c>-methionine</c>, default Y) and Biopython <c>translate(cds=True)</c>.
    /// <see langword="false"/>: stop-to-stop regions exactly as EMBOSS getorf <c>-find 0</c> on a
    /// linear sequence — one region per stop-delimited segment, including the leading segment
    /// that begins at the frame offset and the trailing segment that runs to the last complete
    /// codon without a stop; no start-codon sub-ORFs are added.
    /// </param>
    /// <returns>All found ORFs.</returns>
    public static IEnumerable<OpenReadingFrame> FindOrfs(
        string dnaSequence,
        int minLength = 100,
        bool searchBothStrands = true,
        bool requireStartCodon = true)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            yield break;

        var geneticCode = OrfGeneticCode;

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
        return requireStartCodon
            ? FindStartToStopOrfsInFrame(sequence, frame, geneticCode, minLength, isReverseComplement)
            : FindStopToStopOrfsInFrame(sequence, frame, geneticCode, minLength, isReverseComplement);
    }

    /// <summary>
    /// Start-to-stop ORFs in one frame: every pending start is closed by the next in-frame stop.
    /// </summary>
    private static IEnumerable<OpenReadingFrame> FindStartToStopOrfsInFrame(
        string sequence,
        int frame,
        GeneticCode geneticCode,
        int minLength,
        bool isReverseComplement)
    {
        var currentOrfStarts = new List<int>();

        for (int i = frame; i <= sequence.Length - 3; i += 3)
        {
            string codon = sequence.Substring(i, 3);

            if (StartCodons.Contains(codon))
            {
                currentOrfStarts.Add(i);
            }
            else if (geneticCode.IsStopCodon(codon))
            {
                foreach (int start in currentOrfStarts)
                {
                    int orfAaLength = (i - start) / 3;
                    if (orfAaLength >= minLength)
                    {
                        string orfSeq = sequence.Substring(start, i + 3 - start);
                        yield return new OpenReadingFrame(
                            Start: start,
                            End: i + 3,
                            Frame: frame + 1,
                            IsReverseComplement: isReverseComplement,
                            Sequence: orfSeq,
                            ProteinSequence: TranslateWithInitiatorMet(orfSeq, geneticCode));
                    }
                }

                currentOrfStarts.Clear();
            }
        }
    }

    /// <summary>
    /// Stop-to-stop regions in one frame (EMBOSS getorf <c>-find 0</c>, linear sequence):
    /// the frame is treated as already open at its offset, every stop closes the current
    /// region and opens the next one at the following codon, and the region still open at the
    /// end of the sequence is reported up to the last complete codon. Regions that contain no
    /// sense codon are not ORFs and are skipped.
    /// </summary>
    private static IEnumerable<OpenReadingFrame> FindStopToStopOrfsInFrame(
        string sequence,
        int frame,
        GeneticCode geneticCode,
        int minLength,
        bool isReverseComplement)
    {
        int regionStart = frame;
        int i = frame;

        for (; i <= sequence.Length - 3; i += 3)
        {
            if (!geneticCode.IsStopCodon(sequence.Substring(i, 3)))
                continue;

            int aaLength = (i - regionStart) / 3;
            if (aaLength >= 1 && aaLength >= minLength)
            {
                string orfSeq = sequence.Substring(regionStart, i + 3 - regionStart);
                yield return new OpenReadingFrame(
                    Start: regionStart,
                    End: i + 3,
                    Frame: frame + 1,
                    IsReverseComplement: isReverseComplement,
                    Sequence: orfSeq,
                    ProteinSequence: Translator.Translate(orfSeq, geneticCode).Sequence);
            }

            regionStart = i + 3;
        }

        // Trailing region without a stop: i is now the first codon position that does not fit.
        int endPos = i;
        int trailingAa = (endPos - regionStart) / 3;
        if (endPos > regionStart && trailingAa >= minLength)
        {
            string orfSeq = sequence.Substring(regionStart, endPos - regionStart);
            yield return new OpenReadingFrame(
                Start: regionStart,
                End: endPos,
                Frame: frame + 1,
                IsReverseComplement: isReverseComplement,
                Sequence: orfSeq,
                ProteinSequence: Translator.Translate(orfSeq, geneticCode).Sequence);
        }
    }

    /// <summary>
    /// Translates a start-codon-initiated ORF with the canonical translator, rendering the
    /// initiating codon as Met: an alternative initiator (GTG/TTG) is decoded by initiator
    /// tRNA-fMet, so EMBOSS getorf (<c>-methionine</c>, default Y: "START codons at the
    /// beginning of protein products will usually code for Methionine") and Biopython
    /// <c>Seq.translate(cds=True)</c> both emit M for it (GTGAAAAAAAAATAA → MKKK).
    /// </summary>
    private static string TranslateWithInitiatorMet(string orfSeq, GeneticCode geneticCode)
    {
        string protein = Translator.Translate(orfSeq, geneticCode).Sequence;
        return protein.Length == 0 || protein[0] == 'M' ? protein : string.Concat("M", protein.AsSpan(1));
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
    /// Exact Shine-Dalgarno motifs: every contiguous substring of length ≥ 4 of the consensus
    /// AGGAGG (Shine &amp; Dalgarno 1975, Nature 254:34-38), which pairs with the anti-SD
    /// 3' tail of 16S rRNA (…CCUCCU…). This is the exact-match motif set of Prodigal's
    /// <c>shine_dalgarno_exact()</c> (Hyatt et al. 2010, BMC Bioinformatics 11:119;
    /// <c>sequence.c</c>), which compares a window to AGGAGG and accepts any contiguous
    /// matching sub-motif — its 4-base class is "AGGA/GGAG/GAGG" (<c>gene.c</c> bins 11/12),
    /// so GGAG is included. Prodigal's 3-base class (GGA/GAG/AGG) is not scanned here.
    /// Longest first.
    /// </summary>
    private static readonly string[] ShineDalgarnoMotifs =
        { "AGGAGG", "GGAGG", "AGGAG", "GGAG", "GAGG", "AGGA" };

    /// <summary>
    /// Finds potential Shine-Dalgarno (ribosome binding site) sequences on the FORWARD
    /// strand only. The motif AGGAGG (and shorter variants) is sought upstream of every
    /// forward-strand ORF start codon at an aligned spacing within [minDistance, maxDistance].
    /// </summary>
    /// <remarks>
    /// This overload preserves the original forward-strand-only behaviour. To also report
    /// reverse-strand Shine-Dalgarno hits, use
    /// <see cref="FindRibosomeBindingSitesBothStrands"/>.
    /// </remarks>
    public static IEnumerable<(int position, string sequence, double score)> FindRibosomeBindingSites(
        string dnaSequence,
        int upstreamWindow = 20,
        int minDistance = 4,
        int maxDistance = 15)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            yield break;

        var orfs = FindOrfs(dnaSequence, minLength: 30)
            .Where(o => !o.IsReverseComplement)
            .ToList();

        foreach (var hit in ScanStrandForShineDalgarno(dnaSequence, orfs, upstreamWindow, minDistance, maxDistance))
        {
            yield return (hit.position, hit.sequence, hit.score);
        }
    }

    /// <summary>
    /// Finds potential Shine-Dalgarno (ribosome binding site) sequences on BOTH strands,
    /// reporting the strand of each hit. Forward-strand hits are scanned exactly as in
    /// <see cref="FindRibosomeBindingSites"/>; reverse-strand hits are found by applying the
    /// same scan to the reverse complement (i.e. the mRNA orientation of reverse-strand
    /// genes) and mapping the motif position back to a forward-strand genomic coordinate.
    /// </summary>
    /// <remarks>
    /// For a reverse-strand gene the mRNA is the reverse complement of the forward genomic
    /// strand, so the AGGAGG Shine-Dalgarno motif lies upstream of that gene's start codon
    /// on the reverse strand — which is downstream (higher forward coordinate) of the ORF on
    /// the forward strand. The reported <c>position</c> is the forward-strand index of the
    /// motif's 5' base; <c>sequence</c> is the motif as read 5'→3' on the reverse strand
    /// (so reverse-strand hits read AGGAGG, GGAGG, … just like forward hits); <c>strand</c>
    /// is '+' or '-'.
    /// </remarks>
    public static IEnumerable<(int position, string sequence, double score, char strand)>
        FindRibosomeBindingSitesBothStrands(
        string dnaSequence,
        int upstreamWindow = 20,
        int minDistance = 4,
        int maxDistance = 15)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            yield break;

        var allOrfs = FindOrfs(dnaSequence, minLength: 30).ToList();

        // Forward strand: scan the genomic sequence as-is.
        var forwardOrfs = allOrfs.Where(o => !o.IsReverseComplement).ToList();
        foreach (var hit in ScanStrandForShineDalgarno(dnaSequence, forwardOrfs, upstreamWindow, minDistance, maxDistance))
        {
            yield return (hit.position, hit.sequence, hit.score, '+');
        }

        // Reverse strand: the reverse complement is the mRNA orientation of reverse-strand
        // genes, where the SD motif again reads 5'→3' upstream of the (now forward-facing)
        // start codon. We rebuild the reverse-complement ORFs in that coordinate space and
        // reuse the identical scan, then map each motif's reverse-strand coordinate back to a
        // forward-strand genomic coordinate.
        int len = dnaSequence.Length;
        string revComp = DnaSequence.GetReverseComplementString(dnaSequence);

        var reverseOrfsInRevComp = allOrfs
            .Where(o => o.IsReverseComplement)
            // FindOrfs already mapped Start/End into forward coordinates; undo that mapping to
            // recover the ORF position within the reverse-complement string used for scanning.
            .Select(o => o with { Start = len - o.End, End = len - o.Start })
            .ToList();

        foreach (var hit in ScanStrandForShineDalgarno(revComp, reverseOrfsInRevComp, upstreamWindow, minDistance, maxDistance))
        {
            // hit.position is the 5' base of the motif within the reverse-complement string;
            // the motif occupies revComp[position, position+len) ↔ forward
            // [len - position - motifLen, len - position).
            int forwardPosition = len - hit.position - hit.sequence.Length;
            yield return (forwardPosition, hit.sequence, hit.score, '-');
        }
    }

    /// <summary>
    /// Scans a single strand sequence for Shine-Dalgarno motifs upstream of the supplied ORFs'
    /// start codons. Positions are reported in the coordinate space of <paramref name="sequence"/>.
    /// </summary>
    /// <remarks>
    /// Only maximal motifs are reported: for one start codon, an in-range hit whose span lies
    /// inside another in-range hit (e.g. GGAGG inside AGGAGG) describes the same SD site and is
    /// suppressed, as Prodigal's <c>shine_dalgarno_exact()</c> reports a single maximal
    /// sub-motif of AGGAGG per window (it iterates motif lengths from the longest down). A motif
    /// found upstream of several start codons (nested ORFs) is reported once.
    /// <para>
    /// <b>Score is a heuristic.</b> <c>score = motif.Length / 6</c> ranks motifs by their
    /// complementarity length to the anti-SD only. The published RBS scores are not
    /// reproducible without a genome-specific model: Prodigal scores each (motif, spacer) bin
    /// with log-likelihood weights trained iteratively on the input genome's own gene set
    /// (<c>train_starts_sd()</c>, <c>node.c</c>), and the RBS Calculator (Salis et al. 2009)
    /// sums mRNA/16S rRNA folding free energies with empirically fitted spacing terms.
    /// </para>
    /// </remarks>
    private static IEnumerable<(int position, string sequence, double score)> ScanStrandForShineDalgarno(
        string sequence,
        IReadOnlyList<OpenReadingFrame> orfs,
        int upstreamWindow,
        int minDistance,
        int maxDistance)
    {
        var reported = new HashSet<(int position, string motif)>();

        foreach (var orf in orfs)
        {
            int searchStart = Math.Max(0, orf.Start - upstreamWindow);
            int searchEnd = orf.Start - minDistance;

            if (searchEnd <= searchStart) continue;

            string upstream = sequence.Substring(searchStart, searchEnd - searchStart).ToUpperInvariant();
            var hitsForStart = new List<(int position, string motif)>();

            foreach (string motif in ShineDalgarnoMotifs)
            {
                int pos = upstream.IndexOf(motif, StringComparison.Ordinal);
                while (pos >= 0)
                {
                    int genomicPos = searchStart + pos;
                    int distanceToStart = orf.Start - genomicPos - motif.Length;

                    if (distanceToStart >= minDistance && distanceToStart <= maxDistance
                        && !hitsForStart.Any(h => h.position <= genomicPos
                                                  && genomicPos + motif.Length <= h.position + h.motif.Length))
                    {
                        hitsForStart.Add((genomicPos, motif));
                    }

                    pos = upstream.IndexOf(motif, pos + 1, StringComparison.Ordinal);
                }
            }

            foreach (var hit in hitsForStart)
            {
                if (reported.Add(hit))
                    yield return (hit.position, hit.motif, hit.motif.Length / 6.0);
            }
        }
    }

    /// <summary>
    /// Predicts genes using a simple ORF-based approach: one CDS per stop codon on each strand,
    /// started at the most upstream in-frame ATG/GTG/TTG (the longest ORF for that stop).
    /// </summary>
    /// <remarks>
    /// <see cref="FindOrfs"/> reports every start paired with its stop (nested ORFs sharing a
    /// stop), but a gene is one start–stop pair: Prodigal (Hyatt et al. 2010) builds its gene
    /// set from start/stop nodes and selects a single start per stop, and EMBOSS getorf
    /// <c>-find 1</c> reports one region per stop beginning at the first START codon after the
    /// preceding in-frame STOP. Start choice here follows the longest-ORF (first start) rule;
    /// no RBS- or coding-score-based start refinement is performed.
    /// </remarks>
    public static IEnumerable<GeneAnnotation> PredictGenes(
        string dnaSequence,
        int minOrfLength = 100,
        string prefix = "gene")
    {
        var orfs = FindOrfs(dnaSequence, minOrfLength, searchBothStrands: true, requireStartCodon: true)
            // Stop codon identity: forward stop ends at End; reverse stop sits at forward Start.
            .GroupBy(o => (o.IsReverseComplement, stop: o.IsReverseComplement ? o.Start : o.End))
            .Select(g => g.MaxBy(o => o.End - o.Start))
            .OrderBy(o => o.Start)
            .ThenBy(o => o.IsReverseComplement)
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
    /// Parses the feature lines of a GFF3 document (Sequence Ontology GFF3 Spec v1.26).
    /// <para>
    /// Line handling: blank lines and <c>#</c> comment / <c>##</c> directive lines are skipped; a
    /// trailing CR (CRLF line ending) is not part of the last column; the annotation section ends
    /// at the <c>##FASTA</c> directive or at a line beginning with <c>&gt;</c> ("a GFF line that
    /// begins with the character &gt; creates an implied ##FASTA directive"), and nothing after it
    /// is read as a feature.
    /// </para>
    /// <para>
    /// A feature line is accepted only when it satisfies the column rules of the spec — the same
    /// rules GenomeTools <c>gt gff3validator</c> enforces — and is otherwise skipped (one malformed
    /// line never aborts the parse and does not consume a <c>feature_n</c> number):
    /// exactly nine tab-separated columns; start and end positive integers with start ≤ end;
    /// score "." or a floating-point number; strand one of <c>+ - . ?</c>; phase "." or one of
    /// <c>0 1 2</c>; column 9 "." or ';'-separated tag=value pairs with non-empty tag and value and
    /// no tag repeated. A CDS line with phase "." is accepted (the spec
    /// requires a phase for CDS, but <c>gt gff3validator</c> accepts it too).
    /// </para>
    /// <para>
    /// Decoding: columns 1–3 and every tag and value are RFC 3986 percent-decoded. Column 9 is
    /// split on ';' into tag=value pairs (the first '=' separates tag from value; surrounding
    /// whitespace of a tag is dropped, value text is kept verbatim — "unescaped spaces are allowed
    /// within fields"), and each value list is split on the unescaped ',' before decoding, so
    /// <c>Parent=a,b</c> is two values and <c>Note=a%2Cb</c> is one. Coordinates stay in GFF3
    /// space (1-based, inclusive).
    /// </para>
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="lines"/> is null.</exception>
    public static IEnumerable<GenomicFeature> ParseGff3(IEnumerable<string> lines)
    {
        ArgumentNullException.ThrowIfNull(lines);
        return ParseGff3Lines(lines);
    }

    private static IEnumerable<GenomicFeature> ParseGff3Lines(IEnumerable<string> lines)
    {
        int featureCount = 0;

        foreach (string rawLine in lines)
        {
            if (rawLine is null)
                continue;

            string line = rawLine.EndsWith('\r') ? rawLine[..^1] : rawLine;
            if (string.IsNullOrWhiteSpace(line))
                continue;

            // End of the annotation section: explicit ##FASTA or the implied one ('>' line).
            if (line[0] == '>' || line.TrimEnd() == Gff3FastaDirective)
                yield break;

            if (line[0] == '#')
                continue;

            if (!TryParseGff3FeatureLine(line, featureCount + 1, out var feature))
                continue;

            featureCount++;
            yield return feature;
        }
    }

    private const int Gff3ColumnCount = 9;
    private const string Gff3FastaDirective = "##FASTA";
    private const string Gff3StrandSymbols = "+-.?";

    private static bool TryParseGff3FeatureLine(string line, int featureNumber, out GenomicFeature feature)
    {
        feature = default;

        string[] parts = line.Split('\t');
        if (parts.Length != Gff3ColumnCount)
            return false;

        var invariant = System.Globalization.CultureInfo.InvariantCulture;

        // Columns 4/5: "positive 1-based integer coordinates … Start is always less than or equal to end".
        if (!int.TryParse(parts[3], System.Globalization.NumberStyles.None, invariant, out int start) ||
            !int.TryParse(parts[4], System.Globalization.NumberStyles.None, invariant, out int end) ||
            start < 1 || end < start)
            return false;

        double? score = null;
        if (parts[5] != UndefinedFieldPlaceholder)
        {
            if (!double.TryParse(parts[5], System.Globalization.NumberStyles.Float, invariant, out double parsedScore))
                return false;
            score = parsedScore;
        }

        // Column 7: one of '+', '-', '.', '?'.
        if (parts[6].Length != 1 || !Gff3StrandSymbols.Contains(parts[6][0]))
            return false;
        char strand = parts[6][0];

        // Column 8: "The phase is one of the integers 0, 1, or 2" ('.' when undefined).
        int? phase = null;
        if (parts[7] != UndefinedFieldPlaceholder)
        {
            if (parts[7].Length != 1 || parts[7][0] < '0' || parts[7][0] > '2')
                return false;
            phase = parts[7][0] - '0';
        }

        if (!TryParseGff3Attributes(parts[8], out var attributeValues))
            return false;

        var attributes = new Dictionary<string, string>(attributeValues.Count, StringComparer.Ordinal);
        foreach (var (tag, values) in attributeValues)
            attributes[tag] = string.Join(',', values);

        feature = new GenomicFeature(
            FeatureId: attributes.GetValueOrDefault("ID", $"feature_{featureNumber}"),
            Type: DecodeGff3(parts[2]),
            Start: start,
            End: end,
            Strand: strand,
            Score: score,
            Phase: phase,
            Attributes: attributes,
            Seqid: DecodeGff3(parts[0]),
            Source: DecodeGff3(parts[1]),
            AttributeValues: attributeValues);
        return true;
    }

    /// <summary>
    /// Splits GFF3 column 9 into tag → decoded value list. Returns false — as
    /// <c>gt gff3validator</c> does — for a token that is not tag=value, an empty tag or value, or
    /// a tag that occurs twice (multiple values of one tag must be comma-separated in a single
    /// tag=value pair).
    /// </summary>
    private static bool TryParseGff3Attributes(
        string column,
        out Dictionary<string, IReadOnlyList<string>> attributeValues)
    {
        attributeValues = new Dictionary<string, IReadOnlyList<string>>(StringComparer.Ordinal);
        if (column == UndefinedFieldPlaceholder)
            return true;

        foreach (string pair in column.Split(';'))
        {
            // Empty tokens come from a trailing ';' or "; " separators and carry nothing.
            if (string.IsNullOrWhiteSpace(pair))
                continue;

            // Every token must be tag=value with a non-empty tag and value.
            int eqPos = pair.IndexOf('=');
            if (eqPos < 0 || eqPos == pair.Length - 1)
                return false;

            string tag = DecodeGff3(pair[..eqPos].Trim());
            if (tag.Length == 0)
                return false;

            string[] values = pair[(eqPos + 1)..].Split(',');
            for (int i = 0; i < values.Length; i++)
                values[i] = DecodeGff3(values[i]);

            if (!attributeValues.TryAdd(tag, values))
                return false;
        }

        return true;
    }

    /// <summary>RFC 3986 percent-decoding of one GFF3 field (GFF3 Spec v1.26, "Description of the Format").</summary>
    private static string DecodeGff3(string value) => Uri.UnescapeDataString(value);

    /// <summary>
    /// The GFF3 undefined-field placeholder: "Undefined fields are replaced with the '.' character"
    /// (Sequence Ontology GFF3 Spec v1.26, column-definition preamble).
    /// </summary>
    private const string UndefinedFieldPlaceholder = ".";

    /// <summary>
    /// Exports annotations in GFF3 format: the <c>##gff-version 3</c> directive, then one
    /// nine-column feature line per annotation, in input order.
    /// <para>
    /// Coordinates: the internal 0-based half-open [Start, End) becomes 1-based inclusive
    /// start = Start + 1, end = End. An empty interval [s, s) (s ≥ 1) is written as the spec's
    /// zero-length feature start = end = s ("the implied site is to the right of the indicated
    /// base").
    /// </para>
    /// <para>
    /// Column 2 (source) and column 6 (score) carry the feature's real values when present and
    /// fall back to the "." undefined-field placeholder otherwise. Column 8 (phase) is computed
    /// per the GFF3 spec for CDS features: the number of bases forward from the start of the
    /// current CDS feature to the next codon boundary, accumulated across the ordered CDS
    /// segments of a transcript (grouped by <see cref="GeneAnnotation.GeneId"/>); non-CDS
    /// features emit ".".
    /// </para>
    /// <para>
    /// Column 9: <c>ID</c> = <see cref="GeneAnnotation.GeneId"/> and <c>product</c> =
    /// <see cref="GeneAnnotation.Product"/> first (these two tags are taken only from those
    /// properties, never repeated from the dictionaries), then <see cref="GeneAnnotation.Attributes"/>
    /// (except <c>translation</c>) and <see cref="GeneAnnotation.AttributeValues"/>; a tag with an
    /// empty value is omitted, and "." is written when no tag remains.
    /// </para>
    /// <para>
    /// Encoding (GFF3 Spec v1.26): tab, newline, CR, '%' and control characters are
    /// percent-encoded in every column, "no other characters may be encoded"; in column 9 tags and
    /// values also encode ';', '=', '&amp;' and ','; column 1 (seqid) encodes every character
    /// outside <c>[a-zA-Z0-9.:^*$@!+_?-|]</c> (UTF-8 bytes), so it never contains whitespace or
    /// begins with '&gt;'.
    /// </para>
    /// Source: Sequence Ontology GFF3 Spec v1.26 — column definitions; "The phase is REQUIRED for
    /// all CDS features."
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="annotations"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// <paramref name="seqId"/> is null or empty, or an annotation has an empty type, a strand
    /// outside <c>+ - . ?</c>, a negative start, an end before its start, or the empty interval
    /// [0, 0) (no base to the left of position 1).
    /// </exception>
    public static IEnumerable<string> ToGff3(
        IEnumerable<GeneAnnotation> annotations,
        string seqId = "seq1")
    {
        ArgumentNullException.ThrowIfNull(annotations);
        ArgumentException.ThrowIfNullOrEmpty(seqId);

        // Phase is a per-transcript cumulative quantity, so we must know the order of a
        // transcript's CDS segments before we can emit any of its rows. Materialise once and
        // precompute each CDS row's phase, keyed by the row's position in the input sequence so
        // emission order (and thus output line order) is preserved exactly.
        var annotationList = annotations as IReadOnlyList<GeneAnnotation> ?? annotations.ToList();
        var cdsPhases = ComputeCdsPhases(annotationList);
        string seqidColumn = EncodeGff3Seqid(seqId);
        var invariant = System.Globalization.CultureInfo.InvariantCulture;

        var output = new List<string>(annotationList.Count + 1) { "##gff-version 3" };
        for (int i = 0; i < annotationList.Count; i++)
        {
            var ann = annotationList[i];
            if (string.IsNullOrEmpty(ann.Type))
                throw new ArgumentException($"Annotation {i} has no feature type (GFF3 column 3).", nameof(annotations));
            if (!Gff3StrandSymbols.Contains(ann.Strand))
                throw new ArgumentException($"Annotation {i} has strand '{ann.Strand}'; GFF3 column 7 allows only '+', '-', '.', '?'.", nameof(annotations));
            if (ann.Start < 0 || ann.End < ann.Start || (ann.End == ann.Start && ann.Start == 0))
                throw new ArgumentException($"Annotation {i} has the interval [{ann.Start}, {ann.End}), which has no GFF3 representation (positive 1-based start ≤ end).", nameof(annotations));

            int gffStart = ann.End == ann.Start ? ann.Start : ann.Start + 1;

            string source = string.IsNullOrEmpty(ann.Source) ? UndefinedFieldPlaceholder : EncodeGff3Text(ann.Source, encodeAttributeReserved: false);
            string score = ann.Score.HasValue
                ? ann.Score.Value.ToString(invariant)
                : UndefinedFieldPlaceholder;
            string phase = cdsPhases.TryGetValue(i, out int p)
                ? p.ToString(invariant)
                : UndefinedFieldPlaceholder;
            string type = EncodeGff3Text(ann.Type, encodeAttributeReserved: false);

            output.Add(string.Join('\t',
                seqidColumn, source, type,
                gffStart.ToString(invariant), ann.End.ToString(invariant),
                score, ann.Strand.ToString(), phase, FormatGff3Attributes(ann)));
        }

        return output;
    }

    /// <summary>
    /// Computes the GFF3 column-8 phase for every CDS feature in <paramref name="annotations"/>,
    /// returning a map from the feature's input index to its phase (0, 1, or 2).
    /// <para>
    /// Per the Sequence Ontology GFF3 Spec v1.26, phase is "the number of bases forward from the
    /// start of the current CDS feature the next codon begins", measured relative to the 5' end
    /// (the feature's start on the + strand, its end on the − strand). For a transcript's first
    /// (5'-most) CDS segment the coding frame begins at the segment start, so phase is 0; for each
    /// subsequent segment the phase is the number of bases needed to complete the codon left open
    /// by the preceding segments: <c>(3 − (Σ preceding segment lengths) mod 3) mod 3</c>.
    /// CDS segments are grouped by <see cref="GeneAnnotation.GeneId"/> (the transcript) and ordered
    /// 5'→3': ascending start on the + strand, descending start on the − strand.
    /// </para>
    /// </summary>
    private static Dictionary<int, int> ComputeCdsPhases(IReadOnlyList<GeneAnnotation> annotations)
    {
        var phases = new Dictionary<int, int>();

        // Group the CDS rows of each transcript, tracking their original input index.
        var cdsByTranscript = new Dictionary<string, List<(int Index, GeneAnnotation Ann)>>(StringComparer.Ordinal);
        for (int i = 0; i < annotations.Count; i++)
        {
            var ann = annotations[i];
            if (!string.Equals(ann.Type, "CDS", StringComparison.OrdinalIgnoreCase))
                continue;

            if (!cdsByTranscript.TryGetValue(ann.GeneId, out var segments))
            {
                segments = new List<(int, GeneAnnotation)>();
                cdsByTranscript[ann.GeneId] = segments;
            }
            segments.Add((i, ann));
        }

        foreach (var segments in cdsByTranscript.Values)
        {
            // Order segments 5'→3'. The 5' end is the start on the + strand and the end on the
            // − strand (GFF3 Spec v1.26, column 8), so − strand segments run from the highest
            // genomic coordinate downward.
            bool minusStrand = segments.Count > 0 && segments[0].Ann.Strand == '-';
            var ordered = minusStrand
                ? segments.OrderByDescending(s => s.Ann.Start).ToList()
                : segments.OrderBy(s => s.Ann.Start).ToList();

            int precedingLength = 0;
            foreach (var (index, ann) in ordered)
            {
                // (3 − (preceding mod 3)) mod 3 — bases to remove to reach the next codon.
                phases[index] = (CodonLength - (precedingLength % CodonLength)) % CodonLength;

                // Internal coordinates are 0-based half-open [Start, End); the 1-based inclusive
                // segment length is End − Start.
                precedingLength += ann.End - ann.Start;
            }
        }

        return phases;
    }

    private static string FormatGff3Attributes(GeneAnnotation ann)
    {
        var pairs = new List<string>();

        void Add(string tag, IEnumerable<string?> values)
        {
            var encoded = values
                .Where(v => !string.IsNullOrEmpty(v))
                .Select(v => EncodeGff3Text(v!, encodeAttributeReserved: true))
                .ToList();
            if (tag.Length > 0 && encoded.Count > 0)
                pairs.Add($"{EncodeGff3Text(tag, encodeAttributeReserved: true)}={string.Join(',', encoded)}");
        }

        static bool IsReservedOrSkipped(string tag) =>
            tag is "ID" or "product" or "translation"; // ID/product come from GeneId/Product; translation is too large

        Add("ID", new[] { ann.GeneId });
        Add("product", new[] { ann.Product });

        var multi = ann.AttributeValues;
        if (ann.Attributes is not null)
        {
            foreach (var (tag, value) in ann.Attributes)
            {
                if (!IsReservedOrSkipped(tag) && (multi is null || !multi.ContainsKey(tag)))
                    Add(tag, new[] { value });
            }
        }

        if (multi is not null)
        {
            foreach (var (tag, values) in multi)
            {
                if (!IsReservedOrSkipped(tag) && values is not null)
                    Add(tag, values);
            }
        }

        return pairs.Count == 0 ? UndefinedFieldPlaceholder : string.Join(';', pairs);
    }

    /// <summary>
    /// Percent-encodes a GFF3 field per the GFF3 Specification v1.26: tab (%09), newline (%0A),
    /// CR (%0D), percent (%25) and control characters (%00–%1F, %7F) everywhere; with
    /// <paramref name="encodeAttributeReserved"/> also the column-9 reserved characters
    /// ';' (%3B), '=' (%3D), '&amp;' (%26), ',' (%2C). "No other characters may be encoded" —
    /// in particular spaces stay literal ("unescaped spaces are allowed within fields").
    /// </summary>
    private static string EncodeGff3Text(string value, bool encodeAttributeReserved)
    {
        var sb = new StringBuilder(value.Length);
        foreach (char c in value)
        {
            bool encode = c < 0x20 || c == 0x7F || c == '%'
                || (encodeAttributeReserved && c is ';' or '=' or '&' or ',');
            if (encode)
                sb.Append('%').Append(((int)c).ToString("X2", System.Globalization.CultureInfo.InvariantCulture));
            else
                sb.Append(c);
        }
        return sb.ToString();
    }

    /// <summary>
    /// Encodes GFF3 column 1: "IDs may contain any characters, but must escape any characters not
    /// in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not contain unescaped whitespace
    /// and must not begin with an unescaped '&gt;'" (GFF3 Spec v1.26). Characters outside the set
    /// are written as the percent-encoded bytes of their UTF-8 form (RFC 3986 §2.1).
    /// </summary>
    private static string EncodeGff3Seqid(string seqId)
    {
        var sb = new StringBuilder(seqId.Length);
        Span<byte> utf8 = stackalloc byte[4];
        foreach (Rune rune in seqId.EnumerateRunes())
        {
            if (rune.IsAscii && IsGff3SeqidChar((char)rune.Value))
            {
                sb.Append((char)rune.Value);
                continue;
            }

            int count = rune.EncodeToUtf8(utf8);
            for (int i = 0; i < count; i++)
                sb.Append('%').Append(utf8[i].ToString("X2", System.Globalization.CultureInfo.InvariantCulture));
        }
        return sb.ToString();
    }

    private static bool IsGff3SeqidChar(char c) =>
        char.IsAsciiLetterOrDigit(c) || ".:^*$@!+_?-|".Contains(c);

    /// <summary>
    /// Finds bacterial σ70 promoter boxes: the −35 and −10 hexamers of every −35/−10 pair returned by the canonical
    /// <see cref="Seqeron.Genomics.Analysis.MotifFinder.FindSigma70Promoters"/> (Harley &amp; Reynolds 1987: consensus
    /// TTGACA / TATAAT, spacer 15–21 bp, 17 ± 1 bp in 92 % of promoters), with that method's defaults
    /// (≤ 2 mismatches per box, spacer 15–21, given strand).
    /// </summary>
    /// <remarks>See <see cref="FindPromoterMotifs(string, int, int, int, int)"/>.</remarks>
    /// <param name="dnaSequence">DNA sequence (A/C/G/T, case-insensitive).</param>
    /// <exception cref="ArgumentNullException"><paramref name="dnaSequence"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="dnaSequence"/> contains a non-ACGT character.</exception>
    public static IEnumerable<(int position, string type, string sequence, double score)> FindPromoterMotifs(
        string dnaSequence)
        => FindPromoterMotifs(dnaSequence, maxMismatches35: 2, maxMismatches10: 2, minSpacer: 15, maxSpacer: 21);

    /// <summary>
    /// Finds bacterial σ70 promoter boxes: the −35 and −10 hexamers of every −35/−10 pair returned by the canonical
    /// <see cref="Seqeron.Genomics.Analysis.MotifFinder.FindSigma70Promoters"/> (Harley &amp; Reynolds 1987, NAR
    /// 15:2343: consensus TTGACA / TATAAT, 17 ± 1 bp spacer in 92 % of promoters).
    /// </summary>
    /// <remarks>
    /// <para>A box is reported only when it is part of at least one −35/−10 pair whose hexamers lie within the
    /// mismatch limits of their consensus and whose spacer is <paramref name="minSpacer"/> … <paramref name="maxSpacer"/>
    /// bp; a box shared by several pairs is reported once. Forward strand only (the tuple carries no strand); for
    /// pairs, strand and spacer use <see cref="Seqeron.Genomics.Analysis.MotifFinder.FindSigma70Promoters"/>
    /// directly.</para>
    /// <para><c>score</c> = fraction of the six consensus positions matched, 1 − d/6, where d is the box's Hamming
    /// distance to its consensus as reported by <c>FindSigma70Promoters</c>. It measures consensus identity only, not
    /// promoter strength; for a published quantitative promoter score (free energy, transcription rate) use
    /// <see cref="Seqeron.Genomics.Analysis.MotifFinder.PredictSigma70Promoters"/> (Promoter Calculator v1.0).</para>
    /// <para>Output: 0-based start, type <c>"-35 box"</c> / <c>"-10 box"</c>, the hexamer, score; ordered by
    /// position, then −35 before −10.</para>
    /// </remarks>
    /// <param name="dnaSequence">DNA sequence (A/C/G/T, case-insensitive).</param>
    /// <param name="maxMismatches35">Maximum mismatches of the −35 box to TTGACA (0 … 6).</param>
    /// <param name="maxMismatches10">Maximum mismatches of the −10 box to TATAAT (0 … 6).</param>
    /// <param name="minSpacer">Minimum spacer length (bp).</param>
    /// <param name="maxSpacer">Maximum spacer length (bp).</param>
    /// <exception cref="ArgumentNullException"><paramref name="dnaSequence"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="dnaSequence"/> contains a non-ACGT character.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid mismatch or spacer limits (see <c>FindSigma70Promoters</c>).</exception>
    public static IEnumerable<(int position, string type, string sequence, double score)> FindPromoterMotifs(
        string dnaSequence, int maxMismatches35, int maxMismatches10, int minSpacer, int maxSpacer)
    {
        ArgumentNullException.ThrowIfNull(dnaSequence);
        const int box = 6;

        var candidates = Seqeron.Genomics.Analysis.MotifFinder.FindSigma70Promoters(
            new DnaSequence(dnaSequence), maxMismatches35, maxMismatches10, minSpacer, maxSpacer);

        var boxes = new SortedSet<(int position, int order, string type, string sequence, double score)>();
        foreach (var c in candidates)
        {
            boxes.Add((c.Minus35Start, 0, "-35 box", c.Minus35, (box - c.Mismatches35) / (double)box));
            boxes.Add((c.Minus10Start, 1, "-10 box", c.Minus10, (box - c.Mismatches10) / (double)box));
        }

        return boxes.Select(b => (b.position, b.type, b.sequence, b.score)).ToList();
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
    /// Hexamer score returned when no in-frame hexamer is scorable: CPAT
    /// <c>FrameKmer.kmer_ratio</c> catches the division by a zero count and
    /// returns −1 (CPAT 3.0.5, liguowang/cpat <c>src/cpmodule/FrameKmer.py</c>).
    /// </summary>
    private const double NoScorableHexamerScore = -1.0;

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
    /// sequence is shorter than <paramref name="wordSize"/>, and −1 when the
    /// sequence is long enough but no in-frame hexamer is scorable (every word is
    /// missing from a table or has both frequencies 0) — both exactly as the CPAT
    /// reference <c>kmer_ratio</c>.
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
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(wordSize);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(stepSize);

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
            else if (coding == 0 && noncoding == 0)
                continue; // both-zero: skipped, NOT counted (FrameKmer.kmer_ratio).
            else if (coding == 0 && noncoding > 0)
                sumOfLogRatio -= 1;
            else
                continue; // any other case (e.g. negative values): skipped.

            scoredHexamers++;
        }

        // kmer_ratio: `try: return sum/frame0_count  except: return -1` — when no hexamer
        // was scored (ZeroDivisionError) the reference returns -1, and CPAT feeds that -1
        // into its logistic model as the Hexamer feature (CPAT 3.0.5 cpmodule/FrameKmer.py).
        if (scoredHexamers == 0) return NoScorableHexamerScore;

        return sumOfLogRatio / scoredHexamers;
    }

    #region Fickett TESTCODE (Fickett 1982)

    /// <summary>
    /// Which Fickett TESTCODE look-up tables (and composition convention) to use.
    /// </summary>
    public enum FickettTableSet
    {
        /// <summary>
        /// The published Fickett (1982) tables as distributed with EMBOSS <c>tcode</c>
        /// (<c>Etcode.dat</c>, "Fickett TESTCODE data, NAR 10(17) 5303-5318"): ten position
        /// thresholds 1.9…1.1, 0.0 and ten content thresholds 0.33…0.19, 0.17, 0.00. Base
        /// composition is counted over A/C/G/T only, as <c>tcode_slide</c> does.
        /// </summary>
        Fickett1982,

        /// <summary>
        /// The table set of CPAT 3.0.5 <c>cpmodule/fickett.py</c> (<c>fickett_value</c>): identical
        /// probabilities and weights, but the content thresholds omit 0.19
        /// (0.33…0.21, 0.17, 0), base fractions are divided by the full sequence length, and
        /// sequences shorter than 2 nt score 0. CPAT's prebuilt logistic models were trained on
        /// this feature, so <see cref="PredictCodingProbability"/> uses it.
        /// </summary>
        Cpat,
    }

    /// <summary>Fickett (1982) three-way call for a TESTCODE value.</summary>
    public enum FickettCodingCall
    {
        /// <summary>TESTCODE &lt; 0.74 — predicted non-coding.</summary>
        NonCoding,

        /// <summary>0.74 ≤ TESTCODE &lt; 0.95 — no opinion.</summary>
        NoOpinion,

        /// <summary>TESTCODE ≥ 0.95 — predicted coding.</summary>
        Coding,
    }

    /// <summary>One EMBOSS <c>tcode</c> window: half-open 0-based [Start, End), its TESTCODE and call.</summary>
    public readonly record struct FickettWindow(int Start, int End, double Score, FickettCodingCall Call);

    // Fickett (1982) position-parameter thresholds (max/(min+1) of the per-phase counts),
    // highest first; identical in EMBOSS Etcode.dat and CPAT fickett.py `position_para`.
    private static readonly double[] FickettPositionThresholds =
        { 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0 };

    // Fickett (1982) composition-parameter thresholds (EMBOSS Etcode.dat "Content parameter values").
    private static readonly double[] FickettContentThresholds1982 =
        { 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0.0 };

    // CPAT 3.0.5 fickett.py `content_para`: the 1982 list without 0.19 (9 entries).
    private static readonly double[] FickettContentThresholdsCpat =
        { 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.17, 0.0 };

    // Probability that a region whose parameter falls in each threshold bin is coding,
    // indexed [base A,C,G,T][bin] (Etcode.dat ≡ CPAT position_prob / content_prob).
    private static readonly double[][] FickettPositionProbabilities =
    {
        new[] { 0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22 },
        new[] { 0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23 },
        new[] { 0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08 },
        new[] { 0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09 },
    };

    private static readonly double[][] FickettContentProbabilities =
    {
        new[] { 0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.67, 0.65, 0.81, 0.21 },
        new[] { 0.82, 0.64, 0.51, 0.64, 0.59, 0.59, 0.43, 0.44, 0.39, 0.31 },
        new[] { 0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29 },
        new[] { 0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58 },
    };

    // Per-base weights (A,C,G,T) for the position and composition terms.
    private static readonly double[] FickettPositionWeights = { 0.26, 0.18, 0.31, 0.33 };
    private static readonly double[] FickettContentWeights = { 0.11, 0.12, 0.15, 0.14 };

    // Fickett (1982) decision limits, as applied by EMBOSS tcode_report.
    private const double FickettNonCodingLimit = 0.74;
    private const double FickettCodingLimit = 0.95;

    // EMBOSS tcode.acd: window default 200 (minimum 200), step default 3 (minimum 1).
    private const int FickettMinimumWindow = 200;
    private const int FickettDefaultStep = 3;

    // CPAT fickett_value: `if len(dna) < 2: return 0`.
    private const int CpatFickettMinimumLength = 2;

    private const string FickettBases = "ACGT";

    // Probabilities × weights are products of two-decimal values: exact to 4 decimals.
    private const int FickettScoreDecimals = 4;

    /// <summary>
    /// Computes the Fickett (1982) TESTCODE statistic of a whole sequence with the published
    /// tables (<see cref="FickettTableSet.Fickett1982"/>).
    /// </summary>
    /// <param name="sequence">DNA sequence (case-insensitive).</param>
    /// <returns>The TESTCODE value (≥ 0.95 coding, &lt; 0.74 non-coding; see <see cref="ClassifyFickettScore"/>).</returns>
    public static double CalculateFickettScore(string sequence)
        => CalculateFickettScore(sequence, FickettTableSet.Fickett1982);

    /// <summary>
    /// Computes the Fickett TESTCODE statistic of a whole sequence.
    /// </summary>
    /// <remarks>
    /// For each base b the position parameter is <c>max(b₁,b₂,b₃)/(min(b₁,b₂,b₃)+1)</c> over the
    /// counts of b in codon phases 1–3, and the content parameter is the fraction of b. Each
    /// parameter is mapped to the first threshold bin it reaches (value ≥ threshold), the bin's
    /// coding probability is multiplied by the base's weight, and the eight products are summed
    /// (Fickett 1982, NAR 10(17):5303; EMBOSS <c>tcode_slide</c>; CPAT <c>fickett_value</c>).
    /// <see cref="FickettTableSet.Fickett1982"/> counts composition over A/C/G/T (EMBOSS) and
    /// returns 0 for a sequence with no A/C/G/T (no composition information).
    /// <see cref="FickettTableSet.Cpat"/> reproduces CPAT verbatim (see the enum member).
    /// </remarks>
    /// <param name="sequence">DNA sequence (case-insensitive).</param>
    /// <param name="tables">Table set / convention.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="tables"/> is not a defined value.</exception>
    public static double CalculateFickettScore(string sequence, FickettTableSet tables)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return tables switch
        {
            FickettTableSet.Fickett1982 => FickettScore(sequence.AsSpan(), FickettContentThresholds1982, cpatConvention: false),
            FickettTableSet.Cpat => FickettScore(sequence.AsSpan(), FickettContentThresholdsCpat, cpatConvention: true),
            _ => throw new ArgumentOutOfRangeException(nameof(tables), tables, "Unknown Fickett table set."),
        };
    }

    /// <summary>
    /// Fickett (1982) decision for a TESTCODE value: &lt; 0.74 non-coding, ≥ 0.95 coding,
    /// otherwise no opinion (EMBOSS <c>tcode_report</c>).
    /// </summary>
    public static FickettCodingCall ClassifyFickettScore(double score)
    {
        if (score < FickettNonCodingLimit) return FickettCodingCall.NonCoding;
        if (score >= FickettCodingLimit) return FickettCodingCall.Coding;
        return FickettCodingCall.NoOpinion;
    }

    /// <summary>
    /// Slides a window along the sequence and reports the Fickett (1982) TESTCODE of each
    /// window, as EMBOSS <c>tcode</c>: windows start at 0, 0+step, … while
    /// <c>start + window &lt; length</c> (EMBOSS's strict loop bound, so the window that would end
    /// exactly at the sequence end is not scored), using the <see cref="FickettTableSet.Fickett1982"/> tables.
    /// </summary>
    /// <param name="sequence">DNA sequence (case-insensitive).</param>
    /// <param name="windowSize">Window length; at least 200 (EMBOSS / Fickett minimum).</param>
    /// <param name="stepSize">Window step (default 3, minimum 1).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">window &lt; 200 or step &lt; 1.</exception>
    public static IReadOnlyList<FickettWindow> ScanFickettScore(
        string sequence,
        int windowSize = FickettMinimumWindow,
        int stepSize = FickettDefaultStep)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, FickettMinimumWindow);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(stepSize);

        var windows = new List<FickettWindow>();
        for (int pos = 0; pos + windowSize < sequence.Length; pos += stepSize)
        {
            double score = FickettScore(sequence.AsSpan(pos, windowSize), FickettContentThresholds1982, cpatConvention: false);
            windows.Add(new FickettWindow(pos, pos + windowSize, score, ClassifyFickettScore(score)));
        }

        return windows;
    }

    private static double FickettScore(ReadOnlySpan<char> sequence, double[] contentThresholds, bool cpatConvention)
    {
        if (cpatConvention && sequence.Length < CpatFickettMinimumLength) return 0;

        // counts[phase, base]
        Span<int> counts = stackalloc int[3 * 4];
        for (int i = 0; i < sequence.Length; i++)
        {
            int b = FickettBases.IndexOf(char.ToUpperInvariant(sequence[i]));
            if (b >= 0) counts[(i % 3) * 4 + b]++;
        }

        int acgt = 0;
        for (int k = 0; k < counts.Length; k++) acgt += counts[k];
        if (!cpatConvention && acgt == 0) return 0;

        // CPAT divides by len(dna) (all characters); EMBOSS by the A/C/G/T total.
        double denominator = cpatConvention ? sequence.Length : acgt;

        Span<double> contentTerms = stackalloc double[4];
        Span<double> positionTerms = stackalloc double[4];
        for (int b = 0; b < 4; b++)
        {
            int c0 = counts[b], c1 = counts[4 + b], c2 = counts[8 + b];
            double content = (c0 + c1 + c2) / denominator;
            double position = Math.Max(c0, Math.Max(c1, c2)) / (Math.Min(c0, Math.Min(c1, c2)) + 1.0);

            contentTerms[b] = FickettContentProbabilities[b][FickettBin(content, contentThresholds)] * FickettContentWeights[b];
            positionTerms[b] = FickettPositionProbabilities[b][FickettBin(position, FickettPositionThresholds)] * FickettPositionWeights[b];
        }

        // Summation order of the references (bit-exact): EMBOSS testcode = p1..p4 (position)
        // + p5..p8 (content); CPAT adds the four content terms first, then the position terms.
        Span<double> first = cpatConvention ? contentTerms : positionTerms;
        Span<double> second = cpatConvention ? positionTerms : contentTerms;
        double score = 0;
        for (int b = 0; b < 4; b++) score += first[b];
        for (int b = 0; b < 4; b++) score += second[b];

        // Every probability and weight has two decimals, so the exact TESTCODE is a multiple of
        // 1e-4. The published statistic is returned as that exact value (nearest double), so the
        // 0.74 / 0.95 decision limits are applied without floating-point summation noise. The CPAT
        // feature keeps CPAT's raw float sum (bit-exact with fickett_value).
        return cpatConvention ? score : Math.Round(score, FickettScoreDecimals);
    }

    // First bin whose threshold the value reaches (EMBOSS tcode_index / CPAT look_up_*_prob).
    private static int FickettBin(double value, double[] thresholds)
    {
        for (int i = 0; i < thresholds.Length; i++)
            if (value >= thresholds[i]) return i;
        return thresholds.Length - 1; // unreachable for value ≥ 0 (last threshold is 0)
    }

    #endregion

    #region CPAT coding probability (Wang et al. 2013)

    /// <summary>
    /// A CPAT model: in-frame hexamer frequency tables plus the logistic-regression
    /// coefficients of <c>Label ~ mRNA + ORF + Fickett + Hexamer</c> (CPAT <c>make_logitModel</c>)
    /// and the coding-probability cutoff.
    /// </summary>
    /// <remarks>
    /// CPAT's prebuilt Human/Mouse/Fly/Zebrafish tables and models are distributed under the GNU
    /// GPL (CPAT <c>LICENSE.txt</c>) and are therefore not bundled with this MIT library; load
    /// them with <see cref="ParseHexamerTable"/> and the coefficients printed by
    /// <c>coef(mylogit)</c> (the Human values are listed in the ANNOT-CODING-001 Evidence).
    /// </remarks>
    public sealed class CpatModel
    {
        /// <summary>Creates a model.</summary>
        /// <exception cref="ArgumentNullException">A table is null.</exception>
        /// <exception cref="ArgumentException">A coefficient or the cutoff is not finite, or the cutoff is outside [0, 1].</exception>
        public CpatModel(
            IReadOnlyDictionary<string, double> codingHexamerFrequencies,
            IReadOnlyDictionary<string, double> noncodingHexamerFrequencies,
            double intercept,
            double mrnaLengthCoefficient,
            double orfLengthCoefficient,
            double fickettCoefficient,
            double hexamerCoefficient,
            double codingProbabilityCutoff)
        {
            ArgumentNullException.ThrowIfNull(codingHexamerFrequencies);
            ArgumentNullException.ThrowIfNull(noncodingHexamerFrequencies);
            foreach (double v in new[] { intercept, mrnaLengthCoefficient, orfLengthCoefficient, fickettCoefficient, hexamerCoefficient })
                if (!double.IsFinite(v)) throw new ArgumentException("Logistic coefficients must be finite.");
            if (!double.IsFinite(codingProbabilityCutoff) || codingProbabilityCutoff < 0 || codingProbabilityCutoff > 1)
                throw new ArgumentException("Coding-probability cutoff must lie in [0, 1].", nameof(codingProbabilityCutoff));

            CodingHexamerFrequencies = codingHexamerFrequencies;
            NoncodingHexamerFrequencies = noncodingHexamerFrequencies;
            Intercept = intercept;
            MrnaLengthCoefficient = mrnaLengthCoefficient;
            OrfLengthCoefficient = orfLengthCoefficient;
            FickettCoefficient = fickettCoefficient;
            HexamerCoefficient = hexamerCoefficient;
            CodingProbabilityCutoff = codingProbabilityCutoff;
        }

        /// <summary>In-frame hexamer frequencies of the coding training set.</summary>
        public IReadOnlyDictionary<string, double> CodingHexamerFrequencies { get; }

        /// <summary>In-frame hexamer frequencies of the non-coding training set.</summary>
        public IReadOnlyDictionary<string, double> NoncodingHexamerFrequencies { get; }

        /// <summary>Logistic intercept.</summary>
        public double Intercept { get; }

        /// <summary>Coefficient of the transcript (mRNA) length.</summary>
        public double MrnaLengthCoefficient { get; }

        /// <summary>Coefficient of the ORF length (nt, stop codon included).</summary>
        public double OrfLengthCoefficient { get; }

        /// <summary>Coefficient of the Fickett TESTCODE (CPAT table set).</summary>
        public double FickettCoefficient { get; }

        /// <summary>Coefficient of the hexamer score.</summary>
        public double HexamerCoefficient { get; }

        /// <summary>Coding probability at or above which a transcript is called coding (CPAT TG-ROC cutoff).</summary>
        public double CodingProbabilityCutoff { get; }

        /// <summary>
        /// Logistic response <c>1/(1+e^−η)</c>, η = intercept + Σ coefficient·feature — R
        /// <c>predict(mylogit, newdata, type="response")</c> as called by CPAT <c>coding_prediction</c>.
        /// </summary>
        public double PredictProbability(int mrnaLength, int orfLength, double fickettScore, double hexamerScore)
        {
            double eta = Intercept
                + MrnaLengthCoefficient * mrnaLength
                + OrfLengthCoefficient * orfLength
                + FickettCoefficient * fickettScore
                + HexamerCoefficient * hexamerScore;
            return 1.0 / (1.0 + Math.Exp(-eta));
        }

        /// <summary>
        /// Parses a CPAT hexamer table (<c>make_hexamer_tab</c> output / <c>*_Hexamer.tsv</c>):
        /// whitespace-separated <c>hexamer coding noncoding</c> rows, the header row whose first
        /// field is <c>hexamer</c> skipped, as CPAT's reader. Blank lines are skipped.
        /// </summary>
        /// <exception cref="ArgumentNullException"><paramref name="lines"/> is null.</exception>
        /// <exception cref="FormatException">A row has fewer than three fields or a non-numeric frequency.</exception>
        public static (IReadOnlyDictionary<string, double> Coding, IReadOnlyDictionary<string, double> Noncoding)
            ParseHexamerTable(IEnumerable<string> lines)
        {
            ArgumentNullException.ThrowIfNull(lines);
            var coding = new Dictionary<string, double>(StringComparer.Ordinal);
            var noncoding = new Dictionary<string, double>(StringComparer.Ordinal);
            foreach (string raw in lines)
            {
                string[] fields = raw.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
                if (fields.Length == 0 || fields[0] == "hexamer") continue;
                if (fields.Length < 3)
                    throw new FormatException($"CPAT hexamer row needs 'hexamer coding noncoding': '{raw}'.");
                coding[fields[0]] = double.Parse(fields[1], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture);
                noncoding[fields[0]] = double.Parse(fields[2], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture);
            }

            return (coding, noncoding);
        }
    }

    /// <summary>CPAT ORF-search options (defaults = CPAT 3.0.5 command-line defaults).</summary>
    public sealed record CpatOptions
    {
        /// <summary>Start codons (<c>--start</c>, default ATG). Empty: any non-stop codon opens an ORF.</summary>
        public IReadOnlyList<string> StartCodons { get; init; } = new[] { "ATG" };

        /// <summary>Stop codons (<c>--stop</c>, default TAG, TAA, TGA).</summary>
        public IReadOnlyList<string> StopCodons { get; init; } = new[] { "TAG", "TAA", "TGA" };

        /// <summary>An ORF is kept when its length in nt (stop included) is strictly greater than this (<c>--min-orf</c>, default 75).</summary>
        public int MinOrfLength { get; init; } = 75;

        /// <summary>Number of longest ORF candidates scored (<c>--top-orf</c>, default 5).</summary>
        public int TopOrfCount { get; init; } = 5;

        /// <summary>Also search the reverse complement (<c>--antisense</c>, default false).</summary>
        public bool SearchAntisense { get; init; }
    }

    /// <summary>
    /// One scored CPAT ORF candidate. Coordinates follow CPAT's <c>ORF_info</c> output: 1-based on the
    /// input transcript; for '+' <c>Start</c> ≤ <c>End</c>, for '−' <c>Start</c> &gt; <c>End</c>
    /// (CPAT reports the antisense ORF from its 5' end). <c>Frame</c> is 1–3 on the searched strand.
    /// </summary>
    public sealed record CpatOrfPrediction(
        char Strand,
        int Frame,
        int Start,
        int End,
        int OrfLength,
        string OrfSequence,
        int MrnaLength,
        double FickettScore,
        double HexamerScore,
        double CodingProbability,
        bool IsCoding);

    /// <summary>
    /// CPAT (Wang et al. 2013, CPAT 3.0.5): finds ORF candidates, scores each with ORF length,
    /// transcript length, Fickett TESTCODE (<see cref="FickettTableSet.Cpat"/>) and the hexamer
    /// score (<see cref="CalculateCodingPotential"/>), and returns their logistic coding probabilities.
    /// </summary>
    /// <remarks>
    /// ORF search is CPAT <c>ORFFinder</c>: per frame (0,1,2 of the sense strand, then of the reverse
    /// complement when <see cref="CpatOptions.SearchAntisense"/>), codons are read left to right; a start
    /// codon opens an ORF that runs to and includes the next in-frame stop codon, or to the last complete
    /// codon when none follows; scanning resumes after it, so nested starts are not separate ORFs. ORFs
    /// longer than <see cref="CpatOptions.MinOrfLength"/> nt are ranked by length (stable, sense frames
    /// first) and the top <see cref="CpatOptions.TopOrfCount"/> are scored. This differs from
    /// <see cref="FindOrfs"/> (every nested ATG/GTG/TTG start, stop-terminated ORFs only), so it is a
    /// separate scanner. CPAT: Wang L et al. (2013) NAR 41(6):e74, https://doi.org/10.1093/nar/gkt006.
    /// </remarks>
    /// <param name="transcript">Transcript DNA sequence (case-insensitive; use T, not U).</param>
    /// <param name="model">Hexamer tables + logistic model.</param>
    /// <param name="options">ORF-search options (CPAT defaults when null).</param>
    /// <returns>The scored candidates in CPAT order (longest first); empty when no ORF qualifies.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="transcript"/> or <paramref name="model"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">MinOrfLength &lt; 0 or TopOrfCount &lt; 1.</exception>
    public static IReadOnlyList<CpatOrfPrediction> PredictCodingProbability(
        string transcript,
        CpatModel model,
        CpatOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(transcript);
        ArgumentNullException.ThrowIfNull(model);
        options ??= new CpatOptions();
        ArgumentNullException.ThrowIfNull(options.StartCodons);
        ArgumentNullException.ThrowIfNull(options.StopCodons);
        ArgumentOutOfRangeException.ThrowIfNegative(options.MinOrfLength);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(options.TopOrfCount);

        string seq = transcript.ToUpperInvariant();
        var starts = new HashSet<string>(options.StartCodons.Select(c => c.ToUpperInvariant()), StringComparer.Ordinal);
        var stops = new HashSet<string>(options.StopCodons.Select(c => c.ToUpperInvariant()), StringComparer.Ordinal);

        var candidates = new List<(char Strand, int Frame, int Start, int End, string Sequence)>();
        for (int frame = 0; frame < 3; frame++)
            FindCpatOrfsInFrame(seq, frame, '+', starts, stops, options.MinOrfLength, candidates);
        if (options.SearchAntisense)
        {
            string reverse = DnaSequence.GetReverseComplementString(seq);
            for (int frame = 0; frame < 3; frame++)
                FindCpatOrfsInFrame(reverse, frame, '-', starts, stops, options.MinOrfLength, candidates);
        }

        int mrnaLength = seq.Length;
        var predictions = new List<CpatOrfPrediction>();
        // Python sorted(..., key=length, reverse=True) is stable: OrderByDescending is too.
        foreach (var orf in candidates.OrderByDescending(o => o.Sequence.Length).Take(options.TopOrfCount))
        {
            double fickett = CalculateFickettScore(orf.Sequence, FickettTableSet.Cpat);
            double hexamer = CalculateCodingPotential(orf.Sequence, model.CodingHexamerFrequencies, model.NoncodingHexamerFrequencies);
            int orfLength = orf.Sequence.Length;
            double probability = model.PredictProbability(mrnaLength, orfLength, fickett, hexamer);

            // CPAT reports 1-based coordinates: '+' → (start+1, end); '−' → (L − start, L − end + 1).
            (int start, int end) = orf.Strand == '+'
                ? (orf.Start + 1, orf.End)
                : (mrnaLength - orf.Start, mrnaLength - orf.End + 1);

            predictions.Add(new CpatOrfPrediction(
                orf.Strand, orf.Frame, start, end, orfLength, orf.Sequence, mrnaLength,
                fickett, hexamer, probability, probability >= model.CodingProbabilityCutoff));
        }

        return predictions;
    }

    /// <summary>
    /// CPAT best ORF of a transcript (<c>ORF_prob.best.tsv</c>): among the scored candidates the
    /// one with the highest coding probability (CPAT <c>--best-orf p</c>, default) or the longest
    /// (<c>--best-orf l</c>); ties keep the earlier candidate. Null when no ORF qualifies.
    /// </summary>
    public static CpatOrfPrediction? PredictBestCodingOrf(
        string transcript,
        CpatModel model,
        CpatOptions? options = null,
        bool selectByLength = false)
    {
        CpatOrfPrediction? best = null;
        foreach (var candidate in PredictCodingProbability(transcript, model, options))
        {
            if (best is null
                || (selectByLength ? candidate.OrfLength > best.OrfLength : candidate.CodingProbability > best.CodingProbability))
                best = candidate;
        }

        return best;
    }

    // CPAT find_orfs.ORFFinder.run_one for one frame of one strand.
    private static void FindCpatOrfsInFrame(
        string seq,
        int frame,
        char strand,
        HashSet<string> starts,
        HashSet<string> stops,
        int minOrfLength,
        List<(char Strand, int Frame, int Start, int End, string Sequence)> result)
    {
        int i = frame;
        while (i + CodonLength <= seq.Length)
        {
            string codon = seq.Substring(i, CodonLength);
            int orfStart = i;
            i += CodonLength;

            bool opens = starts.Count > 0 ? starts.Contains(codon) : !stops.Contains(codon);
            if (!opens) continue;

            // Read to the next stop (inclusive) or to the last complete codon.
            int lastCodon = orfStart;
            while (i + CodonLength <= seq.Length)
            {
                lastCodon = i;
                bool isStop = stops.Contains(seq.Substring(i, CodonLength));
                i += CodonLength;
                if (isStop) break;
            }

            int orfEnd = lastCodon + CodonLength;
            if (orfEnd - orfStart > minOrfLength)
                result.Add((strand, frame + 1, orfStart, orfEnd, seq.Substring(orfStart, orfEnd - orfStart)));
        }
    }

    #endregion

    // Tandem-repeat unit upper bound. Microsatellite (STR) motifs are 1-6 bp and
    // minisatellite (VNTR) motifs 10-60 bp per Wikipedia "Tandem repeat" (cites
    // Duitama et al. 2014, Nucleic Acids Res 42(9):5728-5741). 60 bp covers both.
    private const int MaxTandemUnitLength = 60;

    // Inverted-repeat arm upper bound: EMBOSS palindrome -maxpallen default (palindrome.acd) and
    // IUPACpal -M max_len default (Hampson et al. 2021, README) are both 100. Semantics are those of
    // EMBOSS -maxpallen: a longer stem is not reported and is not split into shorter pieces.
    private const int MaxInvertedArmLength = 100;

    // Maximum gap (loop/spacer |G|) between inverted-repeat arms: EMBOSS palindrome -gaplimit default
    // (palindrome.acd) and IUPACpal -g max_gap default (README) are both 100; |G| >= 0 (WGW^R, Hampson 2021).
    private const int DefaultMaxInvertedGap = 100;

    /// <summary>
    /// Finds repetitive elements in a DNA sequence: exact tandem repeats (head-to-tail
    /// consecutive copies of a primitive motif) and exact inverted repeats (a sequence
    /// followed downstream by its reverse complement, form WGW̄ᴿ).
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Tandem repeats</b> are the maximal perfect runs of every period p = 1..60 (Kolpakov &amp; Kucherov 1999;
    /// the MISA / pytrf convention) computed by the canonical
    /// <see cref="Seqeron.Genomics.Analysis.RepeatFinder.FindMicrosatellites(string,int,int,int)"/>: each run is
    /// reported once, at its left end, with its primitive A/C/G/T unit and its complete copies only
    /// (<c>end − start = copies·p</c>, copies ≥ <paramref name="minCopies"/>, span ≥ <paramref name="minRepeatLength"/>);
    /// rotations of a run, non-primitive units and runs of N are not reported. Runs of different periods may overlap.
    /// </para>
    /// <para>
    /// <b>Inverted repeats</b> are the maximal exact stems of
    /// <see cref="Seqeron.Genomics.Analysis.RepeatFinder.FindInvertedRepeats(string,int,int,int,int,int,bool)"/>
    /// with arm ≥ <paramref name="minRepeatLength"/>, arm ≤ max(100, that minimum), loop 0..100 — the result
    /// set of EMBOSS <c>palindrome -minpallen m -maxpallen 100 -gaplimit 100 -nummismatches 0 -overlap Y</c>
    /// (the EMBOSS and IUPACpal default limits). Only A/C/G/T pair (A·T, C·G); sub-stems inside a longer stem are not
    /// reported. For inverted repeats <paramref name="minRepeatLength"/> is the minimum <b>arm</b> length.
    /// </para>
    /// <para>Tandem repeats come first (by period, then start), then inverted repeats (by left-arm start, then
    /// right-arm start).</para>
    /// </remarks>
    /// <param name="dnaSequence">DNA sequence to scan (case-insensitive).</param>
    /// <param name="minRepeatLength">Minimum total span (bp) of a tandem repeat and minimum arm length of an inverted repeat.</param>
    /// <param name="minCopies">Minimum number of adjacent copies for a tandem repeat (>= 2).</param>
    /// <returns>
    /// Tuples of (start, end, type, sequence) where start is 0-based inclusive,
    /// end is exclusive, and type is "tandem_repeat" or "inverted_repeat". The inverted-repeat sequence is
    /// <c>leftArm + rightArm</c> (loop 0) or <c>leftArm[loop]rightArm</c>.
    /// </returns>
    /// <exception cref="ArgumentNullException">dnaSequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">minRepeatLength &lt; 1 or minCopies &lt; 2.</exception>
    public static IEnumerable<(int start, int end, string type, string sequence)> FindRepetitiveElements(
        string dnaSequence,
        int minRepeatLength = 10,
        int minCopies = 2)
    {
        ArgumentNullException.ThrowIfNull(dnaSequence);
        // A tandem repeat requires "two or more" adjacent copies (Benson 1999).
        ArgumentOutOfRangeException.ThrowIfLessThan(minCopies, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(minRepeatLength, 1);

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

        foreach (var run in Seqeron.Genomics.Analysis.RepeatFinder.FindMicrosatellites(
                     sequence, 1, MaxTandemUnitLength, minCopies))
        {
            if (run.TotalLength < minRepeatLength)
                continue;
            yield return (run.Position, run.Position + run.TotalLength, "tandem_repeat",
                sequence.Substring(run.Position, run.TotalLength));
        }

        foreach (var stem in Seqeron.Genomics.Analysis.RepeatFinder.FindInvertedRepeats(
                     sequence, minRepeatLength, DefaultMaxInvertedGap, minLoopLength: 0,
                     maxArmLength: Math.Max(MaxInvertedArmLength, minRepeatLength)))
        {
            string composite = stem.LoopLength == 0
                ? stem.LeftArm + stem.RightArm
                : stem.LeftArm + "[" + stem.Loop + "]" + stem.RightArm;
            yield return (stem.LeftArmStart, stem.RightArmStart + stem.ArmLength, "inverted_repeat", composite);
        }
    }

    /// <summary>
    /// Classifies a repeat sequence against a repeat-element library, mirroring the
    /// RepeatMasker convention of screening the query for occurrences of known
    /// library elements on <b>both strands</b> (RepeatMasker reports matches in "+" or
    /// "C" = complement orientation) and assigning the class of the best (longest) match.
    /// Matching is exact containment of a library element, or of its reverse complement,
    /// within the query; for scored (Smith-Waterman-Gotoh) homology use
    /// <see cref="ClassifyRepeat(string, IReadOnlyDictionary{string, string}, ScoringMatrix, int)"/>.
    /// When no library entry matches, the query falls back to "Simple_repeat" when it is a
    /// perfect tandem repeat (≥ 2 complete copies, whole query) of a primitive A/C/G/T motif
    /// of 1-6 bp, else "Unknown".
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

        // Best library match: the longest known element that occurs within the query on either
        // strand. Matching is one-directional (element ⊆ query): a query that merely is a
        // *substring of* a longer consensus is not an occurrence of that element.
        string? bestClass = null;
        int bestLen = -1;
        foreach (var entry in repeatDb)
        {
            string element = entry.Key.ToUpperInvariant();
            if (element.Length == 0 || element.Length <= bestLen) continue;
            if (query.Contains(element, StringComparison.Ordinal) ||
                query.Contains(DnaSequence.GetReverseComplementString(element), StringComparison.Ordinal))
            {
                bestLen = element.Length;
                bestClass = entry.Value;
            }
        }

        return bestClass ?? SimpleRepeatOrUnknown(query);
    }

    /// <summary>
    /// Classifies a repeat sequence by scored homology to a repeat-element library — the RepeatMasker /
    /// cross_match model: Smith-Waterman-Gotoh local alignment of the query against every library element on
    /// both strands ("+" and "C"), keeping the best-scoring match whose score reaches
    /// <paramref name="minScore"/> (RepeatMasker <c>-cutoff</c>). Alignment is the canonical
    /// <see cref="Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffine(string, string, ScoringMatrix?)"/>
    /// (gap of length k scores <c>GapOpen + k·GapExtend</c>). Ties on score go to the longer element, then to the
    /// first in enumeration order. Without a match the fallback is the same as
    /// <see cref="ClassifyRepeat(string, IReadOnlyDictionary{string, string})"/> (Simple_repeat / Unknown).
    /// </summary>
    /// <remarks>
    /// RepeatMasker's own nucleotide matrices (GC-dependent, transition/transversion-specific) and its
    /// complexity-adjusted score are not reproduced: <see cref="ScoringMatrix"/> is a match/mismatch scheme, and the
    /// scoring scheme and cutoff are the caller's choice (no default is implied).
    /// </remarks>
    /// <param name="sequence">Repeat sequence to classify (case-insensitive).</param>
    /// <param name="repeatDb">Library: element sequence → class label.</param>
    /// <param name="scoring">Match / mismatch / affine gap scores.</param>
    /// <param name="minScore">Minimum local-alignment score for a library match (≥ 1).</param>
    /// <returns>The class label of the best-scoring library entry, or a fallback class.</returns>
    /// <exception cref="ArgumentNullException">sequence, repeatDb or scoring is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">minScore &lt; 1.</exception>
    public static string ClassifyRepeat(
        string sequence,
        IReadOnlyDictionary<string, string> repeatDb,
        ScoringMatrix scoring,
        int minScore)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(repeatDb);
        ArgumentNullException.ThrowIfNull(scoring);
        ArgumentOutOfRangeException.ThrowIfLessThan(minScore, 1);

        if (sequence.Length == 0)
            return SimpleRepeatClass;

        string query = sequence.ToUpperInvariant();

        string? bestClass = null;
        int bestScore = minScore - 1;
        int bestLen = -1;
        foreach (var entry in repeatDb)
        {
            string element = entry.Key.ToUpperInvariant();
            if (element.Length == 0) continue;
            int score = Math.Max(
                Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffine(query, element, scoring).Score,
                Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffine(
                    query, DnaSequence.GetReverseComplementString(element), scoring).Score);
            if (score > bestScore || (score == bestScore && bestClass is not null && element.Length > bestLen))
            {
                bestScore = score;
                bestLen = element.Length;
                bestClass = entry.Value;
            }
        }

        return bestClass ?? SimpleRepeatOrUnknown(query);
    }

    // RepeatMasker output class labels (Smit/Hubley/Green; Repbase Update).
    private const string SimpleRepeatClass = "Simple_repeat";
    private const string UnknownClass = "Unknown";

    // STR/microsatellite motifs are 1-6 bp (Tóth et al. 2000).
    private const int MaxSimpleRepeatMotif = 6;

    /// <summary>
    /// "Simple_repeat" when the whole (upper-case) query is one perfect run of ≥ 2 complete copies of a primitive
    /// A/C/G/T motif of 1-6 bp (canonical <see cref="Seqeron.Genomics.Analysis.RepeatFinder.FindMicrosatellites(string,int,int,int)"/>),
    /// else "Unknown".
    /// </summary>
    private static string SimpleRepeatOrUnknown(string query) =>
        Seqeron.Genomics.Analysis.RepeatFinder.FindMicrosatellites(query, 1, MaxSimpleRepeatMotif, 2)
            .Any(run => run.Position == 0 && run.TotalLength == query.Length)
            ? SimpleRepeatClass
            : UnknownClass;

    /// <summary>
    /// Counts the complete, unambiguous in-frame codons of a coding sequence (frame 0, non-overlapping
    /// triplets). Delegates to the canonical counter <see cref="Translator.CountCodons(string?, int)"/>:
    /// case-insensitive DNA or RNA (U read as T, CodonW <c>ident_codon</c>), keys in upper-case DNA
    /// spelling, triplets containing any other symbol (N, IUPAC codes, gaps) skipped without shifting the
    /// frame, a trailing partial triplet ignored (EMBOSS <c>ajCodSetTripletsS</c>).
    /// </summary>
    /// <param name="dnaSequence">Coding DNA or RNA sequence.</param>
    /// <returns>Codon (upper-case DNA) → occurrence count; codons that do not occur are absent.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="dnaSequence"/> is null.</exception>
    public static IReadOnlyDictionary<string, int> GetCodonUsage(string dnaSequence)
    {
        ArgumentNullException.ThrowIfNull(dnaSequence);
        return Translator.CountCodons(dnaSequence);
    }

    // RSCU is reported over sense codons only; the termination family ('*' in the genetic-code
    // table) is dropped (CodonU `rscu` iterates Biopython forward_table; codonbias ignore_stop=True).
    private const char StopCodonSymbol = '*';

    // A codon is exactly three nucleotides.
    private const int CodonLength = 3;

    /// <summary>
    /// Computes the Relative Synonymous Codon Usage (RSCU) over a set of coding
    /// sequences using the Standard genetic code (NCBI translation table 1).
    /// </summary>
    /// <remarks>
    /// RSCU for codon <c>j</c> of amino acid <c>i</c> is the observed count divided by the count
    /// expected under uniform synonymous usage, <c>RSCU = x_(i,j) / ((1/n_i) Σ_j x_(i,j)) = n_i · x_(i,j) / Σ_j x_(i,j)</c>,
    /// where <c>n_i</c> is the number of synonymous codons for amino acid <c>i</c> (Sharp, Tuohy &amp;
    /// Mosurski 1986, Nucleic Acids Res. 14(13):5125–5143; CodonW <c>rscu_usage_out</c>). An RSCU of
    /// 1.0 indicates no codon bias; values above 1 indicate a preferred codon.
    /// Codon counts are pooled across all input sequences before the RSCU is computed (CodonW
    /// <c>-totals</c>, CodonU <c>rscu</c>). Stop codons are excluded; an unobserved synonymous family
    /// yields RSCU 0.0 for each of its codons (CodonW prints 0.000; the CAI 0.5 pseudocount is not applied).
    /// Computed by the canonical <c>CodonUsageAnalyzer.CalculateRscu(counts, code)</c> (MolTools).
    /// </remarks>
    /// <param name="codingSequences">The coding (CDS) sequences, DNA or RNA, case-insensitive. Each is read
    /// in frame from position 0 in steps of three as by <see cref="GetCodonUsage(string)"/> (U read as T,
    /// triplets with other symbols skipped without shifting the frame, partial trailing codon ignored);
    /// null or empty elements are skipped.</param>
    /// <returns>A dictionary from each sense codon (uppercase DNA) to its RSCU value.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="codingSequences"/> is null.</exception>
    public static IReadOnlyDictionary<string, double> GetCodonUsage(IEnumerable<string> codingSequences)
        => GetCodonUsage(codingSequences, GeneticCode.Standard);

    /// <summary>
    /// Computes the Relative Synonymous Codon Usage (RSCU) over a set of coding
    /// sequences using the supplied genetic code ("RSCU values are genetic code dependent", CodonW).
    /// </summary>
    /// <remarks>
    /// Synonymous families are the amino-acid groups of <see cref="GeneticCode.CodonTable"/>; codons that
    /// are context-dependent stops in NCBI tables 27/28/31 belong to the amino acid they encode
    /// (Biopython <c>forward_table</c>). See <see cref="GetCodonUsage(IEnumerable{string})"/>.
    /// </remarks>
    /// <param name="codingSequences">The coding (CDS) sequences (see the Standard-code overload for
    /// framing and alphabet rules).</param>
    /// <param name="code">The genetic code defining synonymous codon families.</param>
    /// <returns>A dictionary from each sense codon (uppercase DNA) to its RSCU value.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="codingSequences"/> or
    /// <paramref name="code"/> is null.</exception>
    /// <exception cref="OverflowException">The pooled count of one codon exceeds <see cref="int.MaxValue"/>.</exception>
    public static IReadOnlyDictionary<string, double> GetCodonUsage(
        IEnumerable<string> codingSequences,
        GeneticCode code)
    {
        ArgumentNullException.ThrowIfNull(codingSequences);
        ArgumentNullException.ThrowIfNull(code);

        // Pool observed counts across all coding sequences (canonical Core codon counter).
        var counts = new Dictionary<string, int>();
        foreach (string sequence in codingSequences)
        {
            foreach (var (codon, count) in Translator.CountCodons(sequence))
                counts[codon] = checked(counts.GetValueOrDefault(codon) + count);
        }

        // Canonical RSCU (n·x/Σx over the code's synonymous families, 0 for an absent family),
        // restricted to the sense codons.
        var rscu = Seqeron.Genomics.MolTools.CodonUsageAnalyzer.CalculateRscu(counts, code);
        foreach (var (codonRna, aminoAcid) in code.CodonTable)
        {
            if (aminoAcid == StopCodonSymbol)
                rscu.Remove(codonRna.Replace('U', 'T'));
        }

        return rscu;
    }
}
