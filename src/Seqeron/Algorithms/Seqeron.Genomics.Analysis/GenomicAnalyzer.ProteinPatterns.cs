using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Analysis;

/// <summary>A hit of <see cref="GenomicAnalyzer.FindProteinPatternHits(string, string, CancellationToken)"/>: the 0-based start and the length of the matching window.</summary>
/// <param name="Position">0-based start (EMBOSS <c>fuzzpro</c> prints Position + 1).</param>
/// <param name="Length">Window length (a range such as <c>x(2,4)</c> yields hits of several lengths).</param>
public readonly record struct ProteinPatternHit(int Position, int Length)
{
    /// <summary>Exclusive end of the window (0-based).</summary>
    public int End => Position + Length;
}

/// <summary>A hit of <see cref="GenomicAnalyzer.FindTranslatedPatternHits(string, string, string, int)"/> (EMBOSS <c>fuzztran</c>).</summary>
/// <param name="Frame">Translation frame: 1, 2, 3 forward, -1, -2, -3 reverse (EMBOSS transeq convention).</param>
/// <param name="ProteinStart">1-based first residue in the translated frame (the program's <c>*start</c> note).</param>
/// <param name="ProteinEnd">1-based last residue (<c>*end</c>).</param>
/// <param name="Translation">The matched residues (<c>*translation</c>).</param>
/// <param name="Position">0-based start of the codons on the forward strand of the input.</param>
/// <param name="Length">Length in nucleotides (3 per residue).</param>
/// <param name="EmbossStart">Forward frames: the first coordinate the program prints (<c>fuzztran_SourceFeature</c>: <c>1 + frame + 3·start − 4</c>), = <see cref="Position"/> + 1.
/// Null for the reverse frames: EMBOSS 6.6.0 prints meaningless numbers for them (negative, and dependent on state of the program that is not part of the search), which are not reproduced — use <see cref="Position"/>.</param>
/// <param name="EmbossEnd">Forward frames: the second coordinate the program prints (<c>1 + frame + 3·end − 2</c>), = <see cref="Position"/> + <see cref="Length"/>; null for the reverse frames.</param>
public readonly record struct TranslatedPatternHit(
    int Frame, int ProteinStart, int ProteinEnd, string Translation, int Position, int Length, int? EmbossStart, int? EmbossEnd)
{
    /// <summary>True for the reverse frames (-1..-3).</summary>
    public bool IsReverseStrand => Frame < 0;
}

public static partial class GenomicAnalyzer
{
    #region fuzzpro / fuzztran (protein patterns)

    /// <summary>
    /// EMBOSS <c>fuzzpro</c> (6.6.0): every window of a protein sequence that matches a <see cref="ProteinPattern"/> (PROSITE-like: letters, <c>x</c>, <c>[ILV]</c>, <c>{P}</c>,
    /// <c>x(3)</c>, <c>x(2,4)</c>, <c>&lt;</c> / <c>&gt;</c>), exact (the binary has no <c>-pmismatch</c>).
    /// </summary>
    /// <remarks>
    /// Like the binary, every matching window length of a range is reported (overlapping hits included); a start-anchored (<c>&lt;</c>) range pattern reports only its longest match
    /// (the regular-expression search class of <c>embpat.c</c>). The sequence is case-insensitive; a single trailing stop <c>*</c> is not part of it (the binary reads it off). Hits are
    /// ordered by position, then length. Verified against fuzzpro 6.6.0 on 900 random syntax-tree patterns over sequences with X, B, Z, U, J, O, * and lower case
    /// (docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md).
    /// </remarks>
    /// <param name="protein">The protein sequence.</param>
    /// <param name="pattern">The pattern text (<see cref="ProteinPattern.Parse"/>).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The pattern is invalid.</exception>
    public static IReadOnlyList<ProteinPatternHit> FindProteinPatternHits(string protein, string pattern, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(protein);
        return FindProteinPatternHits(protein, ProteinPattern.Parse(pattern), cancellationToken);
    }

    /// <summary><see cref="FindProteinPatternHits(string, string, CancellationToken)"/> with an already parsed pattern.</summary>
    /// <param name="protein">The protein sequence.</param>
    /// <param name="pattern">The parsed pattern.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static IReadOnlyList<ProteinPatternHit> FindProteinPatternHits(string protein, ProteinPattern pattern, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(protein);
        ArgumentNullException.ThrowIfNull(pattern);
        return ProteinPatternHitsCore(protein.EndsWith('*') ? protein[..^1] : protein, pattern, cancellationToken);
    }

    private static List<ProteinPatternHit> ProteinPatternHitsCore(string text, ProteinPattern pattern, CancellationToken cancellationToken)
    {
        var bits = new uint[text.Length];
        for (int i = 0; i < bits.Length; i++)
            bits[i] = ProteinPattern.ResidueBit(text[i]);

        var hits = new HashSet<(int Position, int Length)>();
        foreach (var expansion in pattern.Expand())
        {
            cancellationToken.ThrowIfCancellationRequested();
            int length = expansion.Length;
            for (int pos = 0; pos + length <= bits.Length; pos++)
            {
                if (pattern.AnchoredAtStart && pos != 0)
                    break;
                if (pattern.AnchoredAtEnd && pos + length != bits.Length)
                    continue;
                bool ok = true;
                for (int j = 0; j < length && ok; j++)
                    ok = (bits[pos + j] & expansion[j]) != 0;
                if (ok)
                    hits.Add((pos, length));
            }
        }

        // regex class (embPatMatchFindAll with nterm): a start-anchored range pattern returns only the longest match at the anchor
        if (pattern.AnchoredAtStart && pattern.HasRange && hits.Count > 1)
        {
            int longest = hits.Max(h => h.Length);
            hits.RemoveWhere(h => h.Length != longest);
        }

        return hits.OrderBy(h => h.Position).ThenBy(h => h.Length).Select(h => new ProteinPatternHit(h.Position, h.Length)).ToList();
    }

    /// <summary>
    /// EMBOSS <c>fuzztran</c> (6.6.0): translates a nucleotide sequence in the chosen frames with an EMBOSS genetic-code table and searches every translation with a
    /// <see cref="ProteinPattern"/> (<see cref="FindProteinPatternHits(string, ProteinPattern, CancellationToken)"/>, without reading a trailing stop off).
    /// </summary>
    /// <remarks>
    /// Frames <c>1</c>, <c>2</c>, <c>3</c> (forward), <c>-1</c>, <c>-2</c>, <c>-3</c> (reverse), <c>F</c> (forward three), <c>R</c> (reverse three) or <c>6</c> (all) — the program's
    /// <c>-frame</c> — are translated as <c>transeq</c> does (<see cref="Translator.TranslateSixFrames(string, GeneticCode?, SixFrameConvention)"/> with
    /// <see cref="SixFrameConvention.EmbossTranseq"/>; only the complete codons are searched, the dangling-base residue of transeq is dropped) and the sequence is upper-cased first. <paramref name="table"/> is the program's <c>-table</c> menu (0 =
    /// standard; otherwise the NCBI table number). Hits are listed per frame in the order 1, 2, 3, -1, -2, -3 (those requested), by protein position then length. Each hit carries the
    /// forward-strand nucleotide location (<see cref="TranslatedPatternHit.Position"/>, <see cref="TranslatedPatternHit.Length"/>) and the numbers the binary prints
    /// (<see cref="TranslatedPatternHit.EmbossStart"/> / <see cref="TranslatedPatternHit.EmbossEnd"/>: the program's 1-based coordinates, exact for the forward frames; for the reverse frames EMBOSS 6.6.0
    /// prints meaningless numbers from a defect of its <c>fuzztran_SourceFeature</c> and they are not reproduced). Verified against fuzztran 6.6.0 on 700 random cases (docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md).
    /// </remarks>
    /// <param name="nucleotides">The nucleotide sequence (IUPAC codes allowed).</param>
    /// <param name="pattern">The protein pattern text.</param>
    /// <param name="frames">The frame selection: 1, 2, 3, -1, -2, -3, F, R or 6 (case-insensitive).</param>
    /// <param name="table">EMBOSS genetic code (0 standard, otherwise an NCBI table number).</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The pattern or the frame selection is invalid, or the table is unknown.</exception>
    public static IReadOnlyList<TranslatedPatternHit> FindTranslatedPatternHits(string nucleotides, string pattern, string frames = "1", int table = 0)
    {
        ArgumentNullException.ThrowIfNull(nucleotides);
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(frames);
        var parsed = ProteinPattern.Parse(pattern);
        int[] selected = frames.Trim().ToUpperInvariant() switch
        {
            "F" => [1, 2, 3],
            "R" => [-1, -2, -3],
            "6" => [1, 2, 3, -1, -2, -3],
            "1" => [1],
            "2" => [2],
            "3" => [3],
            "-1" => [-1],
            "-2" => [-2],
            "-3" => [-3],
            _ => throw new ArgumentException("frames must be 1, 2, 3, -1, -2, -3, F, R or 6 (fuzztran -frame).", nameof(frames)),
        };
        var code = GeneticCode.GetByTableNumberOrThrow(table == 0 ? 1 : table, nameof(table));
        string text = nucleotides.ToUpperInvariant();
        var translations = Translator.TranslateSixFrames(text, code, SixFrameConvention.EmbossTranseq);
        int length = text.Length;
        var result = new List<TranslatedPatternHit>();
        foreach (int frame in selected)
        {
            // fuzztran translates only complete codons (the dangling-base residue transeq appends is not searched) and searches the translation as it is (a stop '*' stays)
            string protein = translations[frame].Sequence;
            int complete = Math.Max(0, (length - (Math.Abs(frame) - 1)) / 3);
            if (protein.Length > complete)
                protein = protein[..complete];
            foreach (var hit in ProteinPatternHitsCore(protein, parsed, default))
            {
                int s = hit.Position + 1, e = hit.Position + hit.Length;
                int framenum = Math.Abs(frame);
                int position = frame > 0 ? framenum + 3 * s - 4 : length - (framenum - 1) - 3 * e;
                int? embossStart = frame > 0 ? position + 1 : null;
                int? embossEnd = frame > 0 ? position + 3 * hit.Length : null;
                result.Add(new TranslatedPatternHit(frame, s, e, protein.Substring(hit.Position, hit.Length), position, 3 * hit.Length, embossStart, embossEnd));
            }
        }

        return result;
    }

    #endregion
}
