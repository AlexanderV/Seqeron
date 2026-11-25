// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossMotifs, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.EmbossMotifs.Tools;

public partial class EmbossMotifsTools
{
    [McpServerTool(Name = "find_pattern_hits", Title = "Genome — fuzznuc Pattern Search", ReadOnly = true)]
    [Description("Search a DNA sequence with an EMBOSS fuzznuc pattern (the -pattern language): IUPAC letters, classes [ACG], exclusions {AG}, repeats N(3), length ranges A(2,4) (every matching window length is reported), '-' separators, the anchors '<' (start of the sequence) and '>' (end), a terminating '.'; e.g. [CG](5)TG{A}N(1,5)C. maxMismatches = fuzznuc -pmismatch (Hamming, substitutions only; with a range the extra repeats beyond n of X(n,m) must match exactly, as in fuzznuc); bothStrands = -complement (the reverse-strand pattern is the reversed, complemented one with swapped anchors; positions are forward-strand starts). Give one pattern, or patternFile = the text of a fuzznuc pattern file (one pattern per line, or FASTA '>name' + pattern) for several named patterns. Verified against EMBOSS fuzznuc 6.6.0 on 3300 random patterns. By default a symbol other than A/C/G/T in the sequence matches nothing and costs a mismatch; with sequenceAmbiguity = true the sequence is read as the fuzznuc binary reads it (1350 more oracle rows with ambiguity codes in the sequence): a pattern IUPAC letter matches every sequence symbol whose base set lies within the letter's (pattern R matches A, G, R; sequence R matches pattern R, V, D, N but not A, M, Y), N matches every symbol, a class [ACG] matches its plain bases only, an exclusion {AG} matches every symbol not listed (ambiguity codes included), U reads as T, and the reverse strand is the IUPAC reverse complement of the sequence (the reported patterns are the plain-base normal form). Positions 0-based; hits ordered by position, length, forward before reverse. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static FindPatternHitsResult FindPatternHits(
        [Description("Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive).")] string sequence,
        [Description("One fuzznuc pattern (omit when patternFile is given).")] string? pattern = null,
        [Description("Text of a fuzznuc pattern file: one pattern per line, or FASTA style (>name then the pattern); several patterns are reported under their names (unnamed: pattern1, pattern2, ...).")] string? patternFile = null,
        [Description("Maximum substitutions per hit (fuzznuc -pmismatch, default 0).")] int maxMismatches = 0,
        [Description("Also search the reverse-complement strand (fuzznuc -complement, default false).")] bool bothStrands = false,
        [Description("Read IUPAC ambiguity codes in the SEQUENCE as the fuzznuc binary does (default false: they match nothing).")] bool sequenceAmbiguity = false,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "fuzznuc searches the range of the (reverse-complemented) sequence and reports the ORIGINAL strand: positions are original, the strand is flipped for -sreverse; the sequence is read as fuzznuc reads it (ambiguity codes on)." + " " + EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null,
        [Description("EMBOSS report format (-rformat): any of simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq (diffseq needs 7 report tags: refused as the program does; embl / genbank are nucleotide only, pir / swiss protein only), or the report qualifiers as typed after the files, e.g. '-rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2'. The report FILE of the program (header with the command line, one block per sequence, footer; #Rundate / #!Date is now; sequence name from sequenceName) is returned in the result field report; needs one pattern (pattern, not patternFile). Omit for no report text.")] string? reportFormat = null,
        [Description("Sequence name printed in the report (default sequence).")] string sequenceName = "sequence")
    {
        if (reportFormat is not null)
        {
            if (string.IsNullOrWhiteSpace(pattern) || !string.IsNullOrWhiteSpace(patternFile))
                throw new ArgumentException("reportFormat needs one pattern (the report names the pattern 'pattern')", nameof(reportFormat));
            var plain = FindPatternHits(sequence, pattern, null, maxMismatches, bothStrands, sequenceAmbiguity, cancellationToken, sequenceQualifiers, null, sequenceName);
            string[] cmd = [EmbossCommandLineEcho.Qualifier("pattern", pattern), Qualifier("pmismatch", maxMismatches), EmbossCommandLineEcho.Qualifier("complement", bothStrands ? "Y" : "N", true)];
            var reportHits = plain.Hits.Select(h => new PatternHit(h.Position, h.Length, h.Strand == "reverse", h.Mismatches)).ToList();
            return plain with { Report = EmbossReportText("fuzznuc", sequenceName, sequence, sequenceQualifiers, reportFormat, cmd, r => EmbossMotifReportPrograms.Fuzznuc(r, pattern, maxMismatches, bothStrands, reportHits), out _, cancellationToken) };
        }

        if (sequenceQualifiers is not null)
        {
            var q = EmbossSequenceQualifiers.ParseInputs("fuzznuc", sequenceQualifiers)[0];
            if (maxMismatches < 0)
                throw new ArgumentOutOfRangeException(nameof(maxMismatches), "maxMismatches cannot be negative");
            if (string.IsNullOrWhiteSpace(pattern) == string.IsNullOrWhiteSpace(patternFile))
                throw new ArgumentException("Give exactly one of pattern and patternFile", nameof(pattern));
            var named = string.IsNullOrWhiteSpace(patternFile)
                ? new[] { new NamedPattern("pattern1", FuzznucPattern.Parse(pattern!)) }
                : FuzznucPattern.ParsePatternFile(patternFile, null, maxMismatches).ToArray();
            var qualifiedItems = named.SelectMany(p => EmbossMotifPrograms.FindPatternHits(sequence, q, p.Pattern, maxMismatches, bothStrands, cancellationToken)
                .Select(h => new PatternHitItem(p.Name, h.Position, h.Length, h.IsReverseStrand ? "reverse" : "forward", h.Mismatches))).ToArray();
            return new FindPatternHitsResult(named.ToDictionary(p => p.Name, p => p.Pattern.ToString()), qualifiedItems);
        }

        var nucleotides = RequireEmbossIupacDna(sequence, nameof(sequence));
        if (maxMismatches < 0)
            throw new ArgumentOutOfRangeException(nameof(maxMismatches), "maxMismatches cannot be negative");
        if (string.IsNullOrWhiteSpace(pattern) == string.IsNullOrWhiteSpace(patternFile))
            throw new ArgumentException("Give exactly one of pattern and patternFile", nameof(pattern));

        var patterns = string.IsNullOrWhiteSpace(patternFile)
            ? new[] { new NamedPattern("pattern1", FuzznucPattern.Parse(pattern!)) }
            : FuzznucPattern.ParsePatternFile(patternFile, null, maxMismatches).ToArray();
        var hits = EmbossMotifPrograms.FindPatternHits(nucleotides, patterns, maxMismatches, bothStrands, cancellationToken, sequenceAmbiguity);
        var items = hits
            .SelectMany(kv => kv.Value.Select(h => new PatternHitItem(kv.Key, h.Position, h.Length, h.IsReverseStrand ? "reverse" : "forward", h.Mismatches)))
            .ToArray();
        return new FindPatternHitsResult(patterns.ToDictionary(p => p.Name, p => p.Pattern.ToString()), items);
    }
    [McpServerTool(Name = "generate_emboss_consensus", Title = "Motifs — EMBOSS cons Consensus", ReadOnly = true)]
    [Description("Scoring-matrix plurality consensus of an alignment, identical to EMBOSS 6.6.0 'cons' (EDNAFULL for nucleotides, EBLOSUM62 for proteins; N/X where no residue reaches the plurality; lower case at or below setcase). Gaps '-', '.', '~' allowed. residueType 'auto' decides the type from the first sequence like cons; padRaggedRows pads shorter rows with trailing gaps like cons (ajSeqsetFill). Rows are read as EMBOSS reads an alignment (A4-3): letters and '*' are kept, '.' '~' '-' are gaps, '?' '#' '+' are the unknown residue, every other character (digits, blanks, '_', punctuation) is dropped, and rows that become unequal are padded with trailing gaps.")]
    public static ConsensusResult GenerateEmbossConsensus(
        [Description("At least two aligned sequences (DNA or protein); equal length unless padRaggedRows.")] string[] sequences,
        [Description("Residue type: 'nucleotide' (default, EDNAFULL, no-consensus N), 'protein' (EBLOSUM62, no-consensus X) or 'auto' (decided from the first sequence as the cons program does).")] string residueType = "nucleotide",
        [Description("Minimum positive-match weight for a consensus residue (cons -plurality); default half the total sequence weight.")] double? plurality = null,
        [Description("Required number of identical residues at a position (cons -identity, default 0 = off).")] int identity = 0,
        [Description("Positive-match weight at or below which the residue is written in lower case (cons -setcase); default half the total sequence weight.")] double? setcase = null,
        [Description("Optional per-sequence weights (one per sequence, finite, >= 0; default 1.0 each).")] double[]? weights = null,
        [Description("Pad rows shorter than the longest with trailing '-' gaps (cons ajSeqsetFill) instead of rejecting unequal lengths (default false).")] bool padRaggedRows = false)
    {
        if (sequences is null || sequences.Length < 2)
            throw new ArgumentException("At least two sequences are required", nameof(sequences));
        var type = (residueType ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "nucleotide" => global::Seqeron.Genomics.Analysis.ConsensusResidueType.Nucleotide,
            "protein" => global::Seqeron.Genomics.Analysis.ConsensusResidueType.Protein,
            "auto" => global::Seqeron.Genomics.Analysis.ConsensusResidueType.Auto,
            _ => throw new ArgumentException("residueType must be 'nucleotide', 'protein' or 'auto'", nameof(residueType)),
        };

        var consensus = EmbossConsensus.GenerateEmbossConsensus(
            sequences, padRaggedRows, type, (float?)plurality, identity, (float?)setcase,
            weights?.Select(w => (float)w).ToArray());
        return new ConsensusResult(consensus);
    }
    [McpServerTool(Name = "consambig", Title = "Motifs — EMBOSS consambig Ambiguity Consensus", ReadOnly = true)]
    [Description("Ambiguity-code consensus of an alignment, identical to EMBOSS 6.6.0 'consambig' (1 400 random alignments equal): per column the residue codes are OR-ed (nucleotides A, C, G, T as bits, so R = A|G, V = A|C|G, N = any; amino acids one bit each, B = D|N, Z = E|Q, J = I|L, X = all twenty) and written as the one letter that stands for exactly that set; there is no vote (unlike cons: one deviating base widens the code). A column that contains a gap ('-', '.', '~', or the padding of a shorter row) is written in lower case. An all-gap nucleotide column is 'n'; an amino-acid union that is no letter (A|D) is 'X'; a protein column without any scored residue (all gaps, or only '*') is the NUL character, as in the program. The set type is that of the first sequence (A/C/G/T/U + IUPAC letters = nucleotide) unless residueType 'nucleotide' forces it (residueType 'protein' cannot make a nucleotide-looking set protein, as in the program). Rows are upper-cased, '.'/'~' read as '-', '?' as X; shorter rows are padded with trailing gaps (the program always does). Rows are read as EMBOSS reads an alignment (A4-3): letters and '*' are kept, '.' '~' '-' are gaps, '?' '#' '+' are the unknown residue, every other character (digits, blanks, '_', punctuation) is dropped, and rows that become unequal are padded with trailing gaps.")]
    public static ConsensusResult Consambig(
        [Description("At least two aligned sequences (DNA or protein).")] string[] sequences,
        [Description("'auto' (default: type of the first sequence, as the program), 'nucleotide' (forces nucleotide, -snucleotide) or 'protein'.")] string residueType = "auto",
        CancellationToken cancellationToken = default)
    {
        if (sequences is null || sequences.Length < 2)
            throw new ArgumentException("At least two sequences are required", nameof(sequences));
        return new ConsensusResult(EmbossMotifPrograms.GenerateEmbossAmbiguousConsensus(sequences, ParseConsensusResidueType(residueType), cancellationToken: cancellationToken));
    }

    [McpServerTool(Name = "find_motif", Title = "Genome — Find Motif (Exact)", ReadOnly = true)]
    [Description("All exact occurrences of a motif (case-insensitive, overlapping, ascending 0-based positions, forward strand) in a DNA sequence. IUPAC input (e.g. N runs) is accepted: a symbol other than A/C/G/T never matches, so no window containing one is an occurrence.")]
    public static MotifPositionsResult FindMotif(
        [Description("Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive); non-A/C/G/T symbols match nothing.")] string sequence,
        [Description("Motif to locate.")] string motif)
    {
        var nucleotides = RequireIupacDna(sequence, nameof(sequence));
        if (string.IsNullOrEmpty(motif))
            throw new ArgumentException("Motif cannot be null or empty", nameof(motif));

        var positions = EmbossMotifPrograms.FindMotif(nucleotides, motif);
        return new MotifPositionsResult(positions.ToArray());
    }
    [McpServerTool(Name = "find_known_motifs", Title = "Genome — Search Known Motif Set", ReadOnly = true)]
    [Description("Search for a user-provided set of motifs simultaneously (exact by default; optional IUPAC codes, reverse strand and a mismatch budget). Only motifs with at least one hit are returned; positions ascending, 0-based, forward-strand window starts. mismatches = k > 0 allows up to k substituted positions per window (Hamming distance, EMBOSS fuzznuc -pmismatch; no gaps) and adds matchMismatches / reverseMismatches (the distance of each reported hit, parallel to the position arrays).")]
    public static FindKnownMotifsResult FindKnownMotifs(
        [Description("Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive); a symbol other than A/C/G/T matches no motif symbol (not even N) and costs a mismatch when mismatches > 0.")] string sequence,
        [Description("Set of motifs to search.")] string[] motifs,
        [Description("Also search the reverse strand: reverseMatches lists the windows matching the reverse complement of each motif (default false).")] bool bothStrands = false,
        [Description("Interpret IUPAC ambiguity codes in the motifs, e.g. TATAWAW (default false: codes are literal; without mismatches they match nothing, with mismatches each costs one).")] bool iupac = false,
        [Description("Maximum number of mismatching positions per hit (default 0 = exact; fuzznuc -pmismatch). A bound >= the motif length makes every window a hit.")] int mismatches = 0,
        CancellationToken cancellationToken = default)
    {
        var nucleotides = RequireIupacDna(sequence, nameof(sequence));
        if (motifs is null)
            throw new ArgumentException("Motifs cannot be null", nameof(motifs));
        if (mismatches < 0)
            throw new ArgumentOutOfRangeException(nameof(mismatches), "Mismatches cannot be negative");

        var hits = EmbossMotifPrograms.FindKnownMotifHits(nucleotides, motifs, bothStrands, iupac, mismatches, cancellationToken);
        var matches = new Dictionary<string, int[]>(hits.Count);
        var reverse = new Dictionary<string, int[]>();
        var matchMismatches = mismatches > 0 ? new Dictionary<string, int[]>(hits.Count) : null;
        var reverseMismatches = mismatches > 0 ? new Dictionary<string, int[]>() : null;
        foreach (var kvp in hits)
        {
            var forwardHits = kvp.Value.Where(h => !h.IsReverseStrand).ToArray();
            var reverseHits = kvp.Value.Where(h => h.IsReverseStrand).ToArray();
            if (forwardHits.Length > 0)
            {
                matches[kvp.Key] = forwardHits.Select(h => h.Position).ToArray();
                if (matchMismatches is not null)
                    matchMismatches[kvp.Key] = forwardHits.Select(h => h.Mismatches).ToArray();
            }

            if (reverseHits.Length > 0)
            {
                reverse[kvp.Key] = reverseHits.Select(h => h.Position).ToArray();
                if (reverseMismatches is not null)
                    reverseMismatches[kvp.Key] = reverseHits.Select(h => h.Mismatches).ToArray();
            }
        }

        return new FindKnownMotifsResult(matches, reverse, matchMismatches, reverseMismatches);
    }
}
