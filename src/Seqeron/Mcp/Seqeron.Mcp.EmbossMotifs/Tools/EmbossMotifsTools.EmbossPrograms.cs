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
    [McpServerTool(Name = "fuzzpro", Title = "Protein — EMBOSS fuzzpro pattern search", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 fuzzpro: every window of a protein sequence that matches a PROSITE-like pattern: amino-acid letters (case-insensitive; B = D or N, Z = E or Q, J = I or L; x = any residue), classes [ILV], exclusions {P} (plain letters only: B, Z, J, X inside [] or {} are rejected as in the binary), repeats x(3), length ranges x(2,4) (every matching window length is reported; (0,n) makes the element optional), '-' separators, the anchors '<' (N terminus) and '>' (C terminus), a terminating '.'; e.g. [DE](2)HS{P}x(2)Px(2,4)C. The search is exact unless mismatch > 0 (fuzzpro -pmismatch; verified on 750 random patterns with 1 - 3 mismatches). A start-anchored pattern with a range reports only its longest match (the binary's regular-expression search class); a single trailing stop '*' is not part of the sequence. Positions are 0-based; hits are ordered by position, then length (with mismatch > 0 and a range: the program's order, a window listed once per matching unrolling). Verified against fuzzpro 6.6.0 on 900 random patterns (letters, classes, exclusions, repeats, ranges, anchors; sequences with X, B, Z, U, J, O, * and lower case). For nucleotide patterns see find_pattern_hits (fuzznuc), for a translated search fuzztran. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static FuzzproResult Fuzzpro(
        [Description("Protein sequence (letters; a trailing * is dropped).")] string sequence,
        [Description("The pattern (fuzzpro -pattern).")] string pattern,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "fuzzpro searches the range and counts positions from the start of the original sequence; a protein is never reversed." + " " + EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null,
        [Description("EMBOSS report format (-rformat): any of simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq (diffseq needs 7 report tags: refused as the program does; embl / genbank are nucleotide only, pir / swiss protein only), or the report qualifiers as typed after the files, e.g. '-rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2'. The report FILE of the program (header with the command line, one block per sequence, footer; #Rundate / #!Date is now; sequence name from sequenceName) is returned in the result field report. Omit for no report text.")] string? reportFormat = null,
        [Description("Sequence name printed in the report (default sequence).")] string sequenceName = "sequence",
        [Description("Maximum substituted residues per hit (fuzzpro -pmismatch, default 0 = exact; must be less than the real pattern length for a pattern with a class or an exclusion and no range: the program stops otherwise). Embpat's search engines: Baeza-Yates-Perleberg, Tarhio-Ukkonen-Bleasby (its skip table knows the letters A-Z only) or the recursive brute force, which with a range lets a mismatch fall only on the mandatory part and lists a window once per matching unrolling; each hit then carries its mismatches (the report: Score = length - mismatches and a *mismatch note).")] int mismatch = 0)
    {
        if (reportFormat is not null)
        {
            var plainResult = Fuzzpro(sequence, pattern, cancellationToken, sequenceQualifiers, null, sequenceName, mismatch);          // one search; the report reuses its hits
            var reportHits = plainResult.Hits.Select(h => new ProteinPatternHit(h.Position, h.Length, h.Mismatches)).ToList();
            return plainResult with
            {
                Report = EmbossReportText("fuzzpro", sequenceName, sequence, sequenceQualifiers, reportFormat, [EmbossCommandLineEcho.Qualifier("pattern", pattern), .. (mismatch > 0 ? new[] { Qualifier("pmismatch", mismatch) } : [])], r => EmbossMotifReportPrograms.Fuzzpro(r, pattern, reportHits, mismatch), out _, cancellationToken),
            };
        }

        if (sequenceQualifiers is not null)
        {
            RequireFuzzproInput(sequence, pattern);
            return new FuzzproResult(ProteinPattern.Parse(pattern).ToString(), [.. EmbossMotifPrograms.FindProteinPatternHits(sequence, EmbossSequenceQualifiers.ParseInputs("fuzzpro", sequenceQualifiers)[0], pattern, cancellationToken, mismatch).Select(h => new ProteinPatternHitItem(h.Position, h.Length, h.Mismatches))]);
        }

        var (parsed, hits) = RunFuzzpro(sequence, pattern, cancellationToken, mismatch);
        return new FuzzproResult(parsed.ToString(), [.. hits.Select(h => new ProteinPatternHitItem(h.Position, h.Length, h.Mismatches))]);
    }
    [McpServerTool(Name = "fuzztran", Title = "Genome — EMBOSS fuzztran (protein pattern search after translation)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 fuzztran: translates a nucleotide sequence in the chosen frames (frame 1, 2, 3, -1, -2, -3, F = forward three, R = reverse three, 6 = all; transeq conventions: reverse frame -k is the reverse complement read from its k-th base, dangling bases give X) with an EMBOSS genetic-code table (0 = standard, else the NCBI table number) and searches every translation with a fuzzpro protein pattern (see fuzzpro: letters, x, [ ], { }, repeats, ranges, anchors; exact). Each hit gives the frame, the 1-based protein start / end, the matched residues, the forward-strand nucleotide location (position 0-based, length = 3 per residue) and embossStart / embossEnd, the two numbers the program prints (exact for the forward frames; EMBOSS 6.6.0 prints negative numbers for the reverse frames - a defect of its source-feature mapping, reproduced; use position / length). Hits per frame in the order 1, 2, 3, -1, -2, -3, by protein position, then length. Verified against fuzztran 6.6.0 on 760 random DNA / pattern / frame / table cases (all frame selections, tables 0-11, N / IUPAC codes and lower case). Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static FuzztranResult Fuzztran(
        [Description("Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive).")] string sequence,
        [Description("The protein pattern (fuzztran -pattern).")] string pattern,
        [Description("Frames to translate (fuzztran -frame): 1 (default), 2, 3, -1, -2, -3, F, R or 6.")] string frame = "1",
        [Description("EMBOSS genetic code (fuzztran -table): 0 = standard (default), otherwise an NCBI table number.")] int table = 0,
        CancellationToken cancellationToken = default,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "fuzztran cuts the (reverse-complemented) text at the begin / end AS GIVEN, numbers the hits from the first position of that text and prints negative numbers for the reverse frames (a defect of the program, reproduced)." + " " + EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null,
        [Description("EMBOSS report format (-rformat): any of simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq (diffseq needs 7 report tags: refused as the program does; embl / genbank are nucleotide only, pir / swiss protein only), or the report qualifiers as typed after the files, e.g. '-rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2'. The report FILE of the program (header with the command line, one block per sequence, footer; #Rundate / #!Date is now; sequence name from sequenceName) is returned in the result field report. Omit for no report text.")] string? reportFormat = null,
        [Description("Sequence name printed in the report (default sequence).")] string sequenceName = "sequence",
        [Description("Maximum substituted residues per hit (fuzztran -pmismatch, default 0 = exact; see fuzzpro).")] int mismatch = 0)
    {
        if (reportFormat is not null)
        {
            var plainResult = Fuzztran(sequence, pattern, frame, table, cancellationToken, sequenceQualifiers, null, sequenceName, mismatch);          // one search; the report reuses its hits
            var reportHits = plainResult.Hits.Select(h => new TranslatedPatternHit(h.Frame, h.ProteinStart, h.ProteinEnd, h.Translation, h.Position, h.Length, h.EmbossStart, h.EmbossEnd, h.Mismatches)).ToList();
            return plainResult with
            {
                Report = EmbossReportText("fuzztran", sequenceName, sequence, sequenceQualifiers, reportFormat, [EmbossCommandLineEcho.Qualifier("pattern", pattern), Qualifier("frame", frame), Qualifier("table", table), .. (mismatch > 0 ? new[] { Qualifier("pmismatch", mismatch) } : [])],
                    r => EmbossMotifReportPrograms.Fuzztran(r, pattern, frame, table, reportHits, mismatch), out _, cancellationToken),
            };
        }

        if (sequenceQualifiers is not null)
        {
            if (string.IsNullOrWhiteSpace(pattern))
                throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
            var qualifiedParsed = ProteinPattern.Parse(pattern);
            var qualifiedHits = EmbossMotifPrograms.FindTranslatedPatternHits(sequence, EmbossSequenceQualifiers.ParseInputs("fuzztran", sequenceQualifiers)[0], pattern, frame, table, cancellationToken, mismatch)
                .Select(TranslatedItem).ToArray();
            return new FuzztranResult(qualifiedParsed.ToString(), qualifiedHits);
        }

        var (parsed, hits) = RunFuzztran(sequence, pattern, frame, table, cancellationToken, mismatch);
        return new FuzztranResult(parsed.ToString(), [.. hits.Select(TranslatedItem)]);
    }
}
