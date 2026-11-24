namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// The structure based pair-wise methods of the T-Coffee 13.41 program that call an EXTERNAL structure aligner (audit A2-10):
//   TMalign_pair   = pdb_pair (util_dp_drivers.c) + tc_generic_method.pl -mode=pdb_pair -method=TMalign   (Zhang & Skolnick 2005, `TMalign pdb1 pdb2 > tmp1`)
//   mustang_pair   = pdb_pair + tc_generic_method.pl -mode=pdb_pair -method=mustang                         (Konagurthu et al. 2006, `mustang -i pdb1 pdb2 -F fasta`)
// The structure aligner itself is the external program: this port is T-Coffee's side of the interface, fed with what the aligner wrote.
//   sap_pair (audit A2-12): with SAP NOT installed the method list switches to TMalign (see TcPairJobEntries); with SAP 1.1.3 installed the glue is `sap_pair` of util_dp_drivers.c (L2049):
//       SAP runs on the normalised PDB files (`sap file1 file2 > <template1>_<template2>.sap_results` in the cache), find_token_in_file finds the SECOND token that starts with "Percent"
//       (a token that starts with '*' skips the rest of its line), the lines after that one are read with sscanf " %c %*d %*f %*d %c" after every '*' was removed (lines with "eighted" / "RMSd" skipped,
//       both characters must be letters; the first character is the residue of the SECOND file = sequence s2, the last one that of the first file = s1), so the aligned residues of s2 and s1 are two strings of
//       the same length; they are written as a TC_LIB_FORMAT_01 file (sequence s2 = the first string, sequence s1 = the second string, "#1 2", "a a 100 1 0" for every aligned pair: the weight is the number
//       "100" of the method file) which read_constraint_list_file maps to the real sequences with fix_seq_seq (identical = the same coordinates, else a myers_miller alignment of the compact string to the sequence).
//       A pair that SAP aligned to no residue pair ("score = 0", -nan) adds no entry. A pair for which SAP left no results file at all stops the program (status 1): see TcSapPairEntries.
//       fix_seq_seq types the library-file sequences (get_sequence_type): PROTEIN -> the matrix blosum62mt, else idmat (the sets of aln2seq have no type: idmat).
//
//   fill_P_template (aln_convertion_util.c): a `_P_` template (a PDB file) is attached to a sequence only when the sequence and the PDB sequence are alike enough:
//       align_two_sequences (seq, pdb sequence, "idmat", -3, 0, "fasta_pair_wise")  (lower case, k-tuple 3, terminal gaps free), aln2sim (idmat) > -pdb_min_sim (35)
//       and aln2coverage (of the sequence) > -pdb_min_cov (50). A sequence without a template makes `check_seq_type` fail for the method (SEQ_TYPE P) and the pair is
//       aligned by proba_pair instead ("Method TMalign cannot be applied to [a vs b], proba_pair will be used instead").
//   TMalign: the glue script keeps `tail -4` of the standard output, reverses the lines (l[3] = the first of the four, l[1] = the third) and replaces every character outside
//       [a-zA-Z0-9-] by '-' - the line break included, so each row gets one extra '-': a column of two gaps - and writes ">\nl3\n>\nl1\n".
//   mustang: the file results.afasta (two FASTA records, the first record is the first structure).
//   Either way the alignment goes through aln2constraint_list (weight "sim": 10 x identity of the alignment, see TcPairJobEntries); a structure aligner that wrote nothing adds no entry
//   ("pdb_pair/TMalign failed").
// =====================================================================================================================

public static partial class SequenceAligner
{
    /// <summary>
    /// The alignment the T-Coffee glue script <c>tc_generic_method.pl</c> takes from the standard output of <c>TMalign pdb1.pdb pdb2.pdb</c>: the first and third of the last four lines
    /// (the aligned sequence of the first structure, the aligned sequence of the second one), every character outside <c>[a-zA-Z0-9-]</c> replaced by '-' and the line break turned into a final '-' as well.
    /// </summary>
    /// <param name="tmalignOutput">The standard output of TMalign (only its last four lines matter).</param>
    /// <returns>The two rows, or <c>null</c> when the output has fewer than four lines.</returns>
    public static (string RowA, string RowB)? TCoffeeTMalignAlignment(string tmalignOutput)
    {
        ArgumentNullException.ThrowIfNull(tmalignOutput);
        var lines = tmalignOutput.Split('\n').ToList();
        if (lines.Count > 0 && lines[^1].Length == 0) lines.RemoveAt(lines.Count - 1); // `tail -4` counts the lines ended by a line break
        if (lines.Count < 4) return null;
        var tail = lines.Skip(lines.Count - 4).ToArray();  // [seq1, middle, seq2, empty] when the program ends with an empty line
        static string Clean(string s) => new string(s.Select(c => char.IsAsciiLetterOrDigit(c) || c == '-' ? c : '-').ToArray()) + "-";
        return (Clean(tail[0]), Clean(tail[2]));
    }

    /// <summary>
    /// The text t_coffee prints (standard error) and the status it exits with when the command <c>DaliLite</c> does not exist (DaliLite v5 has only <c>dali.pl</c> and its helpers): the glue <c>tc_generic_method.pl -mode=pdb_pair -method=DaliLite</c> runs
    /// <c>DaliLite -pairwise pdb1.pdb pdb2.pdb &gt;tmp1</c> through its <c>safe_system</c>, which cannot start it and exits with status 1; t_coffee reports "ERROR: UNSPECIFIED UNSPECIFIED" with the glue command line and ends with
    /// "TERMINATION STATUS: FAILURE" (exit status 1). Oracle: 4 x 60 sets (<c>handoff/A2-12b-dali-sap-stderr</c>, runs 3, 4, 5 and the committed run).
    /// </summary>
    internal const string DaliLiteMissingMessage = "ERROR: UNSPECIFIED UNSPECIFIED (COM: <glue>.pl -mode=pdb_pair -method=DaliLite -pdbfile1=<f> -pdbfile2=<f> -outfile=<f> -tmpdir=<d>): no DaliLite command, the glue's `DaliLite -pairwise pdb1.pdb pdb2.pdb` cannot run; exit status 1";

    /// <summary>
    /// The alignment the glue <c>tc_generic_method.pl</c> (<c>seq2pdb_pair</c>, <c>-method=DaliLite</c>) takes from the standard output of <c>DaliLite -pairwise pdb1.pdb pdb2.pdb</c> (DaliLite v3): the first token after every <c>Query</c> is
    /// appended to the first row, the first token after every <c>Sbjct</c> to the second one (Perl <c>/Query\s+(\S+)/</c>, anywhere in the line), every '.' becomes '-'. No <c>result.aln</c> is written when either row stays empty.
    /// </summary>
    /// <param name="daliLiteOutput">The standard output of <c>DaliLite -pairwise</c>.</param>
    /// <returns>The two rows, or <c>null</c> when the output has no <c>Query</c> or no <c>Sbjct</c> token (the pair adds no library entry).</returns>
    public static (string RowA, string RowB)? TCoffeeDaliLiteAlignment(string daliLiteOutput)
    {
        ArgumentNullException.ThrowIfNull(daliLiteOutput);
        var a = new System.Text.StringBuilder();
        var b = new System.Text.StringBuilder();
        foreach (string line in daliLiteOutput.Split('\n'))
        {
            var q = DaliQueryToken.Match(line);
            if (q.Success) a.Append(q.Groups[1].Value);
            else if (DaliSbjctToken.Match(line) is { Success: true } sb) b.Append(sb.Groups[1].Value);
        }

        return a.Length == 0 || b.Length == 0 ? null : (a.ToString().Replace('.', '-'), b.ToString().Replace('.', '-'));
    }

    private static readonly System.Text.RegularExpressions.Regex DaliQueryToken = new(@"Query\s+(\S+)", System.Text.RegularExpressions.RegexOptions.CultureInvariant);
    private static readonly System.Text.RegularExpressions.Regex DaliSbjctToken = new(@"Sbjct\s+(\S+)", System.Text.RegularExpressions.RegexOptions.CultureInvariant);

    /// <summary>
    /// The aligned residue pairs the glue <c>sap_pair</c> of T-Coffee 13.41 (<c>util_dp_drivers.c</c> L2049) reads from the output of SAP 1.1.3 (<c>sap file1 file2</c>, the file <c>&lt;template1&gt;_&lt;template2&gt;.sap_results</c> of the cache):
    /// the characters after the second token that starts with "Percent" (<c>is_sap_file</c> needs two of them) are scanned line by line - every '*' removed, lines with "eighted" or "RMSd" skipped, the rest read as
    /// <c>" %c %*d %*f %*d %c"</c> - and a pair is kept when both characters are letters. <c>First</c> is the residue of the SECOND structure file (the program names it sequence s2), <c>Second</c> that of the first one (s1).
    /// </summary>
    /// <param name="sapResults">The content of the results file.</param>
    /// <returns>The two strings of aligned residues of equal length (empty strings when SAP aligned nothing), or <c>null</c> when the file is not a SAP result (fewer than two "Percent" tokens): "SAP failed to align".</returns>
    public static (string First, string Second)? TCoffeeSapAlignedResidues(string sapResults)
    {
        ArgumentNullException.ThrowIfNull(sapResults);
        // find_token_in_file: whitespace separated tokens; a token that starts with '*' skips the rest of its line; a token whose first characters are the token is found (strncmp)
        int pos = 0, found = 0;
        while (found < 2)
        {
            while (pos < sapResults.Length && char.IsWhiteSpace(sapResults[pos])) pos++;
            if (pos >= sapResults.Length) return null;
            int start = pos;
            while (pos < sapResults.Length && !char.IsWhiteSpace(sapResults[pos])) pos++;
            if (sapResults[start] == '*')
            {
                while (pos < sapResults.Length && sapResults[pos] != '\n') pos++;
            }
            else if (string.CompareOrdinal(sapResults, start, "Percent", 0, 7) == 0) found++;
        }

        while (pos < sapResults.Length && sapResults[pos] != '\n') pos++;
        var first = new System.Text.StringBuilder();
        var second = new System.Text.StringBuilder();
        foreach (string raw in sapResults.Substring(pos).Split('\n'))
        {
            string line = raw.Replace("*", "");
            if (line.Contains("eighted", StringComparison.Ordinal) || line.Contains("RMSd", StringComparison.Ordinal)) continue;
            var m = SapResidueLine.Match(line);
            if (!m.Success) continue;
            char r1 = m.Groups[1].Value[0], r2 = m.Groups[2].Value[0];
            if (!char.IsAsciiLetter(r1) || !char.IsAsciiLetter(r2)) continue;
            first.Append(r1);
            second.Append(r2);
        }

        return (first.ToString(), second.ToString());
    }

    private static readonly System.Text.RegularExpressions.Regex SapResidueLine = new(
        @"^\s*(\S)\s+[+-]?\d+\s+(?:[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|[+-]?(?:nan|inf(?:inity)?))\s+[+-]?\d+\s+(\S)",
        System.Text.RegularExpressions.RegexOptions.IgnoreCase | System.Text.RegularExpressions.RegexOptions.CultureInvariant);

    /// <summary>The entries one pair job (x, y) adds with the REAL SAP installed (<c>sap_pair</c> of util_dp_drivers.c): see the file header. The method-file weight is the number 100 (<c>WEIGHT 100</c>).</summary>
    private static List<(int S1, int R1, int S2, int R2, int W)> TcSapPairEntries(string? sapResults, IReadOnlyList<string> seqs, int x, int y, TcProcessGroups groups)
    {
        var entries = new List<(int, int, int, int, int)>();
        // "SAP failed to align" (no results file, or one without the two "Percent" tokens): the job returns without an entry (util_dp_drivers.c L2144-2150), and in the oracle EVERY run in which a pair of two usable
        // structures has no results file ended with exit status 1 and no alignment ("You have not provided any method or enough Sequences[FATAL]", the only line t_coffee printed with -no_warning): 76 of 76 stopped sets of 520
        // (runs 1, 3, 4, 5 of handoff/A2-12b; the 444 other sets ran through): the program stops where SAP left nothing
        if (sapResults is null || TCoffeeSapAlignedResidues(sapResults) is not var (resS2, resS1))
            throw new TCoffeeAbortedException($"{TCoffeeAbortedException.NoMethodMessage} (SAP 1.1.3 left no result for the pair of structures {x} / {y}: \"SAP failed to align\" is a warning that -no_warning hides, the library of the method ends empty and t_coffee.c L5005 stops with exit status 1)");
        if (resS1.Length == 0) return entries;                                                                // tot = 0: no library file
        // the library file names sequence s2 with the first string and s1 with the second; fix_seq_seq maps both compact strings to the sequences
        // read_seq_in_n_list types the two strings of the library file (get_sequence_type): PROTEIN selects blosum62mt in fix_seq_seq, anything else idmat
        var matrix = TcSeqType(new[] { resS2, resS1 }) == TcSeqKind.Protein ? TcMatrix.ByName("blosum62mt")! : TcMatrix.IdMat;
        int[] map2 = TcFixIndex(resS2, seqs[y], out bool aligned2, matrix);
        int[] map1 = TcFixIndex(resS1, seqs[x], out bool aligned1, matrix);
        groups.SingleLetters ??= aligned2 || aligned1; // aln2sim(B, "idmat") of an unequal pair builds the process-static residue groups
        for (int a = 1; a <= resS1.Length; a++)
            if (map2[a] > 0 && map1[a] > 0)
                entries.Add((y, map2[a], x, map1[a], 100));
        return entries;
    }

    /// <summary>The alignment of the two records of <c>results.afasta</c> written by <c>mustang -i pdb1.pdb pdb2.pdb -F fasta</c> (the first record belongs to the first structure).</summary>
    /// <param name="resultsAfasta">The content of results.afasta.</param>
    /// <returns>The two rows, or <c>null</c> when the file does not hold two records.</returns>
    public static (string RowA, string RowB)? TCoffeeMustangAlignment(string resultsAfasta)
    {
        ArgumentNullException.ThrowIfNull(resultsAfasta);
        var rows = new List<System.Text.StringBuilder>();
        foreach (string raw in resultsAfasta.Split('\n'))
        {
            string line = raw.Trim();
            if (line.Length == 0) continue;
            if (line[0] == '>') rows.Add(new System.Text.StringBuilder());
            else if (rows.Count > 0) rows[^1].Append(line);
        }

        return rows.Count < 2 ? null : (rows[0].ToString(), rows[1].ToString());
    }

    /// <summary>
    /// <c>fill_P_template</c>: whether the PDB sequence is accepted as the structure of the sequence (<c>-pdb_min_sim</c> 35, <c>-pdb_min_cov</c> 50): the identity of their <c>fasta_pair_wise</c> alignment
    /// (lower case, identity matrix, open -3, extend 0, k-tuple 3) must exceed the first and the share of the sequence's residues aligned to a PDB residue the second.
    /// </summary>
    internal static bool TcPdbTemplateUsable(string sequence, string pdbSequence, int minSimilarity = 35, int minCoverage = 50)
    {
        var (a, b) = TcFastaPairRows(sequence.ToLowerInvariant(), pdbSequence.ToLowerInvariant(), TcMatrix.ByName("idmat")!, -3, 0, 3, 0, 0, 0f, 2);
        int sim = TcSeqSimMeasure(a, b, "idmat");
        int npos = 0, covered = 0;
        for (int i = 0; i < a.Length; i++)
        {
            if (a[i] == '-') continue;
            npos++;
            if (b[i] != '-') covered++;
        }

        int cov = npos == 0 ? 0 : 100 * covered / npos;
        return sim > minSimilarity && cov > minCoverage;
    }

    /// <summary>What the pair jobs of the structure methods need: which sequences carry a usable template and what the structure aligner wrote for a pair (indices in the program's order, x &lt; y).</summary>
    internal sealed class TcStructure
    {
        public bool[] Usable = Array.Empty<bool>();
        public Func<TCoffeePairMethod, int, int, string?> Output = (_, _, _) => null;

        /// <summary>SAP is installed (the caller supplied the results of <see cref="TCoffeePairMethod.SapPair"/>, possibly none): <c>sap_pair</c> is the real glue, not the TM-align replacement.</summary>
        public bool SapInstalled;

        /// <summary>A <c>DaliLite</c> command exists (the caller supplied the outputs of <see cref="TCoffeePairMethod.DaliPair"/>, possibly none).</summary>
        public bool DaliInstalled;
    }
}
