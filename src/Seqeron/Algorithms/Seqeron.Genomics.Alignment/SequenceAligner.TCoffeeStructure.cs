namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// The structure based pair-wise methods of the T-Coffee 13.41 program that call an EXTERNAL structure aligner (audit A2-10):
//   TMalign_pair   = pdb_pair (util_dp_drivers.c) + tc_generic_method.pl -mode=pdb_pair -method=TMalign   (Zhang & Skolnick 2005, `TMalign pdb1 pdb2 > tmp1`)
//   mustang_pair   = pdb_pair + tc_generic_method.pl -mode=pdb_pair -method=mustang                         (Konagurthu et al. 2006, `mustang -i pdb1 pdb2 -F fasta`)
// The structure aligner itself is the external program: this port is T-Coffee's side of the interface, fed with what the aligner wrote.
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
    }
}
