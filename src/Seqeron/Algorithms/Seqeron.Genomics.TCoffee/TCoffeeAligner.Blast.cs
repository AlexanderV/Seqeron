// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.TCoffee;

using Seqeron.Genomics.ClustalW;
using static Seqeron.Genomics.ClustalW.ClustalWAligner;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

// =====================================================================================================================
// The BLAST based multiple-alignment methods of the T-Coffee 13.41 program (audit A2-9). `blastp_msa` (ALN_MODE multiple, EXECUTABLE seq_msa, EXECUTABLE2 blastp): the program formats the sequences as a local
// BLAST protein database (`makeblastdb -dbtype prot`), searches every sequence against it (`blastall -p blastp -m7`, NCBI BLAST+ XML via legacy_blast.pl) and turns the high scoring pairs into a library
// (tc_generic_method.pl: blast_msa, blast_xml2profile / ncbi_blast_xml2profile, xml2tag_list):
//     for every HSP of query i against hit j          "# i j" (1-based database order)
//        for every column with a residue in both rows  "\t Qstart-1  Hstart-1  identity"   (1-based residues; identity = %d of Hsp_identity / Hsp_align-len x 100, a double, truncated)
//   the file is read by read_constraint_list_file: entries with weight 0 are skipped, add_entry2list adds the entry and its mirror (weight = max, CONS += 0).
// The BLAST search itself is the external program: this port is T-Coffee's side of the interface, fed with the XML BLAST wrote.
// =====================================================================================================================

public static partial class TCoffeeAligner
{
    /// <summary>
    /// <c>xml2tag_list</c> of tc_generic_method.pl: the bodies of the elements named <paramref name="tag"/> in document order. A body is the concatenation of everything between the opening and the closing tag
    /// (inner tags included, as text), a tag is recognised by the substrings <c>&lt;tag </c> / <c>&lt;tag&gt;</c> (open) and <c>/tag&gt;</c> (close) of a markup token.
    /// </summary>
    private static List<string> TcXmlTagList(string xml, string tag)
    {
        var bodies = new List<string>();
        var current = new System.Text.StringBuilder();
        bool inside = false;
        string open1 = "<" + tag + " ", open2 = "<" + tag + ">", close = "/" + tag + ">";
        int pos = 0;
        while (pos < xml.Length)
        {
            int lt = xml.IndexOf('<', pos);
            if (lt < 0) break;
            int gt = xml.IndexOf('>', lt + 1);
            if (gt < 0) break;
            if (gt == lt + 1) { pos = gt + 1; continue; } // "<>" is no token (the pattern needs a character in between)
            string text = xml[pos..lt];
            if (inside) current.Append(text);
            string markup = xml[lt..(gt + 1)];
            pos = gt + 1;
            if (markup.Contains(open1, StringComparison.Ordinal) || markup.Contains(open2, StringComparison.Ordinal))
            {
                inside = true;
            }
            else if (markup.Contains(close, StringComparison.Ordinal))
            {
                bodies.Add(current.ToString());
                current.Clear();
                inside = false;
            }
            else if (inside) current.Append(markup);
        }

        return bodies;
    }

    /// <summary>The library entries (1-based sequence and residue numbers, weight) of the blast_msa step for <c>blastp_msa</c>, in file order (an entry and its mirror are added by the caller).</summary>
    internal static List<(int S1, int R1, int S2, int R2, int W)> TcBlastMsaEntries(IReadOnlyList<string> names, IReadOnlyList<string> sequences, string xml)
    {
        var entries = new List<(int, int, int, int, int)>();
        var order = new Dictionary<string, int>();
        for (int i = 0; i < names.Count; i++) order.TryAdd(names[i], i);
        var iterations = TcXmlTagList(xml, "Iteration");
        for (int a = 0; a < iterations.Count && a < names.Count; a++)
        {
            int query = order[names[a]];
            foreach (string hit in TcXmlTagList(iterations[a], "Hit"))
            {
                string definition = TcXmlTagList(hit, "Hit_def").FirstOrDefault() ?? string.Empty;
                int target = order.TryGetValue(definition, out int j) ? j : 0; // Perl: an unknown name is order undef + 1 = 1
                foreach (string hsp in TcXmlTagList(hit, "Hsp"))
                {
                    string qseq = TcXmlTagList(hsp, "Hsp_qseq").FirstOrDefault() ?? string.Empty;
                    string hseq = TcXmlTagList(hsp, "Hsp_hseq").FirstOrDefault() ?? string.Empty;
                    int qstart = TcPerlInt(TcXmlTagList(hsp, "Hsp_query-from").FirstOrDefault());
                    int hstart = TcPerlInt(TcXmlTagList(hsp, "Hsp_hit-from").FirstOrDefault());
                    double alignLen = TcPerlInt(TcXmlTagList(hsp, "Hsp_align-len").FirstOrDefault());
                    double identical = TcPerlInt(TcXmlTagList(hsp, "Hsp_identity").FirstOrDefault());
                    int identity = (int)(alignLen == 0 ? 0 : identical / alignLen * 100);
                    for (int c = 0; c < qseq.Length; c++)
                    {
                        int n = 0;
                        if (qseq[c] != '-') { n++; qstart++; }
                        if (c < hseq.Length && hseq[c] != '-') { n++; hstart++; }
                        if (n == 2 && identity > 0 && qstart - 1 >= 1 && hstart - 1 >= 1 && qstart - 1 <= sequences[query].Length && hstart - 1 <= sequences[target].Length)
                            entries.Add((query, qstart - 1, target, hstart - 1, identity));
                    }
                }
            }
        }

        return entries;
    }

    /// <summary>
    /// <c>ncbi_tblastx_xml2lib_file</c> (tc_generic_method.pl): per Iteration (query = <c>Iteration_iter-num</c>) and per Hit of a LATER sequence (<c>Hit_accession + 1 &gt; query</c>) the entries of every nucleotide of the query range of every HSP (reverse
    /// frames: <c>rqe = qlen - qs + 1</c>, <c>rqs = qlen - qe + 1</c>, hits likewise), <c>(e, rhs + d)</c> for <c>e = rqs .. rqe - 1</c>, weight <c>int (identity x 100 / align length)</c>.
    /// </summary>
    internal static List<(int S1, int R1, int S2, int R2, int W)> TcTblastxEntries(IReadOnlyList<string> sequences, string xml)
    {
        var entries = new List<(int, int, int, int, int)>();
        foreach (string iteration in TcXmlTagList(xml, "Iteration"))
        {
            int qindex = TcPerlInt(TcXmlTagList(iteration, "Iteration_iter-num").FirstOrDefault());
            int qlen = TcPerlInt(TcXmlTagList(iteration, "Iteration_query-len").FirstOrDefault());
            foreach (string hit in TcXmlTagList(iteration, "Hit"))
            {
                int hindex = TcPerlInt(TcXmlTagList(hit, "Hit_accession").FirstOrDefault()) + 1;
                if (hindex <= qindex) continue;
                int hlen = TcPerlInt(TcXmlTagList(hit, "Hit_len").FirstOrDefault());
                foreach (string hsp in TcXmlTagList(hit, "Hsp"))
                {
                    int qs = TcPerlInt(TcXmlTagList(hsp, "Hsp_query-from").FirstOrDefault()), qe = TcPerlInt(TcXmlTagList(hsp, "Hsp_query-to").FirstOrDefault()), qf = TcPerlInt(TcXmlTagList(hsp, "Hsp_query-frame").FirstOrDefault());
                    int hs = TcPerlInt(TcXmlTagList(hsp, "Hsp_hit-from").FirstOrDefault()), he = TcPerlInt(TcXmlTagList(hsp, "Hsp_hit-to").FirstOrDefault()), hf = TcPerlInt(TcXmlTagList(hsp, "Hsp_hit-frame").FirstOrDefault());
                    double identical = TcPerlInt(TcXmlTagList(hsp, "Hsp_identity").FirstOrDefault()), alignLen = TcPerlInt(TcXmlTagList(hsp, "Hsp_align-len").FirstOrDefault());
                    int weight = (int)(identical * 100 / alignLen);
                    int rqs = qs, rqe = qe, rhs = hs;
                    if (qf <= 0) { rqe = qlen - qs + 1; rqs = qlen - qe + 1; }
                    if (hf <= 0) rhs = hlen - he + 1;
                    for (int d = 0, e = rqs; e < rqe; e++, d++)
                        if (weight > 0 && e > 0 && rhs + d > 0 && e <= sequences[qindex - 1].Length && rhs + d <= sequences[hindex - 1].Length)
                            entries.Add((qindex - 1, e, hindex - 1, rhs + d, weight));
                }
            }
        }

        return entries;
    }

    /// <summary>
    /// The library <c>t_coffee -method tblastx_msa</c> builds from the XML of the all-against-all <c>tblastx</c> search (A2-9): see <see cref="TCoffeePairMethod.TblastxMsa"/>, printed as <see cref="TCoffeeBlastMsaLibraryEntries"/>.
    /// </summary>
    public static IReadOnlyList<TCoffeeLibraryEntry> TCoffeeTblastxLibraryEntries(IReadOnlyList<string> sequences, string blastXml)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(blastXml);
        foreach (var s in sequences)
            if (string.IsNullOrEmpty(s)) throw new ArgumentException("Empty sequences are not supported.", nameof(sequences));
        var cl = new TcConstraintList(sequences.Select(s => s.Length).ToArray());
        foreach (var (s1, r1, s2, r2, w) in TcTblastxEntries(sequences, blastXml)) cl.AddBoth(s1, r1, s2, r2, w, 0);
        var list = new List<TCoffeeLibraryEntry>();
        for (int s = 0; s < sequences.Count; s++)
            for (int r = 1; r <= cl.Length(s); r++)
                foreach (var e in cl.Entries(s, r))
                    if (s < e.Seq2 && e.Res2 > 0) list.Add(new TCoffeeLibraryEntry(s, r - 1, e.Seq2, e.Res2 - 1, e.Weight, e.Count));
        return list.OrderBy(e => e.SeqA).ThenBy(e => e.SeqB).ThenBy(e => e.PosA).ThenBy(e => e.PosB).ToList();
    }

    /// <summary>What the pair jobs of the BLAST methods need: the sequence names in the program's order and the XML of the search.</summary>
    internal sealed class TcBlastInput
    {
        public string[] Names = Array.Empty<string>();
        public string Xml = string.Empty;
    }

    private static int TcPerlInt(string? text) => int.TryParse(text?.Trim(), System.Globalization.NumberStyles.AllowLeadingSign, System.Globalization.CultureInfo.InvariantCulture, out int v) ? v : 0;

    /// <summary>
    /// The library <c>t_coffee -method blastp_msa</c> builds from the XML of the all-against-all BLAST search of the sequences (A2-9): the residue pairs of the high scoring pairs with the HSP's percent identity (integer, truncated)
    /// as the weight, an entry and its mirror, in the printed format of <see cref="TCoffeePrimaryLibraryEntries"/> (0-based residues, first sequence index smaller than the second; self hits are not printed).
    /// </summary>
    /// <param name="names">The sequence names as they stand in the BLAST database / query (the hits are matched by name).</param>
    /// <param name="sequences">The sequences.</param>
    /// <param name="blastXml">The XML (<c>-m7</c> / <c>-outfmt 5</c>) of the search; one <c>Iteration</c> per sequence in input order.</param>
    public static IReadOnlyList<TCoffeeLibraryEntry> TCoffeeBlastMsaLibraryEntries(IReadOnlyList<string> names, IReadOnlyList<string> sequences, string blastXml)
    {
        ArgumentNullException.ThrowIfNull(names);
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(blastXml);
        if (names.Count != sequences.Count) throw new ArgumentException("names must have one entry per sequence.", nameof(names));
        foreach (var s in sequences)
            if (string.IsNullOrEmpty(s)) throw new ArgumentException("Empty sequences are not supported.", nameof(sequences));
        var cl = new TcConstraintList(sequences.Select(s => s.Length).ToArray());
        foreach (var (s1, r1, s2, r2, w) in TcBlastMsaEntries(names, sequences, blastXml)) cl.AddBoth(s1, r1, s2, r2, w, 0);
        var list = new List<TCoffeeLibraryEntry>();
        for (int s = 0; s < sequences.Count; s++)
            for (int r = 1; r <= cl.Length(s); r++)
                foreach (var e in cl.Entries(s, r))
                    if (s < e.Seq2 && e.Res2 > 0) list.Add(new TCoffeeLibraryEntry(s, r - 1, e.Seq2, e.Res2 - 1, e.Weight, e.Count));
        return list.OrderBy(e => e.SeqA).ThenBy(e => e.SeqB).ThenBy(e => e.PosA).ThenBy(e => e.PosB).ToList();
    }
}
