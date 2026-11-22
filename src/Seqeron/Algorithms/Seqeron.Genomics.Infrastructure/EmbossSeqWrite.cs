#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the writers keep the control flow of ajseqwrite.c

using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>One sequence output format of <c>seqOutFormat</c> (<c>ajseqwrite.c</c>): the name, whether it is an alias, the kinds of sequence it takes, whether the writer keeps the sequences and writes them when the output is closed, and the writer.</summary>
internal sealed record EmbossSeqOutFormat(string Name, bool Alias, bool Save, bool Nucleotide, bool Protein, bool Features, string Writer);

/// <summary>The state of the sequence writers that lasts for a program run: the counters of the names <c>gnl|BL_ORD_ID|N</c> of the ncbi and gifasta writers (a <c>static</c> of each writer).</summary>
internal sealed class EmbossSeqWriteRun
{
    public int NcbiOrdinal, GifastaOrdinal;
}

/// <summary>
/// The sequence writers of EMBOSS 6.6.0 <c>ajseqwrite.c</c> (<c>seqOutFormat</c>, <c>seqSeqFormat</c>, <c>seqWriteSeq</c> and the <c>seqWrite*</c> functions) for the formats that write one entry per sequence from the
/// name, description, accession, type and residues: fasta, gcg, nbrf, codata, strider, staden, text, ig, iguspto, gde, fitch, acedb, mase. A5-8c3h, B09 F252.
/// </summary>
internal static class EmbossSeqWrite
{
    /// <summary>The table <c>seqOutFormat</c> (name, alias, save, nucleotide, protein, feature, writer).</summary>
    internal static readonly EmbossSeqOutFormat[] Table =
    [
        new("unknown", false, false, true, true, false, "Fasta"), new("gcg", false, false, true, true, false, "Gcg"), new("gcg8", true, false, true, true, false, "Gcg"),
        new("embl", false, false, true, false, true, "Emblnew"), new("emold", true, false, true, false, true, "Embl"), new("emblold", true, false, true, false, true, "Embl"), new("em", true, false, true, false, true, "Emblnew"),
        new("emblnew", true, false, true, false, true, "Emblnew"), new("swiss", false, false, false, true, true, "Swissnew"), new("swold", true, false, false, true, true, "Swiss"),
        new("swissold", true, false, false, true, true, "Swiss"), new("swissprotold", true, false, false, true, true, "Swiss"), new("sw", true, false, false, true, true, "Swissnew"),
        new("swissprot", true, false, false, true, true, "Swissnew"), new("swissnew", true, false, false, true, true, "Swissnew"), new("swnew", true, false, false, true, true, "Swissnew"),
        new("swissprotnew", true, false, false, true, true, "Swissnew"), new("uniprot", true, false, false, true, true, "Swissnew"), new("fasta", false, false, true, true, false, "Fasta"),
        new("pearson", true, false, true, true, false, "Fasta"), new("ncbi", false, false, true, true, false, "Ncbi"), new("gifasta", false, false, true, true, false, "Gifasta"), new("gde", false, false, true, true, false, "Gde"),
        new("nib", false, false, true, false, false, "Nibble"), new("nibble", false, false, true, false, false, "Nibble"), new("nbrf", false, false, true, true, true, "Nbrf"), new("pir", true, false, true, true, true, "Nbrf"),
        new("genbank", false, false, true, true, true, "Genbank"), new("gb", true, false, true, false, true, "Genbank"), new("ddbj", true, false, true, false, true, "Genbank"), new("genpept", false, false, true, true, true, "Genpept"),
        new("refseq", false, false, true, false, true, "Refseq"), new("refseqp", false, false, false, true, true, "Refseqp"), new("gff2", false, false, true, true, true, "Gff2"), new("gff3", false, false, true, true, true, "Gff3"),
        new("gff", true, false, true, true, true, "Gff3"), new("ig", false, false, true, true, false, "Ig"), new("iguspto", false, false, true, true, false, "Iguspto"), new("codata", false, false, true, true, false, "Codata"),
        new("strider", false, false, true, true, false, "Strider"), new("acedb", false, false, true, true, false, "Acedb"), new("experiment", false, false, true, true, false, "Experiment"), new("staden", false, false, true, true, false, "Staden"),
        new("text", false, false, true, true, false, "Text"), new("plain", true, false, true, true, false, "Text"), new("raw", true, false, true, true, false, "Text"), new("fitch", false, false, true, true, false, "Fitch"),
        new("msf", false, true, true, true, false, "Msf"), new("clustal", false, true, true, true, false, "Clustal"), new("aln", true, true, true, true, false, "Clustal"), new("selex", false, true, true, true, false, "Selex"),
        new("phylip", false, true, true, true, false, "Phylip"), new("phylipnon", false, true, true, true, false, "Phylipnon"), new("phylip3", true, true, true, true, false, "Phylipnon"), new("asn1", false, false, true, true, false, "Asn1"),
        new("hennig86", false, true, true, true, false, "Hennig86"), new("mega", false, true, true, true, false, "Mega"), new("meganon", false, true, true, true, false, "Meganon"), new("nexus", false, true, true, true, false, "Nexus"),
        new("paup", true, true, true, true, false, "Nexus"), new("nexusnon", false, true, true, true, false, "Nexusnon"), new("paupnon", true, true, true, true, false, "Nexusnon"), new("jackknifer", false, true, true, true, false, "Jackknifer"),
        new("jackknifernon", false, true, true, true, false, "Jackknifernon"), new("treecon", false, true, true, true, false, "Treecon"), new("mase", false, false, true, true, false, "Mase"), new("dasdna", false, false, true, false, false, "Dasdna"),
        new("das", false, false, true, true, false, "Dasseq"), new("fastq-sanger", false, false, true, false, false, "FastqSanger"), new("fastq", true, false, true, false, false, "FastqSanger"),
        new("fastq-illumina", false, false, true, false, false, "FastqIllumina"), new("fastq-solexa", false, false, true, false, false, "FastqSolexa"), new("sam", false, false, true, true, false, "Sam"), new("bam", false, false, true, true, false, "Bam"),
        new("scf", false, false, true, true, false, "Scf"), new("debug", false, false, true, true, false, "Debug"),
    ];

    /// <summary><c>seqoutFindOutFormat</c>: the format of a name (any case); null when there is none.</summary>
    internal static EmbossSeqOutFormat? Find(string name)
    {
        foreach (var f in Table)
        {
            if (string.Equals(f.Name, name, StringComparison.OrdinalIgnoreCase))
                return f;
        }

        return null;
    }

    /// <summary>The writers of this port.</summary>
    internal static bool IsPorted(EmbossSeqOutFormat f) => f.Writer is "Asn1" or "Dasdna" or "Dasseq" or "Nibble" or "Ncbi" or "Gifasta" or "Fasta" or "Gcg" or "Nbrf" or "Codata" or "Strider" or "Staden" or "Text" or "Ig" or "Iguspto" or "Gde" or "Fitch" or "Acedb" or "Mase";

    /// <summary>The <c>SeqSFormat</c> of the library: the layout of <c>seqWriteSeq</c> (the fields the ported formats set; the others stay at the defaults of <c>seqSeqFormat</c>: no name, no right number, no tab, no padding).</summary>
    private sealed class SeqFormat
    {
        public int Spacer, Width = 50, NumWidth, LinePos;
        public bool NumLeft, SkipBefore;
        public string EndStr = "";
    }

    private const int SpaceAll = -9;

    // seqSeqFormat
    private static SeqFormat NewFormat(int seqlen)
    {
        int j = 1;
        for (int i = seqlen; i != 0; i /= 10)
            j++;
        return new SeqFormat { NumWidth = j };
    }

    private static string Pad(string s, int width) => s.Length >= width ? s : new string(' ', width - s.Length) + s;

    // seqWriteSeq
    private static void WriteSeq(StringBuilder o, string seq, SeqFormat sf)
    {
        const int maxSeqWidth = 250;
        int seqlen = seq.Length;
        int width = Math.Min(sf.Width, maxSeqWidth);
        int l1 = sf.LinePos, l = 0, i = 0, ibase = 1;
        var s = new StringBuilder();
        while (i < seqlen)
        {
            if (l1 == 0)
            {
                if (sf.SkipBefore)
                    o.Append('\n');
                if (sf.NumLeft)
                    o.Append(Pad(ibase.ToString(CultureInfo.InvariantCulture), sf.NumWidth));
            }

            l1++;
            if (sf.Spacer == SpaceAll || (sf.Spacer != 0 && (l + 1) % sf.Spacer == 1))
            {
                s.Append(' ');
                l++;
            }

            ibase++;
            s.Append(seq[i++]);
            l++;
            if (l1 == width || i == seqlen)
            {
                l = 0;
                l1 = 0;
                o.Append(s);
                s.Clear();
                if (i == seqlen)
                    o.Append(sf.EndStr);
                o.Append('\n');
            }
        }
    }

    // ajSeqGapS(&seq, '.'): the characters of "-~.? " become the gap character
    private static string GapS(string seq, char gap)
    {
        var sb = new StringBuilder(seq.Length);
        foreach (char c in seq)
            sb.Append(c is '-' or '~' or '.' or '?' or ' ' ? gap : c);
        return sb.ToString();
    }

    // ajSeqoutGetCheckgcg
    private static int CheckGcg(string seq)
    {
        long check = 0, count = 0;
        foreach (char c in seq)
        {
            count++;
            check += count * char.ToUpperInvariant(c);
            if (count == 57)
                count = 0;
        }

        return (int)(check % 10000);
    }

    // seqWriteNibble: the file is rewound before every sequence (ajFileSeek(outseq->File, 0L, SEEK_SET)), so the entry of a later sequence overwrites the start of the file and the rest of a longer earlier entry stays;
    // the magic number 0x6BE9ED3A and the length are written little endian, the residues as 4-bit codes (A 2, C 1, G 3, T and U 0, any other 4), the first of two in the high nibble
    private static void WriteNibble(StringBuilder o, string seq)
    {
        var bytes = new StringBuilder();
        bytes.Append((char)0x3A).Append((char)0xED).Append((char)0xE9).Append((char)0x6B);
        int len = seq.Length;
        bytes.Append((char)(len & 0xFF)).Append((char)((len >> 8) & 0xFF)).Append((char)((len >> 16) & 0xFF)).Append((char)((len >> 24) & 0xFF));
        int current = 0;
        bool fullbyte = true;
        foreach (char c in seq)
        {
            if (c == '\0')
                break;
            current |= c switch { 'A' or 'a' => 2, 'C' or 'c' => 1, 'G' or 'g' => 3, 'T' or 't' or 'U' or 'u' => 0, _ => 4 };
            fullbyte = !fullbyte;
            if (fullbyte)
            {
                bytes.Append((char)(current & 0xFF));
                current = 0;
            }
            else
            {
                current <<= 4;
            }
        }

        if (!fullbyte)
            bytes.Append((char)(current & 0xFF));
        if (o.Length > bytes.Length)
            o.Remove(0, bytes.Length).Insert(0, bytes.ToString());
        else
            o.Clear().Append(bytes);
    }

    /// <summary>The footer of a DAS file (<c>seqCleanDasdna</c> / <c>seqCleanDasseq</c>, called when the file open at the end of the run is closed), empty for the other formats.</summary>
    internal static string Footer(EmbossSeqOutFormat f) => f.Writer switch { "Dasdna" => "</DASDNA>\n", "Dasseq" => "</DASSEQUENCE>\n", _ => "" };

    // seqWriteDasdna (dasseq false) / seqWriteDasseq (dasseq true)
    private static void WriteDas(StringBuilder o, string name, string sv, string type, string seq, int offset, int count, bool dasseq)
    {
        int ilen = seq.Length;
        if (count == 0)
        {
            o.Append("<?xml version=\"1.0\" standalone=\"no\"?>\n");
            o.Append(dasseq ? "<!DOCTYPE DASSEQUENCE SYSTEM \"http://www.biodas.org/dtd/dassequence.dtd\">\n" : "<!DOCTYPE DASDNA SYSTEM \"http://www.biodas.org/dtd/dasdna.dtd\">\n<DASDNA>\n");
        }

        if (dasseq)
        {
            o.Append("<DASSEQUENCE>\n");
            o.Append(CultureInfo.InvariantCulture, $"  <SEQUENCE id=\"{name}\" start=\"{1 + offset}\" stop=\"{ilen + offset}\"\n");
            o.Append(type.Length > 0 && type[0] == 'P' ? "               moltype=\"Protein\"" : "               moltype=\"DNA\"");
            o.Append(sv.Length > 0 ? $" version=\"{sv}\">\n" : " version=\"0.0\">\n");
        }
        else
        {
            o.Append(CultureInfo.InvariantCulture, $"  <SEQUENCE id=\"{name}\" start=\"{1 + offset}\" stop=\"{ilen + offset}\" version=\"{sv}\">\n");
            o.Append(CultureInfo.InvariantCulture, $"    <DNA length=\"{ilen}\">\n");
        }

        for (int i = 0; i < ilen; i += 50)
            o.Append("      ").Append(seq, i, Math.Min(50, ilen - i)).Append('\n');
        if (!dasseq)
            o.Append("    </DNA>\n");
        o.Append("  </SEQUENCE>\n");
    }

    // seqWriteAsn1: the sequence as an ASN.1 value, the lines of the residues of 78 characters (the first one after the 17 characters of `iupacna "`)
    private static void WriteAsn1(StringBuilder o, string desc, string type, string seq)
    {
        o.Append("  seq {\n");
        o.Append("    id { local id 1 },\n");
        o.Append("    descr { title \"").Append(desc).Append("\" },\n");
        o.Append("    inst {\n");
        bool prot = type.Length > 0 && type[0] == 'P';
        if (prot)
            o.Append("      repr raw, mol aa, length ").Append(seq.Length.ToString(CultureInfo.InvariantCulture)).Append(", topology linear,\n {\n");
        else
            o.Append("      repr raw, mol dna, length ").Append(seq.Length.ToString(CultureInfo.InvariantCulture)).Append(", topology linear,\n");
        o.Append("      seq-data\n");
        o.Append(prot ? "        iupacaa \"" : "        iupacna \"");
        WriteSeq(o, seq, new SeqFormat { LinePos = 17, Spacer = 0, Width = 78, EndStr = "\"\n      } } ," });
    }

    private static readonly string[] NcbiDatabases = ["gb", "emb", "dbj", "sp", "ref", "pir", "prf", "tpd", "tpe", "tpg"];

    // seqWriteNcbi / seqWriteGifasta (the writers differ in the GI part and in the databases lcl and bbs)
    private static void WriteNcbi(StringBuilder o, string name, string desc, string acc, string sv, string seq, string gi, string db, string setdb, bool gifasta, EmbossSeqWriteRun? run)
    {
        run ??= new EmbossSeqWriteRun();
        string dbname = setdb.Length > 0 ? setdb : db.Length > 0 ? db : "unk";
        string version = sv.Length > 0 ? sv : acc;
        bool known = Array.IndexOf(NcbiDatabases, dbname) >= 0;
        if (gi.Length > 0 && db.Length == 0 && string.Equals(gi, name, StringComparison.OrdinalIgnoreCase))
        {
            o.Append(">gi|").Append(gi);
        }
        else
        {
            o.Append('>');
            if (gi.Length > 0)
                o.Append("gi|").Append(gi).Append('|');
            else if (gifasta)
                o.Append("gi|000000|");
            if (known)
            {
                o.Append(dbname).Append('|').Append(version).Append('|');
            }
            else if (string.Equals(dbname, "lcl", StringComparison.OrdinalIgnoreCase))
            {
                o.Append(gifasta ? dbname : "lcl").Append('|');
            }
            else if (!gifasta && string.Equals(dbname, "bbs", StringComparison.OrdinalIgnoreCase))
            {
                o.Append("bbs|");
            }
            else if (string.Equals(dbname, "BL_ORD_ID", StringComparison.OrdinalIgnoreCase))
            {
                int n = gifasta ? run.GifastaOrdinal++ : run.NcbiOrdinal++;
                o.Append("gnl|").Append(dbname).Append('|').Append(n).Append(' ');
            }
            else
            {
                o.Append("gnl|").Append(dbname).Append('|');
            }

            if (!string.Equals(version, name, StringComparison.OrdinalIgnoreCase))
                o.Append(name);
        }

        if (version.Length > 0 && !known)
            o.Append(" (").Append(version).Append(')');
        if (desc.Length > 0)
            o.Append(' ').Append(desc);
        o.Append('\n');
        WriteLines60(o, seq);
    }

    private static void WriteLines60(StringBuilder o, string seq)
    {
        for (int i = 0; i < seq.Length; i += 60)
            o.Append(seq, i, Math.Min(60, seq.Length - i)).Append('\n');
    }

    /// <summary>Writes one sequence in a format of this port (the entry of a sequence the writer was given: <paramref name="count"/> is the number of sequences written before it).</summary>
    /// <param name="o">The output file.</param>
    /// <param name="format">The format (<see cref="Find"/>).</param>
    /// <param name="name">The name of the sequence (the default name is set by the caller).</param>
    /// <param name="desc">The description.</param>
    /// <param name="acc">The accession.</param>
    /// <param name="sv">The sequence version.</param>
    /// <param name="type">The type <c>N</c> / <c>P</c> or empty.</param>
    /// <param name="seq">The residues.</param>
    /// <param name="count">The number of sequences written to the file before this one.</param>
    /// <exception cref="InvalidOperationException">The writer is one on which the reference program dies by a fault (iguspto).</exception>
    /// <param name="runDate">The date of the run (the header of mase).</param>
    /// <param name="fullDescMulti">The comment lines of the full description (<c>Fulldesc-&gt;Multi</c>); null when the reader made no full description.</param>
    /// <param name="offset">The offset of the sequence (<c>outseq->Offset</c>: the start of the range of the sequence minus one; the DAS writers print it).</param>
    /// <param name="headerCount">The number of sequences written before this one to ALL the files of the run (<c>outseq->Count</c>, which <c>-ossingle</c> does not reset): the DAS writers write their header when it is 0; -1 = <paramref name="count"/>.</param>
    /// <param name="gi">The GI number of the sequence (the ncbi and gifasta writers).</param>
    /// <param name="db">The database the sequence was read from (<c>Db</c>).</param>
    /// <param name="setdb">The database name of the sequence (<c>Setdb</c>: from the header of an NCBI FASTA entry, else <c>-osdbname</c>).</param>
    /// <param name="run">The state of the writers that lasts for the run (the counter of the <c>BL_ORD_ID</c> names of ncbi and gifasta).</param>
    internal static void Write(StringBuilder o, EmbossSeqOutFormat format, string name, string desc, string acc, string sv, string type, string seq, int count, DateTime runDate, IReadOnlyList<string>? fullDescMulti = null,
        string gi = "", string db = "", string setdb = "", EmbossSeqWriteRun? run = null, int offset = 0, int headerCount = -1)
    {
        bool prot = type.Length > 0 && type[0] == 'P';
        switch (format.Writer)
        {
            case "Asn1":
                WriteAsn1(o, desc, type, seq);
                break;
            case "Dasdna":
                WriteDas(o, name, sv, type, seq, offset, headerCount >= 0 ? headerCount : count, false);
                break;
            case "Dasseq":
                WriteDas(o, name, sv, type, seq, offset, headerCount >= 0 ? headerCount : count, true);
                break;
            case "Nibble":
                WriteNibble(o, seq);
                break;
            case "Ncbi":
                WriteNcbi(o, name, desc, acc, sv, seq, gi, db, setdb, false, run);
                break;
            case "Gifasta":
                WriteNcbi(o, name, desc, acc, sv, seq, gi, db, setdb, true, run);
                break;
            case "Fasta":
                o.Append('>').Append(name);
                if (sv.Length > 0)
                    o.Append(' ').Append(sv);
                else if (acc.Length > 0)
                    o.Append(' ').Append(acc);
                if (desc.Length > 0)
                    o.Append(' ').Append(desc);
                o.Append('\n');
                WriteLines60(o, seq);
                break;
            case "Gcg":
            {
                int ilen = seq.Length;
                char ctype = 'N';
                if (type.Length == 0)
                    o.Append("!!NA_SEQUENCE 1.0\n\n");
                else if (prot)
                {
                    o.Append("!!AA_SEQUENCE 1.0\n\n");
                    ctype = 'P';
                }
                else
                    o.Append("!!NA_SEQUENCE 1.0\n\n");
                string gapped = GapS(seq, '.');
                int check = CheckGcg(gapped);
                if (desc.Length > 0)
                    o.Append(desc).Append("\n\n");
                o.Append(CultureInfo.InvariantCulture, $"{name}  Length: {ilen}  Type: {ctype}  Check: {check,4} ..\n");
                var sf = NewFormat(gapped.Length);
                sf.Spacer = 11;
                sf.NumLeft = true;
                sf.SkipBefore = true;
                sf.EndStr = "\n";
                WriteSeq(o, gapped, sf);
                break;
            }

            case "Nbrf":
            {
                o.Append(type.Length == 0 ? ">D1;" : prot ? ">P1;" : ">D1;").Append(name).Append('\n');
                o.Append(CultureInfo.InvariantCulture, $"{desc}, {seq.Length} bases\n");
                var sf = NewFormat(seq.Length);
                sf.Spacer = 11;
                sf.EndStr = "*\n";
                WriteSeq(o, seq, sf);
                break;
            }

            case "Codata":
            {
                o.Append("ENTRY           ").Append(name).Append(" \n");
                if (desc.Length > 0)
                    o.Append(CultureInfo.InvariantCulture, $"TITLE           {desc}, {seq.Length} bases\n");
                if (acc.Length > 0)
                    o.Append("ACCESSION       ").Append(acc).Append('\n');
                o.Append("SEQUENCE        \n");
                var sf = NewFormat(seq.Length);
                sf.NumWidth = 7;
                sf.Width = 30;
                sf.NumLeft = true;
                sf.Spacer = SpaceAll;
                sf.EndStr = "\n///";
                for (int j = 0; j <= sf.NumWidth; j++)
                    o.Append(' ');
                for (int j = 5; j <= sf.Width; j += 5)
                    o.Append(CultureInfo.InvariantCulture, $"{j,10}");
                o.Append('\n');
                WriteSeq(o, seq, sf);
                break;
            }

            case "Strider":
            {
                o.Append("; ### from DNA Strider ;-)\n");
                o.Append(CultureInfo.InvariantCulture, $"; DNA sequence  {name}, {seq.Length} bases\n;\n");
                var sf = NewFormat(seq.Length);
                sf.EndStr = "\n//";
                WriteSeq(o, seq, sf);
                break;
            }

            case "Staden":
            {
                o.Append('<').Append(name).Append("---->\n");
                var sf = NewFormat(seq.Length);
                sf.Width = 60;
                WriteSeq(o, seq, sf);
                break;
            }

            case "Text":
                WriteSeq(o, seq, NewFormat(seq.Length));
                break;
            case "Ig":
            case "Iguspto":
            {
                if (format.Writer == "Iguspto")
                {
                    // seqWriteIguspto reads outseq->Fulldesc->Multi: only the swiss and iguspto readers make a Fulldesc, for any other sequence the process dies by a segmentation fault
                    if (fullDescMulti is null)
                        throw new InvalidOperationException("The writer of the format iguspto dereferences a null pointer (the sequence has no full description: the reference program dies by a segmentation fault).");
                    if (count > 0)
                        o.Append("\n\f\n");
                }

                if (format.Writer == "Iguspto" && fullDescMulti!.Count > 0)
                {
                    foreach (string line in fullDescMulti)
                        o.Append("; ").Append(line).Append('\n');
                }
                else
                {
                    o.Append(';').Append(desc);
                    o.Append(CultureInfo.InvariantCulture, $", {seq.Length} {(prot ? "aa" : "bases")}\n");
                }

                o.Append(name).Append('\n');
                var sf = NewFormat(seq.Length);
                sf.EndStr = "1";
                if (format.Writer == "Iguspto")
                    sf.Width = 60;
                WriteSeq(o, seq, sf);
                break;
            }

            case "Gde":
                o.Append('#').Append(name).Append('\n');
                WriteLines60(o, seq);
                break;
            case "Fitch":
            {
                o.Append(CultureInfo.InvariantCulture, $"{name}, {seq.Length} bases\n");
                var sf = NewFormat(seq.Length);
                sf.Spacer = 4;
                sf.Width = 60;
                WriteSeq(o, seq, sf);
                break;
            }

            case "Acedb":
            {
                o.Append(prot ? "Peptide : \"" : "DNA : \"").Append(name).Append("\"\n");
                var sf = NewFormat(seq.Length);
                sf.EndStr = "\n";
                WriteSeq(o, seq, sf);
                break;
            }

            case "Mase":
                if (count == 0)
                    o.Append(";;Written by EMBOSS on ").Append(runDate.ToString("ddd ", CultureInfo.InvariantCulture)).Append(runDate.Day.ToString(CultureInfo.InvariantCulture).PadLeft(2)).Append(runDate.ToString(" MMM yyyy HH:mm:ss", CultureInfo.InvariantCulture)).Append('\n');
                o.Append(';').Append(desc).Append('\n').Append(name).Append('\n');
                WriteLines60(o, seq);
                break;
            default:
                throw new NotSupportedException($"The sequence output format '{format.Name}' is not written by this port.");
        }
    }
}
