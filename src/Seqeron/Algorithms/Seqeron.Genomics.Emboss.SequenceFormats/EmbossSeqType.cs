// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854
using System.Text;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>
/// The sequence types of the EMBOSS 6.6.0 AJAX library (<c>ajseqtype.c</c>; B09 WP61, F280): the table <c>seqType</c> (name, gaps, ambiguity codes, nucleotide / protein, the characters converted on input, the description),
/// <c>ajSeqTypeCheckIn</c> (the check of a sequence read by an ACD sequence input: gap characters converted or removed, <c>-snucleotide</c> / <c>-sprotein</c>, the bad characters fixed to N / X, the messages
/// "Sequence is not nucleic" / "Sequence is not a protein" / "Sequence must be ...: found bad character"), <c>ajSeqTypeCheckS</c> (the type of an output sequence, <c>seqClone</c>), <c>ajSeqIsNuc</c> / <c>ajSeqIsProt</c>.
/// </summary>
public static class EmbossSeqType
{
    public sealed record Def(string Name, bool Gaps, bool Ambig, char Kind, string ConvertFrom, string ConvertTo, string Good, string Description)
    {
        /// <summary>The set of the characters a sequence of the type may hold, both cases (<c>ajStrGetfilterCase</c>).</summary>
        public bool[] Filter { get; } = MakeFilter(Good);
    }

    private static bool[] MakeFilter(string chars)
    {
        var f = new bool[256];
        foreach (char c in chars)
        {
            f[char.ToLowerInvariant(c)] = true;
            f[char.ToUpperInvariant(c)] = true;
        }

        return f;
    }

    private const string ProtPure = "ACDEFGHIKLMNPQRSTVWY", ProtAmbig = "BJOUXZ?", ProtStop = "*", NucPure = "ACGTU", NucAmbig = "BDHKMNRSVWXY?", Gap = ".~-", Phylo = "?", GapTest = " .~-";

    private static string Dedup(string s)
    {
        var sb = new StringBuilder();
        foreach (char c in s)
        {
            if (!sb.ToString().Contains(c))
                sb.Append(c);
        }

        return sb.ToString();
    }

    private static readonly string Any = Dedup(ProtPure + ProtAmbig + ProtStop + NucPure + NucAmbig);
    private static readonly string AnyGap = Dedup(ProtPure + ProtAmbig + ProtStop + NucPure + NucAmbig + Gap);
    private static readonly string Nuc = Dedup(NucPure + NucAmbig);
    private static readonly string NucGap = Dedup(NucPure + NucAmbig + Gap);
    private static readonly string NucGapPhylo = Dedup(NucPure + NucAmbig + Phylo + Gap);
    private static readonly string NucPureSet = Dedup(NucPure);
    private static readonly string Prot = Dedup(ProtPure + ProtAmbig);
    private static readonly string ProtGap = Dedup(ProtPure + ProtAmbig + Gap);
    private static readonly string ProtGapPhylo = Dedup(ProtPure + ProtAmbig + ProtStop + Phylo + Gap);
    private static readonly string ProtPureSet = Dedup(ProtPure);
    private static readonly string ProtStopSet = Dedup(ProtPure + ProtAmbig + ProtStop);
    private static readonly string ProtStopGap = Dedup(ProtPure + ProtAmbig + ProtStop + Gap);
    private static readonly string ProtAny = Dedup(ProtPure + ProtAmbig + ProtStop + Phylo + Gap);

    // the table seqType of ajseqtype.c (kind: A any, N nucleotide, P protein)
    private static readonly Def[] Table =
    [
        new("any", false, true, 'A', "?", "X", Any, "any valid sequence"),
        new("gapany", true, true, 'A', "?", "X", AnyGap, "any valid sequence with gaps"),
        new("dna", false, true, 'N', "?XxUu", "NNnTt", Nuc, "DNA sequence"),
        new("puredna", false, false, 'N', "Uu", "Tt", NucPureSet, "DNA sequence, bases ACGT only"),
        new("gapdna", true, true, 'N', "?XxUu", "NNnTt", NucGap, "DNA sequence with gaps"),
        new("gapdnaphylo", true, true, 'N', "Uu", "Tt", NucGapPhylo, "DNA sequence with gaps and queries"),
        new("rna", false, true, 'N', "?XxTt", "NNnUu", Nuc, "RNA sequence"),
        new("purerna", false, false, 'N', "Tt", "Uu", NucPureSet, "RNA sequence, bases ACGU only"),
        new("gaprna", true, true, 'N', "?XxTt", "NNnUu", NucGap, "RNA sequence with gaps"),
        new("gaprnaphylo", true, true, 'N', "Tt", "Uu", NucGapPhylo, "RNA sequence with gaps and queries"),
        new("nucleotide", false, true, 'N', "?Xx", "NNn", Nuc, "nucleotide sequence"),
        new("purenucleotide", false, false, 'N', "", "", NucPureSet, "nucleotide sequence, bases ACGTU only"),
        new("gapnucleotide", true, true, 'N', "?Xx", "NNn", NucGap, "nucleotide sequence with gaps"),
        new("gapnucleotidephylo", true, true, 'N', "", "", NucGapPhylo, "nucleotide sequence with gaps and queries"),
        new("gapnucleotidesimple", true, true, 'N', "BbDdHhKkMmRrSsVvWwXxYy?", "NnNnNnNnNnNnNnNnNnNnNnN", NucGap, "nucleotide sequence with gaps but only N for ambiguity"),
        new("protein", false, true, 'P', "?*", "XX", Prot, "protein sequence"),
        new("pureprotein", false, false, 'P', "", "", ProtPureSet, "protein sequence without BZ U X or *"),
        new("stopprotein", false, true, 'P', "?", "X", ProtStopSet, "protein sequence with possible stops"),
        new("gapprotein", true, true, 'P', "?*", "XX", ProtGap, "protein sequence with gaps"),
        new("gapstopprotein", true, true, 'P', "?", "X", ProtStopGap, "protein sequence with gaps and possible stops"),
        new("gapproteinphylo", true, true, 'P', "", "", ProtGapPhylo, "protein sequence with gaps, stops and queries"),
        new("proteinstandard", false, true, 'P', "?*UuJjOo", "XXXxXxXx", Prot, "protein sequence with no selenocysteine"),
        new("stopproteinstandard", false, true, 'P', "?UuJjOo", "XXxXxXx", ProtStopSet, "protein sequence with a possible stop but no selenocysteine"),
        new("gapproteinstandard", true, true, 'P', "?*UuJjOo", "XXXxXxXx", ProtGap, "protein sequence with gaps but no selenocysteine"),
        new("gapproteinsimple", true, true, 'P', "?*BbZzUuJjOo", "XXXxXxXxXxXx", ProtGap, "protein sequence with gaps but no selenocysteine"),
    ];

    /// <summary><c>seqFindType</c>: the type of a name (case-insensitive) or null.</summary>
    public static Def? Find(string name)
    {
        foreach (var d in Table)
        {
            if (string.Equals(d.Name, name, StringComparison.OrdinalIgnoreCase))
                return d;
        }

        return null;
    }

    // ajStrExchangeSetCC
    private static string ExchangeSet(string s, string from, string to)
    {
        if (from.Length == 0 || s.IndexOfAny(from.ToCharArray()) < 0)
            return s;
        var filter = new char[256];
        for (int i = 0; i < from.Length; i++)
            filter[from[i]] = i < to.Length ? to[i] : to[0];
        var sb = new StringBuilder(s.Length);
        foreach (char c in s)
            sb.Append(c < 256 && from.Contains(c) ? filter[c] : c);
        return sb.ToString();
    }

    /// <summary><c>ajSeqGapS</c> / <c>ajSeqGap</c>: the gap characters (blank . ~ -) become <c>-</c>.</summary>
    public static string GapConvert(string s) => ExchangeSet(s, GapTest, "----");

    /// <summary><c>ajStrRemoveGap</c>: everything that is no letter and no <c>*</c> is removed.</summary>
    public static string RemoveGap(string s)
    {
        var sb = new StringBuilder(s.Length);
        foreach (char c in s)
        {
            if (char.IsAsciiLetter(c) || c == '*')
                sb.Append(c);
        }

        return sb.ToString();
    }

    private static bool IsCharset(string s, bool[] filter)
    {
        if (s.Length == 0)
            return false;
        foreach (char c in s)
        {
            if (c >= 256 || !filter[c])
                return false;
        }

        return true;
    }

    // seqTypeTestS: the first character that is not in the set (0 when there is none or the string is empty)
    private static char TestS(string s, bool[] filter)
    {
        if (s.Length == 0 || IsCharset(s, filter))
            return '\0';
        foreach (char c in s)
        {
            if (c >= 256 || !filter[c])
                return c;
        }

        return '\0';
    }

    private static readonly bool[] FilterNucGap = MakeFilter(NucGap), FilterProtAny = MakeFilter(ProtAny);

    /// <summary><c>ajSeqIsNuc</c> of a sequence with the type <paramref name="type"/> (<c>N</c> / <c>P</c> / empty) and the residues <paramref name="seq"/>.</summary>
    public static bool IsNuc(string type, string seq) => type == "N" || (type != "P" && TestS(seq, FilterNucGap) == '\0');

    /// <summary><c>ajSeqsetIsNuc</c>: the type of the set (the type of its first member) is N, or the residues of the first member are nucleic (whatever its own type says).</summary>
    public static bool IsNucSet(string type, string seq) => type == "N" || TestS(seq, FilterNucGap) == '\0';

    /// <summary><c>ajSeqIsProt</c>.</summary>
    public static bool IsProt(string type, string seq) => type == "P" || (type != "N" && TestS(seq, FilterProtAny) == '\0');

    /// <summary><c>ajSeqTypeCheckS</c>: the type of an output sequence (<c>seqClone</c>): the gaps of a gapped type become <c>-</c>, the characters that are no letters are removed for the other types; when every character is of the type the
    /// characters of the type that are converted on input are converted (an empty type: only the gaps).</summary>
    /// <exception cref="EmbossFatalException">The type is unknown (<c>ajDie</c>).</exception>
    public static string CheckS(string seq, string type)
    {
        if (type.Length == 0)
            return GapConvert(seq);
        var def = Find(type) ?? throw new EmbossFatalException($"Sequence type '{type}' unknown");
        seq = def.Gaps ? GapConvert(seq) : RemoveGap(seq);
        if (IsCharset(seq, def.Filter) && def.ConvertFrom.Length > 0)
            seq = ExchangeSet(seq, def.ConvertFrom, def.ConvertTo);
        return seq;
    }

    // seqTypeStopTrimS: a trailing stop is cut
    private static string StopTrim(string s) => s.Length > 0 && s[^1] == '*' ? s[..^1] : s;

    // seqTypeFixReg: the characters that are not in the set of the type become fixchar
    private static string FixReg(string s, Def def, char fix)
    {
        var sb = new StringBuilder(s.Length);
        foreach (char c in s)
            sb.Append(c < 256 && def.Filter[c] ? c : fix);
        return sb.ToString();
    }

    // seqTypeTestI
    private static bool TestI(ref string seq, ref string typeOfSeq, Def def)
    {
        seq = def.Gaps ? GapConvert(seq) : RemoveGap(seq);
        if (def.Kind == 'P' && !IsProt(typeOfSeq, seq))
            return false;
        if (def.Kind == 'N' && !IsNuc(typeOfSeq, seq))
            return false;
        if (IsCharset(seq, def.Filter))
        {
            if (def.ConvertFrom.Length > 0)
                seq = ExchangeSet(seq, def.ConvertFrom, def.ConvertTo);
            return true;
        }

        return false;
    }

    // seqTypeFix
    private static bool Fix(ref string seq, ref string typeOfSeq, Def def)
    {
        if (!def.Gaps)
            seq = RemoveGap(seq);
        if (def.Name == "pureprotein")
            seq = StopTrim(seq);
        if (def.Ambig)
        {
            switch (def.Kind)
            {
                case 'P':
                    if (def.Name == "protein")
                        seq = StopTrim(seq);
                    seq = FixReg(seq, def, 'X');
                    break;
                case 'N':
                    seq = FixReg(seq, def, 'N');
                    break;
                default:
                    seq = FixReg(seq, def, IsNuc(typeOfSeq, seq) ? 'N' : 'X');
                    break;
            }
        }

        if (def.Name == "pureprotein")
            seq = StopTrim(seq);
        return TestI(ref seq, ref typeOfSeq, def);
    }

    /// <summary>
    /// <c>ajSeqTypeCheckIn</c> of a sequence (type <paramref name="typeOfSeq"/> <c>N</c> / <c>P</c> / empty, residues <paramref name="seq"/>) read by an input of the ACD type <paramref name="inputType"/> (empty: nothing given, only the gaps are converted).
    /// <paramref name="setType"/> is called with <c>N</c> or <c>P</c> when the check makes the sequence nucleic or protein (<c>ajSeqSetNuc</c> / <c>ajSeqSetProt</c>, which change the residues: x becomes n).
    /// </summary>
    /// <param name="error">The message printed on the error stream (<c>ajErr</c>).</param>
    /// <exception cref="EmbossFatalException">The type is unknown.</exception>
    public static bool CheckIn(EmbossSeq t, string inputType, Action<string> error)
    {
        string seq = t.SeqText;
        if (inputType.Length == 0)
        {
            t.Seq.Clear().Append(GapConvert(seq));
            return true;
        }

        var def = Find(inputType) ?? throw new EmbossFatalException($"Sequence type '{inputType}' unknown");
        seq = def.Gaps ? GapConvert(seq) : RemoveGap(seq);
        t.Seq.Clear().Append(seq);
        if (def.Kind == 'P')
        {
            if (IsProt(t.Type, seq))
                t.SetProt();
            else
            {
                error("Sequence is not a protein\n");
                return false;
            }
        }

        if (def.Kind == 'N')
        {
            if (IsNuc(t.Type, seq))
                t.SetNuc();
            else
            {
                error("Sequence is not nucleic\n");
                return false;
            }
        }

        seq = t.SeqText;
        if (IsCharset(seq, def.Filter))
        {
            if (def.ConvertFrom.Length > 0)
                seq = ExchangeSet(seq, def.ConvertFrom, def.ConvertTo);
            t.Seq.Clear().Append(seq);
            return true;
        }

        string type = t.Type;
        bool ok = Fix(ref seq, ref type, def);
        t.Seq.Clear().Append(seq);
        if (ok)
            return true;
        int bad = -1;
        for (int i = 0; i < seq.Length; i++)
        {
            if (seq[i] >= 256 || !def.Filter[seq[i]])
            {
                bad = i;
                break;
            }
        }

        if (bad >= 0)
        {
            error($"ajSeqTypeCheckIn: Sequence must be {def.Description}: found bad character '{seq[bad]}'");
            return false;
        }

        return true;
    }
}
