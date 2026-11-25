// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>
/// <c>ajFeatGetRemoteseq</c> of <c>ajfeat.c</c>: the sequence a feature location names by an identifier of another entry (<c>J00194.1:100..103</c>): the base of the USA of the sequence the features belong to
/// (<c>ajSeqUsaGetBase</c>: the format prefix and the file or database, without the range) plus <c>:</c> and the accession part of the identifier is read as a USA (<c>ajSeqGetFromUsa</c>, no sequence type), set to the range of the
/// feature (<c>ajSeqSetRange</c> / <c>ajSeqSetRangeRev</c>) and trimmed (<c>ajSeqTrim</c>).
/// </summary>
public sealed class EmbossRemoteSeq(IReadOnlyDictionary<string, string> files, IReadOnlyList<string>? definitions, EmbossFeatureReadContext context)
{
    /// <summary>The residues of the remote part <paramref name="gf"/> of a feature of <paramref name="seq"/>; empty when the entry cannot be read.</summary>
    public string Get(EmbossSeqFeature gf, EmbossSeq seq)
    {
        _names ??= definitions is { Count: > 0 } ? EmbossNames.Parse([.. definitions]) : EmbossNames.Parse();
        string remote = gf.Remote ?? "";
        string remoteId = EmbossFastaReader.SeqVersionAccession(remote) ?? remote;
        string usa = UsaGetBase(MakeUsa(seq), _names) + ":" + remoteId;
        var local = new EmbossFeatureReadContext { Program = context.Program };
        var q = EmbossSequenceQualifiers.None with { Nucleotide = !gf.Protein, Protein = gf.Protein };
        var all = EmbossSeqall.Open(context.Program, usa, files, q, false, local, "", true, definitions);
        foreach (string message in local.Messages)
        {
            if (!message.StartsWith("Error: Unable to read sequence '", StringComparison.Ordinal))
                context.Messages.Add(message);
        }

        if (all is null)
            return "";
        var rs = all.Peek();
        if (gf.Strand == '-')
            EmbossSeqRange.SetRangeRev(rs, gf.GetStart(), gf.GetEnd());
        else
            EmbossSeqRange.SetRange(rs, gf.GetStart(), gf.GetEnd());
        EmbossSeqRange.Trim(rs);
        return rs.Seq.ToString();
    }

    private EmbossNames? _names;

    // seqMakeUsa (ajSeqGetUsaS of a sequence with no USA): database-id, asis::TEXT or format::file[:entry], and the range [begin:end] or [-end:-begin:r]
    private static string MakeUsa(EmbossSeq t)
    {
        var sb = new System.Text.StringBuilder();
        if (t.Db.Length > 0)
            sb.Append(t.Db).Append("-id:").Append(t.Entryname);
        else if (t.Formatstr == "asis")
            sb.Append("asis::").Append(t.Reversed ? EmbossNuc.ReverseComplement(t.Seq.ToString()) : t.Seq.ToString());
        else if (t.Entryname.Length > 0)
            sb.Append(t.Formatstr).Append("::").Append(t.Filename).Append(':').Append(t.Entryname);
        else
            sb.Append(t.Formatstr).Append("::").Append(t.Filename);
        var range = new System.Text.StringBuilder("[");
        if (t.Reversed)
        {
            if (t.End != 0)
                range.Append((-t.End).ToString(System.Globalization.CultureInfo.InvariantCulture));
            range.Append(':');
            if (t.Begin != 0)
                range.Append((-t.Begin).ToString(System.Globalization.CultureInfo.InvariantCulture));
            range.Append(":r");
        }
        else
        {
            if (t.Begin != 0)
                range.Append(t.Begin.ToString(System.Globalization.CultureInfo.InvariantCulture));
            range.Append(':');
            if (t.End != 0)
                range.Append(t.End.ToString(System.Globalization.CultureInfo.InvariantCulture));
        }

        range.Append(']');
        if (range.Length > 3)
            sb.Append(range);
        return sb.ToString();
    }

    // ajSeqUsaGetBase
    private static string UsaGetBase(string usa, EmbossNames names)
    {
        string test = usa.Trim(' ', '\t', '\n');
        var range = Regex.Match(test, "^(.*)[[](-?[0-9]*):(-?[0-9]*)(:([Rr])?)?[]]$", RegexOptions.CultureInvariant | RegexOptions.Singleline);
        if (range.Success)
            test = "";                                          // ajRegPre: the text before the match, which starts at the start of the USA
        if (Regex.IsMatch(test, "^[Aa][Ss][Ii][Ss]:+(.+)$", RegexOptions.CultureInvariant) || Regex.IsMatch(test, "^(@|[Ll][Ii][Ss][Tt]:+)(.+)$", RegexOptions.CultureInvariant))
            return "";
        string result = "";
        var fmt = Regex.Match(test, "^([A-Za-z0-9-]*)::(.*)$", RegexOptions.CultureInvariant | RegexOptions.Singleline);
        if (fmt.Success)
        {
            result = fmt.Groups[1].Value + "::";
            test = fmt.Groups[2].Value;
        }

        var db = Regex.Match(test, "^([A-Za-z][A-Za-z0-9_]+)([-]([A-Za-z]+))?([:{]([^}]*)}?)?$", RegexOptions.CultureInvariant | RegexOptions.Singleline);
        if (db.Success)
        {
            string dbname = db.Groups[1].Value;
            if (names.IsDatabase(dbname) || names.AliasDatabase(ref dbname))              // ajNamDatabase
                return result + db.Groups[1].Value;
        }

        var id = Regex.Match(test, "^([^|]+[|]|[^:{%]+)(([:{%])(([^:}]+):)?([^:}]*)}?)?$", RegexOptions.CultureInvariant | RegexOptions.Singleline);
        if (id.Success)
            result += id.Groups[1].Value;
        return result;
    }
}
