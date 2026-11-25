// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

public static partial class EmbossSequencePrograms
{
    // the identifier of an obo term as ajOboinRead sets it: Fullid "EDAM_data:2201" gives the database "EDAM_data" (the text before the last colon) and the identifier "2201" (the text after it)
    private static (string Db, string Id) OboFullId(string entry)
    {
        int pos = 0;
        while (pos < entry.Length)
        {
            int nl = entry.IndexOf('\n', pos);
            string l = (nl < 0 ? entry[pos..] : entry[pos..nl]).TrimEnd('\r');
            pos = nl < 0 ? entry.Length : nl + 1;
            if (l.StartsWith("id:", StringComparison.Ordinal))
            {
                string v = l[3..].Trim();
                int bang = v.IndexOf('!');
                if (bang >= 0)
                    v = v[..bang].Trim();
                int colon = v.LastIndexOf(':');
                return colon < 0 ? ("", v) : (v[..colon], v[(colon + 1)..]);
            }
        }

        return ("", "");
    }

    // ajOboGetTree: the descendants of a term (the terms that name it in is_a, each followed by its own descendants), asked from the database of the term's own prefix (EDAM_data, EDAM_format ...)
    private static void OboTree(string db, string id, List<(string Db, string Id)> list)
    {
        foreach (string entry in EmbossObo.Query($"{db}-isa:{id}", "obo"))
        {
            var child = OboFullId(entry);
            list.Add(child);
            OboTree(child.Db, child.Id, list);
        }
    }

    /// <summary>
    /// <c>drfinddata</c> of EMBOSS 6.6.0 (B09 WP67, F309): the data resources of the catalogue (<c>DRCAT.dat</c>) whose EDAM data terms match the keywords: for every keyword of <c>-query</c> (comma separated) the terms of the ontology <c>edam</c> whose
    /// identifier, name or any word of the name or a synonym equals it (a <c>*</c> is a wildcard; the index words are lower case), with <c>-subclasses</c> all their descendants (<c>is_a</c>), and every catalogue entry that lists such a term as <c>EDAMdat</c>
    /// is written once, in the order the indexes return them, in the format of <c>-oformat</c> (<c>drcat</c> default, <c>basic</c>, <c>list</c>).
    /// </summary>
    internal static EmbossRunResult Drfinddata(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("drfinddata");
        string name = "";
        try
        {
            string query = a.String("query", "");
            string? outfile = a.Raw("outfile");
            bool subclasses = a.Bool("subclasses", true);
            _ = a.Bool("sensitive", false);
            string? oformat = a.RawAssoc("outfile", "oformat");
            string format = (oformat ?? "drcat").ToLowerInvariant();
            if (a.RawAssoc("outfile", "odirectory") is { Length: > 0 })
                throw new InvalidOperationException("drfinddata with an output directory: the reference program dies by a segmentation fault when the directory does not exist and the result is not defined for the port (the directory of the files is not known).");
            name = outfile ?? "drfinddata." + (oformat is null ? "drcat" : oformat.Length > 0 ? oformat : "drfinddata");
            if (outfile is { Length: 0 })
                throw a.Die("outfile", "Error: Output file is required");
            if (name != "stdout" && !EmbossOutPath.DirectoryExists(EmbossReportArgs.ExistingDirectories(files), name))
                throw new InvalidOperationException("drfinddata with an output file in a directory that does not exist: the reference program dies by a segmentation fault.");
            if (format is not ("drcat" or "basic" or "list"))
            {
                var errors = new List<string>();
                if (format.Length > 0)
                    errors.Add($"Error: Unknown output format '{oformat}'");
                errors.Add("Died: Output option -outfile: Format validation failed for type 'outresource'");
                var died = EmbossProgramBase.Died(ctx, new EmbossAcdDieException([.. errors]));
                return died with { Files = name == "stdout" ? new Dictionary<string, string>() : new Dictionary<string, string> { [name] = "" } };
            }

            var obotable = new HashSet<string>(StringComparer.Ordinal);
            var found = new HashSet<string>(StringComparer.Ordinal);
            var sb = new StringBuilder();
            var messages = new List<string>();
            foreach (string q in query.Split(','))
            {
                if (query.Length == 0)
                    break;
                foreach (string field in new[] { "id", "acc", "nam" })
                {
                    foreach (string entry in EmbossObo.Query($"edam-{field}:{q}", "obo", messages))
                    {
                        var first = OboFullId(entry);
                        var list = new List<(string Db, string Id)> { first };
                        if (subclasses)
                            OboTree(first.Db, first.Id, list);
                        foreach (var (_, term) in list)
                        {
                            if (!obotable.Add(term))
                                continue;
                            foreach (string rtext in EmbossObo.Query($"drcat-edat:{term}", "resource"))
                            {
                                var parsed = EmbossDrcat.Parse(rtext, 1);
                                if (parsed.Count == 0 || !found.Add(parsed[0].Id))
                                    continue;
                                var res = parsed[0];
                                sb.Append(format switch { "list" => $"drcat:{res.Id}\n", "basic" => BasicEntry(res), _ => res.Write() });
                            }
                        }
                    }
                }
            }

            string stderr = EmbossRunResult.Join(messages);
            if (name == "stdout")
                return new EmbossRunResult(sb.ToString(), stderr, 0, new Dictionary<string, string>());
            return new EmbossRunResult("", stderr, 0, new Dictionary<string, string> { [name] = sb.ToString() });
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e) with { Files = name.Length > 0 && name != "stdout" ? new Dictionary<string, string> { [name] = "" } : new Dictionary<string, string>() };
        }
    }

    // resourceoutWriteBasic
    private static string BasicEntry(EmbossResource r)
    {
        var sb = new StringBuilder();
        sb.Append("ID      ").Append(r.Id).Append('\n').Append("Name    ").Append(r.Name).Append('\n');
        if (r.Url.Length > 0)
            sb.Append("URL     ").Append(r.Url).Append('\n');
        return sb.Append('\n').ToString();
    }
}
