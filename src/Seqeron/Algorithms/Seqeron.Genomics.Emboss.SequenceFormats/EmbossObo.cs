// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S1121, S3358, S3776, S107, S1854, S1066, S2325, S1172
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>
/// The queries of the drfind programs on the ontology <c>edam</c> and the data resource catalogue <c>drcat</c> (B09 WP67, F309): the databases of <c>emboss.standard</c> (method <c>emboss</c>: the indexes of the installation, read by
/// <see cref="EmbossBtreeCache"/> as <c>textAccessEmboss</c> reads them) answer <c>edam-nam:word</c>, <c>edam-isa:id</c>, <c>drcat-edat:id</c> ...: the texts of the entries in the order of the hits (<c>ajOboinQryS</c> / <c>ajResourceinQryS</c> + <c>ajOboinRead</c> / <c>ajResourceinRead</c>).
/// </summary>
public static class EmbossObo
{
    private static readonly Lazy<EmbossNames> StandardNames = new(() => EmbossNames.ParseWithStandard(null, null));

    private const string IndexDir = "/usr/share/EMBOSS/index/";
    public const string DataDir = "/usr/share/EMBOSS/data/";

    /// <summary>The files of the installation the queries read (the indexes, <c>EDAM.obo</c>, <c>DRCAT.dat</c>), as latin-1 text; null for any other name.</summary>
    internal static string? File(string path)
    {
        if (path.StartsWith(IndexDir, StringComparison.Ordinal) || path.StartsWith(DataDir, StringComparison.Ordinal))
        {
            string name = path[(path.LastIndexOf('/') + 1)..];
            try
            {
                return EmbossDrcat.ReadResource(name + ".gz");
            }
            catch (InvalidOperationException)
            {
                return null;
            }
        }

        return null;
    }

    /// <summary>The texts of the entries that the query <paramref name="query"/> (<c>edam-nam:protein</c>) finds, each from the offset of its entry to the end of the file, in the order of the hits; empty when there is none.</summary>
    /// <exception cref="NotSupportedException">A query with an operator or a character the library parses specially.</exception>
    internal static IReadOnlyList<string> Query(string query, string datatype, List<string>? messages = null)
    {
        var names = StandardNames.Value;
        var env = new EmbossQueryEnv { Files = File, Names = names, IsDatabase = names.IsDatabase, IsServer = names.IsServerName };
        var context = new EmbossFeatureReadContext();
        var textin = new EmbossTextin();
        textin.Query.DataType = datatype;
        string pqry = query;
        var r = EmbossQueryParse.ParseRead(ref pqry, textin, _ => true, env, context);
        messages?.AddRange(context.Messages);
        if (!r.Ok)
            return [];
        return r.HitTexts ?? (r.Text is null ? [] : [r.Text]);
    }
}
