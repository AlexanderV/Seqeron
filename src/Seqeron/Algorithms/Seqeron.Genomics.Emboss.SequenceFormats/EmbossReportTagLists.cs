// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>The report tags (the ACD <c>taglist</c>) of the programs that write feature tables.</summary>
public static class EmbossReportTagLists
{

    /// <summary>The report tags (the ACD <c>taglist</c>) of the programs.</summary>
    public static IReadOnlyList<EmbossReportTag> TagList(string program) => program switch
    {
        "fuzznuc" or "fuzzpro" => [new("str", "pat", "Pattern"), new("int", "mismatch", "Mismatch")],
        "fuzztran" => [new("str", "pat", "Pattern"), new("int", "mismatch", "Mismatch"), new("int", "frame", "Frame"), new("int", "start", "PStart"), new("int", "end", "PEnd"), new("str", "translation", "Translation")],
        "etandem" => [new("int", "rpt_size", "Size"), new("int", "rpt_count", "Count"), new("float", "identity", "Identity"), new("str", "consensus", "Consensus")],
        "equicktandem" => [new("int", "rpt_size", "Size"), new("int", "rpt_count", "Count")],
        "featreport" => [],
        "twofeat" => [new("int", "starta", "starta"), new("int", "enda", "enda"), new("int", "startb", "startb"), new("int", "endb", "endb")],
        _ => throw new ArgumentException($"Unknown program '{program}'.", nameof(program)),
    };
}
