// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossRepeats, licensed under the terms stated in LICENSE (this directory).
using Seqeron.Genomics.Analysis;

namespace Seqeron.Mcp.EmbossRepeats.Tools;

public partial class EmbossRepeatsTools
{
    public static string[] TandemEntries(string kind, int minRepeat, int max, int threshold, bool mismatch, bool uniform)
    {
        var entries = new List<string>();
        if (kind != "equicktandem")
            entries.Add(Qualifier("minrepeat", minRepeat));
        entries.Add(Qualifier("maxrepeat", max));
        entries.Add(Qualifier("threshold", threshold));
        if (kind != "equicktandem" && mismatch)
            entries.Add("-mismatch");
        if (kind != "equicktandem" && uniform)
            entries.Add("-uniform");
        return [.. entries];
    }

    public static EmbossProgramReport RunTandemReport(string kind, EmbossReportRequest r, int minRepeat, int max, int threshold, bool mismatch, bool uniform, IReadOnlyList<EmbossTandemRepeat>? hits)
    {
        if (kind == "equicktandem")
            return EmbossRepeatReportPrograms.Equicktandem(r, max, threshold, hits);
        return EmbossRepeatReportPrograms.Etandem(r, minRepeat, max, threshold, mismatch, uniform, hits);
    }

    public static (string Which, IReadOnlyList<EmbossTandemRepeat> Hits) RunEmbossTandem(string sequence, string? program, int minRepeat, int? maxRepeat, int threshold, bool mismatch, bool uniform, CancellationToken cancellationToken)
    {
        var nucleotides = RequireEmbossIupacDna(sequence, nameof(sequence));
        string which = (program ?? "etandem").Trim().ToLowerInvariant();
        return which switch
        {
            "etandem" => (which, EmbossRepeatPrograms.FindTandemRepeatsEtandem(nucleotides, minRepeat, maxRepeat ?? minRepeat, threshold, mismatch, uniform, cancellationToken)),
            "equicktandem" => (which, EmbossRepeatPrograms.FindTandemRepeatsEquicktandem(nucleotides, maxRepeat ?? 600, threshold, cancellationToken)),
            _ => throw new ArgumentException("program must be etandem or equicktandem", nameof(program)),
        };
    }

    public static FindEmbossTandemRepeatsResult TandemResult(string which, IReadOnlyList<EmbossTandemRepeat> found)
    {
        bool detailed = which == "etandem";
        return new FindEmbossTandemRepeatsResult(which, [.. found.Select(h => new EmbossTandemItem(h.Position, h.Length, h.Score, h.RepeatSize, h.RepeatCount, detailed ? Math.Round((double)h.Identity, 1) : null, h.Consensus))]);
    }
}
