// PROTMOTIF-HMM-001 — faithful killers for the Plan7 HMMER3/f header parsing of the GA gathering
// threshold and the STATS E-value calibration block. These are ordinary under-tested branches (not
// equivalent mutants): the canonical parser test only feeds the bundled profile, whose GA line has
// TWO tokens and whose STATS block is complete, so it never exercises the 1-token GA branch nor the
// "incomplete STATS -> not calibrated" rule. Source: HMMER User's Guide v3.4 (Eddy 2023) file format
// — "GA <ga1> [<ga2>];" and "All three STATS lines or none of them must be present".

using System;
using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

[TestFixture]
[Category("PROTMOTIF-HMM-001")]
public class Plan7ProfileHmm_ParseCalibration_Tests
{
    private static string ReadEmbedded(string fileName)
    {
        var asm = typeof(Plan7ProfileHmm).Assembly;
        string resourceName = $"Seqeron.Genomics.Analysis.Resources.{fileName}";
        using var stream = asm.GetManifestResourceStream(resourceName)
            ?? throw new InvalidOperationException($"Missing embedded resource: {resourceName}");
        using var reader = new System.IO.StreamReader(stream);
        return reader.ReadToEnd();
    }

    // The bundled GA line is "GA    22.9 22.9;" (two thresholds). A single-threshold "GA 22.9;" line is
    // equally valid per the format. The parser keeps GA1 whenever at least ONE token is present
    // (gaToks.Length >= 1); a mutant requiring >= 2 (`> 1`) would drop GA1 on a one-token line.
    [Test]
    public void Parse_SingleTokenGaLine_StillReadsGatheringThreshold()
    {
        var text = ReadEmbedded("PF00018_SH3_1.hmm").Replace("GA    22.9 22.9;", "GA    22.9;");
        var hmm = Plan7ProfileHmm.Parse(text);

        Assert.That(hmm.GatheringThreshold, Is.EqualTo(22.9).Within(1e-12),
            "a one-token GA line must still yield GA1 = 22.9");
    }

    // Calibration requires ALL THREE STATS lines. With the FORWARD line removed, the model is NOT
    // calibrated: Statistics must be null and Parse must not throw. A mutant that ORs the three
    // null-checks (instead of ANDs them) would either build a partial ScoreStatistics or dereference
    // the null FORWARD τ (fwdTau.Value) and throw.
    [Test]
    public void Parse_MissingForwardStatsLine_LeavesModelUncalibrated()
    {
        var text = ReadEmbedded("PF00018_SH3_1.hmm")
            .Replace("STATS LOCAL FORWARD   -4.5735  0.71923\n", string.Empty)
            .Replace("STATS LOCAL FORWARD   -4.5735  0.71923\r\n", string.Empty);

        Plan7ProfileHmm? hmm = null;
        Assert.DoesNotThrow(() => hmm = Plan7ProfileHmm.Parse(text),
            "an incomplete STATS block must parse without throwing");
        Assert.That(hmm!.Statistics, Is.Null,
            "MSV+VITERBI present but FORWARD absent → not calibrated → Statistics is null");
    }

    // Sanity counterpart: the unmodified profile DOES carry all three STATS lines → calibrated.
    [Test]
    public void Parse_AllThreeStatsLines_ModelIsCalibrated()
    {
        var hmm = Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));
        Assert.That(hmm.Statistics, Is.Not.Null, "all three STATS lines present → calibrated");
    }

    // Only STATS LOCAL lines calibrate the model — a GLOBAL distribution must be ignored (the parser
    // gates on sToks[1] == "LOCAL"). With VITERBI changed to GLOBAL, vitMu stays null → not calibrated.
    // A mutant that ORs the length guard ahead of the LOCAL check would accept the GLOBAL line and
    // (wrongly) complete the calibration triple.
    [Test]
    public void Parse_GlobalStatsLine_IsRejected_ModelUncalibrated()
    {
        var text = ReadEmbedded("PF00018_SH3_1.hmm")
            .Replace("STATS LOCAL VITERBI   -8.2932  0.71923", "STATS GLOBAL VITERBI  -8.2932  0.71923");
        var hmm = Plan7ProfileHmm.Parse(text);

        Assert.That(hmm.Statistics, Is.Null,
            "VITERBI given as GLOBAL must be ignored → only MSV+FORWARD (LOCAL) → not calibrated");
    }

    // A profile row with too few numeric fields is rejected with a message reporting how many were
    // actually found: Math.Max(0, toks.Length - skipFirst). Truncate the COMPO line (skipFirst=1,
    // count=20) to 5 numbers → "found only 5". A mutant that ADDs skipFirst would report 7.
    [Test]
    public void Parse_TruncatedCompoRow_ReportsCorrectFieldCount()
    {
        var lines = new List<string>(ReadEmbedded("PF00018_SH3_1.hmm").Split('\n'));
        int idx = lines.FindIndex(l => l.TrimStart().StartsWith("COMPO"));
        Assert.That(idx, Is.GreaterThanOrEqualTo(0), "PF00018 has a COMPO line");
        var toks = lines[idx].Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
        lines[idx] = string.Join("  ", toks[0], toks[1], toks[2], toks[3], toks[4], toks[5]); // COMPO + 5
        var text = string.Join("\n", lines);

        var ex = Assert.Throws<FormatException>(() => Plan7ProfileHmm.Parse(text));
        Assert.That(ex!.Message, Does.Contain("found only 5"),
            "5 numeric fields after the COMPO token must be reported as 'found only 5'");
    }
}
