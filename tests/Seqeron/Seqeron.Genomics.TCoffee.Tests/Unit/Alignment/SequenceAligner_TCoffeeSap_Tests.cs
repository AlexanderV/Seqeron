// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using Seqeron.Genomics.TCoffee;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-12, B09 F298): the glue <c>sap_pair</c> of T-Coffee 13.41 (<c>util_dp_drivers.c</c> L2049) with the REAL SAP 1.1.3 installed.
/// Oracle: the 130 <c>sap_real</c> sets of TestData/Alignment/tcoffee_external_oracle.tsv, run by the user in the Docker package <c>handoff/A2-12-sap-dali</c>
/// (Ubuntu 24.04 linux/amd64, t_coffee 13.41.0.28bdc39 + SAP 1.1.3 built from the GitHub archive, sha256 of the binary in the returned ENV.txt): the rows t_coffee wrote together with the raw
/// <c>&lt;template1&gt;_&lt;template2&gt;.sap_results</c> files of its cache, i.e. SAP's own output for every pair of structures.
/// </summary>
[TestFixture]
public class SequenceAligner_TCoffeeSap_Tests
{
    private static List<SequenceAligner_TCoffeeExternal_Tests.Case> RealSapSets() => SequenceAligner_TCoffeeExternal_Tests.Oracle().Where(c => c.Feature == "sap_real").ToList();

    private static string?[] Pdb(SequenceAligner_TCoffeeExternal_Tests.Case c) =>
        c.Extra.GetProperty("pdb").EnumerateArray().Select(e => e.GetString() is { Length: > 0 } s ? s : null).ToArray();

    /// <summary>The raw SAP files of a set keyed by the pair of structure indices (the cache names are <c>p{i}.pdb_p{j}.pdb.sap_results</c>).</summary>
    private static Dictionary<(int, int), string> SapFiles(SequenceAligner_TCoffeeExternal_Tests.Case c)
    {
        var map = new Dictionary<(int, int), string>();
        if (!c.Extra.TryGetProperty("sap", out var sap)) return map;
        foreach (var f in sap.EnumerateObject())
        {
            var m = System.Text.RegularExpressions.Regex.Match(f.Name, @"^p(\d+)\.pdb_p(\d+)\.pdb\.sap_results$");
            Assert.That(m.Success, Is.True, f.Name);
            map[(int.Parse(m.Groups[1].Value), int.Parse(m.Groups[2].Value))] = f.Value.GetString()!;
        }

        return map;
    }

    private static TCoffeeProgramOptions Options(SequenceAligner_TCoffeeExternal_Tests.Case c, bool realSap)
    {
        var outputs = SequenceAligner_TCoffeeExternal_Tests.Outputs(c);
        if (realSap) outputs[TCoffeePairMethod.SapPair] = SapFiles(c);
        return TCoffeeProgramOptions.FromCommandLine(c.Flags) with { StructureTemplates = Pdb(c), StructureAlignerOutputs = outputs };
    }

    private static string Run(SequenceAligner_TCoffeeExternal_Tests.Case c, bool realSap) =>
        string.Join('|', TCoffeeAligner.MultipleAlignTCoffeeProgram(c.Seqs, Options(c, realSap)).AlignedSequences);

    /// <summary>Every row t_coffee wrote with the real SAP equals the port fed with SAP's raw files; where t_coffee stopped ("FAIL:1") the port raises <see cref="TCoffeeAbortedException"/>.</summary>
    [Test]
    public void RealSap_EveryRowEqualsTheProgram()
    {
        var sets = RealSapSets();
        Assert.That(sets, Has.Count.EqualTo(130));
        var failures = new List<string>();
        int aborted = 0;
        foreach (var c in sets)
        {
            if (c.Rows.StartsWith("FAIL:", StringComparison.Ordinal))
            {
                aborted++;
                Assert.Throws<TCoffeeAbortedException>(() => Run(c, true), $"t_coffee aborts on '{c.Flags}' {string.Join('|', c.Seqs)}");
                continue;
            }

            string got = Run(c, true);
            if (got != c.Rows) failures.Add($"'{c.Flags}' {string.Join('|', c.Seqs)}: got {got} expected {c.Rows}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {sets.Count} sets differ; first: {failures.FirstOrDefault()}");
        Assert.That(aborted, Is.GreaterThan(0), "the sets where the library stays empty and the program stops");
        Assert.That(SapFiles(sets[0]).Count + sets.Sum(s => SapFiles(s).Count), Is.GreaterThan(500), "SAP result files in the oracle (527)");
    }

    /// <summary>The glue matters: the TM-align replacement (SAP not installed) is NOT what the program does with SAP installed.</summary>
    [Test]
    public void RealSap_IsNotTheTMalignReplacement()
    {
        int differs = 0;
        foreach (var c in RealSapSets())
        {
            if (c.Rows.StartsWith("FAIL:", StringComparison.Ordinal)) continue;
            string fallback;
            try { fallback = Run(c, false); }
            catch (TCoffeeAbortedException) { differs++; continue; }
            if (fallback != c.Rows) differs++;
        }

        Assert.That(differs, Is.GreaterThan(10), "sets where the TM-align fallback differs from the real SAP");
    }

    private const string SapOutput =
        " file1 = 1469_1.139.sap_tmp\n prot1->Compound =  CHAIN: A;\n file2 = 1469_2.2513.sap_tmp\n prot2->Compound =  CHAIN: A;\nARNDCQEGHILKMFPSTWYVBZX\nmatrix constant = 8\nCycle 1, 11 residues selected\n\n score = 0.139664\n\n" +
        "Percent sel on aln =   43.75\nPercent aln in sel =   38.89\n" +
        " *L   11   0.0   1 A  \n *A   12   0.0   2 F**\n *Q   13   0.0   3 V* \n  F   14   0.0   4 M* \n" +
        "1469_2.2513.sap_tmp Weighted RMSd =   4.727 over 18 atoms 12.592868 0.834903   5.6 20.316572\n1469_2.2513.sap_tmp Un-weighted RMSd =   2.399 over best 7 atoms ( 0.834903    5.6 )\n";

    /// <summary>The first / last letter of every residue line after the second "Percent" token (the '*' marks removed); SAP's left column is the SECOND structure file.</summary>
    [Test]
    public void SapAlignedResidues_ReadsTheResidueLinesOfSapsOutput()
    {
        var r = TCoffeeAligner.TCoffeeSapAlignedResidues(SapOutput);
        Assert.That(r, Is.EqualTo(("LAQF", "AFVM")));
    }

    /// <summary>SAP's "-nan" result (no residue pair) is still a SAP file (two "Percent" tokens) with no pair; a text without them is "SAP failed to align".</summary>
    [Test]
    public void SapAlignedResidues_NanResultIsEmpty_ForeignTextIsNull()
    {
        const string nan = " score = 0.000000\n\nPercent sel on aln =    0.00\nPercent aln in sel =    -nan\n807_2.8344.sap_tmp Weighted RMSd =    -nan over 0 atoms 0.000000 0.000000  -nan -nan\n";
        Assert.That(TCoffeeAligner.TCoffeeSapAlignedResidues(nan), Is.EqualTo(("", "")));
        Assert.That(TCoffeeAligner.TCoffeeSapAlignedResidues("Segmentation fault"), Is.Null);
        Assert.That(TCoffeeAligner.TCoffeeSapAlignedResidues("Percent sel on aln = 1\n"), Is.Null);
        // a token that starts with '*' hides the rest of its line: the "Percent" after it does not count
        Assert.That(TCoffeeAligner.TCoffeeSapAlignedResidues("*x Percent Percent\nPercent\nPercent a\n A 1 0.0 1 B\n"), Is.EqualTo(("A", "B")));
        // only letters are kept
        Assert.That(TCoffeeAligner.TCoffeeSapAlignedResidues("Percent\nPercent\n 1 2 0.0 3 A\n A 2 0.0 3 -\n B 2 1.5 3 C\n"), Is.EqualTo(("B", "C")));
    }
}
