// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>B09 A5-7: the <c>sequenceQualifiers</c> parameter of fuzznuc / fuzzpro / fuzztran / etandem / equicktandem / seqmatchall / distmat. Values are what the installed EMBOSS 6.6.0 prints.</summary>
[TestFixture]
public class EmbossSequenceQualifiersToolsTests
{

    [Test]
    public void Megamerger_MergesTheTrimmedSequences()
    {
        // megamerger a.fa b.fa -wordsize 6 -sbegin1 37 (installed 6.6.0, oracle row of emboss_seqopts_oracle.tsv)
        const string a = "CAGACCGATGCGCGCCGTCTTCGCCGCTTAGCTTGACGACAGTACCCGGCGGGAATGCTGGCTGGCT";
        const string b = "AGACCGATGCGCGCCGTCTTCGCCGCTTAGCTTGACGACAGTACCCGGCGGGAATGCTGGCTGGCTACACCAGTCGGTCTACAACTTGCTGCGAA";
        var r = EmbossPairwiseAlignmentTools.Megamerger(a, b, wordSize: 6, name1: "a", name2: "b", sequenceQualifiers: "-sbegin1 37");
        Assert.That(r.MergedSequence, Is.EqualTo("agaccgatgcgcgccgtcttcgccgcttagcttgacgacagtacccggcgggaatgctggctggctacaccagtcggtctacaacttgctgcgaa"));
        Assert.That(r.Report, Does.Contain("Using b 1-35 as the initial sequence"));
        Assert.That(r.Report, Does.Contain("Matching region a 1-31 : b 36-66"));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.Megamerger(a, b, wordSize: 6, sequenceQualifiers: "-sprotein1"));
    }}
