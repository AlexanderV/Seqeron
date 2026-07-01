using System;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class GeneratePrimerCandidatesTests
{
    // Same standard 258 bp template used by Seqeron.Genomics.Tests.
    private static string StandardTemplate()
    {
        var sb = new StringBuilder();
        while (sb.Length < 100) sb.Append("GAACTCGT");
        sb.Append(new string('T', 50));
        int revStart = sb.Length;
        while (sb.Length - revStart < 100) sb.Append("TCCGAAGT");
        return sb.ToString();
    }

    [Test]
    public void GeneratePrimerCandidates_Schema_ValidatesCorrectly()
    {
        var t = StandardTemplate();
        Assert.DoesNotThrow(() => MolToolsTools.generate_primer_candidates(t, 0, 60, true));
        Assert.Throws<ArgumentException>(() => MolToolsTools.generate_primer_candidates("", 0, 60));
        Assert.Throws<ArgumentException>(() => MolToolsTools.generate_primer_candidates(null!, 0, 60));
        Assert.Throws<ArgumentException>(() => MolToolsTools.generate_primer_candidates(t, -1, 60));
        Assert.Throws<ArgumentException>(() => MolToolsTools.generate_primer_candidates(t, 0, t.Length + 1));
        Assert.Throws<ArgumentException>(() => MolToolsTools.generate_primer_candidates(t, 60, 10));
    }

    [Test]
    public void GeneratePrimerCandidates_Binding_InvokesSuccessfully()
    {
        var t = StandardTemplate();

        // Default params: lengths 18..25. Region [0,60): start s with s+18<=60 => s in 0..42 (43 starts),
        // and for each s, len 18..min(25, 60-s). Enumerated count is 316 (verified against the algorithm).
        var result = MolToolsTools.generate_primer_candidates(t, 0, 60, true);
        var cands = result.Candidates;

        Assert.Multiple(() =>
        {
            Assert.That(cands.Count, Is.EqualTo(316));
            // Generation order: first candidate is start=0, len=18.
            Assert.That(cands[0].Position, Is.EqualTo(0));
            Assert.That(cands[0].Length, Is.EqualTo(18));
            Assert.That(cands[0].Sequence, Is.EqualTo(t.Substring(0, 18)));
            // All candidates within the default admissible length window.
            Assert.That(cands.All(c => c.Length >= 18 && c.Length <= 25), Is.True);
            // At start position 0, the eight lengths 18..25 are all emitted.
            Assert.That(cands.Where(c => c.Position == 0).Select(c => c.Length),
                Is.EqualTo(new[] { 18, 19, 20, 21, 22, 23, 24, 25 }));
        });
    }

    [Test]
    public void GeneratePrimerCandidates_Reverse_IsReverseComplementOfTemplate()
    {
        // 40 bp AACCGGTT-unit template; each reverse candidate sequence must equal the
        // reverse complement of the corresponding template substring.
        const string template = "AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT";
        var cands = MolToolsTools.generate_primer_candidates(template, 0, 40, forward: false).Candidates;

        Assert.That(cands.Count, Is.GreaterThan(0));
        foreach (var c in cands)
        {
            string sub = template.Substring(c.Position, c.Length);
            string expected = new DnaSequence(sub).ReverseComplement().Sequence;
            Assert.That(c.Sequence, Is.EqualTo(expected));
        }
    }
}
