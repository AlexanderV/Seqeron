// SEQ-AA-COMP-001 — residue counting routed through the canonical Core k-mer counter
// (review 2026-09, FIN-B03 final audit row D1). Evidence: docs/Validation/review-2026-09/B03.md,
// "Duplication sweep (finisher)".
// Oracle: a verbatim copy of the removed private SequenceStatistics.CountAsciiResidues loop
// (ASCII letters only, ASCII-only upper-case fold, Dictionary in first-occurrence order).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_ResidueCountingCanonical_Tests
{
    // Mixed-case ASCII letters, stop/gap/digits/space, U+017F 'ſ', U+0131 'ı', U+212A Kelvin sign,
    // Latin-1 / Greek / CJK letters and a few other BMP characters.
    private const string Alphabet =
        "ACDEFGHIKLMNPQRSTVWYBXZUOJ" + "acdefghiklmnpqrstvwybxzuoj" + "*-0123456789 ." +
        "ſıKéÉßαΑ中Ａａİ";

    // Copy of the removed private loop (SequenceStatistics.cs @ fe26648e, SS:243-256) — the oracle.
    private static Dictionary<char, int> OldCountAsciiResidues(string sequence)
    {
        var counts = new Dictionary<char, int>();
        foreach (char raw in sequence)
        {
            if (char.IsAsciiLetter(raw))
            {
                char ch = char.IsAsciiLetterLower(raw) ? (char)(raw - 32) : raw;
                counts[ch] = counts.GetValueOrDefault(ch) + 1;
            }
        }

        return counts;
    }

    [Test]
    public void Counts_RandomInputs_EqualOldLoopCountsAndKeyOrder()
    {
        var rng = new Random(20261014);
        for (int iter = 0; iter < 2000; iter++)
        {
            // Lengths span both CountKmersSpan buffers (stackalloc ≤ 256, pooled above).
            int length = rng.Next(1, iter % 10 == 0 ? 600 : 80);
            var chars = new char[length];
            for (int i = 0; i < length; i++)
                chars[i] = Alphabet[rng.Next(Alphabet.Length)];
            string sequence = new(chars);

            var expected = OldCountAsciiResidues(sequence).ToList();
            var actual = SequenceStatistics.CalculateAminoAcidComposition(sequence).Counts.ToList();
            Assert.That(actual, Is.EqualTo(expected), $"counts/key order differ for \"{sequence}\"");

            var fractionKeys = SequenceStatistics.CalculateAminoAcidFractions(sequence).Keys.ToList();
            Assert.That(fractionKeys, Is.EqualTo(expected.Select(kv => kv.Key).ToList()),
                $"fraction key order differs for \"{sequence}\"");
        }
    }

    [TestCase("mKſkıK*-1", "M:1,K:3")]
    [TestCase("ſıK", "")]
    [TestCase("wAw", "W:2,A:1")]
    public void Counts_NonAsciiLettersNotFolded_FirstOccurrenceOrder(string sequence, string expected)
    {
        var counts = SequenceStatistics.CalculateAminoAcidComposition(sequence).Counts;
        Assert.That(string.Join(",", counts.Select(kv => $"{kv.Key}:{kv.Value}")), Is.EqualTo(expected));
    }
}
