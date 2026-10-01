using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

/// <summary>
/// Heavy-tier wrapper equivalence for the Sequence-server complexity tools extended by the B04 completeness audit
/// (F35/F36/F53 complexity_mask_low linker / softMask / engine, F37 complexity_linguistic alphabetSize, F54
/// complexity_kmer_entropy correction / normalize): on fixed-seed random inputs every wrapper returns exactly the library
/// value (itself reference-locked in Seqeron.Genomics.Tests).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class B04AuditSequenceMcpEquivalenceTests
{
    private static string Random(Random rng, string alphabet, int length) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    [Test]
    public void ComplexityTools_EqualLibraryCalls()
    {
        var rng = new Random(5005);
        string[] corrections = ["none", "millerMadow", "grassberger"];
        for (int trial = 0; trial < 200; trial++)
        {
            string unit = Random(rng, "ACGT", rng.Next(1, 5));
            string s = Random(rng, "ACGT", rng.Next(1, 40)) + string.Concat(Enumerable.Repeat(unit, rng.Next(1, 30))) + Random(rng, "ACGT", rng.Next(0, 40));

            int k = rng.Next(1, 6);
            string correction = corrections[trial % 3];
            bool normalize = trial % 2 == 0;
            var mode = SequenceComplexity.ParseKmerEntropyCorrection(correction);
            Assert.That(SequenceTools.ComplexityKmerEntropy(s, k, correction, normalize).Entropy,
                Is.EqualTo(SequenceComplexity.CalculateKmerEntropy(s, k, mode, normalize)), s);

            int m = rng.Next(1, 12), a = rng.Next(4, 8);
            Assert.That(SequenceTools.ComplexityLinguistic(s, m, a).Complexity,
                Is.EqualTo(SequenceComplexity.CalculateLinguisticComplexity(s, m, a)), s);
            Assert.That(SequenceTools.ComplexityLinguistic(s, m).Complexity,
                Is.EqualTo(SequenceComplexity.CalculateLinguisticComplexity(s, m)), s);

            string iupac = s.Insert(rng.Next(0, s.Length), new string('N', rng.Next(1, 5)));
            int w = rng.Next(8, 65), level = rng.Next(2, 65), linker = rng.Next(1, 33);
            foreach (string engine in new[] { "sdust", "dustmasker" })
            {
                Assert.That(SequenceTools.ComplexityMaskLow(iupac, w, level / 10.0, 'N', linker, trial % 3 == 0, engine).MaskedSequence,
                    Is.EqualTo(SequenceComplexity.MaskLowComplexity(iupac, w, level / 10.0, 'N', linker, trial % 3 == 0,
                        SequenceComplexity.ParseDustEngine(engine))), iupac);
            }
        }
    }
}
