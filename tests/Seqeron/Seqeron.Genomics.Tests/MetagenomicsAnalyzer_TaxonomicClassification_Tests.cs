using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Test suite for META-CLASS-001: Taxonomic Classification.
/// Covers MetagenomicsAnalyzer.ClassifyReads and BuildKmerDatabase methods.
/// 
/// Evidence sources:
/// - Wikipedia (Metagenomics): K-mer based classification principles
/// - Kraken (CCB JHU): Canonical k-mer algorithm, k=31 default
/// - Wood &amp; Salzberg (2014): Kraken paper
/// </summary>
[TestFixture]
public class MetagenomicsAnalyzer_TaxonomicClassification_Tests
{
    #region BuildKmerDatabase Tests

    [Test]
    [Description("M1: Empty input returns empty database")]
    public void BuildKmerDatabase_EmptyInput_ReturnsEmptyDatabase()
    {
        var references = new List<(string TaxonId, string Sequence)>();

        var database = MetagenomicsAnalyzer.BuildKmerDatabase(references);

        Assert.That(database, Is.Empty);
    }

    [Test]
    [Description("M2: Sequence shorter than k produces no k-mers")]
    public void BuildKmerDatabase_SequenceShorterThanK_ReturnsEmpty()
    {
        var references = new List<(string TaxonId, string Sequence)>
        {
            ("Bacteria|Test", "ATGC") // 4 bp, shorter than k=31
        };

        var database = MetagenomicsAnalyzer.BuildKmerDatabase(references, k: 31);

        Assert.That(database, Is.Empty);
    }

    [Test]
    [Description("M3: Valid reference produces k-mers in database")]
    public void BuildKmerDatabase_ValidReference_ProducesKmers()
    {
        // k=5, "ACGTACTGAC" (10bp) → 6 k-mers, all unique after canonicalization:
        //   i=0: ACGTA (canon: ACGTA, A<T)
        //   i=1: CGTAC (canon: CGTAC, C<G)
        //   i=2: GTACT (canon: AGTAC, A<G)
        //   i=3: TACTG (canon: CAGTA, C<T)
        //   i=4: ACTGA (canon: ACTGA, A<T)
        //   i=5: CTGAC (canon: CTGAC, C<G)
        var references = new List<(string TaxonId, string Sequence)>
        {
            ("Bacteria|Escherichia", "ACGTACTGAC")
        };

        var database = MetagenomicsAnalyzer.BuildKmerDatabase(references, k: 5);

        Assert.Multiple(() =>
        {
            Assert.That(database.Count, Is.EqualTo(6), "10-5+1 = 6 unique canonical k-mers");
            Assert.That(database.Values.All(v => v == "Bacteria|Escherichia"), Is.True,
                "All k-mers should map to the reference taxon");
        });
    }

    [Test]
    [Description("M4: Uses canonical k-mers (lexicographically smaller of kmer or reverse complement)")]
    public void BuildKmerDatabase_UsesCanonicalKmers()
    {
        // "AAAAAAAAAA" is its own reverse complement as "TTTTTTTTTT"
        // Canonical is "AAAAAAAAAA" (lexicographically smaller)
        var references = new List<(string TaxonId, string Sequence)>
        {
            ("Taxon1", "AAAAAAAAAAAAAAAAAAAA") // 20 A's
        };

        var database = MetagenomicsAnalyzer.BuildKmerDatabase(references, k: 10);

        Assert.Multiple(() =>
        {
            Assert.That(database.ContainsKey("AAAAAAAAAA"), Is.True, "Should contain canonical k-mer AAAAAAAAAA");
            Assert.That(database.ContainsKey("TTTTTTTTTT"), Is.False, "Should not contain non-canonical TTTTTTTTTT");
        });
    }

    [Test]
    [Description("M4b: Canonical k-mer uses reverse complement when it is lexicographically smaller")]
    public void BuildKmerDatabase_CanonicalKmer_UsesReverseComplementWhenSmaller()
    {
        // Sequence starting with T's - reverse complement starts with A's
        // T's RC = A's, so A's should be canonical
        var references = new List<(string TaxonId, string Sequence)>
        {
            ("Taxon1", "TTTTTTTTTTTTTTTTTTTT") // 20 T's
        };

        var database = MetagenomicsAnalyzer.BuildKmerDatabase(references, k: 10);

        Assert.Multiple(() =>
        {
            Assert.That(database.ContainsKey("AAAAAAAAAA"), Is.True, "Should contain canonical k-mer AAAAAAAAAA (RC of TTTTTTTTTT)");
            Assert.That(database.ContainsKey("TTTTTTTTTT"), Is.False, "Should not contain non-canonical TTTTTTTTTT");
        });
    }

    [Test]
    [Description("M5: Non-ACGT characters excluded from k-mers")]
    public void BuildKmerDatabase_NonAcgtCharacters_Excluded()
    {
        // Sequence with N in the middle - k-mers containing N should be excluded
        var references = new List<(string TaxonId, string Sequence)>
        {
            ("Taxon1", "ATGCATGCATGCNNNNATGCATGCATGC") // Contains N's
        };

        var database = MetagenomicsAnalyzer.BuildKmerDatabase(references, k: 10);

        // All k-mers in database should only contain ACGT
        foreach (var kmer in database.Keys)
        {
            Assert.That(kmer, Does.Match("^[ACGT]+$"), $"K-mer '{kmer}' contains non-ACGT characters");
        }
    }

    [Test]
    [Description("M6: Database k-mer count follows formula for valid sequences")]
    public void BuildKmerDatabase_KmerCount_FollowsFormula()
    {
        // Non-repeating sequence where all canonical k-mers are unique:
        // k=5, "ACGTACTGAC" (10bp) → 6 raw k-mers, all distinct after canonicalization
        const int k = 5;
        const string sequence = "ACGTACTGAC";
        int expectedKmers = sequence.Length - k + 1; // 6

        var references = new List<(string TaxonId, string Sequence)>
        {
            ("Taxon1", sequence)
        };

        var database = MetagenomicsAnalyzer.BuildKmerDatabase(references, k: k);

        Assert.Multiple(() =>
        {
            Assert.That(database.Count, Is.EqualTo(expectedKmers),
                "For non-repeating sequences with unique canonical k-mers: count = len-k+1 = 6");
            foreach (var kmer in database.Keys)
            {
                Assert.That(kmer.Length, Is.EqualTo(k), $"K-mer '{kmer}' should have length k={k}");
            }
        });
    }

    [Test]
    [Description("S1: Mixed case input handled correctly")]
    public void BuildKmerDatabase_MixedCase_NormalizedToUppercase()
    {
        var references = new List<(string TaxonId, string Sequence)>
        {
            ("Taxon1", "atgcatgcatgcatgcatgc") // lowercase
        };

        var database = MetagenomicsAnalyzer.BuildKmerDatabase(references, k: 10);

        Assert.That(database.Count, Is.GreaterThan(0), "Should process lowercase sequences");

        // All keys should be uppercase
        foreach (var kmer in database.Keys)
        {
            Assert.That(kmer, Is.EqualTo(kmer.ToUpperInvariant()), "K-mer should be uppercase");
        }
    }

    #endregion

    #region ClassifyReads Tests

    [Test]
    [Description("M7: Empty sequence returns Unclassified with Confidence=0")]
    public void ClassifyReads_EmptySequence_ReturnsUnclassified()
    {
        var database = new Dictionary<string, string>();
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "")
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].ReadId, Is.EqualTo("read1"));
            Assert.That(results[0].Kingdom, Is.EqualTo("Unclassified"));
            Assert.That(results[0].Confidence, Is.EqualTo(0));
            Assert.That(results[0].TotalKmers, Is.EqualTo(0));
        });
    }

    [Test]
    [Description("M8: Sequence shorter than k returns Unclassified")]
    public void ClassifyReads_SequenceShorterThanK_ReturnsUnclassified()
    {
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Bacteria|Test" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGC") // 4 bp, shorter than k=14
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].Kingdom, Is.EqualTo("Unclassified"));
            Assert.That(results[0].TotalKmers, Is.EqualTo(0));
        });
    }

    [Test]
    [Description("M9: No database matches returns Unclassified")]
    public void ClassifyReads_NoMatch_ReturnsUnclassified()
    {
        var database = new Dictionary<string, string>
        {
            { "CCCCCCCCCCCCCC", "Bacteria|Other" } // canonical form of GGGG... (C < G)
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCGATCGATCGATCGATCG") // No matching k-mers
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].Kingdom, Is.EqualTo("Unclassified"));
            Assert.That(results[0].MatchedKmers, Is.EqualTo(0));
            Assert.That(results[0].TotalKmers, Is.GreaterThan(0));
        });
    }

    [Test]
    [Description("M10: Matching k-mers classify to correct taxon")]
    public void ClassifyReads_MatchingKmers_ClassifiesCorrectly()
    {
        var database = new Dictionary<string, string>
        {
            { "ATGCGATCGATCGA", "Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia|coli" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCGATCGATCGATCGATCG")
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        // 21bp, k=14 → 8 k-mers; only i=0 (ATGCGATCGATCGA, palindromic) matches the DB entry
        Assert.Multiple(() =>
        {
            Assert.That(results[0].Kingdom, Is.EqualTo("Bacteria"));
            Assert.That(results[0].Phylum, Is.EqualTo("Proteobacteria"));
            Assert.That(results[0].Genus, Is.EqualTo("Escherichia"));
            Assert.That(results[0].Species, Is.EqualTo("coli"));
            Assert.That(results[0].MatchedKmers, Is.EqualTo(1), "Only the first k-mer matches the DB entry");
            Assert.That(results[0].TotalKmers, Is.EqualTo(8), "21-14+1 = 8 non-ambiguous k-mers");
            Assert.That(results[0].Confidence, Is.EqualTo(0.125).Within(0.001), "Confidence = 1/8 = 0.125");
        });
    }

    [Test]
    [Description("M11: Output count equals input read count")]
    public void ClassifyReads_OutputCountEqualsInputCount()
    {
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Bacteria|Test|Genus|species" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCATGCATGCATGC"),
            ("read2", "ATGCATGCATGCATGCATGC"),
            ("read3", "GGGGGGGGGGGGGGGGGGGG"),
            ("read4", ""),
            ("read5", "AT")
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        Assert.That(results.Count, Is.EqualTo(reads.Count), "Output count must equal input count");
    }

    [Test]
    [Description("M12: Confidence = MatchedKmers / TotalKmers")]
    public void ClassifyReads_ConfidenceCalculation_IsCorrect()
    {
        // Create a read where we know exactly how many k-mers will match
        const int k = 10;
        const string kmer = "ATGCATGCAT"; // 10 bp k-mer
        var database = new Dictionary<string, string>
        {
            { kmer, "Bacteria|Test" }
        };

        // Read that contains exactly this k-mer once: ATGCATGCAT + extra
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCATG") // 11 bp → 2 k-mers, first one matches
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: k).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].TotalKmers, Is.EqualTo(2), "Should have 2 k-mers (11-10+1)");
            Assert.That(results[0].MatchedKmers, Is.EqualTo(1), "Should have 1 matching k-mer");
            Assert.That(results[0].Confidence, Is.EqualTo(0.5).Within(0.01), "Confidence should be 1/2 = 0.5");
        });
    }

    [Test]
    [Description("M13: TotalKmers = non-ambiguous k-mers count; equals max(0, len - k + 1) for all-ACGT sequences")]
    public void ClassifyReads_TotalKmers_MatchesFormula()
    {
        const int k = 10;
        var database = new Dictionary<string, string>();
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCATGCATGCATGC"), // 20 bp → 11 k-mers
            ("read2", "ATGCATGCAT"),            // 10 bp → 1 k-mer
            ("read3", "ATGCATGCA"),             // 9 bp → 0 k-mers (shorter than k)
            ("read4", "")                        // 0 bp → 0 k-mers
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: k).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].TotalKmers, Is.EqualTo(11), "20 - 10 + 1 = 11");
            Assert.That(results[1].TotalKmers, Is.EqualTo(1), "10 - 10 + 1 = 1");
            Assert.That(results[2].TotalKmers, Is.EqualTo(0), "9 < 10, so 0");
            Assert.That(results[3].TotalKmers, Is.EqualTo(0), "Empty sequence");
        });
    }

    [Test]
    [Description("M14: MatchedKmers ≤ TotalKmers")]
    public void ClassifyReads_MatchedKmers_BoundedByTotal()
    {
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Bacteria|Test" },
            { "CATGCATGCATGCA", "Bacteria|Test" } // canonical of TGCATGCATGCATG
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCATGCATGCATGCATGCATGCATGC")
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        Assert.That(results[0].MatchedKmers, Is.LessThanOrEqualTo(results[0].TotalKmers),
            "MatchedKmers must be ≤ TotalKmers");
    }

    [Test]
    [Description("M15: Multiple reads all classified")]
    public void ClassifyReads_MultipleReads_AllClassified()
    {
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Bacteria|Genus1|species1" },
            { "GCTAGCTAGCTAGC", "Bacteria|Genus2|species2" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCATGCATGCATGC"),
            ("read2", "GCTAGCTAGCTAGCTAGCTAG"),
            ("read3", "AAAAAAAAAAAAAAAAAAAA") // No match
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results.Count, Is.EqualTo(3));
            Assert.That(results[0].Kingdom, Is.EqualTo("Bacteria"), "Read1 matches Genus1 DB entry");
            Assert.That(results[1].Kingdom, Is.EqualTo("Bacteria"), "Read2 matches Genus2 DB entry");
            Assert.That(results[2].Kingdom, Is.EqualTo("Unclassified"), "Read3 has no DB matches");
        });
    }

    [Test]
    [Description("M16: Taxonomy string parsed correctly (pipe-delimited)")]
    public void ClassifyReads_TaxonomyParsing_PipeDelimited()
    {
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Kingdom|Phylum|Class|Order|Family|Genus|Species" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCATGCATGCATGC")
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].Kingdom, Is.EqualTo("Kingdom"));
            Assert.That(results[0].Phylum, Is.EqualTo("Phylum"));
            Assert.That(results[0].Class, Is.EqualTo("Class"));
            Assert.That(results[0].Order, Is.EqualTo("Order"));
            Assert.That(results[0].Family, Is.EqualTo("Family"));
            Assert.That(results[0].Genus, Is.EqualTo("Genus"));
            Assert.That(results[0].Species, Is.EqualTo("Species"));
        });
    }

    [Test]
    [Description("M16b: Taxonomy string parsed correctly (semicolon-delimited)")]
    public void ClassifyReads_TaxonomyParsing_SemicolonDelimited()
    {
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Kingdom;Phylum;Class;Order;Family;Genus;Species" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCATGCATGCATGC")
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].Kingdom, Is.EqualTo("Kingdom"));
            Assert.That(results[0].Phylum, Is.EqualTo("Phylum"));
            Assert.That(results[0].Genus, Is.EqualTo("Genus"));
            Assert.That(results[0].Species, Is.EqualTo("Species"));
        });
    }

    [Test]
    [Description("S1: Mixed case input handled")]
    public void ClassifyReads_MixedCaseInput_Handled()
    {
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Bacteria|Test" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "atgcatgcatgcatgcatgc") // lowercase
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        // Uppercased read = "ATGCATGCATGCATGCATGC" (20bp, k=14, 7 k-mers)
        // i=0,4: ATGCATGCATGCAT (palindromic) → match; rest don't match
        Assert.Multiple(() =>
        {
            Assert.That(results[0].Kingdom, Is.EqualTo("Bacteria"), "Should classify correctly despite case");
            Assert.That(results[0].MatchedKmers, Is.EqualTo(2), "Should match same k-mers as uppercase input");
        });
    }

    [Test]
    [Description("S2: Multiple taxon matches resolves to highest count")]
    public void ClassifyReads_MultipleTaxonMatches_ResolvesToHighestCount()
    {
        // First k-mer matches Taxon1, next two match Taxon2
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Bacteria|Taxon1" },
            { "CATGCATGCATGCA", "Bacteria|Taxon2" }, // canonical of TGCATGCATGCATG
            { "GCATGCATGCATGC", "Bacteria|Taxon2" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCATGCATGC") // 16 bp → 3 k-mers
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        // Should classify to Taxon2 (2 hits) rather than Taxon1 (1 hit)
        Assert.That(results[0].Phylum, Is.EqualTo("Taxon2"), "Should resolve to taxon with most hits");
    }

    [Test]
    [Description("S3: ClassifyReads uses canonical k-mers for database lookup (non-palindromic k-mer)")]
    public void ClassifyReads_CanonicalKmerLookup_MatchesReverseComplement()
    {
        // DB has canonical form AAAAAAAAAA (A < T, so forward is canonical).
        // Read contains TTTTTTTTTT which canonicalizes to AAAAAAAAAA → should match.
        const int k = 10;
        var database = new Dictionary<string, string>
        {
            { "AAAAAAAAAA", "Bacteria|CanonTest" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "TTTTTTTTTTT") // 11 bp → 2 k-mers, both TTTTTTTTTT → canon AAAAAAAAAA
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: k).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].Kingdom, Is.EqualTo("Bacteria"), "RC k-mer should match canonical DB entry");
            Assert.That(results[0].MatchedKmers, Is.EqualTo(2), "Both k-mers should match via canonicalization");
            Assert.That(results[0].TotalKmers, Is.EqualTo(2));
            Assert.That(results[0].Confidence, Is.EqualTo(1.0).Within(0.001));
        });
    }

    [Test]
    [Description("S4: Ambiguous nucleotides excluded from TotalKmers (per Kraken: Q = non-ambiguous k-mers only)")]
    public void ClassifyReads_AmbiguousNucleotides_ExcludedFromTotalKmers()
    {
        // Read: ATGCATGCA N ATGCATGCAT (20 bp, N at position 9)
        // k=10 → 11 raw k-mers, but 10 of them span the N → only 1 valid k-mer
        // Valid k-mer at i=10: ATGCATGCAT (palindromic)
        const int k = 10;
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCAT", "Bacteria|AmbigTest" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCANATGCATGCAT") // N at position 9
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: k).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].TotalKmers, Is.EqualTo(1),
                "Only 1 of 11 k-mers is non-ambiguous (per Kraken: Q excludes k-mers with non-ACGT)");
            Assert.That(results[0].MatchedKmers, Is.EqualTo(1));
            Assert.That(results[0].Confidence, Is.EqualTo(1.0).Within(0.001));
        });
    }

    [Test]
    [Description("S5: Confidence C/Q where C = k-mers supporting winning taxon only (per Kraken scoring)")]
    public void ClassifyReads_MultiTaxon_ConfidenceUsesWinningTaxonCount()
    {
        // Read: ATGCATGCATGCATGCA (17 bp, k=14 → 4 k-mers)
        // Canonical k-mers:
        //   i=0: ATGCATGCATGCAT (palindromic) → Taxon1
        //   i=1: TGCATGCATGCATG → canon CATGCATGCATGCA → Taxon1
        //   i=2: GCATGCATGCATGC (palindromic) → Taxon2
        //   i=3: CATGCATGCATGCA → Taxon1
        // Taxon1: 3 hits, Taxon2: 1 hit → winner = Taxon1
        // Per Kraken: C = 3 (winning taxon), Q = 4 → Confidence = 0.75
        const int k = 14;
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Bacteria|Taxon1" },
            { "CATGCATGCATGCA", "Bacteria|Taxon1" },
            { "GCATGCATGCATGC", "Bacteria|Taxon2" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCATGCATGCA") // 17 bp → 4 k-mers
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: k).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].Phylum, Is.EqualTo("Taxon1"), "Should classify to taxon with most hits");
            Assert.That(results[0].MatchedKmers, Is.EqualTo(3),
                "C = k-mers supporting winning taxon only (not all 4 matched k-mers)");
            Assert.That(results[0].TotalKmers, Is.EqualTo(4), "Q = all non-ambiguous k-mers");
            Assert.That(results[0].Confidence, Is.EqualTo(0.75).Within(0.001),
                "Confidence = C/Q = 3/4 = 0.75 per Kraken formula");
        });
    }

    #endregion

    #region Invariant Tests

    [Test]
    [Description("Invariant: All outputs have valid confidence range [0, 1]")]
    public void ClassifyReads_AllOutputs_HaveValidConfidence()
    {
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Bacteria|Test" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "ATGCATGCATGCATGCATGC"),
            ("read2", "GGGGGGGGGGGGGGGGGGGG"),
            ("read3", ""),
            ("read4", "AT")
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        foreach (var result in results)
        {
            Assert.That(result.Confidence, Is.InRange(0.0, 1.0),
                $"Confidence for {result.ReadId} should be in [0, 1]");
        }
    }

    [Test]
    [Description("Invariant: Unclassified reads have MatchedKmers = 0")]
    public void ClassifyReads_UnclassifiedReads_HaveZeroMatchedKmers()
    {
        var database = new Dictionary<string, string>
        {
            { "ATGCATGCATGCAT", "Bacteria|Test" }
        };
        var reads = new List<(string Id, string Sequence)>
        {
            ("read1", "GGGGGGGGGGGGGGGGGGGG") // No matches
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].Kingdom, Is.EqualTo("Unclassified"));
            Assert.That(results[0].MatchedKmers, Is.EqualTo(0));
        });
    }

    [Test]
    [Description("Invariant: ReadId preserved in output")]
    public void ClassifyReads_ReadIdPreserved()
    {
        var database = new Dictionary<string, string>();
        var reads = new List<(string Id, string Sequence)>
        {
            ("unique_id_123", "ATGCATGCATGCATGCATGC"),
            ("another-id", "GCTAGCTAGCTAGCTAGCTAG")
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, database, k: 14).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results[0].ReadId, Is.EqualTo("unique_id_123"));
            Assert.That(results[1].ReadId, Is.EqualTo("another-id"));
        });
    }

    #endregion
}
