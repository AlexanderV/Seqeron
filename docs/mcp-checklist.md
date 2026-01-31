# MCP Implementation Checklist v4

## 0. Definition of Done

### Tool DoD
- [ ] MethodId linked
- [ ] toolName/serverName frozen
- [ ] inputSchema defined
- [ ] outputSchema defined
- [ ] Error mapping defined
- [ ] 2 minimal tests (Schema + Binding)
- [ ] `{tool}.mcp.json` created
- [ ] `{tool}.md` created
- [ ] Traceability row added

### Server DoD
- [ ] All tools meet Tool DoD
- [ ] Program.cs configured
- [ ] Server builds
- [ ] Integration test passes
- [ ] README.md created

### Delivery DoD
- [ ] All 12 servers meet Server DoD
- [ ] G1-G5 gates pass
- [ ] Cross-server tests pass
- [ ] Release notes written

---

## 1. Preparation

- [ ] 1.1 Ensure MCP projects are included in `Seqeron.sln`
- [ ] 1.2 Create `SuffixTree.Mcp.Shared` project
- [ ] 1.3 Add `SequenceValidator.cs` helper
- [ ] 1.4 Add `OutputLimiter.cs` helper
- [ ] 1.5 Add `ResponseExtensions.cs`
- [ ] 1.6 Define error code catalog (1000-5999)
- [ ] 1.7 Create `.mcp.json` template
- [ ] 1.8 Create `.md` template
- [ ] 1.9 Set up test project `SuffixTree.Mcp.Tests`

---

## 2. Inventory & Traceability

- [ ] 2.1 Create traceability matrix spreadsheet
- [ ] 2.2 Populate 241 tool rows
- [ ] 2.3 Verify HasDocs status for each tool
- [ ] 2.4 Link DocRef for documented tools
- [ ] 2.5 Mark TBD items for remediation

---

## 3. Per-Server Checklist

### 3.1 Server: Core (12 tools)

- [ ] 3.1.1 Create project `SuffixTree.Mcp.Core`
- [ ] 3.1.2 Add `Program.cs`
- [ ] 3.1.3 Add `SuffixTree.Mcp.Core.csproj`
- [ ] 3.1.4 Create `Tools/SuffixTreeTools.cs`
- [ ] 3.1.5 Create `Tools/ApproximateMatchTools.cs`
- [ ] 3.1.6 Implement all 12 tools
- [ ] 3.1.7 Write 24 tests (12 tools × 2)
- [ ] 3.1.8 Create 12 `.mcp.json` files
- [ ] 3.1.9 Create 12 `.md` files
- [ ] 3.1.10 Build and verify
- [ ] 3.1.11 Create README.md

### 3.2 Server: Sequence (35 tools)

- [ ] 3.2.1 Create project `Seqeron.Mcp.Sequence`
- [ ] 3.2.2 Add `Program.cs`
- [ ] 3.2.3 Create `Tools/DnaTools.cs`
- [ ] 3.2.4 Create `Tools/RnaTools.cs`
- [ ] 3.2.5 Create `Tools/ProteinTools.cs`
- [ ] 3.2.6 Create `Tools/StatisticsTools.cs`
- [ ] 3.2.7 Create `Tools/ComplexityTools.cs`
- [ ] 3.2.8 Create `Tools/KmerTools.cs`
- [ ] 3.2.9 Create `Tools/TranslationTools.cs`
- [ ] 3.2.10 Implement all 35 tools
- [ ] 3.2.11 Write 70 tests (35 tools × 2)
- [ ] 3.2.12 Create 35 `.mcp.json` files
- [ ] 3.2.13 Create 35 `.md` files
- [ ] 3.2.14 Build and verify
- [ ] 3.2.15 Create README.md

### 3.3 Server: Parsers (45 tools)

- [ ] 3.3.1 Create project `Seqeron.Mcp.Parsers`
- [ ] 3.3.2 Add `Program.cs`
- [ ] 3.3.3 Create `Tools/FastaTools.cs`
- [ ] 3.3.4 Create `Tools/FastqTools.cs`
- [ ] 3.3.5 Create `Tools/GenBankTools.cs`
- [ ] 3.3.6 Create `Tools/EmblTools.cs`
- [ ] 3.3.7 Create `Tools/GffTools.cs`
- [ ] 3.3.8 Create `Tools/BedTools.cs`
- [ ] 3.3.9 Create `Tools/VcfTools.cs`
- [ ] 3.3.10 Create `Tools/ConversionTools.cs`
- [ ] 3.3.11 Create `Tools/QualityTools.cs`
- [ ] 3.3.12 Implement all 45 tools
- [ ] 3.3.13 Write 90 tests (45 tools × 2)
- [ ] 3.3.14 Create 45 `.mcp.json` files
- [ ] 3.3.15 Create 45 `.md` files
- [ ] 3.3.16 Build and verify
- [ ] 3.3.17 Create README.md

### 3.4 Server: Alignment (15 tools)

- [ ] 3.4.1 Create project `SuffixTree.Mcp.Alignment`
- [ ] 3.4.2 Add `Program.cs`
- [ ] 3.4.3 Create `Tools/AlignmentTools.cs`
- [ ] 3.4.4 Create `Tools/MotifTools.cs`
- [ ] 3.4.5 Create `Tools/RepeatTools.cs`
- [ ] 3.4.6 Implement all 15 tools
- [ ] 3.4.7 Write 30 tests (15 tools × 2)
- [ ] 3.4.8 Create 15 `.mcp.json` files
- [ ] 3.4.9 Create 15 `.md` files
- [ ] 3.4.10 Build and verify
- [ ] 3.4.11 Create README.md

### 3.5 Server: Variants (22 tools)

- [ ] 3.5.1 Create project `SuffixTree.Mcp.Variants`
- [ ] 3.5.2 Add `Program.cs`
- [ ] 3.5.3 Create `Tools/VariantCallerTools.cs`
- [ ] 3.5.4 Create `Tools/VariantAnnotatorTools.cs`
- [ ] 3.5.5 Create `Tools/QualityTools.cs`
- [ ] 3.5.6 Implement all 22 tools
- [ ] 3.5.7 Write 44 tests (22 tools × 2)
- [ ] 3.5.8 Create 22 `.mcp.json` files
- [ ] 3.5.9 Create 22 `.md` files
- [ ] 3.5.10 Build and verify
- [ ] 3.5.11 Create README.md

### 3.6 Server: MolBio (38 tools)

- [ ] 3.6.1 Create project `SuffixTree.Mcp.MolBio`
- [ ] 3.6.2 Add `Program.cs`
- [ ] 3.6.3 Create `Tools/PrimerTools.cs`
- [ ] 3.6.4 Create `Tools/ProbeTools.cs`
- [ ] 3.6.5 Create `Tools/RestrictionTools.cs`
- [ ] 3.6.6 Create `Tools/CrisprTools.cs`
- [ ] 3.6.7 Create `Tools/CodonTools.cs`
- [ ] 3.6.8 Create `Tools/ThermoTools.cs`
- [ ] 3.6.9 Create `Tools/ReportTools.cs`
- [ ] 3.6.10 Implement all 38 tools
- [ ] 3.6.11 Write 76 tests (38 tools × 2)
- [ ] 3.6.12 Create 38 `.mcp.json` files
- [ ] 3.6.13 Create 38 `.md` files
- [ ] 3.6.14 Build and verify
- [ ] 3.6.15 Create README.md

### 3.7 Server: Assembly (13 tools)

- [ ] 3.7.1 Create project `SuffixTree.Mcp.Assembly`
- [ ] 3.7.2 Add `Program.cs`
- [ ] 3.7.3 Create `Tools/AssemblerTools.cs`
- [ ] 3.7.4 Create `Tools/AssemblyAnalyzerTools.cs`
- [ ] 3.7.5 Create `Tools/PanGenomeTools.cs`
- [ ] 3.7.6 Implement all 13 tools
- [ ] 3.7.7 Write 26 tests (13 tools × 2)
- [ ] 3.7.8 Create 13 `.mcp.json` files
- [ ] 3.7.9 Create 13 `.md` files
- [ ] 3.7.10 Build and verify
- [ ] 3.7.11 Create README.md

### 3.8 Server: Phylogenetics (10 tools)

- [ ] 3.8.1 Create project `SuffixTree.Mcp.Phylogenetics`
- [ ] 3.8.2 Add `Program.cs`
- [ ] 3.8.3 Create `Tools/TreeTools.cs`
- [ ] 3.8.4 Create `Tools/ComparativeTools.cs`
- [ ] 3.8.5 Implement all 10 tools
- [ ] 3.8.6 Write 20 tests (10 tools × 2)
- [ ] 3.8.7 Create 10 `.mcp.json` files
- [ ] 3.8.8 Create 10 `.md` files
- [ ] 3.8.9 Build and verify
- [ ] 3.8.10 Create README.md

### 3.9 Server: Population (15 tools)

- [ ] 3.9.1 Create project `SuffixTree.Mcp.Population`
- [ ] 3.9.2 Add `Program.cs`
- [ ] 3.9.3 Create `Tools/PopulationGeneticsTools.cs`
- [ ] 3.9.4 Create `Tools/MetagenomicsTools.cs`
- [ ] 3.9.5 Implement all 15 tools
- [ ] 3.9.6 Write 30 tests (15 tools × 2)
- [ ] 3.9.7 Create 15 `.mcp.json` files
- [ ] 3.9.8 Create 15 `.md` files
- [ ] 3.9.9 Build and verify
- [ ] 3.9.10 Create README.md

### 3.10 Server: Epigenetics (14 tools)

- [ ] 3.10.1 Create project `SuffixTree.Mcp.Epigenetics`
- [ ] 3.10.2 Add `Program.cs`
- [ ] 3.10.3 Create `Tools/EpigeneticsTools.cs`
- [ ] 3.10.4 Create `Tools/MiRnaTools.cs`
- [ ] 3.10.5 Implement all 14 tools
- [ ] 3.10.6 Write 28 tests (14 tools × 2)
- [ ] 3.10.7 Create 14 `.mcp.json` files
- [ ] 3.10.8 Create 14 `.md` files
- [ ] 3.10.9 Build and verify
- [ ] 3.10.10 Create README.md

### 3.11 Server: Structure (14 tools)

- [ ] 3.11.1 Create project `SuffixTree.Mcp.Structure`
- [ ] 3.11.2 Add `Program.cs`
- [ ] 3.11.3 Create `Tools/RnaStructureTools.cs`
- [ ] 3.11.4 Create `Tools/ProteinMotifTools.cs`
- [ ] 3.11.5 Create `Tools/DisorderTools.cs`
- [ ] 3.11.6 Implement all 14 tools
- [ ] 3.11.7 Write 28 tests (14 tools × 2)
- [ ] 3.11.8 Create 14 `.mcp.json` files
- [ ] 3.11.9 Create 14 `.md` files
- [ ] 3.11.10 Build and verify
- [ ] 3.11.11 Create README.md

### 3.12 Server: Annotation (8 tools)

- [ ] 3.12.1 Create project `SuffixTree.Mcp.Annotation`
- [ ] 3.12.2 Add `Program.cs`
- [ ] 3.12.3 Create `Tools/AnnotationTools.cs`
- [ ] 3.12.4 Create `Tools/SpliceSiteTools.cs`
- [ ] 3.12.5 Create `Tools/ChromosomeTools.cs`
- [ ] 3.12.6 Implement all 8 tools
- [ ] 3.12.7 Write 16 tests (8 tools × 2)
- [ ] 3.12.8 Create 8 `.mcp.json` files
- [ ] 3.12.9 Create 8 `.md` files
- [ ] 3.12.10 Build and verify
- [ ] 3.12.11 Create README.md

---

## 4. Per-Tool Checklist (All 241 Tools)

### Server 1: Core (12 tools)

#### 4.1.1 `suffix_tree_contains` - **Status: Ready**
- **HasDocs**: ✓
- **DocRef**: SuffixTree/SuffixTree.Search.cs#L22
- [x] a) Link MethodId: `SuffixTree.Contains`
- [x] b) Freeze toolName: `suffix_tree_contains`, serverName: `Core`
- [x] c) Define inputSchema: `{ text: string, pattern: string }`
- [x] d) Define outputSchema: `{ found: boolean }`
- [x] e) Define errors: 1001 (empty text), 1002 (null pattern)
- [x] f) Test: `SuffixTreeContains_Schema_ValidatesCorrectly`
- [x] g) Test: `SuffixTreeContains_Binding_InvokesSuccessfully`
- [x] h) Create `suffix_tree_contains.mcp.json`
- [x] i) Create `suffix_tree_contains.md`
- [x] j) Add traceability row

#### 4.1.2 `suffix_tree_count` - **Status: Ready**
- **HasDocs**: ✓
- **DocRef**: SuffixTree/SuffixTree.Search.cs#L132
- [x] a) Link MethodId: `SuffixTree.CountOccurrences`
- [x] b) Freeze toolName: `suffix_tree_count`, serverName: `Core`
- [x] c) Define inputSchema: `{ text: string, pattern: string }`
- [x] d) Define outputSchema: `{ count: integer }`
- [x] e) Define errors: 1001, 1002
- [x] f) Test: `SuffixTreeCount_Schema_ValidatesCorrectly`
- [x] g) Test: `SuffixTreeCount_Binding_InvokesSuccessfully`
- [x] h) Create `suffix_tree_count.mcp.json`
- [x] i) Create `suffix_tree_count.md`
- [x] j) Add traceability row

#### 4.1.3 `suffix_tree_find_all`
- **HasDocs**: ✓
- **DocRef**: SuffixTree.cs#xml
- [ ] a) Link MethodId: `SuffixTree.FindAllOccurrences`
- [ ] b) Freeze toolName: `suffix_tree_find_all`, serverName: `Core`
- [ ] c) Define inputSchema: `{ text: string, pattern: string }`
- [ ] d) Define outputSchema: `{ positions: integer[] }`
- [ ] e) Define errors: 1001, 1002
- [ ] f) Test: `SuffixTreeFindAll_Schema_ValidatesCorrectly`
- [ ] g) Test: `SuffixTreeFindAll_Binding_InvokesSuccessfully`
- [ ] h) Create `suffix_tree_find_all.mcp.json`
- [ ] i) Create `suffix_tree_find_all.md`
- [ ] j) Add traceability row

#### 4.1.4 `suffix_tree_lrs`
- **HasDocs**: ✓
- **DocRef**: SuffixTree.cs#xml
- [ ] a) Link MethodId: `SuffixTree.LongestRepeatedSubstring`
- [ ] b) Freeze toolName: `suffix_tree_lrs`, serverName: `Core`
- [ ] c) Define inputSchema: `{ text: string }`
- [ ] d) Define outputSchema: `{ substring: string, length: integer }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `SuffixTreeLrs_Schema_ValidatesCorrectly`
- [ ] g) Test: `SuffixTreeLrs_Binding_InvokesSuccessfully`
- [ ] h) Create `suffix_tree_lrs.mcp.json`
- [ ] i) Create `suffix_tree_lrs.md`
- [ ] j) Add traceability row

#### 4.1.5 `suffix_tree_lcs`
- **HasDocs**: ✓
- **DocRef**: SuffixTree.cs#xml
- [ ] a) Link MethodId: `SuffixTree.LongestCommonSubstring`
- [ ] b) Freeze toolName: `suffix_tree_lcs`, serverName: `Core`
- [ ] c) Define inputSchema: `{ text1: string, text2: string }`
- [ ] d) Define outputSchema: `{ substring: string, length: integer }`
- [ ] e) Define errors: 1001, 1003 (empty text2)
- [ ] f) Test: `SuffixTreeLcs_Schema_ValidatesCorrectly`
- [ ] g) Test: `SuffixTreeLcs_Binding_InvokesSuccessfully`
- [ ] h) Create `suffix_tree_lcs.mcp.json`
- [ ] i) Create `suffix_tree_lcs.md`
- [ ] j) Add traceability row

#### 4.1.6 `suffix_tree_stats`
- **HasDocs**: ✓
- **DocRef**: SuffixTree.cs#xml
- [ ] a) Link MethodId: `SuffixTree.Properties`
- [ ] b) Freeze toolName: `suffix_tree_stats`, serverName: `Core`
- [ ] c) Define inputSchema: `{ text: string }`
- [ ] d) Define outputSchema: `{ nodeCount: integer, depth: integer, textLength: integer }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `SuffixTreeStats_Schema_ValidatesCorrectly`
- [ ] g) Test: `SuffixTreeStats_Binding_InvokesSuccessfully`
- [ ] h) Create `suffix_tree_stats.mcp.json`
- [ ] i) Create `suffix_tree_stats.md`
- [ ] j) Add traceability row

#### 4.1.7 `find_longest_repeat`
- **HasDocs**: ✓
- **DocRef**: GenomicAnalyzer.cs:L20#xml
- [ ] a) Link MethodId: `GenomicAnalyzer.FindLongestRepeat`
- [ ] b) Freeze toolName: `find_longest_repeat`, serverName: `Core`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ repeat: string, position: integer, length: integer }`
- [ ] e) Define errors: 1001, 2001 (invalid sequence)
- [ ] f) Test: `FindLongestRepeat_Schema_ValidatesCorrectly`
- [ ] g) Test: `FindLongestRepeat_Binding_InvokesSuccessfully`
- [ ] h) Create `find_longest_repeat.mcp.json`
- [ ] i) Create `find_longest_repeat.md`
- [ ] j) Add traceability row

#### 4.1.8 `find_longest_common_region`
- **HasDocs**: ✓
- **DocRef**: GenomicAnalyzer.cs:L178#xml
- [ ] a) Link MethodId: `GenomicAnalyzer.FindLongestCommonRegion`
- [ ] b) Freeze toolName: `find_longest_common_region`, serverName: `Core`
- [ ] c) Define inputSchema: `{ sequence1: string, sequence2: string }`
- [ ] d) Define outputSchema: `{ region: string, position1: integer, position2: integer, length: integer }`
- [ ] e) Define errors: 1001, 1003, 2001
- [ ] f) Test: `FindLongestCommonRegion_Schema_ValidatesCorrectly`
- [ ] g) Test: `FindLongestCommonRegion_Binding_InvokesSuccessfully`
- [ ] h) Create `find_longest_common_region.mcp.json`
- [ ] i) Create `find_longest_common_region.md`
- [ ] j) Add traceability row

#### 4.1.9 `calculate_similarity`
- **HasDocs**: ✓
- **DocRef**: GenomicAnalyzer.cs:L238#xml
- [ ] a) Link MethodId: `GenomicAnalyzer.CalculateSimilarity`
- [ ] b) Freeze toolName: `calculate_similarity`, serverName: `Core`
- [ ] c) Define inputSchema: `{ sequence1: string, sequence2: string, kmerSize?: integer }`
- [ ] d) Define outputSchema: `{ similarity: number }`
- [ ] e) Define errors: 1001, 1003, 2001
- [ ] f) Test: `CalculateSimilarity_Schema_ValidatesCorrectly`
- [ ] g) Test: `CalculateSimilarity_Binding_InvokesSuccessfully`
- [ ] h) Create `calculate_similarity.mcp.json`
- [ ] i) Create `calculate_similarity.md`
- [ ] j) Add traceability row

#### 4.1.10 `hamming_distance`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Pattern_Matching/Approximate_Matching_Hamming.md
- [ ] a) Link MethodId: `ApproximateMatcher.HammingDistance`
- [ ] b) Freeze toolName: `hamming_distance`, serverName: `Core`
- [ ] c) Define inputSchema: `{ sequence1: string, sequence2: string }`
- [ ] d) Define outputSchema: `{ distance: integer }`
- [ ] e) Define errors: 1001, 1003, 1004 (length mismatch)
- [ ] f) Test: `HammingDistance_Schema_ValidatesCorrectly`
- [ ] g) Test: `HammingDistance_Binding_InvokesSuccessfully`
- [ ] h) Create `hamming_distance.mcp.json`
- [ ] i) Create `hamming_distance.md`
- [ ] j) Add traceability row

#### 4.1.11 `edit_distance`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Pattern_Matching/Edit_Distance.md
- [ ] a) Link MethodId: `ApproximateMatcher.EditDistance`
- [ ] b) Freeze toolName: `edit_distance`, serverName: `Core`
- [ ] c) Define inputSchema: `{ sequence1: string, sequence2: string }`
- [ ] d) Define outputSchema: `{ distance: integer }`
- [ ] e) Define errors: 1001, 1003
- [ ] f) Test: `EditDistance_Schema_ValidatesCorrectly`
- [ ] g) Test: `EditDistance_Binding_InvokesSuccessfully`
- [ ] h) Create `edit_distance.mcp.json`
- [ ] i) Create `edit_distance.md`
- [ ] j) Add traceability row

#### 4.1.12 `count_approximate_occurrences`
- **HasDocs**: ✓
- **DocRef**: ApproximateMatcher.cs:L283#xml
- [ ] a) Link MethodId: `ApproximateMatcher.CountApproximateOccurrences`
- [ ] b) Freeze toolName: `count_approximate_occurrences`, serverName: `Core`
- [ ] c) Define inputSchema: `{ sequence: string, pattern: string, maxMismatches: integer }`
- [ ] d) Define outputSchema: `{ count: integer }`
- [ ] e) Define errors: 1001, 1002, 1005 (invalid maxMismatches)
- [ ] f) Test: `CountApproximateOccurrences_Schema_ValidatesCorrectly`
- [ ] g) Test: `CountApproximateOccurrences_Binding_InvokesSuccessfully`
- [ ] h) Create `count_approximate_occurrences.mcp.json`
- [ ] i) Create `count_approximate_occurrences.md`
- [ ] j) Add traceability row

---

### Server 2: Sequence (35 tools)

#### 4.2.1 `dna_validate`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Sequence_Composition/Sequence_Validation.md
- [ ] a) Link MethodId: `DnaSequence.TryCreate`
- [ ] b) Freeze toolName: `dna_validate`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ valid: boolean, length: integer, error?: string }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `DnaValidate_Schema_ValidatesCorrectly`
- [ ] g) Test: `DnaValidate_Binding_InvokesSuccessfully`
- [ ] h) Create `dna_validate.mcp.json`
- [ ] i) Create `dna_validate.md`
- [ ] j) Add traceability row

#### 4.2.2 `dna_reverse_complement`
- **HasDocs**: ✓
- **DocRef**: DnaSequence.cs:L149#xml
- [ ] a) Link MethodId: `DnaSequence.GetReverseComplementString`
- [ ] b) Freeze toolName: `dna_reverse_complement`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ reverseComplement: string }`
- [ ] e) Define errors: 1001, 2001
- [ ] f) Test: `DnaReverseComplement_Schema_ValidatesCorrectly`
- [ ] g) Test: `DnaReverseComplement_Binding_InvokesSuccessfully`
- [ ] h) Create `dna_reverse_complement.mcp.json`
- [ ] i) Create `dna_reverse_complement.md`
- [ ] j) Add traceability row

#### 4.2.3 `rna_validate`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Sequence_Composition/Sequence_Validation.md
- [ ] a) Link MethodId: `RnaSequence.TryCreate`
- [ ] b) Freeze toolName: `rna_validate`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ valid: boolean, length: integer, error?: string }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `RnaValidate_Schema_ValidatesCorrectly`
- [ ] g) Test: `RnaValidate_Binding_InvokesSuccessfully`
- [ ] h) Create `rna_validate.mcp.json`
- [ ] i) Create `rna_validate.md`
- [ ] j) Add traceability row

#### 4.2.4 `rna_from_dna`
- **HasDocs**: ✓
- **DocRef**: RnaSequence.cs:L147#xml
- [ ] a) Link MethodId: `RnaSequence.FromDna`
- [ ] b) Freeze toolName: `rna_from_dna`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ dnaSequence: string }`
- [ ] d) Define outputSchema: `{ rnaSequence: string }`
- [ ] e) Define errors: 1001, 2001
- [ ] f) Test: `RnaFromDna_Schema_ValidatesCorrectly`
- [ ] g) Test: `RnaFromDna_Binding_InvokesSuccessfully`
- [ ] h) Create `rna_from_dna.mcp.json`
- [ ] i) Create `rna_from_dna.md`
- [ ] j) Add traceability row

#### 4.2.5 `protein_validate`
- **HasDocs**: ✓
- **DocRef**: ProteinSequence.cs:L357#xml
- [ ] a) Link MethodId: `ProteinSequence.TryCreate`
- [ ] b) Freeze toolName: `protein_validate`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ valid: boolean, length: integer, error?: string }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `ProteinValidate_Schema_ValidatesCorrectly`
- [ ] g) Test: `ProteinValidate_Binding_InvokesSuccessfully`
- [ ] h) Create `protein_validate.mcp.json`
- [ ] i) Create `protein_validate.md`
- [ ] j) Add traceability row

#### 4.2.6 `nucleotide_composition`
- **HasDocs**: ✓
- **DocRef**: SequenceStatistics.cs:L48#xml
- [ ] a) Link MethodId: `SequenceStatistics.CalculateNucleotideComposition`
- [ ] b) Freeze toolName: `nucleotide_composition`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ A: integer, T: integer, G: integer, C: integer, other: integer }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `NucleotideComposition_Schema_ValidatesCorrectly`
- [ ] g) Test: `NucleotideComposition_Binding_InvokesSuccessfully`
- [ ] h) Create `nucleotide_composition.mcp.json`
- [ ] i) Create `nucleotide_composition.md`
- [ ] j) Add traceability row

#### 4.2.7 `amino_acid_composition`
- **HasDocs**: ✓
- **DocRef**: SequenceStatistics.cs:L98#xml
- [ ] a) Link MethodId: `SequenceStatistics.CalculateAminoAcidComposition`
- [ ] b) Freeze toolName: `amino_acid_composition`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ composition: object }`
- [ ] e) Define errors: 1001, 2002 (invalid protein)
- [ ] f) Test: `AminoAcidComposition_Schema_ValidatesCorrectly`
- [ ] g) Test: `AminoAcidComposition_Binding_InvokesSuccessfully`
- [ ] h) Create `amino_acid_composition.mcp.json`
- [ ] i) Create `amino_acid_composition.md`
- [ ] j) Add traceability row

#### 4.2.8 `molecular_weight_protein`
- **HasDocs**: ✓
- **DocRef**: SequenceStatistics.cs:L159#xml
- [ ] a) Link MethodId: `SequenceStatistics.CalculateMolecularWeight`
- [ ] b) Freeze toolName: `molecular_weight_protein`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ molecularWeight: number, unit: string }`
- [ ] e) Define errors: 1001, 2002
- [ ] f) Test: `MolecularWeightProtein_Schema_ValidatesCorrectly`
- [ ] g) Test: `MolecularWeightProtein_Binding_InvokesSuccessfully`
- [ ] h) Create `molecular_weight_protein.mcp.json`
- [ ] i) Create `molecular_weight_protein.md`
- [ ] j) Add traceability row

#### 4.2.9 `molecular_weight_nucleotide`
- **HasDocs**: ✓
- **DocRef**: SequenceStatistics.cs:L180#xml
- [ ] a) Link MethodId: `SequenceStatistics.CalculateNucleotideMolecularWeight`
- [ ] b) Freeze toolName: `molecular_weight_nucleotide`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, isDna?: boolean }`
- [ ] d) Define outputSchema: `{ molecularWeight: number, unit: string }`
- [ ] e) Define errors: 1001, 2001
- [ ] f) Test: `MolecularWeightNucleotide_Schema_ValidatesCorrectly`
- [ ] g) Test: `MolecularWeightNucleotide_Binding_InvokesSuccessfully`
- [ ] h) Create `molecular_weight_nucleotide.mcp.json`
- [ ] i) Create `molecular_weight_nucleotide.md`
- [ ] j) Add traceability row

#### 4.2.10 `isoelectric_point`
- **HasDocs**: ✓
- **DocRef**: SequenceStatistics.cs:L228#xml
- [ ] a) Link MethodId: `SequenceStatistics.CalculateIsoelectricPoint`
- [ ] b) Freeze toolName: `isoelectric_point`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ pI: number }`
- [ ] e) Define errors: 1001, 2002
- [ ] f) Test: `IsoelectricPoint_Schema_ValidatesCorrectly`
- [ ] g) Test: `IsoelectricPoint_Binding_InvokesSuccessfully`
- [ ] h) Create `isoelectric_point.mcp.json`
- [ ] i) Create `isoelectric_point.md`
- [ ] j) Add traceability row

#### 4.2.11 `hydrophobicity`
- **HasDocs**: ✓
- **DocRef**: SequenceStatistics.cs:L306#xml
- [ ] a) Link MethodId: `SequenceStatistics.CalculateHydrophobicity`
- [ ] b) Freeze toolName: `hydrophobicity`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ gravy: number }`
- [ ] e) Define errors: 1001, 2002
- [ ] f) Test: `Hydrophobicity_Schema_ValidatesCorrectly`
- [ ] g) Test: `Hydrophobicity_Binding_InvokesSuccessfully`
- [ ] h) Create `hydrophobicity.mcp.json`
- [ ] i) Create `hydrophobicity.md`
- [ ] j) Add traceability row

#### 4.2.12 `thermodynamics`
- **HasDocs**: ✓
- **DocRef**: SequenceStatistics.cs:L381#xml
- [ ] a) Link MethodId: `SequenceStatistics.CalculateThermodynamics`
- [ ] b) Freeze toolName: `thermodynamics`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ deltaH: number, deltaS: number, deltaG: number }`
- [ ] e) Define errors: 1001, 2001
- [ ] f) Test: `Thermodynamics_Schema_ValidatesCorrectly`
- [ ] g) Test: `Thermodynamics_Binding_InvokesSuccessfully`
- [ ] h) Create `thermodynamics.mcp.json`
- [ ] i) Create `thermodynamics.md`
- [ ] j) Add traceability row

#### 4.2.13 `melting_temperature`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Molecular_Tools/Melting_Temperature.md
- [ ] a) Link MethodId: `SequenceStatistics.CalculateMeltingTemperature`
- [ ] b) Freeze toolName: `melting_temperature`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, useWallaceRule?: boolean }`
- [ ] d) Define outputSchema: `{ tm: number, unit: string }`
- [ ] e) Define errors: 1001, 2001
- [ ] f) Test: `MeltingTemperature_Schema_ValidatesCorrectly`
- [ ] g) Test: `MeltingTemperature_Binding_InvokesSuccessfully`
- [ ] h) Create `melting_temperature.mcp.json`
- [ ] i) Create `melting_temperature.md`
- [ ] j) Add traceability row

#### 4.2.14 `shannon_entropy`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Sequence_Composition/Shannon_Entropy.md
- [ ] a) Link MethodId: `SequenceStatistics.CalculateShannonEntropy`
- [ ] b) Freeze toolName: `shannon_entropy`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ entropy: number }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `ShannonEntropy_Schema_ValidatesCorrectly`
- [ ] g) Test: `ShannonEntropy_Binding_InvokesSuccessfully`
- [ ] h) Create `shannon_entropy.mcp.json`
- [ ] i) Create `shannon_entropy.md`
- [ ] j) Add traceability row

#### 4.2.15 `linguistic_complexity`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Sequence_Composition/Linguistic_Complexity.md
- [ ] a) Link MethodId: `SequenceStatistics.CalculateLinguisticComplexity`
- [ ] b) Freeze toolName: `linguistic_complexity`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, maxK?: integer }`
- [ ] d) Define outputSchema: `{ complexity: number }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `LinguisticComplexity_Schema_ValidatesCorrectly`
- [ ] g) Test: `LinguisticComplexity_Binding_InvokesSuccessfully`
- [ ] h) Create `linguistic_complexity.mcp.json`
- [ ] i) Create `linguistic_complexity.md`
- [ ] j) Add traceability row

#### 4.2.16 `summarize_sequence`
- **HasDocs**: ✓
- **DocRef**: SequenceStatistics.cs:L775#xml
- [ ] a) Link MethodId: `SequenceStatistics.SummarizeNucleotideSequence`
- [ ] b) Freeze toolName: `summarize_sequence`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ length: integer, gcContent: number, composition: object }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `SummarizeSequence_Schema_ValidatesCorrectly`
- [ ] g) Test: `SummarizeSequence_Binding_InvokesSuccessfully`
- [ ] h) Create `summarize_sequence.mcp.json`
- [ ] i) Create `summarize_sequence.md`
- [ ] j) Add traceability row

#### 4.2.17 `gc_content`
- **HasDocs**: ✓
- **DocRef**: SequenceExtensions.cs:L41#xml
- [ ] a) Link MethodId: `SequenceExtensions.CalculateGcContentFast`
- [ ] b) Freeze toolName: `gc_content`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ gcContent: number, gcCount: integer, totalCount: integer }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `GcContent_Schema_ValidatesCorrectly`
- [ ] g) Test: `GcContent_Binding_InvokesSuccessfully`
- [ ] h) Create `gc_content.mcp.json`
- [ ] i) Create `gc_content.md`
- [ ] j) Add traceability row

#### 4.2.18 `complement_base`
- **HasDocs**: ✓
- **DocRef**: SequenceExtensions.cs:L83#xml
- [ ] a) Link MethodId: `SequenceExtensions.GetComplementBase`
- [ ] b) Freeze toolName: `complement_base`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ nucleotide: string }`
- [ ] d) Define outputSchema: `{ complement: string }`
- [ ] e) Define errors: 1006 (invalid nucleotide)
- [ ] f) Test: `ComplementBase_Schema_ValidatesCorrectly`
- [ ] g) Test: `ComplementBase_Binding_InvokesSuccessfully`
- [ ] h) Create `complement_base.mcp.json`
- [ ] i) Create `complement_base.md`
- [ ] j) Add traceability row

#### 4.2.19 `is_valid_dna`
- **HasDocs**: ✓
- **DocRef**: SequenceExtensions.cs:L210#xml
- [ ] a) Link MethodId: `SequenceExtensions.IsValidDna`
- [ ] b) Freeze toolName: `is_valid_dna`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ valid: boolean }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `IsValidDna_Schema_ValidatesCorrectly`
- [ ] g) Test: `IsValidDna_Binding_InvokesSuccessfully`
- [ ] h) Create `is_valid_dna.mcp.json`
- [ ] i) Create `is_valid_dna.md`
- [ ] j) Add traceability row

#### 4.2.20 `is_valid_rna`
- **HasDocs**: ✓
- **DocRef**: SequenceExtensions.cs:L225#xml
- [ ] a) Link MethodId: `SequenceExtensions.IsValidRna`
- [ ] b) Freeze toolName: `is_valid_rna`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ valid: boolean }`
- [ ] e) Define errors: 1001
- [ ] f) Test: `IsValidRna_Schema_ValidatesCorrectly`
- [ ] g) Test: `IsValidRna_Binding_InvokesSuccessfully`
- [ ] h) Create `is_valid_rna.mcp.json`
- [ ] i) Create `is_valid_rna.md`
- [ ] j) Add traceability row

#### 4.2.21 `complexity_linguistic`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Sequence_Composition/Linguistic_Complexity.md
- [ ] a) Link MethodId: `SequenceComplexity.CalculateLinguisticComplexity`
- [ ] b) Freeze toolName: `complexity_linguistic`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, maxWordLength?: integer }`
- [ ] d) Define outputSchema: `{ complexity: number }`
- [ ] e) Define errors: 1001, 2001
- [ ] f) Test: `ComplexityLinguistic_Schema_ValidatesCorrectly`
- [ ] g) Test: `ComplexityLinguistic_Binding_InvokesSuccessfully`
- [ ] h) Create `complexity_linguistic.mcp.json`
- [ ] i) Create `complexity_linguistic.md`
- [ ] j) Add traceability row

#### 4.2.22 `complexity_shannon`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Sequence_Composition/Shannon_Entropy.md
- [ ] a) Link MethodId: `SequenceComplexity.CalculateShannonEntropy`
- [ ] b) Freeze toolName: `complexity_shannon`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ entropy: number }`
- [ ] e) Define errors: 1001, 2001
- [ ] f) Test: `ComplexityShannon_Schema_ValidatesCorrectly`
- [ ] g) Test: `ComplexityShannon_Binding_InvokesSuccessfully`
- [ ] h) Create `complexity_shannon.mcp.json`
- [ ] i) Create `complexity_shannon.md`
- [ ] j) Add traceability row

#### 4.2.23 `complexity_kmer_entropy`
- **HasDocs**: ✓
- **DocRef**: SequenceComplexity.cs:L128#xml
- [ ] a) Link MethodId: `SequenceComplexity.CalculateKmerEntropy`
- [ ] b) Freeze toolName: `complexity_kmer_entropy`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, k?: integer }`
- [ ] d) Define outputSchema: `{ entropy: number }`
- [ ] e) Define errors: 1001, 2001, 1007 (invalid k)
- [ ] f) Test: `ComplexityKmerEntropy_Schema_ValidatesCorrectly`
- [ ] g) Test: `ComplexityKmerEntropy_Binding_InvokesSuccessfully`
- [ ] h) Create `complexity_kmer_entropy.mcp.json`
- [ ] i) Create `complexity_kmer_entropy.md`
- [ ] j) Add traceability row

#### 4.2.24 `complexity_dust_score`
- **HasDocs**: ✓
- **DocRef**: SequenceComplexity.cs:L296#xml
- [ ] a) Link MethodId: `SequenceComplexity.CalculateDustScore`
- [ ] b) Freeze toolName: `complexity_dust_score`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, wordSize?: integer }`
- [ ] d) Define outputSchema: `{ score: number }`
- [ ] e) Define errors: 1001, 2001
- [ ] f) Test: `ComplexityDustScore_Schema_ValidatesCorrectly`
- [ ] g) Test: `ComplexityDustScore_Binding_InvokesSuccessfully`
- [ ] h) Create `complexity_dust_score.mcp.json`
- [ ] i) Create `complexity_dust_score.md`
- [ ] j) Add traceability row

#### 4.2.25 `complexity_mask_low`
- **HasDocs**: ✓
- **DocRef**: SequenceComplexity.cs:L346#xml
- [ ] a) Link MethodId: `SequenceComplexity.MaskLowComplexity`
- [ ] b) Freeze toolName: `complexity_mask_low`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, threshold?: number, maskChar?: string }`
- [ ] d) Define outputSchema: `{ maskedSequence: string, maskedCount: integer }`
- [ ] e) Define errors: 1001, 2001
- [ ] f) Test: `ComplexityMaskLow_Schema_ValidatesCorrectly`
- [ ] g) Test: `ComplexityMaskLow_Binding_InvokesSuccessfully`
- [ ] h) Create `complexity_mask_low.mcp.json`
- [ ] i) Create `complexity_mask_low.md`
- [ ] j) Add traceability row

#### 4.2.26 `complexity_compression_ratio`
- **HasDocs**: ✓
- **DocRef**: SequenceComplexity.cs:L391#xml
- [ ] a) Link MethodId: `SequenceComplexity.EstimateCompressionRatio`
- [ ] b) Freeze toolName: `complexity_compression_ratio`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string }`
- [ ] d) Define outputSchema: `{ ratio: number }`
- [ ] e) Define errors: 1001, 2001
- [ ] f) Test: `ComplexityCompressionRatio_Schema_ValidatesCorrectly`
- [ ] g) Test: `ComplexityCompressionRatio_Binding_InvokesSuccessfully`
- [ ] h) Create `complexity_compression_ratio.mcp.json`
- [ ] i) Create `complexity_compression_ratio.md`
- [ ] j) Add traceability row

#### 4.2.27 `kmer_count`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/K-mer_Analysis/K-mer_Counting.md
- [ ] a) Link MethodId: `KmerAnalyzer.CountKmers`
- [ ] b) Freeze toolName: `kmer_count`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, k: integer }`
- [ ] d) Define outputSchema: `{ kmers: object, total: integer }`
- [ ] e) Define errors: 1001, 1007
- [ ] f) Test: `KmerCount_Schema_ValidatesCorrectly`
- [ ] g) Test: `KmerCount_Binding_InvokesSuccessfully`
- [ ] h) Create `kmer_count.mcp.json`
- [ ] i) Create `kmer_count.md`
- [ ] j) Add traceability row

#### 4.2.28 `kmer_distance`
- **HasDocs**: ✓
- **DocRef**: KmerAnalyzer.cs:L165#xml
- [ ] a) Link MethodId: `KmerAnalyzer.KmerDistance`
- [ ] b) Freeze toolName: `kmer_distance`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence1: string, sequence2: string, k: integer }`
- [ ] d) Define outputSchema: `{ distance: number }`
- [ ] e) Define errors: 1001, 1003, 1007
- [ ] f) Test: `KmerDistance_Schema_ValidatesCorrectly`
- [ ] g) Test: `KmerDistance_Binding_InvokesSuccessfully`
- [ ] h) Create `kmer_distance.mcp.json`
- [ ] i) Create `kmer_distance.md`
- [ ] j) Add traceability row

#### 4.2.29 `kmer_entropy`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/K-mer_Analysis/K-mer_Frequency_Analysis.md
- [ ] a) Link MethodId: `KmerAnalyzer.CalculateKmerEntropy`
- [ ] b) Freeze toolName: `kmer_entropy`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, k: integer }`
- [ ] d) Define outputSchema: `{ entropy: number }`
- [ ] e) Define errors: 1001, 1007
- [ ] f) Test: `KmerEntropy_Schema_ValidatesCorrectly`
- [ ] g) Test: `KmerEntropy_Binding_InvokesSuccessfully`
- [ ] h) Create `kmer_entropy.mcp.json`
- [ ] i) Create `kmer_entropy.md`
- [ ] j) Add traceability row

#### 4.2.30 `kmer_analyze`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/K-mer_Analysis/K-mer_Search.md
- [ ] a) Link MethodId: `KmerAnalyzer.AnalyzeKmers`
- [ ] b) Freeze toolName: `kmer_analyze`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, k: integer }`
- [ ] d) Define outputSchema: `{ uniqueKmers: integer, totalKmers: integer, entropy: number, mostFrequent: array }`
- [ ] e) Define errors: 1001, 1007
- [ ] f) Test: `KmerAnalyze_Schema_ValidatesCorrectly`
- [ ] g) Test: `KmerAnalyze_Binding_InvokesSuccessfully`
- [ ] h) Create `kmer_analyze.mcp.json`
- [ ] i) Create `kmer_analyze.md`
- [ ] j) Add traceability row

#### 4.2.31 `iupac_code`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Matching.md
- [ ] a) Link MethodId: `ISequence.GetIupacCode`
- [ ] b) Freeze toolName: `iupac_code`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ bases: string[] }`
- [ ] d) Define outputSchema: `{ code: string }`
- [ ] e) Define errors: 1008 (invalid bases)
- [ ] f) Test: `IupacCode_Schema_ValidatesCorrectly`
- [ ] g) Test: `IupacCode_Binding_InvokesSuccessfully`
- [ ] h) Create `iupac_code.mcp.json`
- [ ] i) Create `iupac_code.md`
- [ ] j) Add traceability row

#### 4.2.32 `iupac_match`
- **HasDocs**: ✓
- **DocRef**: docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Matching.md
- [ ] a) Link MethodId: `ISequence.CodesMatch`
- [ ] b) Freeze toolName: `iupac_match`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ code1: string, code2: string }`
- [ ] d) Define outputSchema: `{ match: boolean }`
- [ ] e) Define errors: 1008
- [ ] f) Test: `IupacMatch_Schema_ValidatesCorrectly`
- [ ] g) Test: `IupacMatch_Binding_InvokesSuccessfully`
- [ ] h) Create `iupac_match.mcp.json`
- [ ] i) Create `iupac_match.md`
- [ ] j) Add traceability row

#### 4.2.33 `iupac_matches`
- **HasDocs**: ✓
- **DocRef**: IupacHelper.cs:L14#xml
- [ ] a) Link MethodId: `IupacHelper.MatchesIupac`
- [ ] b) Freeze toolName: `iupac_matches`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ nucleotide: string, iupacCode: string }`
- [ ] d) Define outputSchema: `{ matches: boolean }`
- [ ] e) Define errors: 1006, 1008
- [ ] f) Test: `IupacMatches_Schema_ValidatesCorrectly`
- [ ] g) Test: `IupacMatches_Binding_InvokesSuccessfully`
- [ ] h) Create `iupac_matches.mcp.json`
- [ ] i) Create `iupac_matches.md`
- [ ] j) Add traceability row

#### 4.2.34 `translate_dna`
- **HasDocs**: ✓
- **DocRef**: Translator.cs:L20#xml
- [ ] a) Link MethodId: `Translator.Translate(DnaSequence)`
- [ ] b) Freeze toolName: `translate_dna`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, geneticCodeTable?: integer, readingFrame?: integer }`
- [ ] d) Define outputSchema: `{ protein: string, length: integer }`
- [ ] e) Define errors: 1001, 2001, 1009 (invalid frame)
- [ ] f) Test: `TranslateDna_Schema_ValidatesCorrectly`
- [ ] g) Test: `TranslateDna_Binding_InvokesSuccessfully`
- [ ] h) Create `translate_dna.mcp.json`
- [ ] i) Create `translate_dna.md`
- [ ] j) Add traceability row

#### 4.2.35 `translate_rna`
- **HasDocs**: ✓
- **DocRef**: Translator.cs:L37#xml
- [ ] a) Link MethodId: `Translator.Translate(RnaSequence)`
- [ ] b) Freeze toolName: `translate_rna`, serverName: `Sequence`
- [ ] c) Define inputSchema: `{ sequence: string, geneticCodeTable?: integer, readingFrame?: integer }`
- [ ] d) Define outputSchema: `{ protein: string, length: integer }`
- [ ] e) Define errors: 1001, 2003 (invalid RNA), 1009
- [ ] f) Test: `TranslateRna_Schema_ValidatesCorrectly`
- [ ] g) Test: `TranslateRna_Binding_InvokesSuccessfully`
- [ ] h) Create `translate_rna.mcp.json`
- [ ] i) Create `translate_rna.md`
- [ ] j) Add traceability row

---

### Server 3: Parsers (45 tools)

#### 4.3.1 `fasta_parse`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.2 `fasta_format`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.3 `fasta_write`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.4 `fastq_parse`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.5 `fastq_detect_encoding`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.6 `fastq_encode_quality`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.7 `fastq_phred_to_error`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.8 `fastq_error_to_phred`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.9 `fastq_trim_quality`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.10 `fastq_trim_adapter`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.11 `fastq_statistics`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.12 `fastq_write`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.13 `fastq_format`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.14 `genbank_parse`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.15 `genbank_parse_location`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.16 `genbank_extract_sequence`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.17 `embl_parse`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.18 `embl_parse_location`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.19 `embl_extract_sequence`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.20 `gff_parse`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.21 `gff_statistics`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.22 `gff_write`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.23 `gff_extract_sequence`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.24 `bed_parse`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.25 `bed_block_length`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.26 `bed_statistics`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.27 `bed_write`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.28 `bed_extract_sequence`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.29 `vcf_parse`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.30 `vcf_classify`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.31 `vcf_is_snp`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.32 `vcf_is_indel`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.33 `vcf_variant_length`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.34 `vcf_is_hom_ref`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.35 `vcf_is_hom_alt`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.36 `vcf_is_het`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.37 `vcf_statistics`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.38 `vcf_write`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.39 `vcf_has_flag`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.40 `to_genbank`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.41 `to_bed`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.42 `to_gff`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.43 `to_phylip`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.44 `quality_char_to_phred`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.3.45 `quality_phred_to_char`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

---

### Server 4: Alignment (15 tools)

#### 4.4.1 `align_global`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.2 `align_local`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.3 `align_semi_global`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.4 `align_multiple`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.5 `align_statistics`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.6 `align_format`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.7 `align_global_cancellable`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.8 `motif_find_exact`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.9 `motif_find_degenerate`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.10 `motif_create_pwm`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.11 `motif_scan_pwm`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.12 `motif_consensus`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.13 `repeat_find`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.14 `repeat_summary`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

#### 4.4.15 `find_with_mismatches`
- **HasDocs**: ✓
- [ ] a-j) [Standard checklist items]

---

### Server 5: Variants (22 tools)

#### 4.5.1-4.5.22
All tools listed in plan with HasDocs=✓:
`variant_call`, `variant_classify_mutation`, `variant_titv_ratio`, `variant_statistics`, `variant_predict_effect`, `annotator_classify`, `annotator_normalize`, `annotator_impact_level`, `annotator_pathogenicity`, `annotator_parse_vcf`, `annotator_format_vcf`, `quality_statistics`, `quality_trim`, `quality_sliding_trim`, `quality_expected_errors`, `quality_mask_low`, `quality_consensus`, `quality_detect_encoding`, `quality_phred_to_error`, `quality_error_to_phred`, `quality_phred_string`, `sam_decode_flags`

- [ ] a-j) [Standard checklist items for each]

---

### Server 6: MolBio (38 tools)

#### 4.6.1-4.6.38
All tools listed in plan with HasDocs=✓:
`primer_design`, `primer_evaluate`, `primer_tm`, `primer_tm_salt`, `primer_gc`, `primer_homopolymer`, `primer_dinucleotide_repeat`, `primer_hairpin`, `primer_dimer`, `primer_3prime_stability`, `probe_design_tiling`, `probe_validate`, `probe_specificity`, `probe_molecular_weight`, `probe_extinction`, `probe_concentration`, `restriction_digest_summary`, `restriction_create_map`, `restriction_compatible`, `crispr_system`, `crispr_evaluate_guide`, `crispr_specificity`, `codon_cai`, `codon_enc`, `codon_statistics`, `codon_optimize`, `codon_cai_optimizer`, `codon_remove_restriction`, `codon_reduce_structure`, `codon_compare_usage`, `codon_table_from_sequence`, `genetic_code_get`, `thermo_wallace_tm`, `thermo_marmur_doty_tm`, `thermo_salt_adjusted_tm`, `thermo_salt_correction`, `report_create`, `report_sequence_analysis`

- [ ] a-j) [Standard checklist items for each]

---

### Server 7: Assembly (13 tools)

#### 4.7.1-4.7.13
All tools listed in plan with HasDocs=✓:
`assemble_olc`, `assemble_de_bruijn`, `assemble_identity`, `assemble_stats`, `assemble_merge`, `assemble_consensus`, `assembly_statistics`, `assembly_nx`, `assembly_aun`, `assembly_completeness`, `assembly_compare`, `pangenome_construct`, `pangenome_heaps_law`

- [ ] a-j) [Standard checklist items for each]

---

### Server 8: Phylogenetics (10 tools)

#### 4.8.1-4.8.10
All tools listed in plan with HasDocs=✓:
`tree_build`, `tree_pairwise_distance`, `tree_to_newick`, `tree_parse_newick`, `tree_length`, `tree_depth`, `tree_robinson_foulds`, `tree_patristic_distance`, `comparative_genomes`, `comparative_ani`

- [ ] a-j) [Standard checklist items for each]

---

### Server 9: Population (15 tools)

#### 4.9.1-4.9.15
All tools listed in plan with HasDocs=✓:
`pop_maf`, `pop_nucleotide_diversity`, `pop_watterson_theta`, `pop_tajimas_d`, `pop_diversity_stats`, `pop_hardy_weinberg`, `pop_fst`, `pop_f_statistics`, `pop_ld`, `pop_ihs`, `pop_inbreeding_roh`, `meta_taxonomic_profile`, `meta_alpha_diversity`, `meta_beta_diversity`, `comparative_reversal_distance`

- [ ] a-j) [Standard checklist items for each]

---

### Server 10: Epigenetics (14 tools)

#### 4.10.1-4.10.14
All tools listed in plan with HasDocs=✓:
`epi_cpg_observed_expected`, `epi_bisulfite_simulation`, `epi_methylation_profile`, `epi_chromatin_state`, `epi_epigenetic_age`, `mirna_seed`, `mirna_create`, `mirna_reverse_complement`, `mirna_can_pair`, `mirna_wobble_pair`, `mirna_align_target`, `mirna_site_accessibility`, `mirna_gc_content`, `pangenome_core_alignment`

- [ ] a-j) [Standard checklist items for each]

---

### Server 11: Structure (14 tools)

#### 4.11.1-4.11.14
All tools listed in plan with HasDocs=✓:
`rna_can_pair`, `rna_complement`, `rna_stem_energy`, `rna_hairpin_energy`, `rna_mfe`, `rna_predict_structure`, `rna_validate_dotbracket`, `rna_structure_probability`, `rna_generate_random`, `protein_motif_find`, `protein_prosite_to_regex`, `protein_disorder_predict`, `protein_disorder_propensity`, `protein_disorder_promoting`

- [ ] a-j) [Standard checklist items for each]

---

### Server 12: Annotation (8 tools)

#### 4.12.1-4.12.8
All tools listed in plan with HasDocs=✓:
`annotation_coding_potential`, `splice_predict_structure`, `splice_maxent_score`, `splice_within_coding`, `chromosome_karyotype`, `chromosome_telomeres`, `chromosome_centromere`, `chromosome_arm_ratio`

- [ ] a-j) [Standard checklist items for each]

---

## 5. Final Verification (Quality Gates)

### G1: Coverage Gate
- [ ] 5.1.1 Count all tools = 241
- [ ] 5.1.2 Verify distribution:
  - [ ] Core = 12
  - [ ] Sequence = 35
  - [ ] Parsers = 45
  - [ ] Alignment = 15
  - [ ] Variants = 22
  - [ ] MolBio = 38
  - [ ] Assembly = 13
  - [ ] Phylogenetics = 10
  - [ ] Population = 15
  - [ ] Epigenetics = 14
  - [ ] Structure = 14
  - [ ] Annotation = 8

### G2: Docs-first Gate
- [ ] 5.2.1 All HasDocs=true tools listed first in each server
- [ ] 5.2.2 No HasDocs=false tools without remediation plan

### G3: Traceability Gate
- [ ] 5.3.1 Traceability matrix has 241 rows
- [ ] 5.3.2 All MethodId fields populated or marked TBD
- [ ] 5.3.3 All DocRef fields populated or marked TBD

### G4: Documentation Gate
- [ ] 5.4.1 241 `.mcp.json` files created
- [ ] 5.4.2 241 `.md` files created
- [ ] 5.4.3 All references valid (spot check 10%)
- [ ] 5.4.4 All examples complete (spot check 10%)

### G5: Tests Gate
- [ ] 5.5.1 482 tests created (241 × 2)
- [ ] 5.5.2 All tests passing
- [ ] 5.5.3 No business logic asserts (audit)

---

## Summary Statistics

| Metric | Count |
|--------|-------|
| Total Servers | 12 |
| Total Tools | 241 |
| Total Tests (minimum) | 482 |
| Total `.mcp.json` files | 241 |
| Total `.md` files | 241 |
| HasDocs=true | 241 (100%) |
| HasDocs=false | 0 (0%) |

---

## Notes

1. All 241 tools have XML documentation in source files
2. 31 tools have additional Markdown algorithm documentation
3. Standard checklist items (a-j) apply uniformly to all tools
4. Test count may increase for tools with overloads (justified per-tool)
