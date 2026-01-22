# Algorithms Checklist v2.0

**Дата:** 2026-01-22
**Версия:** 2.1 (100% Coverage)
**Библиотека:** SuffixTree.Genomics

---

## Quick Reference

| Метрика | Значение |
|---------|----------|
| **Total Test Units** | 100 |
| **Completed** | 6 |
| **In Progress** | 0 |
| **Blocked** | 0 |
| **Not Started** | 94 |

---

## Processing Registry

| Status | Test Unit ID | Area | Methods | Evidence | TestSpec | Test File(s) |
|--------|--------------|------|---------|----------|----------|--------------|
| ☑ | SEQ-GC-001 | Composition | 5 | Wikipedia, Biopython | [SEQ-GC-001.md](TestSpecs/SEQ-GC-001.md) | SequenceExtensions_CalculateGcContent_Tests.cs |
| ☑ | SEQ-COMP-001 | Composition | 3 | Wikipedia, Biopython | [SEQ-COMP-001.md](TestSpecs/SEQ-COMP-001.md) | SequenceExtensions_Complement_Tests.cs |
| ☑ | SEQ-REVCOMP-001 | Composition | 4 | Wikipedia, Biopython | [SEQ-REVCOMP-001.md](TestSpecs/SEQ-REVCOMP-001.md) | SequenceExtensions_ReverseComplement_Tests.cs |
| ☑ | SEQ-VALID-001 | Composition | 4 | Wikipedia, IUPAC 1970, Bioinformatics.org | [SEQ-VALID-001.md](TestSpecs/SEQ-VALID-001.md) | SequenceExtensions_SequenceValidation_Tests.cs |
| ☑ | SEQ-COMPLEX-001 | Composition | 4 | Wikipedia, Troyanskaya (2002), Orlov (2004) | [SEQ-COMPLEX-001.md](TestSpecs/SEQ-COMPLEX-001.md) | SequenceComplexityTests.cs |
| ☑ | SEQ-ENTROPY-001 | Composition | 2 | Wikipedia (Entropy, Sequence logo, K-mer), Shannon (1948) | [SEQ-ENTROPY-001.md](TestSpecs/SEQ-ENTROPY-001.md) | SequenceComplexityTests.cs |
| ☐ | SEQ-GCSKEW-001 | Composition | 4 | - | - | - |
| ☐ | PAT-EXACT-001 | Matching | 4 | - | - | - |
| ☐ | PAT-APPROX-001 | Matching | 2 | - | - | - |
| ☐ | PAT-APPROX-002 | Matching | 2 | - | - | - |
| ☐ | PAT-IUPAC-001 | Matching | 2 | - | - | - |
| ☐ | PAT-PWM-001 | Matching | 2 | - | - | - |
| ☐ | REP-STR-001 | Repeats | 4 | - | - | - |
| ☐ | REP-TANDEM-001 | Repeats | 2 | - | - | - |
| ☐ | REP-INV-001 | Repeats | 1 | - | - | - |
| ☐ | REP-DIRECT-001 | Repeats | 1 | - | - | - |
| ☐ | REP-PALIN-001 | Repeats | 2 | - | - | - |
| ☐ | CRISPR-PAM-001 | MolTools | 2 | - | - | - |
| ☐ | CRISPR-GUIDE-001 | MolTools | 2 | - | - | - |
| ☐ | CRISPR-OFF-001 | MolTools | 2 | - | - | - |
| ☐ | PRIMER-TM-001 | MolTools | 2 | - | - | - |
| ☐ | PRIMER-DESIGN-001 | MolTools | 3 | - | - | - |
| ☐ | PRIMER-STRUCT-001 | MolTools | 3 | - | - | - |
| ☐ | PROBE-DESIGN-001 | MolTools | 3 | - | - | - |
| ☐ | PROBE-VALID-001 | MolTools | 2 | - | - | - |
| ☐ | RESTR-FIND-001 | MolTools | 2 | - | - | - |
| ☐ | RESTR-DIGEST-001 | MolTools | 2 | - | - | - |
| ☐ | ANNOT-ORF-001 | Annotation | 3 | - | - | - |
| ☐ | ANNOT-GENE-001 | Annotation | 2 | - | - | - |
| ☐ | ANNOT-PROM-001 | Annotation | 1 | - | - | - |
| ☐ | ANNOT-GFF-001 | Annotation | 2 | - | - | - |
| ☐ | KMER-COUNT-001 | K-mer | 3 | - | - | - |
| ☐ | KMER-FREQ-001 | K-mer | 3 | - | - | - |
| ☐ | KMER-FIND-001 | K-mer | 3 | - | - | - |
| ☐ | ALIGN-GLOBAL-001 | Alignment | 1 | - | - | - |
| ☐ | ALIGN-LOCAL-001 | Alignment | 1 | - | - | - |
| ☐ | ALIGN-SEMI-001 | Alignment | 1 | - | - | - |
| ☐ | ALIGN-MULTI-001 | Alignment | 1 | - | - | - |
| ☐ | PHYLO-DIST-001 | Phylogenetic | 2 | - | - | - |
| ☐ | PHYLO-TREE-001 | Phylogenetic | 1 | - | - | - |
| ☐ | PHYLO-NEWICK-001 | Phylogenetic | 2 | - | - | - |
| ☐ | PHYLO-COMP-001 | Phylogenetic | 3 | - | - | - |
| ☐ | POP-FREQ-001 | PopGen | 3 | - | - | - |
| ☐ | POP-DIV-001 | PopGen | 4 | - | - | - |
| ☐ | POP-HW-001 | PopGen | 1 | - | - | - |
| ☐ | POP-FST-001 | PopGen | 2 | - | - | - |
| ☐ | POP-LD-001 | PopGen | 2 | - | - | - |
| ☐ | CHROM-TELO-001 | Chromosome | 2 | - | - | - |
| ☐ | CHROM-CENT-001 | Chromosome | 1 | - | - | - |
| ☐ | CHROM-KARYO-001 | Chromosome | 2 | - | - | - |
| ☐ | CHROM-ANEU-001 | Chromosome | 2 | - | - | - |
| ☐ | CHROM-SYNT-001 | Chromosome | 2 | - | - | - |
| ☐ | META-CLASS-001 | Metagenomics | 2 | - | - | - |
| ☐ | META-PROF-001 | Metagenomics | 1 | - | - | - |
| ☐ | META-ALPHA-001 | Metagenomics | 1 | - | - | - |
| ☐ | META-BETA-001 | Metagenomics | 1 | - | - | - |
| ☐ | META-BIN-001 | Metagenomics | 1 | - | - | - |
| ☐ | CODON-OPT-001 | Codon | 1 | - | - | - |
| ☐ | CODON-CAI-001 | Codon | 1 | - | - | - |
| ☐ | CODON-RARE-001 | Codon | 1 | - | - | - |
| ☐ | CODON-USAGE-001 | Codon | 2 | - | - | - |
| ☐ | TRANS-CODON-001 | Translation | 3 | - | - | - |
| ☐ | TRANS-PROT-001 | Translation | 1 | - | - | - |
| ☐ | PARSE-FASTA-001 | FileIO | 4 | - | - | - |
| ☐ | PARSE-FASTQ-001 | FileIO | 4 | - | - | - |
| ☐ | PARSE-BED-001 | FileIO | 6 | - | - | - |
| ☐ | PARSE-VCF-001 | FileIO | 4 | - | - | - |
| ☐ | PARSE-GFF-001 | FileIO | 3 | - | - | - |
| ☐ | PARSE-GENBANK-001 | FileIO | 3 | - | - | - |
| ☐ | PARSE-EMBL-001 | FileIO | 2 | - | - | - |
| ☐ | RNA-STRUCT-001 | RnaStructure | 4 | - | - | - |
| ☐ | RNA-STEMLOOP-001 | RnaStructure | 3 | - | - | - |
| ☐ | RNA-ENERGY-001 | RnaStructure | 2 | - | - | - |
| ☐ | MIRNA-SEED-001 | MiRNA | 3 | - | - | - |
| ☐ | MIRNA-TARGET-001 | MiRNA | 2 | - | - | - |
| ☐ | MIRNA-PRECURSOR-001 | MiRNA | 2 | - | - | - |
| ☐ | SPLICE-DONOR-001 | Splicing | 2 | - | - | - |
| ☐ | SPLICE-ACCEPTOR-001 | Splicing | 2 | - | - | - |
| ☐ | SPLICE-PREDICT-001 | Splicing | 3 | - | - | - |
| ☐ | DISORDER-PRED-001 | ProteinPred | 3 | - | - | - |
| ☐ | DISORDER-REGION-001 | ProteinPred | 2 | - | - | - |
| ☐ | PROTMOTIF-FIND-001 | ProteinMotif | 3 | - | - | - |
| ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | 2 | - | - | - |
| ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | 2 | - | - | - |
| ☐ | EPIGEN-CPG-001 | Epigenetics | 3 | - | - | - |
| ☐ | EPIGEN-METHYL-001 | Epigenetics | 3 | - | - | - |
| ☐ | EPIGEN-DMR-001 | Epigenetics | 2 | - | - | - |
| ☐ | VARIANT-CALL-001 | Variants | 3 | - | - | - |
| ☐ | VARIANT-SNP-001 | Variants | 2 | - | - | - |
| ☐ | VARIANT-INDEL-001 | Variants | 2 | - | - | - |
| ☐ | VARIANT-ANNOT-001 | Variants | 2 | - | - | - |
| ☐ | SV-DETECT-001 | StructuralVar | 3 | - | - | - |
| ☐ | SV-BREAKPOINT-001 | StructuralVar | 2 | - | - | - |
| ☐ | SV-CNV-001 | StructuralVar | 2 | - | - | - |
| ☐ | ASSEMBLY-OLC-001 | Assembly | 2 | - | - | - |
| ☐ | ASSEMBLY-DBG-001 | Assembly | 2 | - | - | - |
| ☐ | ASSEMBLY-STATS-001 | Assembly | 4 | - | - | - |
| ☐ | TRANS-EXPR-001 | Transcriptome | 3 | - | - | - |
| ☐ | TRANS-DIFF-001 | Transcriptome | 2 | - | - | - |
| ☐ | TRANS-SPLICE-001 | Transcriptome | 2 | - | - | - |
| ☐ | COMPGEN-SYNTENY-001 | Comparative | 2 | - | - | - |
| ☐ | COMPGEN-ORTHO-001 | Comparative | 2 | - | - | - |
| ☐ | COMPGEN-REARR-001 | Comparative | 2 | - | - | - |
| ☐ | PANGEN-CORE-001 | PanGenome | 2 | - | - | - |
| ☐ | PANGEN-CLUSTER-001 | PanGenome | 2 | - | - | - |
| ☐ | QUALITY-PHRED-001 | Quality | 3 | - | - | - |
| ☐ | QUALITY-STATS-001 | Quality | 2 | - | - | - |

**Статусы:** ☐ Not Started | ⏳ In Progress | ☑ Complete | ⛔ Blocked

---

## Definition of Done (DoD)

### Обязательные критерии

| # | Критерий | Артефакт |
|---|----------|----------|
| 1 | TestSpec создан | `TestSpecs/{TestUnitID}.md` |
| 2 | Тесты написаны | `*.Tests/{Class}_{Method}_Tests.cs` |
| 3 | Покрытие веток ≥ 80% | Coverage report |
| 4 | Edge cases покрыты | null, empty, boundary, error |
| 5 | Тесты проходят | CI green |
| 6 | Evidence задокументирован | PR/commit link в Registry |

### Критерии качества тестов

- [ ] Тесты независимы (не зависят от порядка)
- [ ] Тесты детерминированы (без random без seed)
- [ ] Naming: `Method_Scenario_ExpectedResult`
- [ ] Структура: Arrange-Act-Assert
- [ ] Один assert на логическую проверку

### Для алгоритмов O(n²) и выше

- [ ] Property-based test для инварианта
- [ ] Performance baseline зафиксирован

---

## Test Units by Area

### 1. Sequence Composition (4 units)

#### SEQ-GC-001: GC Content Calculation

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceExtensions.CalculateGcContent(ReadOnlySpan<char>)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ result ≤ 100 (percentage) или 0 ≤ result ≤ 1 (fraction) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateGcContent(ReadOnlySpan<char>)` | SequenceExtensions | Canonical |
| `CalculateGcFraction(ReadOnlySpan<char>)` | SequenceExtensions | Variant (0-1) |
| `CalculateGcContentFast(string)` | SequenceExtensions | Delegate |
| `CalculateGcFractionFast(string)` | SequenceExtensions | Delegate |
| `GcContent` (property) | DnaSequence | Delegate |

**Edge Cases:**
- [ ] Empty sequence → 0
- [ ] All G/C → 100% / 1.0
- [ ] All A/T → 0% / 0.0
- [ ] Mixed case input
- [ ] Non-ACGT characters (N, etc.)

---

#### SEQ-COMP-001: DNA Complement

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceExtensions.GetComplementBase(char)` |
| **Complexity** | O(n) для последовательности |
| **Invariant** | Complement(Complement(x)) = x |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `GetComplementBase(char)` | SequenceExtensions | Canonical |
| `TryGetComplement(ReadOnlySpan, Span)` | SequenceExtensions | Span API |
| `Complement()` | DnaSequence | Instance |

**Edge Cases:**
- [ ] A ↔ T, G ↔ C (both directions)
- [ ] Case insensitivity (a → T)
- [ ] RNA support (U → A)
- [ ] Unknown base → unchanged
- [ ] Destination too small → false

---

#### SEQ-REVCOMP-001: Reverse Complement

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceExtensions.TryGetReverseComplement(ReadOnlySpan, Span)` |
| **Complexity** | O(n) |
| **Invariant** | RevComp(RevComp(x)) = x |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `TryGetReverseComplement(ReadOnlySpan, Span)` | SequenceExtensions | Canonical |
| `ReverseComplement()` | DnaSequence | Instance |
| `GetReverseComplementString(string)` | DnaSequence | Static helper |
| `TryWriteReverseComplement(Span)` | DnaSequence | Span API |

**Edge Cases:**
- [ ] Empty sequence
- [ ] Single nucleotide
- [ ] Palindrome (self-complementary)
- [ ] Destination too small → false

---

#### SEQ-VALID-001: Sequence Validation

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceExtensions.IsValidDna(ReadOnlySpan<char>)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `IsValidDna(ReadOnlySpan<char>)` | SequenceExtensions | Canonical DNA |
| `IsValidRna(ReadOnlySpan<char>)` | SequenceExtensions | Canonical RNA |
| `TryCreate(string, out DnaSequence)` | DnaSequence | Factory |
| `IsValid` (property) | DnaSequence | Instance |

**Edge Cases:**
- [ ] Empty → true (or false? define!)
- [ ] All valid bases
- [ ] Single invalid character
- [ ] Lowercase valid
- [ ] Whitespace handling

---

### 2. Pattern Matching (5 units)

#### PAT-EXACT-001: Exact Pattern Search

| Поле | Значение |
|------|----------|
| **Canonical** | `SuffixTree.FindAllOccurrences(string)` |
| **Complexity** | O(m + k) where m=pattern length, k=occurrences |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindAllOccurrences(string)` | SuffixTree/DnaSequence | Canonical |
| `Contains(string)` | SuffixTree/DnaSequence | Existence check |
| `CountOccurrences(string)` | SuffixTree/DnaSequence | Count only |
| `FindMotif(DnaSequence, string)` | GenomicAnalyzer | Wrapper |
| `FindExactMotif(DnaSequence, string)` | MotifFinder | Wrapper |

**Edge Cases:**
- [ ] Pattern not found → empty
- [ ] Pattern = entire sequence
- [ ] Overlapping occurrences
- [ ] Empty pattern
- [ ] Pattern longer than sequence

---

#### PAT-APPROX-001: Approximate Matching (Hamming)

| Поле | Значение |
|------|----------|
| **Canonical** | `ApproximateMatcher.FindWithMismatches(...)` |
| **Complexity** | O(n × m) |
| **Invariant** | HammingDistance требует equal length |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindWithMismatches(seq, pattern, max)` | ApproximateMatcher | Canonical |
| `HammingDistance(s1, s2)` | ApproximateMatcher | Distance |
| `HammingDistance(span1, span2)` | SequenceExtensions | Span API |

**Edge Cases:**
- [ ] maxMismatches = 0 → exact match
- [ ] maxMismatches ≥ pattern length → all positions
- [ ] Unequal lengths for HammingDistance → exception

---

#### PAT-APPROX-002: Approximate Matching (Edit Distance)

| Поле | Значение |
|------|----------|
| **Canonical** | `ApproximateMatcher.FindWithEdits(...)` |
| **Complexity** | O(n × m²) |
| **Invariant** | EditDistance(a,b) = EditDistance(b,a) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindWithEdits(seq, pattern, maxEdits)` | ApproximateMatcher | Canonical |
| `EditDistance(s1, s2)` | ApproximateMatcher | Distance |

**Edge Cases:**
- [ ] Identical strings → 0
- [ ] One empty string → length of other
- [ ] Single character difference
- [ ] Insertion vs deletion vs substitution

---

#### PAT-IUPAC-001: IUPAC Degenerate Matching

| Поле | Значение |
|------|----------|
| **Canonical** | `MotifFinder.FindDegenerateMotif(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindDegenerateMotif(seq, motif)` | MotifFinder | Canonical |
| `MatchesIupac(nucleotide, iupacCode)` | IupacHelper | Helper |

**Edge Cases:**
- [ ] All IUPAC codes: R, Y, S, W, K, M, B, D, H, V, N
- [ ] Mixed standard + IUPAC
- [ ] N matches everything

---

#### PAT-PWM-001: Position Weight Matrix

| Поле | Значение |
|------|----------|
| **Canonical** | `MotifFinder.ScanWithPwm(...)` |
| **Complexity** | O(n × m) where m=PWM width |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CreatePwm(alignedSequences)` | MotifFinder | Construction |
| `ScanWithPwm(seq, pwm, threshold)` | MotifFinder | Scanning |

**Edge Cases:**
- [ ] Empty alignment
- [ ] Single sequence alignment
- [ ] Threshold at boundary
- [ ] All positions below threshold

---

### 3. Repeat Analysis (5 units)

#### REP-STR-001: Microsatellite Detection (STR)

| Поле | Значение |
|------|----------|
| **Canonical** | `RepeatFinder.FindMicrosatellites(...)` |
| **Complexity** | O(n × U × R) where U=maxUnitLength, R=maxRepeats |
| **Invariant** | Result positions are non-overlapping (or document overlap policy) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindMicrosatellites(DnaSequence, ...)` | RepeatFinder | Canonical |
| `FindMicrosatellites(string, ...)` | RepeatFinder | Overload |
| `FindMicrosatellites(..., CancellationToken)` | RepeatFinder | Cancellable |
| `FindMicrosatellites(..., IProgress)` | RepeatFinder | With progress |

**Edge Cases:**
- [ ] No repeats found
- [ ] Entire sequence is one repeat
- [ ] minRepeats = 2 (minimum)
- [ ] Unit length 1-6 (mono to hexa)
- [ ] Cancellation mid-operation

---

#### REP-TANDEM-001: Tandem Repeat Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindTandemRepeats(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindTandemRepeats(seq, minUnit, maxUnit, minReps)` | GenomicAnalyzer | Canonical |
| `GetTandemRepeatSummary(seq, minRepeats)` | RepeatFinder | Summary |

---

#### REP-INV-001: Inverted Repeat Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `RepeatFinder.FindInvertedRepeats(...)` |
| **Complexity** | O(n² × L) where L=maxLoopLength |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindInvertedRepeats(seq, minArm, maxLoop)` | RepeatFinder | Canonical |

**Edge Cases:**
- [ ] Perfect palindrome (loop = 0)
- [ ] Maximum loop length
- [ ] Arm length at boundaries

---

#### REP-DIRECT-001: Direct Repeat Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `RepeatFinder.FindDirectRepeats(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindDirectRepeats(seq, minLen, maxLen, minSpacing)` | RepeatFinder | Canonical |

---

#### REP-PALIN-001: Palindrome Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `RepeatFinder.FindPalindromes(...)` |
| **Complexity** | O(n²) |
| **Invariant** | Palindrome = reverse complement equals self |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindPalindromes(seq, minLen, maxLen)` | RepeatFinder | Canonical |
| `FindPalindromes(seq, minLen, maxLen)` | GenomicAnalyzer | Alternate |

---

### 4. Molecular Tools (8 units)

#### CRISPR-PAM-001: PAM Site Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `CrisprDesigner.FindPamSites(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindPamSites(seq, systemType)` | CrisprDesigner | Canonical |
| `GetSystem(systemType)` | CrisprDesigner | System info |

**CRISPR Systems:**
- [ ] SpCas9 (NGG)
- [ ] SpCas9-NAG
- [ ] SaCas9 (NNGRRT)
- [ ] Cas12a (TTTV)
- [ ] AsCas12a, LbCas12a
- [ ] CasX

---

#### CRISPR-GUIDE-001: Guide RNA Design

| Поле | Значение |
|------|----------|
| **Canonical** | `CrisprDesigner.DesignGuideRnas(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `DesignGuideRnas(seq, start, end, type)` | CrisprDesigner | Canonical |
| `EvaluateGuideRna(guide, type, params)` | CrisprDesigner | Scoring |

---

#### CRISPR-OFF-001: Off-Target Analysis

| Поле | Значение |
|------|----------|
| **Canonical** | `CrisprDesigner.FindOffTargets(...)` |
| **Complexity** | O(n × m) с maxMismatches; может быть выше |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindOffTargets(guide, genome, maxMismatches)` | CrisprDesigner | Canonical |
| `CalculateSpecificityScore(guide, genome, type)` | CrisprDesigner | Score |

---

#### PRIMER-TM-001: Melting Temperature

| Поле | Значение |
|------|----------|
| **Canonical** | `PrimerDesigner.CalculateMeltingTemperature(...)` |
| **Complexity** | O(n) |
| **Formula** | Wallace rule (<14bp), Marmur-Doty (≥14bp) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateMeltingTemperature(primer)` | PrimerDesigner | Canonical |
| `CalculateMeltingTemperatureWithSalt(primer, Na)` | PrimerDesigner | Salt corrected |

**Constants (ThermoConstants):**
- `WallaceMaxLength` = 14
- `CalculateWallaceTm(at, gc)` = 2×AT + 4×GC
- `CalculateMarmurDotyTm(gc, len)` = 64.9 + 41×(GC-16.4)/len
- `CalculateSaltCorrection(Na)` = 16.6 × log10(Na/1000)

**Edge Cases:**
- [ ] Empty primer → 0
- [ ] Short primer (<14) uses Wallace
- [ ] Long primer (≥14) uses Marmur-Doty
- [ ] Salt concentration variations

---

#### PRIMER-DESIGN-001: Primer Pair Design

| Поле | Значение |
|------|----------|
| **Canonical** | `PrimerDesigner.DesignPrimers(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `DesignPrimers(template, start, end, params)` | PrimerDesigner | Canonical |
| `EvaluatePrimer(seq, pos, isForward, params)` | PrimerDesigner | Single primer |
| `GeneratePrimerCandidates(template, region)` | PrimerDesigner | All candidates |

**Parameters (PrimerParameters):**
- MinLength=18, MaxLength=25, OptimalLength=20
- MinGcContent=40, MaxGcContent=60
- MinTm=55, MaxTm=65, OptimalTm=60
- MaxHomopolymer=4, MaxDinucleotideRepeats=4

---

#### PRIMER-STRUCT-001: Primer Structure Analysis

| Поле | Значение |
|------|----------|
| **Canonical** | `PrimerDesigner.HasHairpinPotential(...)` |
| **Complexity** | O(m²) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `HasHairpinPotential(seq, minStemLength)` | PrimerDesigner | Hairpin |
| `HasPrimerDimer(primer1, primer2, minComp)` | PrimerDesigner | Dimer |
| `Calculate3PrimeStability(seq)` | PrimerDesigner | ΔG calculation |
| `FindLongestHomopolymer(seq)` | PrimerDesigner | Structure |
| `FindLongestDinucleotideRepeat(seq)` | PrimerDesigner | Structure |

---

#### RESTR-FIND-001: Restriction Site Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `RestrictionAnalyzer.FindSites(...)` |
| **Complexity** | O(n × k) where k=enzymes |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindSites(seq, enzymeNames)` | RestrictionAnalyzer | Canonical |
| `FindAllSites(seq)` | RestrictionAnalyzer | All 40+ enzymes |
| `GetEnzyme(name)` | RestrictionAnalyzer | Lookup |

**Enzyme Database:** 40+ enzymes (EcoRI, BamHI, HindIII, NotI, etc.)

---

#### RESTR-DIGEST-001: Digest Simulation

| Поле | Значение |
|------|----------|
| **Canonical** | `RestrictionAnalyzer.Digest(...)` |
| **Complexity** | O(n + k log k) where k=cut sites |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Digest(seq, enzymeNames)` | RestrictionAnalyzer | Canonical |
| `GetDigestSummary(seq, enzymeNames)` | RestrictionAnalyzer | Summary |
| `CreateMap(seq, enzymeNames)` | RestrictionAnalyzer | Full map |

---

### 5. Genome Annotation (4 units)

#### ANNOT-ORF-001: ORF Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `GenomeAnnotator.FindOrfs(...)` |
| **Complexity** | O(n) |
| **Invariant** | ORF starts with start codon, ends with stop codon |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindOrfs(dna, minLen, bothStrands, requireStart)` | GenomeAnnotator | Canonical |
| `FindLongestOrfsPerFrame(dna, bothStrands)` | GenomeAnnotator | Per-frame |
| `FindOpenReadingFrames(seq, minLen)` | GenomicAnalyzer | Alternate |

**Edge Cases:**
- [ ] No ORF found
- [ ] ORF extends to sequence end (no stop)
- [ ] All 6 reading frames
- [ ] Overlapping ORFs

---

#### ANNOT-GENE-001: Gene Prediction

| Поле | Значение |
|------|----------|
| **Canonical** | `GenomeAnnotator.PredictGenes(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `PredictGenes(dna, minOrfLen, prefix)` | GenomeAnnotator | Canonical |
| `FindRibosomeBindingSites(dna, window)` | GenomeAnnotator | RBS/SD |

---

#### ANNOT-PROM-001: Promoter Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `GenomeAnnotator.FindPromoterMotifs(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindPromoterMotifs(dna)` | GenomeAnnotator | Canonical |

**Motifs:** -35 box (TTGACA), -10 box (TATAAT)

---

#### ANNOT-GFF-001: GFF3 I/O

| Поле | Значение |
|------|----------|
| **Canonical** | `GenomeAnnotator.ParseGff3(...)` / `ToGff3(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `ParseGff3(lines)` | GenomeAnnotator | Parse |
| `ToGff3(annotations, seqId)` | GenomeAnnotator | Export |

---

### 6. K-mer Analysis (3 units)

#### KMER-COUNT-001: K-mer Counting

| Поле | Значение |
|------|----------|
| **Canonical** | `KmerAnalyzer.CountKmers(...)` |
| **Complexity** | O(n) |
| **Invariant** | Sum of counts = n - k + 1 |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CountKmers(seq, k)` | KmerAnalyzer | Canonical |
| `CountKmersSpan(seq, k)` | SequenceExtensions | Span API |
| `CountKmersBothStrands(dnaSeq, k)` | KmerAnalyzer | Both strands |

---

#### KMER-FREQ-001: K-mer Frequency Analysis

| Поле | Значение |
|------|----------|
| **Canonical** | `KmerAnalyzer.GetKmerSpectrum(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `GetKmerSpectrum(seq, k)` | KmerAnalyzer | Spectrum |
| `GetKmerFrequencies(seq, k)` | KmerAnalyzer | Normalized |
| `CalculateKmerEntropy(seq, k)` | KmerAnalyzer | Entropy |

---

#### KMER-FIND-001: K-mer Search

| Поле | Значение |
|------|----------|
| **Canonical** | `KmerAnalyzer.FindMostFrequentKmers(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindMostFrequentKmers(seq, k)` | KmerAnalyzer | Most frequent |
| `FindUniqueKmers(seq, k)` | KmerAnalyzer | Count = 1 |
| `FindClumps(seq, k, window, minOcc)` | KmerAnalyzer | Clumps |

---

### 7. Alignment (4 units)

#### ALIGN-GLOBAL-001: Global Alignment (Needleman-Wunsch)

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceAligner.GlobalAlign(...)` |
| **Complexity** | O(n × m) |
| **Invariant** | Optimal global alignment score |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `GlobalAlign(seq1, seq2, scoring)` | SequenceAligner | Canonical |

---

#### ALIGN-LOCAL-001: Local Alignment (Smith-Waterman)

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceAligner.LocalAlign(...)` |
| **Complexity** | O(n × m) |
| **Invariant** | Score ≥ 0, finds best local match |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `LocalAlign(seq1, seq2, scoring)` | SequenceAligner | Canonical |

---

#### ALIGN-SEMI-001: Semi-Global Alignment

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceAligner.SemiGlobalAlign(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `SemiGlobalAlign(seq1, seq2, scoring)` | SequenceAligner | Canonical |

---

#### ALIGN-MULTI-001: Multiple Sequence Alignment

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceAligner.MultipleAlign(...)` |
| **Complexity** | O(n² × m) progressive alignment |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `MultipleAlign(sequences)` | SequenceAligner | Canonical |

---

### 8. Phylogenetics (4 units)

#### PHYLO-DIST-001: Distance Matrix

| Поле | Значение |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.CalculateDistanceMatrix(...)` |
| **Complexity** | O(n² × m) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateDistanceMatrix(seqs, method)` | PhylogeneticAnalyzer | Canonical |
| `CalculatePairwiseDistance(s1, s2, method)` | PhylogeneticAnalyzer | Single pair |

**Distance Methods:** p-distance, Jukes-Cantor, Kimura 2-parameter, Hamming

---

#### PHYLO-TREE-001: Tree Construction

| Поле | Значение |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.BuildTree(...)` |
| **Complexity** | O(n³) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `BuildTree(seqs, distMethod, treeMethod)` | PhylogeneticAnalyzer | Canonical |

**Tree Methods:** UPGMA, Neighbor-Joining

---

#### PHYLO-NEWICK-001: Newick I/O

| Поле | Значение |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.ToNewick(...)` / `ParseNewick(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `ToNewick(treeNode)` | PhylogeneticAnalyzer | Export |
| `ParseNewick(newickString)` | PhylogeneticAnalyzer | Parse |

---

#### PHYLO-COMP-001: Tree Comparison

| Поле | Значение |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.RobinsonFouldsDistance(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `RobinsonFouldsDistance(tree1, tree2)` | PhylogeneticAnalyzer | Topology |
| `FindMRCA(root, taxon1, taxon2)` | PhylogeneticAnalyzer | MRCA |
| `PatristicDistance(root, t1, t2)` | PhylogeneticAnalyzer | Path length |

---

### 9. Population Genetics (5 units)

#### POP-FREQ-001: Allele Frequencies

| Поле | Значение |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(...)` |
| **Complexity** | O(n) |
| **Invariant** | p + q = 1 |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateAlleleFrequencies(hom_maj, het, hom_min)` | PopulationGeneticsAnalyzer | Canonical |
| `CalculateMAF(genotypes)` | PopulationGeneticsAnalyzer | MAF |
| `FilterByMAF(variants, min, max)` | PopulationGeneticsAnalyzer | Filter |

---

#### POP-DIV-001: Diversity Statistics

| Поле | Значение |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateDiversityStatistics(...)` |
| **Complexity** | O(n² × m) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateNucleotideDiversity(seqs)` | PopulationGeneticsAnalyzer | π |
| `CalculateWattersonTheta(segSites, n)` | PopulationGeneticsAnalyzer | θ |
| `CalculateTajimasD(pi, theta, S)` | PopulationGeneticsAnalyzer | D |
| `CalculateDiversityStatistics(seqs)` | PopulationGeneticsAnalyzer | All |

---

#### POP-HW-001: Hardy-Weinberg Test

| Поле | Значение |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.TestHardyWeinberg(...)` |
| **Complexity** | O(1) per variant |
| **Invariant** | Expected: p², 2pq, q² |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `TestHardyWeinberg(variantId, counts)` | PopulationGeneticsAnalyzer | Canonical |

---

#### POP-FST-001: F-Statistics

| Поле | Значение |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateFst(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ Fst ≤ 1 |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateFst(pop1, pop2)` | PopulationGeneticsAnalyzer | Canonical |
| `CalculateFStatistics(variantData)` | PopulationGeneticsAnalyzer | Fis, Fit, Fst |

---

#### POP-LD-001: Linkage Disequilibrium

| Поле | Значение |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateLD(...)` |
| **Complexity** | O(n) per pair |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateLD(var1, var2, genotypes)` | PopulationGeneticsAnalyzer | D', r² |
| `FindHaplotypeBlocks(variants)` | PopulationGeneticsAnalyzer | Blocks |

---

### 10. Chromosome Analysis (5 units)

#### CHROM-TELO-001: Telomere Analysis

| Поле | Значение |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.AnalyzeTelomeres(...)` |
| **Complexity** | O(n) |
| **Constant** | Human telomere repeat: TTAGGG |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `AnalyzeTelomeres(chrName, seq, repeat, ...)` | ChromosomeAnalyzer | Canonical |
| `EstimateTelomereLengthFromTSRatio(tsRatio)` | ChromosomeAnalyzer | qPCR estimate |

---

#### CHROM-CENT-001: Centromere Analysis

| Поле | Значение |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.AnalyzeCentromere(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `AnalyzeCentromere(chrName, seq, windowSize)` | ChromosomeAnalyzer | Canonical |

---

#### CHROM-KARYO-001: Karyotype Analysis

| Поле | Значение |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.AnalyzeKaryotype(...)` |
| **Complexity** | O(k) where k=chromosomes |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `AnalyzeKaryotype(chromosomes, ploidy)` | ChromosomeAnalyzer | Canonical |
| `DetectPloidy(depths, expected)` | ChromosomeAnalyzer | Ploidy |

---

#### CHROM-ANEU-001: Aneuploidy Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.DetectAneuploidy(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `DetectAneuploidy(depthData, medianDepth)` | ChromosomeAnalyzer | Canonical |
| `IdentifyWholeChromosomeAneuploidy(cnStates)` | ChromosomeAnalyzer | Classification |

---

#### CHROM-SYNT-001: Synteny Analysis

| Поле | Значение |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.FindSyntenyBlocks(...)` |
| **Complexity** | O(n log n) — requires verification |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindSyntenyBlocks(orthologPairs, minGenes)` | ChromosomeAnalyzer | Canonical |
| `DetectRearrangements(syntenyBlocks)` | ChromosomeAnalyzer | Rearrangements |

---

### 11. Metagenomics (5 units)

#### META-CLASS-001: Taxonomic Classification

| Поле | Значение |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.ClassifyReads(...)` |
| **Complexity** | O(n × m) where n=reads, m=read length |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `ClassifyReads(reads, kmerDB, k)` | MetagenomicsAnalyzer | Canonical |
| `BuildKmerDatabase(refGenomes, k)` | MetagenomicsAnalyzer | DB construction |

---

#### META-PROF-001: Taxonomic Profile

| Поле | Значение |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.GenerateTaxonomicProfile(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `GenerateTaxonomicProfile(classifications)` | MetagenomicsAnalyzer | Canonical |

---

#### META-ALPHA-001: Alpha Diversity

| Поле | Значение |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.CalculateAlphaDiversity(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateAlphaDiversity(abundances)` | MetagenomicsAnalyzer | Canonical |

**Indices:** Shannon, Simpson, Inverse Simpson, Chao1, Pielou's evenness

---

#### META-BETA-001: Beta Diversity

| Поле | Значение |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.CalculateBetaDiversity(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateBetaDiversity(sample1, sample2)` | MetagenomicsAnalyzer | Canonical |

**Metrics:** Bray-Curtis, Jaccard

---

#### META-BIN-001: Genome Binning

| Поле | Значение |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.BinContigs(...)` |
| **Complexity** | O(n × k × i) where k=bins, i=iterations — needs verification |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `BinContigs(contigs, numBins, minBinSize)` | MetagenomicsAnalyzer | Canonical |

---

### 12. Codon Optimization (4 units)

#### CODON-OPT-001: Sequence Optimization

| Поле | Значение |
|------|----------|
| **Canonical** | `CodonOptimizer.OptimizeSequence(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `OptimizeSequence(seq, organism, strategy)` | CodonOptimizer | Canonical |

**Strategies:** MaximizeCAI, BalancedOptimization, HarmonizeExpression, MinimizeSecondary, AvoidRareCodons

---

#### CODON-CAI-001: CAI Calculation

| Поле | Значение |
|------|----------|
| **Canonical** | `CodonOptimizer.CalculateCAI(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ CAI ≤ 1 |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateCAI(codingSeq, codonTable)` | CodonOptimizer | Canonical |

---

#### CODON-RARE-001: Rare Codon Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `CodonOptimizer.FindRareCodons(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindRareCodons(seq, codonTable, threshold)` | CodonOptimizer | Canonical |

---

#### CODON-USAGE-001: Codon Usage Analysis

| Поле | Значение |
|------|----------|
| **Canonical** | `CodonOptimizer.CalculateCodonUsage(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateCodonUsage(seq)` | CodonOptimizer | Canonical |
| `CompareCodonUsage(seq1, seq2)` | CodonOptimizer | Comparison |

**Organism Tables:** E. coli K12, S. cerevisiae, H. sapiens

---

### 13. Translation (2 units)

#### TRANS-CODON-001: Codon Translation

| Поле | Значение |
|------|----------|
| **Canonical** | `GeneticCode.Translate(...)` |
| **Complexity** | O(1) per codon |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Translate(codon)` | GeneticCode | Canonical |
| `IsStartCodon(codon)` | GeneticCode | Check |
| `IsStopCodon(codon)` | GeneticCode | Check |
| `GetCodonsForAminoAcid(aa)` | GeneticCode | Reverse lookup |

**Genetic Codes:** Standard (1), Vertebrate Mito (2), Yeast Mito (3), Bacterial (11)

---

#### TRANS-PROT-001: Protein Translation

| Поле | Значение |
|------|----------|
| **Canonical** | `Translator.Translate(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Translate(dna, geneticCode)` | Translator | Canonical |

---

### 14. File I/O Parsers (7 units)

#### PARSE-FASTA-001: FASTA Parsing

| Поле | Значение |
|------|----------|
| **Canonical** | `FastaParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | FastaParser |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Parse(fastaContent)` | FastaParser | Canonical |
| `ParseFile(filePath)` | FastaParser | File |
| `ParseFileAsync(filePath)` | FastaParser | Async |
| `ToFasta(entries, lineWidth)` | FastaParser | Export |
| `WriteFile(filePath, entries)` | FastaParser | Write |

**Edge Cases:**
- [ ] Empty content
- [ ] Multi-line sequences
- [ ] Special characters in headers
- [ ] Missing sequence after header

---

#### PARSE-FASTQ-001: FASTQ Parsing

| Поле | Значение |
|------|----------|
| **Canonical** | `FastqParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | FastqParser |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Parse(content, encoding)` | FastqParser | Canonical |
| `ParseFile(filePath, encoding)` | FastqParser | File |
| `CalculateStatistics(records)` | FastqParser | Stats |
| `FilterByQuality(records, minQ)` | FastqParser | Filter |

**Edge Cases:**
- [ ] Phred+33 vs Phred+64 encoding
- [ ] Malformed quality strings
- [ ] Empty records
- [ ] Auto-detect encoding

---

#### PARSE-BED-001: BED File Parsing

| Поле | Значение |
|------|----------|
| **Canonical** | `BedParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | BedParser |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Parse(content, format)` | BedParser | Canonical |
| `ParseFile(filePath, format)` | BedParser | File |
| `FilterByChrom(records, chrom)` | BedParser | Filter |
| `FilterByRegion(records, ...)` | BedParser | Filter |
| `MergeOverlapping(records)` | BedParser | Merge |
| `Intersect(records1, records2)` | BedParser | Set op |

**Edge Cases:**
- [ ] BED3 vs BED6 vs BED12 formats
- [ ] Block structures (BED12)
- [ ] Invalid coordinates

---

#### PARSE-VCF-001: VCF Parsing

| Поле | Значение |
|------|----------|
| **Canonical** | `VcfParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | VcfParser |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Parse(content)` | VcfParser | Canonical |
| `ParseFile(filePath)` | VcfParser | File |
| `ParseWithHeader(content)` | VcfParser | With header |
| `GetVariantType(record)` | VcfParser | Classification |

**Edge Cases:**
- [ ] Multi-allelic variants
- [ ] Missing values (.)
- [ ] Complex INFO fields
- [ ] Sample genotypes

---

#### PARSE-GFF-001: GFF/GTF Parsing

| Поле | Значение |
|------|----------|
| **Canonical** | `GffParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | GffParser |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Parse(content)` | GffParser | Canonical |
| `ParseFile(filePath)` | GffParser | File |
| `ToGff3(features)` | GffParser | Export |

---

#### PARSE-GENBANK-001: GenBank Parsing

| Поле | Значение |
|------|----------|
| **Canonical** | `GenBankParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | GenBankParser |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Parse(content)` | GenBankParser | Canonical |
| `ParseFile(filePath)` | GenBankParser | File |
| `ExtractFeatures(record)` | GenBankParser | Features |

---

#### PARSE-EMBL-001: EMBL Parsing

| Поле | Значение |
|------|----------|
| **Canonical** | `EmblParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | EmblParser |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Parse(content)` | EmblParser | Canonical |
| `ParseFile(filePath)` | EmblParser | File |

---

### 15. Sequence Complexity (3 units)

#### SEQ-COMPLEX-001: Linguistic Complexity

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceComplexity.CalculateLinguisticComplexity(...)` |
| **Complexity** | O(n × k) |
| **Invariant** | 0 ≤ result ≤ 1 |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateLinguisticComplexity(seq, maxLen)` | SequenceComplexity | Canonical |
| `CalculateLinguisticComplexity(string, maxLen)` | SequenceComplexity | String |
| `FindLowComplexityRegions(seq, window, threshold)` | SequenceComplexity | Regions |
| `MaskLowComplexity(seq, threshold)` | SequenceComplexity | Masking |

---

#### SEQ-ENTROPY-001: Shannon Entropy

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceComplexity.CalculateShannonEntropy(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ result ≤ 2 for DNA |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateShannonEntropy(sequence)` | SequenceComplexity | Canonical |
| `CalculateKmerEntropy(seq, k)` | SequenceComplexity | K-mer based |

---

#### SEQ-GCSKEW-001: GC Skew

| Поле | Значение |
|------|----------|
| **Canonical** | `GcSkewCalculator.CalculateGcSkew(...)` |
| **Complexity** | O(n) |
| **Invariant** | -1 ≤ result ≤ 1 |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateGcSkew(sequence)` | GcSkewCalculator | Canonical |
| `CalculateWindowedGcSkew(seq, window, step)` | GcSkewCalculator | Windowed |
| `CalculateCumulativeGcSkew(sequence)` | GcSkewCalculator | Cumulative |
| `FindOriginOfReplication(sequence)` | GcSkewCalculator | Origin detection |

---

### 16. RNA Secondary Structure (3 units)

#### RNA-STRUCT-001: Secondary Structure Prediction

| Поле | Значение |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.Predict(...)` |
| **Complexity** | O(n³) |
| **Class** | RnaSecondaryStructure |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Predict(sequence)` | RnaSecondaryStructure | Canonical |
| `PredictWithConstraints(seq, constraints)` | RnaSecondaryStructure | Constrained |
| `ToDotBracket(structure)` | RnaSecondaryStructure | Notation |
| `FromDotBracket(notation)` | RnaSecondaryStructure | Parse |

---

#### RNA-STEMLOOP-001: Stem-Loop Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.FindStemLoops(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindStemLoops(sequence, minStem, maxLoop)` | RnaSecondaryStructure | Canonical |
| `FindHairpins(sequence, params)` | RnaSecondaryStructure | Hairpins |
| `FindPseudoknots(sequence)` | RnaSecondaryStructure | Pseudoknots |

---

#### RNA-ENERGY-001: Free Energy Calculation

| Поле | Значение |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.CalculateFreeEnergy(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateFreeEnergy(structure)` | RnaSecondaryStructure | Canonical |
| `CalculateStackingEnergy(bp1, bp2)` | RnaSecondaryStructure | Stacking |

---

### 17. MicroRNA Analysis (3 units)

#### MIRNA-SEED-001: Seed Sequence Analysis

| Поле | Значение |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.GetSeedSequence(...)` |
| **Complexity** | O(1) |
| **Class** | MiRnaAnalyzer |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `GetSeedSequence(miRnaSequence)` | MiRnaAnalyzer | Canonical |
| `CreateMiRna(name, sequence)` | MiRnaAnalyzer | Factory |
| `CompareSeedRegions(mirna1, mirna2)` | MiRnaAnalyzer | Compare |

---

#### MIRNA-TARGET-001: Target Site Prediction

| Поле | Значение |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.FindTargetSites(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindTargetSites(mRna, miRna, minScore)` | MiRnaAnalyzer | Canonical |
| `ScoreTargetSite(site)` | MiRnaAnalyzer | Scoring |

**Site Types:**
- [ ] 8mer, 7mer-m8, 7mer-A1, 6mer
- [ ] Supplementary pairing
- [ ] Centered sites

---

#### MIRNA-PRECURSOR-001: Pre-miRNA Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.FindPreMiRnas(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindPreMiRnas(sequence)` | MiRnaAnalyzer | Canonical |
| `ValidateHairpin(structure)` | MiRnaAnalyzer | Validation |

---

### 18. Splice Site Prediction (3 units)

#### SPLICE-DONOR-001: Donor Site Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `SpliceSitePredictor.FindDonorSites(...)` |
| **Complexity** | O(n) |
| **Class** | SpliceSitePredictor |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindDonorSites(sequence, minScore)` | SpliceSitePredictor | Canonical |
| `ScoreDonorSite(context)` | SpliceSitePredictor | Scoring |

---

#### SPLICE-ACCEPTOR-001: Acceptor Site Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `SpliceSitePredictor.FindAcceptorSites(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindAcceptorSites(sequence, minScore)` | SpliceSitePredictor | Canonical |
| `ScoreAcceptorSite(context)` | SpliceSitePredictor | Scoring |

---

#### SPLICE-PREDICT-001: Gene Structure Prediction

| Поле | Значение |
|------|----------|
| **Canonical** | `SpliceSitePredictor.PredictGeneStructure(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `PredictGeneStructure(sequence)` | SpliceSitePredictor | Canonical |
| `FindIntrons(sequence)` | SpliceSitePredictor | Introns |
| `FindExons(sequence)` | SpliceSitePredictor | Exons |

---

### 19. Protein Disorder Prediction (2 units)

#### DISORDER-PRED-001: Disorder Prediction

| Поле | Значение |
|------|----------|
| **Canonical** | `DisorderPredictor.PredictDisorder(...)` |
| **Complexity** | O(n) |
| **Class** | DisorderPredictor |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `PredictDisorder(sequence, windowSize, threshold)` | DisorderPredictor | Canonical |
| `CalculateDisorderScore(window)` | DisorderPredictor | Score |
| `CalculateHydropathy(sequence)` | DisorderPredictor | Hydropathy |

---

#### DISORDER-REGION-001: Disordered Region Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `DisorderPredictor.IdentifyDisorderedRegions(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `IdentifyDisorderedRegions(predictions, minLen)` | DisorderPredictor | Canonical |
| `ClassifyRegionType(region)` | DisorderPredictor | Classification |

---

### 20. Protein Motif Finding (3 units)

#### PROTMOTIF-FIND-001: Motif Search

| Поле | Значение |
|------|----------|
| **Canonical** | `ProteinMotifFinder.FindMotifs(...)` |
| **Complexity** | O(n × m) |
| **Class** | ProteinMotifFinder |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindMotifs(sequence, patterns)` | ProteinMotifFinder | Canonical |
| `FindAllKnownMotifs(sequence)` | ProteinMotifFinder | All patterns |
| `ScanForPattern(sequence, pattern)` | ProteinMotifFinder | Single |

---

#### PROTMOTIF-PROSITE-001: PROSITE Pattern Matching

| Поле | Значение |
|------|----------|
| **Canonical** | `ProteinMotifFinder.MatchPrositePattern(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `MatchPrositePattern(sequence, pattern)` | ProteinMotifFinder | Canonical |
| `ParsePrositePattern(pattern)` | ProteinMotifFinder | Parse |

**Common Patterns:**
- [ ] N-glycosylation (PS00001)
- [ ] Phosphorylation sites
- [ ] Zinc fingers
- [ ] Signal peptides

---

#### PROTMOTIF-DOMAIN-001: Domain Prediction

| Поле | Значение |
|------|----------|
| **Canonical** | `ProteinMotifFinder.PredictDomains(...)` |
| **Complexity** | O(n × d) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `PredictDomains(sequence)` | ProteinMotifFinder | Canonical |
| `PredictSignalPeptide(sequence)` | ProteinMotifFinder | Signal |

---

### 21. Epigenetics Analysis (3 units)

#### EPIGEN-CPG-001: CpG Site Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.FindCpGSites(...)` |
| **Complexity** | O(n) |
| **Class** | EpigeneticsAnalyzer |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindCpGSites(sequence)` | EpigeneticsAnalyzer | Canonical |
| `FindCpGIslands(sequence, params)` | EpigeneticsAnalyzer | Islands |
| `CalculateObservedExpectedCpG(sequence)` | EpigeneticsAnalyzer | O/E ratio |

---

#### EPIGEN-METHYL-001: Methylation Analysis

| Поле | Значение |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.FindMethylationSites(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindMethylationSites(sequence)` | EpigeneticsAnalyzer | Canonical |
| `CalculateMethylationProfile(sites)` | EpigeneticsAnalyzer | Profile |
| `GetMethylationContext(site)` | EpigeneticsAnalyzer | Context (CpG/CHG/CHH) |

---

#### EPIGEN-DMR-001: Differentially Methylated Regions

| Поле | Значение |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.FindDMRs(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindDMRs(profile1, profile2, threshold)` | EpigeneticsAnalyzer | Canonical |
| `AnnotateDMRs(dmrs, annotations)` | EpigeneticsAnalyzer | Annotate |

---

### 22. Variant Calling (4 units)

#### VARIANT-CALL-001: Variant Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `VariantCaller.CallVariants(...)` |
| **Complexity** | O(n × m) |
| **Class** | VariantCaller |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CallVariants(reference, query)` | VariantCaller | Canonical |
| `CallVariantsFromAlignment(aligned1, aligned2)` | VariantCaller | From alignment |
| `ClassifyVariant(variant)` | VariantCaller | Classification |

---

#### VARIANT-SNP-001: SNP Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `VariantCaller.FindSnps(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindSnps(reference, query)` | VariantCaller | Canonical |
| `FindSnpsDirect(ref, query)` | VariantCaller | Direct (no alignment) |

---

#### VARIANT-INDEL-001: Indel Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `VariantCaller.FindInsertions/FindDeletions(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindInsertions(reference, query)` | VariantCaller | Insertions |
| `FindDeletions(reference, query)` | VariantCaller | Deletions |

---

#### VARIANT-ANNOT-001: Variant Annotation

| Поле | Значение |
|------|----------|
| **Canonical** | `VariantAnnotator.Annotate(...)` |
| **Complexity** | O(v × g) |
| **Class** | VariantAnnotator |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `Annotate(variants, annotations)` | VariantAnnotator | Canonical |
| `PredictFunctionalImpact(variant)` | VariantAnnotator | Impact |

---

### 23. Structural Variant Analysis (3 units)

#### SV-DETECT-001: SV Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `StructuralVariantAnalyzer.DetectSVs(...)` |
| **Complexity** | O(n log n) |
| **Class** | StructuralVariantAnalyzer |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `DetectSVs(readPairs, splitReads)` | StructuralVariantAnalyzer | Canonical |
| `FindDiscordantPairs(readPairs, params)` | StructuralVariantAnalyzer | Discordant |
| `ClassifySV(sv)` | StructuralVariantAnalyzer | Classification |

**SV Types:**
- [ ] Deletion, Duplication, Inversion
- [ ] Insertion, Translocation
- [ ] Complex rearrangements

---

#### SV-BREAKPOINT-001: Breakpoint Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `StructuralVariantAnalyzer.FindBreakpoints(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindBreakpoints(splitReads)` | StructuralVariantAnalyzer | Canonical |
| `RefineBreakpoint(region, reads)` | StructuralVariantAnalyzer | Refinement |

---

#### SV-CNV-001: Copy Number Variation

| Поле | Значение |
|------|----------|
| **Canonical** | `StructuralVariantAnalyzer.DetectCNV(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `DetectCNV(depthData, window)` | StructuralVariantAnalyzer | Canonical |
| `SegmentCopyNumber(logRatios)` | StructuralVariantAnalyzer | Segmentation |

---

### 24. Sequence Assembly (3 units)

#### ASSEMBLY-OLC-001: Overlap-Layout-Consensus

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceAssembler.AssembleOLC(...)` |
| **Complexity** | O(n² × m) |
| **Class** | SequenceAssembler |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `AssembleOLC(reads, params)` | SequenceAssembler | Canonical |
| `FindAllOverlaps(reads, minOverlap, minId)` | SequenceAssembler | Overlaps |

---

#### ASSEMBLY-DBG-001: De Bruijn Graph Assembly

| Поле | Значение |
|------|----------|
| **Canonical** | `SequenceAssembler.AssembleDeBruijn(...)` |
| **Complexity** | O(n × k) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `AssembleDeBruijn(reads, params)` | SequenceAssembler | Canonical |
| `BuildDeBruijnGraph(reads, k)` | SequenceAssembler | Graph construction |

---

#### ASSEMBLY-STATS-001: Assembly Statistics

| Поле | Значение |
|------|----------|
| **Canonical** | `GenomeAssemblyAnalyzer.CalculateStatistics(...)` |
| **Complexity** | O(n log n) |
| **Class** | GenomeAssemblyAnalyzer |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateStatistics(contigs)` | GenomeAssemblyAnalyzer | Canonical |
| `CalculateN50(contigs)` | GenomeAssemblyAnalyzer | N50 |
| `CalculateNx(contigs, threshold)` | GenomeAssemblyAnalyzer | Nx/Lx |
| `FindGaps(sequence)` | GenomeAssemblyAnalyzer | Gap detection |

---

### 25. Transcriptome Analysis (3 units)

#### TRANS-EXPR-001: Expression Quantification

| Поле | Значение |
|------|----------|
| **Canonical** | `TranscriptomeAnalyzer.CalculateTPM(...)` |
| **Complexity** | O(n) |
| **Class** | TranscriptomeAnalyzer |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateTPM(geneCounts)` | TranscriptomeAnalyzer | TPM |
| `CalculateFPKM(count, length, total)` | TranscriptomeAnalyzer | FPKM |
| `QuantileNormalize(samples)` | TranscriptomeAnalyzer | Normalization |

---

#### TRANS-DIFF-001: Differential Expression

| Поле | Значение |
|------|----------|
| **Canonical** | `TranscriptomeAnalyzer.FindDifferentiallyExpressed(...)` |
| **Complexity** | O(g × s) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindDifferentiallyExpressed(cond1, cond2, alpha)` | TranscriptomeAnalyzer | Canonical |
| `CalculateFoldChange(expr1, expr2)` | TranscriptomeAnalyzer | Fold change |

---

#### TRANS-SPLICE-001: Alternative Splicing

| Поле | Значение |
|------|----------|
| **Canonical** | `TranscriptomeAnalyzer.DetectAlternativeSplicing(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `DetectAlternativeSplicing(isoforms)` | TranscriptomeAnalyzer | Canonical |
| `CalculatePSI(event, reads)` | TranscriptomeAnalyzer | Percent spliced in |

---

### 26. Comparative Genomics (3 units)

#### COMPGEN-SYNTENY-001: Synteny Detection

| Поле | Значение |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindSyntenicBlocks(...)` |
| **Complexity** | O(n²) |
| **Class** | ComparativeGenomics |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindSyntenicBlocks(genes1, genes2, orthologs)` | ComparativeGenomics | Canonical |
| `VisualizeSynteny(blocks)` | ComparativeGenomics | Visualization |

---

#### COMPGEN-ORTHO-001: Ortholog Identification

| Поле | Значение |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindOrthologs(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `FindOrthologs(genes1, genes2, minIdentity)` | ComparativeGenomics | Canonical |
| `FindParalogs(genes, minIdentity)` | ComparativeGenomics | Paralogs |

---

#### COMPGEN-REARR-001: Genome Rearrangements

| Поле | Значение |
|------|----------|
| **Canonical** | `ComparativeGenomics.DetectRearrangements(...)` |
| **Complexity** | O(n log n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `DetectRearrangements(blocks)` | ComparativeGenomics | Canonical |
| `ClassifyRearrangement(event)` | ComparativeGenomics | Classification |

---

### 27. Pan-Genome Analysis (2 units)

#### PANGEN-CORE-001: Core/Accessory Genome

| Поле | Значение |
|------|----------|
| **Canonical** | `PanGenomeAnalyzer.ConstructPanGenome(...)` |
| **Complexity** | O(g² × s) |
| **Class** | PanGenomeAnalyzer |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `ConstructPanGenome(genomes, idThreshold)` | PanGenomeAnalyzer | Canonical |
| `IdentifyCoreGenes(clusters, threshold)` | PanGenomeAnalyzer | Core genes |

---

#### PANGEN-CLUSTER-001: Gene Clustering

| Поле | Значение |
|------|----------|
| **Canonical** | `PanGenomeAnalyzer.ClusterGenes(...)` |
| **Complexity** | O(g² × s) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `ClusterGenes(genomes, idThreshold)` | PanGenomeAnalyzer | Canonical |
| `GeneratePresenceAbsenceMatrix(clusters)` | PanGenomeAnalyzer | Matrix |

---

### 28. Quality Score Analysis (2 units)

#### QUALITY-PHRED-001: Phred Score Handling

| Поле | Значение |
|------|----------|
| **Canonical** | `QualityScoreAnalyzer.ParseQualityString(...)` |
| **Complexity** | O(n) |
| **Class** | QualityScoreAnalyzer |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `ParseQualityString(qualStr, encoding)` | QualityScoreAnalyzer | Canonical |
| `ToQualityString(scores, encoding)` | QualityScoreAnalyzer | Export |
| `ConvertEncoding(qualStr, from, to)` | QualityScoreAnalyzer | Convert |

---

#### QUALITY-STATS-001: Quality Statistics

| Поле | Значение |
|------|----------|
| **Canonical** | `QualityScoreAnalyzer.CalculateStatistics(...)` |
| **Complexity** | O(n) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `CalculateStatistics(scores)` | QualityScoreAnalyzer | Canonical |
| `CalculateQ30Percentage(scores)` | QualityScoreAnalyzer | Q30 |

---

### 29. Probe Design (2 units)

#### PROBE-DESIGN-001: Hybridization Probe Design

| Поле | Значение |
|------|----------|
| **Canonical** | `ProbeDesigner.DesignProbes(...)` |
| **Complexity** | O(n²) |
| **Class** | ProbeDesigner |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `DesignProbes(template, params)` | ProbeDesigner | Canonical |
| `DesignTilingProbes(template, overlap)` | ProbeDesigner | Tiling |
| `ScoreProbe(sequence, params)` | ProbeDesigner | Scoring |

---

#### PROBE-VALID-001: Probe Validation

| Поле | Значение |
|------|----------|
| **Canonical** | `ProbeDesigner.ValidateProbe(...)` |
| **Complexity** | O(n × g) |

**Methods:**
| Метод | Класс | Тип |
|-------|-------|-----|
| `ValidateProbe(probe, genome)` | ProbeDesigner | Canonical |
| `CheckSpecificity(probe, database)` | ProbeDesigner | Specificity |

---

## Appendix A: Method Index

| Method | Test Unit ID |
|--------|--------------|
| `AssembleDeBruijn` | ASSEMBLY-DBG-001 |
| `AssembleOLC` | ASSEMBLY-OLC-001 |
| `CalculateAlphaDiversity` | META-ALPHA-001 |
| `CalculateAlleleFrequencies` | POP-FREQ-001 |
| `CalculateBetaDiversity` | META-BETA-001 |
| `CalculateCAI` | CODON-CAI-001 |
| `CalculateCodonUsage` | CODON-USAGE-001 |
| `CalculateCumulativeGcSkew` | SEQ-GCSKEW-001 |
| `CalculateDistanceMatrix` | PHYLO-DIST-001 |
| `CalculateDiversityStatistics` | POP-DIV-001 |
| `CalculateFoldChange` | TRANS-DIFF-001 |
| `CalculateFPKM` | TRANS-EXPR-001 |
| `CalculateFreeEnergy` | RNA-ENERGY-001 |
| `CalculateFst` | POP-FST-001 |
| `CalculateGcContent` | SEQ-GC-001 |
| `CalculateGcFraction` | SEQ-GC-001 |
| `CalculateGcSkew` | SEQ-GCSKEW-001 |
| `CalculateKmerEntropy` | KMER-FREQ-001 |
| `CalculateLD` | POP-LD-001 |
| `CalculateLinguisticComplexity` | SEQ-COMPLEX-001 |
| `CalculateMeltingTemperature` | PRIMER-TM-001 |
| `CalculateMethylationProfile` | EPIGEN-METHYL-001 |
| `CalculateN50` | ASSEMBLY-STATS-001 |
| `CalculateNucleotideDiversity` | POP-DIV-001 |
| `CalculateObservedExpectedCpG` | EPIGEN-CPG-001 |
| `CalculatePSI` | TRANS-SPLICE-001 |
| `CalculateQ30Percentage` | QUALITY-STATS-001 |
| `CalculateShannonEntropy` | SEQ-ENTROPY-001 |
| `CalculateSpecificityScore` | CRISPR-OFF-001 |
| `CalculateStatistics` | ASSEMBLY-STATS-001, QUALITY-STATS-001 |
| `CalculateTajimasD` | POP-DIV-001 |
| `CalculateTPM` | TRANS-EXPR-001 |
| `CalculateWattersonTheta` | POP-DIV-001 |
| `CalculateWindowedGcSkew` | SEQ-GCSKEW-001 |
| `CallVariants` | VARIANT-CALL-001 |
| `CallVariantsFromAlignment` | VARIANT-CALL-001 |
| `CheckSpecificity` | PROBE-VALID-001 |
| `ClassifyReads` | META-CLASS-001 |
| `ClusterGenes` | PANGEN-CLUSTER-001 |
| `Complement` | SEQ-COMP-001 |
| `ConstructPanGenome` | PANGEN-CORE-001 |
| `Contains` | PAT-EXACT-001 |
| `ConvertEncoding` | QUALITY-PHRED-001 |
| `CountKmers` | KMER-COUNT-001 |
| `CreateMiRna` | MIRNA-SEED-001 |
| `CreatePwm` | PAT-PWM-001 |
| `DesignGuideRnas` | CRISPR-GUIDE-001 |
| `DesignPrimers` | PRIMER-DESIGN-001 |
| `DesignProbes` | PROBE-DESIGN-001 |
| `DesignTilingProbes` | PROBE-DESIGN-001 |
| `DetectAlternativeSplicing` | TRANS-SPLICE-001 |
| `DetectAneuploidy` | CHROM-ANEU-001 |
| `DetectCNV` | SV-CNV-001 |
| `DetectPloidy` | CHROM-KARYO-001 |
| `DetectRearrangements` | COMPGEN-REARR-001 |
| `DetectSVs` | SV-DETECT-001 |
| `Digest` | RESTR-DIGEST-001 |
| `EditDistance` | PAT-APPROX-002 |
| `EvaluateGuideRna` | CRISPR-GUIDE-001 |
| `EvaluatePrimer` | PRIMER-DESIGN-001 |
| `FilterByQuality` | PARSE-FASTQ-001 |
| `FindAcceptorSites` | SPLICE-ACCEPTOR-001 |
| `FindAllOccurrences` | PAT-EXACT-001 |
| `FindAllOverlaps` | ASSEMBLY-OLC-001 |
| `FindBreakpoints` | SV-BREAKPOINT-001 |
| `FindClumps` | KMER-FIND-001 |
| `FindCpGIslands` | EPIGEN-CPG-001 |
| `FindCpGSites` | EPIGEN-CPG-001 |
| `FindDegenerateMotif` | PAT-IUPAC-001 |
| `FindDifferentiallyExpressed` | TRANS-DIFF-001 |
| `FindDirectRepeats` | REP-DIRECT-001 |
| `FindDiscordantPairs` | SV-DETECT-001 |
| `FindDMRs` | EPIGEN-DMR-001 |
| `FindDonorSites` | SPLICE-DONOR-001 |
| `FindExactMotif` | PAT-EXACT-001 |
| `FindExons` | SPLICE-PREDICT-001 |
| `FindHaplotypeBlocks` | POP-LD-001 |
| `FindIntrons` | SPLICE-PREDICT-001 |
| `FindInvertedRepeats` | REP-INV-001 |
| `FindLowComplexityRegions` | SEQ-COMPLEX-001 |
| `FindMethylationSites` | EPIGEN-METHYL-001 |
| `FindMicrosatellites` | REP-STR-001 |
| `FindMostFrequentKmers` | KMER-FIND-001 |
| `FindMotif` | PAT-EXACT-001 |
| `FindMotifs` | PROTMOTIF-FIND-001 |
| `FindOffTargets` | CRISPR-OFF-001 |
| `FindOrfs` | ANNOT-ORF-001 |
| `FindOriginOfReplication` | SEQ-GCSKEW-001 |
| `FindOrthologs` | COMPGEN-ORTHO-001 |
| `FindPalindromes` | REP-PALIN-001 |
| `FindPamSites` | CRISPR-PAM-001 |
| `FindParalogs` | COMPGEN-ORTHO-001 |
| `FindPreMiRnas` | MIRNA-PRECURSOR-001 |
| `FindPromoterMotifs` | ANNOT-PROM-001 |
| `FindRareCodons` | CODON-RARE-001 |
| `FindRibosomeBindingSites` | ANNOT-GENE-001 |
| `FindSites` | RESTR-FIND-001 |
| `FindStemLoops` | RNA-STEMLOOP-001 |
| `FindSyntenicBlocks` | COMPGEN-SYNTENY-001 |
| `FindSyntenyBlocks` | CHROM-SYNT-001 |
| `FindTandemRepeats` | REP-TANDEM-001 |
| `FindTargetSites` | MIRNA-TARGET-001 |
| `FindUniqueKmers` | KMER-FIND-001 |
| `FindWithEdits` | PAT-APPROX-002 |
| `FindWithMismatches` | PAT-APPROX-001 |
| `GeneratePresenceAbsenceMatrix` | PANGEN-CLUSTER-001 |
| `GenerateTaxonomicProfile` | META-PROF-001 |
| `GetComplementBase` | SEQ-COMP-001 |
| `GetKmerFrequencies` | KMER-FREQ-001 |
| `GetKmerSpectrum` | KMER-FREQ-001 |
| `GetSeedSequence` | MIRNA-SEED-001 |
| `GlobalAlign` | ALIGN-GLOBAL-001 |
| `HammingDistance` | PAT-APPROX-001 |
| `HasHairpinPotential` | PRIMER-STRUCT-001 |
| `HasPrimerDimer` | PRIMER-STRUCT-001 |
| `IdentifyCoreGenes` | PANGEN-CORE-001 |
| `IdentifyDisorderedRegions` | DISORDER-REGION-001 |
| `IsStartCodon` | TRANS-CODON-001 |
| `IsStopCodon` | TRANS-CODON-001 |
| `IsValidDna` | SEQ-VALID-001 |
| `IsValidRna` | SEQ-VALID-001 |
| `LocalAlign` | ALIGN-LOCAL-001 |
| `MaskLowComplexity` | SEQ-COMPLEX-001 |
| `MatchesIupac` | PAT-IUPAC-001 |
| `MatchPrositePattern` | PROTMOTIF-PROSITE-001 |
| `MergeOverlapping` | PARSE-BED-001 |
| `MultipleAlign` | ALIGN-MULTI-001 |
| `OptimizeSequence` | CODON-OPT-001 |
| `Parse` | PARSE-FASTA-001, PARSE-FASTQ-001, PARSE-BED-001, PARSE-VCF-001 |
| `ParseFile` | PARSE-FASTA-001, PARSE-FASTQ-001, PARSE-BED-001, PARSE-VCF-001 |
| `ParseGff3` | ANNOT-GFF-001 |
| `ParseNewick` | PHYLO-NEWICK-001 |
| `ParsePrositePattern` | PROTMOTIF-PROSITE-001 |
| `ParseQualityString` | QUALITY-PHRED-001 |
| `PatristicDistance` | PHYLO-COMP-001 |
| `Predict` | RNA-STRUCT-001 |
| `PredictDisorder` | DISORDER-PRED-001 |
| `PredictDomains` | PROTMOTIF-DOMAIN-001 |
| `PredictGenes` | ANNOT-GENE-001 |
| `PredictGeneStructure` | SPLICE-PREDICT-001 |
| `PredictSignalPeptide` | PROTMOTIF-DOMAIN-001 |
| `QuantileNormalize` | TRANS-EXPR-001 |
| `ReverseComplement` | SEQ-REVCOMP-001 |
| `RobinsonFouldsDistance` | PHYLO-COMP-001 |
| `ScanForPattern` | PROTMOTIF-FIND-001 |
| `ScanWithPwm` | PAT-PWM-001 |
| `ScoreDonorSite` | SPLICE-DONOR-001 |
| `ScoreAcceptorSite` | SPLICE-ACCEPTOR-001 |
| `ScoreProbe` | PROBE-DESIGN-001 |
| `ScoreTargetSite` | MIRNA-TARGET-001 |
| `SegmentCopyNumber` | SV-CNV-001 |
| `SemiGlobalAlign` | ALIGN-SEMI-001 |
| `TestHardyWeinberg` | POP-HW-001 |
| `ToDotBracket` | RNA-STRUCT-001 |
| `ToFasta` | PARSE-FASTA-001 |
| `ToGff3` | ANNOT-GFF-001 |
| `ToNewick` | PHYLO-NEWICK-001 |
| `ToQualityString` | QUALITY-PHRED-001 |
| `Translate` | TRANS-CODON-001, TRANS-PROT-001 |
| `TryGetComplement` | SEQ-COMP-001 |
| `TryGetReverseComplement` | SEQ-REVCOMP-001 |
| `ValidateProbe` | PROBE-VALID-001 |

---

## Appendix B: Complexity Notes

| Test Unit | Claimed | Verified | Notes |
|-----------|---------|----------|-------|
| REP-STR-001 | O(n²) | ⚠️ | Actually O(n × U × R), depends on parameters |
| REP-INV-001 | O(n²) | ⚠️ | O(n² × L) with maxLoopLength |
| CHROM-SYNT-001 | O(n log n) | ⚠️ | Has nested loops, needs verification |
| ALIGN-MULTI-001 | O(n² × m²) | ⚠️ | Progressive is typically O(n² × m) |
| META-BIN-001 | O(n) | ⚠️ | K-means is O(n × k × i), verify iterations |
| CRISPR-OFF-001 | O(n × m) | ⚠️ | May be exponential with high mismatches |
| RNA-STRUCT-001 | O(n³) | ✓ | Standard Nussinov/Zuker |
| ASSEMBLY-OLC-001 | O(n² × m) | ⚠️ | Depends on overlap detection method |
| PANGEN-CORE-001 | O(g² × s) | ⚠️ | All-vs-all comparison |
| SV-CNV-001 | O(n) | ⚠️ | Segmentation may be O(n log n) |

---

## Appendix C: Canonical Implementations

### GC Content (SEQ-GC-001)

```
Canonical: SequenceExtensions.CalculateGcContent(ReadOnlySpan<char>)
         ↑
    ┌────┴────┐
    │         │
CalculateGcContentFast(string)   CalculateGcFraction(ReadOnlySpan)
    │                                      │
    │                              CalculateGcFractionFast(string)
    │                                      │
    ├──────────────────────────────────────┤
    │                                      │
PrimerDesigner.CalculateGcContent   MetagenomicsAnalyzer.CalculateGcContent
ChromosomeAnalyzer (internal)       DnaSequence.GcContent
```

**Test Strategy:** Test canonical `CalculateGcContent(ReadOnlySpan)` thoroughly. Delegates need only smoke tests verifying they call canonical correctly.

---

### DNA Complement (SEQ-COMP-001 + SEQ-REVCOMP-001)

```
Canonical: SequenceExtensions.GetComplementBase(char)
         ↑
    ┌────┴────────────────────┐
    │                         │
TryGetComplement(Span)   TryGetReverseComplement(Span)
    │                         │
DnaSequence.Complement   DnaSequence.ReverseComplement
                              │
                    GetReverseComplementString(string)
```

**Test Strategy:** Test `GetComplementBase` for all nucleotides. Test Span methods for sequence operations. Instance methods are wrappers.

---

## Appendix D: Class Coverage Summary

| Class | Test Units | Status |
|-------|------------|--------|
| ApproximateMatcher | PAT-APPROX-001, PAT-APPROX-002 | ✓ |
| BedParser | PARSE-BED-001 | ✓ |
| ChromosomeAnalyzer | CHROM-TELO-001 to CHROM-SYNT-001 | ✓ |
| CodonOptimizer | CODON-OPT-001 to CODON-USAGE-001 | ✓ |
| ComparativeGenomics | COMPGEN-SYNTENY-001 to COMPGEN-REARR-001 | ✓ |
| CrisprDesigner | CRISPR-PAM-001 to CRISPR-OFF-001 | ✓ |
| DisorderPredictor | DISORDER-PRED-001, DISORDER-REGION-001 | ✓ |
| EmblParser | PARSE-EMBL-001 | ✓ |
| EpigeneticsAnalyzer | EPIGEN-CPG-001 to EPIGEN-DMR-001 | ✓ |
| FastaParser | PARSE-FASTA-001 | ✓ |
| FastqParser | PARSE-FASTQ-001 | ✓ |
| GcSkewCalculator | SEQ-GCSKEW-001 | ✓ |
| GenBankParser | PARSE-GENBANK-001 | ✓ |
| GenomeAnnotator | ANNOT-ORF-001 to ANNOT-GFF-001 | ✓ |
| GenomeAssemblyAnalyzer | ASSEMBLY-STATS-001 | ✓ |
| GenomicAnalyzer | REP-TANDEM-001, REP-PALIN-001 | ✓ |
| GffParser | PARSE-GFF-001 | ✓ |
| GeneticCode | TRANS-CODON-001 | ✓ |
| IupacHelper | PAT-IUPAC-001 | ✓ |
| KmerAnalyzer | KMER-COUNT-001 to KMER-FIND-001 | ✓ |
| MetagenomicsAnalyzer | META-CLASS-001 to META-BIN-001 | ✓ |
| MiRnaAnalyzer | MIRNA-SEED-001 to MIRNA-PRECURSOR-001 | ✓ |
| MotifFinder | PAT-PWM-001, PAT-IUPAC-001 | ✓ |
| PanGenomeAnalyzer | PANGEN-CORE-001, PANGEN-CLUSTER-001 | ✓ |
| PhylogeneticAnalyzer | PHYLO-DIST-001 to PHYLO-COMP-001 | ✓ |
| PopulationGeneticsAnalyzer | POP-FREQ-001 to POP-LD-001 | ✓ |
| PrimerDesigner | PRIMER-TM-001 to PRIMER-STRUCT-001 | ✓ |
| ProbeDesigner | PROBE-DESIGN-001, PROBE-VALID-001 | ✓ |
| ProteinMotifFinder | PROTMOTIF-FIND-001 to PROTMOTIF-DOMAIN-001 | ✓ |
| QualityScoreAnalyzer | QUALITY-PHRED-001, QUALITY-STATS-001 | ✓ |
| RepeatFinder | REP-STR-001 to REP-PALIN-001 | ✓ |
| RestrictionAnalyzer | RESTR-FIND-001, RESTR-DIGEST-001 | ✓ |
| RnaSecondaryStructure | RNA-STRUCT-001 to RNA-ENERGY-001 | ✓ |
| SequenceAligner | ALIGN-GLOBAL-001 to ALIGN-MULTI-001 | ✓ |
| SequenceAssembler | ASSEMBLY-OLC-001, ASSEMBLY-DBG-001 | ✓ |
| SequenceComplexity | SEQ-COMPLEX-001, SEQ-ENTROPY-001 | ✓ |
| SequenceExtensions | SEQ-GC-001 to SEQ-VALID-001 | ✓ |
| SpliceSitePredictor | SPLICE-DONOR-001 to SPLICE-PREDICT-001 | ✓ |
| StructuralVariantAnalyzer | SV-DETECT-001 to SV-CNV-001 | ✓ |
| TranscriptomeAnalyzer | TRANS-EXPR-001 to TRANS-SPLICE-001 | ✓ |
| Translator | TRANS-PROT-001 | ✓ |
| VariantAnnotator | VARIANT-ANNOT-001 | ✓ |
| VariantCaller | VARIANT-CALL-001 to VARIANT-INDEL-001 | ✓ |
| VcfParser | PARSE-VCF-001 | ✓ |

**Total Classes Covered: 42/42 (100%)**

---

*Сгенерировано: 2026-01-22*
*Версия чек-листа: 2.1 (100% Coverage)*
