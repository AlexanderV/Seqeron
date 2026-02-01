# Algorithm: Guide RNA Design

## Algorithm Group
**MolTools** (Molecular Biology Tools)

## Test Unit ID
CRISPR-GUIDE-001

## Overview
Guide RNA (gRNA/sgRNA) design algorithms for CRISPR-Cas genome editing. These methods evaluate and design optimal single guide RNAs based on sequence composition metrics that correlate with editing efficiency and specificity.

## Canonical Methods
- `CrisprDesigner.DesignGuideRnas(DnaSequence, int, int, CrisprSystemType, GuideRnaParameters?)` - Design guide RNAs in a target region
- `CrisprDesigner.EvaluateGuideRna(string, CrisprSystemType, GuideRnaParameters?)` - Evaluate a guide RNA sequence for quality

## Theoretical Foundation

### Guide RNA Structure
A single guide RNA (sgRNA) consists of:
1. **Spacer Sequence**: ~20 nucleotides at the 5' end that determine targeting specificity
2. **Scaffold Sequence**: Constant region required for Cas9 binding (tracrRNA-derived)

**Source**: Addgene CRISPR Guide - "The gRNA is a short synthetic RNA composed of a scaffold sequence necessary for Cas-binding and a user-defined ∼20-nucleotide spacer that defines the genomic target to be modified."

### Key Quality Metrics

#### 1. GC Content (40-70% optimal)
Higher GC content enhances RNA-DNA duplex stability and reduces off-target hybridization. However, extremely high GC can cause secondary structures.

**Source**: Wikipedia Guide RNA - "The GC content of sgRNA should optimally be over 50% to improve the efficiency of targeting"

**Implementation**: Scoring penalizes guides outside 40-70% GC range.

#### 2. Poly-T Sequences (TTTT/UUUU)
Poly-T stretches of 4+ bases act as Pol III termination signals and must be avoided.

**Source**: Addgene - RNA polymerase III terminates at poly-T sequences, prematurely ending gRNA transcription.

**Implementation**: Detection of TTTT results in score penalty and issue flagging.

#### 3. Seed Region GC Content
The seed sequence (last 8-12 nucleotides at 3' end) initiates target annealing. Seed region GC affects binding initiation.

**Source**: Addgene CRISPR Guide - "the seed sequence (8–10 bases at the 3′ end of the gRNA targeting sequence) will begin to anneal to the target DNA"

**Implementation**: Calculates GC% of last 12 nucleotides separately.

#### 4. Self-Complementarity
Internal secondary structures reduce guide efficacy by sequestering the spacer sequence.

**Implementation**: Checks for potential stem-loop formation within the guide.

#### 5. Restriction Sites (Optional)
Common restriction sites in the guide can complicate cloning workflows.

**Implementation**: Optional detection of common 6-cutter recognition sequences.

### PAM Requirement
Guides are designed upstream of PAM sites. The DesignGuideRnas method identifies PAM sites in the target region and extracts the 20bp upstream sequence as guide candidates.

**Source**: Addgene - "The genomic target of the gRNA can be any ~20 nucleotide sequence, provided it meets two conditions: The sequence is unique compared to the rest of the genome [and] The target is present immediately adjacent to a Protospacer Adjacent Motif (PAM)"

## Algorithm Complexity
- `DesignGuideRnas`: O(n) where n is the target region length
- `EvaluateGuideRna`: O(k²) where k is guide length (due to self-complementarity check)

## Input Parameters

### EvaluateGuideRna
| Parameter | Type | Description |
|-----------|------|-------------|
| guide | string | 20bp guide sequence to evaluate |
| systemType | CrisprSystemType | CRISPR system (SpCas9, SaCas9, etc.) |
| parameters | GuideRnaParameters? | Optional custom parameters |

### DesignGuideRnas
| Parameter | Type | Description |
|-----------|------|-------------|
| sequence | DnaSequence | Target DNA sequence |
| regionStart | int | Start of target region (0-indexed) |
| regionEnd | int | End of target region (0-indexed, exclusive) |
| systemType | CrisprSystemType | CRISPR system type |
| parameters | GuideRnaParameters? | Optional custom parameters |

### GuideRnaParameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| MinGcContent | 40 | Minimum acceptable GC% |
| MaxGcContent | 70 | Maximum acceptable GC% |
| MinScore | 50 | Minimum score threshold |
| AvoidPolyT | true | Penalize TTTT sequences |
| CheckSelfComplementarity | true | Check for internal structures |

## Output

### GuideRnaCandidate Record
| Field | Type | Description |
|-------|------|-------------|
| Sequence | string | 20bp guide sequence |
| Position | int | Position in source sequence |
| GcContent | double | Overall GC percentage |
| SeedGcContent | double | GC% of seed region (last 12bp) |
| HasPolyT | bool | Contains TTTT |
| Score | double | Composite quality score (0-100) |
| Issues | IReadOnlyList<string> | Quality concerns |
| FullGuideRna | string | Guide + scaffold sequence |

## Scoring Algorithm

### Base Score: 100 points

### Deductions:
1. **GC Content Outside Range**:
   - Below MinGcContent: -20 points
   - Above MaxGcContent: -20 points

2. **Poly-T Detection**:
   - Contains TTTT: -15 points

3. **Self-Complementarity**:
   - High complementarity detected: -10 points

4. **Seed Region GC**:
   - Seed GC outside 30-80%: -5 points

5. **Restriction Sites**:
   - Each detected site: -5 points (max -15)

## Edge Cases

### Documented in Evidence
1. **Empty guide sequence**: Throws ArgumentNullException
2. **Null sequence for DesignGuideRnas**: Throws ArgumentNullException
3. **Region start < 0**: Throws ArgumentOutOfRangeException
4. **Region end > sequence length**: Throws ArgumentOutOfRangeException
5. **Guide length ≠ 20**: May produce unexpected results (ASSUMPTION - verify standard length handling)

### Biological Edge Cases
1. **All-A guide (0% GC)**: Valid but low score, issues reported
2. **All-G/C guide (100% GC)**: Valid but low score, issues reported
3. **Guide with TTTT in middle**: Detected, penalized
4. **Guide with TTTT at boundary**: Detected if 4+ consecutive T's

## Evidence Sources

### Primary Sources
1. **Addgene CRISPR Guide** (https://www.addgene.org/guides/crispr/)
   - Authoritative resource from nonprofit plasmid repository
   - Details on gRNA structure, PAM requirements, seed sequence

2. **Wikipedia: Guide RNA** (https://en.wikipedia.org/wiki/Guide_RNA)
   - GC content >50% optimal
   - Length typically 17-24bp with 20bp standard

3. **Wikipedia: CRISPR gene editing** (https://en.wikipedia.org/wiki/CRISPR_gene_editing)
   - SpCas9 PAM is 5'-NGG-3'
   - Target sequence is 20 bases long

4. **Wikipedia: Protospacer adjacent motif** (https://en.wikipedia.org/wiki/Protospacer_adjacent_motif)
   - PAM sequences for different Cas proteins
   - PAM is required for Cas9 cleavage

### Academic References (from Addgene)
- Doench et al. (2014) "Rational design of highly active sgRNAs" - Nature Biotechnology
- Hsu et al. (2013) "DNA targeting specificity of RNA-guided Cas9 nucleases" - Nature Biotechnology

## Implementation Notes

### Current Implementation Accuracy
The implementation aligns well with documented guidelines:
- ✅ 40-70% GC range matches commonly used thresholds
- ✅ Poly-T detection (TTTT) correctly identifies termination signals
- ✅ Seed region (12bp) slightly larger than minimum (8-10bp), more conservative
- ✅ 20bp standard guide length

### Known Gaps
- Position-weighted scoring (Doench rules) not implemented
- Machine learning-based prediction not available
- No consideration of chromatin accessibility

## Related Test Units
- **CRISPR-PAM-001**: PAM Site Detection (prerequisite)
- **CRISPR-OFF-001**: Off-Target Analysis (related)

## Change History
| Date | Version | Changes |
|------|---------|---------|
| 2025-01-08 | 1.0 | Initial documentation for CRISPR-GUIDE-001 |
