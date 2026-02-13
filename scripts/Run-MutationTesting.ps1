<#
.SYNOPSIS
    Runs Stryker mutation testing per project, targeting only files with processed algorithms.
.DESCRIPTION
    Each project run targets specific source files from ALGORITHMS_CHECKLIST_V2.md (☑ items only).
    Heavy algorithms (O(n²)+) are excluded to prevent hangs.
    Results are collected in a summary at the end.
.NOTES
    Excluded heavy algorithms:
      - SequenceAligner.cs       (976 LOC)  — NW O(n×m), SW O(n×m), MSA O(n²×k)
      - ApproximateMatcher.cs    (415 LOC)  — EditDistance O(n×m)
      - RnaSecondaryStructure.cs (845 LOC)  — Nussinov O(n³), MFE O(n³)
      - PhylogeneticAnalyzer.cs  (768 LOC)  — Distance matrix O(n²), UPGMA/NJ
      - ProteinMotifFinder.cs   (1022 LOC)  — Domain prediction heavy scans
      - ChromosomeAnalyzer.cs    (889 LOC)  — Synteny analysis, karyotype
      - MetagenomicsAnalyzer.cs  (701 LOC)  — Genome binning with permutations
#>
param(
    [string]$TestProjectDir = "d:\Prototype\SuffixTree\tests\Seqeron\Seqeron.Genomics.Tests",
    [switch]$DryRun
)

$ErrorActionPreference = 'Continue'

# Define projects and their processed source files (heavy algorithms excluded)
$projects = @(
    # 1. Core — remaining processed files (SequenceExtensions already done — 100%)
    @{
        Project = "Seqeron.Genomics.Core.csproj"
        Files = @(
            "GeneticCode.cs",
            "Translator.cs",
            "DnaSequence.cs",
            "RnaSequence.cs",
            "ProteinSequence.cs",
            "IupacHelper.cs"
        )
    },

    # 2. Analysis — EXCLUDED: RnaSecondaryStructure.cs (O(n³)), ProteinMotifFinder.cs (heavy scans)
    @{
        Project = "Seqeron.Genomics.Analysis.csproj"
        Files = @(
            "GcSkewCalculator.cs",
            "GenomicAnalyzer.cs",
            "SequenceComplexity.cs",
            "KmerAnalyzer.cs",
            "MotifFinder.cs",
            "RepeatFinder.cs",
            "DisorderPredictor.cs"
        )
    },

    # 3. MolTools — all lightweight
    @{
        Project = "Seqeron.Genomics.MolTools.csproj"
        Files = @(
            "CrisprDesigner.cs",
            "PrimerDesigner.cs",
            "ProbeDesigner.cs",
            "RestrictionAnalyzer.cs",
            "CodonOptimizer.cs",
            "CodonUsageAnalyzer.cs"
        )
    },

    # 4. Annotation — all processed files are manageable
    @{
        Project = "Seqeron.Genomics.Annotation.csproj"
        Files = @(
            "GenomeAnnotator.cs",
            "EpigeneticsAnalyzer.cs",
            "MiRnaAnalyzer.cs",
            "SpliceSitePredictor.cs"
        )
    },

    # 5. Alignment — EXCLUDED ENTIRELY: SequenceAligner.cs O(n×m), ApproximateMatcher.cs O(n×m)

    # 6. Phylogenetics — EXCLUDED ENTIRELY: PhylogeneticAnalyzer.cs O(n²)

    # 7. Population — PopulationGeneticsAnalyzer (lightweight stats)
    @{
        Project = "Seqeron.Genomics.Population.csproj"
        Files = @("PopulationGeneticsAnalyzer.cs")
    },

    # 8. Chromosome — EXCLUDED ENTIRELY: ChromosomeAnalyzer.cs (synteny heavy)

    # 9. Metagenomics — EXCLUDED ENTIRELY: MetagenomicsAnalyzer.cs (permutations)

    # 10. IO — All parsers (linear I/O, lightweight)
    @{
        Project = "Seqeron.Genomics.IO.csproj"
        Files = @(
            "FastaParser.cs",
            "FastqParser.cs",
            "BedParser.cs",
            "VcfParser.cs",
            "GffParser.cs",
            "GenBankParser.cs",
            "EmblParser.cs",
            "SequenceIO.cs"
        )
    }
)

$results = @()
$startTime = Get-Date

foreach ($p in $projects) {
    $projectName = $p.Project -replace '\.csproj$', ''
    $mutateArgs = ($p.Files | ForEach-Object { "--mutate `"**/$_`"" }) -join ' '
    
    Write-Host "`n========================================" -ForegroundColor Cyan
    Write-Host "  PROJECT: $projectName" -ForegroundColor Cyan
    Write-Host "  Files: $($p.Files -join ', ')" -ForegroundColor DarkCyan
    Write-Host "========================================" -ForegroundColor Cyan

    $cmd = "dotnet stryker --project `"$($p.Project)`" $mutateArgs -r json -r cleartext"
    
    if ($DryRun) {
        Write-Host "[DRY RUN] $cmd" -ForegroundColor Yellow
        continue
    }

    Push-Location $TestProjectDir
    try {
        $output = Invoke-Expression $cmd 2>&1 | Out-String
        Write-Host $output

        # Extract mutation score from output
        if ($output -match 'final mutation score is (\d+\.\d+) %') {
            $score = [double]$Matches[1]
        } elseif ($output -match 'unable to calculate') {
            $score = -1
        } else {
            $score = -2
        }

        $results += [PSCustomObject]@{
            Project  = $projectName
            Score    = if ($score -ge 0) { "$score%" } else { "N/A" }
            RawScore = $score
        }
    }
    finally {
        Pop-Location
    }
}

# Summary
$elapsed = (Get-Date) - $startTime
Write-Host "`n`n========================================" -ForegroundColor Green
Write-Host "  MUTATION TESTING SUMMARY" -ForegroundColor Green
Write-Host "  Elapsed: $($elapsed.ToString('hh\:mm\:ss'))" -ForegroundColor Green
Write-Host "========================================`n" -ForegroundColor Green
$results | Format-Table -AutoSize
