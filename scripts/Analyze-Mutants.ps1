#!/usr/bin/env pwsh
# Analyze survived mutants from Stryker mutation testing reports

$ErrorActionPreference = 'Stop'
$root = "d:\Prototype\SuffixTree\tests\Seqeron\Seqeron.Genomics.Tests\StrykerOutput"

$reports = Get-ChildItem -Path $root -Recurse -Filter "mutation-report.json" |
    Sort-Object FullName |
    Select-Object -Last 6 -ExpandProperty FullName

Write-Host "Found $($reports.Count) batch reports" -ForegroundColor Cyan

# Parse all reports
$survived = @()
$scores = @()

foreach ($rp in $reports) {
    $json = Get-Content $rp -Raw | ConvertFrom-Json
    foreach ($prop in $json.files.PSObject.Properties) {
        $fileName = ($prop.Name -split '[\\/]' | Select-Object -Last 1)
        $mutants = $prop.Value.mutants
        $killed = @($mutants | Where-Object { $_.status -eq 'Killed' }).Count
        $surv = @($mutants | Where-Object { $_.status -eq 'Survived' }).Count
        $noCov = @($mutants | Where-Object { $_.status -eq 'NoCoverage' }).Count
        $total = $mutants.Count
        $score = if (($killed + $surv) -gt 0) { [math]::Round($killed / ($killed + $surv) * 100, 1) } else { 100 }

        $scores += [PSCustomObject]@{
            FileName = $fileName
            Score    = $score
            Killed   = $killed
            Survived = $surv
            NoCov    = $noCov
            Total    = $total
        }

        foreach ($m in @($mutants | Where-Object { $_.status -eq 'Survived' })) {
            $survived += [PSCustomObject]@{
                FileName    = $fileName
                LineNum     = $m.location.start.line
                MutatorName = $m.mutatorName
                Repl        = if ($m.replacement) { $m.replacement.ToString() } else { '' }
            }
        }
    }
}

Write-Host "`nTotal survived mutants: $($survived.Count)" -ForegroundColor Yellow
Write-Host "Total files with mutations: $($scores.Count)" -ForegroundColor Yellow

# === SECTION 1: Scores by file ===
Write-Host "`n========================================" -ForegroundColor Green
Write-Host "MUTATION SCORES BY FILE (ascending)" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
$scores | Where-Object { $_.Survived -gt 0 -or $_.Killed -gt 0 } | Sort-Object Score | Format-Table -AutoSize

# === SECTION 2: By mutator type ===
Write-Host "`n========================================" -ForegroundColor Green
Write-Host "SURVIVED BY MUTATOR TYPE" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
$survived | Group-Object MutatorName | Sort-Object Count -Descending | Select-Object Count, Name | Format-Table -AutoSize

# === SECTION 3: Top-10 worst files - detailed hot lines ===
Write-Host "`n========================================" -ForegroundColor Green
Write-Host "TOP-10 WORST FILES - HOT LINE ANALYSIS" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green

$top10 = $survived | Group-Object FileName | Sort-Object Count -Descending | Select-Object -First 10
foreach ($grp in $top10) {
    $fileScore = ($scores | Where-Object { $_.FileName -eq $grp.Name }).Score
    Write-Host "`n--- $($grp.Name) ($($grp.Count) survived, score=$fileScore%) ---" -ForegroundColor Yellow
    
    # Mutator breakdown
    $grp.Group | Group-Object MutatorName | Sort-Object Count -Descending | ForEach-Object {
        Write-Host "  $($_.Count.ToString().PadLeft(4)) $($_.Name)"
    }
    
    # Top hot lines
    Write-Host "  Hot lines:" -ForegroundColor DarkCyan
    $hotLines = $grp.Group | Group-Object LineNum | Sort-Object Count -Descending | Select-Object -First 8
    foreach ($hl in $hotLines) {
        $sample = $hl.Group[0]
        $replSnippet = if ($sample.Repl.Length -gt 75) { $sample.Repl.Substring(0, 75) + "..." } else { $sample.Repl }
        Write-Host "    L$($hl.Name) ($($hl.Count)x): $replSnippet"
    }
}

# === SECTION 4: Pattern classification ===
Write-Host "`n========================================" -ForegroundColor Green
Write-Host "WEAKNESS PATTERN CLASSIFICATION" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green

$eqCount = @($survived | Where-Object { $_.MutatorName -eq 'Equality mutation' }).Count
$arithCount = @($survived | Where-Object { $_.MutatorName -eq 'Arithmetic mutation' }).Count
$logicalCount = @($survived | Where-Object { $_.MutatorName -eq 'Logical mutation' }).Count
$blockCount = @($survived | Where-Object { $_.MutatorName -eq 'Block removal mutation' }).Count
$nullCount = @($survived | Where-Object { $_.MutatorName -match 'Null coalescing' }).Count
$bitCount = @($survived | Where-Object { $_.MutatorName -eq 'Bitwise mutation' }).Count

Write-Host @"

Pattern 1: BOUNDARY CONDITIONS ($eqCount mutations, 54%)
  - < vs <=, > vs >= swaps survive  
  - Tests do not exercise exact boundary values
  - Fix: Add boundary-value tests (off-by-one: n-1, n, n+1)

Pattern 2: ARITHMETIC PRECISION ($arithCount mutations, 25%)
  - +/-, */divide swaps survive in formulas
  - Tests use approximate assertions or skip edge cases
  - Fix: Assert exact values with known inputs, test with 0/1/negative

Pattern 3: BOOLEAN LOGIC ($logicalCount mutations, 12%)
  - && to || and ! negation survive
  - Compound conditions not fully branch-tested
  - Fix: Test each clause independently (true/false matrix)

Pattern 4: DEAD CODE / BLOCK REMOVAL ($blockCount mutations, 6%)
  - Entire if-blocks can be removed without failing tests
  - Missing test coverage for specific code paths
  - Fix: Add tests that require the removed block's effect

Pattern 5: NULL HANDLING ($nullCount mutations, 3%)
  - ?? operator mutations survive
  - Tests never pass null arguments
  - Fix: Add null-input tests

Pattern 6: BITWISE ($bitCount mutations, <1%)
  - & vs | bit operator swaps
  - Affects SequenceIO.cs quality encoding
  - Fix: Add tests with specific bit patterns
"@

# === SECTION 5: Priority list for mutation-killing tests ===
Write-Host "`n========================================" -ForegroundColor Green
Write-Host "PRIORITY FILES FOR MUTATION-KILLING TESTS" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green

$priority = $scores | Where-Object { $_.Survived -gt 0 } | Sort-Object Score |
    Select-Object FileName, Score, Survived, @{N='Priority';E={
        if ($_.Score -lt 50) { 'CRITICAL' }
        elseif ($_.Score -lt 60) { 'HIGH' }
        elseif ($_.Score -lt 70) { 'MEDIUM' }
        else { 'LOW' }
    }}
$priority | Format-Table -AutoSize

Write-Host "CRITICAL: $(@($priority | Where-Object { $_.Priority -eq 'CRITICAL' }).Count) files" -ForegroundColor Red
Write-Host "HIGH:     $(@($priority | Where-Object { $_.Priority -eq 'HIGH' }).Count) files" -ForegroundColor DarkYellow
Write-Host "MEDIUM:   $(@($priority | Where-Object { $_.Priority -eq 'MEDIUM' }).Count) files" -ForegroundColor Yellow
Write-Host "LOW:      $(@($priority | Where-Object { $_.Priority -eq 'LOW' }).Count) files" -ForegroundColor Green
