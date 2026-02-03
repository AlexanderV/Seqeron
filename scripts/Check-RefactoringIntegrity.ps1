# Скрипт проверки целостности рефакторинга
# Запускать до и после рефакторинга для сравнения

param(
    [string]$SourcePath = "src/Seqeron",
    [string]$TestsPath = "tests",
    [switch]$SaveBaseline = $false,
    [string]$BaselineFile = "refactoring_baseline.json"
)

$ErrorActionPreference = "Stop"

function Get-ProjectMetrics {
    param([string]$Path)
    
    $metrics = @{
        Timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
        SourceFiles = @()
        TotalSourceFiles = 0
        TotalSourceLines = 0
        TestFiles = @()
        TotalTestFiles = 0
    }
    
    # Исходные файлы
    $sourceFiles = Get-ChildItem -Path $Path -Recurse -Filter "*.cs" |
        Where-Object { $_.FullName -notmatch "\\obj\\" -and $_.FullName -notmatch "\\bin\\" }
    
    foreach ($file in $sourceFiles) {
        $lines = (Get-Content $file.FullName | Measure-Object -Line).Lines
        $metrics.SourceFiles += @{
            Name = $file.Name
            Path = $file.FullName.Replace((Get-Location).Path + "\", "")
            Lines = $lines
        }
        $metrics.TotalSourceLines += $lines
    }
    $metrics.TotalSourceFiles = $sourceFiles.Count
    
    return $metrics
}

function Get-TestMetrics {
    param([string]$Path)
    
    $testFiles = Get-ChildItem -Path $Path -Recurse -Filter "*Tests.cs" |
        Where-Object { $_.FullName -notmatch "\\obj\\" -and $_.FullName -notmatch "\\bin\\" }
    
    return @{
        TotalTestFiles = $testFiles.Count
        TestFiles = $testFiles | ForEach-Object { $_.Name }
    }
}

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Refactoring Integrity Check" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Получаем метрики
Write-Host "Collecting source file metrics..." -ForegroundColor Yellow
$sourceMetrics = Get-ProjectMetrics -Path $SourcePath

Write-Host "Collecting test file metrics..." -ForegroundColor Yellow
$testMetrics = Get-TestMetrics -Path $TestsPath

# Объединяем
$allMetrics = @{
    Timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    Source = $sourceMetrics
    Tests = $testMetrics
}

# Выводим результаты
Write-Host ""
Write-Host "SOURCE FILES:" -ForegroundColor Green
Write-Host "  Total files: $($sourceMetrics.TotalSourceFiles)"
Write-Host "  Total lines: $($sourceMetrics.TotalSourceLines)"
Write-Host ""

# Группируем по папкам
$byFolder = $sourceMetrics.SourceFiles | Group-Object { Split-Path (Split-Path $_.Path) -Leaf }
foreach ($folder in $byFolder | Sort-Object Name) {
    $folderLines = ($folder.Group | Measure-Object -Property Lines -Sum).Sum
    Write-Host "  $($folder.Name): $($folder.Count) files, $folderLines lines" -ForegroundColor Gray
}

Write-Host ""
Write-Host "TEST FILES:" -ForegroundColor Green
Write-Host "  Total test files: $($testMetrics.TotalTestFiles)"

# Сохраняем baseline если нужно
if ($SaveBaseline) {
    $allMetrics | ConvertTo-Json -Depth 10 | Set-Content -Path $BaselineFile
    Write-Host ""
    Write-Host "Baseline saved to: $BaselineFile" -ForegroundColor Yellow
}

# Сравниваем с baseline если существует
if (Test-Path $BaselineFile) {
    Write-Host ""
    Write-Host "COMPARISON WITH BASELINE:" -ForegroundColor Cyan
    
    $baseline = Get-Content $BaselineFile | ConvertFrom-Json
    
    $filesDiff = $sourceMetrics.TotalSourceFiles - $baseline.Source.TotalSourceFiles
    $linesDiff = $sourceMetrics.TotalSourceLines - $baseline.Source.TotalSourceLines
    $testsDiff = $testMetrics.TotalTestFiles - $baseline.Tests.TotalTestFiles
    
    $filesColor = if ($filesDiff -eq 0) { "Green" } else { "Red" }
    $linesColor = if ([Math]::Abs($linesDiff) -lt 10) { "Green" } else { "Yellow" }
    $testsColor = if ($testsDiff -eq 0) { "Green" } else { "Red" }
    
    Write-Host "  Source files: $($sourceMetrics.TotalSourceFiles) (baseline: $($baseline.Source.TotalSourceFiles), diff: $filesDiff)" -ForegroundColor $filesColor
    Write-Host "  Source lines: $($sourceMetrics.TotalSourceLines) (baseline: $($baseline.Source.TotalSourceLines), diff: $linesDiff)" -ForegroundColor $linesColor
    Write-Host "  Test files:   $($testMetrics.TotalTestFiles) (baseline: $($baseline.Tests.TotalTestFiles), diff: $testsDiff)" -ForegroundColor $testsColor
    
    # Проверяем, что ничего не потерялось
    if ($filesDiff -ne 0) {
        Write-Host ""
        Write-Host "WARNING: Source file count changed!" -ForegroundColor Red
        
        $baselineNames = $baseline.Source.SourceFiles | ForEach-Object { $_.Name }
        $currentNames = $sourceMetrics.SourceFiles | ForEach-Object { $_.Name }
        
        $missing = $baselineNames | Where-Object { $_ -notin $currentNames }
        $added = $currentNames | Where-Object { $_ -notin $baselineNames }
        
        if ($missing) {
            Write-Host "  Missing files:" -ForegroundColor Red
            $missing | ForEach-Object { Write-Host "    - $_" -ForegroundColor Red }
        }
        
        if ($added) {
            Write-Host "  Added files:" -ForegroundColor Yellow
            $added | ForEach-Object { Write-Host "    + $_" -ForegroundColor Yellow }
        }
    }
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
