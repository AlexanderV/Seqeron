# Скрипт миграции using-директив в тестах
# Запускать после переноса всех файлов в новые пакеты

param(
    [string]$TestsPath = "tests",
    [switch]$DryRun = $false,
    [switch]$Verbose = $false
)

$ErrorActionPreference = "Stop"

# Маппинг старых using на новые
$usingMappings = @"
using Seqeron.Genomics.Core;
using Seqeron.Genomics.IO;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Annotation;
using Seqeron.Genomics.Phylogenetics;
using Seqeron.Genomics.Population;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Genomics.MolTools;
using Seqeron.Genomics.Reports;
"@

function Update-UsingDirectives {
    param(
        [string]$FilePath,
        [bool]$DryRun
    )
    
    $content = Get-Content $FilePath -Raw
    
    # Проверяем, есть ли using Seqeron.Genomics;
    if ($content -match "using Seqeron\.Genomics;") {
        
        # Заменяем using Seqeron.Genomics; на новые using
        $newContent = $content -replace "using Seqeron\.Genomics;", $usingMappings
        
        if ($DryRun) {
            Write-Host "[DRY-RUN] Would update: $FilePath" -ForegroundColor Yellow
        } else {
            Set-Content -Path $FilePath -Value $newContent -NoNewline
            Write-Host "[UPDATED] $FilePath" -ForegroundColor Green
        }
        
        return $true
    }
    
    return $false
}

# Найти все .cs файлы в тестах
$testFiles = Get-ChildItem -Path $TestsPath -Recurse -Filter "*.cs" |
    Where-Object { $_.FullName -notmatch "\\obj\\" -and $_.FullName -notmatch "\\bin\\" }

Write-Host "Found $($testFiles.Count) test files" -ForegroundColor Cyan
Write-Host ""

$updatedCount = 0
$skippedCount = 0

foreach ($file in $testFiles) {
    $updated = Update-UsingDirectives -FilePath $file.FullName -DryRun $DryRun
    
    if ($updated) {
        $updatedCount++
    } else {
        $skippedCount++
        if ($Verbose) {
            Write-Host "[SKIPPED] $($file.FullName)" -ForegroundColor Gray
        }
    }
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Summary:" -ForegroundColor Cyan
Write-Host "  Updated: $updatedCount files" -ForegroundColor Green
Write-Host "  Skipped: $skippedCount files" -ForegroundColor Gray
Write-Host "========================================" -ForegroundColor Cyan

if ($DryRun) {
    Write-Host ""
    Write-Host "This was a DRY RUN. No files were modified." -ForegroundColor Yellow
    Write-Host "Run without -DryRun to apply changes." -ForegroundColor Yellow
}
