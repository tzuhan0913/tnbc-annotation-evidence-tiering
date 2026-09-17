param(
  [string[]]$Names = @("SeuratObject_TNBC.rds"),
  [int]$MaxAttempts = 40
)

$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$outDir = Join-Path $root "data\figshare_17058077"
$manifestPath = Join-Path $outDir "relevant_files.tsv"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$manifest = Import-Csv -Delimiter "`t" -Path $manifestPath
$selected = $manifest | Where-Object { $Names -contains $_.name }
if ($selected.Count -eq 0) {
  throw "No matching Figshare files in $manifestPath"
}

foreach ($row in $selected) {
  $url = "https://ndownloader.figshare.com/files/$($row.file_id)"
  $outFile = Join-Path $outDir $row.name
  $expected = [int64]$row.size_bytes

  for ($attempt = 1; $attempt -le $MaxAttempts; $attempt++) {
    $current = 0
    if (Test-Path -LiteralPath $outFile) {
      $current = (Get-Item -LiteralPath $outFile).Length
    }

    if ($current -ge $expected) {
      Write-Host "Complete: $($row.name) ($current bytes)"
      break
    }

    Write-Host "Attempt $attempt/$MaxAttempts for $($row.name): $current / $expected bytes"
    curl.exe -L -C - --ssl-no-revoke --retry 5 --retry-delay 5 --connect-timeout 60 --speed-limit 2048 --speed-time 120 --output $outFile $url
    if ($LASTEXITCODE -ne 0) {
      Write-Host "curl exit code $LASTEXITCODE; will retry if attempts remain."
      Start-Sleep -Seconds 5
    }
  }

  $final = 0
  if (Test-Path -LiteralPath $outFile) {
    $final = (Get-Item -LiteralPath $outFile).Length
  }
  if ($final -lt $expected) {
    throw "Incomplete download for $($row.name): $final / $expected bytes"
  }
}

Write-Host "Figshare download step finished."
