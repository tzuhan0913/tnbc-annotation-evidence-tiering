param(
  [switch]$DownloadRawTar,
  [switch]$ExtractPriority1,
  [switch]$ExtractAllManifestSamples
)

$ErrorActionPreference = "Stop"

$baseUrl = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161529/suppl"
$root = Split-Path -Parent $PSScriptRoot
$dataDir = Join-Path $root "data\GSE161529"
$rawDir = Join-Path $dataDir "raw"
$manifestPath = Join-Path $dataDir "reanalysis_sample_manifest.tsv"
$rawTarPath = Join-Path $dataDir "GSE161529_RAW.tar"

New-Item -ItemType Directory -Force -Path $dataDir, $rawDir | Out-Null

function Download-IfMissing {
  param([string]$Url, [string]$OutFile)
  if (Test-Path -LiteralPath $OutFile) {
    Write-Host "Exists: $OutFile"
  } else {
    Write-Host "Downloading: $Url"
    Invoke-WebRequest -Uri $Url -OutFile $OutFile
  }
}

function Download-LargeFile {
  param([string]$Url, [string]$OutFile)
  Write-Host "Downloading with resume support: $Url"
  curl.exe -L -C - --ssl-no-revoke --retry 10 --retry-delay 10 --connect-timeout 60 --output $OutFile $Url
  if ($LASTEXITCODE -ne 0) {
    throw "curl.exe failed with exit code $LASTEXITCODE"
  }
}

Download-IfMissing "$baseUrl/filelist.txt" (Join-Path $dataDir "filelist.txt")
Download-IfMissing "$baseUrl/GSE161529_features.tsv.gz" (Join-Path $dataDir "GSE161529_features.tsv.gz")

if ($DownloadRawTar) {
  $expectedBytes = 2359685120
  $needsDownload = $true
  if (Test-Path -LiteralPath $rawTarPath) {
    $currentBytes = (Get-Item -LiteralPath $rawTarPath).Length
    if ($currentBytes -ge $expectedBytes) {
      Write-Host "Exists: $rawTarPath"
      $needsDownload = $false
    } else {
      Write-Host "Partial file found: $currentBytes / $expectedBytes bytes"
    }
  }
  if ($needsDownload) {
    Download-LargeFile "$baseUrl/GSE161529_RAW.tar" $rawTarPath
  }
}

if ($ExtractPriority1 -or $ExtractAllManifestSamples) {
  if (-not (Test-Path -LiteralPath $rawTarPath)) {
    throw "Missing $rawTarPath. Re-run with -DownloadRawTar first."
  }

  $manifest = Import-Csv -Delimiter "`t" -Path $manifestPath
  if ($ExtractPriority1) {
    $manifest = $manifest | Where-Object { $_.priority -eq "1" }
  }

  $wanted = @()
  foreach ($row in $manifest) {
    $wanted += "$($row.file_prefix)-barcodes.tsv.gz"
    $wanted += "$($row.file_prefix)-matrix.mtx.gz"
  }

  $wanted = $wanted | Sort-Object -Unique
  Write-Host "Extracting $($wanted.Count) files to $rawDir"
  tar -xf $rawTarPath -C $rawDir $wanted
}

Write-Host "Done."
