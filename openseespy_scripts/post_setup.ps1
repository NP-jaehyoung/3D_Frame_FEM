$envPath = $env:CONDA_PREFIX
New-Item -ItemType Directory -Force -Path (Join-Path $envPath "DLLs") | Out-Null

Copy-Item (Join-Path $envPath "Library\bin\tcl86t.dll") (Join-Path $envPath "DLLs\tcl86t.dll") -Force

if (Test-Path (Join-Path $envPath "Library\bin\tk86t.dll")) {
  Copy-Item (Join-Path $envPath "Library\bin\tk86t.dll") (Join-Path $envPath "DLLs\tk86t.dll") -Force
}

Write-Host "Done: copied tcl/tk DLLs into $envPath\DLLs"
