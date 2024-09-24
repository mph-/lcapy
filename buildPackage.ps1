param (
[string]$pythonPath
)
if($pythonPath){
    Set-Alias pythonPath $pythonPath
    $pytestPath = Split-Path -Path pythonPath
    $pytestPath = Join-Path -Path $pytestPath -ChildPath "pytest.exe"
    Set-Alias pytest $pytestPath
}
else{
    Set-Alias pythonPath python
    Set-Alias pytest pytest
}
$output = [string]::Concat("executing with: ", (Get-Alias pythonPath).Definition, "`npython refers to the standard python installation or the current aktiv venv")
Write-Output $output

$startDir = Get-Location

# Make sure the test circuits can be executed before building the Package
Write-Host "Execute tests" -ForegroundColor Green
Start-Sleep -Milliseconds 500
Set-Location $PSScriptRoot/NonLcapyFiles/tests

C:\Users\yannick.wieland\AppData\Local\Programs\Python\lcapy1_24\Scripts\pytest.exe
if ($LASTEXITCODE -ne 0)
{
    Write-Host "An error occured while simplifiing the test circuits" -ForegroundColor Red
    Set-Location $startDir
    return
}

# build the python package
Set-Location $PSScriptRoot
Write-Host "Building package" -ForegroundColor Green
Start-Sleep -Milliseconds 1000

try{
    pythonPath setup.py sdist bdist_wheel
}
catch{
    Write-Host "An error occured while building the package" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Red
    Set-Location $startDir
    return
}

Set-Location $PSScriptRoot
$version = pythonPath setup.py --version

$content = Get-Content -Path "NonLcapyFiles\solve.py"
if($content[0][0] -eq "#"){
    $content[0] = [string]::Concat("# for lcapy version: ", $version)
    Set-Content -Path "NonLcapyFiles\solve.py" -Value $content
}
else{
    $newContent = ,[string]::Concat("# for lcapy version: ", $version)
    $newContent += $content
    $content = $newContent
    Set-Content -Path "NonLcapyFiles\solve.py" -Value $content
}

Set-Location $PSScriptRoot
Write-Host "Successfully tested and build package" -ForegroundColor Green
Write-Host ([string]::Concat("Version is: ", $version)) -ForegroundColor Green

$oldPackages = Get-ChildItem -Path "..\Pyodide\Packages" -Filter "*lcapy*" -Recurse
foreach($item in $oldPackages){
    $path = [string]::Concat("..\Pyodide\Packages\", $item)
    Remove-Item -Path $path
    Write-Output ([string]::Concat("Removed: ", $path))
}


try {
    Copy-Item -Path "NonLcapyFiles\solve.py" -Destination "..\Pyodide\"
}
catch{
    Write-Host "could not copy solve.py" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Red
    Set-Location $startDir
    return
}
Write-Output "Copied solve.py to: ..\Pyodide\"


try{
    $newPackage = Get-ChildItem -Path "dist" -Filter "*.whl" -File | Sort-Object -Property LastWriteTime -Descending | Select-Object -First 1
}
catch {
    Write-Host "couldn't find new package in .\dist" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Red
    Set-Location $startDir
    return
}


try{
    Copy-Item -Path ([string]::Concat("dist\", $newPackage)) -Destination "..\Pyodide\Packages\"
}
catch {
    Write-Host "could not copy new package" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Red
    Set-Location $startDir
    return
}
Write-Output ([string]::Concat("Copied ", $newPackage, " to: ..\Pyodide\Packages\"))


Write-Host "Successfully updated solve.py and lcapy package in Pyodide distribution" -ForegroundColor Green
Set-Location $startDir