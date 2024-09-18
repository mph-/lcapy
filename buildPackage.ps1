param (
[string]$pythonPath
)
if($pythonPath){
    Set-Alias pythonPath $pythonPath
}
else{
    Set-Alias pythonPath python
}
$output = [string]::Concat("executing with: ", (Get-Alias pythonPath).Definition, "`npython refers to the standard python installation or the current aktiv venv")
Write-Output $output

$startDir = Get-Location

# Make sure the test circuits can be executed before building the Package
Write-Host "Execute test circuits" -ForegroundColor Green
Start-Sleep -Milliseconds 500
Set-Location $PSScriptRoot/Non-Lcapy-Files


pythonPath TestWithStandardCircuits.py
if ($LASTEXITCODE -ne 0){
    Write-Host "An error occured while simplifiing the test circuits" -ForegroundColor Red
    Set-Location $startDir
    return
}

# build the python package
Set-Location $PSScriptRoot
Write-Host "Building package" -ForegroundColor Green
Start-Sleep -Milliseconds 1000

pythonPath setup.py sdist bdist_wheel
if ($LASTEXITCODE -ne 0){
    Write-Host "An error occured while building the package" -ForegroundColor Red
    Set-Location $startDir
    return
}

Set-Location $PSScriptRoot
$version = pythonPath setup.py --version

$content = Get-Content -Path "Non-Lcapy-Files\solve.py"
if($content[0][0] -eq "#"){
    $content[0] = [string]::Concat("# for lcapy version: ", $version)
    Set-Content -Path "Non-Lcapy-Files\solve.py" -Value $content
}
else{
    $newContent = ,[string]::Concat("# for lcapy version: ", $version)
    $newContent += $content
    $content = $newContent
    Set-Content -Path "Non-Lcapy-Files\solve.py" -Value $content
}

Set-Location $startDir
Write-Host "Successfully tested and build package" -ForegroundColor Green
Write-Host ([string]::Concat("Version is: ", $version)) -ForegroundColor Green

$oldPackages = Get-ChildItem -Path "..\Pyodide\Packages" -Filter "*lcapy*" -Recurse
foreach($item in $oldPackages){
    $path = [string]::Concat("..\Pyodide\Packages\", $item)
    Remove-Item -Path $path
    Write-Output ([string]::Concat("Removed: ", $path))
}

Copy-Item -Path "Non-Lcapy-Files\solve.py" -Destination "..\Pyodide\"
if ($LASTEXITCODE -ne 0){
    Write-Host "could not copy solve.py" -ForegroundColor Red
    return
}
Write-Output "Copied solve.py to: ..\Pyodide\"

$newPackage = Get-ChildItem -Path "dist" -Filter "*.whl" -File | Sort-Object -Property LastWriteTime -Descending | Select-Object -First 1
if ($LASTEXITCODE -ne 0){
    Write-Host "could find new package in .\dist" -ForegroundColor Red
    return
}

Copy-Item -Path ([string]::Concat("dist\", $newPackage)) -Destination "..\Pyodide\Packages\"
if ($LASTEXITCODE -ne 0){
    Write-Host "could not copy new package" -ForegroundColor Red
    return
}
Write-Output ([string]::Concat("Copied ", $newPackage, " to: ..\Pyodide\Packages\"))


Write-Host "Successfully updated solve.py and lcapy package in Pyodide distribution" -ForegroundColor Green