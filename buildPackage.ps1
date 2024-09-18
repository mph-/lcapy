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

$curDir = Get-Location
# Make sure the test circuits can be executed before building the Package
Write-Host "Execute test circuits" -ForegroundColor Green
Start-Sleep -Milliseconds 500
Set-Location $PSScriptRoot/Non-Lcapy-Files


pythonPath TestWithStandardCircuits.py
if ($LASTEXITCODE -ne 0){
    Write-Host "An error occured while simplifiing the test circuits" -ForegroundColor Red
    Set-Location $curDir
    return
}

# build the python package
Set-Location $PSScriptRoot
Write-Host "Building package" -ForegroundColor Green
Start-Sleep -Milliseconds 1000

pythonPath setup.py sdist bdist_wheel
if ($LASTEXITCODE -ne 0){
    Write-Host "An error occured while building the package" -ForegroundColor Red
    Set-Location $curDir
    return
}

Set-Location $curDir
Write-Host "Successfully tested and build package" -ForegroundColor Green