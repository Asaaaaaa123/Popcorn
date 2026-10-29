<#
Start Bone Spatial Studio on Windows.

  .\start.ps1                          (if scripts are blocked: powershell -ExecutionPolicy Bypass -File .\start.ps1)
  .\start.ps1 --port 8800 --no-browser   every argument is passed on to "python -m webapp"

The first start creates a virtual environment in .venv (with uv when it is installed, else with Python 3.12 from the
"py" launcher, else with python/python3) and installs requirements.txt; requirements.txt is installed again whenever
it changes.
#>
$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
$venv = Join-Path $root '.venv'
$python = Join-Path $venv 'Scripts\python.exe'
$requirements = Join-Path $root 'requirements.txt'
$uv = Get-Command uv -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1

function Invoke-Native([string]$exe, [string[]]$arguments) {
    # pip, uv and the server log to stderr; under 'Stop', PowerShell 5.1 would turn that into a terminating error
    # whenever the output is redirected. Success is judged by $LASTEXITCODE instead.
    $ErrorActionPreference = 'Continue'
    & $exe @arguments
}

function Find-Python {
    # The first Python 3.10+ among the candidates, as @(executable, arguments...); $null when there is none.
    foreach ($candidate in @(@('py', '-3.12'), @('py', '-3'), @('python'), @('python3'))) {
        $command = Get-Command $candidate[0] -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
        if (-not $command) { continue }
        $extra = @($candidate | Select-Object -Skip 1)
        try {
            # the Microsoft Store placeholder "python.exe" only prints a message and fails
            $version = & $command.Source @extra -c "import sys; print('%d.%d' % sys.version_info[:2])" 2>$null
        } catch {
            continue
        }
        if ($LASTEXITCODE -eq 0 -and "$version" -match '^3\.(\d+)$' -and [int]$Matches[1] -ge 10) {
            return ,(@($command.Source) + $extra)
        }
    }
    return $null
}

function Fail([string]$message) {
    Write-Host $message -ForegroundColor Red
    exit 1
}

Push-Location -LiteralPath $root
try {
    if (-not (Test-Path -LiteralPath $python)) {
        if ($uv) {
            Write-Host 'Creating the virtual environment .venv with uv (Python 3.12) ...'
            Invoke-Native $uv.Source @('venv', '--seed', '--python', '3.12', $venv)
        } else {
            $base = Find-Python
            if (-not $base) {
                Fail ("Python 3.10 or newer (3.12 recommended) is required but was not found.`n" +
                      "Install it from https://www.python.org/downloads/ (tick 'Add python.exe to PATH') " +
                      "or run:  winget install Python.Python.3.12`n" +
                      'Then open a new terminal and run .\start.ps1 again.')
            }
            Write-Host "Creating the virtual environment .venv with $($base -join ' ') ..."
            Invoke-Native $base[0] (@($base | Select-Object -Skip 1) + @('-m', 'venv', $venv))
        }
        if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $python)) {
            Fail 'Creating the virtual environment failed (see above). Delete the .venv folder and try again.'
        }
    }

    $stamp = Join-Path $venv '.requirements.sha256'
    $hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $requirements).Hash.ToLowerInvariant()
    $installed = if (Test-Path -LiteralPath $stamp) { (Get-Content -LiteralPath $stamp -Raw).Trim() } else { '' }
    if ($hash -ne $installed) {
        Write-Host 'Installing the requirements (the first time this takes several minutes) ...'
        if ($uv) {
            Invoke-Native $uv.Source @('pip', 'install', '--python', $python, '-r', $requirements)
        } else {
            Invoke-Native $python @('-m', 'pip', 'install', '--disable-pip-version-check', '-r', $requirements)
        }
        if ($LASTEXITCODE -ne 0) {
            Fail 'Installing the requirements failed (see above). Check the internet connection and run .\start.ps1 again.'
        }
        Set-Content -LiteralPath $stamp -Value $hash -Encoding Ascii
    }

    Invoke-Native $python (@('-m', 'webapp') + @($args))
    $code = $LASTEXITCODE
} finally {
    Pop-Location
}
exit $code
