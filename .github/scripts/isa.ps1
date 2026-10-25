# Reports the vector tier the .NET runtime actually accelerates under the current DOTNET_* knobs,
# appends it to the job summary, and optionally writes "<label>|<vector512>|<cpu>" to -OutFile.
param([string]$Label = '', [string]$OutFile = '')

$probe = Join-Path $env:RUNNER_TEMP 'isa-probe.cs'
@(
    'using System.Runtime.InteropServices;'
    'using System.Runtime.Intrinsics;'
    'Console.WriteLine($"{RuntimeInformation.ProcessArchitecture}|{Vector512.IsHardwareAccelerated}|{Vector256.IsHardwareAccelerated}|{Vector128.IsHardwareAccelerated}");'
) | Set-Content $probe
$out = dotnet run --file $probe | Select-Object -Last 1
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
$arch, $v512, $v256, $v128 = $out -split '\|'

$cpu = if ($IsLinux) { ((Get-Content /proc/cpuinfo | Select-String '^model name' | Select-Object -First 1).Line -split ':\s*', 2)[1] }
       elseif ($IsMacOS) { sysctl -n machdep.cpu.brand_string }
       else { (Get-CimInstance Win32_Processor | Select-Object -First 1).Name }
$cpu = if ($cpu) { "$cpu".Trim() } else { 'unknown CPU' }

$line = "${arch}: Vector512=$v512 Vector256=$v256 Vector128=$v128 ($cpu)"
Write-Host $line
"**$Label**: $line" | Out-File -Append -Encoding utf8 $env:GITHUB_STEP_SUMMARY
if ($OutFile) { "$Label|$v512|$cpu" | Set-Content $OutFile }
