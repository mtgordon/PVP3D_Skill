#!/bin/bash
# usage: stop_runs.sh name1 name2 ...  -- kills the febio4.exe whose command line has <name>.feb
for n in "$@"; do
  powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='febio4.exe'\" | Where-Object { \$_.CommandLine -match '$n\.feb' } | ForEach-Object { Stop-Process -Id \$_.ProcessId -Force; 'stopped $n pid ' + \$_.ProcessId }"
done
