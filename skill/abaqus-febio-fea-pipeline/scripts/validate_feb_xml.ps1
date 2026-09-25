# Quick well-formedness check for a FEBio .feb file before spending time
# running the solver on it. This catches "forgot to close a tag" / mismatched
# block boundaries from hand-editing, but does NOT catch FEBio's own stricter
# semantic rules (see references/febio-xml-format.md) -- a file can pass this
# and still fail to parse in febio4. Treat it as a fast pre-filter, not full
# validation.
#
# Usage: powershell -File validate_feb_xml.ps1 "C:\path\to\model.feb"

param([Parameter(Mandatory=$true)][string]$Path)

try {
  [xml]$xml = Get-Content -Raw $Path
  Write-Output "XML VALID: $Path"
} catch {
  Write-Output "XML ERROR in $Path : $($_.Exception.Message)"
  exit 1
}
