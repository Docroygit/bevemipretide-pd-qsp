param([Parameter(Mandatory = $true)][string]$WorkDir)
# Renders the tracked documents with Word: PDF with markup, clean .docx with all changes accepted, and clean PDF.
# Word is driven from local copies in $WorkDir (a local, non-synced folder). If a hidden Word instance is ever
# force-closed, clear its entries under HKCU:\Software\Microsoft\Office\16.0\Word\Resiliency\DocumentRecovery,
# otherwise the next hidden instance waits on an invisible recovery prompt.
$project = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force $WorkDir | Out-Null
foreach ($name in @('Manuscript_bjp', 'Supplementary_Tables')) {
  $w = New-Object -ComObject Word.Application
  $w.Visible = $false
  $w.DisplayAlerts = 0
  $local = Join-Path $WorkDir ($name + '_in.docx')
  Copy-Item -LiteralPath (Join-Path $project "documents\$($name)_tracked.docx") -Destination $local -Force
  $d = $w.Documents.Open($local, $false, $true, $false)
  $revisions = $d.Revisions.Count
  $d.ExportAsFixedFormat((Join-Path $WorkDir "$($name)_tracked.pdf"), 17, $false, 0, 0, 1, 1, 7)
  $d.Revisions.AcceptAll()
  $d.SaveAs2((Join-Path $WorkDir "$($name)_clean.docx"), 16)
  $d.ExportAsFixedFormat((Join-Path $WorkDir "$($name)_clean.pdf"), 17)
  "$($name): $revisions revisions, $($d.ComputeStatistics(2)) pages after acceptance"
  $d.Close(0)
  $w.Quit()
}
