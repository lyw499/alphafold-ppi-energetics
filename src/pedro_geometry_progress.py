"""Cheap checkpoint-only progress report; never opens coordinates."""
import json
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'results/local_interface_geometry/pedro'
records=[json.loads(f.read_text()) for f in (p/'checkpoints').glob('a???_j????.done.json')]
failed=list((p/'checkpoints').glob('a???_j????.failed.json'))
print(f"Completed jobs: {len(records)}/2156; interfaces: {sum(r['interfaces'] for r in records)}/170644; contacts: {sum(r['contacts'] for r in records)}/6671664; failed jobs: {len(failed)}")
for f in failed:print(f.name,json.loads(f.read_text())['reason'])
for name in ['validation/status.json','run_receipt.txt','report_complete.json']:
 f=p/name
 if f.exists():print(name+':',f.read_text().strip())
