"""Check the distribution hashes before running analyses that create outputs."""
from pathlib import Path
import hashlib
root=Path(__file__).resolve().parents[1]
rows=(root/'SHA256SUMS').read_text().splitlines()
for row in rows:
 expected,name=row.split('  ',1)
 path=root/name
 if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
  raise SystemExit('Hash mismatch or missing file: '+name)
print(f'PASS: {len(rows)} distribution file hashes match.')
