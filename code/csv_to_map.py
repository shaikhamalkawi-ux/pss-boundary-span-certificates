#!/usr/bin/env python3
"""Convert an engineer's binary occupancy CSV (1 free, 0 blocked) to strict MovingAI .map."""
import argparse,csv,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
a=p.parse_args()
with Path(a.input).open(newline='') as f:rows=list(csv.reader(f))
if not rows or not rows[0] or any(len(x)!=len(rows[0]) or any(v not in ('0','1') for v in x) for x in rows):
    raise SystemExit('INPUT HOLD: rectangular CSV with 0/1 ONLY; no header, missing cells or other symbols')
body=[''.join('.' if v=='1' else '@' for v in x) for x in rows]
raw=f'type octile\nheight {len(body)}\nwidth {len(body[0])}\nmap\n'+'\n'.join(body)+'\n'
p=Path(a.output);p.parent.mkdir(exist_ok=True,parents=True);p.write_bytes(raw.encode('ascii'))
print(json.dumps({'output':str(p),'width':len(body[0]),'height':len(body),'sha256':hashlib.sha256(raw.encode()).hexdigest(),'source_class':'AUTHOR-SUPPLIED occupancy; source fidelity to field data is NOT automatically certified'},sort_keys=True))
