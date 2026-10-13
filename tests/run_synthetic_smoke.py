#!/usr/bin/env python3
"""Source-free PSS regression: ring topology and both exact proof modes."""
import json,subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='pss_smoke_') as d:
    d=Path(d)
    source=d/'ring.map'
    result=subprocess.run([sys.executable,str(root/'code'/'csv_to_map.py'),
        '--input',str(root/'tests'/'synthetic_8x8_ring.csv'),'--output',str(source)],
        capture_output=True,text=True,check=True)
    print('CSV to map:',result.stdout.strip())
    out=d/'generated'
    result=subprocess.run([sys.executable,str(root/'code'/'pss_engineer.py'),
        'window','--map',str(source),'--allow-other','--x','0','--y','0',
        '--factor','1','--window','8','--clearance-cells','0','--output',str(out)],
        capture_output=True,text=True)
    print('Engineering:',result.stdout.strip())
    if result.returncode: raise RuntimeError(result.stderr)
    summary=json.loads((out/'engineering_summary.json').read_text())
    result=json.loads(result.stdout.strip().splitlines()[-1])
    assert (result['beta0'],result['beta1'])==(1,1)
    assert summary['certificate_status']=='PASS'
    verification=json.loads((out/'independent_certificate_verification.json').read_text())
    assert verification['full']['status']=='PASS',verification
    assert verification['rank']['status']=='PASS',verification
    assert (out/'certificate_full.json').is_file() and (out/'certificate_rank.json').is_file()
    print('SYNTHETIC_SMOKE_PASS: exact geometry beta0=1 beta1=1, full/rank certificates verified')
