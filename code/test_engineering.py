#!/usr/bin/env python3
"""Independent sanity, geometry, corruption and source adapter tests for PSS R8."""
from pathlib import Path
import subprocess,sys,json,csv,copy,random,itertools,hashlib
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from pss_engineer import (read_map,topology,make_nerve,certificates_for_nerve,clearance_square,coarsen_all_passable,get_components)
from batched_core import check as check_batch
from exact_core import CertificateError
OUT=ROOT/'qa';OUT.mkdir(exist_ok=True,parents=True)
tests=[]
def test(name,condition,extra=None):
    tests.append({'test':name,'passed':bool(condition),'details':extra})
    if not condition:raise AssertionError((name,extra))
# Simple known closed-cubical limiting cases.
base={
 'empty':([[0,0],[0,0]],(0,0)),
 'solid':([[1,1],[1,1]],(1,0)),
 'diagonal_corner_only':([[1,0],[0,1]],(1,0)),
 'separated_two':([[1,0,1]],(2,0)),
 'single_hole':([[1,1,1],[1,0,1],[1,1,1]],(1,1)),
 'open_ring_no_hole':([[1,0,1],[1,0,1],[1,1,1]],(1,0)),
}
for name,(g,expected) in base.items():
 q=topology(g);test('topology_'+name,(q['beta0_closed'],q['beta1_closed'])==expected)
# Corner-only case physically dubious for positive-width robot, must flag.
q=topology(base['diagonal_corner_only'][0]);test('diagonal_corner_warning',q['corner_touch_discrepancy'] and q['beta0_4neighbor_gridgraph']==2)
# Strict unknown-source error and bad characters.
locked=ROOT/'data/Berlin_0_256.map';g,meta=read_map(locked)
test('source_sha_locked',meta['source_sha256']==hashlib.sha256(locked.read_bytes()).hexdigest())
try:read_map(ROOT/'data/not_here.map');bad=False
except FileNotFoundError:bad=True
test('missing_source_fails',bad)
try:coarsen_all_passable(g,3);bad=False
except ValueError:bad=True
test('invalid_coarsening_fails',bad)
# Check independent all-sources full-grid Euler==background verification for all clearance states.
for name in ('Berlin_0_256','Berlin_1_256','Boston_0_256','Paris_1_256'):
 a,_=read_map(ROOT/'data'/f'{name}.map')
 old=-1
 for r in range(4):
  z=topology(clearance_square(a,r))
  test('full_grid_dual_euler_'+name+'_r'+str(r),z['dual_and_euler_agree'])
  if old>=0:test('free_space_nonincreasing_'+name+'_r'+str(r),z['free_cells']<=old)
  old=z['free_cells']
# Proof corruption checks on exact V23 example regenerated without new theory.
p=ROOT/'examples/Berlin_0_256_x98_y38_f2_clear0'
system=json.loads((p/'system.json').read_text());cert=json.loads((p/'certificate_span.json').read_text());
test('untampered_certificate_verifies',check_batch(system,cert)['status']=='PASS')
mods=[('claim_beta1',lambda c:c['claimed_pss1'].__setitem__(4,c['claimed_pss1'][4]+1)),
      ('pair_removed',lambda c:c['pairs'].pop()),
      ('fake_batch_support',lambda c:c['positive_batches'][0].__setitem__('ids',[0,1,2,3,99])),
      ('fake_positive_point',lambda c:c['positive_batches'][0].__setitem__('point',['99999','99999'])),
      ('fake_system_digest',lambda c:c.update(system_sha256='0'*64))]
for name,fn in mods:
    z=copy.deepcopy(cert);fn(z)
    try:check_batch(system,z);rejected=False
    except (CertificateError,ValueError,KeyError,IndexError,TypeError):rejected=True
    test('malformed_certificate_rejected_'+name,rejected)
# Independent oracle: GEOS/Shapely checks triangle-pair and candidate triple geometry.
try:
 from shapely.geometry import Polygon
 selections=[[[int(z) for z in row] for row in csv.reader((ROOT/'examples'/name/'admitted_window.csv').open())] for name in (
    'Berlin_0_256_x98_y38_f2_clear0','Boston_0_256_x108_y30_f2_clear0','Paris_1_256_x208_y72_f1_clear0')]
 total_pairs=total_triples=0
 for idx,patch in enumerate(selections):
  m=make_nerve(patch);s=[Polygon(t) for t in m['triangles']];positive=m['pair_positive']
  for a,b in itertools.combinations(range(len(s)),2):
   test('shapely_pair_case'+str(idx)+'_'+str(a)+'_'+str(b),((a,b) in positive)==s[a].intersects(s[b]));total_pairs+=1
  for t in m['candidates']:
   true=s[t[0]].intersection(s[t[1]]).intersects(s[t[2]])
   test('shapely_triple_case'+str(idx)+'_'+str(t),((t in m['pos'])==true));total_triples+=1
 geom={'geometry_engine':'Shapely 2.1.2 / GEOS','crosscheck_pairs':total_pairs,'crosscheck_candidates_triples':total_triples,'mismatches':0}
except ImportError:
 geom={'geometry_engine':'not installed','crosscheck_pairs':0,'crosscheck_candidates_triples':0,'status':'SKIPPED; optional dependency'}
# Check generic customer CSV -> .map adapter -> same Betti proof on a real derived V23 window.
csvsrc=ROOT/'examples/Berlin_0_256_x98_y38_f2_clear0/admitted_window.csv'
converted=OUT/'customer_like_input_from_real_Berlin_window.map'
cmd=[sys.executable,str(ROOT/'code/csv_to_map.py'),'--input',str(csvsrc),'--output',str(converted)]
subprocess.run(cmd,check=True,capture_output=True,text=True)
z,meta=read_map(converted,allow_other=True)
t=topology(z);test('user_csv_converter_topology',(t['beta0_closed'],t['beta1_closed'])==(2,1))
m=make_nerve(z);r=certificates_for_nerve(m,2,1)
test('user_csv_converter_certificates_both_pass',r['validation']['full']['status']==r['validation']['rank']['status']=='PASS')
res={'test_assertions':len(tests),'failures':sum(not q['passed'] for q in tests),'all_pass':all(t['passed'] for t in tests),'crosscheck':geom,
     'customer_supplied_csv_roundtrip':{'converted_source_sha256':meta['source_sha256'],'beta0':2,'beta1':1,'both_certificates_pass':True},
     'negative_control_rejections':len(mods),'edge_cases':len(base)}
(OUT/'engineering_validation.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
