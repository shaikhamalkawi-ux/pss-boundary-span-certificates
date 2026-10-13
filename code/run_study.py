#!/usr/bin/env python3
"""Deterministic R8 engineering evaluation using V23's nine source-derived windows.
Not independent sampling; no hypothesis tests or performance claim are made.
"""
from __future__ import annotations
import sys,csv,json,time,hashlib
from pathlib import Path
from pss_engineer import (ROOT,read_map,clearance_square,coarsen_all_passable,topology,make_nerve,certificates_for_nerve,draw_window,save_json,compact)

cases=json.loads((ROOT/'data/V23_selected_cases.json').read_text())
old=json.loads((ROOT/'source_history/V23_summary_original.json').read_text())
old={q['case_id']:q for q in old['all_requested_baseline_results']}
locked_proof_hashes=json.loads((ROOT/'data/V23_original_proof_hashes.json').read_text())
root=ROOT/'results';root.mkdir(parents=True,exist_ok=True)
rows=[];valid=[]
for k,c in enumerate(cases):
    path=ROOT/'data'/f"{c['map']}.map"; original,source=read_map(path)
    for r in (0,1,2):
        cleared=clearance_square(original,r)
        coarse=coarsen_all_passable(cleared,c['coarse_factor'])
        x,y=c['x'],c['y']
        p=[line[x:x+8] for line in coarse[y:y+8]]
        ts=topology(p)
        baseline=c['case_id'];n=ts['free_cells'];model=make_nerve(p,x0=c['coarse_factor']*x,y0=c['coarse_factor']*y,stride=c['coarse_factor'])
        # Original nine cases: exact independent re-creation and byte-for-byte proof parity against V23.
        # New clearance=1 case: full+span certificates via existing V23 verifier.
        # Clearance=2: direct cubical computation; any nonempty case may be certified by separate CLI.
        cert=None
        if model and r<=1:
            cert=certificates_for_nerve(model,ts['beta0_closed'],ts['beta1_closed'])
            if r==0:
                exact=['branch_count','pair_predicate_queries','clique_triple_queries','true_positive_triples','true_negative_triples','rank_positive','full_batched_bytes','span_batched_bytes']
                got=[cert['stats'][key] for key in ['branches','pairs_tested_all','candidate_triples_tested_all','positive_triples','negative_triples','boundary_rank','proof_bytes_full','proof_bytes_rank']]
                expect=[old[baseline][v] for v in exact]
                if got!=expect:raise AssertionError((baseline,dict(zip(exact,zip(got,expect)))))
                if (ts['beta0_closed'],ts['beta1_closed'])!=(old[baseline]['beta0'],old[baseline]['beta1']):raise AssertionError('baseline beta mismatch')
                hashes={fname:hashlib.sha256(compact(cert[key])).hexdigest() for fname,key in [('system.json','system'),('full_batched.json','full'),('span_batched.json','rank')]}
                if any(hashes[fname]!=locked_proof_hashes[baseline][fname]['sha256'] for fname in hashes):
                    raise AssertionError(('original_exact_payload_identity_mismatch',baseline,hashes))
                valid.append({'case':baseline,'original_file_sha256_parity':True,'original_proof_bytes_parity':True,'original_pss_parity':True,'frozen_core_verifier':True})
            if k in (0,4,8):
                out=ROOT/'examples'/f"{c['map']}_x{x}_y{y}_f{c['coarse_factor']}_clear{r}"
                out.mkdir(parents=True,exist_ok=True)
                for fname,v in [('system',cert['system']),('certificate_full',cert['full']),('certificate_span',cert['rank']),('proof_verifier',cert['validation'])]:save_json(out/(fname+'.json'),v)
                with (out/'admitted_window.csv').open('w',newline='') as f:wr=csv.writer(f);wr.writerows([[int(q) for q in pp] for pp in p])
                draw_window(out/'figure.png',original,cleared,coarse,p,x,y,c['coarse_factor'],ts,cert['stats'])
                save_json(out/'case_result.json',{'source':source,'case':baseline,'clearance_cells':r,'topology':ts,'proof':cert['stats'],'verified':cert['validation']['rank']['status']})
        item={'case_id':baseline,'map':c['map'],'factor':c['coarse_factor'],'coarse_x':x,'coarse_y':y,
              'clearance_cells':r,'free_cells':n,'beta0_closed':ts['beta0_closed'],'beta1_closed':ts['beta1_closed'],
              'beta0_4grid':ts['beta0_4neighbor_gridgraph'],'corner_touch_warning':ts['corner_touch_discrepancy'],
              'certificate_generated':bool(cert),'certificate_pass':cert['validation']['rank']['status'] if cert else 'not_requested',
              'branches':cert['stats']['branches'] if cert else 2*n,
              'offline_bytes_saved':cert['stats']['proof_bytes_saved'] if cert else None}
        rows.append(item)
        print(json.dumps({v:item[v] for v in ('case_id','clearance_cells','free_cells','beta0_closed','beta1_closed','certificate_pass')},sort_keys=True),flush=True)
with (root/'nine_V23_cases_x_3_clearance_scenarios.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
save_json(root/'nine_case_source_parity.json',{'comparison_basis':'V23 source and original selected cases; same frozen proof checker','nine_case_checks':valid,
 'nine_case_exact_parity':len(valid)==9,'clearance_scenarios':[0,1,2],'rows':rows})
print('NINE_CASE_BASELINE_REPRODUCED',len(valid),'ADDITIONAL_CELL_RADIUS_1_CERTIFICATES',sum(x['clearance_cells']==1 and x['certificate_generated'] for x in rows))
