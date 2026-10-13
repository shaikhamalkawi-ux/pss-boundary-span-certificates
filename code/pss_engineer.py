#!/usr/bin/env python3
"""PSS engineering application: exact topology & independently checked rational certificates.

Data contract: MovingAI octile occupancy grid, '.' = admitted square, '@'=blocked.
This demonstrates mathematical *closed-square* topology, not robot collision-free paths.
Every model/assumption is emitted alongside numerical results.
Copyright: new demonstration logic prepared for PSS; core proof verifier inherited verbatim
from frozen V23 (its files and hashes are preserved in CORE_ORIGIN_SHA256.json).
"""
from __future__ import annotations
import argparse,csv,hashlib,itertools,json,sys,time
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F

sys.path.insert(0,str(Path(__file__).resolve().parent/'core'))
from exact_core import BinaryBasis, boundary_vector, canonical_hash
from batched_core import check as check_batch
from geometry2d import sparse_farkas

ROOT=Path(__file__).resolve().parents[1]
LOCK={x['name']:x for x in json.loads((ROOT/'data/map_source_lock.json').read_text())}


def compact(x):return (json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode('utf-8')
def sha(x):return hashlib.sha256(x).hexdigest()
def save_json(p,x):p=Path(p);p.parent.mkdir(exist_ok=True,parents=True);p.write_bytes(compact(x))

def read_map(path:Path,allow_other=False):
    raw=path.read_bytes();digest=sha(raw)
    if not allow_other:
        if path.stem not in LOCK:raise ValueError('Source basename is not locked; use --allow-other explicitly')
        meta=LOCK[path.stem]
        if digest!=meta['raw_sha256']:raise ValueError('Raw map SHA-256 mismatch: source model is not authenticated')
    ls=raw.decode('ascii').splitlines()
    if len(ls)<5 or not ls[0].startswith('type ') or not ls[1].startswith('height ') or not ls[2].startswith('width ') or ls[3]!='map':
        raise ValueError('Expected MovingAI octile format with first four header lines')
    h=int(ls[1].split()[1]);w=int(ls[2].split()[1])
    if len(ls)!=h+4 or any(len(x)!=w or not set(x)<={'.','@'} for x in ls[4:]):
        raise ValueError('Bad map dimensions or unknown symbols; only . and @ are admitted')
    return [[ch=='.' for ch in line] for line in ls[4:]],{'name':path.name,'source_sha256':digest,'source_bytes':len(raw),'width':w,'height':h,'free_symbol':'.','blocked_symbol':'@'}

def clearance_square(grid,r):
    if not isinstance(r,int) or r<0:raise ValueError('clearance is a nonnegative *cell* integer')
    if r==0:return grid
    H=len(grid);W=len(grid[0]);result=[[False]*W for _ in range(H)]
    # Exact Chebyshev-square clearance around cell centers, applied on FULL raw map first.
    for y in range(r,H-r):
        for x in range(r,W-r):
            result[y][x]=all(grid[yy][xx] for yy in range(y-r,y+r+1) for xx in range(x-r,x+r+1))
    return result

def coarsen_all_passable(g,f):
    if f<1:raise ValueError('factor must be >=1')
    H=len(g);W=len(g[0]);
    if H%f or W%f:raise ValueError('coarse factor must divide both full-map dimensions exactly')
    return [[all(g[yy][xx] for yy in range(y,y+f) for xx in range(x,x+f))
            for x in range(0,W,f)] for y in range(0,H,f)]

def count_components(g,foreground=True,diag=False):
    H=len(g);W=len(g[0]); seen=set();n=0;touches=[]
    offsets=[(1,0),(-1,0),(0,1),(0,-1)]
    if diag:offsets += [(dx,dy) for dx in (-1,1) for dy in (-1,1)]
    for y in range(H):
        for x in range(W):
            if bool(g[y][x])!=foreground or (x,y) in seen:continue
            n+=1;edge=(x in (0,W-1) or y in (0,H-1));seen.add((x,y));queue=[(x,y)]
            for ax,ay in queue:
                for dx,dy in offsets:
                    xx,yy=ax+dx,ay+dy
                    if 0<=xx<W and 0<=yy<H and bool(g[yy][xx])==foreground and (xx,yy) not in seen:
                        seen.add((xx,yy));queue.append((xx,yy));edge |= (xx in (0,W-1) or yy in (0,H-1))
            touches.append(bool(edge))
    return n,touches

def topology(g):
    """Two *independent* first-homology calculations for finite unions of CLOSED unit squares."""
    if not g or not g[0] or any(len(row)!=len(g[0]) for row in g):raise ValueError('invalid grid')
    H=len(g);W=len(g[0]);free=sum(map(sum,g))
    b0,_=count_components(g,True,True)
    b0_4,_=count_components(g,True,False)
    bg_comps,borders=count_components(g,False,False)
    holes=sum(not t for t in borders)
    # Independent Euler-cell realization. All 4 vertices and 4 edges of every closed free square admitted.
    vertices=set();edges=set()
    for y,row in enumerate(g):
        for x,freecell in enumerate(row):
            if not freecell:continue
            vs=((x,y),(x+1,y),(x+1,y+1),(x,y+1))
            vertices.update(vs)
            edges.update(tuple(sorted((vs[k],vs[(k+1)%4]))) for k in range(4))
    chi=len(vertices)-len(edges)+free
    beta1_euler=b0-chi
    if beta1_euler!=holes:raise AssertionError(f'Independent cubical audit mismatch: Euler={beta1_euler} dual_4_bg={holes}')
    return {'free_cells':free,'total_cells':H*W,'free_fraction':free/(H*W),'beta0_closed':b0,'beta1_closed':holes,
            'beta0_4neighbor_gridgraph':b0_4,'chi':chi,'vertices':len(vertices),'edges':len(edges),
            'faces':free,'background_4connected_components':bg_comps,'dual_and_euler_agree':True,
            'corner_touch_discrepancy':bool(b0!=b0_4)}

def make_nerve(patch,x0=0,y0=0,stride=1):
    triangles=[];rows=[];support=defaultdict(set)
    for py,line in enumerate(patch):
        for px,active in enumerate(line):
            if not active:continue
            xx=x0+stride*px;yy=y0+stride*py;f=stride
            for verts,rr in [
                (((xx,yy),(xx+f,yy),(xx+f,yy+f)),((1,0,xx+f),(0,-1,-yy),(-1,1,yy-xx))),
                (((xx,yy),(xx,yy+f),(xx+f,yy+f)),((-1,0,-xx),(0,1,yy+f),(1,-1,xx-yy)))]:
                i=len(triangles);triangles.append(verts);rows.append(rr)
                for pt in verts:support[pt].add(i)
    m=len(triangles)
    if m==0:return None
    edges=sorted({p for star in support.values() for p in itertools.combinations(sorted(star),2)})
    eidx={e:i for i,e in enumerate(edges)}
    pair_positive=set(edges)
    adjacent=[set() for _ in range(m)]
    for a,b in edges:adjacent[a].add(b);adjacent[b].add(a)
    candidates=[]
    for a in range(m):
        for b in sorted(z for z in adjacent[a] if z>a):
            for c in sorted(z for z in adjacent[a].intersection(adjacent[b]) if z>b):candidates.append((a,b,c))
    truepos={t for ids in support.values() for t in itertools.combinations(sorted(ids),3)}
    stars=[{'kind':'batch_point','point':[str(v) for v in pt],'ids':sorted(ids)}
           for pt,ids in sorted(support.items()) if len(ids)>=3]
    return {'triangles':triangles,'rows':rows,'edges':edges,'pair_positive':pair_positive,'candidates':candidates,'pos':truepos,
            'stars':stars,'n':m,'eidx':eidx,'support':support}

def certificates_for_nerve(model,expected_b0,expected_b1):
    """Complete source-truth discovery, THEN offline full/rank certificate retention.
    Does NOT claim online-oracle query saving. All pairs/triples are evaluated/discovered.
    """
    rows=model['rows'];tri=model['triangles'];n=model['n'];points=model['support'];edges=model['edges'];eidx=model['eidx']
    system={'schema':'PSS-H-union-1','ambient_dimension':2,
      'branches':[{'id':i,'rows':[[str(a),str(b),str(c)] for a,b,c in rr]} for i,rr in enumerate(rows)]}
    branch_evidence=[{'id':i,'kind':'dimension','dimension':2,'lower_points':[[str(u),str(v)] for u,v in tri[i]],'equalities':[]} for i in range(n)]
    pair_records=[]
    pos=model['pair_positive']
    for a in range(n):
        pa=frozenset(tri[a])
        for b in range(a+1,n):
            pb=frozenset(tri[b])
            if (a,b) in pos:
                p=min(pa.intersection(pb))
                pair_records.append({'kind':'point','ids':[a,b],'point':[str(v) for v in p]})
            else:
                lam=sparse_farkas(list(rows[a])+list(rows[b]))
                pair_records.append({'kind':'farkas','ids':[a,b],'lambda':[str(v) for v in lam]})
    basis=BinaryBasis(); neg=[];essential=[]
    for t in sorted(model['pos']):basis.add(boundary_vector(t,eidx))
    for t in model['candidates']:
        if t in model['pos']:continue
        lam=sparse_farkas([rr for i in t for rr in rows[i]])
        rec={'kind':'farkas','ids':list(t),'lambda':[str(v) for v in lam]}
        neg.append(rec)
        if basis.reduce(boundary_vector(t,eidx))!=0:essential.append(rec)
    b0=get_components(n,edges)
    b1=len(edges)-n+b0-basis.rank
    if (b0,b1)!=(expected_b0,expected_b1):raise AssertionError(f'nerve/source mismatch {(b0,b1)} vs {(expected_b0,expected_b1)}')
    keys={'schema':'PSS-batched-boundary-1','system_sha256':canonical_hash(system),
          'branches':branch_evidence,'pairs':pair_records,
          'positive_batches':model['stars'],'claimed_pss1':[1,2,0,b0,b1]}
    full={**keys,'mode':'full','negative_triples':neg}
    rank={**keys,'mode':'rank','negative_triples':essential}
    v1=check_batch(system,full);v2=check_batch(system,rank)
    assert v1['status']==v2['status']=='PASS' and v1['pss1']==v2['pss1']==[1,2,0,b0,b1]
    return {'system':system,'full':full,'rank':rank,'validation':{'full':v1,'rank':v2},'stats':{
       'branches':n,'pairs_tested_all':n*(n-1)//2,'pairs_positive':len(edges),'pairs_negative':len(pair_records)-len(edges),
       'candidate_triples_tested_all':len(model['candidates']),'positive_triples':len(model['pos']),
       'negative_triples':len(neg),'negative_essential':len(essential),'boundary_rank':basis.rank,
       'positive_point_batches':len(model['stars']),'proof_bytes_full':len(compact(full)),
       'proof_bytes_rank':len(compact(rank)),
       'proof_bytes_saved':len(compact(full))-len(compact(rank)),
       'system_bytes':len(compact(system)),
       'mode':'offline proof retention; not query reduction'}}

def get_components(n,edges):
    root=list(range(n))
    def find(x):
        while root[x]!=x:x=root[x]
        return x
    for a,b in edges:
        root[find(a)]=find(b)
    return len({find(i) for i in range(n)})

def draw_window(out,full,after,coarse,patch,x,y,f,info,proof_stats):
    try:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        import numpy as np
    except ImportError:
        return False
    fig,axes=plt.subplots(1,3,figsize=(11.8,4.2),layout='constrained')
    def image(g,ax,title):
        ax.imshow(np.asarray(g,dtype=np.uint8),cmap='Greys_r',vmin=0,vmax=1,origin='upper',interpolation='nearest')
        ax.set_title(title,size=10);ax.set_xticks([]);ax.set_yticks([])
    image(full,axes[0],'Raw benchmark occupancy (256 x 256)' if len(full)==256 else 'Raw occupancy')
    axes[0].add_patch(Rectangle((f*x-0.5,f*y-0.5),f*len(patch[0]),f*len(patch),fill=False,edgecolor='#f16e50',linewidth=2.5))
    image(patch,axes[1],f'Admitted {len(patch[0])}x{len(patch)} cells; clearance + coarsening')
    image(patch,axes[2],'Closed-square union and fixed diagonal split')
    axes[2].set_xticks(np.arange(-.5,len(patch[0]),1),minor=True)
    axes[2].set_yticks(np.arange(-.5,len(patch),1),minor=True)
    axes[2].grid(which='minor',color='#b1b1b1',linewidth=.6)
    for yy,row in enumerate(patch):
        for xx,on in enumerate(row):
            if on:axes[2].plot([xx-.5,xx+.5],[yy-.5,yy+.5],color='#178576',linewidth=.9)
    if proof_stats is not None:
        txt=f"$\\beta_0={info['beta0_closed']}~;~\\beta_1={info['beta1_closed']}$\n{proof_stats['branches']} branches; {proof_stats['candidate_triples_tested_all']} triple candidates"
    else:
        txt=f"$\\beta_0={info['beta0_closed']}~;~\\beta_1={info['beta1_closed']}$\nDirect cubical check only"
    axes[2].text(.5,-.18,txt,ha='center',va='center',transform=axes[2].transAxes,size=10)
    fig.suptitle('PSS Engineering Adapter | source-derived closed-cell topology (NOT route safety)',size=11)
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out,dpi=180,bbox_inches='tight');plt.close(fig)
    return True

def run_window(args):
    p=Path(args.map)
    original,source=read_map(p,allow_other=args.allow_other)
    eroded=clearance_square(original,args.clearance_cells)
    coarse=coarsen_all_passable(eroded,args.factor)
    x,y=args.x,args.y;w=args.window
    if min(x,y)<0 or x+w>len(coarse[0]) or y+w>len(coarse):raise ValueError('window outside coarsened map')
    patch=[row[x:x+w] for row in coarse[y:y+w]]
    topo=topology(patch)
    out=Path(args.output);out.mkdir(exist_ok=True,parents=True)
    with (out/'admitted_window.csv').open('w',newline='') as f:
        q=csv.writer(f);q.writerows([[int(v) for v in row] for row in patch])
    # Large full-source or operator-only sweeps must not construct a massive nerve.
    model=make_nerve(patch,x0=x*args.factor,y0=y*args.factor,stride=args.factor) if args.certificate and topo['free_cells'] else None
    proof=None
    if args.certificate and model:
        if model['n']>args.max_certificate_branches:raise ValueError('More branches than declared engineering CLI certificate budget; use --no-certificate')
        proof=certificates_for_nerve(model,topo['beta0_closed'],topo['beta1_closed'])
        save_json(out/'system.json',proof['system'])
        save_json(out/'certificate_full.json',proof['full'])
        save_json(out/'certificate_rank.json',proof['rank'])
        save_json(out/'independent_certificate_verification.json',proof['validation'])
    # The --no-certificate or empty-region path is intentionally direct-cubical only.
    metadata={'source':source,'transform':{'clearance_cells':args.clearance_cells,'clearance_norm':'Chebyshev / axis-aligned square',
          'clearance_applied_to':'entire source map BEFORE subwindow and coarsening','coarse_factor':args.factor,
          'coarse_passable_rule':'ALL original cells in coarse footprint are passable','window_origin_coarse_xy':[x,y],
          'window_size_cells':w,'raw_window_origin_xy':[x*args.factor,y*args.factor],
          'square_topology':'union of CLOSED cells (diagonal corner contact counts for beta0)','triangular_cover':'NW-SE split of every admitted square, exact union'},
          'topology':topo,'pss1':[int(bool(topo['free_cells'])), 2 if topo['free_cells'] else None,0,topo['beta0_closed'],topo['beta1_closed']],
          'certificate_status':'PASS' if proof else ('NOT_APPLICABLE_EMPTY_REGION' if topo['free_cells']==0 else 'NOT_REQUESTED'),
          'certificate_pass': bool(proof and proof['validation']['rank']['status']=='PASS'),
          'certificate':proof['stats'] if proof else None,
          'engineering_claim_boundary':[
            'beta0/beta1 are properties of modeled closed-square geometry, not navigability or travel-time guarantees.',
            'Diagonal-only cell contact can connect the closed geometric union but may be blocked for finite-width robots.',
            'No physical distance, map scale, obstacle height, response time, or safe clearance can be inferred from benchmark pixels.',
            'Chebyshev clearance in map cells is a declared scenario, NOT field-measured robot radius.',
            'All pair and candidate triple predicates are discovered before offline proof retention; query speedup NOT claimed.',
            'If complete map is trusted, direct cubical topology is a stronger and often cheaper comparator.'
          ]}
    save_json(out/'engineering_summary.json',metadata)
    draw_window(out/'engineering_figure.png',original,eroded,coarse,patch,x,y,args.factor,topo,proof['stats'] if proof else None)
    print(json.dumps({'name':p.name,'window':[x,y,w],'factor':args.factor,'clearance':args.clearance_cells,
                      'beta0':topo['beta0_closed'],'beta1':topo['beta1_closed'],'free':topo['free_cells'],
                      'corner_touch_warning':topo['corner_touch_discrepancy'],
                      'certified':bool(proof),'output':str(out)},sort_keys=True))
    return metadata

def run_full(args):
    g,source=read_map(Path(args.map),args.allow_other)
    output=Path(args.output);output.mkdir(exist_ok=True,parents=True)
    rows=[]
    for r in args.clearance_cells:
        a=clearance_square(g,r)
        tr=topology(a)
        rows.append({'source':source['name'],'sha256':source['source_sha256'],'clearance_cells':r,**tr})
    with (output/'full_map_clearance_topology.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    save_json(output/'full_map_source_contract.json',{'source':source,'clearance_units':'native pixels only; NOT metres','method':'direct cubical Euler and independent dual-background flood fill','results':rows})
    for x in rows: print(json.dumps({'map':x['source'],'clearance':x['clearance_cells'],'free':x['free_cells'],'beta0':x['beta0_closed'],'beta1':x['beta1_closed'],'beta0_4':x['beta0_4neighbor_gridgraph']}))
    return rows

def main():
    a=argparse.ArgumentParser(description='PSS engineering tool: authenticated grid -> topology -> optional verified exact proof')
    sub=a.add_subparsers(dest='task',required=True)
    b=sub.add_parser('window',help='compute source-derived window + complete certificate if requested')
    b.add_argument('--map',required=True);b.add_argument('--x',required=True,type=int);b.add_argument('--y',required=True,type=int)
    b.add_argument('--window',type=int,default=8);b.add_argument('--factor',type=int,default=2)
    b.add_argument('--clearance-cells',type=int,default=0);b.add_argument('--output',required=True)
    b.add_argument('--allow-other',action='store_true',help='accept an unlocked external MovingAI .map with its SHA recorded')
    b.add_argument('--max-certificate-branches',type=int,default=128)
    b.add_argument('--certificate',dest='certificate',action='store_true',default=True)
    b.add_argument('--no-certificate',dest='certificate',action='store_false')
    b.set_defaults(func=run_window)
    f=sub.add_parser('full',help='all-source exact topology sweep without enormous triangle proof')
    f.add_argument('--map',required=True);f.add_argument('--clearance-cells',type=int,nargs='+',default=[0,1,2,3]);f.add_argument('--output',required=True)
    f.add_argument('--allow-other',action='store_true');f.set_defaults(func=run_full)
    z=a.parse_args();z.func(z)
if __name__=='__main__':main()
