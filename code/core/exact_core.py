"""Small solver-free PSS_1 checker, exact rational geometry + GF(2) homology.
The checker authenticates the ENCODED branch union, not an external model adapter.
All pairs are required. Triple omissions are accepted only by boundary-span proof.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import combinations
import hashlib, json

class CertificateError(ValueError):
    pass

def require(ok, message):
    if not ok: raise CertificateError(message)

def frac(x):
    if isinstance(x, bool) or isinstance(x, float):
        raise CertificateError('Rational coefficients must be integers or strings, not float/bool')
    return F(x)

def canonical_hash(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def rr_rank(rows):
    a=[[frac(x) for x in r] for r in rows]
    if not a: return 0
    require(all(len(r)==len(a[0]) for r in a),'ragged rank matrix')
    k=0
    for c in range(len(a[0])):
        p=next((i for i in range(k,len(a)) if a[i][c]),None)
        if p is None:continue
        a[k],a[p]=a[p],a[k]
        s=a[k][c];a[k]=[v/s for v in a[k]]
        for i in range(k+1,len(a)):
            s=a[i][c]
            if s:a[i]=[v-s*w for v,w in zip(a[i],a[k])]
        k+=1
        if k==len(a):break
    return k

class BinaryBasis:
    """Column basis in F_2; integer bitsets represent chain vectors."""
    def __init__(self):self.pivots={}
    def reduce(self,v):
        while v:
            p=v.bit_length()-1
            if p not in self.pivots:break
            v^=self.pivots[p]
        return v
    def add(self,v):
        v=self.reduce(v)
        if not v:return False
        self.pivots[v.bit_length()-1]=v
        return True
    @property
    def rank(self):return len(self.pivots)

def gf2_rank(vectors):
    b=BinaryBasis()
    for v in vectors:b.add(v)
    return b.rank

def components(vertices, edges):
    parent={i:i for i in vertices}
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in edges:
        x,y=root(i),root(j)
        if x!=y:parent[x]=y
    cc={}
    for i in vertices:cc.setdefault(root(i),[]).append(i)
    return list(cc.values())

def triangle_universe(vertices,edges):
    e=set(edges)
    return [t for t in combinations(vertices,3) if all(x in e for x in combinations(t,2))]

def boundary_vector(t,edge_index):
    v=0
    for e in combinations(t,2):v^=1<<edge_index[e]
    return v

def parse_rows(raw,n):
    require(isinstance(raw,list),'rows not a list')
    out=[]
    for row in raw:
        require(len(row)==n+1,'row dimension mismatch')
        out.append(tuple(frac(x) for x in row))
    return out

def feasible(point,rows,n):
    if len(point)!=n:return False
    return all(sum(a*x for a,x in zip(r[:n],point))<=r[n] for r in rows)

def combo(rows,raw,n):
    require(len(raw)==len(rows),'multiplier length mismatch')
    lam=[frac(x) for x in raw]
    require(all(x>=0 for x in lam),'negative inequality multiplier')
    return [sum(l*r[j] for l,r in zip(lam,rows)) for j in range(n+1)]

def check_geometry(record,rows,n):
    require(isinstance(record,dict),'invalid geometry record')
    if record.get('kind')=='point':
        p=[frac(x) for x in record['point']]
        require(feasible(p,rows,n),'infeasible witness')
        return True
    if record.get('kind')=='farkas':
        v=combo(rows,record['lambda'],n)
        require(all(x==0 for x in v[:n]) and v[n]<0,'invalid strict Farkas contradiction')
        return False
    raise CertificateError('unknown evidence kind')

def check_dimension(record,rows,n):
    d=record.get('dimension')
    require(isinstance(d,int) and not isinstance(d,bool) and 0<=d<=n,'invalid dimension')
    pts=[[frac(x) for x in p] for p in record['lower_points']]
    require(len(pts)==d+1 and all(feasible(p,rows,n) for p in pts),'invalid dimension lower witnesses')
    require(rr_rank([[x-y for x,y in zip(p,pts[0])] for p in pts[1:]])==d,'dimension lower rank mismatch')
    eqs=record['equalities'];require(len(eqs)==n-d,'wrong number of affine equalities')
    normals=[]
    for eq in eqs:
        a=[frac(x) for x in eq['a']];b=frac(eq['b'])
        require(len(a)==n,'affine equality dimension')
        plus=combo(rows,eq['upper'],n);minus=combo(rows,eq['lower'],n)
        require(plus[:n]==a and plus[n]<=b,'invalid upper affine proof')
        require(minus[:n]==[-x for x in a] and minus[n]<=-b,'invalid lower affine proof')
        normals.append(a)
    require(rr_rank(normals)==n-d,'dependent affine hull equalities')
    return d

def check(system, cert):
    require(system.get('schema')=='PSS-H-union-1','unsupported system schema')
    require(cert.get('schema')=='PSS-boundary-span-1','unsupported certificate schema')
    require(cert.get('system_sha256')==canonical_hash(system),'system hash mismatch')
    n=system['ambient_dimension'];require(isinstance(n,int) and n>=1,'ambient dimension')
    rawbranches=system['branches'];m=len(rawbranches)
    rows=[parse_rows(b['rows'],n) for b in rawbranches]
    require([b['id'] for b in rawbranches]==list(range(m)),'branch IDs not canonical')
    br=cert['branches'];require(len(br)==m,'branch coverage')
    dims={};vertices=[]
    for i,r in enumerate(br):
        require(r['id']==i,'branch record ID')
        if r.get('kind')=='farkas':
            require(not check_geometry(r,rows[i],n),'unexpected feasible branch')
        else:
            require(r.get('kind')=='dimension','invalid branch record')
            dims[i]=check_dimension(r,rows[i],n);vertices.append(i)
    pairs=cert['pairs'];expected=set(combinations(vertices,2));got=set();edges=[]
    for r in pairs:
        k=tuple(r['ids']);require(k in expected and k not in got,'duplicate/invalid pair')
        got.add(k)
        if check_geometry(r,[a for i in k for a in rows[i]],n):edges.append(k)
    require(got==expected,'incomplete pair universe')
    edges.sort();index={e:i for i,e in enumerate(edges)}
    candidates=triangle_universe(vertices,edges);candset=set(candidates)
    found={};basis=BinaryBasis();positive=0;negative=0
    for r in cert['triples']:
        t=tuple(r['ids']);require(t in candset and t not in found,'duplicate/invalid triple')
        ok=check_geometry(r,[a for i in t for a in rows[i]],n);found[t]=ok
        if ok:
            basis.add(boundary_vector(t,index));positive+=1
        else:negative+=1
    unresolved=[t for t in candidates if t not in found]
    for t in unresolved:
        require(basis.reduce(boundary_vector(t,index))==0,'unresolved rank-relevant triple')
    cc=components(vertices,edges)
    z=len(edges)-len(vertices)+len(cc)
    eta=int(bool(vertices));d=max(dims.values()) if dims else None
    niso=sum(max(dims[i] for i in c)==0 for c in cc)
    pss=[eta,d,niso,len(cc),z-basis.rank]
    require(cert['claimed_pss1']==pss,'claimed signature mismatch')
    return {'pss1':pss,'branches':m,'active_branches':len(vertices),'pairs':len(pairs),'edges':len(edges),
            'triangle_candidates':len(candidates),'positive_triples':positive,'negative_triples':negative,
            'unqueried_span_redundant':len(unresolved),'boundary_rank':basis.rank,'status':'PASS'}
