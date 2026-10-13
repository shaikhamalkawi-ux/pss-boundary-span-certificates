"""Independent certificate GENERATOR utilities for bounded rational 2-D polyhedra.
No floating-point predicate is used. The exact_core checker does not import this module.
"""
from fractions import Fraction as F
from itertools import combinations
from exact_core import feasible, rr_rank

def dedupe(poly):
    out=[]
    for p in poly:
        p=tuple(p)
        if not out or p!=out[-1]:out.append(p)
    if len(out)>1 and out[0]==out[-1]:out.pop()
    return out

def vertices(rows):
    # Every admitted generator input explicitly contains |x|,|y|<=B, B<=1000.
    poly=[(F(-1000),F(-1000)),(F(1000),F(-1000)),(F(1000),F(1000)),(F(-1000),F(1000))]
    for a,b,c in rows:
        if not poly:break
        out=[]
        for u,v in zip(poly,poly[1:]+poly[:1]):
            fu=a*u[0]+b*u[1]-c;fv=a*v[0]+b*v[1]-c
            if fu<=0:out.append(u)
            if (fu<0 and fv>0) or (fu>0 and fv<0):
                t=fu/(fu-fv);out.append(tuple(x+t*(y-x) for x,y in zip(u,v)))
        poly=dedupe(out)
    return poly

def cross(u,v):return u[0]*v[1]-u[1]*v[0]

def _sparse_farkas(rows):
    n=len(rows)
    for i,r in enumerate(rows):
        if r[0]==r[1]==0 and r[2]<0:
            ans=[F(0)]*n;ans[i]=F(1);return ans
    for i,j in combinations(range(n),2):
        u,v=rows[i],rows[j]
        if cross(u,v)!=0:continue
        h=next((k for k in (0,1) if v[k]),None)
        if h is None:continue
        l=-u[h]/v[h]
        if l>=0 and all(u[k]+l*v[k]==0 for k in (0,1)) and u[2]+l*v[2]<0:
            ans=[F(0)]*n;ans[i]=F(1);ans[j]=l;return ans
    for i,j,k in combinations(range(n),3):
        u,v,w=rows[i],rows[j],rows[k]
        l=[cross(v,w),cross(w,u),cross(u,v)]
        if all(x<=0 for x in l):l=[-x for x in l]
        if all(x>=0 for x in l) and any(l) and sum(x*r[2] for x,r in zip(l,[u,v,w]))<0:
            ans=[F(0)]*n
            for a,b in zip([i,j,k],l):ans[a]=b
            return ans
    raise ValueError('No 2-D Farkas certificate found; input may be outside bounded generator contract')

def sparse_farkas(rows):
    # Keep strongest parallel rows and expand the proof to original row indices.
    chosen={}
    for i,r in enumerate(rows):
        norm=abs(r[0]) if r[0] else abs(r[1])
        if not norm:
            if r[2]<0:
                out=[F(0)]*len(rows);out[i]=F(1);return out
            continue
        key=(r[0]/norm,r[1]/norm);rhs=r[2]/norm
        if key not in chosen or rhs<chosen[key][0]:chosen[key]=(rhs,i)
    indices=sorted(v[1] for v in chosen.values())
    small=_sparse_farkas([rows[i] for i in indices]);out=[F(0)]*len(rows)
    for i,l in zip(indices,small):out[i]=l
    return out


def implication(rows,target):
    a,b,c=target;n=len(rows)
    if a==b==0 and c>=0:return [F(0)]*n
    for i,r in enumerate(rows):
        h=next((k for k in (0,1) if r[k]),None)
        if h is None:continue
        l=target[h]/r[h]
        if l>=0 and l*r[0]==a and l*r[1]==b and l*r[2]<=c:
            out=[F(0)]*n;out[i]=l;return out
    for i,j in combinations(range(n),2):
        u,v=rows[i],rows[j];det=cross(u,v)
        if not det:continue
        l=cross((a,b),v)/det;s=cross(u,(a,b))/det
        if l>=0 and s>=0 and l*u[2]+s*v[2]<=c:
            out=[F(0)]*n;out[i]=l;out[j]=s;return out
    raise ValueError('No valid inequality implication found')

def evidence(rows,known_vertices=None):
    pts=vertices(rows) if known_vertices is None else known_vertices
    if pts:
        p=tuple(sum(v[j] for v in pts)/len(pts) for j in range(2))
        assert feasible(p,rows,2)
        return {'kind':'point','point':[str(x) for x in p]}
    return {'kind':'farkas','lambda':[str(x) for x in sparse_farkas(rows)]}

def dimension_certificate(rows):
    pts=vertices(rows)
    if not pts:return evidence(rows,pts)
    chosen=[pts[0]];diff=[]
    for p in pts[1:]:
        v=[p[j]-pts[0][j] for j in range(2)]
        if rr_rank(diff+[v])>len(diff):diff.append(v);chosen.append(p)
    d=len(diff);eqs=[];normals=[]
    for r in rows:
        if all(r[0]*p[0]+r[1]*p[1]==r[2] for p in pts) and rr_rank(normals+[r[:2]])>len(normals):
            eqs.append({'a':list(map(str,r[:2])),'b':str(r[2]),'upper':list(map(str,implication(rows,r))),
                        'lower':list(map(str,implication(rows,tuple(-x for x in r))))})
            normals.append(r[:2])
        if len(normals)==2-d:break
    assert len(normals)==2-d
    return {'kind':'dimension','dimension':d,'lower_points':[list(map(str,p)) for p in chosen],'equalities':eqs}

def box(x0,x1,y0,y1):
    return [tuple(map(F,r)) for r in [(1,0,x1),(-1,0,-x0),(0,1,y1),(0,-1,-y0)]]

def bounded(rows,b=100):return box(-b,b,-b,b)+[tuple(map(F,r)) for r in rows]
