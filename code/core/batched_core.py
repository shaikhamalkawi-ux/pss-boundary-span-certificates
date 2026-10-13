"""Solver-free checker for batched rational witnesses and boundary-span omissions.

A positive batch is ONE point checked in every named branch. It certifies the
entire simplex on that support. Star boundaries generate its boundary image;
no absent triangle is silently filled. All branch pairs remain fully certified.
This module imports no geometry generator, optimizer, network, or external package.
"""
from __future__ import annotations
from itertools import combinations
from exact_core import (require, frac, canonical_hash, parse_rows, feasible,
    check_geometry, check_dimension, BinaryBasis, components,
    triangle_universe, boundary_vector, CertificateError)


def check_lower(system: dict, cert: dict) -> dict:
    require(system.get('schema') == 'PSS-H-union-1', 'unsupported system schema')
    require(cert.get('system_sha256') == canonical_hash(system), 'system hash mismatch')
    n = system.get('ambient_dimension')
    require(type(n) is int and n >= 1, 'invalid ambient dimension')
    raw = system['branches']
    require(isinstance(raw, list), 'branch array missing')
    m = len(raw)
    require(all(type(b.get('id')) is int for b in raw), 'noninteger branch ID')
    require([b['id'] for b in raw] == list(range(m)), 'branch IDs not canonical')
    rows = [parse_rows(b['rows'], n) for b in raw]
    br = cert['branches']
    require(len(br) == m, 'branch coverage')
    dims, vertices = {}, []
    for i, r in enumerate(br):
        require(type(r['id']) is int and r['id'] == i, 'branch record ID')
        if r.get('kind') == 'farkas':
            require(not check_geometry(r, rows[i], n), 'unexpected feasible branch')
        else:
            require(r.get('kind') == 'dimension', 'invalid branch evidence')
            dims[i] = check_dimension(r, rows[i], n)
            vertices.append(i)
    expected = set(combinations(vertices, 2))
    seen, edges = set(), []
    for r in cert['pairs']:
        ids = tuple(r['ids'])
        require(all(type(i) is int for i in ids), 'noninteger pair ID')
        require(ids in expected and ids not in seen, 'duplicate or invalid pair')
        seen.add(ids)
        if check_geometry(r, [a for i in ids for a in rows[i]], n):
            edges.append(ids)
    require(seen == expected, 'incomplete pair universe')
    edges.sort()
    cc = components(vertices, edges)
    return {'n': n, 'rows': rows, 'vertices': vertices, 'dims': dims,
            'edges': edges, 'index': {e: i for i, e in enumerate(edges)},
            'cc': cc, 'z': len(edges)-len(vertices)+len(cc)}


def star_columns(ids: tuple[int, ...], edge_index: dict) -> list[int]:
    """Anchor at the smallest support ID; span equals all triangle boundaries."""
    if len(ids) < 3:
        return []
    a = ids[0]
    return [boundary_vector((a, b, c), edge_index)
            for b, c in combinations(ids[1:], 2)]


def check(system: dict, cert: dict) -> dict:
    require(cert.get('schema') == 'PSS-batched-boundary-1', 'unsupported batch schema')
    mode = cert.get('mode')
    require(mode in ('full', 'rank'), 'unsupported proof mode')
    lower = check_lower(system, cert)
    n, rows, vs, index = lower['n'], lower['rows'], lower['vertices'], lower['index']
    candidates = triangle_universe(vs, lower['edges'])
    candset = set(candidates)
    basis = BinaryBasis()
    supports, seen_records = [], set()
    star_count = incidence = inequality_checks = 0
    for r in cert['positive_batches']:
        require(r.get('kind') == 'batch_point', 'invalid positive batch kind')
        ids = tuple(r['ids'])
        require(all(type(i) is int for i in ids), 'noninteger batch ID')
        require(len(ids) >= 3 and tuple(sorted(set(ids))) == ids, 'noncanonical batch support')
        require(set(ids) <= set(vs), 'inactive batch branch')
        p = tuple(frac(x) for x in r['point'])
        key = (ids, p)
        require(key not in seen_records, 'duplicate positive batch')
        seen_records.add(key)
        require(all(feasible(p, rows[i], n) for i in ids), 'invalid batched common point')
        require(all(e in index for e in combinations(ids, 2)), 'batch contradicts pair graph')
        columns = star_columns(ids, index)
        for v in columns:
            basis.add(v)
        supports.append(frozenset(ids))
        star_count += len(columns)
        incidence += len(ids)
        inequality_checks += sum(len(rows[i]) for i in ids)
    negatives = set()
    for r in cert['negative_triples']:
        t = tuple(r['ids'])
        require(all(type(i) is int for i in t), 'noninteger negative ID')
        require(t in candset and t not in negatives, 'duplicate or invalid negative triple')
        require(r.get('kind') == 'farkas', 'negative record is not Farkas')
        require(not check_geometry(r, [a for i in t for a in rows[i]], n), 'invalid negative')
        require(not any(set(t) <= s for s in supports), 'conflicting batch and negative')
        negatives.add(t)
    positive_count = unresolved = 0
    for t in candidates:
        if t in negatives:
            continue
        positive = any(set(t) <= s for s in supports)
        positive_count += int(positive)
        unresolved += int(not positive)
        if mode == 'full':
            require(positive, 'full mode missing triangle membership')
        else:
            require(basis.reduce(boundary_vector(t, index)) == 0,
                    'unresolved rank-relevant triangle')
    dims, cc = lower['dims'], lower['cc']
    pss = [int(bool(vs)), max(dims.values()) if dims else None,
           sum(max(dims[i] for i in c) == 0 for c in cc), len(cc), lower['z']-basis.rank]
    require(cert['claimed_pss1'] == pss, 'claimed signature mismatch')
    return {'status': 'PASS', 'pss1': pss, 'active_branches': len(vs),
            'pairs': len(cert['pairs']), 'edges': len(lower['edges']),
            'triangle_candidates': len(candidates), 'boundary_rank': basis.rank,
            'positive_batches': len(supports), 'positive_branch_incidences': incidence,
            'positive_inequality_checks': inequality_checks,
            'star_columns_generated': star_count,
            'negative_triples': len(negatives), 'inferred_positive_triangles': positive_count,
            'unresolved_span_redundant': unresolved}


def expand_for_original_checker(system: dict, cert: dict) -> dict:
    """Expand each batch to original per-triangle point evidence for a cross-check.
    The old checker independently validates every expanded geometry record.
    """
    lower = check_lower(system, cert)
    triples = list(cert['negative_triples'])
    negative = {tuple(r['ids']) for r in triples}
    for t in triangle_universe(lower['vertices'], lower['edges']):
        if t in negative:
            continue
        group = next((r for r in cert['positive_batches'] if set(t) <= set(r['ids'])), None)
        if group is not None:
            triples.append({'kind': 'point', 'ids': list(t), 'point': group['point']})
    return {'schema': 'PSS-boundary-span-1', 'system_sha256': canonical_hash(system),
            'branches': cert['branches'], 'pairs': cert['pairs'],
            'triples': triples, 'claimed_pss1': cert['claimed_pss1']}
