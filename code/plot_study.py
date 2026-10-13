#!/usr/bin/env python3
"""Render PSS R8 real-source descriptive figures; no inferential error bars."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((ROOT/'results/full_maps_x_4_clearance_scenarios.csv').open()))
series={mapname:sorted([q for q in rows if q['map']==mapname],key=lambda z:int(z['clearance_cells'])) for mapname in sorted(set(z['map'] for z in rows))}
figout=ROOT/'figures';figout.mkdir(exist_ok=True,parents=True)
for field,caption,slug in [('beta0_closed','Connected components of the closed-cell union, beta0','beta0'),
                           ('beta1_closed','Independent holes in the closed-cell union, beta1','beta1'),
                           ('free_cells','Admitted free cells after native-pixel square clearance','freecells')]:
 fig,ax=plt.subplots(figsize=(7.4,4.4),layout='constrained')
 for name,ds in series.items():
  xs=[int(x['clearance_cells']) for x in ds]
  ys=[int(x[field]) for x in ds]
  ax.plot(xs,ys,marker='o',linewidth=1.9,label=name.replace('_256.map',''))
 ax.set(xlabel='Square-Chebyshev clearance (native map cells, NOT metres)',ylabel=caption)
 ax.set_xticks([0,1,2,3]);ax.grid(True,alpha=.2);ax.legend(title='MovingAI source map',loc='best',fontsize=8)
 fig.savefig(figout/('R8_fullmap_'+slug+'.png'),dpi=220,bbox_inches='tight')
 plt.close(fig)
print('FIGURES',sorted(q.name for q in figout.iterdir()))
# A single legible source-specific operational-screening figure.
import numpy as np
from matplotlib.patches import Rectangle
windows=[]
for r in (0,1):
 src=ROOT/'examples'/f'Boston_0_256_x108_y30_f2_clear{r}'
 with (src/'admitted_window.csv').open() as f:
  data=np.asarray([[int(x) for x in row] for row in csv.reader(f)],dtype=int)
 info=__import__('json').loads((src/'case_result.json').read_text())['topology']
 windows.append((r,data,info))
fig,axs=plt.subplots(1,2,figsize=(8.4,4.05),layout='constrained')
for ax,(r,arr,info) in zip(axs,windows):
 ax.imshow(arr,cmap='Greys_r',origin='upper',vmin=0,vmax=1,interpolation='none')
 ax.set_xticks(np.arange(-.5,8,1),minor=True);ax.set_yticks(np.arange(-.5,8,1),minor=True)
 ax.grid(which='minor',color='#ababab',linewidth=.65);ax.tick_params(which='both',left=False,bottom=False,labelleft=False,labelbottom=False)
 ax.set_title(f'Scenario r = {r} native cell(s)\n{info["free_cells"]} admitted coarse cells, β₀={info["beta0_closed"]}, β₁={info["beta1_closed"]}\n4-neighbour β₀={info["beta0_4neighbor_gridgraph"]}',size=11)
fig.text(.5,.008,'Same Boston source, same 8 × 8 window (coarse origin x=108, y=30); scenario-only erosion before cropping.',ha='center',fontsize=9)
fig.savefig(figout/'R8_Boston_window_clearance_comparison.png',dpi=220,bbox_inches='tight')
plt.close(fig)
