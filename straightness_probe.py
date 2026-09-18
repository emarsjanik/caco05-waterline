import csv, glob, sys, numpy as np
from collections import defaultdict

def load(path):
    c,r,s=[],[],[]
    for x in csv.DictReader(open(path)):
        c.append(int(x['Column']))
        r.append(float(x.get('Row_Precise') or x['Row']))
        s.append(int(float(x.get('Has_Signal') or 1)))
    return np.array(c), np.array(r), np.array(s, dtype=bool)

def residual_rms(rows, win):
    """RMS deviation from a local straight-line fit, per column."""
    n=len(rows); out=np.full(n, np.nan)
    half=win//2
    x=np.arange(win, dtype=float)
    for i in range(half, n-half):
        seg=rows[i-half:i-half+win]
        if len(seg)<win: continue
        A=np.vstack([x, np.ones(win)]).T
        coef,_,_,_=np.linalg.lstsq(A, seg, rcond=None)
        out[i]=np.sqrt(np.mean((seg-(A@coef))**2))
    return out

WIN=int(sys.argv[2]) if len(sys.argv)>2 else 61
files=sorted(glob.glob(sys.argv[1]))[:40]
good=[]; bad=[]
for f in files:
    c,r,s=load(f)
    if len(c)<200: continue
    res=residual_rms(r, WIN)
    ok=np.isfinite(res)
    good.append(res[ok & s]); bad.append(res[ok & ~s])
good=np.concatenate([g for g in good if len(g)]) if good else np.array([])
bad=np.concatenate([b for b in bad if len(b)]) if bad else np.array([])
print('window = %d columns, files = %d' % (WIN, len(files)))
print()
for name,v in [('Has_Signal=1 (kept)',good), ('Has_Signal=0 (flagged)',bad)]:
    if len(v)==0: print('%-24s none'%name); continue
    print('%-24s n=%7d  median %7.3f  p10 %7.3f  p90 %7.3f  frac<0.05: %5.1f%%'
          % (name, len(v), np.median(v), np.percentile(v,10), np.percentile(v,90),
             100*np.mean(v<0.05)))
