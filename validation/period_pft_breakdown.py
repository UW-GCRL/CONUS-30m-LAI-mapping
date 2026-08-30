"""
Period x PFT x class accuracy from lai_val_samples.csv. Shows whether any biome degrades in
the non-training periods (2000-2005, 2019-2022) vs training (2006-2018). Streams the CSV.
"""
import numpy as np, pandas as pd, collections
CSV = r"K:\Hangkai\CONUS_LAI\lai_val_samples.csv"
SCALE = 0.001
PFT = {1:'decid',2:'evergreen',3:'mixed',4:'shrub',5:'grass',6:'crop',7:'woodywet',8:'herbwet'}
PERIODS = ('2000-2005','2006-2018','2019-2022')
acc = collections.defaultdict(lambda: np.zeros(7))     # n,Sa,Sb,Saa,Sbb,Sab,Sabsd

def add(key, a, b):
    if a.size == 0: return
    v = acc[key]; v[0]+=a.size; v[1]+=a.sum(); v[2]+=b.sum()
    v[3]+=(a*a).sum(); v[4]+=(b*b).sum(); v[5]+=(a*b).sum(); v[6]+=np.abs(a-b).sum()

for ch in pd.read_csv(CSV, chunksize=2_000_000):
    mod = ch['modis'].to_numpy(float)*SCALE; good = mod >= 0
    per = ch['period'].to_numpy(); pft = ch['pft'].to_numpy(int); qa = ch['qa'].to_numpy(int)
    ret = ch['retrieved'].to_numpy(float)*SCALE; ov = ch['overall'].to_numpy(float)*SCALE
    isret = qa==0; isfill = (qa>=1)&(qa<=4)
    for p in PERIODS:
        pm = good & (per==p)
        for k in PFT:
            base = pm & (pft==k)
            add((p,k,'R'), ret[base & isret], mod[base & isret])
            add((p,k,'F'), ov[base & isfill], mod[base & isfill])

def stats(v):
    n,Sa,Sb,Saa,Sbb,Sab,Sabsd = v
    if n < 20: return (int(n), float('nan'), float('nan'), float('nan'))
    bias=(Sa-Sb)/n; rmse=np.sqrt(max(0,(Saa-2*Sab+Sbb)/n))
    num=(n*Sab-Sa*Sb)**2; den=(n*Saa-Sa*Sa)*(n*Sbb-Sb*Sb); r2=num/den if den>0 else float('nan')
    return (int(n), r2, rmse, bias)

for cls,lbl in [('R','RETRIEVED (QA=0)'),('F','GAP-FILLED (QA 1-4)')]:
    print(f"\n===== {lbl}: R2 (RMSE) by PFT x period =====")
    print(f"{'PFT':<10} " + " ".join(f"{p:>18}" for p in PERIODS))
    for k in PFT:
        cells=[]
        for p in PERIODS:
            n,r2,rm,bi = stats(acc[(p,k,cls)])
            cells.append(f"{r2:4.2f} ({rm:4.2f}) n{n//1000}k" if n>=20 else "   -   ")
        print(f"{PFT[k]:<10} " + " ".join(f"{c:>18}" for c in cells))
print("\n(each cell: R2 (RMSE)  n=thousands;  compare across the 3 periods per PFT)")
