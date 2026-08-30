"""
Period-stratified accuracy from lai_val_samples.csv (answers Reviewer 1: accuracy in the
non-training periods 2000-2005 and 2019-2022 vs the 2006-2018 training window).
Streams the CSV in chunks (69M rows) and accumulates moments per (period, class) -> R2/RMSE/MAE/bias.
class: RETRIEVED (qa==0, value=retrieved) | GAP-FILLED (qa 1-4, value=overall).
"""
import numpy as np, pandas as pd, collections
CSV = r"K:\Hangkai\CONUS_LAI\lai_val_samples.csv"
SCALE = 0.001
acc = collections.defaultdict(lambda: np.zeros(7))   # n, Sa, Sb, Saa, Sbb, Sab, Sabsd

def add(key, a, b):
    if a.size == 0: return
    v = acc[key]
    v[0]+=a.size; v[1]+=a.sum(); v[2]+=b.sum(); v[3]+=(a*a).sum()
    v[4]+=(b*b).sum(); v[5]+=(a*b).sum(); v[6]+=np.abs(a-b).sum()

for ch in pd.read_csv(CSV, chunksize=2_000_000):
    mod = ch['modis'].to_numpy(float)*SCALE
    good = mod >= 0
    per = ch['period'].to_numpy(); qa = ch['qa'].to_numpy(int)
    ret = ch['retrieved'].to_numpy(float)*SCALE; ov = ch['overall'].to_numpy(float)*SCALE
    for p in ('2000-2005','2006-2018','2019-2022'):
        pm = good & (per==p)
        add((p,'RETRIEVED'), ret[pm & (qa==0)],            mod[pm & (qa==0)])
        add((p,'GAP-FILLED'), ov[pm & (qa>=1) & (qa<=4)],  mod[pm & (qa>=1) & (qa<=4)])

def stats(v):
    n,Sa,Sb,Saa,Sbb,Sab,Sabsd = v
    if n < 5: return (int(n), float('nan'), float('nan'), float('nan'), float('nan'))
    bias=(Sa-Sb)/n; rmse=np.sqrt(max(0,(Saa-2*Sab+Sbb)/n)); mae=Sabsd/n
    num=(n*Sab-Sa*Sb)**2; den=(n*Saa-Sa*Sa)*(n*Sbb-Sb*Sb)
    r2=num/den if den>0 else float('nan')
    return (int(n), r2, rmse, mae, bias)

print(f"{'period':<12}{'class':<12}{'N':>12}{'R2':>7}{'RMSE':>7}{'MAE':>7}{'bias':>7}")
for p in ('2000-2005','2006-2018','2019-2022'):
    for c in ('RETRIEVED','GAP-FILLED'):
        n,r2,rm,ma,bi = stats(acc[(p,c)])
        print(f"{p:<12}{c:<12}{n:>12,}{r2:>7.2f}{rm:>7.2f}{ma:>7.2f}{bi:>7.2f}")
    print()
print("Reviewer 1: compare 2000-2005 & 2019-2022 (outside RF training) vs 2006-2018 (training).")
