"""
SCALE-MATCHED validation: aggregate 30 m LAI to the 500 m MODIS grid (mean of the 30 m pixels
whose centers fall in each MODIS cell), then compare cell-mean vs MODIS-median. Handles the
non-integer 30m:500m ratio by exact coordinate binning. Reports per PFT and per period, for
RETRIEVED (observed 30m only, qa==0) and GAP-FILLED (all valid gf), plus a homogeneous-cell
variant (within-cell CV of 30m LAI below a threshold).
Run: python agg_validate.py
"""
import os, csv, numpy as np, rasterio, collections
from rasterio.windows import Window
from affine import Affine

SRC=r"K:\Hangkai\CONUS_LAI"; GF=r"K:\Hangkai\CONUS_LAI_gapfilled"; VAL=r"K:\Hangkai\CONUS_LAI_gapfilled\MODIS"
MON=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
YEARS=range(2000,2023); NODATA=-32768; SCALE=0.001
WIN5=48; WINS_PER_MONTH=12; MINPIX=50; CV_PURE=0.30
PFT={1:'decid',2:'evergreen',3:'mixed',4:'shrub',5:'grass',6:'crop',7:'woodywet',8:'herbwet'}
rng=np.random.default_rng(0)
def period(y): return '2000-2005' if y<=2005 else ('2019-2022' if y>=2019 else '2006-2018')
acc=collections.defaultdict(lambda: np.zeros(7))   # (period,pft,class,pure) -> n,Sa,Sb,Saa,Sbb,Sab,Sabsd
def add(key,a,b):
    if a.size==0: return
    v=acc[key]; v[0]+=a.size; v[1]+=a.sum(); v[2]+=b.sum(); v[3]+=(a*a).sum(); v[4]+=(b*b).sum(); v[5]+=(a*b).sum(); v[6]+=np.abs(a-b).sum()

with rasterio.open(os.path.join(SRC,"2015_Jul_LAI.tif")) as r: tr30=r.transform; W30=r.width; H30=r.height

def cellmeans(gd,qd,sd,mi,mds,c5,r5,w5,h5):
    """aggregate 30m -> the 500m window (c5,r5,w5,h5); returns cell-mean gf, ret, cv, and mask."""
    t5=mds.transform
    # geo extent of the 500m window -> 30m pixel bounds
    xs=[t5.c+c5*t5.a, t5.c+(c5+w5)*t5.a]; ys=[t5.f+r5*t5.e, t5.f+(r5+h5)*t5.e]
    inv30=~tr30
    cols=[(inv30*(x,y))[0] for x in xs for y in ys]; rows=[(inv30*(x,y))[1] for x in xs for y in ys]
    cc0=max(0,int(np.floor(min(cols)))); cc1=min(W30,int(np.ceil(max(cols))))
    rr0=max(0,int(np.floor(min(rows)))); rr1=min(H30,int(np.ceil(max(rows))))
    if cc1<=cc0 or rr1<=rr0: return None
    w30=Window(cc0,rr0,cc1-cc0,rr1-rr0)
    gf=gd.read(1,window=w30); qa=qd.read(1,window=w30); ret=sd.read(1,window=w30)
    jj,ii=np.meshgrid(np.arange(cc1-cc0),np.arange(rr1-rr0))
    x=tr30.c+(cc0+jj+0.5)*tr30.a; y=tr30.f+(rr0+ii+0.5)*tr30.e
    inv5=~t5
    c5f=np.floor(inv5.a*x+inv5.b*y+inv5.c).astype(int)-c5
    r5f=np.floor(inv5.d*x+inv5.e*y+inv5.f).astype(int)-r5
    inside=(c5f>=0)&(c5f<w5)&(r5f>=0)&(r5f<h5)
    lin=(r5f*w5+c5f); N=w5*h5
    def means(mask,val):
        m=inside&mask
        s=np.bincount(lin[m],weights=val[m].astype(float),minlength=N)
        c=np.bincount(lin[m],minlength=N)
        return s,c
    s_gf,c_gf=means(gf>0,gf); ss_gf,_=means(gf>0,gf.astype(float)**2)
    s_rt,c_rt=means((qa==0)&(ret>0),ret)
    mean_gf=np.where(c_gf>=MINPIX,s_gf/np.maximum(c_gf,1),np.nan)
    mean_rt=np.where(c_rt>=MINPIX,s_rt/np.maximum(c_rt,1),np.nan)
    var=np.where(c_gf>=MINPIX,ss_gf/np.maximum(c_gf,1)-(s_gf/np.maximum(c_gf,1))**2,np.nan)
    cv=np.sqrt(np.maximum(var,0))/np.maximum(mean_gf,1e-6)
    return mean_gf.reshape(h5,w5),mean_rt.reshape(h5,w5),cv.reshape(h5,w5)

def main():
    for y in YEARS:
        mp=os.path.join(VAL,f"MODISmed_{y}.tif"); pp=os.path.join(VAL,f"PFT_{y}.tif")
        if not(os.path.exists(mp) and os.path.exists(pp)): print(f"{y}: skip"); continue
        mds=rasterio.open(mp); pds=rasterio.open(pp); W5=mds.width; H5=mds.height
        for mi,mon in enumerate(MON):
            rp=os.path.join(SRC,f"{y}_{mon}_LAI.tif"); gp=os.path.join(GF,f"{y}_{mon}_LAI_gf.tif"); qp=os.path.join(GF,f"{y}_{mon}_LAI_qa.tif")
            if not(os.path.exists(rp) and os.path.exists(gp) and os.path.exists(qp)): continue
            sd=rasterio.open(rp); gd=rasterio.open(gp); qd=rasterio.open(qp)
            for _ in range(WINS_PER_MONTH):
                c5=int(rng.integers(0,W5-WIN5)); r5=int(rng.integers(0,H5-WIN5))
                res=cellmeans(gd,qd,sd,mi,mds,c5,r5,WIN5,WIN5)
                if res is None: continue
                mgf,mrt,cv=res
                mod=mds.read(mi+1,window=Window(c5,r5,WIN5,WIN5)).astype(float)*SCALE
                pf=pds.read(1,window=Window(c5,r5,WIN5,WIN5))
                vmod=mod>=0
                for k in PFT:
                    kk=(pf==k)&vmod
                    gmask=kk&np.isfinite(mgf); rmask=kk&np.isfinite(mrt)
                    add((period(y),k,'F',0), mgf[gmask]*SCALE, mod[gmask])
                    add((period(y),k,'R',0), mrt[rmask]*SCALE, mod[rmask])
                    pure=gmask&(cv<CV_PURE)
                    add((period(y),k,'F',1), mgf[pure]*SCALE, mod[pure])
            sd.close(); gd.close(); qd.close()
        mds.close(); pds.close()
        print(f"{y} done", flush=True)
    report()

def stats(v):
    n,Sa,Sb,Saa,Sbb,Sab,Sd=v
    if n<20: return (int(n),float('nan'),float('nan'),float('nan'),float('nan'))
    bias=(Sa-Sb)/n; rmse=np.sqrt(max(0,(Saa-2*Sab+Sbb)/n)); mae=Sd/n
    num=(n*Sab-Sa*Sb)**2; den=(n*Saa-Sa*Sa)*(n*Sbb-Sb*Sb); r2=num/den if den>0 else float('nan')
    return (int(n),r2,rmse,mae,bias)

def pool(**filt):
    tot=np.zeros(7)
    for (p,k,c,pu),v in acc.items():
        if 'period' in filt and p!=filt['period']: continue
        if 'pft' in filt and k!=filt['pft']: continue
        if c!=filt['cls']: continue
        if pu!=filt.get('pure',0): continue
        tot+=v
    return tot

def row(name,v):
    n,r2,rm,ma,bi=stats(v); print(f"{name:<12}{n:>10,}{r2:>7.2f}{rm:>7.2f}{ma:>7.2f}{bi:>7.2f}")

def report():
    print("\n===== 500m-AGGREGATED vs MODIS-median (all periods pooled) =====")
    for cls,lbl in [('R','RETRIEVED (obs 30m cell-mean)'),('F','GAP-FILLED (all 30m cell-mean)')]:
        print(f"\n--- {lbl} ---"); print(f"{'PFT':<12}{'N':>10}{'R2':>7}{'RMSE':>7}{'MAE':>7}{'bias':>7}")
        for k in PFT: row(PFT[k], pool(pft=k,cls=cls))
        row('ALL', pool(cls=cls))
    print("\n===== GAP-FILLED by period (pooled PFT) =====")
    print(f"{'period':<12}{'N':>10}{'R2':>7}{'RMSE':>7}{'MAE':>7}{'bias':>7}")
    for p in ('2000-2005','2006-2018','2019-2022'): row(p, pool(period=p,cls='F'))
    print("\n===== HOMOGENEOUS cells (within-cell CV<0.30), GAP-FILLED =====")
    print(f"{'PFT':<12}{'N':>10}{'R2':>7}{'RMSE':>7}{'MAE':>7}{'bias':>7}")
    for k in PFT: row(PFT[k], pool(pft=k,cls='F',pure=1))
    row('ALL', pool(cls='F',pure=1))

if __name__=="__main__": main()
