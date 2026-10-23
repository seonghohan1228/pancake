import pandas as pd, numpy as np
f="/mnt/user-data/uploads/Bearing Test Rig/logs/live-daq_20261020_213732_047904700.csv"
d=pd.read_csv(f); d=d[(d.pressure_available==1)&(d.load_cell_available==1)].copy()
d=d[['elapsed_s','timestamp_utc','pressure_bar_abs','ctt_v','force_n','tc0_c','load_raw_u16']].reset_index(drop=True)
print(len(d), d.elapsed_s.diff().describe())
# atmospheric reference
print(d[d.elapsed_s<55][['pressure_bar_abs','ctt_v']].describe().loc[['mean','std','min','max']])
print(d[d.elapsed_s>440][['pressure_bar_abs','ctt_v']].describe().loc[['mean','std','min','max']])
print(d[(d.elapsed_s>385)&(d.elapsed_s<430)][['pressure_bar_abs','ctt_v']].describe().loc[['mean','std','min','max']])
# rolling stability
w=15
d['p_sd']=d.pressure_bar_abs.rolling(w,center=True).std(); d['v_sd']=d.ctt_v.rolling(w,center=True).std()
st=(d.p_sd<0.02)&(d.v_sd<0.004)
d['seg']=(st!=st.shift()).cumsum()
rows=[]
for s,g in d[st].groupby('seg'):
    if g.elapsed_s.iloc[-1]-g.elapsed_s.iloc[0]<3: continue
    rows.append(dict(t0=g.elapsed_s.iloc[0],t1=g.elapsed_s.iloc[-1],n=len(g),p=g.pressure_bar_abs.mean(),psd=g.pressure_bar_abs.std(),v=g.ctt_v.mean(),vsd=g.ctt_v.std(),F_old=g.force_n.mean()))
P=pd.DataFrame(rows); pd.set_option('display.width',200); print(P.round(4))
P.to_csv('plateaus_raw.csv',index=False); d.to_csv('cal_clean.csv',index=False)
import numpy as np, json
prod=1.0027  # rod side: transducer reading at rest (user: 1 bar abs)
A=np.pi/4*0.100**2   # m2, 100 mm bore, full-bore extension
# loading envelope (virgin rise): slow ramp + peak-of-rise points (p_abs, V)
up=np.array([[1.50,0.0],[1.597,0.042],[1.731,0.078],[1.825,0.124],[1.934,0.172],[2.000,0.202],[2.131,0.254],
             [4.784,1.296],[6.349,1.814],[7.106,2.036],[7.993,2.307],[8.222,2.403],[9.440,2.748]])
# unloading envelope (trough of each fall)
dn=np.array([[1.0027,0.138],[1.081,0.272],[1.909,0.741],[2.687,1.064],[3.597,1.411],[4.849,1.866],[6.228,2.311],[6.668,2.412],[9.440,2.748]])
for a in (up,dn): a[:,0]=(a[:,0]-prod)*1e5*A   # -> A*dp in N
# up branch: F <= A dp (upper bound); dn: F >= A dp (lower bound)
V=np.linspace(0.15,2.74,60)
Fu=np.interp(V,up[:,1],up[:,0]); Fd=np.interp(V,dn[:,1],dn[:,0])
Fm=(Fu+Fd)/2; ff=(Fu-Fd)/2
for v in [0.15,0.3,0.5,1.0,1.44,2.0,2.5,2.74]:
    i=np.argmin(abs(V-v)); print(f'V={V[i]:.2f}  F_up(<=)={Fu[i]:7.0f}  F_dn(>=)={Fd[i]:7.0f}  mid={Fm[i]:7.0f}  friction={ff[i]:5.0f} N  K_mid={Fm[i]/V[i]:6.0f}  old={2752.246*(V[i]-0.013031):7.0f}')
# fit mid curve through origin: F=a V + b V^2
X=np.c_[V,V**2]; c,*_=np.linalg.lstsq(X,Fm,rcond=None); print('mid fit a,b',c, 'rms',np.sqrt(np.mean((X@c-Fm)**2)))
cl,*_=np.linalg.lstsq(V[:,None],Fm,rcond=None); print('linear through 0 K',cl,'rms',np.sqrt(np.mean((V*cl[0]-Fm)**2)))
cu,*_=np.linalg.lstsq(X,Fu,rcond=None); cd,*_=np.linalg.lstsq(X,Fd,rcond=None); print('bounds fits',cu,cd)
json.dump(dict(A_m2=A,p_rod_bar_abs=prod,up=up.tolist(),dn=dn.tolist(),mid=c.tolist(),upper=cu.tolist(),lower=cd.tolist(),lin_K=cl[0]),open('cal3.json','w'),indent=1)
import numpy as np, json
c=json.load(open('cal3.json')); up=np.array(c['up']); dn=np.array(c['dn'])
# exclude breakaway (V<0.3) and reversal zone at top (V>2.45 on dn)
u=up[(up[:,1]>0.25)]; d=dn[(dn[:,1]>0.3)&(dn[:,1]<2.45)]
Kmax=(u[:,0]/u[:,1]); Kmin=(d[:,0]/d[:,1])
print('up F/V',np.c_[u[:,1],Kmax].round(0)); print('dn F/V',np.c_[d[:,1],Kmin].round(0))
kmax=Kmax.min(); kmin=Kmin.max(); kmid=(kmax+kmin)/2
print('K range',kmin,kmax,'mid',kmid,'±%',(kmax-kmin)/2/kmid*100)
# implied seal friction at mid K
print('friction up',(u[:,0]-kmid*u[:,1]).round(0),'dn',(kmid*d[:,1]-d[:,0]).round(0))
# Sep 2 check: 1.4399 V at 'dp=5 bar'
print('Sep2 point implied K if dp=5:',3926.99/1.4399, ' dp needed on today loading branch for 1.44 V:', np.interp(1.4399,up[:,1],up[:,0])/785.398)
c.update(K_min=kmin,K_max=kmax,K_mid=kmid,V0=0.0); json.dump(c,open('cal3.json','w'),indent=1)
import numpy as np, pandas as pd, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral','Times New Roman'],'mathtext.fontset':'stix','font.size':10.5})
c=json.load(open('cal3.json')); Kold,V0old=2752.246065784855,0.013031204699778745
Ks={'mid':c['K_mid'],'lo':c['K_min'],'hi':c['K_max']}
R,cl=0.0215,42.5e-6
RB=np.array([[10,2*np.pi**2*10],[1.33,26.4],[0.631,12.8],[0.264,5.79],[0.121,3.22],[0.0446,1.70],[0.0188,1.05],[0.00474,0.514]])[::-1]
g=lambda S: np.exp(np.interp(np.log(S),np.log(RB[:,0]),np.log(RB[:,1])))
W=lambda Wo,K: K*(Wo/Kold+V0old)
P=pd.read_pickle('P.pkl'); C=pd.read_pickle('C.pkl'); M=pd.read_pickle('M.pkl')
for k,K in Ks.items():
    r=W(C.load_n,K)/C.load_n; C['r_'+k]=r
    C['tauCFD_'+k]=C.cfd_torque_nm*r*g(C.S_TC0/r)/g(C.S_TC0)
    C['err_'+k]=100*(C['tauCFD_'+k]/C.new_exp_torque_nm-1)
C['exp_f_mid']=C.new_exp_f/C.r_mid; C['S2_mid']=C.S_TC2/C.r_mid; C['cfd_f_mid']=C.tauCFD_mid/(C.load_n*C.r_mid*R)
tab=C[['case_id','rpm','load_n','r_mid','r_lo','r_hi','S_TC2','S2_mid','new_exp_f','exp_f_mid','new_exp_torque_nm','cfd_torque_nm','tauCFD_mid','new_relative_error_percent','err_mid','err_lo','err_hi']]
pd.set_option('display.width',250); print(tab.round(4).to_string())
tab.to_csv('LGE-2026-10-20-recal-case-errors.csv',index=False)
# ---------- Figure 1: calibration
up=np.array(c['up']); dn=np.array(c['dn']); d=pd.read_csv('cal_clean.csv'); A=c['A_m2']
fig,ax=plt.subplots(figsize=(6.4,5.2))
ax.plot(d.ctt_v,(d.pressure_bar_abs-c['p_rod_bar_abs'])*1e5*A,'.',color='0.75',ms=2,label=r'Record 21:37, all samples ($A\Delta p$)')
ax.plot(up[:,1],up[:,0],'^',color='C3',ms=6,label=r'Loading envelope, $F\leq A\Delta p$')
ax.plot(dn[:,1],dn[:,0],'v',color='C0',ms=6,label=r'Unloading envelope, $F\geq A\Delta p$')
V=np.linspace(0,2.8,10)
ax.plot(V,Kold*(V-V0old),'k--',lw=1.3,label=r'Old: $2752\,(V-0.013)$')
ax.plot(V,Ks['mid']*V,'g-',lw=1.8,label=rf'Proposed: ${Ks["mid"]:.0f}\,V$')
ax.fill_between(V,Ks['lo']*V,Ks['hi']*V,color='g',alpha=0.15,label=rf'Bounds ${Ks["lo"]:.0f}$–${Ks["hi"]:.0f}\ \mathrm{{N/V}}$')
ax.set_xlabel(r'Amplifier output $V_\mathrm{CTT}$ (V)'); ax.set_ylabel(r'Force (N)'); ax.set_xlim(0,2.85); ax.set_ylim(-200,8000)
ax.grid(alpha=0.3); ax.legend(fontsize=8.5,loc='upper left'); ax.set_title(r'Load-cell recalibration, Ø100 mm piston, rod side 1 bar abs',fontsize=10.5)
plt.tight_layout(); plt.savefig('LGE-2026-10-20-loadcell-recal.png',dpi=200)
# ---------- Figure 2: Stribeck + error
fig,ax=plt.subplots(1,2,figsize=(12,5.0),gridspec_kw={'width_ratios':[1.35,1]})
a=ax[0]; s=P.S_TC2; f=P.new_f
a.loglog(s,f,'o',mfc='none',mec='0.6',ms=4,label='FW68D exp., old load scale')
a.loglog(s/P.r_mid,f/P.r_mid,'o',color='C0',ms=4,alpha=0.8,label='FW68D exp., recalibrated')
Ms=M.S;Mf=M.f
a.loglog(Ms,Mf,'D',mfc='none',mec='C1',ms=5,label='FW68D–R410A exp., old')
a.loglog(Ms/M.r_mid,Mf/M.r_mid,'D',color='C1',ms=5,label='FW68D–R410A exp., recalibrated')
a.loglog(C.S_TC2,C.cfd_f,'s',mfc='none',mec='k',ms=8,label='CFD, old load')
a.loglog(C.S2_mid,C.cfd_f_mid,'s',color='k',ms=6,label='CFD, recalibrated load (R–B scaling)')
for _,r in C.iterrows():
    a.annotate('',xy=(r.S2_mid,r.exp_f_mid),xytext=(r.S_TC2,r.new_exp_f),arrowprops=dict(arrowstyle='->',color='C3',lw=1.2))
    a.plot([r.S_TC2/r.r_lo,r.S_TC2/r.r_hi],[r.new_exp_f/r.r_lo,r.new_exp_f/r.r_hi],'-',color='C3',lw=3,alpha=0.5)
    a.annotate(r.case_id,(r.S2_mid,r.exp_f_mid),xytext=(5,-11),textcoords='offset points',fontsize=10,color='C3')
SS=np.logspace(-2.2,0.4,50); a.loglog(SS,2*np.pi**2*SS/(R/cl),':',color='0.4',label='Petroff')
a.set_xlabel(r'Sommerfeld number $S=(\mu N/P)(R/c)^2$ (–), TC2 viscosity'); a.set_ylabel(r'Friction coefficient $f=\tau/(WR)$ (–)')
a.set_xlim(0.008,2.2); a.set_ylim(0.002,0.06); a.grid(which='both',alpha=0.25); a.legend(fontsize=8,loc='upper left')
a.set_title('(a) Stribeck plane: points move along slope-one lines',fontsize=10.5)
b=ax[1]; x=np.arange(6); 
b.bar(x-0.19,C.new_relative_error_percent,0.36,color='0.7',label='Old load scale')
b.bar(x+0.19,C.err_mid,0.36,color='C0',label='Recalibrated (estimate)')
b.errorbar(x+0.19,C.err_mid,yerr=[C.err_mid-np.minimum(C.err_lo,C.err_hi),np.maximum(C.err_lo,C.err_hi)-C.err_mid],fmt='none',ecolor='k',capsize=3)
for i,(o,n) in enumerate(zip(C.new_relative_error_percent,C.err_mid)):
    b.text(i+0.19,n+(2 if n>0 else -2),f'{n:+.0f}',ha='center',va='bottom' if n>0 else 'top',fontsize=8.5)
    b.text(i-0.19,o+(2 if o>0 else -2),f'{o:+.0f}',ha='center',va='bottom' if o>0 else 'top',fontsize=8.5,color='0.35')
b.axhline(0,color='k',lw=0.8); b.set_xticks(x); b.set_xticklabels([f'{r.case_id}\n{r.rpm:.0f} rpm' for _,r in C.iterrows()])
b.set_ylabel(r'Torque error $100(\tau_\mathrm{CFD}/\tau_\mathrm{exp}-1)$ (%)'); b.set_ylim(-62,50); b.grid(axis='y',alpha=0.3); b.legend(fontsize=8.5,loc='upper left')
b.set_title('(b) CFD torque error, no new CFD',fontsize=10.5)
plt.tight_layout(); plt.savefig('LGE-2026-10-20-recal-stribeck-cfd-error.png',dpi=200)
json.dump(dict(K_mid=Ks['mid'],K_lo=Ks['lo'],K_hi=Ks['hi'],ratio_points=[P.r_mid.min(),P.r_mid.max()],ratio_mix=[M.r_mid.min(),M.r_mid.max()]),open('summary.json','w'))
