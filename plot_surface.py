import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from calibrate_ssvi import calibrate_ssvi, ssvi_total_variance

q=pd.read_csv("data/assumed_option_quotes.csv")
fit=calibrate_ssvi(q)
Ts=np.linspace(q["T"].min(),q["T"].max(),60)
ks=np.linspace(-.30,.30,100)
nodes=np.array(sorted(fit["theta_by_T"]))
ths=np.array([fit["theta_by_T"][x] for x in nodes])
theta=np.interp(Ts,nodes,ths)

W=np.array([ssvi_total_variance(ks,t,fit["rho"],fit["eta"],fit["gamma"])
            for t in theta])
IV=np.sqrt(W/Ts[:,None])

fig=plt.figure(figsize=(10,7))
ax=fig.add_subplot(111,projection="3d")
K,TT=np.meshgrid(ks,Ts)
ax.plot_surface(K,TT,IV,alpha=.85)
ax.set_xlabel("log-moneyness k=log(K/F)")
ax.set_ylabel("T")
ax.set_zlabel("Implied volatility")
ax.set_title("Calibrated SSVI Surface")
plt.tight_layout(); plt.savefig("ssvi_surface.png",dpi=180)

fig,ax=plt.subplots(figsize=(10,6))
for T in nodes:
    s=q[q["T"]==T]
    ax.scatter(s.log_moneyness,s.true_iv,s=25,label=f"T={T:g}")
    kfine=np.linspace(-.3,.3,250)
    w=ssvi_total_variance(kfine,fit["theta_by_T"][T],
                          fit["rho"],fit["eta"],fit["gamma"])
    ax.plot(kfine,np.sqrt(w/T))
ax.set_xlabel("log-moneyness"); ax.set_ylabel("implied volatility")
ax.set_title("SSVI Smile Fits"); ax.legend(ncol=2,fontsize=8)
ax.grid(alpha=.25); plt.tight_layout(); plt.savefig("ssvi_smiles.png",dpi=180)
