import numpy as np
import pandas as pd
from scipy.optimize import least_squares

def phi_power(theta, eta, gamma):
    return eta/(theta**gamma*(1+theta)**(1-gamma))

def ssvi_total_variance(k, theta, rho, eta, gamma):
    phi=phi_power(theta,eta,gamma)
    return .5*theta*(1+rho*phi*k+np.sqrt((phi*k+rho)**2+1-rho**2))

def calibrate_ssvi(quotes):
    Tvals=np.sort(quotes["T"].unique()); n=len(Tvals)
    theta0=[]
    for T in Tvals:
        s=quotes[quotes.T==T]
        j=(s.log_moneyness.abs()).idxmin()
        theta0.append(float(s.loc[j,"true_iv"])**2*T)
    theta0=np.array(theta0)
    k=quotes.log_moneyness.to_numpy()
    T=quotes.T.to_numpy()
    wobs=quotes.true_iv.to_numpy()**2*T
    idx=np.searchsorted(Tvals,T)

    # Transforms enforce theta>0, -1<rho<1, eta>0, 0<gamma<1.
    def unpack(x):
        theta=np.exp(x[:n])
        rho=np.tanh(x[n])
        eta=np.exp(x[n+1])
        gamma=1/(1+np.exp(-x[n+2]))
        return theta,rho,eta,gamma

    def resid(x):
        theta,rho,eta,gamma=unpack(x)
        return ssvi_total_variance(k,theta[idx],rho,eta,gamma)-wobs

    x0=np.r_[np.log(theta0),np.arctanh(-.5),np.log(.6),0.]
    fit=least_squares(resid,x0,max_nfev=5000,
                      xtol=1e-13,ftol=1e-13,gtol=1e-13)
    theta,rho,eta,gamma=unpack(fit.x)
    return {"theta_by_T":dict(zip(Tvals,theta)),
            "rho":rho,"eta":eta,"gamma":gamma,
            "rmse_total_variance":float(np.sqrt(np.mean(fit.fun**2)))}

if __name__=="__main__":
    q=pd.read_csv("data/assumed_option_quotes.csv")
    x=calibrate_ssvi(q)
    print(f"rho={x['rho']:.6f}, eta={x['eta']:.6f}, gamma={x['gamma']:.6f}")
    print(f"RMSE(w)={x['rmse_total_variance']:.8f}")
    for T,th in x["theta_by_T"].items():
        print(f"T={T:.2f}: theta={th:.8f}, ATM vol={(th/T)**.5:.4%}")
