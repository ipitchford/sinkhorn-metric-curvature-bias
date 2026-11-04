"""Independent spectral check; output location made portable for this bundle."""
from pathlib import Path
import numpy as np
from scipy.special import i0
from scipy.linalg import eigh
import json
n=768;theta=2*np.pi*np.arange(n)/n;kap=1.3;rho=np.exp(kap*np.cos(theta))/(2*np.pi*i0(kap));mu=rho*2*np.pi/n
phi=np.sin(2*theta)+.3*np.cos(3*theta);dphi=2*np.cos(2*theta)-.9*np.sin(3*theta)
h=4*np.sin(2*theta)+2.7*np.cos(3*theta)+kap*np.sin(theta)*dphi
h-=np.dot(mu,h);g0=np.dot(mu,dphi*dphi)
c=2*(1-np.cos(theta[:,None]-theta[None,:]));output=[]
for eps in [.1,.05,.025,.0125,.00625]:
 T=np.exp(-c/eps)/np.sqrt(np.pi*eps);v=1/np.sqrt(rho)
 for it in range(100):
  residual=v*(T@(mu*v))-1
  if max(abs(residual))<3e-14:break
  v/=np.sqrt(1+residual)
 S=(np.sqrt(mu)*v)[:,None]*T*(np.sqrt(mu)*v)[None,:]
 vals,V=eigh(S,check_finite=False);p=V.T@(np.sqrt(mu)*h);vals=vals[:-1];p=p[:-1]
 g=eps/2*np.sum(p*p*vals/(1-vals*vals))
 row={'eps':eps,'g':float(g),'g0':float(g0),'linear_scaled':float((g/g0-1)/eps),'remainder_scaled':float((g-g0+eps*g0/4)/eps**2),'scaling_residual':float(max(abs(residual))),'iterations':it}
 output.append(row);print(row,flush=True)
open(Path(__file__).resolve().parents[1]/'results'/('independent_circle_check.json'),'w').write(json.dumps(output,indent=2))
