"""Independent spectral check; output location made portable for this bundle."""
from pathlib import Path
import numpy as np,json
from scipy.linalg import eigh
n=1024;theta=2*np.pi*np.arange(n)/n;kap=.8;a=1.4;b=.7;speed=np.sqrt(a*a*np.sin(theta)**2+b*b*np.cos(theta)**2);speed_ratio=(a*a-b*b)*np.sin(theta)*np.cos(theta)/(speed*speed)
rho=np.exp(kap*np.cos(theta));rho/=np.sum(rho*speed)*2*np.pi/n;mu=rho*speed*2*np.pi/n
phi=np.sin(2*theta)+.3*np.cos(3*theta);dphi=2*np.cos(2*theta)-.9*np.sin(3*theta);ddphi=-4*np.sin(2*theta)-2.7*np.cos(3*theta)
h=(-ddphi+(kap*np.sin(theta)+speed_ratio)*dphi)/(speed*speed);h-=np.dot(mu,h)
g0=np.dot(mu,(dphi/speed)**2);curv=a*b/speed**3;bias=-.25*np.dot(mu,curv**2*(dphi/speed)**2)
x=a*np.cos(theta);y=b*np.sin(theta);c=(x[:,None]-x[None,:])**2+(y[:,None]-y[None,:])**2;output=[]
for eps in [.1,.05,.025,.0125,.00625]:
 T=np.exp(-c/eps)/np.sqrt(np.pi*eps);v=1/np.sqrt(rho)
 for it in range(100):
  residual=v*(T@(mu*v))-1
  if max(abs(residual))<3e-14:break
  v/=np.sqrt(1+residual)
 S=(np.sqrt(mu)*v)[:,None]*T*(np.sqrt(mu)*v)[None,:]
 vals,V=eigh(S,check_finite=False);p=V.T@(np.sqrt(mu)*h);vals=vals[:-1];p=p[:-1]
 g=eps/2*np.sum(p*p*vals/(1-vals*vals))
 row={'eps':eps,'g':float(g),'g0':float(g0),'predicted_bias':float(bias),'measured_bias':float((g-g0)/eps),'remainder_scaled':float((g-g0-eps*bias)/eps**2),'scaling_residual':float(max(abs(residual)))}
 output.append(row);print(row,flush=True)
open(Path(__file__).resolve().parents[1]/'results'/('independent_ellipse_check.json'),'w').write(json.dumps(output,indent=2))
