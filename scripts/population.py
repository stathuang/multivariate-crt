"""Compute the numerical population approximations in Supplement E.5 and E.8.
B and tau are Monte Carlo expectations using exactly the original 500 backgrounds.
"""
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS']: os.environ[k]='1'
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/e5-max-coefficient-t-20261103/src'))
from inputs import reconstruct,old_records,conditional_parameters

def scan(theta,L):
    j=np.eye(len(theta))-theta/np.diag(theta)[:,None]
    forward=np.linalg.solve(np.eye(len(j))-np.tril(j,-1),np.triu(j,1))
    reverse=np.linalg.solve(np.eye(len(j))-np.triu(j,1),np.tril(j,-1))
    return np.linalg.matrix_power(reverse,L)@np.linalg.matrix_power(forward,L)

out={}
for design in ['A','B']:
    records=[r for r in old_records() if r['design']==design and r['pattern']=='dense' and r['theta']==.75]
    assert len(records)==500
    B=0.;tau=0.
    for r in records:
        x,y,t,s,omega,graph,key=reconstruct(r)
        S=x[:,s]-x[:,s].mean(0); n=len(x)
        H=np.linalg.solve(S@S.T+.2*n*np.eye(n),.2*n*np.eye(n))
        # Full centered representation has one extra constant direction; remove it from trace.
        Q=H@H
        B=B+S.T@Q@S/n/500
        tau=tau+(np.trace(Q)-1)/n/500
    _,C,V=conditional_parameters(omega,t)
    direction=np.ones(len(t))/np.sqrt(len(t));direction/=np.sqrt(direction@V@direction)
    beta=.75*direction;gamma=np.zeros(len(s));gamma[:4]=[1,-.7,.5,.3]
    # gamma is the coefficient of S in the conditional representation Y=S gamma+E beta+epsilon.
    m=C@B@gamma;u=tau*V@beta;D=tau*V+C@B@C.T+np.eye(len(t))/n
    F=scan(omega[np.ix_(t,t)],3)
    target=float(u@np.linalg.solve(D,u)-(F@u)@np.linalg.solve(D,F@u))
    cross=float(2*m@np.linalg.solve(D,(np.eye(len(t))-F)@u))
    effects=np.sqrt(n)*np.abs(beta)/np.sqrt(np.diag(omega)[t])
    out[design]=dict(target=target,background_interaction=cross,total=target+cross,
                    conditional_effect_min=float(effects.min()),conditional_effect_max=float(effects.max()),
                    group_signal_variance=float(beta@V@beta),backgrounds=500)
    assert abs(target+cross-{'A':.040760044950,'B':.094472168835}[design])<1e-10,out[design]
theta=np.eye(8)
for j in range(7):theta[j,j+1]=theta[j+1,j]=-.2
F=scan(theta,1)
w,U=np.linalg.eigh(theta);sqrt=(U*np.sqrt(w))@U.T;inverse=(U*(1/np.sqrt(w)))@U.T
norm=float(np.linalg.norm(sqrt@F@inverse,2));assert abs(norm-.0598)<.00005
out['connected_scan_theta_norm']=norm
path=ROOT/'outputs/rebuilt/population.json';path.parent.mkdir(parents=True,exist_ok=True)
path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
