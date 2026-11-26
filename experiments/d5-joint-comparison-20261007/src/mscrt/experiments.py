"""Pre-specified paired experiments; raw replicate records are always retained."""
from __future__ import annotations
import argparse, hashlib, json, platform, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import beta, f, norm
from .core import (RidgeStatistic, conditional_parameters, coordinate_pvalues,
                   dense_pvalue, g_copies, mvn_copies, oracle_copies, rank_pvalue,
                   ridge_spectrum_moments, spherical_pvalue, residual_basis)


def generator(seed: int, *keys: object) -> np.random.Generator:
    """Stable, method-specific streams, invariant to method order and job sharding."""
    digest = hashlib.sha256(json.dumps([seed, *keys], sort_keys=True).encode()).digest()
    return np.random.default_rng(np.frombuffer(digest[:16], dtype='<u4'))


def graph_design(p: int, t: int, kind: str, bandwidth: int, offdiag: float):
    omega = np.eye(p)
    for j in range(1, bandwidth + 1):
        omega += offdiag * (np.eye(p, k=j) + np.eye(p, k=-j))
    positions = np.arange(t) if kind == 'connected' else np.linspace(0, p - 1, t, dtype=int)
    if kind == 'dispersed' and np.any(np.diff(positions) <= bandwidth):
        raise ValueError('Dispersed targets must have no within-target edges.')
    order = np.r_[positions, np.setdiff1d(np.arange(p), positions)]
    omega = omega[np.ix_(order, order)]
    if np.linalg.eigvalsh(omega).min() <= 0:
        raise ValueError('Precision is not positive definite.')
    sigma = np.linalg.inv(omega)
    scales = np.sqrt(np.diag(sigma))
    # Standardize each marginal to variance one.
    omega = scales[:, None] * omega * scales[None, :]
    sigma = np.linalg.inv(omega)
    graph = abs(omega) > 1e-12; np.fill_diagonal(graph, False)
    return omega, sigma, graph


def wilson(k: int, r: int, z: float = 1.95996398454):
    phat = k / r; den = 1 + z*z/r
    center = (phat + z*z/(2*r)) / den
    half = z*np.sqrt(phat*(1-phat)/r + z*z/(4*r*r))/den
    return center-half, center+half


def summarize(raw: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    rows = []
    for idx, g in raw.groupby(keys, dropna=False, sort=True):
        idx = (idx,) if not isinstance(idx, tuple) else idx
        r = len(g); k = int(g['reject'].sum()); ph = k/r; lo, hi = wilson(k, r)
        row = dict(zip(keys, idx)); row.update(reps=r, rejections=k, rate=ph,
                mcse=np.sqrt(ph*(1-ph)/r), wilson_low=lo, wilson_high=hi)
        if 'seconds' in g:
            row['mean_seconds'] = g.seconds.mean()
        rows.append(row)
    return pd.DataFrame(rows)


def matched(cfg: dict, out: Path, rep_start: int, rep_stop: int | None):
    n,p,t,m = (cfg[k] for k in ('n','p','t','M'))
    reps=cfg['reps']; stop=reps if rep_stop is None else min(reps, rep_stop)
    target=np.arange(t); raw=[]; infos=[]
    cells=[('null',0.)] + [(pat,float(a)) for pat in ('dense','sparse') for a in cfg['signals'] if a>0]
    for kind in cfg['graphs']:
        om, sig, graph=graph_design(p,t,kind,cfg['bandwidth'],cfg['offdiag'])
        chol=np.linalg.cholesky(sig); s,b,v=conditional_parameters(om,target)
        gamma=np.zeros(len(s)); gamma[:min(4,len(s))]=np.array([1.,-.7,.5,.3])[:min(4,len(s))]
        np.savez_compressed(out/f'design_{kind}.npz',omega=om,graph=graph,gamma=gamma)
        for pattern,theta in cells:
            direction=np.ones(t)/np.sqrt(t) if pattern!='sparse' else np.eye(t)[0]
            direction=direction/np.sqrt(direction@v@direction)
            for rep in range(rep_start,stop):
                key=(kind,pattern,theta,rep)
                rng=generator(cfg['seed'],'data',*key)
                x=rng.normal(size=(n,p))@chol.T
                eps=rng.normal(size=n)
                mu=x[:,s]@b.T
                y=x[:,s]@gamma + theta*((x[:,target]-mu)@direction) + eps
                statistic=RidgeStatistic.fit(y,x[:,s],n*cfg['nuisance_penalty'],cfg['target_penalty'])
                t0=statistic(x[:,target])
                start=time.perf_counter()
                xc=oracle_copies(x,target,om,m,generator(cfg['seed'],'oracle',*key))
                por=rank_pvalue(t0,statistic(xc))
                raw.append(dict(graph=kind,pattern=pattern,theta=theta,rep=rep,method='Oracle-group',pvalue=por,
                                reject=por<=cfg['alpha'],seconds=time.perf_counter()-start,local_updates=0,stat_evals=m+1))
                for sweeps in cfg['sweeps']:
                    start=time.perf_counter()
                    copies=g_copies(x,target,graph,m,sweeps,generator(cfg['seed'],'g',sweeps,*key))
                    pg=rank_pvalue(t0,statistic(copies[:,:,target]))
                    raw.append(dict(graph=kind,pattern=pattern,theta=theta,rep=rep,method=f'G-group-L{sweeps}',pvalue=pg,
                                    reject=pg<=cfg['alpha'],seconds=time.perf_counter()-start,
                                    local_updates=(m+1)*sweeps*t,stat_evals=m+1))
                    if theta==0:
                        # For the Gaussian reference path only, with TRUE nuisance
                        # mean and direction. Not an estimated diagnostic for real data.
                        contrasts=np.einsum('n,mnt,t->m',eps,x[None,:,target]-copies[:,:,target],direction)
                        estimate=float(.5*np.mean(contrasts**2))
                        infos.append(dict(graph=kind,rep=rep,L=sweeps,retained_info=estimate,
                                          ratio=estimate/n,raw_info=n))
                start=time.perf_counter()
                pc,poc=coordinate_pvalues(x,target,graph,statistic,om)
                elapsed=time.perf_counter()-start
                for name,pv in [('G-coordinate-analytic',pc),('Oracle-coordinate-analytic',poc)]:
                    raw.append(dict(graph=kind,pattern=pattern,theta=theta,rep=rep,method=name,pvalue=pv,
                                    reject=pv<=cfg['alpha'],seconds=elapsed/2,local_updates=0,stat_evals=t))
            print(f'matched {kind} {pattern} theta={theta}: {stop-rep_start} reps',flush=True)
    raw=pd.DataFrame(raw); raw.to_csv(out/'raw.csv',index=False)
    summarize(raw,['graph','pattern','theta','method']).to_csv(out/'summary.csv',index=False)
    # Paired differences are computed at the replicate level, including within-data
    # dependence. They are not obtained by summing marginal binomial variances.
    pairs=[]
    for (kind,pat,th),g in raw.groupby(['graph','pattern','theta']):
        wide=g.pivot(index='rep',columns='method',values='reject').astype(float)
        for a in wide.columns:
            if a=='Oracle-group': continue
            delta=wide[a]-wide['Oracle-group']
            se=float(delta.std(ddof=1)/np.sqrt(len(delta))) if len(delta)>1 else np.nan
            pairs.append(dict(graph=kind,pattern=pat,theta=th,method=a,reference='Oracle-group',
                              mean_difference=delta.mean(),paired_mcse=se,reps=len(delta)))
    pd.DataFrame(pairs).to_csv(out/'paired_differences.csv',index=False)
    info=pd.DataFrame(infos); info.to_csv(out/'information_raw.csv',index=False)
    if len(info):
        si=info.groupby(['graph','L']).ratio.agg(['mean','std','count']).reset_index()
        si['mcse']=si['std']/np.sqrt(si['count'])
        # Exact corollary only applies to the dispersed, independent target set.
        om, sig, gr=graph_design(p,t,'dispersed',cfg['bandwidth'],cfg['offdiag'])
        _,_,vv=conditional_parameters(om,target); bb=np.ones(t)/np.sqrt(t); bb/=np.sqrt(bb@vv@bb)
        degrees=gr[target].sum(axis=1)
        exact=float(np.sum(bb*bb*np.diag(vv)*(n-degrees-1))/n)
        si['independent_target_exact_ratio']=np.where(si.graph=='dispersed',exact,np.nan)
        si.to_csv(out/'information_summary.csv',index=False)


def local_dimension(cfg: dict, out: Path, rep_start: int, rep_stop: int | None):
    n,p,t,m=(cfg[k] for k in ('n','p','t','M')); stop=cfg['reps'] if rep_stop is None else min(cfg['reps'],rep_stop)
    raw=[]; info=[]; target=np.arange(t)
    # Every target depends only on the same two nuisance covariates. Target
    # coordinates are independent conditional on XS. Extra neighbors are fixed.
    for q in cfg['degrees']:
        graph=np.zeros((p,p),bool)
        for j in target:
            graph[j,t:t+q]=True; graph[t:t+q,j]=True
        for theta in cfg['signals']:
            for rep in range(rep_start,stop):
                rng=generator(cfg['seed'],'dimension-data',theta,rep)
                xs=rng.normal(size=(n,p-t)); et=rng.normal(size=(n,t)); eps=rng.normal(size=n)
                xt=.7*xs[:,[0]]-.4*xs[:,[1]]+et
                x=np.column_stack((xt,xs)); direction=np.ones(t)/np.sqrt(t)
                y=xs[:,0]-.7*xs[:,1]+theta*(et@direction)+eps
                stat=RidgeStatistic.fit(y,xs,n*cfg['nuisance_penalty'],1.)
                copies=g_copies(x,target,graph,m,1,generator(cfg['seed'],'dimension-copies',q,theta,rep))
                pv=rank_pvalue(stat(xt),stat(copies[:,:,target]))
                raw.append(dict(degree=q,df=max(n-q-1,0),theta=theta,rep=rep,reject=pv<=cfg['alpha'],pvalue=pv))
                if theta==0:
                    ds=np.einsum('n,mnt,t->m',eps,xt[None]-copies[:,:,target],direction)
                    info.append(dict(degree=q,df=max(n-q-1,0),rep=rep,ratio=.5*np.mean(ds**2)/n,
                                     exact_ratio=max(n-q-1,0)/n))
    raw=pd.DataFrame(raw);raw.to_csv(out/'raw.csv',index=False)
    summarize(raw,['degree','df','theta']).to_csv(out/'summary.csv',index=False)
    info=pd.DataFrame(info);info.to_csv(out/'information_raw.csv',index=False)
    z=info.groupby(['degree','df','exact_ratio']).ratio.agg(['mean','std','count']).reset_index()
    z['mcse']=z['std']/np.sqrt(z['count']);z.to_csv(out/'information_summary.csv',index=False)


def omission(cfg: dict, out: Path, rep_start: int, rep_stop: int | None):
    raw=[];stop=cfg['reps'] if rep_stop is None else min(cfg['reps'],rep_stop)
    for n in cfg['sample_sizes']:
        for rho in cfg['correlations']:
            for rep in range(rep_start,stop):
                rng=generator(cfg['seed'],'omission',n,rho,rep)
                xs=rng.normal(size=n); x=rho*xs+np.sqrt(1-rho*rho)*rng.normal(size=n)
                y=xs+rng.normal(size=n)
                d=np.column_stack((np.ones(n),xs))
                yres=y-d@np.linalg.lstsq(d,y,rcond=None)[0]
                for name,dd,yy in [('Correct-graph',d,y),('Omitted-edge-raw',np.ones((n,1)),y),
                                   ('Omitted-edge-distilled',np.ones((n,1)),yres)]:
                    pv=spherical_pvalue(x,dd,yy)
                    raw.append(dict(n=n,rho=rho,rep=rep,method=name,pvalue=pv,reject=pv<=cfg['alpha']))
    raw=pd.DataFrame(raw);raw.to_csv(out/'raw.csv',index=False)
    summ=summarize(raw,['n','rho','method'])
    summ['tv_upper_bound']=np.minimum(1.,cfg['alpha']+.5*np.sqrt(-summ.n*np.log1p(-summ.rho**2)))
    summ.to_csv(out/'summary.csv',index=False)


def algebra(cfg: dict, out: Path, rep_start: int=0, rep_stop: int | None=None):
    rng=generator(cfg['seed'],'algebra'); n,s,t=cfg['n'],cfg['s'],cfg['t']; m=cfg['draws']
    xs=rng.normal(size=(n,s));d=np.column_stack((np.ones(n),xs))
    xt=xs@rng.normal(size=(s,t))+rng.normal(size=(n,t))
    # Arbitrary nonlinear null response, independent target errors.
    y=xs[:,0]**2+np.sin(xs[:,1])+rng.normal(size=n)
    pa,t0,df=dense_pvalue(y,xt,d)
    pf=f.sf((df-t)*t0/(t*(1-t0)),t,df-t)
    b=residual_basis(d);r=b.T@y
    copies=mvn_copies(xt,d,m,rng);z=b.T@copies
    q,_=np.linalg.qr(z,mode='reduced')
    tc=np.sum((q.swapaxes(-1,-2)@r)**2,axis=1)/(r@r)
    w,mean,var=ridge_spectrum_moments(b.T@xt,1.)
    u=rng.normal(size=(m,df));u/=np.linalg.norm(u,axis=1)[:,None]
    tr=np.sum(u*u*w,axis=1)
    summary=dict(n=n,s=s,t=t,d=df,draws=m,analytic_beta_p=pa,analytic_F_p=float(pf),
                 absolute_beta_F_error=abs(pa-pf),copy_beta_empirical_mean=float(tc.mean()),
                 copy_beta_theory_mean=t/df,copy_beta_empirical_variance=float(tc.var(ddof=1)),
                 copy_beta_theory_variance=2*t*(df-t)/(df*df*(df+2)),
                 ridge_mean_empirical=float(tr.mean()),ridge_mean_theory=mean,
                 ridge_variance_empirical=float(tr.var(ddof=1)),ridge_variance_theory=var)
    # Exact saturation: ZZ' = c I_d gives a constant ridge statistic even with
    # nontrivial Haar rotations. This is a constructed numerical identity, not a
    # claim that an iid Gaussian design is exactly spherical.
    zs=2*np.column_stack((np.eye(8),np.zeros((8,4))))
    ws,ms,vs=ridge_spectrum_moments(zs,1.)
    summary.update(saturated_ridge_mean=ms,saturated_ridge_variance=vs)
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    pd.DataFrame({'dense_copy_statistic':tc,'ridge_weighted_sphere_statistic':tr}).to_csv(out/'draws.csv',index=False)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--out',required=True)
    ap.add_argument('--rep-start',type=int,default=0);ap.add_argument('--rep-stop',type=int)
    args=ap.parse_args();cfg=json.loads(Path(args.config).read_text());out=Path(args.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f'Refusing to overwrite nonempty output directory: {out}')
    out.mkdir(parents=True,exist_ok=True)
    (out/'config.json').write_text(json.dumps(cfg,indent=2));started=time.time()
    globals()[cfg['experiment']](cfg,out,args.rep_start,args.rep_stop)
    import scipy,pandas,sklearn
    meta=dict(elapsed_seconds=time.time()-started,python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
              pandas=pandas.__version__,sklearn=sklearn.__version__,platform=platform.platform(),
              rep_start=args.rep_start,rep_stop=args.rep_stop,config_sha256=hashlib.sha256(Path(args.config).read_bytes()).hexdigest(),
              scope='NEW independently implemented experiment; not a reproduction of original Tables 1--5')
    (out/'run_metadata.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2),flush=True)

if __name__=='__main__': main()
