#!/usr/bin/env python3
"""Recompute manuscript summaries from saved replication-level results."""
from pathlib import Path
import json
import re
import math
import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / os.environ.get('REPLICATION_REFERENCE','reference')
OUT = ROOT / 'outputs/rebuilt'
OUT.mkdir(parents=True, exist_ok=True)
CHECKS = []
TABLES = {}


def save(name, frame):
    frame.to_csv(OUT / (name + '.csv'), index=False, float_format='%.17g')


def reject_summary(frame, keys, p='pvalue', inclusive=False, expected=None, sample_se=True):
    result = []
    for key, group in frame.groupby(keys, sort=True, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        valid = np.isfinite(group[p])
        if 'status' in group:
            valid &= group.status.isin(['ok','numerically_certified'])
        values = group.loc[valid, p].to_numpy()
        assert np.all((values >= 0) & (values <= 1))
        indicators = values <= .05 if inclusive else values < .05
        n, count = len(values), int(indicators.sum())
        assert n > 0 and (expected is None or len(group) == expected), (key, len(group), expected)
        rate = count/n
        result.append(dict(zip(keys, key), attempts=len(group), valid=n,
                           failed=len(group)-n, rejections=count, estimate=rate,
                           se=math.sqrt(rate*(1-rate)/(n-1 if sample_se else n)),
                           ci_lower=0 if count == 0 else beta.ppf(.025,count,n-count+1),
                           ci_upper=1 if count == n else beta.ppf(.975,count+1,n-count)))
    return pd.DataFrame(result)


def vec(frame, where, field='estimate'):
    z=frame
    for k,v in where.items():
        z=z.loc[z[k].eq(v)]
    assert len(z)==1, (where,len(z))
    return float(z.iloc[0][field])


def paired(frame, keys, methods, value='reject', rep='replicate'):
    result=[]
    for key,g in frame.groupby(keys,sort=True):
        key=key if isinstance(key,tuple) else (key,)
        pivot=g.pivot(index=rep,columns='method',values=value)
        assert not pivot.isna().any().any()
        for a,b in methods:
            d=pivot[a].astype(float)-pivot[b].astype(float)
            result.append(dict(zip(keys,key),method_a=a,method_b=b,pairs=len(d),
                               difference=d.mean(),paired_mcse=d.std(ddof=1)/math.sqrt(len(d))))
    return pd.DataFrame(result)


def table(name, rows):
    TABLES[name]=rows
    path=OUT/(Path(name).stem+'.tex')
    lines=['% Generated numerical rows; column definitions are in RESULT_MAP.md.',
           '\\begin{tabular}{'+'r'*len(rows[0])+'}']
    lines.extend(' & '.join(f'{v:.10g}' for v in row)+' \\\\' for row in rows)
    path.write_text('\n'.join(lines+['\\end{tabular}'])+'\n')


def check_table(name, selector, skip=0):
    text=(ROOT/'expected'/name).read_text()
    lines=[line for line in text.splitlines() if selector(line)]
    computed=TABLES[name]
    assert len(lines)==len(computed), (name,len(lines),len(computed))
    count=0
    for line,values in zip(lines,computed):
        cells=line.split('&')[skip:]
        tokens=[]
        for cell in cells:
            # Data rows only. Formatting dimensions and row indices are excluded.
            cell=re.sub(r'\\(?:shortstack|hspace)(?:\[[^]]*\])?', '', cell)
            tokens.extend(re.findall(r'<?-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?',cell))
        assert len(tokens)==len(values), (name,tokens,values)
        for printed,value in zip(tokens,values):
            if printed.startswith('<'):
                assert value < float(printed[1:]), (name,printed,value)
            else:
                decimal=len(printed.split('.')[1]) if '.' in printed else 0
                tolerance=.5000001*10**(-decimal) if decimal else 1e-10
                assert abs(float(printed)-value)<=tolerance, (name,printed,value,tolerance)
            count+=1
    CHECKS.append(dict(check=name,items=count,status='passed',basis='Each displayed value versus replication-level recomputation; tolerance is half its printed unit.'))


def main():
    methods=['Oracle CRT','Plug-in MLE CRT','Plug-in CRT (d.f. adjustment)','MVN-CRT','LW plug-in CRT (D.10)']
    setting_names=['LD-tight','LD-moderate','p>n target']
    size_raw=pd.read_csv(REF/'gaussian/E01/results.csv')
    size=reject_summary(size_raw,['setting','statistic','method'],'p_value',True,5000,False)
    power_raw=pd.read_csv(REF/'gaussian/E02/results.csv')
    power=reject_summary(power_raw,['setting','theta','method'],'p_value',True,400,False)
    save('gaussian-size',size); save('gaussian-power',power)
    rows=[]
    for setting,dims in zip(setting_names,[(60,45,10,14),(80,50,10,29),(80,60,40,19)]):
        row=list(dims)
        for m in methods:
            where=dict(setting=setting,statistic='ridge gain',method=m)
            row += [100*vec(size,where),100*vec(size,where,'se')]
        rows.append(row)
    name='Simulation/CRT/oracle_plugin_benchmark_size.tex'
    table(name,rows);check_table(name,lambda l: bool(re.match(r'^\$\d',l)))
    gdiff=power.pivot(index=['setting','theta'],columns='method',values='estimate')
    power_gap=float((gdiff['Oracle CRT']-gdiff['MVN-CRT']).abs().max())
    assert abs(power_gap-.0225)<1e-12
    CHECKS.append(dict(check='Section 7.1 maximum oracle/MVN difference',value=power_gap,manuscript=.023,status='passed'))

    graphical=[]
    configs=json.loads((ROOT/'experiments/reproducible-study/config/simulations.json').read_text())
    configs=configs['scenarios'] if isinstance(configs,dict) else configs
    for path in sorted((REF/'graphical').glob('*.csv')):
        if path.name=='counts.csv': continue
        raw=pd.read_csv(path)
        if 'pvalue' not in raw: continue
        cfg=next(c for c in configs if c['id']==path.stem)
        z=reject_summary(raw,['scenario','theta','graph','method'],expected=cfg['repetitions'])
        assert not raw.duplicated(['scenario','replicate','theta','graph','method']).any()
        graphical.append(z)
    graphical=pd.concat(graphical,ignore_index=True)
    old=pd.read_csv(REF/'graphical/counts.csv')
    old=old.loc[old.scenario.isin(graphical.scenario.unique()) & old.metric.eq('rejection')]
    joined=graphical.merge(old,on=['scenario','theta','graph','method'],suffixes=('','_reference'),validate='one_to_one')
    assert len(joined)==len(graphical)
    for col in ['estimate','se','attempts','valid','failed','rejections']:
        assert np.allclose(joined[col],joined[col+'_reference'],atol=1e-12,rtol=0),col
    assert graphical.failed.sum()==2
    save('graphical',graphical)
    CHECKS.append(dict(check='Figures 3-9 plotted estimates and SEs',items=len(graphical),failed_method_rows_preserved=2,status='passed'))

    e5=[json.loads(s) for s in (REF/'maxcoeft/raw.jsonl').read_text().splitlines()]
    assert len(e5)==7000
    matched=REF/'matched/raw.jsonl'
    if not matched.exists(): matched=ROOT/'experiments/d5-joint-comparison-20261007/results/complete/raw.jsonl'
    old=[json.loads(s) for s in matched.read_text().splitlines()]
    identity=lambda r: tuple(r[k] for k in ['design','arrangement','pattern','theta','rep'])
    old={identity(r):r for r in old}
    rows=[]
    for r in e5:
        baseline=old[identity(r)]
        assert r['reused_pvalues']==baseline['pvalues']
        assert r['data_sha256']==baseline['data_sha256']
        pv=baseline['pvalues']|r['pvalues']
        for sampler in ['G','Oracle']:
            assert pv[sampler+'RidgeMaxCoefTBonf']==min(1,2*min(pv[sampler+'Ridge'],pv[sampler+'MaxCoefT']))
        for m in ['CoordG','GRidge','GMaxCoefT','GRidgeMaxCoefTBonf','OracleRidge','OracleMaxCoefT','OracleRidgeMaxCoefTBonf']:
            rows.append(dict(design=r['design'],pattern=r['pattern'],theta=r['theta'],replicate=r['rep'],method=m,pvalue=pv[m]))
    e5_raw=pd.DataFrame(rows)
    e5_summary=reject_summary(e5_raw,['design','pattern','theta','method'],inclusive=True,expected=500)
    save('matched',e5_summary)
    e5_raw['reject']=e5_raw.pvalue<=.05
    e5_pairs=paired(e5_raw,['design','pattern','theta'],[('GRidgeMaxCoefTBonf','GRidge'),('GRidgeMaxCoefTBonf','GMaxCoefT')])
    save('matched-paired',e5_pairs)
    for pattern,method,difference,se in [('dense','GRidge',-.058,.011),('sparse','GMaxCoefT',.058,.014)]:
        where=dict(design='B',pattern=pattern,theta=.75,method_b=method)
        assert abs(vec(e5_pairs,where,'difference')-difference)<1e-12
        assert abs(vec(e5_pairs,where,'paired_mcse')-se)<.0005
    CHECKS.append(dict(check='E.5 paired Bonferroni comparison changes and MCSEs',status='passed'))
    rows=[]
    for design in ['A','B']:
        for pattern,theta in [('null',0)]+[(p,t) for p in ['dense','sparse'] for t in [.5,.75,1]]:
            rows.append([theta]+[100*vec(e5_summary,dict(design=design,pattern=pattern,theta=theta,method=m)) for m in
                                  ['CoordG','GRidge','GMaxCoefT','GRidgeMaxCoefTBonf','OracleRidge','OracleMaxCoefT','OracleRidgeMaxCoefTBonf']])
    name='Simulation/CRT/joint_coordinatewise_oracle_comparison.tex'
    table(name,rows);check_table(name,lambda l: bool(re.match('^(Null|Dense|Sparse) &',l)),1)
    gaps=[]
    for (d,p,t),grp in e5_summary.groupby(['design','pattern','theta']):
        for m in ['Ridge','MaxCoefT','RidgeMaxCoefTBonf']:
            gaps.append(abs(vec(grp,{'method':'G'+m})-vec(grp,{'method':'Oracle'+m})))
    assert abs(max(gaps)-.044)<1e-12
    CHECKS.append(dict(check='E.5 maximum G-CRT/oracle difference',value=max(gaps),manuscript=.044,status='passed'))

    non_raw=pd.concat([pd.read_csv(f) for f in sorted((REF/'nongaussian').glob('cell-*.csv'))])
    assert len(non_raw)==12000 and non_raw.status.eq('ok').all()
    levels=np.arange(1,21)/200
    rows=[]
    for (dist,n),g in non_raw.groupby(['distribution','n']):
        assert len(g)==1000 and not g.repetition.duplicated().any()
        for alpha in levels:
            k=int((g.p_value<=alpha).sum());estimate=k/1000
            z=1.959963984540054;den=1+z*z/1000
            mid=(estimate+z*z/2000)/den
            half=z*math.sqrt(estimate*(1-estimate)/1000+z*z/(4*1000**2))/den
            rows.append(dict(distribution=dist,n=n,alpha=alpha,rejections=k,estimate=estimate,lower=mid-half,upper=mid+half))
    non=pd.DataFrame(rows);save('nongaussian',non)
    assert len(non)==240
    expected_non=pd.read_csv(REF/'nongaussian/rejection-rate-intervals.csv')
    check=non.merge(expected_non,on=['distribution','n','alpha'],validate='one_to_one',suffixes=('','_old'))
    assert len(check)==240
    for a,b in [('estimate','rejection_rate'),('lower','wilson_lower'),('upper','wilson_upper')]:
        assert np.allclose(check[a],check[b],atol=1e-14,rtol=0)
    rows=[]
    for n in [50,100,200,400]:
        row=[n]
        for d in ['Gaussian','t5','t3']:
            row += [100*vec(non,dict(n=n,distribution=d,alpha=.05),f) for f in ['estimate','lower','upper']]
        rows.append(row)
    name='inline/nongaussian-size.tex'
    table(name,rows);check_table(name,lambda l:bool(re.match(r'^\$\d',l)))
    CHECKS.append(dict(check='Figure 10 and Table 4: all Wilson intervals recomputed',items=240,status='passed'))

    direction=pd.read_csv(REF/'direction/direction-tests.csv')
    direction['case']=direction['case'].fillna('null')
    dsum=reject_summary(direction,['n','case'],inclusive=True,expected=500)
    rows=[]
    for n in [40,120,360]:
        for case in ['null','positive','negative']:
            g=direction.loc[(direction.n==n)&(direction['case']==case)]
            w=dict(n=n,case=case)
            rows.append([n,g.T0.mean(),g.Tbar.mean(),g.gap.mean(),g.gap.std()/math.sqrt(500),
                         100*vec(dsum,w),100*vec(dsum,w,'ci_lower'),100*vec(dsum,w,'ci_upper')])
    name='Simulation/CRT/direction_mechanism.tex'
    table(name,rows)
    # The Response column contains text only and contributes no number.
    check_table(name,lambda l: bool(re.match(r'^\$\d',l)))
    save('direction',direction.groupby(['n','case']).agg(mean_T0=('T0','mean'),mean_Tbar=('Tbar','mean'),mean_gap=('gap','mean')).reset_index().merge(dsum,on=['n','case']))

    methods=['original','centered_pool','centered_identity','scalar_twosided']
    center=pd.read_csv(REF/'centering/tests.csv')
    center['case']=center['case'].fillna('null')
    cs=reject_summary(center,['n','case','method'],inclusive=True,expected=100)
    cp=paired(center,['n','case'],[(m,'original') for m in methods[1:]])
    save('centering',cs);save('centering-paired',cp)
    rows=[];diff=[]
    for case in ['null','positive','negative','orthogonal']:
        for n in [40,120,360]:
            rows.append([n]+[100*vec(cs,dict(n=n,case=case,method=m)) for m in methods])
            row=[n]
            for m in methods[1:]:
                row.extend([100*vec(cp,dict(n=n,case=case,method_a=m),'difference'),100*vec(cp,dict(n=n,case=case,method_a=m),'paired_mcse')])
            diff.append(row)
    for fname,values in [('centering_rates',rows),('centering_paired',diff)]:
        name=f'Simulation/CRT/{fname}.tex';table(name,values);check_table(name,lambda l: bool(re.match(r'^\$\d',l)))

    connected=pd.read_csv(REF/'connected/tests.csv')
    assert len(connected)==100000
    qs=reject_summary(connected,['b','r','method'],inclusive=True,expected=2500,sample_se=False)
    methods=['original','uncentered_pool','centered_pool','uncentered_identity','centered_identity']
    qp=paired(connected,['b','r'],[('centered_pool','uncentered_pool'),('centered_identity','uncentered_identity')])
    sign=[]
    for r,g in connected.groupby('r'):
        for m in methods:
            p=g.loc[g.method==m].pivot(index='replicate',columns='b',values='reject')
            d=p[1].astype(float)-p[-1].astype(float)
            sign.append(dict(r=r,method=m,difference=d.mean(),paired_mcse=d.std()/50))
    sign=pd.DataFrame(sign)
    save('connected',qs);save('connected-paired',qp);save('connected-sign',sign)
    rates=[];differences=[]
    for b in [-1,1]:
        for r in [0,.15,.3,.6]:
            rates.append([r]+[100*vec(qs,dict(b=b,r=r,method=m)) for m in methods])
            row=[r]
            for m in ['centered_pool','centered_identity']:
                row.extend([100*vec(qp,dict(b=b,r=r,method_a=m),'difference'),100*vec(qp,dict(b=b,r=r,method_a=m),'paired_mcse')])
            differences.append(row)
    for r in [0,.15,.3,.6]:
        row=[r]
        for m in methods:
            row.extend([100*vec(sign,dict(r=r,method=m),'difference'),100*vec(sign,dict(r=r,method=m),'paired_mcse')])
        differences.append(row)
    for fname,values in [('connected_rates',rates),('connected_paired',differences)]:
        name=f'Simulation/CRT/{fname}.tex';table(name,values);check_table(name,lambda l: bool(re.match(r'^\$\d',l)))

    stock=pd.concat([pd.read_csv(p) for p in sorted((REF/'stock').glob('stock_*.csv'))])
    assert len(stock)==14012 and stock.status.eq('ok').all()
    null=stock.loc[stock.scenario.str.startswith('stock_typeI')]
    ns=reject_summary(null,['scenario','method'],expected=400)
    noise=stock.loc[stock.scenario.str.startswith('stock_noise')]
    nsum=reject_summary(noise,['theta','method'],expected=400)
    obs=stock.loc[stock.scenario.str.startswith('stock_three')]
    save('stock-null',ns);save('stock-noise',nsum);save('stock-observed',obs)
    rows=[]
    for group in ['T1','T2']:
        row=[]
        for m in ['LM-L1-R-SSR','RF-RR','De-Lasso','dCRT (L1)','dCRT (RF)']:
            row.extend([100*vec(ns,dict(scenario='stock_typeI_'+group,method=m)),100*vec(ns,dict(scenario='stock_typeI_'+group,method=m),'se')])
        rows.append(row)
    name='Application/Stock/Stock_Type-I-error.tex';table(name,rows)
    check_table(name,lambda l:l.startswith('$\\mc T_'),1)
    rows=[]
    for group in ['T1','T2']:
        rows.append([vec(obs,dict(scenario='stock_three_stat_replay_'+group,method=m),'pvalue') for m in ['LM-L1-R-SSR','RF-RR','dCRT (L1)','dCRT (RF)']])
    name='Application/Stock/Stock_pvalue.tex';table(name,rows);check_table(name,lambda l:l.startswith('$\\mc{T}_'),1)
    for method,expected in [('LM-L1-R-SSR',.9425),('dCRT (L1)',.7275)]:
        assert abs(vec(nsum,dict(theta=.1,method=method))-expected)<1e-12
    CHECKS.append(dict(check='Section 8.1 94.3% versus 72.8% and all stock null/observed numbers',status='passed'))
    gof=pd.read_csv(REF/'stock/Stock_GoF.csv')
    name='Application/Stock/Stock_GoF.tex'
    table(name,[[vec(gof,dict(method=m),'pvalue') for m in ['VV','DP','F_sum']]])
    check_table(name,lambda l:'p-value' in l,1)
    bcr=pd.read_csv(REF/'bcr/UNRESTRICTED_RESULTS.csv')
    bs=reject_summary(bcr,['method'],expected=400)
    save('bcr-null',bs)
    row=[]
    for m in ['LRT','RF','dCRT (L1)','dCRT (RF)']:
        row.extend([100*vec(bs,dict(method=m)),100*vec(bs,dict(method=m),'se')])
    name='Application/BCR/BCR_Type-I-error.tex';table(name,[row]);check_table(name,lambda l:l.startswith('\\shortstack[l]'),1)
    bo=pd.read_csv(REF/'bcr/observed.csv');save('bcr-observed',bo)
    cert=pd.read_csv(REF/'bcr/likelihood-certificates.csv')
    assert len(cert)==401 and cert.replicate.nunique()==401
    assert cert[['full_gap','reduced_gap']].to_numpy().max()<=1e-5
    for col in ['full_exact_face','reduced_exact_face','full_exact_dual','reduced_exact_dual']:
        assert cert[col].all()
    for _,r in cert.iterrows():
        actual=vec(bo,dict(method='LRT'),'pvalue') if r.replicate==0 else vec(bcr,dict(method='LRT',replicate=r.replicate),'pvalue')
        assert abs(actual-r.pvalue)<1e-10
    CHECKS.append(dict(check='BCR 401 saved objective-bound summaries, p-values, and certificate flags',status='passed',maximum_gap=float(cert[['full_gap','reduced_gap']].to_numpy().max())))
    expected={'LRT':9.58e-6,'RF':.018,'dCRT (L1)':.697,'dCRT (RF)':.407}
    for m,v in expected.items():
        actual=vec(bo,dict(method=m),'pvalue')
        assert abs(actual-v)<(5e-9 if m=='LRT' else .000501),(m,actual,v)
    CHECKS.append(dict(check='Section 8.2 all four observed BCR p-values and 160/29/13/13 null counts',status='passed'))
    (OUT/'table-values.json').write_text(json.dumps(TABLES,indent=2)+'\n')
    (OUT/'verification.json').write_text(json.dumps(dict(status='passed',checks=CHECKS,scope='Saved replication summaries; no model fitting in this command.'),indent=2)+'\n')
    print(f'PASS: {len(CHECKS)} verification groups; {len(TABLES)} numerical tables rebuilt.')


if __name__=='__main__':
    main()
