"""Event-level timing metrics on validated canonical ROOT exports (ns throughout)."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from plot_cfd_per_trigger_locallinear_peak import _build_reconstructed_waveform, _toa_local_linear_from_xy


def covariance_report(pairs):
    if len(pairs)<2:return {'n':len(pairs),'status':'insufficient'}
    a=np.asarray(pairs,float);v=np.cov(a.T,ddof=1);observed=float(np.var(a[:,0]-a[:,1],ddof=1))
    predicted=float(v[0,0]+v[1,1]-2*v[0,1])
    if not np.isclose(observed,predicted,rtol=1e-10,atol=1e-12):raise ValueError('variance identity failed')
    return dict(n=len(a),var_t1_ns2=float(v[0,0]),var_t2_ns2=float(v[1,1]),cov_ns2=float(v[0,1]),
                var_delta_ns2=observed,identity_residual_ns2=observed-predicted,
                single_detector_resolution=None,reason='pair timing only; no sqrt(2) assumption')


def width(values):
    return float(np.diff(np.quantile(values,[.16,.84]))[0]/2) if len(values)>=2 else None


def event_metrics(event,dt=.01,fwhm=3.,fraction=.3,half=2):
    out={'event':event['event'],'N':event['detected'],'primary':event['primary'],'triggers':[]}
    for n,bins in zip(event['detected'],event['time_bins']):
        assert sum(b[1] for b in bins)==n
        tr={'N':n,'toa_ns':None,'quality':{'method':'no_data'},'arrival':None}
        if n:
            centers=np.array([(b[2]+b[3])/2 for b in bins]);counts=np.array([b[1] for b in bins])
            # Counts are integer photoelectrons; exact empirical event quantiles on bin centers.
            photons=np.repeat(centers,counts);q10,q50,q90=np.quantile(photons,[.1,.5,.9])
            tr['arrival']=dict(q10_ns=float(q10),q50_ns=float(q50),q90_minus_q10_ns=float(q90-q10),late_fraction=float(np.mean(photons>q10+10)))
            if any(b[3]>=99999 or b[2]<0 for b in bins):
                tr['quality']={'method':'out_of_time_range'}
            else:
                t,y=_build_reconstructed_waveform([b[2] for b in bins],[b[3] for b in bins],counts,
                 sample_step_ns=dt,response_rise_ns=1.3,response_fwhm_ns=fwhm,
                 gamma_shape=1.915604733026,gamma_tau_ns=None,transit_time_ns=14.)
                tr['toa_ns'],tr['quality']=_toa_local_linear_from_xy(t,y,fraction,half)
        out['triggers'].append(tr)
    times=[t['toa_ns'] for t in out['triggers']]
    out['delta_ns']=times[0]-times[1] if all(t is not None for t in times) else None
    return out


def summarize(events):
    total=len(events);valid=[e for e in events if e['delta_ns'] is not None]
    delta=[e['delta_ns'] for e in valid];coinc=sum(all(n>0 for n in e['N']) for e in events)
    w=width(delta);med=float(np.median(delta)) if delta else None
    failures={};low_r2=0
    for e in events:
        for t in e['triggers']:
            method=t['quality']['method'];failures[method]=failures.get(method,0)+1
            low_r2+= t['quality'].get('r2') is not None and t['quality']['r2']<.90
    result=dict(events=total,coincidence=coinc,coincidence_fraction=coinc/total if total else None,
      timing_success=len(valid),success_all=len(valid)/total if total else None,
      success_given_coincidence=len(valid)/coinc if coinc else None,
      N_mean=np.mean([e['N'] for e in events],axis=0).tolist() if total else [],
      N_zero_fraction=[sum(e['N'][i]==0 for e in events)/total for i in [0,1]] if total else [],
      width68_ns=w,std_delta_ns=float(np.std(delta,ddof=1)) if len(delta)>1 else None,
      median_delta_ns=med,tail_fraction=float(np.mean(np.abs(np.array(delta)-med)>3*w)) if w and med is not None else None,
      covariance=covariance_report([[e['triggers'][0]['toa_ns'],e['triggers'][1]['toa_ns']] for e in valid]),
      extraction_methods=failures,low_r2_trigger_count=int(low_r2),sqrt2_conversion=False)
    for i in [0,1]:
        arrivals=[e['triggers'][i]['arrival'] for e in events if e['triggers'][i]['arrival'] is not None]
        result[f'arrival_T{i+1}']={'events':len(arrivals),**{key:float(np.mean([a[key] for a in arrivals])) for key in ['q10_ns','q50_ns','q90_minus_q10_ns','late_fraction']} } if arrivals else {'events':0}
        slopes=[e['triggers'][i]['quality']['slope_per_ns'] for e in events if 'slope_per_ns' in e['triggers'][i]['quality']]
        result[f'mean_slope_T{i+1}']=float(np.mean(slopes)) if slopes else None
    return result


def analyze(paths,**kwargs):
    events=[];inputs={}
    for path in paths:
        data=Path(path).read_bytes();inputs[str(path)]=hashlib.sha256(data).hexdigest()
        for e in json.loads(data):
            row=event_metrics(e,**kwargs);row['run']=str(Path(path).parent);events.append(row)
    return dict(inputs_sha256=inputs,settings=kwargs,summary=summarize(events),events=events)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('inputs',type=Path,nargs='+');p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=analyze(a.inputs)
    with a.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
