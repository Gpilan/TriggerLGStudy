"""Fixed-threshold LED selection; CFD remains the timing estimator.

k multiplies the analytic peak-normalized SPE amplitude (1 a.u.), never an
individual event's peak. Times and coincidence half-window are in ns.
"""
import math
import numpy as np
from plot_cfd_per_trigger_locallinear_peak import _build_reconstructed_waveform, _toa_local_linear_from_xy


def first_crossing(time, amplitude, threshold):
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError('LED threshold must be finite and positive')
    x, y = np.asarray(time, float), np.asarray(amplitude, float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError('Waveform arrays must be equal-length vectors')
    if not len(x):
        return dict(time_ns=None, status='no_signal')
    if not np.isfinite(x).all() or not np.isfinite(y).all() or np.any(np.diff(x) <= 0):
        raise ValueError('Waveform must be finite and strictly increasing in time')
    if y[0] >= threshold:
        return dict(time_ns=None, status='already_above_at_start')
    hits = np.flatnonzero((y[:-1] < threshold) & (y[1:] >= threshold))
    if not len(hits):
        return dict(time_ns=None, status='below_threshold')
    i = int(hits[0])
    t = x[i] + (threshold-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i])
    return dict(time_ns=float(t), status='crossed')


def select_pair(waveforms, k, window_ns=None, offset_ns=0., spe_peak=1.):
    if len(waveforms) != 2:
        raise ValueError('Exactly two channels are required')
    if not math.isfinite(k) or k <= 0 or not math.isfinite(spe_peak) or spe_peak <= 0:
        raise ValueError('k and SPE reference peak must be positive and finite')
    if not math.isfinite(offset_ns) or (window_ns is not None and (not math.isfinite(window_ns) or window_ns < 0)):
        raise ValueError('Invalid LED window or offset')
    channels = [first_crossing(t,y,k*spe_peak) for t,y in waveforms]
    both = all(c['time_ns'] is not None for c in channels)
    delta = channels[0]['time_ns']-channels[1]['time_ns'] if both else None
    timed = both and window_ns is not None and abs(delta-offset_ns) <= window_ns
    return dict(k=k,threshold=k*spe_peak,spe_peak=spe_peak,channels=channels,
                both_above_threshold=both,delta_LED_ns=delta,window_half_width_ns=window_ns,
                offset_ns=offset_ns,timed_coincidence=timed if window_ns is not None else None,
                selected=timed if window_ns is not None else both)


def analyze_event(event, ks, windows=(None,), offset_ns=0., sample_ns=.01):
    """Build each channel once, compute CFD once, then scan LED masks.

Censored waveforms remain in the incident denominator with an explicit status;
no timing is invented from incomplete time bins.
"""
    waveforms=[];cfd=[];quality=[]
    if len(event['detected']) != 2 or len(event['time_bins']) != 2:
        raise ValueError('Exactly two channels are required')
    for n,bins in zip(event['detected'],event['time_bins']):
        if sum(b[1] for b in bins) != n:
            raise ValueError('Detected count mismatch')
        if any(b[3] >= 99999 or b[2] < 0 for b in bins):
            waveforms.append(([],[]));cfd.append(None);quality.append('out_of_time_range');continue
        t,y=_build_reconstructed_waveform([b[2] for b in bins],[b[3] for b in bins],[b[1] for b in bins],
            sample_step_ns=sample_ns,response_rise_ns=1.3,response_fwhm_ns=3.,
            gamma_shape=1.915604733026,gamma_tau_ns=None,transit_time_ns=14.)
        toa,q=_toa_local_linear_from_xy(t,y,.3,2)
        waveforms.append((t,y));cfd.append(toa);quality.append(q['method'])
    delta=cfd[0]-cfd[1] if all(t is not None for t in cfd) else None
    selections=[]
    for k in ks:
        for window in windows:
            r=select_pair(waveforms,k,window,offset_ns)
            r['waveform_complete']='out_of_time_range' not in quality
            r['analysis_selected']=r['selected'] and delta is not None and r['waveform_complete']
            if not r['waveform_complete']:
                r.update(selected=None,both_above_threshold=None,timed_coincidence=None)
            selections.append(r)
    return dict(event=event['event'],CFD_times_ns=cfd,delta_CFD_ns=delta,CFD_quality=quality,LED=selections)


def analyze_file(path, ks, windows=(None,), offset_ns=0.):
    import json,hashlib
    from pathlib import Path
    source=Path(path);raw=source.read_bytes();events=json.loads(raw)
    if len({e['event'] for e in events})!=len(events):raise ValueError('Duplicate event IDs within input')
    rows=[analyze_event(e,ks,windows,offset_ns) for e in events]
    for row,event in zip(rows,events):
        row['primary']=event['primary'];row['primary_entered_tiles']=event.get('primary_entered_tiles')
    summaries=[]
    for i,(k,w) in enumerate((k,w) for k in ks for w in windows):
        sel=[r['LED'][i] for r in rows]
        count=sum(s['selected'] is True for s in sel);valid=sum(s['analysis_selected'] for s in sel)
        summaries.append(dict(k=k,window_half_width_ns=w,incident=len(rows),LED_selected=count,
            waveform_unknown=sum(not s['waveform_complete'] for s in sel),CFD_valid_selected=valid,
            LED_efficiency=count/len(rows) if rows else None,
            CFD_success_given_LED=valid/count if count else None,
            LED_efficiency_upper_bound=(count+sum(not s['waveform_complete'] for s in sel))/len(rows) if rows else None,
            delta_CFD_ns=[r['delta_CFD_ns'] for r in rows if r['LED'][i]['analysis_selected']]))
    return dict(input_sha256=hashlib.sha256(raw).hexdigest(),input=str(source),
        settings=dict(k=ks,window_half_width_ns=list(windows),offset_ns=offset_ns,SPE_peak=1.,CFD_fraction=.3,sample_ns=.01),
        incident_events=len(rows),CFD_before_LED_window_ns=[r['delta_CFD_ns'] for r in rows if r['delta_CFD_ns'] is not None],
        summary=summaries,events=rows)


if __name__=='__main__':
    import argparse,json
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--k',nargs='+',type=float,required=True)
    parser.add_argument('--window-ns',nargs='+',default=['none'],help='Coincidence half-width(s); none means thresholds only')
    parser.add_argument('--offset-ns',type=float,default=0.)
    args=parser.parse_args();windows=[None if s.lower()=='none' else float(s) for s in args.window_ns]
    if args.output.exists():raise FileExistsError(args.output)
    result=analyze_file(args.input,args.k,windows,args.offset_ns)
    with args.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
