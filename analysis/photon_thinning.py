"""Uniform detected-photon thinning without replacement, preserving event clusters."""
import copy
import numpy as np


def thin_event(event, target, rng):
    if target<=0 or any(n<target for n in event['detected']):
        raise ValueError('target must be positive and supported by both triggers')
    result=copy.deepcopy(event)
    for trig,bins in enumerate(event['time_bins']):
        counts=np.array([b[1] for b in bins],dtype=int)
        if int(counts.sum())!=event['detected'][trig]:raise ValueError('photon count mismatch')
        indices=rng.choice(int(counts.sum()),size=target,replace=False)
        chosen=np.bincount(np.searchsorted(np.cumsum(counts),indices,side='right'),minlength=len(counts))
        assert int(chosen.sum())==target and np.all(chosen<=counts)
        result['time_bins'][trig]=[[b[0],int(c),b[2],b[3]] for b,c in zip(bins,chosen) if c]
    result['detected']=[target,target]
    return result


def cluster_width_difference(a,b,rng,replicates=2000):
    """Rows are independent original events; columns are thinning realizations."""
    a=np.asarray(a,float);b=np.asarray(b,float)
    if len(a)<20 or len(b)<20:return {'status':'insufficient original events'}
    x=a[rng.integers(0,len(a),(replicates,len(a))),:].reshape(replicates,-1)
    y=b[rng.integers(0,len(b),(replicates,len(b))),:].reshape(replicates,-1)
    delta=(np.diff(np.nanquantile(x,[.16,.84],axis=1),axis=0)[0]-np.diff(np.nanquantile(y,[.16,.84],axis=1),axis=0)[0])/2
    value=(np.diff(np.nanquantile(a,[.16,.84]))[0]-np.diff(np.nanquantile(b,[.16,.84]))[0])/2
    return {'LG_minus_noLG_width68_ns':float(value),'cluster_bootstrap95_ns':np.quantile(delta,[.025,.975]).tolist(),
            'original_events':[len(a),len(b)],'replicates_per_event':a.shape[1],
            'method':'resample original event rows; keep all thinning realizations together'}
