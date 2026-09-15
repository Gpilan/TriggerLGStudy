"""W01: exercise production CFD/MPV synthesis without ROOT, plus measured widths."""
import json
import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import spe_response as spe
import plot_cfd_per_trigger_locallinear_peak as cfd
import plot_event_optical_mpv_matrix_root as mpv


def widths(t, y):
    p = max(range(len(y)), key=y.__getitem__)
    def cross(level, falling=False):
        target = level * y[p]
        for i in (range(p + 1, len(y)) if falling else range(1, p + 1)):
            if min(y[i-1], y[i]) <= target <= max(y[i-1], y[i]):
                return t[i-1] + (target-y[i-1])/(y[i]-y[i-1])*(t[i]-t[i-1])
        raise AssertionError('missing crossing')
    return cross(.9)-cross(.1), cross(.5, True)-cross(.5)


def wave(rise=1.3, width=3., dt=.005, tail=1e-10, counts=None, shift=0.):
    counts = counts or [1]
    centers = [shift + i * .5 for i in range(len(counts))]
    lo = [t-.005 for t in centers]; hi = [t+.005 for t in centers]
    kw = dict(sample_step_ns=dt, response_rise_ns=rise, response_fwhm_ns=width,
              gamma_shape=1.915604733026, gamma_tau_ns=None,
              transit_time_ns=14., kernel_tail_level=tail)
    t,y = cfd._build_reconstructed_waveform(lo,hi,counts,**kw)
    pulse = mpv._build_reconstructed_waveform(mpv.EventSummary(0,sum(counts),lo,hi,counts),peak_frac=.3,**kw)
    assert t == pulse.time and y == pulse.voltage
    return t,y


class ResponseTests(unittest.TestCase):
    def test_width_grid_tail(self):
        rows=[]
        for rise in [.65,1.3,2.6]:
            for width in [1.,3.,6.]:
                try:
                    response=spe.resolve_response(rise,width)
                except ValueError:
                    rows.append(dict(rise=rise,fwhm=width,status='unsupported'))
                    continue
                measured=[]
                for dt in [.05,.01,.005]:
                    t,y=wave(rise,width,dt)
                    r,w=widths(t,y); measured.append((r,w))
                    tol=.02 if dt==.05 else .001
                    self.assertLessEqual(abs(r-rise),tol)
                    self.assertLessEqual(abs(w-width),tol)
                    rows.append(dict(rise=rise,fwhm=width,dt=dt,measured_rise=r,measured_fwhm=w,status='pass'))
                self.assertLessEqual(max(abs(a-b) for a,b in zip(measured[-1],measured[-2])),.001)
                t2,y2=wave(rise,width,.005,1e-12)
                self.assertLessEqual(max(abs(a-b) for a,b in zip(widths(t2,y2),measured[-1])),1e-9)
                a=cfd._toa_local_linear_from_xy(t,y,.3,2)[0]
                b=cfd._toa_local_linear_from_xy(t2,y2,.3,2)[0]
                self.assertLessEqual(abs(a-b),1e-9)
        self.assertEqual(sum(r['status']=='pass' for r in rows),21)
        self.assertEqual([(r['rise'],r['fwhm']) for r in rows if r['status']=='unsupported'],[(1.3,1.),(2.6,1.)])
        print('WIDTH_TABLE='+json.dumps(rows))

    def test_multi_photon_and_cfd(self):
        t,y=wave(counts=[2,3,1])
        ta,ya=wave(counts=[2,0,1]); tb,yb=wave(counts=[0,3,0])
        self.assertEqual(t,ta)
        # Sum on a common physical time grid (builders use first/last nonzero hit).
        by_time={round(x,8):v for x,v in zip(tb,yb)}
        self.assertLessEqual(max(abs(v-a-by_time.get(round(x,8),0.)) for x,v,a in zip(t,y,ya)),1e-12)
        original=cfd._toa_local_linear_from_xy(t,y,.3,2)[0]
        ts,ys=wave(counts=[2,3,1],shift=2.)
        shifted=cfd._toa_local_linear_from_xy(ts,ys,.3,2)[0]
        scaled=cfd._toa_local_linear_from_xy(t,[4*v for v in y],.3,2)[0]
        self.assertLessEqual(abs(shifted-original-2.),.001)
        self.assertLessEqual(abs(scaled-original),.001)
        # Delayed photons make falling-branch width affect the rising sum/CFD.
        values=[]
        for width in [2.,3.,6.]:
            tx,yx=wave(width=width,counts=[1,0,0,0,1])
            values.append(cfd._toa_local_linear_from_xy(tx,yx,.3,2)[0])
        self.assertGreater(max(values)-min(values),.01)
        print('MULTI_PHOTON_CFD='+json.dumps(values))

    def test_guards_and_gamma_limit(self):
        for kw in [dict(rise_ns=float('nan'),fwhm_ns=3.),dict(rise_ns=1.3,fwhm_ns=float('inf')),
                   dict(rise_ns=1.3,fwhm_ns=.1),dict(rise_ns=1.3,fwhm_ns=3.,tau_ns=1.),
                   dict(rise_ns=1.3,fwhm_ns=3.,tail_level=0.),dict(rise_ns=1.3,fwhm_ns=3.,shape=0.)]:
            with self.assertRaises(ValueError): spe.resolve_response(**kw)
        response=spe.resolve_response(1.3,3.)
        for dt in [0.,-1.,float('nan'),float('inf')]:
            with self.assertRaises(ValueError): spe.sample_kernel(response,dt)
        shape=response.shape; tau=response.rise_tau_ns
        natural=(spe._fall_root(shape,.5)-spe._gamma_root(shape,.5,0.,shape,True))*tau
        same=spe.resolve_response(1.3,natural,tau_ns=tau)
        self.assertAlmostEqual(same.fall_tau_ns,tau,places=12)
        for i in range(1000):
            t=i*.01
            self.assertAlmostEqual(same.value(t),spe._gamma_response(shape,t/tau),places=12)
        self.assertEqual(cfd._build_reconstructed_waveform([],[],[],sample_step_ns=.05,response_rise_ns=1.3,response_fwhm_ns=3.,gamma_shape=shape,gamma_tau_ns=None,transit_time_ns=14.),([],[]))

if __name__=='__main__': unittest.main()
