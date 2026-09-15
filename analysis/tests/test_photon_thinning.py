import unittest,sys,copy
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from photon_thinning import thin_event,cluster_width_difference
class ThinningTests(unittest.TestCase):
 def test_counts_and_distribution(self):
  bins=[[0,50,0.,.01],[1,150,1.,1.01],[2,300,2.,2.01]]
  ev=dict(event=0,detected=[500,500],time_bins=[bins,copy.deepcopy(bins)]);original=copy.deepcopy(ev)
  rng=np.random.default_rng(42);draws=[]
  for i in range(1000):
   r=thin_event(ev,100,rng);self.assertEqual(r['detected'],[100,100])
   self.assertEqual(sum(b[1] for b in r['time_bins'][0]),100)
   counts={b[0]:b[1] for b in r['time_bins'][0]};draws.append([counts.get(i,0) for i in range(3)])
  observed=np.mean(draws,axis=0);p=np.array([.1,.3,.6]);error=6*np.sqrt(100*p*(1-p)*(400/499)/1000)
  self.assertTrue(np.all(np.abs(observed-100*p)<error));self.assertEqual(ev,original)
  self.assertEqual(thin_event(ev,100,np.random.default_rng(9)),thin_event(ev,100,np.random.default_rng(9)))
  with self.assertRaises(ValueError):thin_event(ev,501,rng)
 def test_cluster_bootstrap(self):
  r=cluster_width_difference(np.zeros((20,10)),np.zeros((20,10)),np.random.default_rng(1))
  self.assertEqual(r['cluster_bootstrap95_ns'],[0.,0.]);self.assertEqual(r['original_events'],[20,20])
  self.assertEqual(cluster_width_difference(np.zeros((2,10)),np.zeros((20,10)),np.random.default_rng(1))['status'],'insufficient original events')
if __name__=='__main__':unittest.main()
