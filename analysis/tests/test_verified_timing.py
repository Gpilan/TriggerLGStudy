import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from verified_timing import event_metrics,summarize,covariance_report
from spe_response import convolve_impulses
class MetricsTests(unittest.TestCase):
 def test_empty_kept(self):
  ev=dict(event=0,detected=[0,0],primary=[60000,0,0,0,0,0,1],time_bins=[[],[]])
  s=summarize([event_metrics(ev)])
  self.assertEqual(s['events'],1);self.assertEqual(s['coincidence'],0);self.assertEqual(s['success_all'],0)
  self.assertIsNone(s['width68_ns']);self.assertFalse(s['sqrt2_conversion'])
 def test_count_mismatch_rejected(self):
  with self.assertRaises(AssertionError):event_metrics(dict(event=0,primary=[],detected=[1,0],time_bins=[[],[]]))
 def test_covariance(self):
  self.assertEqual(covariance_report([(1,1),(2,2),(3,3)])['var_delta_ns2'],0)
  self.assertEqual(covariance_report([(1,3),(2,2),(3,1)])['var_delta_ns2'],4)
 def test_linear_convolution(self):
  impulse=[0,2,0,1,0,0];kernel=[0,.1,1,.5]
  expected=[sum(impulse[j]*kernel[i-j] for j in range(len(impulse)) if 0<=i-j<len(kernel)) for i in range(len(impulse))]
  for a,b in zip(convolve_impulses(impulse,kernel),expected):self.assertAlmostEqual(a,b,places=14)
if __name__=='__main__':unittest.main()
