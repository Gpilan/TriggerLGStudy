import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from LED import first_crossing, select_pair, analyze_event
from verified_timing import event_metrics
from spe_response import resolve_response

class LEDTests(unittest.TestCase):
 def test_absolute_not_fraction_and_first_crossing(self):
  self.assertEqual(first_crossing([0,1,2,3,4],[0,2,0,8,0],1)['time_ns'],.5)
  self.assertEqual(first_crossing([0,1,2],[0,4,0],1)['time_ns'],.25)
  self.assertIsNone(first_crossing([0,1],[0,.5],1)['time_ns'])
 def test_window(self):
  w=[([0,1,2],[0,2,0]),([1,2,3],[0,2,0])]
  self.assertFalse(select_pair(w,1,.9)['selected'])
  self.assertTrue(select_pair(w,1,1)['selected'])
  self.assertTrue(select_pair(w,1,0,offset_ns=-1)['selected'])
  self.assertIsNone(select_pair(w,1)['timed_coincidence'])
 def test_missing_invalid(self):
  self.assertFalse(select_pair([([],[]),([0,1],[0,2])],1)['selected'])
  for k in [0,-1,float('nan')]:
   with self.assertRaises(ValueError):first_crossing([0,1],[0,1],k)
  with self.assertRaises(ValueError):first_crossing([1,0],[0,1],1)
 def test_CFD_independent_and_censor(self):
  e=dict(event=0,primary=[60000,0,0,0,0,0,1],detected=[5,4],time_bins=[[[0,5,1,1.01]],[[0,4,1.1,1.11]]])
  a=analyze_event(e,[.5,2,100],[None,1]);b=event_metrics(e)
  self.assertEqual(a['delta_CFD_ns'],b['delta_ns'])
  self.assertTrue(a['LED'][0]['analysis_selected']);self.assertFalse(a['LED'][-1]['analysis_selected'])
  e['time_bins'][0][0][3]=99999
  a=analyze_event(e,[1]);self.assertIsNone(a['LED'][0]['selected']);self.assertFalse(a['LED'][0]['analysis_selected'])
 def test_spe_peak(self):
  r=resolve_response(1.3,3.);self.assertAlmostEqual(r.value(r.peak_ns),1.,places=12)
if __name__=='__main__':unittest.main()
