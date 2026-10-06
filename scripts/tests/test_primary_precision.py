import sys,struct,math,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from validate_frozen_run import primary_matches, audited_primary_matches
class PrecisionTests(unittest.TestCase):
 def test_double_audit(self):
  request=[60000.,0.,0.,.1,0.,0.,1.]
  actual=[60000.,2.345,-1.23,.1,0.,0.,1.]
  self.assertTrue(audited_primary_matches(actual,request))
  actual[3]=.2
  self.assertFalse(audited_primary_matches(actual,request))
 def test_float_storage(self):
  norm=math.sqrt(1+.1**2+.05**2);request=[60000.,1.23456789,-2.34567891,0.,.1/norm,.05/norm,1/norm]
  stored=[struct.unpack('f',struct.pack('f',v))[0] for v in request]
  self.assertTrue(primary_matches(stored,request))
  self.assertFalse(primary_matches([stored[0],stored[1]+.01,*stored[2:]],request))
  self.assertFalse(primary_matches([*stored[:4],0,0,1],request))
  self.assertFalse(primary_matches([float('nan'),*stored[1:]],request))
if __name__=='__main__':unittest.main()
