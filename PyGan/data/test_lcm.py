#
# test_lcm: non regression testing for lcm class
#
import lcm
from assertS import *
import numpy as np
my_lcm=lcm.new('LCM_INP','nonfuel')
my_lcm._impx=3
my_lcm.lib()
my_lcm.keys()
sign=my_lcm['SIGNATURE']
print('object signature=', sign)
daughter=my_lcm['REFL']
daughter.lib()
o2=daughter.copy('new_branch_of_dictionary')
state=o2['STATE-VECTOR']
print('state vector=', state)
o3=daughter['MIXTURES']
ia=np.array([8, 7, 8, 4, 9, 1, 0, 4], dtype='i')
ra=np.array([8.0,6.0,5.0,2.0,1.0], dtype='f')
da=np.array([8.0,6.0,5.0,2.0,1.0], dtype='d')
o2['key1']='new comments for this record'
o2['key2']=ia
o2['key3']=ra
o2['key4']=da
print('key2=',o2['key2'])
print('np type of key2=',o2['key2'].dtype)
print('np type of key3=',o2['key3'].dtype)
print('np type of key4=',o2['key4'].dtype)
o4=o2.rep('key5')
o2.lib()
print('o2 object name=',o2._name)
o5=o3[0]['CALCULATIONS'][0]
o5.lib()
print('ISOTOPESUSED=',o5['ISOTOPESUSED'])
print('ISOTOPESTEMP=',o5['ISOTOPESTEMP'])
print('ISOTOPESMIX=',o5['ISOTOPESMIX'])
assertS(o5,'ISOTOPESTEMP',0,5.23150024E+02)
del o5
print("test test_lcm completed")
