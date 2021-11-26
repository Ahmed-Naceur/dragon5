#
# test_lifo: non regression testing for lifo class
#
import lifo
import lcm
from assertS import *
import numpy as np
x=lifo.new()
x._impx=1
x.push('tata')
x.push(int)
x.push(12345678)
x.push(2.5)
x.push(3.5)
x.push(float)
x.pushEmpty("my_new_LCM_object","LCM")
x.push('this is very looong text')
val=x.node(3)
print("val(3)=",val)
val=x.pop()
print("val(pop)=",val)
ilen=x.getMax()
print("stack length=",ilen)
my_lcm=lcm.new('LCM_INP','nonfuel')
x.push(my_lcm)
x.push('more text')
new_lcm=x.node('nonfuel')
print("stack content before popping the lcm node")
x.lib()
x.pop()
x.pop()
print("stack content after popping twice the lcm node")
x.lib()
new_lcm.lib()
microlib=new_lcm['REFL']['MIXTURES'][0]['CALCULATIONS'][0]
assertS(microlib,'ISOTOPESTEMP',0,5.23150024E+02)
print("test test_lifo completed")
