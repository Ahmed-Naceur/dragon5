# -*- coding: utf-8 -*-
"""
Created on 14/02/2022
@author: BENEDET-BEN01
Python version : 3.6.9
DRAGON version : v5bev2068

This script allows to test 2 C2M procedures
FloatList: create a list of 1 floats
StringList: create a list of 1 character
"""

import lifo, cle2000
from assertS import *

stack1=lifo.new()
stack1.pushEmpty("my_LCM", "LCM")
stack2=lifo.new()
stack2.pushEmpty("my_LCM", "LCM")

#Store a double
my_float_list = cle2000.new('FloatList',stack1,impx=1)
my_float_list.exec()
stack1.lib()
stack1.node(0).lib()
value = stack1.node(0)["Float_List"][0]
print("value=",value)

#Store a character :
my_string_list = cle2000.new('StringList',stack2,impx=1)
my_string_list.exec()
stack2.lib()
stack2.node(0).lib()

assertS(stack1.node(0),'Float_List',0,3.141592654)
print("test Exemple_jdd completed")
