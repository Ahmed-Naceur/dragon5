#!/bin/env python3
""" generation of pylib variable for use in Makefiles """

import os
from distutils.sysconfig import get_config_var
pylib = os.path.basename(get_config_var("LIBDIR"))
print(pylib)
