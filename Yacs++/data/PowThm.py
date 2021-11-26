import os
from module_generator import Generator,Module,Service
from module_generator import CPPComponent, PYComponent

MACH=os.uname()[0]+'_'+os.uname()[4]
HOME=os.getenv("HOME")
VERSION5=os.getenv("VERSION5")
KERNEL=os.getenv("SALOME_KERNEL")
PREREQUISITES=os.getenv("SALOME_KERNEL")+"/../env_build.sh"
print "YACSGEN: install POW/THM components on "+MACH
print "         Salome kernel: "+KERNEL

context={"update":1,
	"makeflags":"-j2",
	"prerequisites":HOME+"/salome_6.6.0/env_build.sh",
	"kernel":KERNEL
        }

#--------------
# Composant POW
#--------------
POW=CPPComponent("Composant_POW",
	services=[
		Service("runPOW",
			defs="""
#include "POW.hxx"
			""",
			inport=[("Power", "double"),],
			instream=[("fuel_temperature","CALCIUM_real", "I"),
				  ("water_temperature","CALCIUM_real", "I"),
				  ("water_density","CALCIUM_real", "I"),
				  ("continueLoop","CALCIUM_integer", "I")],
			outstream=[("powerTab","CALCIUM_real", "I"),
				  ("irradiationTab","CALCIUM_real", "I")],
			body="""
POW pow_;
pow_.initialize(Power, component);
pow_.run();
			""",
		),
	],
        libs="-L"+VERSION5+"/Yacs++/lib/"+MACH+" -lYacs++" \
            " -L"+VERSION5+"/Skin++/lib/"+MACH+" -lSkin++" \
            " -L"+VERSION5+"/Ganlib/lib/"+MACH+" -lGanlib" \
            " -L"+VERSION5+"/Donjon/lib/"+MACH+" -lDonjon" \
            " -L"+VERSION5+"/Dragon/lib/"+MACH+" -lDragon" \
            " -L"+VERSION5+"/Trivac/lib/"+MACH+" -lTrivac" \
            " -L"+VERSION5+"/Utilib/lib/"+MACH+" -lUtilib" \
            " -L"+VERSION5+"/Ganlib/lib/"+MACH+" -lGanlib" \
            " -lgfortran",
        includes="-I"+VERSION5+"/Yacs++/src" \
                " -I"+VERSION5+"/Skin++/src" \
                " -I"+VERSION5+"/Ganlib/src"
)

#--------------
# Composant THM
#--------------
THM=CPPComponent("Composant_THM",
	services=[
		Service("runTHM",
			defs="""
#include "THM.hxx"
			""",
			instream=[("powerTab","CALCIUM_real", "I"),
				  ("irradiationTab","CALCIUM_real", "I")],
			outstream=[("fuel_temperature","CALCIUM_real", "I"),
				  ("water_temperature","CALCIUM_real", "I"),
				  ("water_density","CALCIUM_real", "I"),
				  ("continueLoop","CALCIUM_integer", "I")],
			body="""
THM thm_;
thm_.initialize(component);
thm_.run();
			""",
		),
	],
        libs="-L"+VERSION5+"/Yacs++/lib/"+MACH+" -lYacs++" \
            " -L"+VERSION5+"/Skin++/lib/"+MACH+" -lSkin++" \
            " -L"+VERSION5+"/Ganlib/lib/"+MACH+" -lGanlib" \
            " -L"+VERSION5+"/Donjon/lib/"+MACH+" -lDonjon" \
            " -L"+VERSION5+"/Dragon/lib/"+MACH+" -lDragon" \
            " -L"+VERSION5+"/Trivac/lib/"+MACH+" -lTrivac" \
            " -L"+VERSION5+"/Utilib/lib/"+MACH+" -lUtilib" \
            " -L"+VERSION5+"/Ganlib/lib/"+MACH+" -lGanlib" \
            " -lgfortran",
        includes="-I"+VERSION5+"/Yacs++/src" \
                " -I"+VERSION5+"/Skin++/src" \
                " -I"+VERSION5+"/Ganlib/src"
)

g=Generator(Module("couplage_module",components=[POW, THM],prefix="./DonjonTHM_install"),context)
g.generate()
g.bootstrap()
g.configure()
g.make()
g.install()
g.make_appli("PowThm_appli",
             restrict=["KERNEL","GUI","COMPONENT","YACS","CALCULATOR"])
