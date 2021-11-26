import os
from module_generator import Generator,Module,Service
from module_generator import CPPComponent, PYComponent

MACH=os.uname()[0]+'_'+os.uname()[4]
HOME=os.getenv("HOME")
VERSION5=os.getenv("VERSION5")
KERNEL=os.getenv("SALOME_KERNEL")
PREREQUISITES=os.getenv("SALOME_KERNEL")+"/../env_build.sh"
print "YACSGEN: install FACT component on "+MACH
print "         Salome kernel: "+KERNEL

context={"update":1,
	"makeflags":"-j2",
	"prerequisites":PREREQUISITES,
	"kernel":KERNEL
        }

#---------------
# Composant FACT
#---------------
FACT=CPPComponent("Composant_FACT",
	compodefs="""
#include "FACT.hxx"
static FACT fact_;
	""",
	services=[
		Service("runFACT",
			inport=[("a", "long"),],
			outport=[("b", "long"),],
			body="""
long bb;
fact_.run(a, bb);
b=bb;
			""",
		),
	],
        libs="-L"+VERSION5+"/Yacs++/lib/"+MACH+" -lYacs++" \
            " -L"+VERSION5+"/Skin++/lib/"+MACH+" -lSkin++" \
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

g=Generator(Module("couplage_module_fact",components=[FACT],prefix="./FACT_install"),context)
g.generate()
g.bootstrap()
g.configure()
g.make()
g.install()
g.make_appli("FACT_appli",
             restrict=["KERNEL","GUI","COMPONENT","YACS","CALCULATOR"])
