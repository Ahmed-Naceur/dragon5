#---------------------------------------------------------------------------
#
#  Makefile for executing the Dragon non-regression tests
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
OS = $(shell uname -s | cut -d"_" -f1)
ifneq (,$(filter $(OS),SunOS AIX))
  MAKE = gmake
endif
all :
	$(MAKE) -C src
clean :
	$(MAKE) clean -C src
tests :
	./rdragon iaea2d.x2m -quiet
	./rdragon g2s_prestation.x2m -quiet
	./rdragon salmacro.x2m -quiet
	./rdragon tmacro.x2m  -quiet
	./rdragon tdraglib.x2m -quiet
	./rdragon twimsE.x2m -quiet
	./rdragon twlup.x2m -quiet
	./rdragon tmatxs2.x2m -quiet
	./rdragon trowland_shem295_jeff3.1.x2m -quiet
	./rdragon sens.x2m -quiet
	./rdragon testDuo.x2m -quiet
	./rdragon testDuo_B1.x2m -quiet
	./rdragon uo2_evo_xsm.x2m -quiet
	./rdragon C2D20.x2m -quiet
	./rdragon CFC-CELL.x2m -quiet
	./rdragon RegtestCNG_mccg.x2m -quiet
	./rdragon RegtestLZC_mccg.x2m -quiet
