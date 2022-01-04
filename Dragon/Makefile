#---------------------------------------------------------------------------
#
#  Makefile for executing the Dragon non-regression tests
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
#
OS = $(shell uname -s | cut -d"_" -f1)
ifneq (,$(filter $(OS),SunOS AIX))
  MAKE = gmake
endif
ifeq ($(intel),1)
  fcompiler = intel
else
  fcompiler =
endif
all :
	$(MAKE) -C src
clean :
	$(MAKE) clean -C src
tests :
	./rdragon iaea2d.x2m $(fcompiler) -quiet
	./rdragon g2s_prestation.x2m $(fcompiler) -quiet
	./rdragon salmacro.x2m $(fcompiler) -quiet
	./rdragon tmacro.x2m  $(fcompiler) -quiet
	./rdragon tdraglib.x2m $(fcompiler) -quiet
	./rdragon twimsE.x2m $(fcompiler) -quiet
	./rdragon twlup.x2m $(fcompiler) -quiet
	./rdragon tmatxs2.x2m $(fcompiler) -quiet
	./rdragon trowland_shem295_jeff3.1.x2m $(fcompiler) -quiet
	./rdragon sens.x2m $(fcompiler) -quiet
	./rdragon testDuo.x2m $(fcompiler) -quiet
	./rdragon testDuo_B1.x2m $(fcompiler) -quiet
	./rdragon uo2_evo_xsm.x2m $(fcompiler) -quiet
	./rdragon C2D20.x2m $(fcompiler) -quiet
	./rdragon CFC-CELL.x2m $(fcompiler) -quiet
	./rdragon ErmBeavrsPwrRefl.x2m $(fcompiler) -quiet
	./rdragon RegtestCNG_mccg.x2m $(fcompiler) -quiet
	./rdragon RegtestLZC_mccg.x2m $(fcompiler) -quiet
