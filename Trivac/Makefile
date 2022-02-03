#---------------------------------------------------------------------------
#
#  Makefile for executing the Trivac non-regression tests
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
	./rtrivac DIFtst.x2m $(fcompiler) -quiet
	./rtrivac Ktests.x2m $(fcompiler) -quiet
	./rtrivac SPNtst.x2m $(fcompiler) -quiet
	./rtrivac iaea3d.x2m $(fcompiler) -quiet
	./rtrivac monju3D_spn.x2m $(fcompiler) -quiet
	./rtrivac multigroup_albedo_2d.x2m $(fcompiler) -quiet
	./rtrivac multigroup_albedo_3d.x2m $(fcompiler) -quiet
	./rtrivac takedaM4_spn.x2m $(fcompiler) -quiet
	./rtrivac prob5p12_nem.x2m $(fcompiler) -quiet
