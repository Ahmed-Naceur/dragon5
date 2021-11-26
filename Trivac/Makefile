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
all :
	$(MAKE) -C src
clean :
	$(MAKE) clean -C src
tests :
	./rtrivac DIFtst.x2m -quiet
	./rtrivac Ktests.x2m -quiet
	./rtrivac SPNtst.x2m -quiet
	./rtrivac iaea3d.x2m -quiet
	./rtrivac monju3D_spn.x2m -quiet
	./rtrivac multigroup_albedo_2d.x2m -quiet
	./rtrivac multigroup_albedo_3d.x2m -quiet
	./rtrivac takedaM4_spn.x2m -quiet
