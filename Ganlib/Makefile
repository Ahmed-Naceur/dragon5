#---------------------------------------------------------------------------
#
#  Makefile for executing the Ganlib non-regression tests
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
	./rganlib testgan1.x2m $(fcompiler) -quiet
	./rganlib testgan2.x2m $(fcompiler) -quiet
	./rganlib testgan3.x2m $(fcompiler) -quiet
