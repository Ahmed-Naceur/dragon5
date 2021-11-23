#---------------------------------------------------------------------------
#
#  Makefile for executing the Skin++ non-regression test
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
	./rskin -quiet
