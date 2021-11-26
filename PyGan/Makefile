#---------------------------------------------------------------------------
#
#  Makefile for executing the Pygan non-regression tests
#  Author : A. Hebert (2020-7-3)
#
#---------------------------------------------------------------------------
#
OS = $(shell uname -s | cut -d"_" -f1)
ifneq (,$(filter $(OS),SunOS AIX))
  MAKE = gmake
endif
all :
	$(MAKE) donjon -C src
ganlib :
	$(MAKE) ganlib -C src
trivac :
	$(MAKE) trivac -C src
dragon :
	$(MAKE) dragon -C src
donjon :
	$(MAKE) donjon -C src
clean :
	$(MAKE) clean -C src
tests :
	./rpython test_lcm.py -quiet
	./rpython test_lifo.py -quiet
	./rpython fact.py -quiet
	./rpython simplePOW.py -quiet
