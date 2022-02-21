#---------------------------------------------------------------------------
#
#  Makefile for building the PyGan extension modules
#  Author : A. Hebert (2020-7-3)
#
#---------------------------------------------------------------------------
#
PYTHON3_OK := $(shell command -v python3)
ifeq ('$(PYTHON3_OK)','')
    $(error package 'python3' not found)
endif
ifeq ($(FORTRANPATH),)
  $(error FORTRANPATH is not set)
endif

DIRNAME = $(shell uname -sm | sed 's/[ ]/_/')
OS = $(shell uname -s | cut -d"_" -f1)
lib = ../lib/$(DIRNAME)
libGan = ../../Ganlib/lib/$(DIRNAME)
pylib = $(shell python3 ../../script/find_pylib.py)

all:
	$(MAKE) donjon
ganlib: clean sub-make-ganlib pygan-ganlib
ifdef openmp
	@echo 'pygan_ganlib: openmp is defined'
endif
ifdef hdf5
	@echo 'pygan_ganlib: hdf5 is defined'
endif
trivac: clean sub-make-trivac pygan-trivac
ifdef openmp
	@echo 'pygan_trivac: openmp is defined'
endif
ifdef hdf5
	@echo 'pygan_trivac: hdf5 is defined'
endif
dragon: clean sub-make-dragon pygan-dragon
ifdef openmp
	@echo 'pygan_dragon: openmp is defined'
endif
ifdef hdf5
	@echo 'pygan_dragon: hdf5 is defined'
endif
donjon: clean sub-make-donjon pygan-donjon
ifdef openmp
	@echo 'pygan_donjon: openmp is defined'
endif
ifdef hdf5
	@echo 'pygan_donjon: hdf5 is defined'
endif
sub-make-ganlib:
	$(MAKE) openmp=$(openmp) -C ../../Utilib/src
	$(MAKE) openmp=$(openmp) hdf5=$(hdf5) -C ../../Ganlib/src
sub-make-trivac: sub-make-ganlib
	$(MAKE) openmp=$(openmp) hdf5=$(hdf5) -C ../../Trivac/src
sub-make-dragon: sub-make-trivac
	$(MAKE) openmp=$(openmp) hdf5=$(hdf5) -C ../../Dragon/src
sub-make-donjon: sub-make-dragon
	$(MAKE) openmp=$(openmp) hdf5=$(hdf5) -C ../../Donjon/src
libGanlib.a: $(lib)/ sub-make-ganlib
	cp $(libGan)/libGanlib.a .
	ar -d libGanlib.a xabort_c.o
	mv libGanlib.a $(lib)
$(lib)/:
	mkdir -p $(lib)/
pygan-ganlib: libGanlib.a sub-make-ganlib
	mkdir -p $(DIRNAME)
	cp *.[ch] $(DIRNAME)
	cp *.py $(DIRNAME)
	cd $(DIRNAME); python3 setup_lcm.py install --home=.
	cd $(DIRNAME); python3 setup_lifo.py install --home=.
ifdef openmp
	export CODE_EMBEDDED=GANLIB_OMP; cd $(DIRNAME); python3 setup_cle2000.py install --home=.
else
	export CODE_EMBEDDED=GANLIB; cd $(DIRNAME); python3 setup_cle2000.py install --home=.
endif
	mv $(DIRNAME)/$(pylib)/* $(lib)/
	/bin/rm -r $(DIRNAME)
pygan-trivac: libGanlib.a sub-make-trivac
	mkdir -p $(DIRNAME)
	cp *.[ch] $(DIRNAME)
	cp *.py $(DIRNAME)
	cd $(DIRNAME); python3 setup_lcm.py install --home=.
	cd $(DIRNAME); python3 setup_lifo.py install --home=.
ifdef openmp
	export CODE_EMBEDDED=TRIVAC_OMP; cd $(DIRNAME); python3 setup_cle2000.py install --home=.
else
	export CODE_EMBEDDED=TRIVAC; cd $(DIRNAME); python3 setup_cle2000.py install --home=.
endif
	mv $(DIRNAME)/$(pylib)/* $(lib)/
	/bin/rm -r $(DIRNAME)
pygan-dragon: libGanlib.a sub-make-dragon
	mkdir -p $(DIRNAME)
	cp *.[ch] $(DIRNAME)
	cp *.py $(DIRNAME)
	cd $(DIRNAME); python3 setup_lcm.py install --home=.
	cd $(DIRNAME); python3 setup_lifo.py install --home=.
ifdef openmp
	export CODE_EMBEDDED=DRAGON_OMP; cd $(DIRNAME); python3 setup_cle2000.py install --home=.
else
	export CODE_EMBEDDED=DRAGON; cd $(DIRNAME); python3 setup_cle2000.py install --home=.
endif
	mv $(DIRNAME)/$(pylib)/* $(lib)/
	/bin/rm -r $(DIRNAME)
pygan-donjon: libGanlib.a sub-make-donjon
	mkdir -p $(DIRNAME)
	cp *.[ch] $(DIRNAME)
	cp *.py $(DIRNAME)
	cd $(DIRNAME); python3 setup_lcm.py install --home=.
	cd $(DIRNAME); python3 setup_lifo.py install --home=.
ifdef openmp
	export CODE_EMBEDDED=DONJON_OMP; cd $(DIRNAME); python3 setup_cle2000.py install --home=.
else
	export CODE_EMBEDDED=DONJON; cd $(DIRNAME); python3 setup_cle2000.py install --home=.
endif
	mv $(DIRNAME)/$(pylib)/* $(lib)/
	/bin/rm -r $(DIRNAME)
clean:
	/bin/rm -r -f ../lib
	/bin/rm -f sub-make*
