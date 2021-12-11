#---------------------------------------------------------------------------
#
#  Makefile for building the Dragon library and load module
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
#
nbit = -m32
ARCH = $(shell uname -m)
ifeq (,$(filter $(ARCH),i386 i686))
  nbit = -m64
endif

DIRNAME = $(shell uname -sm | sed 's/[ ]/_/')
OS = $(shell uname -s | cut -d"_" -f1)
ccompiler = gcc
opt = -O
PREPRO = cpp
ifeq ($(openmp),1)
  COMP = -fopenmp
  PREPRO = cpp -D_OPENMP
  ccompiler = gcc-11
else
  COMP =
endif

ifeq ($(intel),1)
  fcompiler = ifort
else
  fcompiler = gfortran
endif

ifeq ($(OS),AIX)
  python_version_major := 2
else
  python_version_full := $(wordlist 2,4,$(subst ., ,$(shell python --version 2>&1)))
  python_version_major := $(word 1,${python_version_full})
endif

ifeq ($(OS),Darwin)
  F90 = $(fcompiler)
  C = $(ccompiler)
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
  FFLAGS = -Wall $(nbit) -fPIC -ffpe-summary=none
  FFLAG77 = -Wall $(nbit) -fPIC -ffpe-summary=none -frecord-marker=4
  LFLAGS = $(nbit)
else
ifeq ($(OS),Linux)
  F90 = $(fcompiler)
  C = $(ccompiler)
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
  FFLAGS = -Wall $(nbit) -fPIC
  FFLAG77 = -Wall $(nbit) -fPIC -frecord-marker=4
  LFLAGS = $(nbit)
else
ifeq ($(OS),CYGWIN)
  F90 = $(fcompiler)
  C = $(ccompiler)
  FLAGS = -DLinux -DUnix
  CFLAGS = -Wall $(nbit) -fPIC
  FFLAGS = -Wall $(nbit) -fPIC
  FFLAG77 = -Wall $(nbit) -fPIC -frecord-marker=4
  LFLAGS = $(nbit)
else
ifeq ($(OS),SunOS)
  MAKE = gmake
  F90 = f90
  C = cc
  PREPRO = /usr/lib/cpp
  FLAGS = -DSunOS -DUnix
  CFLAGS = $(nbit)
  FFLAGS = $(nbit) -s -ftrap=%none
  FFLAG77 = $(nbit) -s -ftrap=%none
  LFLAGS = $(nbit)
else
ifeq ($(OS),AIX)
  opt = -O4
  MAKE = gmake
  DIRNAME = AIX
  F90 = xlf90
  C = cc
  FLAGS = -DAIX -DUnix
  CFLAGS = -qstrict
  FFLAGS = -qstrict -qmaxmem=-1 -qsuffix=f=f90
  FFLAG77 = -qstrict -qmaxmem=-1 -qxlf77=leadzero -qfixed
  LFLAGS = -qstrict -bmaxdata:0x80000000 -qipa
else
  $(error $(OS) is not a valid OS)
endif
endif
endif
endif
endif

ifeq ($(intel),1)
  FFLAGS = -fPIC
  FFLAG77 = -fPIC
  lib = ../lib/$(DIRNAME)_intel
  libUtl = ../../Utilib/lib/$(DIRNAME)_intel
  libGan = ../../Ganlib/lib/$(DIRNAME)_intel
  libTri = ../../Trivac/lib/$(DIRNAME)_intel
  bin = ../bin/$(DIRNAME)_intel
  lib_module = ../lib/$(DIRNAME)_intel/modules
  INCLUDE = -I../../Ganlib/lib/$(DIRNAME)_intel/modules/
else
  lib = ../lib/$(DIRNAME)
  libUtl = ../../Utilib/lib/$(DIRNAME)
  libGan = ../../Ganlib/lib/$(DIRNAME)
  libTri = ../../Trivac/lib/$(DIRNAME)
  bin = ../bin/$(DIRNAME)
  lib_module = ../lib/$(DIRNAME)/modules
  INCLUDE = -I../../Ganlib/lib/$(DIRNAME)/modules/
endif

ifeq ($(hdf5),1)
  LFLAGS += -L${HDF5_DIR}/lib -lhdf5
endif

SRCC = $(shell ls *.c)
SRC77 = $(shell ls *.f)
SRCF77 = $(shell ls *.F)
ifeq ($(python_version_major),2)
  SRC90 = $(shell python ../../script/make_depend.py ls *.f90)
else
  SRC90 = $(shell python3 ../../script/make_depend_py3.py ls *.f90)
endif
OBJC = $(SRCC:.c=.o)
OBJ90 = $(SRC90:.f90=.o)
OBJ77 = $(SRC77:.f=.o)
OBJF77 = $(SRCF77:.F=.o)
all : sub-make Dragon
ifeq ($(openmp),1)
	@echo 'Dragon: openmp is defined'
endif
ifeq ($(intel),1)
	@echo 'Dragon: intel is defined'
endif
ifeq ($(hdf5),1)
	@echo 'Dragon: hdf5 is defined'
endif
sub-make:
	$(MAKE) openmp=$(openmp) intel=$(intel) hdf5=$(hdf5) -C ../../Trivac/src
%.o : %.c
	$(C) $(CFLAGS) $(opt) $(COMP) -c $< -o $@
%.o : %.f90
	$(F90) $(FFLAGS) $(opt) $(COMP) $(INCLUDE) -c $< -o $@
%.o : %.f
	$(F90) $(FFLAG77) $(opt) $(COMP) $(INCLUDE) -c $< -o $@
%.o : %.F
	$(PREPRO) -P -W -traditional  $(FLAGS) $< temp.f
	$(F90) $(FFLAG77) $(opt) $(COMP) $(INCLUDE) -c temp.f -o $@
	/bin/rm temp.f
$(lib_module)/:
	mkdir -p $(lib_module)/
libDragon.a: $(OBJC) $(OBJ90) $(OBJ77) $(OBJF77) $(lib_module)/
	ar r $@ $(OBJC) $(OBJ90) $(OBJ77) $(OBJF77)
	cp $@ $(lib)/$@
	cp *.mod $(lib_module)
$(bin)/:
	mkdir -p $(bin)/
Dragon: libDragon.a DRAGON.o $(bin)/ sub-make
	$(F90) $(LFLAGS) $(opt) $(COMP) DRAGON.o $(lib)/libDragon.a $(libTri)/libTrivac.a \
	$(libUtl)/libUtilib.a $(libGan)/libGanlib.a -o Dragon
	cp $@ $(bin)/$@
clean:
	$(MAKE) -C ../../Trivac/src clean
	/bin/rm -f *.o *.mod *.a sub-make temp.* Dragon
