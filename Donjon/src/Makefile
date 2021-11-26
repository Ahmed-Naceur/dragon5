#---------------------------------------------------------------------------
#
#  Makefile for building the Donjon library and load module
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
fcompiler = gfortran
ccompiler = gcc
opt = -O
ifeq ($(openmp),1)
  COMP = -fopenmp
  ccompiler = gcc-11
else
  COMP =
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
lib = ../lib/$(DIRNAME)
libUtl = ../../Utilib/lib/$(DIRNAME)
libGan = ../../Ganlib/lib/$(DIRNAME)
libTri = ../../Trivac/lib/$(DIRNAME)
libDra = ../../Dragon/lib/$(DIRNAME)
bin = ../bin/$(DIRNAME)
INCLUDE = -I../../Ganlib/lib/$(DIRNAME)/modules/

ifeq ($(intel),1)
  fcompiler = ifort
  FFLAGS = -fPIC
  FFLAG77 = -fPIC
else
  fcompiler = gfortran
endif

SRCC = $(shell ls *.c)
SRC77 = $(shell ls *.f)
ifeq ($(python_version_major),2)
  SRC90 = $(shell python ../../script/make_depend.py ls *.f90)
else
  SRC90 = $(shell python3 ../../script/make_depend_py3.py ls *.f90)
endif
OBJC = $(SRCC:.c=.o)
all : sub-make Donjon
ifeq ($(openmp),1)
	@echo 'Donjon: openmp is defined'
endif
ifeq ($(intel),1)
	@echo 'Donjon: intel is defined'
endif
sub-make:
	$(MAKE) openmp=$(openmp) intel=$(intel) -C ../../Dragon/src
%.o : %.c
	$(C) $(CFLAGS) $(opt) $(COMP) -c $< -o $@
OBJ90 = $(SRC90:.f90=.o)
%.o : %.f90
	$(F90) $(FFLAGS) $(opt) $(COMP) $(INCLUDE) -c $< -o $@
OBJ77 = $(SRC77:.f=.o)
%.o : %.f
	$(F90) $(FFLAG77) $(opt) $(COMP) $(INCLUDE) -c $< -o $@
$(lib)/:
	mkdir -p $(lib)/
libDonjon.a: $(OBJC) $(OBJ90) $(OBJ77) $(lib)/
	ar r $@ $(OBJC) $(OBJ90) $(OBJ77)
	cp $@ $(lib)/$@
$(bin)/:
	mkdir -p $(bin)/
Donjon: libDonjon.a DONJON.o $(bin)/ sub-make
	$(F90) $(LFLAGS) $(opt) $(COMP) DONJON.o $(lib)/libDonjon.a $(libDra)/libDragon.a \
	$(libTri)/libTrivac.a $(libUtl)/libUtilib.a $(libGan)/libGanlib.a -o Donjon
	cp $@ $(bin)/$@
clean:
	$(MAKE) -C ../../Dragon/src clean
	/bin/rm -f *.o *.mod *.a Donjon
