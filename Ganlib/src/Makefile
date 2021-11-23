#---------------------------------------------------------------------------
#
#  Makefile for building the Ganlib library and load module
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
opt = -O
PREPRO = cpp
ifeq ($(mpi),1)
  FMPI = -DMPI
  fcompiler = mpif90
  ccompiler = mpicc
else
  FMPI =
  ccompiler = gcc
endif
ifeq ($(openmp),1)
  COMP = -fopenmp
  ccompiler = gcc-11
  PREPRO = cpp -DOPENMP
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
  DIRNAME = AIX
  F90 = xlf90
  C = cc
  FLAGS = -DF90 -DAIX -DUnix
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
bin = ../bin/$(DIRNAME)
lib_module = ../lib/$(DIRNAME)/modules

ifeq ($(intel),1)
  fcompiler = ifort
  FFLAGS = -fPIC
  FFLAG77 = -fPIC
else
  fcompiler = gfortran
endif

SRCC = $(shell ls *.c)
SRC77 = $(shell ls *.f)
SRCF77 = $(shell ls *.F)
ifeq ($(python_version_major),2)
  SRC90 = $(shell python ../../script/make_depend.py ls *.f90)
else
  SRC90 = $(shell python3 ../../script/make_depend_py3.py ls *.f90)
endif
SRCF90 = $(shell ls *.F90)
OBJC = $(SRCC:.c=.o)
all : Ganlib
ifeq ($(mpi),1)
	@echo 'Ganlib: mpi is defined'
endif
ifeq ($(openmp),1)
	@echo 'Ganlib: openmp is defined'
endif
ifeq ($(intel),1)
	@echo 'Ganlib: intel is defined'
endif
%.o : %.c
	$(C) $(CFLAGS) $(FLAGS) $(opt) $(COMP) -c $< -o $@
OBJ90 = $(SRC90:.f90=.o)
%.o : %.f90
	$(F90) $(FFLAGS) $(opt) $(COMP) -c $< -o $@
OBJF90 = $(SRCF90:.F90=.o)
%.o : %.F90
	$(PREPRO) -P -W -traditional $(FLAGS) $< temp.f90
	$(F90) $(FFLAGS) $(opt) $(COMP) -c temp.f90 -o $@
	/bin/rm temp.f90
OBJ77 = $(SRC77:.f=.o)
%.o : %.f
	$(F90) $(FFLAG77) $(opt) $(COMP) -c $< -o $@
OBJF77 = $(SRCF77:.F=.o)
%.o : %.F
	$(PREPRO) -P -W -traditional $(FLAGS) $(FMPI) $< temp.f
	$(F90) $(FFLAG77) $(opt) $(COMP) -c temp.f -o $@
	/bin/rm temp.f
$(lib_module)/:
	mkdir -p $(lib_module)/
libGanlib.a: $(OBJC) $(OBJ90) $(OBJF90) $(OBJ77) $(OBJF77) $(lib_module)/
	ar r $@ $(OBJC) $(OBJ90) $(OBJF90) $(OBJ77) $(OBJF77)
	cp $@ $(lib)/$@
	cp *.mod $(lib_module)
$(bin)/:
	mkdir -p $(bin)/
Ganlib: libGanlib.a GANMAIN.o $(bin)/
	$(F90) $(LFLAGS) $(opt) $(COMP) GANMAIN.o $(lib)/libGanlib.a -o Ganlib
	cp $@ $(bin)/$@
clean:
	/bin/rm -f *.o *.mod *.a Ganlib
