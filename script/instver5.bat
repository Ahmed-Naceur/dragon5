rem *********************************************************************************
rem   Script to compile and link Version4 components on MS-DOS
rem   compile with intel Fortran and intel C++
rem   run the script in Dos windows from Visual studio command prompt
rem   execute script from Version4 directory as .\script\instver5
rem   base on a script created by E. Varin for Dragon 3.06
rem *********************************************************************************
rem ----compile utilib----
cd Utilib
md lib\ms-dos
copy src\*.* lib\ms-dos\
cd lib\ms-dos\
ifort /c /assume:byterecl /warn:nofileopt *.f 1>OUT 2>&1
rem make utilib library
lib /nologo /out:utilib.lib *.obj
del *.obj *.f
cd ..\..\..
rem ----------------------
rem ----compile ganlib----
cd Ganlib
md lib\ms-dos
copy src\*.* lib\ms-dos\
cd lib\ms-dos\
rem ren KDRCPU.F KDRCPU.fpp
rem ren KDROPN.F KDROPN.fpp
rem ren KDRSTD.F KDRSTD.fpp
ren DRVMPI.F DRVMPI.fpp
ren SNDMPI.F SNDMPI.fpp
ifort /c /DMSDOS /assume:byterecl /warn:nofileopt FILMODx.f90  1>XOUT0 2>&1
ifort /c /DMSDOS /assume:byterecl /warn:nofileopt GANLIBx.f90  1>XOUT1 2>&1
ifort /c /DMSDOS /assume:byterecl /warn:nofileopt *.f *.fpp *.f90  1>XOUT2 2>&1
cl /c /DMSDOS *.c
rem make ganlib library
move ganmain.obj ..\
lib /nologo /out:ganlib.lib *.obj
del *.obj *.f *.f90 *.fpp *.c *.h
move ..\ganmain.obj .
cd ..\..
rem make ganlib executable
md bin\ms-dos
cd bin\ms-dos\
rem
ifort -o ganlib.exe ..\..\lib\ms-dos\ganmain.obj ..\..\lib\ms-dos\ganlib.lib ^
   ..\..\..\Utilib\lib\ms-dos\utilib.lib 1>XOUT3 2>&1
cd ..\..\..
rem ----------------------
rem ----compile trivac----
cd Trivac
md lib\ms-dos
copy src\*.* lib\ms-dos\
cd lib\ms-dos\
copy ..\..\..\ganlib\src\filmod.f90 filmod.f90
copy ..\..\..\ganlib\src\ganlib.f90 ganlib.f90
ifort /c /DMSDOS /assume:byterecl /warn:nofileopt filmod.f90  1>OUT0 2>&1
ifort /c /DMSDOS /assume:byterecl /warn:nofileopt ganlib.f90  1>OUT1 2>&1
ifort /c /assume:byterecl /warn:nofileopt *.f90 *.f 1>OUT2 2>&1
rem make trivac library
move trivac.obj ..\
lib /nologo /out:trivac.lib *.obj >OUT3
del *.obj *.f *.f90
move ..\trivac.obj .
cd ..\..
rem make trivac executable
md bin\ms-dos
cd bin\ms-dos\
ifort -o trivac.exe ..\..\lib\ms-dos\trivac.obj ..\..\lib\ms-dos\trivac.lib ^
   ..\..\..\Ganlib\lib\ms-dos\ganlib.lib ..\..\..\Utilib\lib\ms-dos\utilib.lib ^
   1>OUT3 2>&1
cd ..\..\..
rem ----------------------
rem ----compile dragon----
cd Dragon
md lib\ms-dos
copy src\*.* lib\ms-dos\
cd lib\ms-dos\
copy ..\..\..\ganlib\src\filmod.f90 filmod.f90
copy ..\..\..\ganlib\src\ganlib.f90 ganlib.f90
rem ren DRAGON.F DRAGON.fpp
ifort /c /DMSDOS /assume:byterecl /warn:nofileopt filmod.f90  1>OUT0 2>&1
ifort /c /DMSDOS /assume:byterecl /warn:nofileopt ganlib.f90  1>OUT1 2>&1
ifort /c /assume:byterecl /warn:nofileopt *.f90 *.f 1>OUT2 2>&1
cl /c /DMSDOS *.c   1>OUT3 2>&1
rem make dragon library
move dragon.obj ..\
lib /nologo /out:dragon.lib *.obj  1>OUT4 2>&1
del *.obj *.f90 *.f *.c *.h
move ..\dragon.obj .
cd ..\..
rem make dragon executable
md bin\ms-dos
cd bin\ms-dos\
ifort -o dragon.exe ..\..\lib\ms-dos\dragon.obj ..\..\lib\ms-dos\dragon.lib ^
   ..\..\..\Trivac\lib\ms-dos\trivac.lib ..\..\..\Ganlib\lib\ms-dos\ganlib.lib ^
   ..\..\..\Utilib\lib\ms-dos\utilib.lib  1>OUT4 2>&1
cd ..\..\..
rem ----------------------
rem ----compile donjon----
cd Donjon
md lib\ms-dos
copy src\*.* lib\ms-dos\
cd lib\ms-dos\
rem ren DONJON.F DONJON.fpp
copy ..\..\..\ganlib\src\filmod.f90 filmod.f90
copy ..\..\..\ganlib\src\ganlib.f90 ganlib.f90
ifort /c /DMSDOS /assume:byterecl /warn:nofileopt filmod.f90  1>OUT0 2>&1
ifort /c /DMSDOS /assume:byterecl /warn:nofileopt ganlib.f90  1>OUT1 2>&1
ifort /c /assume:byterecl /warn:nofileopt *.f90 *.f 1>OUT2 2>&1
cl /c /DMSDOS *.c  1>OUT3 2>&1
rem make donjon library
move donjon.obj ..\
lib /nologo /out:donjon.lib *.obj  1>OUT4 2>&1
del *.obj *.f90 *.f *.c *.h
move ..\donjon.obj .
cd ..\..
rem make donjon executable
md bin\ms-dos
cd bin\ms-dos\
ifort -o donjon.exe ..\..\lib\ms-dos\donjon.obj ..\..\lib\ms-dos\donjon.lib ^
   ..\..\..\Dragon\lib\ms-dos\dragon.lib ..\..\..\Trivac\lib\ms-dos\trivac.lib ^
   ..\..\..\Ganlib\lib\ms-dos\ganlib.lib ..\..\..\Utilib\lib\ms-dos\utilib.lib ^
    1>OUT4 2>&1
cd ..\..\..
rem ---------------------
rem ----compile optex----
rem cd Optex
rem md lib\ms-dos
rem copy src\*.* lib\ms-dos\
rem cd lib\ms-dos\
rem ifort /c /assume:byterecl /warn:nofileopt *.f
rem rem make optex library
rem move optex.obj ..
rem lib /nologo /out:optex.lib *.obj
rem del *.obj *.f
rem move ..\optex.obj .
rem cd ..\..
rem rem make optex executable
rem md bin\ms-dos
rem cd bin\ms-dos\
rem ifort -o optex.exe ..\..\lib\ms-dos\optex.obj ..\..\lib\ms-dos\optex.lib ^
rem    ..\..\..\Donjon\lib\ms-dos\donjon.lib  ..\..\..\Dragon\lib\ms-dos\dragon.lib ^
rem    ..\..\..\Trivac\lib\ms-dos\trivac.lib ..\..\..\Ganlib\lib\ms-dos\ganlib.lib ^
rem    ..\..\..\Utilib\lib\ms-dos\utilib.lib
rem cd ..\..\..
rem *********************************************************************************
