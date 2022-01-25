*--------------------------------------------------------------------- 
*  Electron Source on 1D Slab Geometry
*  CEPXS-BFP library
*  Author: Ahmed Naceur, Sep.2021
*---------------------------------------------------------------------
*  Define STRUCTURES and MODULES used
*---------------------------------------------------------------------

MODULE       GEO: FMAC: MAC: SNT: ASM: FLU: OUT: DELETE: END: ;
PROCEDURE    assertS ;

LINKED_LIST  MACRO OLDMACRO GEOM TRACK SYSTEM FLUX EEDITS1 GWT0 ;
SEQ_ASCII    FMAC_M :: FILE 
'/home/ahnac/dragon-radiotherapy/cepxs-dragon/FMAC_e_3.0MeV_0g_80e_P8_NC.txt' ;


*---------------------------------------------------------------------
*  Set the GEOMETRY
*---------------------------------------------------------------------
* lx=2  : number of subdivisions along x
* VOID  : bc-zero re-entrant flux
* MIX   : mixture by subdivisions (=lx)
* MESHX : spatial mesh along x (in cm)
* SPLITX: mesh splitting 
*---------------------------------------------------------------------
GEOM := GEO: :: CAR1D 2
           X- VOID  X+ VOID
           MIX  1 2
           MESHX 0.0 0.3 4.0
           SPLITX 30 70
           ;

*---------------------------------------------------------------------
*  Track the GEOMETRY
*---------------------------------------------------------------------
* MAXR  :
* DIAM  : spatial approximation order
*          =0: (constant) HODD
*          =1: (linear)   Disc. Galerkin
*          =2: (parabolic) 
*          =3: (cubic) Disc. Galerkin
* NLIVO : disable Livolant acceleration of scattering iteration
* SN    : angular approximation order of the flux
* SCAT  : limit of anistropy of the scattering source (=SN)
* EPSI  : convergence criterion of inner iterations (def. =1.0E-5)
* QUAD  : type of the angular quadrature (1,2,3,4,5,10)
* MAXI  : maximum number of inner iterations (def.=100)
* BFPL  : solution of the BFP with Przybylski and Ligou energy  
*         propagation factors
*---------------------------------------------------------------------          
TRACK := SNT: GEOM ::
      TITLE 'HOMOGENIZED W_6MeV_e_CSD BENCHMARK WITH FMAC-M'
      EDIT 2 MAXR 1000 DIAM 1 NLIVO SN 64 SCAT 64 EPSI 1.E-4
      QUAD 1 MAXI 300 BFPL ;


*---------------------------------------------------------------------
*  Access LIBRARY
*---------------------------------------------------------------------
* FMAC: extract macroscopic cross section data from FMAC-M ASCII file
*   FMAC_M  : .txt (character*12): ASCII file containing FMAC-M information
*   PARTICLE: type of particle irradiating (NEUT,GAMA,BETA,POSITR,PROT)
*
* MAC:
*   MACRO is the MACROLIB updated
*   MIX to be updated (2 is created, while 1 is the mixture in the old
*                      macrolib mixture)
*   OLDL says take cross-sections associated with 1 from OLDMACRO
*
* MAC: 
*   READ  to specify the input file format
*   INPUT to specify that mixture cross sections will be read on the input
*         stream
*   MIX   to specify that the macro. cross.sect. associated with a new 
*         mixture are to be read.
*      1: identifier for the next mixture to be read
*   FIXE  to specify fixed particle source density
*         <<qsour>> 0.0 ... 
*          multigroup array of fixed particle source densiy in [cm-3 s-1]
*---------------------------------------------------------------------


MACRO    := FMAC: FMAC_M :: EDIT 2 PARTICLE B ;


OLDMACRO := MACRO ;

MACRO    := MAC: MACRO OLDMACRO :: MIX 2 1 OLDL ;

OLDMACRO := DELETE: OLDMACRO ;

REAL qsour := 1.E10 2.0 / $Pi_R / ;
ECHO "qsour=" qsour ;

MACRO := MAC: MACRO ::
     READ INPUT
     MIX 1 FIXE
     <<qsour>> 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
     0.0 0.0 0.0 0.0 0.0
     ;

*---------------------------------------------------------------------
*  Compute the FLUX
*---------------------------------------------------------------------
* ASM: assembly module for matrices required by the flux solution module
*   ARM: assembly calculation without building the full collision probability matrices
* FLU: solve the BFP equation
*   TYPE S: fixed source problem
*   EXTE:   control parameters for the external iteration
*      100: maximum number of external iterations
*---------------------------------------------------------------------

SYSTEM := ASM: MACRO TRACK :: EDIT 2 ARM ;
FLUX   := FLU: MACRO TRACK SYSTEM :: EDIT 1 TYPE S EXTE 100 ;


EEDITS1 := OUT: FLUX TRACK MACRO GEOM ::
    EDIT 2 COND INTG NONE ! 1/(\ro_{W} g/cm^{3})
    ;


*----
*  Asserts
*----
* GWT0 := EEDITS1 :: STEP UP GROUP STEP AT 1 ;
* assertS GWT0 :: 'FLUX-INTG' 1 1.649578E-4 ;
* assertS GWT0 :: 'FLUX-INTG' 5 3.512070E-4 ;
* assertS GWT0 :: 'FLUX-INTG' 10 4.093370E-4 ;
* assertS GWT0 :: 'FLUX-INTG' 15 4.156484E-4 ;
* assertS GWT0 :: 'FLUX-INTG' 20 4.124526E-4 ;
* assertS GWT0 :: 'FLUX-INTG' 25 3.760881E-4 ;
* assertS GWT0 :: 'FLUX-INTG' 30 2.358792E-4 ;
*

END: ;


