*DECK SEN
      SUBROUTINE SEN(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To create sensitivity profiles to cross-section on the reactivity
* using first order perturbation method based on the 
* adjoint calculation.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Laville, G. Marleau
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1) create or modification type(L_SENS);
*         HENTRY(2) read-only type(L_MACROLIB or L_LIBRARY);
*         HENTRY(3) read-only type(L_TRACK);
*         HENTRY(4) read-only type(L_FLUX);
*         HENTRY(5) read only type(L_AFLUX).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*Comments:
*  Call format:
*    File.sdf := SENS: Flux Adjoint Biblio Track
*                      :: [ EDIT iprint ANIS nanis ] ;
*  with
*    File.sdf = sdf (SEQ_ASCII) file in creation mode     
*    Flux     = Flux (LINKED_LIST or XSM_FILE) in read only mode     
*    Adjoint  = Adjoint (LINKED_LIST or XSM_FILE) in read only mode    
*    Biblio   = Biblio (LINKED_LIST or XSM_FILE) in read only mode    
*    Track    = Track (LINKED_LIST or XSM_FILE) in read only mode   
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT      NONE
*----
*  Routine arguments
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  Parameters
*----
      INTEGER      NCAR,NSTATE,IOUT
      CHARACTER    NAMSBR*6,HSIGN*12,SDF*12
      PARAMETER   (NSTATE=40,IOUT=6,NAMSBR='SEN   ')
      PARAMETER   (NCAR=3)
*----
*  Local variables
*----
      TYPE(C_PTR) IPLIB,IPTRACK,IPFLUX,IPAFLUX
      INTEGER     IPSENS,NR,NU,NM,NI,NG,NGD,NGA,NUD,NUA,NMT,NL,IFMT
      INTEGER     I,IEN,ISDF,IADJ,ITYPE,ISTATE(NSTATE),IPRINT,
     >            NANIS,NLTERM
*----
*  Verify if call format is adequate
*----
      IF(NENTRY .NE. 5) CALL XABORT(NAMSBR//
     > ': FIVE data structure EXPECTED.')
*----
*  First data structure .sdf file
*----
      IEN=1
      IF(IENTRY(IEN) .NE. 4 ) CALL XABORT(NAMSBR//
     >  ': SEQ_ASCII format expected for .sdf file')
      IF(JENTRY(IEN) .NE. 0 ) CALL XABORT(NAMSBR//
     >  ': .sdf file must be in creation mode')
      SDF=HENTRY(IEN)
      ISDF=0
      DO I=1,9
         IF(SDF(I:I+3).EQ.'.sdf') ISDF=1
      ENDDO
      IF(ISDF.NE.1) CALL XABORT(NAMSBR//
     >  ': The extension of the first structure has be ".sdf"')
      IPSENS=FILUNIT(KENTRY(IEN))
*----
*  Process the other 4 data structures (arbitrary order)
*----
      IPLIB=C_NULL_PTR
      IPTRACK=C_NULL_PTR
      IPFLUX=C_NULL_PTR
      IPAFLUX=C_NULL_PTR
      NUD=0
      NUA=0
      NMT=0
      NGD=0
      NGA=0
      DO IEN=2,5
        IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2))
     >  CALL XABORT(NAMSBR//
     >  ': LINKED_LIST or XSM_FILE expected')
        IF(JENTRY(IEN).NE.2) CALL XABORT(NAMSBR//
     >  ': data structure must be in READ_ONLY mode')
        CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.EQ.'L_FLUX') THEN
          CALL LCMGET(KENTRY(IEN),'STATE-VECTOR',ISTATE)
          ITYPE=ISTATE(3)        
          IF((ITYPE.NE.1).AND.(ITYPE.NE.10)) CALL XABORT(NAMSBR//
     >    ': Keff problem required')
          IADJ=MOD(ISTATE(3)/10,10)
          IF(IADJ .EQ. 1) THEN
            IPAFLUX=KENTRY(IEN)
            NGA=ISTATE(1)
            NUA=ISTATE(2)
          ELSE
            IPFLUX=KENTRY(IEN)
            NGD=ISTATE(1)
            NUD=ISTATE(2)
          ENDIF
        ELSE IF(HSIGN.EQ.'L_TRACK') THEN
          CALL LCMGET(KENTRY(IEN),'STATE-VECTOR',ISTATE)
          IPTRACK=KENTRY(IEN)
          NR=ISTATE(1)
          NU=ISTATE(2)
          NMT=ISTATE(4)
        ELSE IF(HSIGN.EQ.'L_LIBRARY') THEN
          CALL LCMGET(KENTRY(IEN),'STATE-VECTOR',ISTATE)
          IPLIB=KENTRY(IEN)
          NM=ISTATE(1)
          NI=ISTATE(2)
          NG=ISTATE(3)
          NL=ISTATE(4)
          IFMT=ISTATE(5)
        ELSE
          CALL XABORT(NAMSBR//': '//HSIGN//' is an invalid signature ')
        ENDIF
      ENDDO 
*----
*  Test if all data structures required are available
*----
      IF(.NOT.C_ASSOCIATED(IPLIB)) CALL XABORT(NAMSBR//
     > ': No microlib data structure found') 
      IF(.NOT.C_ASSOCIATED(IPTRACK)) CALL XABORT(NAMSBR//
     > ': No tracking data structure found') 
      IF(.NOT.C_ASSOCIATED(IPFLUX)) CALL XABORT(NAMSBR//
     > ': No direct flux data structure found') 
      IF(.NOT.C_ASSOCIATED(IPAFLUX)) CALL XABORT(NAMSBR//
     > ': No adjoint flux data structure found') 
*----
*  Test if parameters are compatibles
*  NR Number of region in Tracking object.
*  NU Number of unkwnow in Tracking/Flux objects.
*  NM Number of mixture in Library object.
*  NI Number of isotopes in Library object.
*  NG Number of energy group in Library object.
*----
      IF(NGD .NE. NG) CALL XABORT(NAMSBR//
     > ': Number of groups in flux and microlib not identical')
      IF(NGA .NE. NG) CALL XABORT(NAMSBR//
     > ': Number of groups in adjoint and microlib not identical')
      IF(NUD .NE. NU) CALL XABORT(NAMSBR//
     > ': Number of unknowns in flux and tracking not identical')
      IF(NUA .NE. NU) CALL XABORT(NAMSBR//
     > ': Number of unknowns in adjoint and tracking not identical')
      IF(NMT .GT. NM) CALL XABORT(NAMSBR//
     > ': Number of mixtures in tracking larger that microlib')
*----
*  Read input parameters
*----
      NANIS=1
      CALL SENGET(IPRINT,NL,NANIS)
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6000) NG,NR,NU,NM,NI,NL,NANIS
      ENDIF     
*----
*  Launch sensitivity analysis main routine.
*----
*                            ! A 1D calculation have to use NANIS=1
*      IF(NDIM .EQ. 3) THEN  ! It is necessary to introduce the parameter NDIM
*        NLTERM=NANIS*NANIS  ! 3D calculation
*      ELSEIF(NDIM .EQ. 2) THEN
        NLTERM=(NANIS*(NANIS+1))/2 ! 2D calculation
*      ELSE
*        NLTERM=NANIS        ! 1D calculation
*      ENDIF
      CALL SENDRV(IPSENS,IPTRACK,IPLIB,IPFLUX,IPAFLUX,IPRINT,
     >            NR,NU,NI,NG,NANIS,NLTERM)
      RETURN
*----
*  Format
*----
 6000 FORMAT(' Number of groups           =',I10/
     >       ' Number of regions          =',I10/
     >       ' Number of unknowns         =',I10/
     >       ' Maximum number of mixtures =',I10/
     >       ' Number of isotopes         =',I10/
     >       ' Number anisotropy order    =',I10/
     >       ' Anisotropy order kept      =',I10)
      END
