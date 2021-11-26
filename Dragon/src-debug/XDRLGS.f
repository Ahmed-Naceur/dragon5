*DECK XDRLGS
      SUBROUTINE XDRLGS(IPLIB ,IGS   ,IPRINT,MINLEG,MAXLEG,IORD  ,
     >                  NGROUP,XSREC ,SCAT  ,ITYPRO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Get/save Legendre-dependent cross section data from/on IPLIB.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input
* IPLIB   pointer to the internal library.
* IGS     get or save flag:
*         >0 save;
*         <0 get.
* IPRINT  Print level (cross sections printed if IPRINT>99).
* MINLEG  mimimum Legendre order to process for scattering.
* MAXLEG  maximum Legendre order to process for scattering.
* IORD    cross section order:
*         =1 constant;
*         =2 linear;
*         =3 quadratic.
* NGROUP  number of energy groups.
*
*Parameters: input/output
* XSREC   cross section records (scattering cross section of
*         order MINLEG to MAXLEG for IRPROC=MINLEG+1,MAXLEG+1).
* SCAT    complete scattering matrix (SCAT(JG,IG) is from IG to JG
*         for order MINLEG to MAXLEG).
*
*Parameters: output
* ITYPRO  vector for cross section processed indices:
*         =0 absent  (not processed);
*         >0 present (processed).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPLIB
      INTEGER          IGS,IPRINT,MINLEG,MAXLEG,IORD,NGROUP,
     >                 ITYPRO(MAXLEG-MINLEG+1)
      REAL             XSREC(NGROUP,MAXLEG-MINLEG+1),
     >                 SCAT(NGROUP,NGROUP,MAXLEG-MINLEG+1)
*----
*  LOCAL PARAMETERS
*  SCATTERING CROSS SECTIONS START AT MINLEG+1 WITH
*  NAME NAMSCT='SIGS'//NAMLEG AND NAMSCT='SCAT'//NAMLEG
*  WITH NAMLEG DEFINED BY
*  WRITE(NAMLEG ,'(I2.2)') ILEG
*  FOR ILEG=MINLEG+1 TO MAXLEG+1
*----
      INTEGER          IOUT,MAXGAR
      PARAMETER       (IOUT=6,MAXGAR=100)
      INTEGER          NPROC,IGAR(MAXGAR),IODIV,LONG,ITYP,LONG2,ILEG,
     >                 IXSR,IXSTN,IG,JG,NXSCMP,IGTO,IGMIN,IGMAX,IGFROM
      CHARACTER*12     NAMXS
      CHARACTER        NAMLEG*2,NORD*6,HCM(0:10)*2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,IJJ
      REAL, ALLOCATABLE, DIMENSION(:) :: XSSCMP
      DATA HCM         /'00','01','02','03','04','05','06','07','08',
     >                  '09','10'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NJJ(NGROUP),IJJ(NGROUP),XSSCMP(NGROUP*NGROUP))
*
      IODIV=0
      IF(IORD.EQ.1) THEN
        NORD='      '
        IODIV=1
      ELSE IF(IORD.EQ.2) THEN
        NORD='   LIN'
        IODIV=2
      ELSE IF(IORD.EQ.3) THEN
        NORD='   QUA'
        IODIV=4
      ENDIF
      NPROC=MAXLEG-MINLEG+1
*----
*  READ/INITIALIZE STATE VECTORS
*----
      IF(MAXLEG+1.GT.MAXGAR) THEN
        WRITE(IOUT,9000) 'SCAT-SAVED',MAXGAR,MAXLEG+1
        CALL XABORT('XDRLGS: INVALID VALUE FOR MAXLEG')
      ELSE IF(MAXLEG.LT.MINLEG) THEN
        CALL XABORT('XDRLGS: MAXLEG.LT.MINLEG')
      ENDIF
*
      CALL XDISET(ITYPRO,NPROC,0)
      CALL LCMLEN(IPLIB,'SCAT-SAVED',LONG,ITYP)
      LONG2=MAX(LONG,MAXLEG+1)
      CALL XDISET(IGAR,LONG2,0)
      IF(LONG.NE.0) THEN
        CALL LCMGET(IPLIB,'SCAT-SAVED',IGAR)
        DO 20 ILEG=MINLEG+1,MIN(LONG,MAXLEG+1)
          ITYPRO(ILEG-MINLEG)=IGAR(ILEG)
  20    CONTINUE
      ELSE IF(IGS.LT.0) THEN
        CALL XABORT('XDRLGS: NO SCAT-SAVED RECORD AVAILABLE')
      ENDIF
*
      IF(IGS.GT.0) THEN
*----
*  SAVE LOCAL SCATTERING MATRIX
*----
        IXSR=0
        DO 120 ILEG=MINLEG+1,MAXLEG+1
          IXSR=IXSR+1
          IXSTN=MOD(ITYPRO(IXSR)/IODIV,2)
*----
*  FIND IF SCATTERING XS NOT ALL 0.0
*----
          IF((ILEG.EQ.1).AND.(IODIV.EQ.1)) THEN
            IF(IXSTN.EQ.0) THEN
              ITYPRO(IXSR)=ITYPRO(IXSR)+IODIV
              IGAR(ILEG)=IGAR(ILEG)+IODIV
              IXSTN=1
            ENDIF
          ELSE
            DO 130 IG=1,NGROUP
              DO 131 JG=1,NGROUP
                IF(SCAT(IG,JG,IXSR).NE.0.0) THEN
                  IF(IXSTN.EQ.0) THEN
                    ITYPRO(IXSR)=ITYPRO(IXSR)+IODIV
                    IGAR(ILEG)=IGAR(ILEG)+IODIV
                    IXSTN=1
                  ENDIF
                  GO TO 135
                ENDIF
 131          CONTINUE
 130        CONTINUE
          ENDIF
 135      IF(IXSTN.NE.0) THEN
            IF(ILEG.LE.11) THEN
               NAMLEG=HCM(ILEG-1)
            ELSE
               WRITE(NAMLEG,'(I2.2)') ILEG-1
            ENDIF
            CALL LCMPUT(IPLIB,'SIGS'//NAMLEG//NORD,NGROUP,2,
     >      XSREC(1,IXSR))
*----
*  COMPRESS SCATTERING MATRIX
*  SCAT(IGTO,IGFROM) REPRESENT SCATTERING CROSS SECTION
*    FROM GROUP "IGFROM" TO GROUP "IGTO"
*  IJJ(IGTO) IS MAXIMUM GROUP NUMBER
*    WITH SCATTERING TO "IGTO" GROUP
*  NJJ(IGTO) IS NUMBER OF GROUPS
*    WITH SCATTERING TO "IGTO" GROUP
*  XSSCMP(IX) IS COMPRESSED SCATTERING MATRIX
*  IX CAN BE LOCALIZED IN SCAT(IGTO,IGFROM) USING
*    IF(IGTO=1) THEN
*      IPOSD=1
*    ELSE
*      IPOSD=1+SUM( NJJ(IGF) , IGF=1,IGTO-1)
*    ENDIF
*    IF(IGFROM.GT.IJJ(IGTO)) THEN
*      XSSCMP NOT STORED
*    ELSE IF(IGFROM.LT.IJJ(IGTO)-NJJ(IGTO)+1) THEN
*      XSSCMP NOT STORED
*    ELSE
*      IX=IPOSD+IJJ(IGTO)-IGFROM
*      XSSCMP(IX)=SCAT(IGTO,IGFROM)
*  ENDIF
*----
            NXSCMP=0
            DO 140 IGTO=1,NGROUP
              IGMIN=IGTO
              IGMAX=IGTO
              DO 150 IGFROM=1,NGROUP
                IF(SCAT(IGTO,IGFROM,IXSR).NE.0.0) THEN
                  IGMIN=MIN(IGMIN,IGFROM)
                  IGMAX=MAX(IGMAX,IGFROM)
                ENDIF
 150          CONTINUE
              IJJ(IGTO)=IGMAX
              NJJ(IGTO)=IGMAX-IGMIN+1
              DO 160 IGFROM=IGMAX,IGMIN,-1
                NXSCMP=NXSCMP+1
                XSSCMP(NXSCMP)=SCAT(IGTO,IGFROM,IXSR)
 160          CONTINUE
 140        CONTINUE
            CALL LCMPUT(IPLIB,'NJJS'//NAMLEG//NORD,NGROUP,1,NJJ)
            CALL LCMPUT(IPLIB,'IJJS'//NAMLEG//NORD,NGROUP,1,IJJ)
            CALL LCMPUT(IPLIB,'SCAT'//NAMLEG//NORD,NXSCMP,2,XSSCMP)
          ENDIF
 120    CONTINUE
        CALL LCMPUT(IPLIB,'SCAT-SAVED',LONG2,1,IGAR)
      ELSE
*----
*  GET LOCAL SCATTERING MATRIX
*----
        IXSR=0
        DO 220 ILEG=MINLEG+1,MAXLEG+1
          IXSR=IXSR+1
*----
*  READ IF IXSTN = 1
*  INITIALIZE TO 0.0 IF IXSTN = 0
*----
          CALL XDRSET(XSREC(1,IXSR),NGROUP,0.0)
          CALL XDRSET(SCAT(1,1,IXSR),NGROUP*NGROUP,0.0)
          IXSTN=MOD(ITYPRO(IXSR)/IODIV,2)
          IF(IXSTN.EQ.1) THEN
            IF(ILEG.LE.11) THEN
               NAMLEG=HCM(ILEG-1)
            ELSE
               WRITE(NAMLEG,'(I2.2)') ILEG-1
            ENDIF
            CALL LCMGET(IPLIB,'SIGS'//NAMLEG//NORD,XSREC(1,IXSR))
            CALL LCMGET(IPLIB,'NJJS'//NAMLEG//NORD,NJJ)
            CALL LCMGET(IPLIB,'IJJS'//NAMLEG//NORD,IJJ)
            CALL LCMGET(IPLIB,'SCAT'//NAMLEG//NORD,XSSCMP)
*----
*  DECOMPRESS SCATTERING MATRIX
*  SCAT(IGTO,IGFROM) REPRESENT SCATTERING CROSS SECTION
*    FROM GROUP "IGFROM" TO GROUP "IGTO"
*  IJJ(IGTO) IS MAXIMUM GROUP NUMBER
*    WITH SCATTERING TO "IGTO" GROUP
*  NJJ(IGTO) IS NUMBER OF GROUPS
*    WITH SCATTERING TO "IGTO" GROUP
*  XSSCMP(IX) IS COMPRESSED SCATTERING MATRIX
*  SCAT(IGTO,IGFROM) CAN BE LOCALIZED IN XSSCMP(IX) USING
*    IF(IGTO=1) THEN
*      IPOSD=1
*    ELSE
*      IPOSD=1+SUM( NJJ(IGF) , IGF=1,IGTO-1)
*    ENDIF
*    IF(IGFROM.GT.IJJ(IGTO)) THEN
*      SCAT(IGTO,IGFROM)=0.0
*    ELSE IF(IGFROM.LT.IJJ(IGTO)-NJJ(IGTO)+1) THEN
*      SCAT(IGTO,IGFROM)=0.0
*    ELSE
*      SCAT(IGTO,IGFROM)=XSSCMP(IPOSD+IJJ(IGTO)-IGFROM)
*    ENDIF
*----
            NXSCMP=0
            DO 240 IGTO=1,NGROUP
              IGMAX=IJJ(IGTO)
              IGMIN=IGMAX-NJJ(IGTO)+1
              DO 250 IGFROM=IGMAX,IGMIN,-1
                NXSCMP=NXSCMP+1
                SCAT(IGTO,IGFROM,IXSR)=XSSCMP(NXSCMP)
 250          CONTINUE
 240        CONTINUE
          ENDIF
 220    CONTINUE
      ENDIF
      IF(IPRINT .GE. 100) THEN
*----
*  Print XS
*----
        IXSR=0
        DO ILEG=MINLEG+1,MAXLEG+1
          IXSR=IXSR+1
          IXSTN=MOD(ITYPRO(ILEG)/IODIV,2)
          IF(IXSTN.NE.0) THEN
            WRITE(NAMXS,'(A4,I2.2,A6)') 'SIGS',ILEG-1,NORD
            WRITE(IOUT,6000) NAMXS
            WRITE(IOUT,6010) (XSREC(IG,IXSR),IG=1,NGROUP)
*----
*  SCAT(IGTO,IGFROM) REPRESENT SCATTERING CROSS SECTION
*  FROM GROUP "IGFROM" TO GROUP "IGTO"
*----
            WRITE(NAMXS,'(A4,I2.2,A6)') 'SCAT',ILEG-1,NORD
            WRITE(IOUT,6000) NAMXS
            DO IGFROM=1,NGROUP
              WRITE(IOUT,6001) IGFROM
              WRITE(IOUT,6010) (SCAT(IGTO,IGFROM,IXSR),IGTO=1,NGROUP)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XSSCMP,IJJ,NJJ)
      RETURN
*----
*  ABORT FORMAT
*----
 6000 FORMAT(/' CROSS SECTION TYPE    = ',A12)
 6001 FORMAT(/' SCATTERING FROM GROUP = ',I10)
 6010 FORMAT(1P,5E16.7)
 9000 FORMAT(/' XDRLGS: ****** ABORT ******'/
     >        ' INVALID LENGTH OF RECORD ',A10/
     >        ' STORAGE SPACE          = ',I10/
     >        ' LENGTH OF RECORD LONG  = ',I10/
     >        ' ***************************')
      END
