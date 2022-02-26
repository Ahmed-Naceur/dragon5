*DECK NXTAGM
      SUBROUTINE NXTAGM(IPRINT,NFSUR ,NFREG ,NEREG ,NESUR ,
     >                  KEYMRG,MATALB,MATRT ,SURVOL,
     >                  KEYFLX,MATCOD,MATRTN,VOLUME)
*
*----------
*
*Purpose:
* To apply general merge vector to geometry
* and to create the L_TRACK data structure MATCOD VOLUME
* and KEYFLX vectors.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* NFSUR   final number of surfaces.
* NFREG   final number of regions.
* KEYMRG  global merging vector.
* MATALB  global mixture/albedo identification vector (including HMIX).
* MATRT   global BC-REFL+TRAN.
* SURVOL  global surface volume vector.
*
*Parameters: output
* NEREG   final number of regions after MERGE.
* NESUR   final number of surfaces after MERGE.
* KEYFLX  final flux index vector after MERGE.
* MATCOD  final mixture vector after MERGE (including HMIX).
* MATRTN  final BC-REFL+TRAN.
* VOLUME  final volume vector after MERGE.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NFSUR,NFREG,NEREG,NESUR
      INTEGER          KEYMRG(-NFSUR:NFREG),MATALB(-NFSUR:NFREG,2),
     >                 MATRT(NFSUR)
      DOUBLE PRECISION SURVOL(-NFSUR:NFREG)
      INTEGER          KEYFLX(NFREG),MATCOD(NFREG,2),
     >                 MATRTN(NFSUR)
      REAL             VOLUME(NFREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTAGM')
*----
*  Local variables
*----
      INTEGER          IREG,JREG,IMIX,ITST,KSUR,LSUR
      INTEGER          MIMIX,MTST
      DOUBLE PRECISION DVR
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,'(A16,2X,I10)') 'Surface Merge   ',NFSUR
        WRITE(IOUT,'(3I10,E20.10)')
     > (JREG,KEYMRG(JREG),MATALB(JREG,1),SURVOL(JREG),JREG=-1,-NFSUR,-1)
        WRITE(IOUT,'(A16,2X,I10)') 'Region Merge    ',NFREG
        WRITE(IOUT,'(3I10,E20.10)')
     > (JREG,KEYMRG(JREG),MATALB(JREG,1),SURVOL(JREG),JREG=1,NFREG)
      ENDIF
*----
*  Determine number of merge regions
*----
      DVR=0.0D0
      IMIX=0
      MIMIX=0
      KSUR=0
      NEREG=0
      DO JREG=1,NFREG
        NEREG=MAX(NEREG,KEYMRG(JREG))
      ENDDO
      NESUR=0
      DO JREG=-1,-NFSUR,-1
        NESUR=MIN(NESUR,KEYMRG(JREG))
      ENDDO
      NESUR=-NESUR
      DO IREG=1,NEREG
        ITST=-1
        MTST=-1
        DO JREG=1,NFREG
          IF(KEYMRG(JREG) .EQ. IREG) THEN
            IF(ITST .EQ. -1) THEN
              IMIX=MATALB(JREG,1)
              DVR=SURVOL(JREG)
              ITST=1
            ELSE
              IF(IMIX .NE. MATALB(JREG,1) ) CALL XABORT(NAMSBR//
     >': Merging region with different mixtures not permitted')
              DVR=DVR+SURVOL(JREG)
            ENDIF
            IF(MTST .EQ. -1) THEN
              MIMIX=MATALB(JREG,2)
              MTST=1
            ELSE
              IF(MIMIX .NE. MATALB(JREG,2) ) CALL XABORT(NAMSBR//
     >': Merging region with different mixtures not permitted')
            ENDIF
          ENDIF
        ENDDO
        IF(ITST .EQ. -1) CALL XABORT(NAMSBR//
     >': One merge region not defined')
        VOLUME(IREG)=REAL(DVR)
        KEYFLX(IREG)=IREG
        MATCOD(IREG,1)=IMIX
        MATCOD(IREG,2)=MIMIX
      ENDDO
*----
*  Compress MATRT to MATRTN
*----
      CALL XDISET(MATRTN,NESUR,0)
      DO IREG=1,NFSUR
        KSUR=-KEYMRG(-IREG)
        LSUR=-KEYMRG(-MATRT(IREG))
        IF(MATRTN(KSUR) .EQ. 0) THEN
          MATRTN(KSUR)=LSUR
        ELSE
          IF(MATRTN(KSUR) .NE. LSUR) CALL XABORT(NAMSBR//
     >': Merging BC-REFL+TRAN with different surface coupling '//
     >'not permitted')
        ENDIF
      ENDDO
      DO IREG=1,NESUR
         IF(MATRTN(KSUR) .EQ. 0) CALL XABORT(NAMSBR//
     >': Some surfaces in BC-REFL+TRAN have no coupling ')
      ENDDO
*----
*  Print output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
