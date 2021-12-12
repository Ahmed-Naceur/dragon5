*DECK KELRNG
      FUNCTION KELRNG(  IPRT,   NDIM, NEXTGE,   NCPC,  MINDO,  MAXDO,
     >                ICORDO,  NSURO,  NVOLO, IDLGEO,
     >                  MAXC, RMESHO, MATGEO,  VOLSO, INDEXO )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Renumber all zones and surfaces for a block by the coordinate 
* (rect/cyl) values.
*
*Copyright:
* Copyright (C) 1990 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPRT    intermediate printing level for output.      
* NDIM    number of dimensions.                              
* NEXTGE  rectangular(0)/circular(1) boundary.         
* NCPC    number of cylinders in a type + 3.
* MINDO   min index values for all axes (rect/cyl).    
* MAXDO   max index values for all axes (rect/cyl).    
* ICORDO  principal axis directions (X/Y/Z) meshes.    
* NSURO   number of surfaces for a specific geometry.       
* NVOLO   number of zones for a specific geometry.          
* IDLGEO  specific position for a geometry.            
* MAXC    dimension of rmesho.                       
* RMESHO  real mesh values (rect/cyl).                 
* MATGEO  material numbers corresponding to geometries.     
* VOLSO   volumes and surfaces for each geometry.        
*
*Parameters: output
* INDEXO  coordinates for zones & surfaces of a cell.
* KELRNG  number of surfaces and zones renumbered.
*
*-----------------------------------------------------------------------
*
      IMPLICIT        NONE
*
      INTEGER         KELRNG, IPRT, NDIM, NEXTGE, NCPC, MAXC,
     >                NSURO, NVOLO, IDLGEO
      INTEGER         ICUR(4), MINDO(NCPC), MAXDO(NCPC),
     >                ICORDO(NCPC), INDEXO(4,*), MATGEO(*), IXYZ(3),
     >                MINT(3), MAXT(3), INCT(3), JINT(3), JAXT(3)
      REAL            VOLSO(*), RMESHO(MAXC)
      DOUBLE PRECISION  RECT(2,3), RAC(2), CEC(2), RAYMAX, RAYXY
      LOGICAL         LELCRN
*
      INTEGER         NSU, NVO, IDLGE, NCP, I, J, KSUR, ISUX, IVOX,
     >                INOX, ICX, ICY, ICZ, IDX, IDY, IDZ, JX, JY, JZ,
     >                IMAT, IMATX, IMATY, IMATZ, IMATYZ, IMATR, NBEXT,
     >                JCP, IX, IY, JRAY, NO, NESUR, NEVOL
*
      INTEGER         IOUT, IND
      PARAMETER     ( IOUT=6 )
*
      IND(I)= IDLGE + I
*
      NSU   =  NSURO
      NVO   =  NVOLO
      IDLGE = IDLGEO
      NCP   =  NCPC 
*
*     INITIALISATION OF INDEX AND VARIOUS THINGS
      DO 20 I= NSU, NVO
         MATGEO(IND(I))=0
         VOLSO(IND(I))=0.0
         DO 10 J= 1, 4
            INDEXO(J,IND(I))= 0
   10    CONTINUE
   20 CONTINUE
      KSUR= MOD(NDIM+1,3)
      DO 25 I= 1, 3
         IXYZ(I)= ABS(ICORDO(I))
         JINT(I)= MINDO(IXYZ(I))
         JAXT(I)= MAXDO(IXYZ(I))+1
         IF( ICORDO(I).GT.0 )THEN
            IF( I.EQ.3 )THEN
               MINT(I)= MINDO(IXYZ(I))+1-KSUR
               MAXT(I)= MAXDO(IXYZ(I))+KSUR
            ELSE
               MINT(I)= MINDO(IXYZ(I))
               MAXT(I)= MAXDO(IXYZ(I))+1
            ENDIF
            INCT(I)= +1
         ELSE
            IF( I.EQ.3 )THEN
               MINT(I)= MAXDO(IXYZ(I))+KSUR
               MAXT(I)= MINDO(IXYZ(I))+1-KSUR
            ELSE
               MINT(I)= MAXDO(IXYZ(I))+1
               MAXT(I)= MINDO(IXYZ(I))
            ENDIF
            INCT(I)= -1
         ENDIF
   25 CONTINUE
*
      KELRNG= 0
      ISUX= 0
      IVOX= 0
      INOX= 0
*
*     NUMBER ZONES & SURFACES
      IF( NCP.LT.4 )THEN
*        THERE ARE NO CYLINDER AT ALL
         J= 3
         ICUR(4)= 0
         ICZ= 3
      ELSE
         J= 4
         CEC(1)= DBLE(RMESHO(MINDO(J)-2))
         CEC(2)= DBLE(RMESHO(MINDO(J)-1))
         ICZ= ICORDO(J)
      ENDIF
*
*     AXIS ORDER IN TRUE GEOMETRY
      ICX= MOD(ICZ  , 3) + 1
      ICY= MOD(ICZ+1, 3) + 1
*
*     AXIS ORDER FOR NUMBERING PROCESS
      IDX= IXYZ(ICX)
      IDY= IXYZ(ICY)
      IDZ= IXYZ(ICZ)
*
*     LOOP OVER ALL "ICZ,ICY,ICX" ZONES, THEN RADIUS
      DO 260 JZ= MINT(ICZ), MAXT(ICZ), INCT(ICZ)
         ICUR(IDZ)= JZ-1
         IF( JZ.NE.JINT(ICZ).AND.JZ.NE.JAXT(ICZ) )THEN
            IMATZ= 0
         ELSE
            IMATZ=   - 2*ICZ
            IF(    (INCT(ICZ).EQ.+1.AND.JZ.EQ.MINT(ICZ))
     >         .OR.(INCT(ICZ).EQ.-1.AND.JZ.EQ.MAXT(ICZ)) )
     >              IMATZ= IMATZ+1
         ENDIF
         DO 250 JY= MINT(ICY), MAXT(ICY), INCT(ICY)
            RECT(1,IDY)= DBLE(RMESHO(MAX(JINT(ICY)  ,JY-1)))
            RECT(2,IDY)= DBLE(RMESHO(MIN(JAXT(ICY)-1,JY  )))
            ICUR(IDY)= JY-1
            IF( JY.NE.JINT(ICY).AND.JY.NE.JAXT(ICY) )THEN
               IMATY= 0
            ELSE
               IMATY= -2*IDY
               IF(    (INCT(ICY).EQ.+1.AND.JY.EQ.MINT(ICY))
     >            .OR.(INCT(ICY).EQ.-1.AND.JY.EQ.MAXT(ICY)) )
     >                 IMATY= IMATY+1
            ENDIF
*
*           TO EXCLUDE LINES
            IF( IMATY*IMATZ .NE. 0 ) GO TO 250
            IMATYZ= IMATY + IMATZ
            DO 240 JX= MINT(ICX), MAXT(ICX), INCT(ICX)
               RECT(1,IDX)= DBLE(RMESHO(MAX(JINT(ICX)  ,JX-1)))
               RECT(2,IDX)= DBLE(RMESHO(MIN(JAXT(ICX)-1,JX  )))
               ICUR(IDX)= JX-1
               IF( JX.NE.JINT(ICX).AND.JX.NE.JAXT(ICX) )THEN
                  IMATX= 0
               ELSE
                  IMATX= -2*IDX
                  IF(    (INCT(ICX).EQ.+1.AND.JX.EQ.MINT(ICX))
     >               .OR.(INCT(ICX).EQ.-1.AND.JX.EQ.MAXT(ICX)) )
     >                    IMATX= IMATX+1
               ENDIF
*
*              TO EXCLUDE SINGLE POINTS
               IF( IMATYZ*IMATX .NE. 0 ) GO TO 240
               IMAT= IMATYZ + IMATX
               NBEXT=1
               IF( NCP.GT.3 )THEN
                  IMATR= IMAT
                  RAC(1)= 0.0D0
                  DO 230 JRAY= MINDO(J), MAXDO(J)
                     RAC(2)= DBLE(RMESHO(JRAY))
                     ICUR(4)= JRAY-1
                     IF(LELCRN(CEC,RAC,RECT(1,ICX),RECT(1,ICY)))THEN
                        IF( IMAT.EQ.0 )THEN
*                          ZONE NUMBERING
                           IVOX= IVOX + 1
                           INOX= INOX + 1
                           NO=INOX
                           IMATR= IVOX
                        ELSE
*                          SURFACE NUMBERING
                           ISUX= ISUX - 1
                           NO= ISUX
                        ENDIF
*                       IDENTIFY FACE AND CHARGE THE ZONE OR SURFACE NO
                        MATGEO(IND(NO))= IMATR
                        DO 220 JCP= 1, 4
                           INDEXO(JCP,IND(NO))= ICUR(JCP)
  220                   CONTINUE
                     ELSE
                        IF( IMAT.EQ.0 )THEN
*                          ZONE NUMBERING
                           INOX= INOX + 1
*                       IDENTIFY FACE AND CHARGE THE ZONE OR SURFACE NO
                           MATGEO(IND(INOX))= -1
*                        ELSE
*                          ISUX=ISUX-1
                        ENDIF
                     ENDIF
                     RAC(1)= RAC(2)
  230             CONTINUE
                  ICUR(4)= MAXDO(J)
                  RAYMAX= DBLE(RMESHO(MAXDO(J)))
                  NBEXT=0
                  DO 233 IX= 1, 2
                  DO 232 IY= 1, 2
                  RAYXY= (RECT(IX,ICX)-CEC(1))*(RECT(IX,ICX)-CEC(1))
     >                 + (RECT(IY,ICY)-CEC(2))*(RECT(IY,ICY)-CEC(2))
                    IF( RAYXY.GE.RAYMAX ) NBEXT= NBEXT + 1
  232             CONTINUE
  233             CONTINUE
               ENDIF
               IF( NBEXT.EQ.0 )THEN
*
*                 NUMBER 'INSIDE' OF CYLINDER
                  IF( NEXTGE.EQ.0 )THEN
*
*                    CONSIDER ONLY FOR OVERALL CARTESIAN GEOMETRY
*                    SET IMAT TO -1 TO IDENTIFY REGION EXTRACTED
                     IF( IMAT.EQ.0 )THEN
*
*                       ZONE NUMBERING
                        INOX= INOX + 1
                        MATGEO(IND(INOX))= -1
                     ENDIF
                  ENDIF
               ELSE
*
*                 NUMBER 'OUTSIDE' OF CYLINDER
                  IF( IMAT.EQ.0 )THEN
*
*                    ZONE NUMBERING
                     IVOX= IVOX + 1
                     INOX= INOX + 1
                     IMAT= IVOX
                     NO  = INOX
                  ELSE
*
*                    SURFACE NUMBERING
                     ISUX= ISUX - 1
                     NO= ISUX
                  ENDIF
*
*                 IDENTIFY FACE AND CHARGE THE ZONE OR SURFACE NO
                  MATGEO(IND(NO))= IMAT
                  DO 235 JCP= 1, 4
                     INDEXO(JCP,IND(NO))= ICUR(JCP)
  235             CONTINUE
               ENDIF
  240       CONTINUE
  250    CONTINUE
  260 CONTINUE
*
      KELRNG= IVOX - ISUX + 1
*
      IF( IPRT.GT.5 )THEN
         NESUR= 0
         NEVOL= 0
         DO 549 I= NSU, NVO
            IF( I.GT.0.AND.MATGEO(IND(I)).LT.0 ) NEVOL= NEVOL+1
            IF( I.LT.0.AND.MATGEO(IND(I)).EQ.0 ) NESUR= NESUR-1
  549    CONTINUE
         WRITE(IOUT,'(/13H   NUMBERING ,I8,13H VOLUMES AND ,'//
     >           'I8,10H SURFACES.)') NVO-NEVOL,-NSU+NESUR
         WRITE(IOUT,'(17X,7HMINDIM=,10I8)') (MINDO(J),J=1,NCP)
         WRITE(IOUT,'(17X,7HMAXDIM=,10I8)') (MAXDO(J),J=1,NCP)
*
         DO 550 I= NSU-NESUR, NVO
             WRITE(IOUT,'(8H MATGEO(,I8,2H)=,I6,7H INDEX=,4I8)')
     >                    I, MATGEO(IND(I)), (INDEXO(J,IND(I)),J=1,4)
  550    CONTINUE
      ENDIF
*
      RETURN
      END
