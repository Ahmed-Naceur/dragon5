*DECK KELSYM
      FUNCTION KELSYM(   IPRT,   NDIM,   MAXDO,  NSURO,  NVOLO,
     >                 IDLGEO,  INDEXO, MATGEO, KEYSYM )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Generate the vector KEYSYM for a block. 
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
* IPRT    intermediate printing level.                
* NDIM    number of dimensions (2 or 3).                    
* MAXDO   max index values for all axes (rect/cyl).   
* NSURO   number of surfaces for a specific geometry.      
* NVOLO   number of zones for a specific geometry.         
* IDLGEO  specific position for a geometry.           
* INDEXO  coordinates for zones & surfaces of a cell. 
* MATGEO  material numbers corresponding to geometries.    
*
*Parameters: output
* KEYSYM  symmetry numbers corresponding to geometries.
* KELSYM  number of surfaces and zones renumbered.
*
*-----------------------------------------------------------------------
*
      IMPLICIT        NONE
*
      INTEGER         KELSYM, IPRT, NDIM, NSURO, NVOLO, IDLGEO
      INTEGER         MAXDO(*),INDEXO(4,*),KEYSYM(*),MATGEO(*)
*
      INTEGER         ICUR(4), I, J, IVS, MAXPRC, MAXSUI, ISYM, IND
      LOGICAL         SWITCH
      INTEGER         IOUT
      PARAMETER     ( IOUT=6 )
*
      IND(I)= IDLGEO + I
*
      DO 5 IVS= 0, NVOLO
         KEYSYM(IND(IVS))= 0
    5 CONTINUE
      KELSYM= 0
*
*     LOCATES THE SYMMETRIC SURFACE TO EACH SURFACE
      DO 50 IVS = NSURO, -1
         IF( MATGEO(IND(IVS)).EQ.0 )GO TO 51
         MAXPRC= 0
         DO 10 J = 1, 4
            ICUR(J)= INDEXO(J,IND(IVS))
*
*           FIND THE SYMMETRIC SURFACE BY CHANGING END-FACE
            IF( J.LE.NDIM )THEN
               MAXSUI= MAXDO(J)
               IF( ICUR(J).EQ.MAXPRC)THEN
                  ICUR(J)= MAXSUI
               ELSEIF( ICUR(J).EQ.MAXSUI)THEN
                  ICUR(J)= MAXPRC
               ENDIF
               MAXPRC= MAXSUI
            ENDIF
*
*           THE SENTINEL VALUE IS IVS=0
            INDEXO(J,IND(0))= ICUR(J)
   10    CONTINUE
         ISYM= NSURO
   20       SWITCH= .TRUE.
            DO 30 J    = 1, 4
               SWITCH= SWITCH .AND. ICUR(J).EQ.INDEXO(J,IND(ISYM))
   30       CONTINUE
            IF( SWITCH )GO TO 40
            ISYM= ISYM + 1
         GO TO 20
   40    KEYSYM(IND(IVS))= ISYM
         IF( IPRT.GE.10 )THEN
            WRITE(IOUT,'(22H SURFACE SYMMETRIC TO ,I6,4H IS ,I6)')
     >                                         -IVS,     -ISYM
         ENDIF
         IF( ISYM.NE.0 ) KELSYM=KELSYM-1
   51    CONTINUE
   50 CONTINUE
*
*     RESET SENTINEL INDEXO(J,IND(0)) FOR SUBSEQUENT USES
      DO 60 J= 1, 4
         INDEXO(J,IND(0))= 0
   60 CONTINUE
*
      RETURN
      END
