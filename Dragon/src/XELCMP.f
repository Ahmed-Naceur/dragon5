*DECK XELCMP
      SUBROUTINE XELCMP(     NS,     NV,  VOLIN,  MATIN,  MRGIN,
     >                    NSOUT,  NVOUT, VOLOUT, MATOUT,  ITGEO,  ICODE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Merge volumes and surfaces and recompute the number of surfaces 
* and volumes.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NS      number of surfaces before merging.      
* NV      number of zones before merging.      
* VOLIN   volumes and surfaces before merging.      
* MATIN   numbering of sufaces and zones before merging.      
* MRGIN   merging index.                               
* ITGEO   kind of geometry(0,1,2,3).                   
* ICODE   index of boundary conditions.
*
*Parameters: output
* NSOUT   number of surfaces after merging.
* NVOUT   number of zones after merging.
* VOLOUT  volumes and surfaces after merging.
* MATOUT  numbering of sufaces and zones after merging.
*
*-----------------------------------------------------------------------
*
      IMPLICIT      NONE
*
      INTEGER       NS,NV,NSOUT,NVOUT,ITGEO,IVS,IMR,ICNT,I0,IOUT,
     >              IR,JR,MATMRG,LESOIR,
     >              MATIN(-NS:NV),MRGIN(-NS:NV),MATOUT(*),ICODE(6)
      REAL          VOLIN(-NS:NV),VOLOUT(*),ZERO
      CHARACTER*4   CORIEN(0:3,-6:0)
      PARAMETER    ( ZERO= 0.0, IOUT=6 )
      DATA         ((CORIEN(JR,IR),IR=-6,0),JR=0,3)
     >       / ' O6 ',' O5 ',' O4 ',' O3 ',' O2 ',' O1 ','    ',
     >         ' Z+ ',' Z- ','****','****',' R+ ','****','    ',
     >         ' Z+ ',' Z- ','****','****','****','HBC ','    ',
     >         ' Z+ ',' Z- ',' Y+ ',' Y- ',' X+ ',' X- ','    '/
*
*     FIND NSOUT AND NVOUT & INITIALIZE VOLOUT AND MATOUT
      NSOUT= 0
      NVOUT= 0
      DO 10 IVS= -NS, NV
         VOLOUT(IVS+NS+1)= ZERO
         MATOUT(IVS+NS+1)= 0
         IF( IVS.GT.0 )THEN
            IF( MRGIN(IVS).LT.0 )THEN
               CALL XABORT( 'XELCMP: 1.INCOMPATIBLE MERGE INDEX' )
            ENDIF
         ELSEIF( IVS.LT.0 )THEN
            IF( MRGIN(IVS).GT.0 )THEN
               CALL XABORT( 'XELCMP: 2.INCOMPATIBLE MERGE INDEX' )
            ENDIF
         ELSE
            IF( MRGIN(IVS).NE.0 )THEN
               WRITE(IOUT,*) 'XELCMP: *KEYMRG* VECTOR IS:', MRGIN
               CALL XABORT( 'XELCMP: 3.INCOMPATIBLE MERGE INDEX' )
            ENDIF
            IF( VOLIN(IVS).NE.0.0 )THEN
               WRITE(IOUT,*) 'XELCMP: *VOLSUR* VECTOR IS:', VOLIN
               CALL XABORT( 'XELCMP: 4. VOLSUR(0).NE.0 ON TRACK-FILE' )
            ENDIF
            IF( MATIN(IVS).NE.0 )THEN
               WRITE(IOUT,*) 'XELCMP: *MATALB* VECTOR IS:', MATIN
               CALL XABORT( 'XELCMP: 5. MATALB(0).NE.0 ON TRACK-FILE' )
            ENDIF
         ENDIF
         NSOUT= MIN(NSOUT,MRGIN(IVS))
         NVOUT= MAX(NVOUT,MRGIN(IVS))
   10 CONTINUE
      NSOUT= -NSOUT
*
*     ALL VALUES MUST BE PRESENT BETWEEN -NSOUT AND NVOUT IN MRGIN(*)
*     BUT WITH THE SAME MATIN(*) NUMBER FOR MERGED ZONES.
*     NEW(97/11): 0 MEANS REGION IS REMOVED
      DO 30 IMR= -NSOUT, NVOUT
         ICNT= 0
         DO 20 IVS= -NS, NV
            IF( ICNT.EQ.0 ) MATMRG= MATIN(IVS)
            IF( MRGIN(IVS).EQ.IMR )THEN
               ICNT= ICNT+1
               IF( MATMRG.NE.MATIN(IVS) )THEN
                  LESOIR= MATIN(IVS)
                  IF( IVS.GE.0 )THEN
*
*                    FOR MERGING ZONES, ABORT IF NOT SAME *MATALB*
                     WRITE(IOUT,*) '*** ABORT *** ATTEMPT TO MERGE '//
     >                          'MIX ',MATMRG,' WITH MIX ',
     >                                 LESOIR,' IN ZONE #',IVS
                     CALL XABORT( 'XELCMP: 6.INCOMPATIBLE MERGE INDEX' )
                  ELSE
*
*                    FOR MERGING FACES, ABORT IF NOT SAME *ICODE*
                     IF( ICODE(-MATMRG).NE.ICODE(-LESOIR) )THEN
                        WRITE(IOUT,*) '*** ABORT *** ATTEMPT TO MERGE ',
     >                                ' FACE ',-IVS,
     >                             '( ',CORIEN(ITGEO,MATMRG),',ICODE=',
     >                                  ICODE(-MATMRG),') WITH A FACE ',
     >                             '( ',CORIEN(ITGEO,LESOIR),',ICODE=',
     >                                  ICODE(-LESOIR),'). '
                     CALL XABORT( 'XELCMP: 7.INCOMPATIBLE MERGE INDEX' )
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
   20    CONTINUE
         IF( ICNT.EQ.0 )THEN
            CALL XABORT( 'XELCMP: 8.MISSING VALUES IN THE MERGE INDEX' )
         ENDIF
   30 CONTINUE
*
*     COMPUTE VOLOUT AND MATOUT VALUES
      I0= 1 + NSOUT
      DO 40 IVS= -NS, NV
         VOLOUT(I0+MRGIN(IVS))= VOLOUT(I0+MRGIN(IVS))+VOLIN(IVS)
         MATOUT(I0+MRGIN(IVS))= MATIN(IVS)
   40 CONTINUE
*
      RETURN
      END
