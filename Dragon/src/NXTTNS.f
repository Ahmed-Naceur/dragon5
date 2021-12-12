*DECK NXTTNS
      SUBROUTINE NXTTNS(IFTRK ,IFTEMP,IPRINT,RENO  ,NFSUR ,NFREG ,
     >                  NDIM  ,MAXSUB, MAXSGL,NTLINE,NBDR ,IFMT  ,
     >                  KEYMRG,DVNOR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To normalize tracking lines and save track volume 
* normalisation factors on tracking data structure.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IFTRK   pointer to the TRACKING file in
*         creation mode.
* IFTEMP  pointer to a temporary TRACKING data structure in
*         update or creation mode.
* IPRINT  print level.
* RENO    track normalisation option. A value RENO=-1 implies
*         a direction dependent normalization of the tracks
*         for the volume while a value RENO=0, implies
*         a global normalisation.
* NFSUR   number of surfaces.
* NFREG   number of regions.
* NDIM    problem dimensions.
* MAXSUB  maximum number of subtracks in a track.
* MAXSGL  maximum number of segments in a track.
* NTLINE  number of track generated.
* NBDR    number of directions for track normalization.
* IFMT    track format: =0 short; =1 long.
* KEYMRG  index array for surface and volume renumbering.
* DVNOR   track volume normalisation factors.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*  \\\\
*  Based on the XELTI2 and XELTI3 routines.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IFTRK,IFTEMP
      INTEGER          IPRINT,RENO,NFSUR,NFREG,NDIM,MAXSUB,
     >                 MAXSGL,NTLINE,NBDR,IFMT
      INTEGER          KEYMRG(-NFSUR:NFREG)
      DOUBLE PRECISION DVNOR(NFREG,NBDR)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTNS')
*----
*  Local variables
*----
      INTEGER          IRLINE,NSUB,NBSEG,ISEG,IREG,JSEG,JREG,II
      DOUBLE PRECISION WEIGHT
      INTEGER          IRA,IADD(4),INREG,JNREG,ITDIR,IND
      LOGICAL          LNEW
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NUMERO,IANGL
      DOUBLE PRECISION , ALLOCATABLE, DIMENSION(:) :: LENGTH
      DOUBLE PRECISION , ALLOCATABLE, DIMENSION(:,:) :: DADD
*----
*  Scratch storage allocation
*   NUMERO  region/surface identification number for segment.
*   LENGTH  segment length.
*----
      ALLOCATE(NUMERO(MAXSGL),LENGTH(MAXSGL),IANGL(MAXSUB),
     > DADD(NDIM,MAXSUB))
*----
*  Processing starts:
*  print routine opening output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      IRLINE=0
 100  CONTINUE
        IF(IFMT .EQ. 1) THEN
          READ(IFTEMP,END=105) NSUB,NBSEG,WEIGHT,
     >                         (IANGL(II),II=1,NSUB),
     >                         (NUMERO(ISEG),ISEG=1,NBSEG),
     >                         (LENGTH(ISEG),ISEG=1,NBSEG),
     >                         (IADD(IRA),IRA=1,4),
     >                         ((DADD(IRA,II),IRA=1,NDIM),II=1,NSUB)
        ELSE
          READ(IFTEMP,END=105) NSUB,NBSEG,WEIGHT,
     >                         (IANGL(II),II=1,NSUB),
     >                         (NUMERO(ISEG),ISEG=1,NBSEG),
     >                         (LENGTH(ISEG),ISEG=1,NBSEG)
        ENDIF
        IRLINE=IRLINE+1
*----
*  Normalize track LENGTH globally
*----
        IF((RENO .EQ. -1) .AND. (NSUB .GT. 1)) THEN
*         Angular-dependent normalization of a cyclic multi-track
          IND=0
          LNEW=.TRUE.
          DO ISEG=1,NBSEG
            IREG=NUMERO(ISEG)
            IF(IREG .GT. NFREG) THEN
              WRITE(IOUT,9001) NAMSBR,(NUMERO(JSEG),JSEG=1,NBSEG)
              CALL XABORT(NAMSBR//
     >        ': Region number larger than maximum permitted')
            ELSE IF(IREG .GT. 0) THEN
              IF(LNEW) THEN
                IND=IND+1
                IF(IND.GT.NSUB) CALL XABORT(NAMSBR//': NSUB overflow')
                LNEW=.FALSE.
              ENDIF
              ITDIR=IANGL(IND)
              LENGTH(ISEG)=LENGTH(ISEG)*DVNOR(IREG,ITDIR+1)
            ELSE
              LNEW=.TRUE.
            ENDIF
          ENDDO
          IF(IND.NE.NSUB) CALL XABORT(NAMSBR//': Algorithm failure')
        ELSE IF(RENO .LE. 0) THEN
          DO ISEG=1,NBSEG
            IREG=NUMERO(ISEG)
            IF(IREG .GT. NFREG) THEN
              WRITE(IOUT,9001) NAMSBR,(NUMERO(JSEG),JSEG=1,NBSEG)
              CALL XABORT(NAMSBR//
     >        ': Region number larger than maximum permitted')
            ELSE IF(IREG .GT. 0) THEN
              IF(RENO .EQ. -1) THEN
                ITDIR=IANGL(1)
                LENGTH(ISEG)=LENGTH(ISEG)*DVNOR(IREG,ITDIR+1)
              ELSE
                LENGTH(ISEG)=LENGTH(ISEG)*DVNOR(IREG,1)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
*----
*  Change region and surface numbering and
*  compress track line for successive segment with same region
*----
        JSEG=1
        JREG=NUMERO(1)
        JNREG=KEYMRG(JREG)
        NUMERO(1)=JNREG
        DO ISEG=2,NBSEG
          IREG=NUMERO(ISEG)
          INREG=KEYMRG(IREG)
          NUMERO(ISEG)=INREG
          IF(INREG .LT. 0 .OR. INREG .NE. JNREG) THEN
            JSEG=JSEG+1
            NUMERO(JSEG)=NUMERO(ISEG)
            LENGTH(JSEG)=LENGTH(ISEG)
            JNREG=INREG
          ELSE
            LENGTH(JSEG)=LENGTH(JSEG)+LENGTH(ISEG)
          ENDIF
        ENDDO
        NBSEG=JSEG
        IF(IFMT .EQ. 1) THEN
          WRITE(IFTRK) NSUB,NBSEG,WEIGHT,
     >                (IANGL(II),II=1,NSUB),
     >                (NUMERO(ISEG),ISEG=1,NBSEG),
     >                (LENGTH(ISEG),ISEG=1,NBSEG),
     >                (IADD(IRA),IRA=1,4),
     >                ((DADD(IRA,II),IRA=1,NDIM),II=1,NSUB)
        ELSE
          WRITE(IFTRK) NSUB,NBSEG,WEIGHT,
     >                (IANGL(II),II=1,NSUB),
     >                (NUMERO(ISEG),ISEG=1,NBSEG),
     >                (LENGTH(ISEG),ISEG=1,NBSEG)
        ENDIF
        GO TO 100
 105  CONTINUE
      IF(IRLINE .NE. NTLINE) THEN
        WRITE(IOUT,9000) NAMSBR,IRLINE,NTLINE
        CALL XABORT(NAMSBR//
     >': Problem with number of lines on tracking file')
      ENDIF
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(DADD,IANGL,LENGTH,NUMERO)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 9000 FORMAT(' ***** Error in ',A6,' *****'/,
     >       '       Number of lines : ',I10,' and ',I10)
 9001 FORMAT(' ***** Error in ',A6,' *****'/,
     >       ' Regions crossed by line segment :'/10I10)
      END
