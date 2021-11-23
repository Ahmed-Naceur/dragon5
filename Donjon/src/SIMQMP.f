*DECK SIMQMP
      SUBROUTINE SIMQMP(LX,LY,LXMIN,LYMIN,HHX,IHY,CYCLE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Unfold the quarter shuffling map to full shuffling map, using
* rotations around the center.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
*
*Author(s): 
* V. Salino
*
*Parameters: input
* LX      number of assemblies along the X axis.
* LY      number of assemblies along the Y axis.
* LXMIN   coordinates on X axis of the first assembly.
* LYMIN   coordinates on Y axis of the first assembly.
* HHX     naval battle indices along X axis.
* IHY     naval battle indices along Y axis.
*
*Parameters: input/output
* CYCLE   shuffling matrix for refuelling given as a quarter map,
*         and returned as a full reconstructed matrix
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LX,LY,LXMIN,LYMIN,IHY(LY)
      CHARACTER HHX(LX)*1,CYCLE(LX,LY)*4
*----
*  LOCAL VARIABLES
*   ROTMAT  counter-clockwise rotation matrices, with an Y-axis directed
*           downward.
*           ROTMAT(x,x,1) <  90 degrees rotation matrix
*           ROTMAT(x,x,2) < 180 degrees rotation matrix
*           ROTMAT(x,x,3) < 270 degrees rotation matrix
*----
      INTEGER INTG2,XPOS,YPOS,Q
      REAL XCENTER,YCENTER,ROTX(3),ROTY(3),IROT(3),JROT(3),
     >     ROTMAT(2,2,3)
      CHARACTER TEXT4*4,TEXT1*1,RECONS(3)*4
*
      DATA ROTMAT(1,1,1), ROTMAT(1,2,1)/+0.0, +1.0/
      DATA ROTMAT(2,1,1), ROTMAT(2,2,1)/-1.0, +0.0/
*
      DATA ROTMAT(1,1,2), ROTMAT(1,2,2)/-1.0, +0.0/
      DATA ROTMAT(2,1,2), ROTMAT(2,2,2)/+0.0, -1.0/
*
      DATA ROTMAT(1,1,3), ROTMAT(1,2,3)/+0.0, -1.0/
      DATA ROTMAT(2,1,3), ROTMAT(2,2,3)/+1.0, +0.0/
*
      IF(LX.NE.LY) CALL XABORT('@SIMQMP: QMAP KEYWORD IS NOT
     > COMPATIBLE WITH A NON-SQUARE REFUELLING SCHEME.')
      XCENTER=(REAL(LX)+1)/2
      YCENTER=(REAL(LY)+1)/2
      DO J=LYMIN,LY
        DO I=LXMIN,LX
*         Excluding potential central assembly from reconstruction
          IF(.NOT.(MOD(LX,2).EQ.1.AND.I.EQ.LXMIN.AND.J.EQ.LYMIN)) THEN
            TEXT4=CYCLE(I,J)
            DO Q=1,3
              IF((TEXT4.NE.'NEW').AND.(TEXT4.NE.'|').AND.(TEXT4.NE.'-')
     >        .AND.(TEXT4.NE.'-|-').AND.(TEXT4.NE.'SPC').AND.
     >        (TEXT4(4:).NE.'@')) THEN
                READ(TEXT4,'(A1,I2)') TEXT1,INTG2
                XPOS=0
                DO K=1,LX
                  IF(HHX(K).EQ.TEXT1) THEN
                    IF(XPOS.NE.0)CALL XABORT('@SIMQMP: X-AXIS HAS '
     >              //'MULTIPLE TIMES THE SAME COORDINATES. CHECK '
     >              //'YOUR RESINI: CALL.')
                    XPOS=K
                  ENDIF
                ENDDO
                IF(XPOS.EQ.0) CALL XABORT('@SIMQMP: UNABLE TO FIND XPO'
     >          //'S(1).')
                YPOS=0
                DO K=1,LY
                  IF(IHY(K).EQ.INTG2) THEN
                    IF(YPOS.NE.0)CALL XABORT('@SIMQMP: Y-AXIS HAS '
     >              //'MULTIPLE TIMES THE SAME COORDINATES. CHECK '
     >              //'YOUR RESINI: CALL.')
                    YPOS=K
                  ENDIF
                ENDDO
                IF(YPOS.EQ.0) CALL XABORT('@SIMQMP: UNABLE TO FIND YPO'
     >          //'S(2).')
*               Reconstruction of an element of the matrix
                ROTX(Q)=ROTMAT(1,1,Q)*(REAL(XPOS)-XCENTER)
     >                 +ROTMAT(1,2,Q)*(REAL(YPOS)-YCENTER)+XCENTER
                ROTY(Q)=ROTMAT(2,1,Q)*(REAL(XPOS)-XCENTER)
     >                 +ROTMAT(2,2,Q)*(REAL(YPOS)-YCENTER)+YCENTER
                WRITE(RECONS(Q),'(A1,I2.2)') HHX(INT(ROTX(Q))),
     >                                       IHY(INT(ROTY(Q)))
              ELSE
                RECONS(Q)=TEXT4
              ENDIF
*             Coordinates of the assembly to be filled with
*             reconstructed information
              IROT(Q)=ROTMAT(1,1,Q)*(REAL(I)-XCENTER)
     >               +ROTMAT(1,2,Q)*(REAL(J)-YCENTER)+XCENTER
              JROT(Q)=ROTMAT(2,1,Q)*(REAL(I)-XCENTER)
     >               +ROTMAT(2,2,Q)*(REAL(J)-YCENTER)+YCENTER
            ENDDO
*         
            IF((J.EQ.LYMIN).AND.(MOD(LX,2).EQ.1)) THEN
              IF(RECONS(3).NE.CYCLE(INT(IROT(3)),INT(JROT(3)))) THEN
                WRITE(6,10)
                WRITE(6,20) HHX(I),IHY(J),CYCLE(I,J),RECONS(3),
     >          HHX(INT(IROT(3))),IHY(INT(JROT(3))),
     >          CYCLE(INT(IROT(3)),INT(JROT(3)))
                CALL XABORT('@SIMQMP: CHECK FOR AN ERROR IN THE QUARTE'
     >          //'R-MAP RELOADING PATTERN OR SWITCH TO MAP KEYWORD.')
              ENDIF
            ENDIF
*         
            DO Q=1,3
              CYCLE(INT(IROT(Q)),INT(JROT(Q)))=RECONS(Q)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      RETURN
*
  10  FORMAT('@SIMQMP: INCONSISTENCY IN REDUNDANT DATA. THE ',
     > 'QUARTER-MAP RELOADING PATTERN IS NOT QUARTER-SYMETRIC.')
  20  FORMAT('CONTENT OF ',A1,I2.2,' (',A4,') IS SUPPOSED TO LEAD TO "'
     > ,A4,'" IN ',A1,I2.2,', BUT "',A4,'" HAS BEEN SPECIFIED ',
     > 'INSTEAD.')
      END
