*DECK THMAVG
      SUBROUTINE THMAVG(IPMAP,IMPX,NX,NY,NZ,NCH,TCOMB,TSURF,DCOOL,
     > TCOOL,PCOOL,HCOOL,POW,NSIMS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print averaged thermalhydraulics properties over the core map.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal.
*
*Author(s): 
* M. Cordiez
*
*Parameters: input
* IPMAP   pointer to the fuelmap object.
* IMPX    printing index (=0 for no print).
* NX      number of meshes along X direction.
* NY      number of meshes along Y direction.
* NZ      number of meshes along Z direction (channel direction).
* NCH     number of fuel channels in the axial plane.
* TCOMB   averaged fuel temperature distribution in K.
* TSURF   surface fuel temperature distribution in K.
* DCOOL   coolant density distribution in g/cc.
* TCOOL   coolant temperature distribution in K.
* PCOOL   coolant pressure distribution in Pa.
* HCOOL   coolant enthalpty distribution in J/kg.
* POW     power distribution in W.
* NSIMS   flag greater than zero to activate axial averaging of
*         thermohydraulics information.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXHHX
      PARAMETER(MAXHHX=30)
      TYPE(C_PTR) IPMAP
      INTEGER IMPX,NX,NY,NZ,NCH,NSIMS
      REAL TCOMB(NZ,NX,NY),TSURF(NZ,NX,NY),DCOOL(NZ,NX,NY),
     > TCOOL(NZ,NX,NY),PCOOL(NZ,NX,NY),HCOOL(NZ,NX,NY),POW(NZ,NX,NY)
*----
*  LOCAL VARIABLES
*----
*     Variables for an averaged fuel bundle
      INTEGER NBLEVELCOMB,IHY(MAXHHX)
      REAL TCOMBAVGAVG, TSURFAVGAVG, DCOOLAVGAVG, TCOOLAVGAVG,
     > PCOOLAVGAVG, HCOOLAVGAVG, POWAVGAVG, POWRELAVGAVG
      REAL TCOMBAVG(NZ), TSURFAVG(NZ), DCOOLAVG(NZ), TCOOLAVG(NZ),
     > PCOOLAVG(NZ), HCOOLAVG(NZ), POWERAVG(NZ), POWRELAVG(NZ)
*     --> POWRELAVG : relative power by axial plane
*     Variables for axially averaged to draw a core map
      REAL TCOMBCM(NX,NY),TSURFCM(NX,NY),DCOOLCM(NX,NY),TCOOLCM(NX,NY),
     > PCOOLCM(NX,NY),HCOOLCM(NX,NY),POWERCM(NX,NY),POWRELCM(NX,NY)
      REAL POWAVGCM, POWRELAVGCM
      CHARACTER HHX(MAXHHX)*1,TEXT1*1,TEXT1B*1,TEXT4*4
*----
*  ALLOCATABLE ARRAYS
*----
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: HZONE
*----
*  RECOVER NAVAL BATTLE COORDINATES OF THE MAP
*----
      IF(NSIMS.GT.0) THEN
        LX=NSIMS/100
        LY=MOD(NSIMS,100)
        ALLOCATE(HZONE(NCH))
        CALL LCMGTC(IPMAP,'S-ZONE',4,NCH,HZONE)
        TEXT4=HZONE(1)
        READ(TEXT4,'(A1,I2)') TEXT1,INTG2
        L=0
        DO K=1,NCH
          TEXT4=HZONE(K)
          READ(TEXT4,'(A1,I2)') TEXT1B,INTG2B
          IF(TEXT1B.EQ.TEXT1) THEN
            L=L+1
            IF(L.GT.MAXHHX)CALL XABORT('@THMAVG: MAXHHX OVERFLOW.(1)')
            IF(L.GT.LY)CALL XABORT('@THMAVG: INCOHERENCE IN BASIC '
     >      //'ASSEMBLY LAYOUT GIVEN IN RESINI: (1).')
            IHY(L)=INTG2B
          ENDIF
        ENDDO
        L=L+1
        IF(L.GT.MAXHHX)CALL XABORT('@THMAVG: MAXHHX OVERFLOW.(2)')
        IHY(L)=0
        L=0
        DO K=1,NCH
          TEXT4=HZONE(K)
          READ(TEXT4,'(A1,I2)') TEXT1B,INTG2B
          IF(INTG2B.EQ.IHY((LY+1)/2)) THEN
            L=L+1
            IF(L.GT.MAXHHX)CALL XABORT('@THMAVG: MAXHHX OVERFLOW.(3)')
            IF(L.GT.LX)CALL XABORT('@THMAVG: INCOHERENCE IN BASIC '
     >      //'ASSEMBLY LAYOUT GIVEN IN RESINI: (2).')
            HHX(L)=TEXT1B
          ENDIF
        ENDDO
        DEALLOCATE(HZONE)
      ENDIF
*----
*  VARIABLES INITIALIZATION
*----
*     Variables for an average fuel bundle
      TCOMBAVGAVG  = 0
      TSURFAVGAVG  = 0
      DCOOLAVGAVG  = 0
      TCOOLAVGAVG  = 0
      PCOOLAVGAVG  = 0
      HCOOLAVGAVG  = 0
      POWAVGAVG    = 0
      POWRELAVGAVG = 0
      NBLEVELCOMB = 0
      DO L=1,NZ
        TCOMBAVG(L) = 0.0
        TSURFAVG(L) = 0.0
        DCOOLAVG(L) = 0.0
        TCOOLAVG(L) = 0.0
        PCOOLAVG(L) = 0.0
        HCOOLAVG(L) = 0.0
        POWERAVG(L) = 0.0
        POWRELAVG(L) = 0.0
      ENDDO
*     Variables for an averaged core layer (map of values)
      POWAVGCM     = 0
      POWRELAVGCM  = 0
*----
*  SUM THE VALUES FOR A EVERY FUEL BUNDLE TO AVERAGE THEM
*----
      NBASS=0
      DO 95 I=1,NX
      DO 90 J=1,NY
      TCOMBCM(I,J) = 0.0
      TSURFCM(I,J) = 0.0
      DCOOLCM(I,J) = 0.0
      TCOOLCM(I,J) = 0.0
      PCOOLCM(I,J) = 0.0
      HCOOLCM(I,J) = 0.0
      POWERCM(I,J) = 0.0
      POWRELCM(I,J) = 0.0
      IF(POW((NZ+1)/2,I,J).GT.0.0) THEN
*        We do not average on the reflectors whose values equal 0
         NBASS=NBASS+1
         TCOMBAVG=TCOMBAVG+TCOMB(:,I,J)
         TSURFAVG=TSURFAVG+TSURF(:,I,J)
         DCOOLAVG=DCOOLAVG+DCOOL(:,I,J)
         TCOOLAVG=TCOOLAVG+TCOOL(:,I,J)
         PCOOLAVG=PCOOLAVG+PCOOL(:,I,J)
         HCOOLAVG=HCOOLAVG+HCOOL(:,I,J)
         POWERAVG=POWERAVG+POW(:,I,J)
      ENDIF
   90 CONTINUE
   95 CONTINUE
*----
*  COMPUTE THE AVERAGED VALUES FOR A GENERIC FUEL BUNDLE
*----
      IF(NSIMS.GT.0) THEN
         TCOMBAVG=TCOMBAVG/REAL(NBASS)
         TSURFAVG=TSURFAVG/REAL(NBASS)
         DCOOLAVG=DCOOLAVG/REAL(NBASS)
         TCOOLAVG=TCOOLAVG/REAL(NBASS)
         PCOOLAVG=PCOOLAVG/REAL(NBASS)
         HCOOLAVG=HCOOLAVG/REAL(NBASS)
         POWERAVG=POWERAVG/REAL(NBASS)
*
*        Computation of the relative power by axial plane and
*        computation of the averaged-on-z-axis values of an average
*        fuel bundle
         DO L=1,NZ
           TCOMBAVGAVG=TCOMBAVGAVG+TCOMBAVG(L)
           TSURFAVGAVG=TSURFAVGAVG+TSURFAVG(L)
           DCOOLAVGAVG=DCOOLAVGAVG+DCOOLAVG(L)
           TCOOLAVGAVG=TCOOLAVGAVG+TCOOLAVG(L)
           PCOOLAVGAVG=PCOOLAVGAVG+PCOOLAVG(L)
           HCOOLAVGAVG=HCOOLAVGAVG+HCOOLAVG(L)
           POWAVGAVG=POWAVGAVG+POWERAVG(L)
           IF(POWERAVG(L).NE.0) NBLEVELCOMB=NBLEVELCOMB+1
         ENDDO
         TCOMBAVGAVG=TCOMBAVGAVG/REAL(NBLEVELCOMB)
         TSURFAVGAVG=TSURFAVGAVG/REAL(NBLEVELCOMB)
         DCOOLAVGAVG=DCOOLAVGAVG/REAL(NBLEVELCOMB)
         TCOOLAVGAVG=TCOOLAVGAVG/REAL(NBLEVELCOMB)
         PCOOLAVGAVG=PCOOLAVGAVG/REAL(NBLEVELCOMB)
         HCOOLAVGAVG=HCOOLAVGAVG/REAL(NBLEVELCOMB)
         POWAVGAVG=POWAVGAVG/REAL(NBLEVELCOMB)
         POWRELAVG=POWERAVG/POWAVGAVG
*
*        Computation of the average relative power by axial plane
*        (it must be equal to 1)
         DO L=1,NZ
           POWRELAVGAVG=POWRELAVGAVG+POWRELAVG(L)
         ENDDO
         POWRELAVGAVG=POWRELAVGAVG/REAL(NBLEVELCOMB)
*
*        There is no use in computing them if the user does not want them
         IF(IMPX.GT.2) THEN
           WRITE(6,'(/28H THMAVG: AVERAGE FUEL BUNDLE/1X,27(1H-))')
           WRITE(6,210) ' ___________________________________________',
     >          '_____________________________________________________',
     >          '___________________'
           WRITE(6,210) '|     |    TFUEL    |    TSURF    |    DCOOL ',
     >          '   |    TCOOL    |    PCOOL    |    HCOOL    |    ',
     >          'POWER    |  POW REL  |'
           WRITE(6,230) '| AVG |',TCOMBAVGAVG,' |',TSURFAVGAVG,' |',
     >          DCOOLAVGAVG,' |',TCOOLAVGAVG,' |',PCOOLAVGAVG,' |',
     >          HCOOLAVGAVG,' |',POWAVGAVG,' |',POWRELAVGAVG,' |'
           WRITE(6,210) '|_____|_____________|_____________|__________',
     >          '___|_____________|_____________|_____________|_______',
     >          '______|___________|'
           DO L=NZ,1,-1
             IF(L.EQ.1) THEN
               WRITE(6,230) '| BOT |',TCOMBAVG(L),' |',TSURFAVG(L),
     >            ' |',DCOOLAVG(L),' |',TCOOLAVG(L),' |',PCOOLAVG(L),
     >            ' |',HCOOLAVG(L),' |',POWERAVG(L),' |',
     >            POWRELAVG(L),' |'
             ELSEIF(L.EQ.NZ) THEN
               WRITE(6,230) '| TOP |',TCOMBAVG(L),' |',TSURFAVG(L),
     >            ' |',DCOOLAVG(L),' |',TCOOLAVG(L),' |',PCOOLAVG(L),
     >            ' |',HCOOLAVG(L),' |',POWERAVG(L),' |',
     >            POWRELAVG(L),' |'
             ELSE
               WRITE(6,235) '| ',L,' |',TCOMBAVG(L),' |',TSURFAVG(L),
     >            ' |',DCOOLAVG(L),' |',TCOOLAVG(L),' |',PCOOLAVG(L),
     >            ' |',HCOOLAVG(L),' |',POWERAVG(L),' |',
     >            POWRELAVG(L),' |'
             ENDIF
           ENDDO
           WRITE(6,210) '|_____|_____________|_____________|_________',
     >          '____|_____________|_____________|_____________|______',
     >          '_______|___________|'
         ENDIF
*----
*  COMPUTE THE AVERAGED VALUES ON THE CORE MAP
*----
*       We do not average on the reflectors whose values equal 0
        DO K=1,NZ
            TCOMBCM(:,:)=TCOMBCM(:,:)+TCOMB(K,:,:)
            TSURFCM(:,:)=TSURFCM(:,:)+TSURF(K,:,:)
            DCOOLCM(:,:)=DCOOLCM(:,:)+DCOOL(K,:,:)
            TCOOLCM(:,:)=TCOOLCM(:,:)+TCOOL(K,:,:)
            PCOOLCM(:,:)=PCOOLCM(:,:)+PCOOL(K,:,:)
            HCOOLCM(:,:)=HCOOLCM(:,:)+HCOOL(K,:,:)
            POWERCM(:,:)=POWERCM(:,:)+POW(K,:,:)
        ENDDO
        TCOMBCM=TCOMBCM/REAL(NBLEVELCOMB)
        TSURFCM=TSURFCM/REAL(NBLEVELCOMB)
        DCOOLCM=DCOOLCM/REAL(NBLEVELCOMB)
        TCOOLCM=TCOOLCM/REAL(NBLEVELCOMB)
        PCOOLCM=PCOOLCM/REAL(NBLEVELCOMB)
        HCOOLCM=HCOOLCM/REAL(NBLEVELCOMB)
        POWERCM=POWERCM/REAL(NBLEVELCOMB)
*       Calculation of the relative power distribution (avg = 1)
        DO 106 I=1,NX
        DO 105 J=1,NY
          POWAVGCM=POWAVGCM+POWERCM(I,J)
  105   CONTINUE
  106   CONTINUE
        POWAVGCM=POWAVGCM/REAL(NBASS)
        POWRELCM=POWERCM/POWAVGCM
        DO 108 I=1,NX
        DO 107 J=1,NY
          POWRELAVGCM=POWRELAVGCM+POWRELCM(I,J)
  107   CONTINUE
  108   CONTINUE
        POWRELAVGCM=POWRELAVGCM/REAL(NBASS)
*
*        There is no use in computing them if the user does not want them
         IF(IMPX.GT.2) THEN
           IDEB=1
           JDEB=1
*          We do not draw the reflector
           I=1
           DO WHILE (POW((NZ+1)/2,I,(NY+1)/2).EQ.0)
             IDEB=IDEB+1
             I=I+1
           END DO
           J=1
           DO WHILE (POW((NZ+1)/2,(NX+1)/2,J).EQ.0)
             JDEB=JDEB+1
             J=J+1
           END DO
*          ********************
*          We write the results
*          ********************
*          Fuel temperature
           WRITE(6,'(/25H THMAVG: AVERAGE CORE MAP/1X,24(1H-))')
           WRITE(6,202) 'TCOMB'
           WRITE(6,203,ADVANCE='NO') (HHX(I),I=1,NX-2*IDEB+2)
           DO 111 J=JDEB,NY
             WRITE(6,204,ADVANCE='NO') IHY(J-JDEB+1)
             ENDLINE=0
           DO 110 I=IDEB,NX
             IF(POW((NZ+1)/2,I,J).GT.0) THEN
               WRITE(6,205,ADVANCE='NO') TCOMBCM(I,J)
               ENDLINE=1
             ELSE IF(ENDLINE.EQ.0) THEN
               WRITE(6,208,ADVANCE='NO')
             ENDIF
  110      CONTINUE
  111      CONTINUE
*          Surface temperature
           WRITE(6,202) 'TSURF'
           WRITE(6,203,ADVANCE='NO') (HHX(I),I=1,NX-2*IDEB+2)
           DO 113 J=JDEB,NY
             WRITE(6,204,ADVANCE='NO') IHY(J-JDEB+1)
             ENDLINE=0
           DO 112 I=IDEB,NX
             IF(POW((NZ+1)/2,I,J).GT.0) THEN
               WRITE(6,205,ADVANCE='NO') TSURFCM(I,J)
               ENDLINE=1
             ELSE IF(ENDLINE.EQ.0) THEN
               WRITE(6,208,ADVANCE='NO')
             ENDIF
  112      CONTINUE
  113      CONTINUE
*          Coolant density
           WRITE(6,202) 'DCOOL'
           WRITE(6,203,ADVANCE='NO') (HHX(I),I=1,NX-2*IDEB+2)
           DO 115 J=JDEB,NY
             WRITE(6,204,ADVANCE='NO') IHY(J-JDEB+1)
             ENDLINE=0
           DO 114 I=IDEB,NX
             IF(POW((NZ+1)/2,I,J).GT.0) THEN
               WRITE(6,206,ADVANCE='NO') DCOOLCM(I,J)
               ENDLINE=1
             ELSE IF(ENDLINE.EQ.0) THEN
               WRITE(6,208,ADVANCE='NO')
             ENDIF
  114      CONTINUE
  115      CONTINUE
*          Coolant temperature
           WRITE(6,202) 'TCOOL'
           WRITE(6,203,ADVANCE='NO') (HHX(I),I=1,NX-2*IDEB+2)
           DO 117 J=JDEB,NY
             WRITE(6,204,ADVANCE='NO') IHY(J-JDEB+1)
             ENDLINE=0
           DO 116 I=IDEB,NX
             IF(POW((NZ+1)/2,I,J).GT.0) THEN
               WRITE(6,205,ADVANCE='NO') TCOOLCM(I,J)
               ENDLINE=1
             ELSE IF(ENDLINE.EQ.0) THEN
               WRITE(6,208,ADVANCE='NO')
             ENDIF
  116      CONTINUE
  117      CONTINUE
*          Coolant pressure
           WRITE(6,202) 'PCOOL'
           WRITE(6,203,ADVANCE='NO') (HHX(I),I=1,NX-2*IDEB+2)
           DO 119 J=JDEB,NY
             WRITE(6,204,ADVANCE='NO') IHY(J-JDEB+1)
             ENDLINE=0
           DO 118 I=IDEB,NX
             IF(POW((NZ+1)/2,I,J).GT.0) THEN
               WRITE(6,207,ADVANCE='NO') PCOOLCM(I,J)
               ENDLINE=1
             ELSE IF(ENDLINE.EQ.0) THEN
               WRITE(6,208,ADVANCE='NO')
             ENDIF
  118      CONTINUE
  119      CONTINUE
*          Coolant enthalpy
           WRITE(6,202) 'HCOOL'
           WRITE(6,203,ADVANCE='NO') (HHX(I),I=1,NX-2*IDEB+2)
           DO 121 J=JDEB,NY
             WRITE(6,204,ADVANCE='NO') IHY(J-JDEB+1)
             ENDLINE=0
           DO 120 I=IDEB,NX
             IF(POW((NZ+1)/2,I,J).GT.0) THEN
               WRITE(6,207,ADVANCE='NO') HCOOLCM(I,J)
               ENDLINE=1
             ELSE IF(ENDLINE.EQ.0) THEN
               WRITE(6,208,ADVANCE='NO')
             ENDIF
  120      CONTINUE
  121      CONTINUE
*          Power
           WRITE(6,202) 'POWER'
           WRITE(6,203,ADVANCE='NO') (HHX(I),I=1,NX-2*IDEB+2)
           DO 123 J=JDEB,NY
             WRITE(6,204,ADVANCE='NO') IHY(J-JDEB+1)
             ENDLINE=0
           DO 122 I=IDEB,NX
             IF(POW((NZ+1)/2,I,J).GT.0) THEN
               WRITE(6,207,ADVANCE='NO') POWERCM(I,J)
               ENDLINE=1
             ELSE IF(ENDLINE.EQ.0) THEN
               WRITE(6,208,ADVANCE='NO')
             ENDIF
  122      CONTINUE
  123      CONTINUE
*          Power
           WRITE(6,209) 'RELATIVE POWER, REFLECTORS EXCLUDED (AVG:',
     1                  POWRELAVGCM,')'
           WRITE(6,203,ADVANCE='NO') (HHX(I),I=1,NX-2*IDEB+2)
           DO 125 J=JDEB,NY
             WRITE(6,204,ADVANCE='NO') IHY(J-JDEB+1)
             ENDLINE=0
           DO 124 I=IDEB,NX
             IF(POW((NZ+1)/2,I,J).GT.0) THEN
               WRITE(6,206,ADVANCE='NO') POWRELCM(I,J)
               ENDLINE=1
             ELSE IF(ENDLINE.EQ.0) THEN
               WRITE(6,208,ADVANCE='NO')
             ENDIF
  124      CONTINUE
  125      CONTINUE
         ENDIF
      ENDIF
      RETURN
*
  202 FORMAT(/1X,A)
  203 FORMAT(1X,20(8X,1A1))
  204 FORMAT(/1X,I2)
  205 FORMAT(F9.1)
  206 FORMAT(F9.3)
  207 FORMAT(1P,E9.2)
  208 FORMAT(9X)
  209 FORMAT(/1X,A,F7.4,A)
  210 FORMAT(1X,A,A,A,A)
  230 FORMAT(1X,A,F12.2,A,F12.2,A,F12.4,A,F12.2,A,3P,E12.4,
     >       A,1P,E12.4,A,1P,E12.4,A,0P,F10.4,A)
  235 FORMAT(1X,A,I3,A,F12.2,A,F12.2,A,F12.4,A,F12.2,A,3P,E12.4,
     >       A,1P,E12.4,A,1P,E12.4,A,0P,F10.4,A)
      END
