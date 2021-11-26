MODULE PCREAD
!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran support module for PMAXS reading.
!
!Copyright:
! Copyright (C) 2019 Ecole Polytechnique de Montreal
!
!Author(s): A. Hebert
!
!-----------------------------------------------------------------------
!
    use PCRDATA

    IMPLICIT NONE

    TYPE(PMAXS_WISE_TYPE),target :: PMAXO
    CHARACTER(80), ALLOCATABLE ::PMAXS_F_name(:)
    INTEGER(4),   ALLOCATABLE ::Bran_struct(:)

    INTEGER(4) :: ntbase1,nhinc
    TYPE(BRANCH_WISE_TYPE), DIMENSION(:,:), POINTER :: incbase !(ntbase1,nhinc)
    INTEGER(4), DIMENSION(:), allocatable :: bset

    INTEGER(4) :: n_hist_type, hist_type(Nallvar)

CONTAINS

!---------------------------------------------------------------------
   SUBROUTINE AllocateBranch(Bran)
!---------------------------------------------------------------------
!
    IMPLICIT NONE   
   
    TYPE(BRANCH_WISE_TYPE),target :: Bran
    INTEGER(4) :: k,NBURN
    NBURN=PMAX%Bset(Bran%ibset)%NBURN
    allocate(Bran%XS(NBURN))
    do k=1,NBURN
        XS=>Bran%XS(k)
        call AllocateXSBlock
    enddo
   END SUBROUTINE AllocateBranch

!---------------------------------------------------------------------
   SUBROUTINE ClearBranch(Bran)
!---------------------------------------------------------------------
!
    IMPLICIT NONE   
   
    TYPE(BRANCH_WISE_TYPE),target :: Bran
    INTEGER(4) :: k,NBURN

    NBURN=PMAX%Bset(bran%ibset)%NBURN
    do k=1,NBURN
        XS=>bran%XS(k)
        CALL Clear_XS
    enddo
    deallocate(Bran%XS)
   END SUBROUTINE ClearBranch

!---------------------------------------------------------------------
   SUBROUTINE AllocateTIVB(TIVB)
!---------------------------------------------------------------------
!
    IMPLICIT NONE
   
    TYPE(HIST_TIV_TYPE),target :: TIVB
    INTEGER(4) :: k,NBURN

    NBURN=PMAX%Bset(TIVB%ibset)%NBURN
    allocate(TIVB%TIV(NBURN))
    do k=1,NBURN
        TIV=>TIVB%TIV(k)
        call Allocate_TIV
    enddo
   END SUBROUTINE AllocateTIVB

!---------------------------------------------------------------------
   SUBROUTINE ClearTIVB(TIVB)
!---------------------------------------------------------------------
!
    IMPLICIT NONE
   
    TYPE(HIST_TIV_TYPE),target :: TIVB
    INTEGER(4) :: k,NBURN

    NBURN=PMAX%Bset(TIVB%ibset)%NBURN
    do k=1,NBURN
        TIV=>TIVB%TIV(k)
        call Deallocate_TIV
    enddo
    deallocate(TIVB%TIV)
   END SUBROUTINE ClearTIVB

!---------------------------------------------------------------------
   SUBROUTINE read_TIV(PMAXS_unit)
!---------------------------------------------------------------------
      IMPLICIT NONE
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TIVtemp
      INTEGER(4) :: PMAXS_unit
      INTEGER(4) i,j,k  !,m
      REAL(8) dump
99    format(8E12.5)

!1)   fission spectrum inverse velocity and detector xs
        if(iXSTI .GT. 0)then
          allocate(TIVtemp(NGROUP,4))
          read(PMAXS_unit,99)((TIVtemp(i,j),i=1,NGROUP),j=1,iXSTI)
          do j=1,iXSTI
             k=iTIV(j)
             if(k .GT. 0)then
                do i=1,NGROUP
                   TIV%sig(i,k)=TIVtemp(i,j)
                enddo
             endif
          enddo
          deallocate(TIVtemp)
       endif
!2)   yiled
    if(pyld)then
       if(lyld)then
          read(PMAXS_unit,99)TIV%YLD(:)
       else
          read(PMAXS_unit,99)
       endif
    endif

!cdf
       IF(tcdf)THEN
          READ(PMAXS_unit,99)((dump,i=1,NGROUP),j=1,NCD)
       ENDIF
! gff
       if(tgff.and.NRODS .GT. 0)then
             read(PMAXS_unit,99)((dump,i=1,NGROUP),j=1,NRODS)
       endif

!3)   BETA of Delayed neutron data
    if(pbet)then
        if(lbet)then
           read(PMAXS_unit,99)TIV%kinp(BBET:EBET)
        else
           read(PMAXS_unit,99)(dump,i=1,NDLAY)
        endif
     endif
!4)lambda of Delayed neutron data
    if(pamb)then
        if(lamb)then
           read(PMAXS_unit,99)TIV%kinp(BLAM:ELAM)
        else
           read(PMAXS_unit,99)(dump,i=1,NDLAY)
        endif
     endif
!5)   Decay heat data
    if(pdec)then
        if(ldec)then
              read(PMAXS_unit,99)TIV%kinp(BDHB:EDHB)
              read(PMAXS_unit,99)TIV%kinp(BDHL:EDHL)
        else
          read(PMAXS_unit,99)(dump,i=1,NDCAY)
          read(PMAXS_unit,99)(dump,i=1,NDCAY)
        endif
     endif
     return
   END SUBROUTINE read_TIV

!---------------------------------------------------------------------
   SUBROUTINE read_XS_Block(PMAXS_unit)
!---------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER(4) :: PMAXS_unit
       REAL(8) LPFtemp(8)
       REAL(8) dump
       INTEGER(4) i,j  !,k
       read(PMAXS_unit,99)((XS%sig(i,j),i=1,NGROUP),j=1,4)

       if(pxes)then
          if(pdet)then
             if(lxes)then
                if(ldet)then
                   read(PMAXS_unit,99)((XS%sig(i,j),i=1,NGROUP),j=5,7),XS%det
                else
                   read(PMAXS_unit,99)((XS%sig(i,j),i=1,NGROUP),j=5,7),(dump,i=1,NGROUP)
                endif
             else
                if(ldet)then
                   read(PMAXS_unit,99)((dump,i=1,NGROUP),j=1,3),XS%det
                else
                   read(PMAXS_unit,99)((dump,i=1,NGROUP),j=1,4)
                endif
             endif
          else
             if(lxes)then
                   read(PMAXS_unit,99)((XS%sig(i,j),i=1,NGROUP),j=5,7)
                   XS%sig(1,5) = XS%sig(1,5) * 1E24
                   XS%sig(2,5) = XS%sig(2,5) * 1E24
                   XS%sig(1,6) = XS%sig(1,6) * 1E24
                   XS%sig(2,6) = XS%sig(2,6) * 1E24
             else
                   read(PMAXS_unit,99)((dump,i=1,NGROUP),j=1,3)
             endif
          endif
       else
          if(pdet)then
                if(ldet)then
                   read(PMAXS_unit,99)XS%det
                else
                   read(PMAXS_unit,99)(dump,i=1,NGROUP)
                endif
          endif
       endif
! sct (scattering cross sections)
       read(PMAXS_unit,99)XS%sct
! adf
            IF(padf)THEN
               IF(ladf)THEN
                  READ(PMAXS_unit,99)((XS%adf(i,j),i=1,NGROUP),j=1,NAD)
                  if(NAD .LT. NADF)then
                     if(NAD .EQ. 1)then
                        do i=1,NGROUP
                            XS%adf(i,:)=XS%adf(i,1)
                        enddo
                     elseif(NAD .EQ. 2)then
                        if(NADF .EQ. 3) then
                           call XABORT('read_XS_Block: Error - Please Use Same NADF In All PMAXS Files')
                        elseif(NADF .EQ. 4)then
                            do i=1,NGROUP
                               XS%adf(i,3)=XS%adf(i,2)
                               XS%adf(i,4)=XS%adf(i,1)
                            enddo
                        else
                            do i=1,NGROUP
                               XS%adf(i,3)=XS%adf(i,1)
                               XS%adf(i,4)=XS%adf(i,2)
                               XS%adf(i,5)=XS%adf(i,1)
                               XS%adf(i,6)=XS%adf(i,2)
                            enddo
                        endif
                     else
                            do i=1,NGROUP
                               XS%adf(i,4)=XS%adf(i,1)
                               XS%adf(i,5)=XS%adf(i,2)
                               XS%adf(i,6)=XS%adf(i,3)
                            enddo
                     endif
                  endif
               ELSE
                  READ(PMAXS_unit,99)((dump,i=1,NGROUP),j=1,NAD)
               ENDIF
            ENDIF
! lpf
  99   format(8e12.5)
       if(iLPF .GT. 0)then
          read(PMAXS_unit,99)(LPFtemp(j),j=1,iLPF)
          if(pded.and.lded)XS%LPF(1:4)=LPFtemp(1:4)
          if(pj1f.and.lj1f)XS%LPF(xlpk:xj1c)=LPFtemp(ilpk:ij1c)
       endif
! cdf
            IF(pcdf)THEN
               IF(lcdf)THEN
                  READ(PMAXS_unit,99)((XS%cdf(i,j),i=1,NGROUP),j=1,NCD)
                  if(NCD .LT. NCDF)then
                     if(NCD .EQ. 1)then
                        do i=1,NGROUP
                            XS%cdf(i,:)=XS%cdf(i,1)
                        enddo
                     elseif(NCD .EQ. 2)then
                               if(NCDF .EQ. 5)then
                                  do i=1,NGROUP
                                     XS%cdf(i,5)=XS%cdf(i,2)
                                     XS%cdf(i,4)=XS%cdf(i,2)
                                     XS%cdf(i,3)=XS%cdf(i,1)
                                     XS%cdf(i,2)=XS%cdf(i,1)
                                  enddo
                               elseif(NCDF .EQ. 6)then
                                  do i=1,NGROUP
                                     XS%cdf(i,3)=XS%cdf(i,1)
                                     XS%cdf(i,4)=XS%cdf(i,2)
                                     XS%cdf(i,5)=XS%cdf(i,1)
                                     XS%cdf(i,6)=XS%cdf(i,2)
                                  enddo
                               elseif(NCDF .EQ. 8)then
                                  do i=1,NGROUP
                                     XS%cdf(i,8)=XS%cdf(i,2)
                                     XS%cdf(i,7)=XS%cdf(i,2)
                                     XS%cdf(i,6)=XS%cdf(i,2)
                                     XS%cdf(i,5)=XS%cdf(i,2)
                                     XS%cdf(i,4)=XS%cdf(i,1)
                                     XS%cdf(i,3)=XS%cdf(i,1)
                                     XS%cdf(i,2)=XS%cdf(i,1)
                                  enddo
                               else
                                  call XABORT('read_XS_Block: Error - Please Use Same NCDF In All PMAXS Files')
                               endif
                     elseif(NCD .EQ. 3)then
                               if(NCDF .EQ. 4)then
                                  do i=1,NGROUP
                                     XS%cdf(i,4)=XS%cdf(i,2)
                                  enddo
                               elseif(NCDF .EQ. 5)then
                                  do i=1,NGROUP
                                     XS%cdf(i,4)=XS%adf(i,1)
                                     XS%cdf(i,5)=XS%adf(i,2)
                                  enddo
                               elseif(NCDF .EQ. 6)then
                                  do i=1,NGROUP
                                     XS%cdf(i,4)=XS%cdf(i,1)
                                     XS%cdf(i,5)=XS%cdf(i,2)
                                     XS%cdf(i,6)=XS%cdf(i,3)
                                  enddo
                               elseif(NCDF .EQ. 8)then
                                  do i=1,NGROUP
                                     XS%cdf(i,4)=XS%cdf(i,2)
                                     XS%cdf(i,5)=XS%adf(i,1)
                                     XS%cdf(i,6)=XS%adf(i,2)
                                     XS%cdf(i,7)=XS%adf(i,2)
                                     XS%cdf(i,8)=XS%adf(i,1)
                                  enddo
                               endif
                     elseif(NCD .EQ. 4)then
                              if(NCDF .EQ. 8)then
                                  do i=1,NGROUP
                                     XS%cdf(i,5)=XS%adf(i,1)
                                     XS%cdf(i,6)=XS%adf(i,2)
                                     XS%cdf(i,7)=XS%adf(i,3)
                                     XS%cdf(i,8)=XS%adf(i,4)
                                  enddo
                              else
                                  call XABORT('read_XS_Block: Error - Please Use Same NCDF In All PMAXS Files')
                              endif
                     elseif(NCD .EQ. 5)then
                                  do i=1,NGROUP
                                     XS%cdf(i,8)=XS%cdf(i,4)
                                     XS%cdf(i,7)=XS%cdf(i,5)
                                     XS%cdf(i,6)=XS%cdf(i,5)
                                     XS%cdf(i,5)=XS%cdf(i,4)
                                     XS%cdf(i,4)=XS%cdf(i,2)
                                  enddo
                     else
                                  call XABORT('read_XS_Block: Please use same NCDF in all PMAXS files')
                     endif
                  endif
               ELSE
                  READ(PMAXS_unit,99)((dump,i=1,NGROUP),j=1,NCD)
               ENDIF
            ENDIF
! gff
       if(pgff.and.NRODS .GT. 0)then
          if(lgff)then
             read(PMAXS_unit,99)XS%gff
          else
             read(PMAXS_unit,99)((dump,i=1,NGROUP),j=1,NRODS)
          endif
       endif
       return
   END SUBROUTINE read_XS_Block

!---------------------------------------------------------------------
   SUBROUTINE det_var_position
!---------------------------------------------------------------------
!          determine variable position in input PMAXS file
         
         IMPLICIT NONE
         
         INTEGER(4) i  !, i_pri,i_adf,i_lpf
! XS block
! LPF
         if(pded)then
            i=4
         else
            i=0
         endif
         if(pj1f)then
            ilpk=i+1
            ij1c=i+4
            iLPF=ij1c
         else
            iLPF=i
         endif

! TIV block
         if(pchi)then
            i=1
            iTIV(i)=xchi
         else
            i=0
         endif
         if(pchd)then
            i=i+1
            iTIV(i)=xchd
         endif
         if(pvel)then
            i=i+1
            iTIV(i)=xinv
         endif
         iXSTI=i
   END SUBROUTINE det_var_position

!---------------------------------------------------------------------
   SUBROUTINE set_var_position
!---------------------------------------------------------------------
!          set variable position in memory and output PMAXS
         
         IMPLICIT NONE
         INTEGER(4) :: i
         formng='(1P008E12.5)'
         if(NGROUP .GT. 4)then
            write(formng(4:6),'(I3.3)')NGROUP
         elseif(NGROUP .EQ. 3)then
            formng='(1P006E12.5)'
         endif

         if(NADF .EQ. 0)ladf=.false.
         if(NCDF .EQ. 0)lcdf=.false.
         if(NRODS .EQ. 0)lgff=.false.

!    ded
         if(lded)then
            i=4
         else
            i=0
         endif
         if(lj1f)then
            xlpk=i+1
            xj1i=i+2
            xj1s=i+3
            xj1c=i+4
            NLPF=xj1c
         else
            NLPF=i
         endif

! TIV block
         if(lchi)then
            i=1
             TIVname(i)='Chi'
         else
            i=0
         endif
         if(lchd)then
            i=i+1
            TIVname(i)='Chd'
         endif
         xchd=i
         if(linv)then
            i=i+1
            TIVname(i)='inV'
         endif
         xinv=i

! beta and lambda
         EBET=NDLAY
         BLAM=EBET+1
         ELAM=EBET+NDLAY
! decay heat
         BDHB=ELAM+1
         EDHB=ELAM+NDCAY
         BDHL=EDHB+1
         EDHL=EDHB+NDCAY

! format
         formng='(1P008E12.5)'
         if(NGROUP .GT. 4)then
             write(formng(4:6),'(I3.3)')NGROUP
         elseif(NGROUP .EQ. 3)then
             formng='(1P006E12.5)'
         endif
   END SUBROUTINE set_var_position

!---------------------------------------------------------------------
   SUBROUTINE read_PMAXS_file(iPMAX,kread,PMAXS_unit)
!---------------------------------------------------------------------
      use PCRDATA

      IMPLICIT NONE

      INTEGER(4) :: iPMAX,kread,PMAXS_unit
      INTEGER(4) :: itemp,i_s
      CHARACTER(8) :: tit
      CHARACTER(80) :: oneline

      read(PMAXS_unit,'(A80)',end=101)oneline
      if(oneline(1:8).NE.'GLOBAL_V') call XABORT('dep_read_pmaxs_file: GLOBAL_V expected.')
!1) global variables
      if(oneline(64:64).eq.' ')then
         read(oneline,*)tit,NHST,NGR,NDL,NDC,NAD,NCD,NRODS,NCOLA, &
            padf,pxes,pded,pj1f,pchi,pchd,pvel,pdet,pyld,pcdf,pgff,pbet,pamb,pdec
            derivatives=.true.
            pzdf = .FALSE.
      else if(oneline(66:66).eq.' ')then
         read(oneline,*)tit,NHST,NGR,NDL,NDC,NAD,NCD,NRODS,NCOLA, &
            padf,pxes,pded,pj1f,pchi,pchd,pvel,pdet,pyld,pcdf,pgff,pbet,pamb,pdec,derivatives
            pzdf = .FALSE.
      else
         read(oneline,*)tit,NHST,NGR,NDL,NDC,NAD,NCD,NRODS,NCOLA, &
            padf,pxes,pded,pj1f,pchi,pchd,pvel,pdet,pyld,pcdf,pgff,pbet,pamb,pdec,pzdf,derivatives
      endif
      if(kread.LE.0)THEN
        if(kread .EQ. -1)THEN
          NGROUP=NGR
          NDLAY =NDL
          NDCAY =NDC
          NADF  =NAD
          NCDF  =NCD
          MHST  =NHST
          MRODS =NRODS
          MCOLA =NCOLA
          if(MCOLA .LT. NROWA)MCOLA=NROWA
          MBset=1
          MBRA=1
          MBCR=0
        else
          if(NGROUP.NE.NGR) then
            call XABORT('read_PMAXS_file: Error - NGROUP must be the same in all PMAXS files')
          endif
          if( NDLAY.NE.NDL)THEN
            if(NDLAY .EQ. 0)THEN
              NDLAY=NDL
            ELSEif(NDL .GT. 0 .AND. pbet .AND. pamb)THEN
              call XABORT('read_PMAXS_file: Error - NDLAY must be the same in all PMAXS files')
            ENDIF
          endif
          if( NDCAY.NE.NDC)THEN
            if(NDCAY .EQ. 0)THEN
               NDCAY=NDC
            elseif(NDC .GT. 0 .AND. pdec)THEN
              call XABORT('read_PMAXS_file: Error - NDCAY must be same in all PMAXS files')
            endif
          endif
          if(  NADF .LT. NAD)NADF=NAD
          if(  NCDF .LT. NCD)NCDF=NCD
          if( MHST   .LT. NHST ) MHST  =NHST
          if( MRODS  .LT. NRODS) MRODS =NRODS
          if( MCOLA  .LT. NCOLA) MCOLA =NCOLA
          if( MCOLA  .LT. NROWA) MCOLA =NROWA
        endif
      endif
      call set_var_position

      read(PMAXS_unit,'(A80)') hcomment(1)
      read(PMAXS_unit,'(A80)') hcomment(2)
      read(PMAXS_unit,'(A80)') hcomment(3)
      lxes=.false.
      NXST=4
      if(INDEX(hcomment(3),"xe,sm" ) /= 0) THEN
        lxes=.true.
        NXST=7
      endif
      if(INDEX(hcomment(3),"det" ) /= 0) THEN
        lxes=.true.
        NXST=8
      endif
      tcdf=.false.
      tgff=.false.
      if(pcdf)then
         if(INDEX(hcomment(3),"CDF" ) /= 0) THEN
           tcdf=.true.
           pcdf=.false.
         ENDIF
      endif
      if(pgff)then
         if(INDEX(hcomment(3),"GFF" ) /= 0) THEN
           tgff=.true.
           pgff=.false.
         ENDIF
      endif
      read(PMAXS_unit,'(A80)') hcomment(4)
      read(PMAXS_unit,'(A80)') hcomment(5)
      read(PMAXS_unit,'(A80)') hcomment(6)

      call read_pmax_head(iPMAX, PMAXS_unit)
!4) XS Set identification
      do
         read(PMAXS_unit,*,end=101)tit
         if(tit .EQ. 'XS_SET')exit
      enddo

      backspace(PMAXS_unit)
      read(PMAXS_unit,*)tit,itemp,i_s,itemp,itemp,NCOLA,NROWA,NPART,PITCH,XBE,YBE,iHMD,Dsat,ARWatR,ARByPa,ARConR

      CALL test_pinpower

      call AllocatePMAXS
      call DEP_read_main(PMAXS_unit)
      return
      101 call XABORT('read_PMAXS_file: Error - Reached The End Of PMAXS File')
   END SUBROUTINE read_PMAXS_file

!---------------------------------------------------------------------
   SUBROUTINE test_pinpower 
!---------------------------------------------------------------------
     IMPLICIT NONE

     computed_part: SELECT CASE (NPART)
     CASE (0)
         NCOL=NCOLA
         NROW=NROWA
         NRODS=NCOL*NROW
     CASE (1)
         NCOL=NCOLA
         NROW=NROWA
         if(NCOL.ne.NROW) THEN
            call XABORT('test_pinpower: Error - Assembly Must Be Square For NPART=1')
         END IF
         NRODS=NCOL*(NCOL+1)/2
     CASE (2)
         NCOL=(NCOLA+1)/2
         NROW=(NROWA+1)/2
         NRODS=NCOL*NROW
     CASE (3)
         NCOL=(NCOLA+1)/2
         NROW=(NROWA+1)/2
         if(NCOL.ne.NROW) THEN
            call XABORT('test_pinpower: Error - Assembly Must Be Square For NPART=3')
         END IF
         NRODS=NCOL*(NCOL+1)/2
     END SELECT computed_part
   END SUBROUTINE test_pinpower

!---------------------------------------------------------------------
   SUBROUTINE read_pmax_head(iPMAX, PMAXS_unit)
!---------------------------------------------------------------------
     use PCRDATA
     IMPLICIT NONE
     INTEGER(4) :: iPMAX, PMAXS_unit
     INTEGER(4) :: i,ibra,itemp,inb,j
     CHARACTER(8) :: tit

     if(NDL .EQ. 0)then
       pbet=.false.
       pamb=.false.
     endif
     if(NDC .EQ. 0)pdec=.false.
     if(NAD .EQ. 0)padf=.false.
     if(NCD .EQ. 0)pcdf=.false.
     if(NRODS .EQ. 0)pgff=.false.

     call det_var_position

     bran_i=>Bran_info(iPMAX)
     if(bran_i%NOT_assigned)then
       bran_i%NOT_assigned=.false.
!2) States data
       do
          read(PMAXS_unit,*,end=101)tit
          if(tit .EQ. 'STA_VAR') then
             backspace(PMAXS_unit)
             read(PMAXS_unit,*)tit,Nstat_var
             bran_i%Nstat_var=Nstat_var
             allocate(bran_i%var_ind(Nstat_var),bran_i%var_nam(Nstat_var))
             var_ind=>bran_i%var_ind
             backspace(PMAXS_unit)
             read(PMAXS_unit,*)tit,Nstat_var,bran_i%var_nam(1:Nstat_var)
             ktf=0
             inb=1
             do i=1,Nstat_var
                validname=.false.
                do j=inb,Nallvar
                     if(bran_i%var_nam(i).eq.all_var_nam(j))then
                      validname=.true.
                      var_ind(i)=j
                      inb=j+1
                      exit
                      endif
                enddo
                if(validname)then
                   if(inb .EQ. 5)ktf=i
                else
                   call XABORT('read_pmax_head: Error - State Variable Name Invalid')
                endif
             enddo
             exit
          endif
          if(tit .EQ. 'BRANCHES'.or.tit .EQ. 'BURNUPS'.or.tit .EQ. 'XS_SET') then
             backspace(PMAXS_unit)
             Nstat_var=5
             bran_i%Nstat_var=Nstat_var
             allocate(bran_i%var_ind(Nstat_var),bran_i%var_nam(Nstat_var))
             var_ind=>bran_i%var_ind
             ktf=4
             do i=1,Nstat_var
                var_ind(i)=i
                bran_i%var_nam(i)=all_var_nam(i)
             enddo
             exit
          endif
       enddo

!2) States data
       allocate(bran_i%NBR(Nstat_var))
       NBR=>bran_i%NBR
       NBRA=1
       do
          read(PMAXS_unit,*,end=101)tit
          if(tit .EQ. 'BRANCHES') then
             backspace(PMAXS_unit)
             read(PMAXS_unit,*)tit,itemp,NBR
             do i=1,Nstat_var
                NBRA=NBRA+NBR(i)
             enddo
             allocate(bran_i%state(Nstat_var,NBRA),bran_i%state_nam(NBRA))
             if(NBRA .GT. 1)then
                state=>bran_i%state
                 if(ktf .GT. 0)then
                   do ibra=1,NBRA
                      read(PMAXS_unit,*,end=101)bran_i%state_nam(ibra),itemp,state(:,ibra)
                      state(ktf,ibra)=dsqrt(state(ktf,ibra))
                   enddo
                else
                   do ibra=1,NBRA
                      read(PMAXS_unit,*,end=101)bran_i%state_nam(ibra),itemp,state(:,ibra)
                   enddo
                endif
             else
                bran_i%state=0
             endif
             exit
          endif
          if(tit .EQ. 'BURNUPS'.or.tit .EQ. 'XS_SET')   then
             backspace(PMAXS_unit)
             allocate(bran_i%state(Nstat_var,NBRA))
             bran_i%state=0
             NBR=0
             exit
          endif
       enddo
       bran_i%NBRA=NBRA
       bran_i%ktf=ktf
       if(MBRA .LT. NBRA)MBRA=NBRA
       if(var_ind(1) .EQ. 1)then
          if(MBCR .LT. NBR(1))MBCR=NBR(1)
       endif
     else
       Nstat_var=bran_i%Nstat_var
       NBRA=bran_i%NBRA
       ktf=bran_i%ktf
     endif

!3) Burnup  information
     do
        read(PMAXS_unit,*,end=101)tit
        if(tit .EQ. 'BURNUPS') then
            backspace(PMAXS_unit)
            read(PMAXS_unit,*)tit,PMAX%NBset
            if(MBset .LT. PMAX%NBset)MBset=PMAX%NBset
            allocate(PMAX%Bset(PMAX%NBset))
            do i=1,PMAX%NBset
                read(PMAXS_unit,*)itemp,itemp
                allocate(PMAX%Bset(i)%burns(itemp))
                backspace(PMAXS_unit)
                read(PMAXS_unit,*)itemp,PMAX%Bset(i)%NBURN,PMAX%Bset(i)%burns
            enddo
            exit
        endif
        if(tit .EQ. 'XS_SET')   then
            backspace(PMAXS_unit)
            PMAX%NBset=1
            allocate(PMAX%Bset(PMAX%NBset))
            allocate(PMAX%Bset(1)%burns(1))
            PMAX%Bset(1)%NBURN=1
            PMAX%Bset(1)%burns(1)=0
            exit
        endif
     enddo
     return
     101  call XABORT('read_pmax_head: Error - Reached The End Of PMAXS File')
     STOP
   END SUBROUTINE read_pmax_head

!---------------------------------------------------------------------
   SUBROUTINE DEP_read_main(PMAXS_unit)
!---------------------------------------------------------------------
     use PCRDATA
     IMPLICIT NONE
     INTEGER(4) :: PMAXS_unit
     INTEGER(4) :: i,ihst,ibra,itemp,iBset,NBURN
     CHARACTER(4) :: tit4
     CHARACTER(8) :: tit
!     History case wise data
     do ihst=1,NHST
!6) History case identification
       do
         read(PMAXS_unit,*,end=101)tit
         if(tit .EQ. 'HST_CASE')then
            backspace(PMAXS_unit)
            read(PMAXS_unit,*)tit,PMAX%history(:,ihst)
            PMAX%TIVB(ihst)%ibset=1
            exit
         endif
         if(tit .EQ. 'HISTORYC')then
            backspace(PMAXS_unit)
            read(PMAXS_unit,*)tit,PMAX%TIVB(ihst)%ibset,PMAX%history(:,ihst)
            exit
         endif
       enddo
       if(ktf .GT. 0)PMAX%history(ktf,ihst)=sqrt(PMAX%history(ktf,ihst))
       call AllocateTIVB(PMAX%TIVB(ihst))
       NBURN=PMAX%Bset(PMAX%TIVB(ihst)%ibset)%NBURN
       do i=1,NBURN
         TIV=>PMAX%TIVB(ihst)%TIV(i)
         call read_TIV(PMAXS_unit)
       enddo
!branch wise data
!7) State identification      Always
       do ibra=1,NBRA
         read(PMAXS_unit,'(A4,2I4)')tit4,itemp,iBset
         PMAX%branch(ibra,ihst)%iBset=iBset
         NBURN=PMAX%bset(iBset)%NBURN
         call read_branches(NBURN,PMAXS_unit,PMAX%branch(ibra,ihst))
       enddo !ibra
     enddo !ihst
     return
     101  call XABORT('DEP_read_main: Error - Reached The End Of PMAXS File')
   END SUBROUTINE DEP_read_main

!------------------------------------------------------------------
   SUBROUTINE read_branches(NBURN,PMAXS_unit,bran)
!------------------------------------------------------------------
       TYPE(BRANCH_WISE_TYPE) :: bran
       INTEGER(4) :: NBURN,PMAXS_unit,iburn
       call AllocateBranch(bran)
       do iburn=1,NBURN
          XS=>bran%XS(iburn)
          call read_XS_Block(PMAXS_unit)
       enddo
   END SUBROUTINE read_branches

!---------------------------------------------------------------------
   SUBROUTINE AllocatePMAXS
!---------------------------------------------------------------------
!
      PMAX%NCOL=NCOL
      PMAX%NRODS=NRODS
      PMAX%NHST=NHST
      PMAX%NROW=NROW
      PMAX%NPART=NROW
      PMAX%NROWA=NROWA
      PMAX%NCOLA=NCOLA
      PMAX%iHMD=iHMD
      PMAX%Dsat=Dsat
      PMAX%ARWatR=ARWatR
      PMAX%ARByPa=ARByPa
      PMAX%ARConR=ARConR
      PMAX%PITCH=PITCH
      PMAX%XBE=XBE
      PMAX%YBE=YBE
      PMAX%derivatives=derivatives
      allocate(PMAX%TIVB(NHST))
      allocate(PMAX%branch(NBRA,NHST))
      allocate(PMAX%history(Nstat_var,NHST))
      allocate(PMAX%base(NHST))
      allocate(PMAX%invdiff(NHST))
   END SUBROUTINE AllocatePMAXS

!---------------------------------------------------------------------
   SUBROUTINE Clear_PMAXS_file(iPMAX)
!---------------------------------------------------------------------
      use PCRDATA
      IMPLICIT NONE
      INTEGER(4) :: iPMAX, i, ihst, ibra
      bran_i=>Bran_info(iPMAX)
      if(Nstat_var > 0) then
         deallocate(bran_i%var_ind,bran_i%var_nam,bran_i%NBR)
         if(NBRA.GT.0) deallocate(bran_i%state,bran_i%state_nam)
      endif
      do ihst=1,NHST
        call ClearTIVB(PMAX%TIVB(ihst))
        do ibra=1,NBRA
          call ClearBranch(PMAX%branch(ibra,ihst))
        enddo !ibra
      enddo
      if(PMAX%NBset > 0) then
         do i=1,PMAX%NBset
            deallocate(PMAX%Bset(i)%burns)
         enddo
         deallocate(PMAX%Bset)
      endif
      if(NHST > 0) then
        deallocate(PMAX%TIVB)
        if(NBRA.GT.0) deallocate(PMAX%branch)
        if(Nstat_var > 0) deallocate(PMAX%history)
        deallocate(PMAX%base)
        deallocate(PMAX%invdiff)
      endif
      return
   END SUBROUTINE Clear_PMAXS_file
END MODULE PCREAD
