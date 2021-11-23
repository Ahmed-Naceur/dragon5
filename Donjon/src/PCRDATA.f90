MODULE PCRDATA
!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran module containing PMAXS file information.
!
!Copyright:
! Copyright (C) 2019 Ecole Polytechnique de Montreal
!
!Author(s): A. Hebert
!
!-----------------------------------------------------------------------
!
    IMPLICIT NONE      
    INTEGER(4), PARAMETER  :: Nallvar=12
    REAL(8)             :: state_value(Nallvar)
    
    CHARACTER(2)        :: all_var_nam(Nallvar)
    DATA                   all_var_nam/'CR','DC','PC','TF','TC','IC','DM','PM','TM','IM','DN','BN'/
    INTEGER(4), PARAMETER   :: svCR=1,svDC=2,svPC=3,svTF=4,svTC=5
    INTEGER(4), PARAMETER   :: svIC=6,svDM=7,svPM=8,svTM=9,svIM=10
    INTEGER(4), PARAMETER   :: svDN=11,svBN=12
    
    REAL(8)             :: Sref(Nallvar)
    DATA                Sref/0.0, 0.71, 600.0, 28.3, 580.0, 0.0, 0.71, 600.0, 580.0, 0.0, 0.0, 0.0 /

    LOGICAL :: validname
    LOGICAL :: lHST(Nallvar), lSTT(Nallvar)

    CHARACTER(12) :: formng
    CHARACTER(3) :: TIVname(4)

    INTEGER(4), PARAMETER :: xtr=1,xab=2,xnf=3,xkf=4,xxe=5,xsm=6,xfi=7
    INTEGER(4), PARAMETER :: xdcl=1,xdwr=2,xdbp=3,xdcr=4,xchi=1,BBET=1
    INTEGER(4) :: xchd,xinv, EBET,BLAM,ELAM,BDHB,EDHB,BDHL,EDHL
    INTEGER(4) :: xlpk,xj1i,xj1s,xj1c
    INTEGER(4) :: NXST,NLPF,iLPF,iXSTI
    INTEGER(4) :: iTIV(4),ilpk,ij1c,iread_xs

    INTEGER(4) :: NGR,NDL,NDC,NAD,NCD
    INTEGER(4) :: NGROUP,NDLAY,NDCAY,NADF,NZDF,NCDF,iups,Nset,NTDF
    INTEGER(4) :: NHST,NBRA,NBCR,NBset,NRODS,NCOL,NROW,NPART,NROWA,NCOLA,NXSB
    INTEGER(4) :: MHST,MBRA,MBCR,MBset,MRODS,MCOLA
    INTEGER(4) :: N_Bran_struct,Nstat_var,ktf
    INTEGER(4), DIMENSION(:), POINTER :: var_ind,NBR
    REAL(8) :: iHMD,Dsat,ARWatR,ARByPa,ARConR,PITCH,XBE,YBE,minw,maxw,maxws,minws
    REAL(8), DIMENSION(:,:), POINTER :: state

    LOGICAL :: ladf,lxes,lded,lj1f,lchi,lchd,linv,ldet,lyld,lcdf,lgff,lbet,lamb,ldec,lzdf
    LOGICAL :: tcdf,tgff
    LOGICAL :: padf,pxes,pded,pj1f,pchi,pchd,pvel,pdet,pyld,pcdf,pgff,pbet,pamb,pdec,pzdf
    LOGICAL :: lcrp,lppm,lxesm,derivatives, outrange

    CHARACTER(80),DIMENSION(6) :: hcomment

    TYPE Branches_info
       INTEGER(4) :: Nstat_var, ktf, NBRA
       INTEGER(4), DIMENSION(:),   POINTER   :: var_ind     !(Nstat_var)
       INTEGER(4), DIMENSION(:),   POINTER   :: NBR         !(Nstat_var)
       Character(2), DIMENSION(:), POINTER   :: state_nam   !(NBRA)
       REAL(8), DIMENSION(:,:), POINTER   :: state          !(Nstat_var,NBRA)
       logical :: NOT_assigned
       Character(2), DIMENSION(:), POINTER :: var_nam       !(Nstat_var)
    END TYPE Branches_info
   
    TYPE XSBLOCK_TYPE
         REAL(8), DIMENSION(:,:), POINTER ::sig  !(NGROUP,NXST)
!                                           1   2   3   4   5   6   7
!                                          xtr,xab,xnf,xkf,xfi,xxe,xsm
         REAL(8), DIMENSION(:,:), POINTER ::sct,adf,cdf,gff,zdf
         REAL(8), DIMENSION(:),   POINTER ::LPF,det
         REAL(8) ::    kinf,B2,kinfB,kinfL
! Average assembly flux
         REAL(8), DIMENSION(:),   POINTER :: flux
! Axial surface flux, (g,bottom->top)         
         REAL(8), DIMENSION(:,:), POINTER :: zflx
! Radial surface flux, (g,W-S-E-N) if cart, if hex(g,NW-W-SW-SE-E-NE)
         REAL(8), DIMENSION(:,:), POINTER :: rflx
! Axial current in, (g,bottom->top,i-o-n)         
         REAL(8), DIMENSION(:,:,:), POINTER :: zcur
! Radial surface flux, (g,W-S-E-N,i-o-n) if cart, if hex(g,NW-W-SW-SE-E-NE)
         REAL(8), DIMENSION(:,:,:), POINTER :: rcur  
! Groupwise yields
         REAL(8), DIMENSION(:),   POINTER :: yldI, yldXe, yldPm
! Xe, Sm, I, Pm Number Densities         
         REAL(8) :: NDXE,NDSM,NDI,NDPM
    END TYPE XSBLOCK_TYPE

    TYPE BRANCH_WISE_TYPE
        INTEGER(4)                        :: iBset
        TYPE(XSBLOCK_TYPE),dimension(:),pointer:: XS(:) !(NBURN)
    END TYPE BRANCH_WISE_TYPE

    TYPE PMAXS_WISE_TYPE
         logical   derivatives
         INTEGER(4) :: NCOL,NRODS,NROW,NPART,NROWA,NCOLA
         INTEGER(4) :: NHST,NBset
         REAL(8),  DIMENSION(:,:), POINTER         :: history      !(Nstat_var,NHST)
         REAL(8),  DIMENSION(:), POINTER           :: invdiff      !(NHST)
         INTEGER(4),  DIMENSION(:), POINTER           :: base         !(NHST)
         REAL(8):: iHMD,Dsat,ARWatR,ARByPa,ARConR,PITCH,XBE,YBE

         TYPE(BRANCH_WISE_TYPE), DIMENSION(:,:), POINTER :: branch !(NBRA,NHST)
         TYPE(HIST_TIV_TYPE),     DIMENSION(:), POINTER :: TIVB   !(NHST)
         TYPE(Burnup_info),DIMENSION(:), POINTER  :: Bset !(NBset)
    END TYPE PMAXS_WISE_TYPE

    TYPE HIST_TIV_TYPE
        INTEGER(4)                        :: iBset
        TYPE(TH_INDEP_VAR),dimension(:),pointer:: TIV(:) !(NBURN)
    END TYPE HIST_TIV_TYPE

    TYPE Burnup_info
        INTEGER(4) :: NBURN
        REAL(8),DIMENSION(:), POINTER  :: burns !(NBURN)
    end type Burnup_info

    TYPE TH_INDEP_VAR
       REAL(8), DIMENSION(:,:),POINTER :: sig   !(NGROUP,xinv)
!                                                xchi,xchd,xinv,
       REAL(8), DIMENSION(:),POINTER :: kinp !bet,lam,dhb,dhl
       REAL(8) :: YLD(3) !YID,YXE,YPM
       REAL(8) :: NDXE,NDSM,NDI,NDPM
       REAL(8) :: POWER
       REAL(8) :: DAYS
       REAL(8) :: BURNUP
    END TYPE TH_INDEP_VAR

    TYPE XSBLOCK_ITEM
       INTEGER :: IBURN ! burnup step -- added by EPM
       REAL(8) :: DELTA ! delta local variable -- added by EPM
       TYPE(XSBLOCK_TYPE), POINTER :: XS
       TYPE(TH_INDEP_VAR), POINTER :: TIV
    END TYPE XSBLOCK_ITEM

    TYPE(PMAXS_WISE_TYPE),DIMENSION(:), POINTER :: PMAXS !(XS_F_NUM)
    TYPE(PMAXS_WISE_TYPE),pointer:: PMAX
    TYPE(Branches_info), DIMENSION(:), target,allocatable ::Bran_info
    TYPE(Branches_info), POINTER ::bran_i
    TYPE(XSBLOCK_TYPE), target,allocatable :: XSCR(:)
    TYPE(XSBLOCK_TYPE), target :: XSND
    TYPE(XSBLOCK_TYPE), pointer:: XS
    TYPE(TH_INDEP_VAR), pointer:: TIV

contains

!---------------------------------------------------------------------
   SUBROUTINE AllocateXSBlock
!---------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER(4) :: ireg
      
      IF (NADF .GT. 4) THEN
          ireg = NADF
      ELSE
          ireg = 4
      END IF

      allocate(XS%sig(NGROUP,NXST))
      allocate(XS%sct(NGROUP,NGROUP))
      allocate(XS%flux(NGROUP))
      allocate(XS%rflx(NGROUP,ireg))
      allocate(XS%zflx(NGROUP,2))
      allocate(XS%rcur(NGROUP,ireg,3))
      allocate(XS%zcur(NGROUP,2,3))
      IF (lyld) THEN
          allocate(XS%yldI(NGROUP))
          allocate(XS%yldXe(NGROUP))
          allocate(XS%yldPm(NGROUP))
      ELSE
          allocate(XS%yldI(1))
          allocate(XS%yldXe(1))
          allocate(XS%yldPm(1))
      END IF
      if(ladf)then
         allocate(XS%adf(NGROUP,NADF))
      else
         allocate(XS%adf(1,1))
      endif
      if(lzdf)then
         allocate(XS%ZDF(NGROUP,NZDF))
         NTDF = NADF + NZDF
      else
         allocate(XS%ZDF(1,1))
         NTDF = NADF
      endif
      if(NLPF .GT. 0)then
         allocate(XS%LPF(NLPF))
      else
         allocate(XS%LPF(1))
      endif
      if(ldet)then
         allocate(XS%det(NGROUP))
      else
         allocate(XS%det(1))
      endif
      if(lcdf)then
         allocate(XS%cdf(NGROUP,NCDF))
      else
         allocate(XS%cdf(1,1))
      endif
      if(lgff.and.NRODS .GT. 0)then
         allocate(XS%gff(NGROUP,NRODS))
      else
         allocate(XS%gff(1,1))
      endif
      
      CALL Default_XS
   END SUBROUTINE AllocateXSBlock

!---------------------------------------------------------------------
   SUBROUTINE DeallocateXSBlock
!---------------------------------------------------------------------
      deallocate(XS%sig)
      deallocate(XS%sct)
      deallocate(XS%adf)
      deallocate(XS%zdf)
      deallocate(XS%cdf)      
      deallocate(XS%LPF)
      deallocate(XS%det)
      deallocate(XS%gff)
      deallocate(XS%flux)
      deallocate(XS%rflx)
      deallocate(XS%zflx)
      deallocate(XS%rcur)
      deallocate(XS%zcur)      
   END SUBROUTINE DeallocateXSBlock

!---------------------------------------------------------------------
   SUBROUTINE Clear_XS
!---------------------------------------------------------------------
      XS%sig=0
      XS%sct=0
      XS%adf=0
      XS%zdf=0
      XS%cdf=0
      XS%LPF=0
      XS%det=0
      XS%gff=0
      XS%flux=0
      XS%yldI=0
      XS%yldXe=0
      XS%yldPm=0
      XS%rflx=0
      XS%zflx=0
      XS%rcur=0
      XS%zcur=0           
   END SUBROUTINE Clear_XS

!---------------------------------------------------------------------
   SUBROUTINE Default_XS
!---------------------------------------------------------------------
      CALL Clear_XS

      XS%kinf  = 0.0
      XS%kinfB = 1.0
      XS%kinfL = 1.0
      XS%B2    = 0.0
                 
      XS%ndxe  = 0.0
      XS%ndsm  = 0.0
      XS%ndi   = 0.0
      XS%ndpm  = 0.0      

      XS%adf   = 1.0
      XS%zdf   = 1.0
      XS%cdf   = 1.0      
      XS%gff   = 1.0      

      XS%yldi   = 0.0
      XS%yldxe  = 0.0
      XS%yldpm  = 0.0            
   END SUBROUTINE Default_XS   
   
!---------------------------------------------------------------------
   SUBROUTINE Allocate_TIV
!---------------------------------------------------------------------
     if(xinv .GT. 0)then
       allocate(TIV%sig(NGROUP,xinv))
     else
       allocate(TIV%sig(1,1))
     endif
     if(EDHL .GT. 0)then
       allocate(TIV%kinp(EDHL))
     else
       allocate(TIV%kinp(1))
     endif
     TIV%sig=0
     TIV%kinp=0
     TIV%yld=0
     TIV%power=0.0
     TIV%days=0.0
     TIV%burnup=0.0      
     TIV%ndxe=0.0
     TIV%ndsm=0.0
     TIV%ndi =0.0
     TIV%ndpm=0.0      
   END SUBROUTINE Allocate_TIV

!---------------------------------------------------------------------
   SUBROUTINE Deallocate_TIV
!---------------------------------------------------------------------
     deallocate(TIV%sig)
     deallocate(TIV%kinp)
   END SUBROUTINE Deallocate_TIV
END MODULE PCRDATA
