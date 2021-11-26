!
!-----------------------------------------------------------------------
!
!Purpose:
! Generate data relative to gigognes originating from nodes and generate
! tracking indices assigned to them.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s):
! G. Civario (CS-SI)
!
!Comments:
! Ce module possede trois fonctions:
!  - generateTrack : fonction d'entree du module
!  - calculFinalMerge : calcul de la numerotation type dragon
!  - ltVec : fonction d'ordre specifique pour des vecteurs d'entiers de tailles
!            differentes
!
module track

    use cellulePlaced
    use segArc
    use ptNodes
    implicit none

contains
  subroutine generateTrack(szP,szSa,nbNode,lgMaxGig,gig,merg)
    integer,intent(in) :: szP,szSa,nbNode,lgMaxGig
    integer,dimension(lgMaxGig*nbNode),intent(out) :: gig
    integer,dimension(nbNode),intent(inout) :: merg
    
    integer :: i,lgMaxMrg,s,d,lgMaxGigTest
    integer,dimension(:,:),allocatable :: mrgMat
    
    lgMaxGigTest = 0
    lgMaxMrg = 0
    do i = 1,szP
       lgMaxGigTest = max(lgMaxGigTest,size(tabCellulePlaced(i)%gig))
       lgMaxMrg = max(lgMaxMrg,size(tabCellulePlaced(i)%mrg))
    end do
    if(lgMaxGigTest /= lgMaxGig) call XABORT('g2s_generatingTrack: lgMax error')
    lgMaxMrg = lgMaxMrg + 1
    allocate(mrgMat(lgMaxMrg,nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: generateTrack => allocation pb")
    gig = 0
    mrgMat = 0
    do i = 1,szSa
       if (tabSegArc(i)%nodeg>0.and.tabSegArc(i)%indCellPg>0) then
          s = size(tabCellulePlaced(tabSegArc(i)%indCellPg)%gig)
          d = (tabSegArc(i)%nodeg-1)*lgMaxGig+1
          gig(d:d+s-1) = tabCellulePlaced(tabSegArc(i)%indCellPg)%gig(1:s)
          s = size(tabCellulePlaced(tabSegArc(i)%indCellPg)%mrg)
          d = tabSegArc(i)%nodeg
          mrgMat(1:s,d) = tabCellulePlaced(tabSegArc(i)%indCellPg)%mrg(1:s)
          mrgMat(s+1,d) = merg(d)
       end if
       if (tabSegArc(i)%noded>0.and.tabSegArc(i)%indCellPd>0) then
          s = size(tabCellulePlaced(tabSegArc(i)%indCellPd)%gig)
          d = (tabSegArc(i)%noded-1)*lgMaxGig+1
          gig(d:d+s-1) = tabCellulePlaced(tabSegArc(i)%indCellPd)%gig(1:s)
          s = size(tabCellulePlaced(tabSegArc(i)%indCellPd)%mrg)
          d = tabSegArc(i)%noded
          mrgMat(1:s,d) = tabCellulePlaced(tabSegArc(i)%indCellPd)%mrg(1:s)
          mrgMat(s+1,d) = merg(d)
       end if
    end do
    call calculFinalMerge(mrgMat,merg)
    deallocate(mrgMat)
  end subroutine generateTrack

  subroutine calculFinalMerge(inMat,outVec)
    integer,dimension(:,:),intent(in) :: inMat
    integer,dimension(:),intent(out)  :: outVec
    
    integer :: i,j,d1,d2,maxD2
    logical :: found,sorted
    integer,dimension(:),allocatable   :: tmpVec
    integer,dimension(:,:),allocatable :: workMat
    
    d1 = size(inMat,1) !profondeur max de gigogne + 1
    d2 = size(inMat,2) !nombre de nodes
    maxD2 = 0
    allocate(workMat(d1,d2))
    allocate(tmpVec(d1),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: calculFinalMerge => allocation pb")
    workMat(:d1,:d2) = 0
    !remplissage de workMat a l'aide d'occurences uniques de lignes de inMat
    do i = 1,d2
       found = .false.
       do j = 1,maxD2
          if (all(workMat(:d1,j)==inMat(:d1,i))) then
             found = .true.
             exit
          end if
       end do
       if (.not. found) then
          maxD2 = maxD2 + 1
          workMat(:d1,maxD2) = inMat(:d1,i)
       end if
    end do
    !classement des lignes de workMat en ordre croissant par bubble-sort
    do j = maxD2,2,-1
       sorted = .true.
       do i = 1,j-1
          if (ltVec(workMat(:d1,i+1),workMat(:d1,i))) then
             tmpVec(:d1)      = workMat(:d1,i+1)
             workMat(:d1,i+1) = workMat(:d1,i)
             workMat(:d1,i)   = tmpVec(:d1)
             sorted = .false.
          end if
       end do
       if (sorted) exit
    end do
    !remplissage de outVec en fonction de l'egalite entre les lignes de workMat
    !et inMat, apres le classement
    do i = 1,d2
       do j = 1,maxD2
          if (all(workMat(:d1,j)==inMat(:d1,i))) then
             outVec(i) = j
             exit
          end if
       end do
    end do
    
    deallocate(workMat,tmpVec)
  end subroutine calculFinalMerge
  
  function ltVec(v1,v2)
    integer,dimension(:),intent(in) :: v1,v2
    logical :: ltVec
    integer :: i
    
    do i = 1,size(v1)
       if (v1(i) < v2(i)) then
          ltVec = .true.
          return
       else if (v1(i) > v2(i)) then
          ltVec = .false.
          return
       end if
    end do
    ltVec = .false.
  end function ltVec
end module track
