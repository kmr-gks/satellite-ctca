#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine ccinit
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   C C I N I T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer,parameter :: PRECI=4
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ngx,ngy,ngz, ngt
  integer(kind=4) :: hdat(10)=0, ndat=10


!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)
      ngx = xu - xl; ngy = yu - yl; ngz = zu - zl
      if(ivplane.eq.1) then
        ngt = ngx*ngy
        allocate(datopn(0:ngx-1,0:ngy-1,10))
      else if(ivplane.eq.2) then
        ngt = ngy*ngz
        allocate(datopn(0:ngy-1,0:ngz-1,10))
      else if(ivplane.eq.3) then
        ngt = ngx*ngz
        allocate(datopn(0:ngx-1,0:ngz-1,10))
      end if

!-------------------- 
      allocate(phi_ctca(phi_area_size))
      call CTCAR_regarea_real4(datopn(:,:,1), ngt, dareaid)
      call CTCAR_regarea_int(wflag, 10, iareaid)
      !call CTCAR_regarea_real8(phi(1,1:phi_area_size,0,0,1),phi_area_size,phi_areaid)
      call CTCAR_regarea_real8(phi_ctca,phi_area_size,phi_areaid)
      

!-------------------- 
      datopn(:,:,1) = 0.0d0
      wflag(:) = 0
      wflag(2) = myid


!-------------------- 
      if(myid.eq.0) then
        hdat(1) = nstep/ifdiag + 1
        hdat(2) = ivplane
        hdat(3) = ivcut
        hdat(4:6) = (/nx,ny,nz/)
        hdat(7:9) = nodes(1:3)
        call CTCAR_sendreq(hdat,ndat)
      end if


  return
  end subroutine ccinit


!
  subroutine wrtarea
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   W R T A R E A
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer,parameter :: PRECI=4
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ngx,ngy,ngz


!--------------------
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)
      ngx = xu - xl; ngy = yu - yl; ngz = zu - zl


!--------------------
      do while(.true.)
        if(wflag(1).ne.0) cycle
!
        ivsnap = ivsnap + 1
        if(ivplane.eq.1) then
          if(ivcut.ge.zl.and.ivcut.lt.zu) then
            datopn(0:ngx-1,0:ngy-1,1) = &
           &  real(abs(rhoav(1,0:nxsd-1,0:nysd-1,ivcut-zl)/renrho/rho0),PRECI)
            wflag(1) = ivsnap
          end if
        else if(ivplane.eq.2) then
          if(ivcut.ge.xl.and.ivcut.lt.xu) then
            datopn(0:ngy-1,0:ngz-1,1) = &
           &  real(abs(rhoav(1,ivcut-xl,0:nysd-1,0:nzsd-1)/renrho/rho0),PRECI)
            wflag(1) = ivsnap
          end if
        else if(ivplane.eq.3) then
          if(ivcut.ge.yl.and.ivcut.lt.yu) then
            datopn(0:ngx-1,0:ngz-1,1) = &
           &  real(abs(rhoav(1,0:nxsd-1,ivcut-yl,0:nzsd-1)/renrho/rho0),PRECI)
            wflag(1) = ivsnap
          end if
        end if
!
        exit
      end do


  return
  end subroutine wrtarea
