program coupler_test
  use mpi
  use ctca
  implicit none
!
  integer(kind=4) :: ierr, myrank, nprocs
  integer(kind=4) :: dareaid, iareaid, phi_areaid
  integer(kind=4) :: reqinf(4)
  integer(kind=4) :: frmrank, progid, ndat=10
  integer(kind=4) :: datint(10)
!
  call CTCAC_init()
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
!
! 領域確保 ハンドルを取得する感じ
  call CTCAC_regarea_real4(dareaid)
  call CTCAC_regarea_int(iareaid)
  call CTCAC_regarea_real8(phi_areaid)
!
  do while (.true.)
    call CTCAC_pollreq(reqinf,frmrank,datint,ndat)
    if( CTCAC_isfin() ) exit
    progid = 0
    call CTCAC_enqreq(reqinf,progid,datint,ndat)
  end do
!
  call CTCAC_finalize()
!
  stop
end program coupler_test
