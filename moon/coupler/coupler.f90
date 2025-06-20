program coupler
  use mpi
  use ctca
  implicit none
!
  integer(kind=4) :: ierr, myrank, nprocs
  !エリアID
  integer(kind=4) :: reqinf(4)
  integer(kind=4) :: frmrank, progid
  integer(kind=4) :: req_params(10)
  real(kind=8) :: req_params_real(10)
  integer :: flag_id,num_par_id,num_par_v_id
!
  call CTCAC_init()
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
!
! エリアIDを取得
  call CTCAC_regarea_int(flag_id)
  call CTCAC_regarea_int(num_par_id)
  call CTCAC_regarea_int(num_par_v_id)
!
  do while (.true.)
    call CTCAC_pollreq_withreal8(reqinf,frmrank,req_params,size(req_params),req_params_real,size(req_params_real))
    if( CTCAC_isfin() ) exit
    progid = 0
    call CTCAC_enqreq_withreal8(reqinf,progid,req_params,size(req_params),req_params_real,size(req_params_real))
  end do
!
  call CTCAC_finalize()
!
  stop
end program coupler
