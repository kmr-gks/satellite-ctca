program sample_coupler
  use mpi
  use ctca
  implicit none

  integer(kind=4) :: dataintnum, myrank, nprocs, ierr, areaid(2), progid, fromrank
  integer(kind=4),dimension(4) :: reqinfo
  integer(kind=4),dimension(2) :: dataint
  real(kind=8),dimension(6,400) :: datareal8
  logical :: CTCAC_verbose = .false.

  call CTCAC_init()
  if(CTCAC_verbose) print *, "coupler init done"

  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)

  call CTCAC_regarea_real8(areaid(1))
  call CTCAC_regarea_real8(areaid(2))

  do while (.true.)
    call CTCAC_pollreq(reqinfo, fromrank, dataint, 2)
    if(CTCAC_verbose) print *, "coupler poll req done"

    if (CTCAC_isfin()) then
      if(CTCAC_verbose) print *, "coupler got fin"
      exit
    end if

!    call CTCAC_readarea_real8(areaid, fromrank, (dataint(2)-1)*6*400, 6*400, datareal8)
!    if(CTCAC_verbose) print *, "coupler read data"

    if (fromrank >= 0) then
      if(CTCAC_verbose) print *, "coupler got req"
      ! decide progid
      !progid = dataint(1)
      progid = 1

      ! enqueue request with data
      call CTCAC_enqreq_withreal8(reqinfo, progid, dataint, 2, datareal8, 6*400)
      if(CTCAC_verbose) print *, "coupler enqreq done"

    end if
  end do

  call CTCAC_finalize()

end program

