program sample_requester
  use mpi
  use ctca
  implicit none

  integer(kind=4),parameter :: nstep = 100
  integer(kind=4) :: dataintnum, myrank, nprocs, ierr, areaid(2), progid, fromrank
  integer(kind=4),dimension(4) :: reqinfo
  integer(kind=4),dimension(2) :: dataint
  integer(kind=4) :: istep
  real(kind=8),dimension(6,400,10) :: datareal8
  real(kind=8),dimension(1),volatile :: rtowdat=0, wtordat=0
  real(kind=8) :: tick_rtow = 0.0d0
  integer(kind=4) :: i, j, k
  integer(kind=4), parameter :: prognum = 1

  call CTCAR_init()
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
  if (myrank == 0) print *, "requester init done"

  call CTCAR_regarea_real8(rtowdat, 1, areaid(1))
  call CTCAR_regarea_real8(wtordat, 1, areaid(2))

  progid = 1
  dataint(1) = progid
  dataint(2) = nstep
  tick_rtow = 1.0d0/nstep

  if (myrank == 0) then
    print *, "requester sendreq : ", i
    call CTCAR_sendreq(dataint, 2)
    print *, "requester sendreq done"
  end if

  if (myrank == 0) wtordat(1) = -1.0d0
  do istep = 1, nstep
    if (myrank == 0) then
      print *, "req-step:", istep

      ! R1
      do while (.true.)
        if (rtowdat(1) < 0.0d0) exit ! R1C
      end do

      ! R2
      rtowdat(1) = istep*tick_rtow
      print *, "requester put data", rtowdat(1)

      ! R3
      do while (.true.)
        if (wtordat(1) >= 0.0d0) exit
      end do

      ! R4
      print *, "requester got data from w1", wtordat(1)

      ! R5
      wtordat(1) = -1.0d0

!      call sleep(1)
    end if

!   computation body processed within each iteration
    call MPI_Barrier(CTCA_subcomm, ierr)

  end do

  call CTCAR_finalize()

end program

