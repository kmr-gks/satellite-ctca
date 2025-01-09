program sample_worker
  use mpi
  use ctca
  implicit none

  integer(kind=4) :: nstep
  integer(kind=4) :: dataintnum, myrank, nprocs, ierr, areaid(2), progid, fromrank
  integer(kind=4) :: world_myrank
  integer(kind=4),dimension(2) :: dataint
  integer(kind=4) :: istep
  real(kind=8),dimension(6,400) :: datareal8
  real(kind=8),dimension(1) :: rtowd, wtord
  real(kind=8) :: tick_wtor = 1.0d0
  character(len=256) :: fname
  integer(kind=4) :: i, j, c
  
  call CTCAW_init(1, 8)
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, world_myrank, ierr)
  if(myrank.eq.0) print *, "worker1 init done"

  call CTCAW_regarea_real8(areaid(1))
  call CTCAW_regarea_real8(areaid(2))

  do while (.true.)
    call CTCAW_pollreq_withreal8(fromrank, dataint, 2, datareal8, 6*400)
    nstep = dataint(2)

    if (CTCAW_isfin()) then
      exit
    end if

    if (myrank == 0) print *, "worker1 poll req done", fromrank

    rtowd(1) = -1.0d0
    if (myrank == 0) then
      print *, "wrk1 first write start"
      call CTCAW_writearea_real8(areaid(1), fromrank, 0, 1, rtowd)
      print *, "wrk1 first write end"
    end if
    do istep = 1, nstep
      if (myrank == 0) then
        print *, "wrk1-step:", istep

        ! W1
        do while (.true.)
          call CTCAW_readarea_real8(areaid(1), fromrank, 0, 1, rtowd)
          if (rtowd(1) >= 0.0d0) exit ! W1C
        end do

        ! W2
        print *, "worker1 got data from req.", rtowd(1)
 
        ! W3
        rtowd(1) = -1.0d0
        call CTCAW_writearea_real8(areaid(1), fromrank, 0, 1, rtowd)

        ! W4
        do while (.true.)
          call CTCAW_readarea_real8(areaid(2), fromrank, 0, 1, wtord)
          if (wtord(1) < 0.0d0) exit ! W4C
        end do

        ! W5
        wtord(1) = istep*tick_wtor
        call CTCAW_writearea_real8(areaid(2), fromrank, 0, 1, wtord)
        print *, "worker1 put data", wtord(1)

!        call sleep(1)
      end if

!     computation body processed within each iteration
      call MPI_Barrier(CTCA_subcomm, ierr)

    end do

    ! end this work
    call CTCAW_complete()
    if (myrank == 0) print *, world_myrank, "worker1 complete done"
  end do

  call CTCAW_finalize()

end program sample_worker

