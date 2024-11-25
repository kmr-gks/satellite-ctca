program dout
  use mpi
  use ctca
  implicit none
!
  type pinfo
    integer(kind=4) :: rank, pgrid(2), offset(2), dsize(2), stat
  end type
  integer(kind=4) :: ierr, myrank, nprocs, len
  integer(kind=4) :: phi_areaid,phi_area_size=35
  integer(kind=4) :: frmrank, hdat(10), ndat=10
  integer(kind=4) :: nsdom, lsdom(2)
  integer(kind=4) :: wflag(10) = 0
  integer(kind=4) :: nread, iread, jread
  integer(kind=4) :: szx, szy, szt, osx, osy
  real(kind=4),allocatable :: odat(:,:), sdat(:,:)
  type(pinfo),allocatable :: ptable(:)
  character(128) :: hostname
  character(7) :: ofname
  logical :: vb=.false.
  integer(kind=4) :: i, j, k, m
  integer :: socket, ret, comlen
  real*8,allocatable    :: read_data(:)
!
  call CTCAW_init(0, 1)
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
  call MPI_Get_processor_name(hostname, len, ierr)
!
  print*, "dout: ", myrank, " / ", nprocs, " : ", hostname
!
  call CTCAW_regarea_real8(phi_areaid)

!
  allocate(read_data(phi_area_size))
  do while( .true. )
    call CTCAW_pollreq(frmrank,hdat,ndat)
    if( CTCAW_isfin() ) exit
    !read_areaしてポテンシャルを取得する
    call CTCAW_readarea_real8(phi_areaid,51,0,phi_area_size,read_data)
    print*, "CTCAworker: read_data="!, read_data
    do i = 1, size(read_data)
      write(*, "(A,E)", advance="no") ",", read_data(i)
    end do

    call CTCAW_complete()
  end do

  print*, "dout is finalizing..."
  call CTCAW_finalize()
!
  stop
end program dout
