program worker
  use mpi
  use ctca
  implicit none
!
  integer :: ierr, myrank, nprocs
  !area id
  integer :: pbuf_id
  integer :: i,step=0
  integer :: from_rank
  !data for request
  integer ::req_params(10)
  real*8,allocatable    :: energy(:)
  integer :: pbuf_size, pbuf_mem, energy_size=0
  real(kind=8) :: satellite_pos, grid_length=0.5
  !flag of completion
  integer :: flag_id,flag_size=1,flag(1)
  !output file of energy
  character(len=100) :: output_file_name
  integer :: output_file_unit=10,particle_per_rank(130)
    integer date_time(8)
    character(len=100) :: date_str(3)
!
  call CTCAW_init(0, 1)
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
!
  print*, "worker: ", myrank, " / ", nprocs
  call get_environment_variable("OUTPUT_FILE_NAME",output_file_name)
! get area id
  call CTCAW_regarea_real8(pbuf_id)
  call CTCAW_regarea_int(flag_id)
  !open the output file
  open(unit=output_file_unit,file=output_file_name, status='replace', action='write')
  write(output_file_unit,'(A)') "step,energy"
!
  do while( .true. )
    !polling request
    call CTCAW_pollreq(from_rank,req_params,size(req_params))
    if( CTCAW_isfin() ) exit
    !get information from request
    from_rank = req_params(1)
    pbuf_size = req_params(2)
    pbuf_mem = req_params(3)
    step = req_params(4)
    satellite_pos = step*grid_length
    energy_size = req_params(5)
    !allocate energy array
    if(.not.allocated(energy)) then
      allocate(energy(pbuf_size))
    end if

    !read energy from pbuf
    call CTCAW_readarea_real8(pbuf_id,from_rank,0,energy_size,energy)
    call CTCAW_readarea_int(flag_id,from_rank,0,flag_size,flag)
    print*,"from_rank=",from_rank,"flag=",flag
    
    write(output_file_unit, '( *(F,",",F,/) )') (satellite_pos, energy(i), i=1, energy_size)

    particle_per_rank(from_rank) = particle_per_rank(from_rank) + energy_size

    call CTCAW_complete()
  end do
  print*,particle_per_rank
  print*, "worker is finalizing..."
  call date_and_time(date_str(1),date_str(2),date_str(3),date_time)
   print '("worker: ",I4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)',date_time(1),date_time(2),date_time(3),date_time(5),date_time(6),date_time(7)
  call CTCAW_finalize()
  close(output_file_unit)
!
  stop
end program worker
