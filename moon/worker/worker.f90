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
  integer :: flag_id,flag_size,flag(6)
  !output file of energy
  character(len=100) :: output_file_name
  integer :: output_file_unit=10,particle_per_rank(130)
    integer date_time(8)
    character(len=100) :: date_str(3)
  integer finished_rank,waiting_rank
!
  call CTCAW_init(0, 1)
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
  flag_size = size(flag)
!
  print*, "worker: ", myrank, " / ", nprocs
  call get_environment_variable("OUTPUT_FILE_NAME",output_file_name)
! get area id
  call CTCAW_regarea_real8(pbuf_id)
  call CTCAW_regarea_int(flag_id)
  !open the output file
  open(unit=output_file_unit,file=output_file_name, status='replace', action='write')
  write(output_file_unit,'(A)') "step,energy"
  !polling request
  call CTCAW_pollreq(from_rank,req_params,size(req_params))
  !allocate energy array
  if(.not.allocated(energy)) then
    allocate(energy(req_params(1)))
  end if
  call CTCAW_complete()
!
  do while( .true. )
    if( CTCAW_isfin() ) exit
    
    finished_rank=0
    waiting_rank=0
    do from_rank=0,127
      call CTCAW_readarea_int(flag_id,from_rank,0,flag_size,flag)
      !print*,"from_rank=",from_rank,"flag=",flag(1)
      if (flag(1).eq.0) then
        !write to csv
        print*,"from_rank=",from_rank,"step=",flag(5)
        !read energy from pbuf
        step=flag(5)
        energy_size=flag(6)
        call CTCAW_readarea_real8(pbuf_id,from_rank,0,energy_size,energy)
        !print*,"CTCAW_readarea_real8 pbuf_id"
        satellite_pos = step*grid_length
        write(output_file_unit, '( *(F,",",F,/) )') (satellite_pos, energy(i), i=1, energy_size)
        particle_per_rank(from_rank) = particle_per_rank(from_rank) + energy_size
        !set flag
        flag(1)=1
        !time up
        call CTCAW_writearea_int(flag_id,from_rank,0,flag_size,flag)
        !print*,"CTCAW_writearea_int flag_id"
      else if (flag(1).eq.1) then
        !wait for requester
        print*, "waiting for requester...", from_rank
        waiting_rank=waiting_rank+1
      else if (flag(1).eq.2) then
        finished_rank=finished_rank+1
      end if
    end do
    if (finished_rank.eq.128) then
      exit
    end if
    if (waiting_rank.gt.0) then
      call sleep(1)
    end if
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
