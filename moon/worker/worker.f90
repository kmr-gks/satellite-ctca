module common_module
    use iso_c_binding
    implicit none

    ! usleep関数のインターフェース
    interface
        subroutine usleep(microseconds) bind(C, name='usleep')
            import :: C_INT
            integer(C_INT), value :: microseconds
        end subroutine usleep
    end interface
contains
    subroutine mysleep(seconds)
        real, intent(in) :: seconds
        call usleep(int(seconds*1000000, C_INT))
    end subroutine mysleep
end module common_module
program worker
  use mpi
  use ctca
  use common_module
  implicit none
!
  integer :: ierr, myrank, nprocs
  !area id
  integer :: i,j,k,l,step=0
  integer :: from_rank
  !data for request
  integer ::req_params(10)
  real(kind=8) :: req_params_real(10)
  integer :: pbuf_size, pbuf_mem, real_par_num_per_sup_par
  !time is real unit
  real(kind=8) :: grid_length=0.5,time_ratio,time
  !flag of completion
  integer :: flag_id,flag_size,flag(10)
  !output file of energy
  character(len=100) :: output_file_name
  integer :: output_file_unit=10
    integer date_time(8)
    character(len=100) :: date_str(3)
  integer finished_rank,waiting_rank
  !number of super particles per time(sec), energy(log10 eV), and species(1or2)
  integer,allocatable    :: num_par(:,:),num_par_total(:,:,:),num_par_v(:,:,:),num_par_v_total(:,:,:,:)
  integer num_par_id,num_par_v_id,energy_bin,spec_num,v_dim,nstep
  character(len=100) :: format_string
!
  call CTCAW_init(0, 1)
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
  flag_size = size(flag)
!
  print*, "worker: ", myrank, " / ", nprocs
  call get_environment_variable("OUTPUT_FILE_NAME",output_file_name)
! get area id
  call CTCAW_regarea_int(flag_id)
  call CTCAW_regarea_int(num_par_id)
  call CTCAW_regarea_int(num_par_v_id)
  !open the output file
  open(unit=output_file_unit,file=output_file_name, status='replace', action='write',buffered='yes')
  write(output_file_unit,'(A)') "time,species,energy(10*log10eV),par-count,vx-count,vy-count,vz-count"
  !polling request
  call ctcaw_pollreq_withreal8(from_rank,req_params,size(req_params),req_params_real,size(req_params_real))
  print*,"req_params_real(1)=",req_params_real(1)
  time_ratio=req_params_real(1)
  !allocate energy array
   energy_bin=req_params(3)
  spec_num=req_params(4)
  nstep=req_params(5)
  real_par_num_per_sup_par=req_params(6)
  v_dim=req_params(7)
  allocate(num_par(-energy_bin:energy_bin,spec_num))
  allocate(num_par_total(-energy_bin:energy_bin,spec_num,nstep))
  num_par_total=0
  allocate(num_par_v(-energy_bin:energy_bin,v_dim,spec_num))
  allocate(num_par_v_total(-energy_bin:energy_bin,v_dim,spec_num,nstep))
  num_par_v_total=0
  call CTCAW_complete()
!
  do while( .true. )
    if( CTCAW_isfin() ) exit
    
    finished_rank=0
    waiting_rank=0
    do from_rank=0,127
      call CTCAW_readarea_int(flag_id,from_rank,0,flag_size,flag)
      if (flag(1).eq.0) then
        !write to csv
        !print*,"from_rank=",from_rank,"step=",flag(5)
        !read energy from pbuf
        step=flag(5)
        time=step/time_ratio

        call CTCAW_readarea_int(num_par_id,from_rank,0,size(num_par),num_par)
        call CTCAW_readarea_int(num_par_v_id,from_rank,0,size(num_par_v),num_par_v)
        !set flag
        flag(1)=1
        !time up
        call CTCAW_writearea_int(flag_id,from_rank,0,flag_size,flag)
        num_par_total(:,:,step)=num_par_total(:,:,step)+num_par(:,:)
        num_par_v_total(:,:,:,step)=num_par_v_total(:,:,:,step)+num_par_v(:,:,:)
      else if (flag(1).eq.1) then
        !wait for requester
        waiting_rank=waiting_rank+1
      else if (flag(1).eq.2) then
        finished_rank=finished_rank+1
      end if
    end do
    if (finished_rank.gt.64) then
      !if half of the ranks are finished, worker will exit
      exit
    end if
    if (waiting_rank.eq.128) then
      call mysleep(0.01)
    end if
  end do
  print*, "worker is writing data to file"
  call system("date")
  num_par_total=num_par_total*real_par_num_per_sup_par
  num_par_v_total=num_par_v_total*real_par_num_per_sup_par
  ! create dynamic format string
  format_string = '( *(G0, ",", I4, ",", I4, ",", I,' // repeat('",", I,', v_dim) // '/) )'
  write(output_file_unit, format_string) &
    (((real(i)/time_ratio, j, k, num_par_total(k, j, i), &
        (num_par_v_total(k, l, j, i), l = 1, v_dim), &
        k = -energy_bin, energy_bin), j = 1, spec_num), i = 1, nstep)

  call system("date")
  call CTCAW_finalize()
  close(output_file_unit)
!
  stop
end program worker
