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
  integer :: ierr, myrank, nprocs, nprocs_reqester
  !area id
  integer :: i,j,k,l,step=0
  integer :: from_rank
  !data for request
  integer ::req_params(10)
  real(kind=8) :: req_params_real(10)
  integer :: pbuf_size, pbuf_mem, real_par_num_per_sup_par
  !time is real unit
  real(kind=8) :: grid_length,time_ratio,time,neighbour_vol_real, energy_extent
  !flag of completion
  integer :: flag_id,flag_size,flag(10)
  !output file of energy
  character(len=100) :: output_file_name
  integer :: output_file_unit=10
    integer date_time(8)
    character(len=100) :: date_str(3)
  integer finished_rank,waiting_rank
  !number of super particles per time(sec), energy(log10 eV), and species(1or2)
  integer,allocatable    :: num_par(:,:),num_par_v(:,:,:)
  real(kind=8),allocatable :: num_par_total(:,:,:),num_par_v_total(:,:,:,:)
  !step_csv time-step number for csv output
  integer :: num_par_id,num_par_v_id,energy_bin,spec_num,v_dim,nstep,step_csv=100,step_from,step_to
  character(len=1000) :: format_string
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
  write(output_file_unit,'(A)') "time,species,energy(10*log10eV),par-count,vx-count,vy-count,vz-count,vx-p,vy-p,vz-p,vx-n,vy-n,vz-n"
  !polling request
  call ctcaw_pollreq_withreal8(from_rank,req_params,size(req_params),req_params_real,size(req_params_real))
  print*,"req_params_real(1)=",req_params_real(1)
  time_ratio=req_params_real(1)
  neighbour_vol_real=req_params_real(2)
  grid_length=req_params_real(3)
  print*, "neighbour_vol_real=",neighbour_vol_real,"[m^3]", "grid_length=",grid_length,"[m]"
  !allocate energy array
  energy_bin=req_params(3)
  spec_num=req_params(4)
  nstep=req_params(5)
  real_par_num_per_sup_par=req_params(6)
  v_dim=req_params(7)
  nprocs_reqester=req_params(8)
  step_from=req_params(9)
  step_to=req_params(10)
  allocate(num_par(-energy_bin:energy_bin,spec_num))
  allocate(num_par_total(-energy_bin:energy_bin,spec_num,step_csv))
  num_par_total=0
  allocate(num_par_v(-energy_bin:energy_bin,v_dim,spec_num))
  allocate(num_par_v_total(-energy_bin:energy_bin,v_dim,spec_num,step_csv))
  num_par_v_total=0
  call CTCAW_complete()
!
  do while( .true. )
    if( CTCAW_isfin() ) exit
    
    finished_rank=0
    waiting_rank=0
    do from_rank=0,nprocs_reqester-1
      call CTCAW_readarea_int(flag_id,from_rank,0,flag_size,flag)
      if (flag(1).eq.0) then
        !write to csv
        !print*,"from_rank=",from_rank,"step=",flag(5)
        !read energy from pbuf
        step=flag(5)
        time=step/time_ratio
        if (step_from.lt.step.and.step.le.step_to) then
          call CTCAW_readarea_int(num_par_id,from_rank,0,size(num_par),num_par)
          call CTCAW_readarea_int(num_par_v_id,from_rank,0,size(num_par_v),num_par_v)
        end if
        !set flag
        flag(1)=1
        !time up
        call CTCAW_writearea_int(flag_id,from_rank,0,flag_size,flag)
        if (step_from.lt.step.and.step.le.step_to) then
          if (flag(2).eq.0) then
            print*,"from_rank=",flag(2),"step=",step,"index=",ceiling(real(step_csv)*(step-step_from)/(step_to-step_from+1))
          end if
          step=ceiling(real(step_csv)*(step-step_from)/(step_to-step_from+1))
          ! Measures for floating point calculation errors
          step=max(step,1)
          step=min(step,step_csv)
          num_par_total(:,:,step)=num_par_total(:,:,step)+num_par(:,:)
          num_par_v_total(:,:,:,step)=num_par_v_total(:,:,:,step)+num_par_v(:,:,:)
        end if
      else if (flag(1).eq.1) then
        !wait for requester
        waiting_rank=waiting_rank+1
      else if (flag(1).eq.2) then
        finished_rank=finished_rank+1
      end if
    end do
    if (finished_rank.gt.nprocs_reqester/2) then
      !if half of the ranks are finished, worker will exit
      exit
    end if
  end do
  call system("date")
  !transform number of super particles to number of real particles
  num_par_total=num_par_total*real_par_num_per_sup_par
  num_par_v_total=num_par_v_total*real_par_num_per_sup_par
  !take average by step (time)
  !配列の時間の次元について、複数のステップのデータが1行に対応するときは対応ステップ数で割り平均を取る
  num_par_total=num_par_total/(nstep/step_csv)
  num_par_v_total=num_par_v_total/(nstep/step_csv)
  !devide by volume of observation range
  num_par_total=num_par_total/neighbour_vol_real
  !normalize
  num_par_v_total = num_par_v_total / (spread(sum(num_par_v_total, dim=1), dim=1, ncopies=size(num_par_v_total, 1))+1)

  !devide by extent of energy bin in log10 scale
  do i=-energy_bin,energy_bin
    energy_extent=10**((i+1)/10.0)-10**(i/10.0)
    !print*,"energy_extent=10^",(i+1)/10.0-"10^",i/10.0,"=",10**((i+1)/10.0),"-",10**(i/10.0),"=",energy_extent
    num_par_total(i,:,:) = num_par_total(i,:,:)/energy_extent
  end do
  ! create dynamic format string
  format_string = '( *(G0, ",", I4, ",", I4, ",", G0,' // repeat('",", G0,', v_dim) // '/) )'
  write(output_file_unit, format_string) &
    ((((real(i)/step_csv*(step_to-step_from+1)+step_from)/time_ratio, j, k, num_par_total(k, j, i), &
        (num_par_v_total(k, l, j, i), l = 1, v_dim), &
        k = -energy_bin, energy_bin), j = 1, spec_num), i = 1, step_csv)

  call system("date")
  call CTCAW_finalize()
  close(output_file_unit)
!
  stop
end program worker
