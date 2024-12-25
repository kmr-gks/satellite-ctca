program worker
  use mpi
  use ctca
  implicit none
!
  integer :: ierr, myrank, nprocs
  !エリアID
  integer :: pbuf_id
  integer :: i,step=0
  integer :: from_rank
  !リクエストを受け取るときのデータ
  integer ::req_params(10)
  !受け取るデータとそのサイズ
  real*8,allocatable    :: energy(:)
  integer :: pbuf_size, pbuf_mem, energy_size=0
  real(kind=8) :: satellite_pos, grid_length=0.5
  !output file of energy
  character(len=100) :: output_file_name
  integer :: output_file_unit=10,particle_per_rank(130)
!
  call CTCAW_init(0, 1)
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
!
  print*, "worker: ", myrank, " / ", nprocs
  call get_environment_variable("OUTPUT_FILE_NAME",output_file_name)
! エリアIDを取得
  call CTCAW_regarea_real8(pbuf_id)
  !open the output file
  open(unit=output_file_unit,file=output_file_name, status='replace', action='write')
  write(output_file_unit,'(A)') "step,energy"
!
  do while( .true. )
    !リクエストを受けとる
    call CTCAW_pollreq(from_rank,req_params,size(req_params))
    if( CTCAW_isfin() ) exit
    !リクエスト時のデータをもとに受け取るデータの情報をみる
    from_rank = req_params(1)
    pbuf_size = req_params(2)
    pbuf_mem = req_params(3)
    step = req_params(4)
    satellite_pos = step*grid_length
    energy_size = req_params(5)
    !初回のみ領域確保
    if(.not.allocated(energy)) then
      allocate(energy(pbuf_size))
    end if

    !read_dataにデータを読み込む
    call CTCAW_readarea_real8(pbuf_id,from_rank,0,energy_size,energy)
    
    do i=1, energy_size
      write(output_file_unit,'(F,",",F)') satellite_pos,energy(i)
    end do
    particle_per_rank(from_rank) = particle_per_rank(from_rank) + energy_size

    call CTCAW_complete()
  end do
  print*,particle_per_rank
  print*, "worker is finalizing..."
  call CTCAW_finalize()
  close(output_file_unit)
!
  stop
end program worker
