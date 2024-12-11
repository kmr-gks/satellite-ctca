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
  real*8,allocatable    :: pbuf_data(:,:)
  integer :: pbuf_size, pbuf_mem
  !position of satellite
  real(kind=8) :: shipx,shipy,shipz
  real(kind=8) :: dist,mass=1,energy
  !output file of energy
  character(len=100) :: output_file_name="output.csv"
  integer :: output_file_unit=10
!
  call CTCAW_init(0, 1)
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
!
  print*, "worker: ", myrank, " / ", nprocs
!
! エリアIDを取得
  call CTCAW_regarea_real8(pbuf_id)
  !open the output file
  open(unit=output_file_unit,file=output_file_name, status='replace', action='write')
  write(output_file_unit,'(A)') "step,energy"
!
  do while( .true. )
    step = step + 1
    !リクエストを受けとる
    call CTCAW_pollreq(from_rank,req_params,size(req_params))
    if( CTCAW_isfin() ) exit
    !リクエスト時のデータをもとに受け取るデータの情報をみる
    from_rank = req_params(1)
    pbuf_size = req_params(2)
    pbuf_mem = req_params(3)
    !初回のみ領域確保
    if(.not.allocated(pbuf_data)) then
      allocate(pbuf_data(pbuf_size,pbuf_mem))
    end if

    !read_dataにデータを読み込む
    call CTCAW_readarea_real8(pbuf_id,from_rank,0,pbuf_size*pbuf_mem,pbuf_data)

    !position of satellite
    shipx=step
    shipy=5
    shipz=145
    
    do i=1, pbuf_size
      !check the distance between satellite and the object
      dist=sqrt((pbuf_data(i,1)-shipx)**2+(pbuf_data(i,2)-shipy)**2+(pbuf_data(i,3)-shipz)**2)
      if (dist.lt.80) then
        !calculate the energy
        energy=mass*0.5*(pbuf_data(i,4)**2+pbuf_data(i,5)**2+pbuf_data(i,6)**2)
        write(output_file_unit,'(I,",",F)') step,energy
      end if
    end do

    !print*, "CTCAworker: pbuf_data="
    !do i = 1, 10
    !  print*,pbuf_data(i,:)
    !end do

    call CTCAW_complete()
  end do

  print*, "worker is finalizing..."
  call CTCAW_finalize()
  close(output_file_unit)
!
  stop
end program worker
