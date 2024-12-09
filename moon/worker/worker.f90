program worker
  use mpi
  use ctca
  implicit none
!
  integer :: ierr, myrank, nprocs
  !エリアID
  integer :: pbuf_id
  integer :: i
  integer :: from_rank
  !リクエストを受け取るときのデータ
  integer ::req_params(10)
  !受け取るデータとそのサイズ
  real*8,allocatable    :: pbuf_data(:,:)
  integer :: pbuf_size, pbuf_mem
!
  call CTCAW_init(0, 1)
  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
!
  print*, "worker: ", myrank, " / ", nprocs
!
! エリアIDを取得
  call CTCAW_regarea_real8(pbuf_id)
!
  do while( .true. )
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
    !pbuf_velを速度単位に変換
    print*, "CTCAworker: pbuf_data="
    do i = 1, 10
      print*,pbuf_data(i,:)
    end do

    call CTCAW_complete()
  end do

  print*, "worker is finalizing..."
  call CTCAW_finalize()
!
  stop
end program worker
