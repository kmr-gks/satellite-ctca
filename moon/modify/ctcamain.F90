module m_ctcamain
    use oh_type
    use paramt
    use allcom
#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
    implicit none
    private
    public cotocoa_init, cotocoa_mainstep, cotocoa_finalize

    !共有するphiの配列
    real(kind=8),allocatable :: pbuf_data(:,:)
    !共有する配列のサイズ
    integer(kind=8)       :: pbuf_size, pbuf_mem=6
    !共有領域のエリアID
    integer               :: pbuf_id

contains

    subroutine cotocoa_init
        pbuf_size=size(pbuf)
        !共有する配列の領域を確保
        allocate(pbuf_data(pbuf_size,pbuf_mem))
        !エリアIDを取得
        call CTCAR_regarea_real8(pbuf_data,pbuf_size*pbuf_mem,pbuf_id)
    end subroutine cotocoa_init

    subroutine cotocoa_mainstep
        implicit none
    
        !リクエストを送るときのデータ
        integer(kind=4) ::req_params(10)
        integer ::from_rank=10

        if(myid.eq.from_rank) then
            !pbuf_posには位置情報を格納
            pbuf_data(:,1)=pbuf(:)%x
            pbuf_data(:,2)=pbuf(:)%y
            pbuf_data(:,3)=pbuf(:)%z
            pbuf_data(:,4)=pbuf(:)%vx
            pbuf_data(:,5)=pbuf(:)%vy
            pbuf_data(:,6)=pbuf(:)%vz
            !print*, "requester: pbuf_vel=", pbuf_data(1:10,:)
            !リクエスト時のデータを設定
            req_params(1)=from_rank
            req_params(2)=pbuf_size
            req_params(3)=pbuf_mem
            !リクエスト時にデータを送ることができる
            call CTCAR_sendreq(req_params,size(req_params))
        end if
    
    end subroutine cotocoa_mainstep

    subroutine cotocoa_finalize
    
    end subroutine cotocoa_finalize

end module m_ctcamain
