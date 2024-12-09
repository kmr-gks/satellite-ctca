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
    real*8,allocatable    :: phi_data(:)
    !共有する配列のサイズ
    integer(kind=8)       :: phi_data_size=35
    !共有領域のエリアID
    integer               :: phi_areaid

contains

    subroutine cotocoa_init
        !共有する配列の領域を確保
        allocate(phi_data(phi_data_size))
        !エリアIDを取得
        call CTCAR_regarea_real8(phi_data,phi_data_size,phi_areaid)
    end subroutine cotocoa_init

    subroutine cotocoa_mainstep
        implicit none
    
        !リクエストを送るときのデータ
        integer(kind=4) ::req_params(10)
        integer(kind=4) ::x,y,z,phi_shape(5),i
        integer ::from_rank=10

        if(myid.eq.from_rank) then
            phi_shape=shape(phi)
            x=phi_shape(2)
            y=phi_shape(3)
            z=phi_shape(4)
            phi_data = phi(1,1:phi_data_size,y/2,z/2,1)
            print*, "requester: phi_data=", phi_data
            !リクエスト時のデータを設定
            req_params(1)=from_rank
            req_params(2)=phi_data_size
            !リクエスト時にデータを送ることができる
            call CTCAR_sendreq(req_params,size(req_params))
        end if
    
    end subroutine cotocoa_mainstep

    subroutine cotocoa_finalize
    
    end subroutine cotocoa_finalize

end module m_ctcamain
