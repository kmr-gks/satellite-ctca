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
    real(kind=8),allocatable :: dist(:),energy(:)
    !共有する配列のサイズ
    integer(kind=8)       :: pbuf_size, pbuf_mem=6
    !共有領域のエリアID
    integer               :: pbuf_id,energy_size,i
    !position of satellite
    real(kind=8) :: shipx,shipy,shipz
    real(kind=8) :: neighbour_thr,mass=1
    !environment variables
    character(len=100) :: env_shipy,env_shipz,env_neighbour_thr

contains

    subroutine cotocoa_init
        pbuf_size=size(pbuf)
        !共有する配列の領域を確保
        allocate(dist(pbuf_size))
        allocate(energy(pbuf_size))
        !エリアIDを取得
        call CTCAR_regarea_real8(energy,pbuf_size,pbuf_id)

        ! set parameters from environment variables
        call get_environment_variable("SHIPY",env_shipy)
        read(env_shipy,*) shipy
        call get_environment_variable("SHIPZ",env_shipz)
        read(env_shipz,*) shipz
        call get_environment_variable("NEIGHBOUR_THR",        env_neighbour_thr)
        read(env_neighbour_thr,*) neighbour_thr
    end subroutine cotocoa_init

    subroutine cotocoa_mainstep
        implicit none

        !リクエストを送るときのデータ
        integer(kind=4) ::req_params(10)

        shipx=istep
        dist(:)=sqrt((pbuf(:)%x-shipx)**2+(pbuf(:)%y-shipy)**2+(pbuf(:)%z-shipz)**2)
        energy_size=0
        do i=1, pbuf_size
            if (dist(i).lt.neighbour_thr) then
                energy_size=energy_size+1
                energy(energy_size)=mass*(pbuf(i)%vx**2+pbuf(i)%vy**2+pbuf(i)%vz**2)/2
            end if
        end do

        !リクエスト時のデータを設定
        req_params(1)=myid
        req_params(2)=pbuf_size
        req_params(3)=pbuf_mem
        req_params(4)=istep
        req_params(5)=energy_size
        !リクエスト時にデータを送ることができる
        call CTCAR_sendreq(req_params,size(req_params))
    
    end subroutine cotocoa_mainstep

    subroutine cotocoa_finalize
    
    end subroutine cotocoa_finalize

end module m_ctcamain
