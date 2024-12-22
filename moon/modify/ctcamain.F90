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
        implicit none
        !plasma frequency for ion (real unit)
        real(kind=8) :: wp_ion_emses,wp_ion_real
        real(kind=8) :: len_ratio,vel_ratio,time_ratio,freq_ratio,grid_length=0.5,ion_density
        real(kind=8) :: ion_charge=1.60217662d-19,ion_mass=1.6726219d-27,permittivity=8.854187817d-12
        real(kind=8) :: simu_vol,sup_par_mass
        integer sup_par_num

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
        call get_environment_variable("NEIGHBOUR_THR",env_neighbour_thr)
        read(env_neighbour_thr,*) neighbour_thr

        ! set plasma frequency for ion
        wp_ion_emses=wp(2)
        len_ratio=1/grid_length
        vel_ratio=cv*1e3/2.99792458d8
        time_ratio=len_ratio/vel_ratio
        freq_ratio=1/time_ratio
        wp_ion_real=wp_ion_emses/freq_ratio
        ion_density=wp_ion_real**2*ion_mass*permittivity/ion_charge**2/1e6
        simu_vol=nx*ny*nz*grid_length**3
        !super particle number
        sup_par_num=nodes(1)*nodes(2)*nodes(3)*pbuf_size
        sup_par_mass=ion_mass*(simu_vol*1d6*ion_density/sup_par_num)
        if (myid.eq.0) then
            print *, "wp_ion_emses=",wp_ion_emses
            print *, "wp_ion_real=",wp_ion_real,"[Hz]"
            print *, "ion_density=",ion_density,"[/cc]"
            print *, "simu_vol=",simu_vol,"[m^3]"
            print *, "sup_par_num=",sup_par_num
            print *, "sup_par_mass=",sup_par_mass,"[kg]"
        end if

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
