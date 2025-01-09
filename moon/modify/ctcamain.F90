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
module m_ctcamain
    use oh_type
    use paramt
    use allcom
    use common_module
#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
    implicit none
    private
    public cotocoa_init, cotocoa_mainstep, cotocoa_finalize
    !elementary charge[C], ion mass[kg], electron mass[kg], permittivity of vacuum[F/m]
    real(kind=8) :: ion_charge=1.6021766d-19,ion_mass=1.6726219d-27,electron_mass=9.109383d-31,permittivity=8.85418781d-12
    !ratio of emses to real unit, ion density[/cc]
    real(kind=8) :: len_ratio,vel_ratio,time_ratio,freq_ratio,ion_density

    !distance between satellite and particles, energy density
    real(kind=8),allocatable :: dist(:),energy(:)
    integer,allocatable :: species(:)
    !size of pbuf,members of pbuf
    integer(kind=8)       :: pbuf_size, pbuf_mem=6
    !area id of pbuf, size of energy
    integer               :: pbuf_id,species_id,energy_size,i
    !flag of completion
    integer :: flag_id,flag_size,flag(6)
    !position of satellite (emses unit)
    real(kind=8) :: shipx,shipy,shipz
    !neighbour threshold, super particle mass, grid length, neighbour volume
    real(kind=8) :: neighbour_thr,sup_par_mass,grid_length=0.5,neighbour_vol
    !environment variables
    character(len=100) :: env_shipy,env_shipz,env_neighbour_thr
    !data for request
    integer(kind=4) ::req_params(10)
    real(kind=8) :: req_params_real(10)
    !number of super particles per energy(10*log10eV), and species(1or2)
    integer :: num_par(-100:100,2),num_par_id

contains

    subroutine cotocoa_init
        implicit none
        !plasma frequency for ion (emses, real unit)
        real(kind=8) :: wp_ion_emses,wp_ion_real
        !simulation volume
        real(kind=8) :: simu_vol
        !number of super particles
        integer sup_par_num,real_par_num_per_sup_par

        pbuf_size=size(pbuf)
        flag_size = size(flag)
        allocate(dist(pbuf_size))
        allocate(energy(pbuf_size))
        allocate(species(pbuf_size))
        call CTCAR_regarea_real8(energy,pbuf_size,pbuf_id)
        call CTCAR_regarea_int(species,pbuf_size,species_id)
        call CTCAR_regarea_int(flag,flag_size,flag_id)
        call CTCAR_regarea_int(num_par,size(num_par),num_par_id)

        ! set parameters from environment variables
        call get_environment_variable("SHIPY",env_shipy)
        read(env_shipy,*) shipy
        call get_environment_variable("SHIPZ",env_shipz)
        read(env_shipz,*) shipz
        call get_environment_variable("NEIGHBOUR_THR",env_neighbour_thr)
        read(env_neighbour_thr,*) neighbour_thr

        ! get plasma frequency for ion
        wp_ion_emses=wp(2)

        len_ratio=1/grid_length
        vel_ratio=cv*1e3/2.99792458d8
        time_ratio=len_ratio/vel_ratio
        freq_ratio=1/time_ratio

        wp_ion_real=wp_ion_emses/freq_ratio
        ion_density=wp_ion_real**2*ion_mass*permittivity/ion_charge**2/1e6
        simu_vol=nx*ny*nz*grid_length**3

        sup_par_num=nodes(1)*nodes(2)*nodes(3)*pbuf_size
        real_par_num_per_sup_par=simu_vol*1d6*ion_density/sup_par_num
        sup_par_mass=ion_mass*real_par_num_per_sup_par
        neighbour_vol=4/3*pi*(neighbour_thr*grid_length)**3
        if (myid.eq.0) then
            print *, "wp_ion_emses=",wp_ion_emses
            print *, "wp_ion_real=",wp_ion_real,"[Hz]"
            print *, "ion_density=",ion_density,"[/cc]"
            print *, "simu_vol=",simu_vol,"[m^3]"
            print *, "sup_par_num=",sup_par_num
            print *, "real_par_num_per_sup_par=",real_par_num_per_sup_par,"[1]"
            print *, "sup_par_mass=",sup_par_mass,"[kg]"
        end if
        !send request
        if (myid.eq.0) then
            req_params(1)=pbuf_size
            req_params_real(1)=time_ratio
            print*,"req_params_real(1)=",req_params_real(1)
            call CTCAR_sendreq_withreal8(req_params,size(req_params),req_params_real,size(req_params_real))
        end if
        flag(1)=1

    end subroutine cotocoa_init

    subroutine cotocoa_mainstep
        implicit none
        logical status1
        integer status2(MPI_STATUS_SIZE),ierr
        integer energy_index,j
        
        num_par(:,:)=0

        !set position of satellite
        shipx=istep/10

        dist(:)=sqrt((pbuf(:)%x-shipx)**2+(pbuf(:)%y-shipy)**2+(pbuf(:)%z-shipz)**2)
        
        !check flag
        do while (.true.)
            if (flag(1).eq.0) then
                !wait for worker
                !print *, "wait for worker",myid
                call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status1, status2, ierr)
                call mysleep(0.01)
            else
                exit
            end if
        end do

        !calculate energy
        energy_size=0
        do i=1, pbuf_size
            if (dist(i).lt.neighbour_thr) then
                energy_size=energy_size+1
                species(energy_size)=pbuf(i)%spec
                !energy(eV)
                if (species(energy_size).eq.1) then
                    energy(energy_size)=electron_mass*(pbuf(i)%vx**2+pbuf(i)%vy**2+pbuf(i)%vz**2)/(vel_ratio**2)/2/ion_charge
                else
                    energy(energy_size)=ion_mass*(pbuf(i)%vx**2+pbuf(i)%vy**2+pbuf(i)%vz**2)/(vel_ratio**2)/2/ion_charge
                end if
                if (energy(energy_size).gt.0) then
                    energy_index=int(10*log10(energy(energy_size)))
                    num_par(energy_index,species(energy_size))=num_par(energy_index,species(energy_size))+1
                end if
            end if
        end do

        !set flag
        flag(1)=0
        flag(2)=myid
        flag(3)=pbuf_size
        flag(4)=pbuf_mem
        flag(5)=istep
        flag(6)=energy_size
    end subroutine cotocoa_mainstep

    subroutine cotocoa_finalize
    integer date_time(8)
    character(len=100) :: date_str(3)
    flag(1)=2
    if (myid.eq.0) print *, "requester finished:",myid
    end subroutine cotocoa_finalize

end module m_ctcamain
