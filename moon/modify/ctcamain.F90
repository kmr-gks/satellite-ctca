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
    real(kind=8) :: len_ratio,vel_ratio,time_ratio,freq_ratio,ion_density,vel_ratio_k

    !distance between satellite and particles, energy density
    real(kind=8),allocatable :: dist(:)
    !size of pbuf,members of pbuf
    integer(kind=8)       :: pbuf_size, pbuf_mem=6
    !area id of pbuf, size of energy
    integer               :: pbuf_id,species_id,energy_size,i
    !flag of completion
    integer :: flag_id,flag_size,flag(10)
    !position of satellite (emses unit)
    real(kind=8) :: ship_x_from,ship_x_to,ship_y_from,ship_y_to,ship_z_from,ship_z_to,shipx,shipy,shipz
    !index of env_buffer
    integer :: p1, p2, p3, p4, p5
    !neighbour threshold, super particle mass, grid length, neighbour volume
    real(kind=8) :: neighbour_thr,sup_par_mass,grid_length,neighbour_vol_real
    integer correct_by_bin_width,step_from,step_to
    !environment variables
    character(len=100) :: env_buffer
    !data for request
    integer(kind=4) ::req_params(10)
    real(kind=8) :: req_params_real(10)
    !number of super particles per energy(10*log10eV), and species(1or2)
    integer,allocatable :: num_par(:,:),num_par_v(:,:,:)
    integer :: num_par_id,num_par_v_id,energy_bin=100,spec_num=2,v_dim=9

contains

    subroutine cotocoa_init
        implicit none
        !plasma frequency for ion (emses, real unit)
        real(kind=8) :: wp_ion_emses,wp_ion_real
        !simulation volume
        real(kind=8) :: simu_vol
        !number of super particles
        integer sup_par_num,real_par_num_per_sup_par
        integer nprocess, ierr

        call MPI_Comm_size(CTCA_subcomm, nprocess, ierr)

        pbuf_size=size(pbuf)
        flag_size = size(flag)
        allocate(dist(pbuf_size))
        allocate(num_par(-energy_bin:energy_bin,spec_num))
        allocate(num_par_v(-energy_bin:energy_bin,v_dim,spec_num))
        call CTCAR_regarea_int(flag,flag_size,flag_id)
        call CTCAR_regarea_int(num_par,size(num_par),num_par_id)
        call CTCAR_regarea_int(num_par_v,size(num_par_v),num_par_v_id)

        ! set parameters from environment variables
        call get_environment_variable("GRID_LENGTH",env_buffer)
        read(env_buffer,*) grid_length
        call get_environment_variable("NEIGHBOUR_THR",env_buffer)
        read(env_buffer,*) neighbour_thr
        call get_environment_variable("CORRECT_BY_BIN_WIDTH",env_buffer)
        read(env_buffer,*) correct_by_bin_width
        call get_environment_variable("STEP_FROM",env_buffer)
        read(env_buffer,*) step_from
        call get_environment_variable("STEP_TO",env_buffer)
        read(env_buffer,*) step_to

        call get_environment_variable("SHIP_COORD", env_buffer)
        !format: "(x1,y1,z1)-(x2,y2,z2)"
        env_buffer = adjustl(env_buffer)  ! remove leading spaces
        call replace_chars(env_buffer, '(', ' ')
        call replace_chars(env_buffer, ')', ' ')
        call replace_chars(env_buffer, '-', ' ')
        call replace_chars(env_buffer, ',', ' ')
        read(env_buffer, *) ship_x_from, ship_y_from, ship_z_from, ship_x_to, ship_y_to, ship_z_to
        if (myid.eq.0) then
            print *, "Ship Coordinates:"
            print *, "From: (", ship_x_from, ",", ship_y_from, ",", ship_z_from, ")"
            print *, "To:   (", ship_x_to, ",", ship_y_to, ",", ship_z_to, ")"
        end if

        ! get plasma frequency for ion
        wp_ion_emses=wp(2)

        !plasma.inpでcv=10000としてもcv=10として実行されてしまう
        len_ratio=1/grid_length
        vel_ratio=cv*1e3/2.99792458d8
        vel_ratio_k=cv/2.99792458d8
        time_ratio=len_ratio/vel_ratio
        freq_ratio=1/time_ratio

        wp_ion_real=wp_ion_emses/freq_ratio
        ion_density=wp_ion_real**2*ion_mass*permittivity/ion_charge**2/1e6
        simu_vol=nx*ny*nz*grid_length**3

        sup_par_num=nodes(1)*nodes(2)*nodes(3)*pbuf_size
        real_par_num_per_sup_par=simu_vol*1d6*ion_density/sup_par_num
        sup_par_mass=ion_mass*real_par_num_per_sup_par
        neighbour_vol_real=4.0/3*pi*(neighbour_thr)**3
        if (myid.eq.0) then
            print *, "cv=",cv
            print *, "wp=",wp
            print *, "grid_length=",grid_length,"[m]"
            print *, "ion_charge=",ion_charge,"[C]"
            print *, "ion_mass=",ion_mass,"[kg]"
            print *, "electron_mass=",electron_mass,"[kg]"
            print *, "len_ratio=",len_ratio
            print *, "vel_ratio=",vel_ratio
            print *, "time_ratio=",time_ratio
            print *, "freq_ratio=",freq_ratio
            print *, "wp_ion_emses=",wp_ion_emses
            print *, "wp_ion_real=",wp_ion_real,"[Hz]"
            print *, "ion_density=",ion_density,"[/cc]"
            print *, "simu_vol=",simu_vol,"[m^3]"
            print *, "sup_par_num=",sup_par_num
            print *, "real_par_num_per_sup_par=",real_par_num_per_sup_par,"[1]"
            print *, "sup_par_mass=",sup_par_mass,"[kg]"
        end if
        !convert ship position to emses unit
        ship_x_from=ship_x_from*len_ratio
        ship_x_to=ship_x_to*len_ratio
        ship_y_from=ship_y_from*len_ratio
        ship_y_to=ship_y_to*len_ratio
        ship_z_from=ship_z_from*len_ratio
        ship_z_to=ship_z_to*len_ratio
        neighbour_thr=neighbour_thr*len_ratio
        !send request
        if (myid.eq.0) then
            req_params(1)=pbuf_size
            req_params(2)=-energy_bin
            req_params(3)=energy_bin
            req_params(4)=spec_num
            req_params(5)=nstep
            req_params(6)=real_par_num_per_sup_par
            req_params(7)=v_dim
            req_params(8)=nprocess
            req_params(9)=step_from
            req_params(10)=step_to
            req_params_real(1)=time_ratio
            req_params_real(2)=neighbour_vol_real
            req_params_real(3)=grid_length
            call CTCAR_sendreq_withreal8(req_params,size(req_params),req_params_real,size(req_params_real))
        end if
        flag(1)=1

    end subroutine cotocoa_init

    subroutine cotocoa_mainstep
        implicit none
        logical status1
        integer status2(MPI_STATUS_SIZE),ierr
        integer energy_index,j,species,pbuf_valid_size,v_sum_size(spec_num)
        real(kind=8) :: v(v_dim),energy,v_sum(spec_num,3),v_avg(spec_num,3)
        
        num_par(:,:)=0
        num_par_v(:,:,:)=0
        v_sum(:,:)=0
        v_sum_size(:)=0

        !set position of satellite
        shipx=ship_x_from+(ship_x_to-ship_x_from)*float(istep-step_from)/float(step_to-step_from)
        shipy=ship_y_from+(ship_y_to-ship_y_from)*float(istep-step_from)/float(step_to-step_from)
        shipz=ship_z_from+(ship_z_to-ship_z_from)*float(istep-step_from)/float(step_to-step_from)
        if (myid.eq.0.and.step_from.lt.istep.and.istep.lt.step_to) then
            print *,"ship position=",shipx,shipy,shipz
        end if

        dist(:)=sqrt((pbuf(:)%x-shipx)**2+(pbuf(:)%y-shipy)**2+(pbuf(:)%z-shipz)**2)
        
        !check flag
        do while (.true.)
            if (flag(1).eq.0) then
                !wait for worker
                !print *, "wait for worker",myid
                call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status1, status2, ierr)
            else
                exit
            end if
        end do

        !calculate energy
        pbuf_valid_size=sum(totalp)
        energy_size=0
        do i=1, pbuf_valid_size
            if (dist(i).lt.neighbour_thr.and.step_from.lt.istep.and.istep.lt.step_to) then
                energy_size=energy_size+1
                species=pbuf(i)%spec
                !energy(eV)
                if (species.eq.1) then
                    energy=electron_mass*(pbuf(i)%vx**2+pbuf(i)%vy**2+pbuf(i)%vz**2)/(vel_ratio_k**2)/2/ion_charge
                else
                    energy=ion_mass*(pbuf(i)%vx**2+pbuf(i)%vy**2+pbuf(i)%vz**2)/(vel_ratio_k**2)/2/ion_charge
                end if
                !velocity(m/s)
                v(4)=pbuf(i)%vx
                v(5)=pbuf(i)%vy
                v(6)=pbuf(i)%vz
                v(1:3)=abs(v(4:6))
                v(7:9)=-v(4:6)
                v(:)=v(:)/vel_ratio_k
                v_sum(species,:)=v_sum(species,:)+v(1:3)
                v_sum_size(species)=v_sum_size(species)+1
                !check energy
                if (energy.gt.0) then
                    energy_index=int(10*log10(energy))
                    !check boundary
                    energy_index=max(energy_index,lbound(num_par, 1))
                    energy_index=min(energy_index,ubound(num_par, 1))
                    !switch to change unit of histogram
                    if (correct_by_bin_width.eq.1) then
                        ![1/m^3/eV]
                        num_par(energy_index,species)=num_par(energy_index,species)+1
                    else
                        ![eV/m^3/eV]
                        num_par(energy_index,species)=num_par(energy_index,species)+energy
                    end if
                end if
                do j=1,v_dim
                    if (v(j).gt.0) then
                        energy_index=int(10*log10(v(j)))
                        !check boundary
                        energy_index=max(energy_index,lbound(num_par_v, 1))
                        energy_index=min(energy_index,ubound(num_par_v, 1))
                        num_par_v(energy_index,j,species)=num_par_v(energy_index,j,species)+1
                    end if
                end do
            end if
        end do
        
        !set flag
        flag(1)=0
        flag(2)=myid
        flag(3)=pbuf_size
        flag(4)=pbuf_mem
        flag(5)=istep
        flag(6)=energy_size
        !ランクが保有する粒子のうち衛星に近いものの速度成分の平均
        if (v_sum_size(1).gt.0.and.v_sum_size(2).gt.0) then
            v_avg(1,:)=v_sum(1,:)/v_sum_size(1)
            v_avg(2,:)=v_sum(2,:)/v_sum_size(2)
            !print *,"step",istep,"rank",myid,"avg_v=",v_avg
        end if
    end subroutine cotocoa_mainstep

    subroutine cotocoa_finalize
    integer date_time(8)
    character(len=100) :: date_str(3)
    flag(1)=2
    if (myid.eq.0) print *, "requester finished:",myid
    end subroutine cotocoa_finalize

    subroutine replace_chars(str, old, new)
        character(len=*), intent(inout) :: str
        character(len=1), intent(in) :: old, new
        integer :: i
        do i = 1, len_trim(str)
            if (str(i:i) == old) str(i:i) = new
        end do
    end subroutine replace_chars

end module m_ctcamain
