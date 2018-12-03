program mcpolar

    use mpi

    use iso_fortran_env, only : int64

    !shared data
    use utils, only : green, blue, bold, colour, red, white_b, black, str
    use constants
    use photon_vars
    use iarray
    use opt_prop

    !subroutines
    use subs
    use gridset_mod
    use sourceph_mod
    use inttau2
    use ch_opt
    use stokes_mod
    use writer_mod

    implicit none

    integer(kind=int64) :: nphotons, j
    integer :: iseed, xcell, ycell, zcell, holdseed
    logical :: tflag
    real    :: nscatt, raxi, dtoskin, binwid
    real    :: delta, start, finish, start2, finish2, tmp

    !mpi variables
    integer :: id, error, numproc
    real    :: nscattGLOBAL,ran2

    call MPI_init(error)

    call MPI_Comm_size(MPI_COMM_WORLD, numproc, error)

    call MPI_Comm_rank(MPI_COMM_WORLD, id, error)

    !set directory paths
    call directory()

    !allocate and set arrays to 0
    call alloc_array(numproc, id)
    call zarray()


    !**** Read in parameters from the file input.params
    open(10,file=trim(resdir)//'input.params',status='old')
    read(10,*) nphotons
    read(10,*) xmax
    read(10,*) ymax
    read(10,*) zmax
    read(10,*) n1
    read(10,*) n2
    read(10,*) beam
    read(10,*) raxi
    read(10,*) n
    read(10,*) dtoskin
    close(10)

    iseed = -876535443 + id
    iseed = -abs(iseed+id)  ! Random number seed must be negative for ran2

    holdseed = iseed

    call init_opt4
    wavelength = 1435.e-9
    fact = twopi/wavelength
    tana=tan(5.d0*pi/180.d0)

    binwid = 2.*xmax / real(nxg)
    if(id == 0)then
        print*, ''      
        print*,'# of photons to run',nphotons*int(numproc,kind=int64)
    end if

    !***** Set up density grid *******************************************
    call gridset(id)

    !***** Set small distance for use in optical depth integration routines 
    !***** for roundoff effects when crossing cell walls
    delta = 1.e-8*(2.*zmax/nzg)
    nscatt = 0
    nscattGLOBAL = 0

    call cpu_time(start)
    call cpu_time(start2)
    !loop over photons 
    call MPI_Barrier(MPI_COMM_WORLD, error)
    print*,'Photons now running on core: ',colour(id, green)
    do j = 1 , nphotons

        call init_opt4

        tflag=.FALSE.

        if(j == 100000 .and. id == 0)then
            call cpu_time(finish2)
            print*,' '
            tmp = (finish2-start2)/100000.*real(nphotons)
            if(tmp >= 60.)then
                tmp = tmp / 60.
                if(tmp > 60)then
                    tmp = tmp / 60.
                    print*,str(tmp),' hrs'
                else
                    print*,str(tmp),' mins'
                end if
            else
                print*,str(tmp),' s'
            end if
            print*,' '
        end if

        if(mod(j,10000000_int64) == 0)then
            print *, colour(j, blue, bold),' scattered photons completed on core: ',colour(id, str(30+mod(id,7)), bold)
        end if

        !***** Release photon from point source *******************************
        call sourceph(xcell,ycell,zcell,raxi, dtoskin, iseed)

        ! image(xcell, ycell) = image(xcell, ycell) + cmplx(cos((phase * fact)), sin(phase * fact))


        !****** Find scattering location

        call tauint1(xcell,ycell,zcell,tflag,iseed,delta)

        !******** Photon scatters in grid until it exits (tflag=TRUE) 
        do while(tflag.eqv..FALSE.)
            ! ran = ran2(iseed)

            if(ran2(iseed) < 0.)then!interacts with tissue
            !       ! call stokes(iseed)
            !       ! nscatt = nscatt + 1        
               else
                  tflag=.true.
                  exit
            end if

            !************ Find next scattering location

            call tauint1(xcell,ycell,zcell,tflag,iseed,delta)
            ! if(.not. tflag)call peeling(xcell,ycell,zcell,delta)

        end do

    end do      ! end loop over nph photons

    ! call mpi_reduce(image, imageGLOBAL, size(image), mpi_double_complex, mpi_sum, 0, mpi_comm_world, error)
    call mpi_reduce(phasor, phasorGLOBAL, size(phasor), mpi_double_complex, mpi_sum, 0, mpi_comm_world, error)

    call cpu_time(finish)
    if(finish-start.ge.60.)then
        print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
    else
        print*, 'time taken ~ ',colour(floor(finish-start/60.),red, bold),'s'
    end if


    if(id == 0)then
        print*,'Average # of scatters per photon:',nscattGLOBAL/(nphotons*numproc)
        !write out files

        call writer(nphotons, numproc)
        print*,colour('write done',black,white_b,bold)
    end if

    call MPI_Finalize(error)
end program mcpolar