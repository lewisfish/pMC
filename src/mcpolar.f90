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
    integer :: iseed, xcell, ycell, zcell, holdseed, u, idx, idy
    logical :: tflag
    real    :: nscatt, raxi, dtoskin, binwid, ran2
    real    :: delta, start, finish, start2, finish2, tmp, phiim, thetaim

    real    :: pixres

    !mpi variables
    integer :: id, error, numproc
    real    :: nscattGLOBAL

    call MPI_init(error)
    call MPI_Comm_size(MPI_COMM_WORLD, numproc, error)
    call MPI_Comm_rank(MPI_COMM_WORLD, id, error)

    !set directory paths
    call directory()

    !allocate and set arrays to 0
    call alloc_array(numproc, id)
    call zarray()


    !**** Read in parameters from the file input.params
    open(newunit=u,file=trim(resdir)//'input.params',status='old')
    read(u,*) nphotons
    read(u,*) xmax
    read(u,*) ymax
    read(u,*) zmax
    read(u,*) n1
    read(u,*) n2
    read(u,*) beam
    read(u,*) raxi
    read(u,*) n
    read(u,*) dtoskin
    read(u,*) l
    read(u,*) imgsize
    read(u,*) pixres

    close(u)

    iseed = -123456879 + id
    iseed = -abs(iseed+id)  ! Random number seed must be negative for ran2

    holdseed = iseed

    phiim   = 0. * pi/180.
    thetaim = 180. * pi/180.


    !image postion vector
    !angle for vector
    costim = cos(thetaim)
    sintim = sin(thetaim)
    sinpim = sin(phiim)
    cospim = cos(phiim)

    !vector
    v(1) = sintim * cospim     
    v(2) = sintim * sinpim
    v(3) = costim      

    call init_opt4
    wavelength = 780.d-9!587.6d-9!7.6e-9
    fact = twopi/wavelength
    tana=tan(5.d0*pi/180.d0)

    binwid = 2.d-3/1000.!pixres
    pixels = imgsize/binwid
    allocate(imageb(1000, 1000))
    allocate(imagebGLOBAL(1000, 1000))

    !2.*xmax / real(nxg)
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
    print*,'Photons now running on core: ',colour(id, str(30+mod(id,7)), bold)
    do j = 1 , nphotons

        call init_opt4
        phase = 0.d0
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

        if(mod(j,1000000_int64) == 0)then
            print *, colour(j, blue, bold),' scattered photons completed on core: ',colour(id, str(30+mod(id,7)), bold)
        end if

        !***** Release photon from point source *******************************
        call sourceph(xcell,ycell,zcell,raxi, dtoskin, iseed)



        ! imaget(xcell, ycell) = imaget(xcell, ycell) + cmplx(cos((phase * fact)), sin(phase * fact))


        !****** Find scattering location
        call tauint1(xcell,ycell,zcell,tflag,iseed,delta)

        !******** Photon scatters in grid until it exits (tflag=TRUE) 
        do while(tflag.eqv..FALSE.)

            ! if(ran2(iseed) < albedo)then!interacts with tissue
            !       call stokes(iseed)
            !       nscatt = nscatt + 1        
            !    else
                  tflag = .true.
                  exit
            ! end if

            ! !************ Find next scattering location

            ! call tauint1(xcell,ycell,zcell,tflag,iseed,delta)
            ! if(.not. tflag)call peeling(xcell,ycell,zcell,delta)

        end do
        ! stop
        if(xcell /= -1 .and. ycell /= -1 .and. tflag)then
            ! if(abs(xp) <= imgsize/2. .and. abs(yp) <= imgsize/2.)then
                ! print*,xp+xmax,yp+xmax,xcell,ycell
                ! idx = floor(xcell/binwid) + 1
                ! idy = floor(((yp + ymax) - (2.*xmax-imgsize)/2.)/binwid) + 1
                ! print*,idx,idy
                ! print*,xcell,ycell
                ! imageb(xcell, ycell) = imageb(xcell, ycell) + cmplx(cos((phase * fact)), sin(phase * fact))
            ! end if
        end if

    end do      ! end loop over nph photons

    ! call mpi_reduce(imaget, imagetGLOBAL, size(imaget), mpi_double_complex, mpi_sum, 0, mpi_comm_world, error)
    ! call mpi_reduce(imageb, imagebGLOBAL, size(imageb), mpi_double_complex, mpi_sum, 0, mpi_comm_world, error)
    call mpi_reduce(phasor, phasorGLOBAL, size(phasor), mpi_double_complex, mpi_sum, 0, mpi_comm_world, error)

    call mpi_reduce(nscatt, nscattGLOBAL, 1, mpi_double, mpi_sum, 0, mpi_comm_world, error)

    call cpu_time(finish)
    if(finish-start.ge.60.)then
        print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
    else
        print*, 'time taken ~ ',colour(floor(finish-start/60.),red, bold),'s'
    end if


    if(id == 0)then
        print*,'Average # of scatters per photon:',nscattGLOBAL/(nphotons*numproc)
        !write out files

        ! open(newunit=u,file="bessel-l"//str(int(l))//"-top-int-test.dat",access="stream",form="unformatted",status="replace")
        ! write(u)abs(imagetGLOBAL)**2
        ! close(u)

        ! open(newunit=u,file="bessel-l"//str(int(l))//"-bot-int-test.dat",access="stream",form="unformatted",status="replace")
        ! write(u)abs(imagebGLOBAL)**2
        ! close(u)

        ! open(newunit=u,file="bessel-l"//str(int(l))//"-phase.dat",access="stream",form="unformatted",status="replace")
        ! write(u)real(imageGLOBAL)
        ! close(u)

        call writer(nphotons, numproc)
        print*,colour('write done',black,white_b,bold)
    end if

    call MPI_Finalize(error)
end program mcpolar