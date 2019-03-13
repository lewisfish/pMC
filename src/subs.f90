module subs

implicit none

    contains

        subroutine directory
        !  subroutine defines vars to hold paths to various folders   
        !   
        !   
            use constants, only : cwd, homedir, fileplace, resdir

            implicit none

            !get current working directory

            call get_environment_variable('PWD', cwd)

            ! get 'home' dir from cwd
            homedir = trim(cwd(1:len(trim(cwd))-3))
            ! get data dir
            fileplace = trim(homedir)//'data/'
            ! get res dir
            resdir = trim(homedir)//'res/'

        end subroutine directory


        subroutine zarray

            use iarray

            !sets all arrays to zero
            implicit none
            rhokap = 0.
            xface = 0.
            yface = 0.
            zface = 0.
            phasor= cmplx(0.,0.)
            phasorGLOBAL = cmplx(0.,0.)
            imageb = cmplx(0.d0, 0.d0)
            imagebGLOBAL = cmplx(0.d0, 0.d0)

        end subroutine zarray


        subroutine alloc_array(numproc, id)
        !  subroutine allocates allocatable arrays
        !   
        !   
            use iso_fortran_env, only : int64
            use utils,           only : str, mem_free
            use constants,       only : nxg, nyg, nzg
            use iarray

            implicit none

            integer, intent(IN) :: numproc, id

            integer(int64) :: limit, cnt, i


            limit = mem_free()
            cnt = 0_int64

            allocate(xface(nxg+1))
            inquire(iolength=i)xface(:)
            call chck_mem(cnt, i, limit, 'xface', numproc)

            allocate(yface(nyg+1))
            inquire(iolength=i)yface
            call chck_mem(cnt, i, limit, 'yface', numproc)

            allocate(zface(nzg+1))
            inquire(iolength=i)zface
            call chck_mem(cnt, i, limit, 'zface', numproc)

            allocate(imageb(-50:50, -50:50))
            allocate(imagebGLOBAL(-50:50, -50:50))

            allocate(rhokap(nxg, nyg, nzg))

            allocate(phasor(nxg,nyg,nzg))
            allocate(phasorGLOBAL(nxg,nyg,nzg))


            if(id == 0)print'(A,1X,F5.2,A)','allocated:',dble(cnt)/dble(limit)*100.d0,' % of total RAM'
        end subroutine alloc_array


        subroutine chck_mem(cur, new, limit, name, numproc)
        !routine to check if the system has enough RAM available in order to run the simulation
        !cur: current memory assigned, new: new memory to be assigned
        !limit: the limit of RAM available, name: name of array to be assigned, numproc: processor #

            use iso_fortran_env, only : int64
            use utils,           only : str

            implicit none

            integer(int64), intent(IN)    :: new, limit
            integer(int64), intent(INOUT) :: cur 
            integer,        intent(IN)    :: numproc
            character(*),   intent(IN)    :: name

            integer :: error

            cur = cur + new * numproc
            if(cur > limit)then
                print*,'Need '//str(cur-limit)//' more memory to run. '//name
                call mpi_finalize(error)
                stop
            end if
        end subroutine chck_mem
end MODULE subs