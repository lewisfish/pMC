module subs

implicit none

save

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
            resdir=trim(homedir)//'res/'

        end subroutine directory


        subroutine zarray

            use iarray

            !sets all arrays to zero
            implicit none

            xface = 0.
            yface = 0.
            zface = 0.
            ! rhokap = 0.
            ! jmean = 0.
            ! intensity = 0.
            ! image = 0.
            ! imaget = 0.
            ! imagep = 0.
            ! imagethg = 0.
            phasor = 0.

        end subroutine zarray


        subroutine alloc_array(numproc, id)
        !  subroutine allocates allocatable arrays
        !   
        !   
            use iso_fortran_env, only : int64
            use utils,           only : str, mem_free
            use constants,       only : nxg, nyg, nzg, nbins
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

            ! allocate(rhokap(nxg,nyg,nzg))
            ! inquire(iolength=i)rhokap
            ! call chck_mem(cnt, i, limit, 'rhokap', numproc)

            ! allocate(intensity(nxg,nyg,nzg))
            ! inquire(iolength=i)intensity
            ! call chck_mem(cnt, i, limit, 'intensity', numproc)

            ! allocate(image(-((Nbins-1)/2):((Nbins-1)/2), -((Nbins-1)/2):((Nbins-1)/2)))
            ! inquire(iolength=i)image
            ! call chck_mem(cnt, i, limit, 'image', numproc)

            ! allocate(imaget(-((Nbins-1)/2):((Nbins-1)/2), -((Nbins-1)/2):((Nbins-1)/2)))
            ! inquire(iolength=i)imaget
            ! call chck_mem(cnt, i, limit, 'imaget', numproc)

            ! allocate(imagethg(-((Nbins-1)/2):((Nbins-1)/2), -((Nbins-1)/2):((Nbins-1)/2)))
            ! inquire(iolength=i)imagethg
            ! call chck_mem(cnt, i, limit, 'imagethg', numproc)

            !       allocate(imagep(-((Nbins-1)/2):((Nbins-1)/2), -((Nbins-1)/2):((Nbins-1)/2)))
            ! inquire(iolength=i)imagep
            ! call chck_mem(cnt, i, limit, 'imagep', numproc)

            allocate(phasor(nbins,nbins))
            inquire(iolength=i)phasor
            call chck_mem(cnt, i, limit, 'phasor', numproc)

            ! allocate(jmean(nxg,nyg,nzg))
            ! inquire(iolength=i)jmean
            ! call chck_mem(cnt, i, limit, 'jmean', numproc)

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