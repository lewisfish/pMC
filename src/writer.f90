module writer_mod

    implicit none

    private
    public :: writer

    interface write_binary
        module procedure write_binaryR2
        module procedure write_binaryR3
        module procedure write_binaryR4
        module procedure write_binaryI3
    end interface

    contains
        subroutine writer(nphotons, numproc)

            ! use iarray,          only : phasorGLOBAL
            use utils,           only : str
            use iso_fortran_env, only : int64

            implicit none

            integer,             intent(IN) :: numproc
            integer(kind=int64), intent(IN) :: nphotons

            character(len=256)  :: filename


            ! filename = 'jmean/testp.raw'
            
            ! call write_binary(trim(filename), real(phasorGLOBAL))
            ! print*,trim(filename)

            ! filename = 'jmean/test.raw'
            
            ! call write_binary(trim(filename), abs(phasorGLOBAL)**2)
            ! print*,trim(filename)

        end subroutine writer


        subroutine write_binaryR2(filename, array)

            use constants, only : fileplace

            implicit none

            character(len=*), intent(IN) :: filename
            real,             intent(IN) :: array(:,:)
            
            integer(kind=8) :: u, i

            inquire(iolength=i)array
            open(newunit=u,file=trim(fileplace)//trim(filename),access='stream',status='REPLACE',form='unformatted')
            write(u) array
            close(u)

        end subroutine write_binaryR2


        subroutine write_binaryR3(filename, array)

            use constants, only : fileplace

            implicit none

            character(len=*), intent(IN) :: filename
            real,             intent(IN) :: array(:,:,:)
            
            integer(kind=8) :: u, i

            inquire(iolength=i)array
            open(newunit=u,file=trim(fileplace)//trim(filename),access='stream',status='REPLACE',form='unformatted')
            write(u) array
            close(u)

        end subroutine write_binaryR3


        subroutine write_binaryR4(filename, array)

            use constants, only : fileplace

            implicit none

            character(len=*), intent(IN) :: filename
            real,             intent(IN) :: array(:,:,:,:)
            
            integer :: u, i

            inquire(iolength=i)array
            open(newunit=u,file=trim(fileplace)//trim(filename),access='direct',status='REPLACE',form='unformatted',&
            recl=i)
            write(u,rec=1) array
            close(u)

        end subroutine write_binaryR4


        subroutine write_binaryI3(filename, array)

            use constants, only : fileplace

            implicit none

            character(len=*), intent(IN) :: filename
            integer,          intent(IN) :: array(:,:,:)
            
            integer :: u, i

            inquire(iolength=i)array
            open(newunit=u,file=trim(fileplace)//trim(filename),access='direct',status='REPLACE',form='unformatted',&
            recl=i)
            write(u,rec=1) array
            close(u)

        end subroutine write_binaryI3

end module writer_mod
