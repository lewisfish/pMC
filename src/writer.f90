MODULE writer_mod

implicit none
save

private
public :: writer

interface write_binary
    module procedure write_binaryR2
    module procedure write_binaryR3
    module procedure write_binaryR4
    module procedure write_binaryI3
end interface

CONTAINS
    subroutine writer(nphotons, numproc, iseed)

        use constants,       only : nxg, nyg, nzg, xmax, ymax, zmax, beam, nbins
        use iarray,          only : phasorGLOBAL!, jmeanGLOBAL!,, imageGLOBAL, imagetGLOBAL, imagepGLOBAL,  
        use iso_fortran_env, only : int64
        use utils,           only : str

        implicit none

        integer(kind=int64), intent(IN) :: nphotons
        integer, intent(IN) :: iseed
        integer :: numproc
        integer :: u, i
        character(len=256) :: filename


        filename = 'jmean/'//trim(beam)//'-tester-phase'//str(iseed)//'.raw'
        
        call write_binary(trim(filename), phasorGLOBAL)
        print*,trim(filename)

        filename = 'jmean/'//trim(beam)//'-tester-intesity'//str(iseed)//'.raw'
        phasorGLOBAL = phasorGLOBAL **2
        call write_binary(trim(filename), phasorGLOBAL)
        print*,trim(filename)


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

end MODULE writer_mod
