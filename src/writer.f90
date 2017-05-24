MODULE writer_mod

implicit none
save

private
public :: writer

interface write_binary
    module procedure write_binaryR3
    module procedure write_binaryR4
    module procedure write_binaryI3
end interface

CONTAINS
    subroutine writer(nphotons, numproc)

        use constants, only : nxg, nyg, nzg, xmax, ymax, zmax, beam
        use iarray,    only : phasorGLOBAL, jmeanGLOBAL, intensityGLOBAL

        implicit none

        integer, intent(IN) :: nphotons, numproc
        character(len=256) :: filename

        jmeanGLOBAL = jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

        print*,

        print*,'Written to:'
        filename = 'jmean/'//trim(beam)//'small-jmean.dat'
        print*,trim(filename)
        call write_binary(trim(filename), jmeanGLOBAL)


        filename = 'jmean/'//trim(beam)//'small-phase.dat'
        print*,trim(filename)
        call write_binary(trim(filename), phasorGLOBAL)


        filename = 'jmean/'//trim(beam)//'small-intesity.dat'
        print*,trim(filename)
        call write_binary(trim(filename), intensityGLOBAL)

        ! filename = 'jmean/'//trim(beam)//'small-absenergy.dat'
        ! print*,trim(filename)
        ! tmp = jmeanGLOBAL*mua
        ! call write_binary(trim(filename), tmp)

    end subroutine writer


    subroutine write_binaryR3(filename, array)

        use constants, only : fileplace

        implicit none

        character(len=*), intent(IN) :: filename
        real,             intent(IN) :: array(:,:,:)
        
        integer :: u, i

        inquire(iolength=i)array
        open(newunit=u,file=trim(fileplace)//trim(filename),access='direct',status='REPLACE',form='unformatted',&
        recl=i)
        write(u,rec=1) array
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
