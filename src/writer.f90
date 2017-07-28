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

        use constants, only : nxg, nyg, nzg, xmax, ymax, zmax, beam, nbins
        use iarray,    only : phasorGLOBAL!, jmeanGLOBAL!,, imageGLOBAL, imagetGLOBAL, imagepGLOBAL,  

        implicit none

        integer, intent(IN) :: nphotons, numproc
        integer :: u, i, j
        character(len=256) :: filename

        ! jmeanGLOBAL = jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

        print*,' '

        ! print*,'Written to:'
        ! filename = 'jmean/'//trim(beam)//'-test-jmean.dat'
        ! print*,trim(filename)
        ! call write_binary(trim(filename), jmeanGLOBAL)


        filename = 'jmean/'//trim(beam)//'-test-phase.raw'
        call write_binary(trim(filename), phasorGLOBAL)
        print*,trim(filename)

        filename = 'jmean/'//trim(beam)//'-test-intesity.raw'
        call write_binary(trim(filename), phasorGLOBAL**2.)
        print*,trim(filename)

        ! imageGLOBAL = imageGLOBAL**2

        ! filename = 'im/image-'//trim(beam)//'-test.dat'
        ! open(newunit=u,file=trim(filename),status='replace')
        ! do i = -((Nbins-1)/2), ((Nbins-1)/2)
        !     write(u,*) (imageGLOBAL(j,i),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        ! end do
        ! close(u)
        ! print*,trim(filename)

        ! filename = 'im/image-'//trim(beam)//'-tau.dat'
        ! open(newunit=u,file='im/image-'//trim(beam)//'-tau.dat',status='replace')
        ! do i = -((Nbins-1)/2), ((Nbins-1)/2)
        !     write(u,*) (imagetGLOBAL(j,i),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        ! end do
        ! close(u)
        ! print*,trim(filename)

        ! filename = 'im/image-'//trim(beam)//'-tau-hg.dat'
        ! open(newunit=u,file='im/image-'//trim(beam)//'-tau-hg.dat',status='replace')
        ! do i = -((Nbins-1)/2), ((Nbins-1)/2)
        !     write(u,*) (imagethgGLOBAL(j,i),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        ! end do
        ! close(u)
        ! print*,trim(filename)
        ! imagepGLOBAL = imagepGLOBAL**2

        ! filename = 'im/image-'//trim(beam)//'-test-phase.dat'
        ! open(newunit=u,file=trim(filename),status='replace')
        ! do i = -((Nbins-1)/2), ((Nbins-1)/2)
        !     write(u,*) (imagepGLOBAL(j,i),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        ! end do
        ! close(u)
        ! print*,trim(filename)


        ! filename = 'jmean/'//trim(beam)//'-test-intesity.dat'
        ! phasorGLOBAL = phasorGLOBAL**2
        ! call write_binary(trim(filename), phasorGLOBAL)
        ! print*,trim(filename)

        ! filename = 'jmean/'//trim(beam)//'-cur-absenergy.dat'
        ! print*,trim(filename)
        ! tmp = jmeanGLOBAL*mua
        ! call write_binary(trim(filename), tmp)

    end subroutine writer


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
