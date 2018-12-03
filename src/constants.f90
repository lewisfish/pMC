module constants
!
! module containing constants:
!         PI,TWOPI, the number of grid elements in each direction n%g,
!         and the vars for the filepaths.

    implicit none

    integer, parameter :: nxg=100, nyg=100, nzg=250, nbins=100
    real,    parameter :: PI=4.d0*atan(1.d0), TWOPI=2.d0*PI
    real               :: xmax, ymax, zmax, fact, tana
    character(len=10)  :: beam
    character(len=255) :: cwd, homedir, fileplace, resdir

end module constants
