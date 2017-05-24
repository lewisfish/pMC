MODULE constants
!
! Module containing constants:
!         PI,TWOPI, the number of grid elements in each direction n%g,
!         and the vars for the filepaths.

    implicit none

    integer, parameter :: nxg=200, nyg=200, nzg=200
    real,    parameter :: PI=4.*atan(1.), TWOPI=2.*PI
    real               :: xmax, ymax, zmax
    character(len=10)  :: beam
    character(len=255) :: cwd, homedir, fileplace, resdir

end MODULE constants
