MODULE constants
!
! Module containing constants:
!         PI,TWOPI, the number of grid elements in each direction n%g,
!         and the vars for the filepaths.

    implicit none

    integer, parameter :: nxg=500, nyg=500, nzg=500, nbins=1
    real,    parameter :: PI=4.*atan(1.), TWOPI=2.*PI
    real               :: xmax, ymax, zmax, costim, sintim, sinpim, cospim, v(3), phiim, thetaim, bin_wid
    character(len=10)  :: beam
    character(len=255) :: cwd, homedir, fileplace, resdir

end MODULE constants
