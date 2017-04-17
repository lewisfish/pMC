MODULE constants
!
! Module containing constants:
!         PI,TWOPI, the number of grid elements in each direction n%g,
!         and the vars for the filepaths.

implicit none
save

integer, parameter :: nxg=200, nyg=200, nzg=500
real,    parameter :: PI=3.141592, TWOPI=6.283185
real               :: xmax, ymax, zmax
character(len=10)  :: beam
character(len=255) :: cwd, homedir, fileplace, resdir

end MODULE constants
