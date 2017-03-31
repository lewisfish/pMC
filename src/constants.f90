MODULE constants
!
! Module containing constants:
!         PI,TWOPI, the number of grid elements in each direction n%g,
!         the various bin # params,
!         and the vars for the filepaths.
!

implicit none
save

integer, parameter :: nxg=250, nyg=250, nzg=250
real,    parameter :: PI = 3.141592, TWOPI = 6.283185
character(len=255) :: cwd, homedir, fileplace, resdir

end MODULE constants
