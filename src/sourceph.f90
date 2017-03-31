MODULE sourceph_mod

implicit none
save

CONTAINS
   subroutine sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed)

   use constants, only : nxg,nyg,nzg,twopi
   use photon_vars

   implicit none


   integer, intent(OUT)   :: xcell, ycell, zcell
   integer, intent(INOUT) :: iseed
   real,    intent(IN)    :: xmax, ymax, zmax
   real                   :: ran2, theta

   zp = zmax - (1.e-5*(2.*zmax/nzg))
   xp = 0.
   if(ran2(iseed) < 0.5)then
      yp = 0.
   else
      yp = 0.0000625
   end if
   
   phi = twopi * ran2(iseed)
   cosp = cos(phi)
   sinp = sin(phi)
   theta = acos(ran2(iseed) - 1.)        
   sint = sin(theta)
   cost = cos(theta)
 
   nxp = sint * cosp  
   nyp = sint * sinp
   nzp = cost

   phase = 0.

   !*************** Linear Grid *************************
   xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
   ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
   zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
   !*****************************************************
   end subroutine sourceph
end MODULE sourceph_mod
