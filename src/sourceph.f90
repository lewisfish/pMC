MODULE sourceph_mod

implicit none
save

CONTAINS
   subroutine sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed)

   use constants, only : nxg,nyg,nzg,twopi,pi
   use photon_vars

   implicit none


   integer, intent(OUT)   :: xcell, ycell, zcell
   integer, intent(INOUT) :: iseed
   real,    intent(IN)    :: xmax, ymax, zmax
   real                   :: ran2, theta, angle, R_axicon, r_pos, n, wavelength

   zp = zmax - (1.e-5*(2.*zmax/nzg))
   xp = 0.
   yp = 0.
   do
      xp = -.00005+ran2(iseed)*0.0001
      yp = -.00005+ran2(iseed)*0.0001
      if(xp**2. + yp **2. <= .00005**2)exit
   end do
      angle = 5.

theta = (180-angle)*pi/180.
cost = cos(theta)
sint = sin(theta)
   phi = atan2(xp,yp)
   cosp = cos(phi)
   sinp = sin(phi)


   ! phi = twopi * ran2(iseed)
   ! cosp = cos(phi)
   ! sinp = sin(phi)
   ! theta = acos(ran2(iseed) - 1.)        
   ! sint = sin(theta)
   ! cost = cos(theta)
 
   nxp = sint * cosp  
   nyp = sint * sinp
   nzp = cost

   n = 1.55
   wavelength = 1435.e-9
   R_axicon = 25.4
   r_pos = sqrt(xp**2 + yp**2)

   ! phase = (twopi * n / wavelength) * (R_axicon - r_pos)*tan(angle*pi/180.)
   ! phase = 0.
   phase = -twopi*r_pos/wavelength * sin(angle*pi/180.)
   !*************** Linear Grid *************************
   xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
   ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
   zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
   !*****************************************************
   end subroutine sourceph
end MODULE sourceph_mod
