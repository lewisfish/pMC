MODULE sourceph_mod

implicit none
save

CONTAINS
   subroutine sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed)

   use constants, only : nxg,nyg,nzg,twopi,pi
   use opt_prop,  only : wavelength
   use photon_vars

   implicit none


   integer, intent(OUT)   :: xcell, ycell, zcell
   integer, intent(INOUT) :: iseed
   real,    intent(IN)    :: xmax, ymax, zmax
   real                   :: ran2, theta, R_axicon, r_pos, n, phigauss, r1, w

   zp = zmax - (1.e-5*(2.*zmax/nzg))
   xp = 0.

   w = 0.00002
   r1 = w * sqrt(-log(ran2(iseed)))
   phigauss = twopi * ran2(iseed)
   xp = r1 * cos(phigauss)
   yp = r1 * sin(phigauss)
   ! print*,xp, yp
   ! call exit(0)
   ! if(ran2(iseed) < 0.5)then
   !    ycell = 50
   ! else
   !    ycell = 150
   ! end if
   ! yp = (real(ycell -.5)/nyg * 2.*ymax )- ymax 
   ! do
   !    xp = -.00005+ran2(iseed)*0.001
   !    yp = -.00005+ran2(iseed)*0.001
   !    if(xp**2. + yp **2. <= .00005**2)exit
   ! end do
      
   angle = 170.

   theta = angle*pi/180.
   cost = cos(theta)
   sint = sin(theta)


! < < is top left
! > > is bottom right
! > < is bottom left
! < > is top right


if(xp >0. .and. yp >0.)then
   phi = atan(xp/yp)   
elseif(xp >0. .and. yp <0.)then
   phi = pi/2. + atan(yp/(-xp))
elseif(xp <0. .and. yp <0.)then
   phi = pi + atan(xp/yp)
elseif(xp<0. .and. yp > 0.)then
   phi = 3.*pi/2. + atan((-yp)/xp)
end if

   phi = atan2(xp,yp)

   cosp = cos(phi)
   sinp = sin(phi)


   ! phi = twopi * ran2(iseed)
   ! cosp = cos(phi)
   ! sinp = sin(phi)

   nxp = sint * cosp  
   nyp = sint * sinp
   nzp = cost

   n = 1.55
   R_axicon = .00254
   r_pos = sqrt(xp**2 + yp**2)

   phase = (twopi * n / wavelength) * (R_axicon - r_pos)*tan(angle*pi/180.)
   ! phase = -twopi*r_pos/wavelength * sin(angle*pi/180.)
   !*************** Linear Grid *************************
   xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
   ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
   zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
   !*****************************************************
   end subroutine sourceph
end MODULE sourceph_mod
