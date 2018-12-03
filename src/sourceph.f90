MODULE sourceph_mod

    implicit none
    save

    private
    public :: sourceph

    CONTAINS
        subroutine sourceph(xcell, ycell, zcell,raxi, dtoskin, iseed)

            use constants,   only : nxg, nyg, nzg, xmax, ymax, zmax
            use photon_vars, only : xp, yp, zp

            implicit none

            integer,          intent(OUT)   :: xcell, ycell, zcell
            integer,          intent(INOUT) :: iseed
            real, intent(IN) :: raxi, dtoskin

            call bessel(raxi, dtoskin, iseed)
           

            xcell = int(nxg * (xp + xmax) / (2. * xmax)) + 1
            ycell = int(nyg * (yp + ymax) / (2. * ymax)) + 1
            zcell = int(nzg * (zp + zmax) / (2. * zmax)) + 1
        end subroutine sourceph


        subroutine bessel(Raxi, d, iseed)

            use constants,   only : nzg, xmax, zmax, ymax,pi
            use opt_prop,    only : n
            use photon_vars, only : xp, yp, zp, sint, cost, sinp, cosp, phi, phase, nxp, nyp, nzp

            implicit none

            integer, intent(INOUT) :: iseed
            real, intent(IN) :: raxi, d

            real :: r_pos, tana, x0, y0, z0, dist

            ! axi_thickness = .11d0 ! 1.1mm
            ! Raxi = 2.54d0/2.d0 !25.4/2mm
            ! d = 1.11d0! 11.1mm

            phase = 0.d0
            tana = tan(5.d0*pi/180.d0)
            !top of axicon
            call rang(xp, yp, 0.d0, 1.d-3/4.d0, iseed)
            ! y = qrang(0.d0, 1.d-3/4.d0, iseed)

            r_pos = sqrt(xp**2 + yp**2)
            zp = (raxi - r_pos) * tana

            phase = n*zp

            x0 = ranu(-xmax, xmax, iseed)!2.*sqrt(2.)*xmax * ran*sin(37.*pi/180.)!ranu(-xmax, xmax, iseed)
            y0 = ranu(-ymax, ymax, iseed)!2.*sqrt(2.)*xmax * ran*cos(37.*pi/180.)!ranu(-ymax, ymax, iseed)

            z0 = (r_pos* tana) + zp + d !11.1mm

            dist = sqrt((x0 - xp)**2 + (y0 - yp)**2 + (z0 - zp)**2)

            phase = phase + dist

            nxp = (x0 - xp) / dist
            nyp = (y0 - yp) / dist
            nzp = -(z0 - zp) / dist !-ive due to way z pos is defined

            cost = nzp
            sint = sqrt(1.d0 - cost**2)

            phi = atan2(nyp, nxp)

            cosp = cos(phi)
            sinp = sin(phi)

            xp = x0
            yp = y0
            zp = zmax - (1.e-5*(2.*zmax/nzg))

         end subroutine bessel


        ! subroutine gaussian_phase(iseed)
        ! ! f => focal length of lens
        ! ! zr => rayleigh length
        ! ! zf => z location of of the focus


        !     use constants,   only : pi, nzg, zmax, twopi, xmax
        !     use photon_vars, only : nxp, nyp, nzp, xp,yp,zp, phase,zr
        !     use opt_prop,    only : wavelength

        !     implicit none

        !     integer, intent(INOUT) :: iseed
        !     real :: ran2,w
        !     real :: xo, yo, zo, zf, f, fact, rz,  D, wo, r1, phigauss, r_pos


        !     w = 0.0000000015
        !     zf = zmax-.01
        !     zo = zmax - (1.e-5*(2.*zmax/nzg))
        !     xo = 0.

        !     f = zmax/2000.
        !     D = 6.e-8
        !     wo = (2./pi) * (wavelength*f)/D
        !     zr = (pi * wo**2)/wavelength
        !     ! print*,wo,zr
        !     !gaussian beam via box-muller method
        !     do
        !         r1 = w*sqrt(-2.*log(ran2(iseed)))
        !         phigauss = twopi * ran2(iseed)
        !         xo = r1 * cos(phigauss)
        !         yo = r1 * sin(phigauss)
        !         if(xo**2 + yo**2 < xmax**2.)exit
        !     end do

        !     xp = -xo*((zo - zf)/sqrt(xo**2 + yo**2 + f**2)) 
        !     yp = -yo*((zo - zf)/sqrt(xo**2 + yo**2 + f**2)) 
        !     zp = zf + f*((zo - zf)/sqrt(xo**2 + yo**2 + f**2)) 

        !     rz = -(zp - zf) * (1 + (zr / (zp - zf))**2)

        !     fact = 1./sqrt(1 + ((xp**2 + yp**2) / rz**2))

        !     nxp = fact * (xp/rz)
        !     nyp = fact * (yp/rz)
        !     nzp = -1. * fact

        !     r_pos = sqrt(xp**2 + yp**2)
        !     phase = (twopi*(abs(zp-zo)+r_pos)/wavelength)
        !     zp = zo

        ! end subroutine gaussian_phase


        ! subroutine circular(radius, iseed)

        !     use constants, only : nzg, zmax
        !     use photon_vars, only : xp, yp, zp, sint, cost, sinp, cosp, phi, phase

        !     implicit none

        !     real,    intent(IN)    :: radius
        !     integer, intent(INOUT) :: iseed

        !     real :: ran2

        !     zp = zmax - (1.e-5*(2.*zmax/nzg))

        !     do
        !        xp = 2. * radius * ran2(iseed) - radius
        !        yp = 2. * radius * ran2(iseed) - radius
        !        if(xp**2 + yp **2 <= radius**2)exit
        !     end do

        !     cost = -1.
        !     sint = 1. - cost**2.

        !     phi = 0.
        !     cosp = cos(phi)
        !     sinp = sin(phi)

        !     phase = 0.

        ! end subroutine circular


        ! subroutine point(x, y , iseed)

        !     use constants,   only : nzg, zmax, twopi
        !     use photon_vars, only : xp, yp, zp, phi, cosp, sinp, cost, sint, phase

        !     implicit none

        !     real,    intent(IN)    :: x, y
        !     integer, intent(INOUT) :: iseed

        !     real :: ran2

        !     zp = zmax - (1.e-5*(2.*zmax/nzg))

        !     xp = x
        !     yp = y

        !     cost = 2. * ran2(iseed) - 1.
        !     sint = sqrt(1. - cost * cost) 

        !     phi = twopi * ran2(iseed)
        !     cosp = cos(phi)
        !     sinp = sin(phi)

        !     phase = 0.

        ! end subroutine point


        ! subroutine uniform(iseed)

        !     use photon_vars, only : phi, phase, cost, sint, cosp, sinp, xp, yp, zp
        !     use constants,   only : xmax, ymax, zmax, nzg

        !     implicit none

        !     integer, intent(INOUT) :: iseed
        !     real :: ran2

        !     zp = zmax - (1.e-5*(2.*zmax/nzg))
        !     xp = 2. * xmax * ran2(iseed) - xmax
        !     yp = 2. * ymax * ran2(iseed) - xmax

        !     cost = -1.
        !     sint = 1. - cost**2.

        !     phi = 0.
        !     cosp = cos(phi)
        !     sinp = sin(phi)

        !     phase = 0.
        ! end subroutine uniform


        subroutine rang(x, y, avg, sigma, iseed)

            implicit none

            real,    intent(IN) :: avg, sigma
            integer, intent(INOUT) :: iseed
            real, intent(OUT):: x,y
            
            real :: s, tmp

            s = 1.

            do while(s >= 1.)
                x = ranu(-1., 1., iseed)
                y = ranu(-1., 1., iseed)
                s = y**2 + x**2
            end do

            tmp = x*sqrt(-2.*log(s)/s)
            x = avg + sigma*tmp

           tmp = y*sqrt(-2.*log(s)/s)
            y = avg + sigma*tmp

        end subroutine rang


        real function ranu(a, b, iseed)

            implicit none


            real, intent(IN)       :: a, b
            integer, intent(INOUT) :: iseed

            real :: ran2

            ranu = a + ran2(iseed) * (b - a)

        end function ranu

end MODULE sourceph_mod
