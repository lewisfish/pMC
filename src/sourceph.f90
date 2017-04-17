MODULE sourceph_mod

    implicit none
    save

    private
    public :: sourceph

    CONTAINS
        subroutine sourceph(xcell, ycell, zcell, iseed)

            use constants,   only : nxg, nyg, nzg, beam, xmax, ymax, zmax
            use photon_vars, only : xp, yp, zp, nxp, nyp, nzp, sint, cost, sinp, cosp


            implicit none

            integer,          intent(OUT)   :: xcell, ycell, zcell
            integer,          intent(INOUT) :: iseed

            real    :: w, x, y
            integer :: error

            w = 0.00005
            x = 0.
            y = 0.

            select case (beam)
            case('bessel')
                call bessel(w, iseed)
            case('gaussian')
                call gaussian(w, iseed)
            case('point')
                call point(x, y, iseed)
            case('circular')
                call circular(w, iseed)
            case('uniform')
                call uniform(iseed)
            case default
                print*,'Error, unrecognised beam type!'
                call exit(0)
                call mpi_finalize(error)
            end select

            nxp = sint * cosp  
            nyp = sint * sinp
            nzp = cost

            xcell = int(nxg * (xp + xmax) / (2. * xmax)) + 1
            ycell = int(nyg * (yp + ymax) / (2. * ymax)) + 1
            zcell = int(nzg * (zp + zmax) / (2. * zmax)) + 1
        end subroutine sourceph


        subroutine bessel(w, iseed)

            use constants,   only : twopi, pi, nzg, xmax, zmax
            use opt_prop,    only : wavelength
            use photon_vars, only : xp, yp, zp, sint, cost, sinp, cosp, phi, phase, angle

            implicit none

            real, intent(IN)       :: w
            integer, intent(INOUT) :: iseed

            real :: r1, phigauss, ran2, theta, r_pos

            zp = zmax - (1.e-5*(2.*zmax/nzg))
            xp = 0.

            !gaussian beam via box-muller method
            do
                r1 = w*sqrt(-2.*log(ran2(iseed)))
                phigauss = twopi * ran2(iseed)
                xp = r1 * cos(phigauss)
                yp = r1 * sin(phigauss)
                if(xp**2 + yp**2 < xmax**2.)exit
            end do

            !axicon lens angle
            angle = 182.

            theta = angle*pi/180.
            cost = cos(theta)
            sint = sin(theta)

            phi = atan2(yp,xp)
            cosp = cos(phi)
            sinp = sin(phi)

            !initial phase
            r_pos = sqrt(xp**2 + yp**2)
            phase = twopi*r_pos/wavelength * sin((angle-180)*pi/180.)
        end subroutine bessel


        subroutine gaussian(w, iseed)

            use constants,   only : twopi, nzg, xmax, zmax
            use photon_vars, only : xp, yp, zp, sint, cost, sinp, cosp, phi, phase

            implicit none

            real, intent(IN)       :: w
            integer, intent(INOUT) :: iseed

            real :: phigauss, r1, ran2

            zp = zmax - (1.e-5*(2.*zmax/nzg))
            xp = 0.

            !gaussian beam via box-muller method
            do
                r1 = w*sqrt(-2.*log(ran2(iseed)))
                phigauss = twopi * ran2(iseed)
                xp = r1 * cos(phigauss)
                yp = r1 * sin(phigauss)
                if(xp**2 + yp**2 < xmax**2.)exit
            end do

            cost = -1.
            sint = 1. - cost**2.

            phi = 0.
            cosp = cos(phi)
            sinp = sin(phi)

            !initial phase
            phase = 0.
        end subroutine gaussian


        subroutine circular(radius, iseed)

            use constants, only : nzg, zmax
            use photon_vars, only : xp, yp, zp, sint, cost, sinp, cosp, phi, phase

            implicit none

            real,    intent(IN)    :: radius
            integer, intent(INOUT) :: iseed

            real :: ran2

            zp = zmax - (1.e-5*(2.*zmax/nzg))

            do
               xp = 2. * radius * ran2(iseed) - radius
               yp = 2. * radius * ran2(iseed) - radius
               if(xp**2 + yp **2 <= radius**2)exit
            end do

            cost = -1.
            sint = 1. - cost**2.

            phi = 0.
            cosp = cos(phi)
            sinp = sin(phi)

            phase = 0.

        end subroutine circular


        subroutine point(x, y , iseed)

            use constants,   only : nzg, zmax, twopi
            use photon_vars, only : xp, yp, zp, phi, cosp, sinp, cost, sint, phase

            implicit none

            real,    intent(IN)    :: x, y
            integer, intent(INOUT) :: iseed

            real :: ran2

            zp = zmax - (1.e-5*(2.*zmax/nzg))

            xp = x
            yp = y

            cost = 2. * ran2(iseed) - 1.
            sint = 1. - cost * cost 

            phi = twopi * ran2(iseed)
            cosp = cos(phi)
            sinp = sin(phi)

            phase = 0.

        end subroutine point

        subroutine uniform(iseed)

            use photon_vars, only : phi, phase, cost, sint, cosp, sinp, xp, yp, zp
            use constants,   only : xmax, ymax, zmax, nzg

            implicit none

            integer, intent(INOUT) :: iseed
            real :: ran2

            zp = zmax - (1.e-5*(2.*zmax/nzg))
            xp = 2. * xmax * ran2(iseed) - xmax
            yp = 2. * ymax * ran2(iseed) - xmax

            cost = -1.
            sint = 1. - cost**2.

            phi = 0.
            cosp = cos(phi)
            sinp = sin(phi)

            phase = 0.
        end subroutine uniform

end MODULE sourceph_mod
