MODULE sourceph_mod

    implicit none
    save

    private
    public :: sourceph

    CONTAINS
        subroutine sourceph(xcell, ycell, zcell, iseed)

            use constants,   only : nxg, nyg, nzg, beam, xmax, ymax, zmax
            use photon_vars, only : xp, yp, zp, nxp, nyp, nzp, sint, cost, sinp, cosp, phi


            implicit none

            integer,          intent(OUT)   :: xcell, ycell, zcell
            integer,          intent(INOUT) :: iseed

            real    :: w, x, y, ran2
            integer :: error

            w = 0.0000075

            select case (beam)
            case('bessel')
                call bessel(w, iseed)
                nxp = sint * cosp  
                nyp = sint * sinp
                nzp = cost
            case('gaussian')
                call gaussian_phase(iseed)
                phi = atan2(nyp,nxp)
                cosp = cos(phi)
                sinp = sin(phi)

                cost = nzp
                sint = sqrt(1. - cost**2.)
            case('point')
                x = 0.
                y = 0.
                if(ran2(iseed) < .5)x = x - xmax/10.
                call point(x, y, iseed)
                nxp = sint * cosp  
                nyp = sint * sinp
                nzp = cost
            case('circular')
                call circular(w, iseed)
                nxp = sint * cosp  
                nyp = sint * sinp
                nzp = cost
            case('uniform')
                call uniform(iseed)
                nxp = sint * cosp  
                nyp = sint * sinp
                nzp = cost
            case default
                call mpi_finalize(error)
                error stop 'Error, unrecognised beam type!'
            end select




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
            angle = 190.

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


        subroutine gaussian_phase(iseed)
        ! f => focal length of lens
        ! zr => rayleigh length
        ! zf => z location of of the focus


            use constants,   only : pi, nzg, zmax, twopi, xmax
            use photon_vars, only : nxp, nyp, nzp, xp,yp,zp, phase,zr
            use opt_prop,    only : wavelength

            implicit none

            integer, intent(INOUT) :: iseed
            real :: ran2
            real :: xo, yo, zo, zf, f, fact, rz,  D, wo, r1, phigauss, r_pos

            zo = zmax - (1.e-5*(2.*zmax/nzg))
            xo = 0.

            f = zmax/8.
            D = 1.e-6
            wo = (2./pi) * (wavelength*f)/D
            zr = (pi * wo**2)/wavelength

            !gaussian beam via box-muller method
            do
                r1 = d*sqrt(-2.*log(ran2(iseed)))
                phigauss = twopi * ran2(iseed)
                xo = r1 * cos(phigauss)
                yo = r1 * sin(phigauss)
                if(xo**2 + yo**2 < xmax**2.)exit
            end do

            xp = -xo*((zo - zf)/sqrt(xo**2 + yo**2 + f**2)) 
            yp = -yo*((zo - zf)/sqrt(xo**2 + yo**2 + f**2)) 
            zp = zf + f*((zo - zf)/sqrt(xo**2 + yo**2 + f**2)) 

            rz = -(zp - zf) * (1 + (zr / (zp - zf))**2)

            fact = 1./sqrt(1 + ((xp**2 + yp**2) / rz**2))

            nxp = fact * (xp/rz)
            nyp = fact * (yp/rz)
            nzp = -1. * fact

            r_pos = sqrt(xp**2 + yp**2)
            phase = cos(twopi*(abs(zp-zo)+r_pos)/wavelength)
            zp = zo

        end subroutine gaussian_phase


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
