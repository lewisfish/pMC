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

            w = 0.00005

            select case (beam)
            case('bessel')
                call bessel(iseed)
            case('gaussian')
                call gaussian_phase(iseed)
                phi = atan2(nyp,nxp)
                cosp = cos(phi)
                sinp = sin(phi)

                cost = nzp
                sint = sqrt(1. - cost**2.)
            case('point')
                x = xmax/4.!need large $\lambda$ small seperation for a few maxima
                y = 0.
                if(ran2(iseed) < .5)x = -xmax/4.
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


        subroutine bessel(iseed)

            use constants,   only : twopi, pi, nzg, xmax, zmax, nxg, nyg, ymax
            use opt_prop,    only : wavelength
            use photon_vars, only : xp, yp, zp, sint, cost, sinp, cosp, phi, phase, angle, nxp, nyp, nzp

            implicit none

            integer, intent(INOUT) :: iseed

            real :: r1, phigauss, ran2, theta, r_pos, alpha, beta, rpos_skin
            real :: x0, y0, z0, t, t1
            integer,save :: count


            zp = zmax - (1.e-5*(2.*zmax/nzg))
            xp = 0.
            !gaussian beam via margasalia method
            do
                zp = zmax - (1.e-5*(2.*zmax/nzg))
                !sigma = 1/e^2 diameter of beam / 4
                xp = rang(0., .1d0/4.d0, iseed)
                yp = rang(0., .1d0/4.d0, iseed)


                !axicon deflection angle, beta
                angle = 180. + 2.314!18.!2.314could be 5...

                !adjust photon from on axicon to just in skin, using trig see bessel_geo_skin.svg
                alpha = 5.*pi/180.
                beta = angle*pi/180.

                r_pos = sqrt(xp**2 + yp**2)
                rpos_skin = (1. - tan(alpha)*tan(beta))*r_pos
                phi = atan2(yp, xp)

                xp = rpos_skin * cos(phi)
                yp = rpos_skin * sin(phi)

                cost = cos(beta)
                sint = sin(beta)

                cosp = cos(phi)
                sinp = sin(phi)

                nxp = sint * cosp  
                nyp = sint * sinp
                nzp = cost

                x0 = xp
                y0 = yp
                z0 =  zmax + 1.

                t = (zp - z0)/nzp
                yp = nyp * t + y0
                xp = nxp * t + x0


                if(abs(xp) > xmax .or. abs(yp) > ymax)then
                    continue
                else
                    exit
                end if

            end do
            !initial phase
            r_pos = sqrt(xp**2 + yp**2)
            phase = cos(twopi*r_pos/wavelength * sin((angle-180)*pi/180.) + twopi*1.0/wavelength)

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
            real :: ran2,w
            real :: xo, yo, zo, zf, f, fact, rz,  D, wo, r1, phigauss, r_pos


            w = 0.0000000015
            zf = zmax-.01
            zo = zmax - (1.e-5*(2.*zmax/nzg))
            xo = 0.

            f = zmax/2000.
            D = 6.e-8
            wo = (2./pi) * (wavelength*f)/D
            zr = (pi * wo**2)/wavelength
            ! print*,wo,zr
            !gaussian beam via box-muller method
            do
                r1 = w*sqrt(-2.*log(ran2(iseed)))
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
            phase = (twopi*(abs(zp-zo)+r_pos)/wavelength)
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
            sint = sqrt(1. - cost * cost) 

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


        real function rang(avg, sigma, iseed)

            implicit none

            real,    intent(IN) :: avg, sigma
            integer, intent(INOUT) :: iseed
            
            real :: s, u, tmp

            s = 1.

            do while(s >= 1.)
                u = ranu(-1., 1., iseed)
                s = ranu(-1., 1., iseed)
                s = s**2 + u**2
            end do

            tmp = u*sqrt(-2.*log(s)/s)
            rang = avg + sigma*tmp

        end function rang


        real function ranu(a, b, iseed)

            implicit none


            real, intent(IN)       :: a, b
            integer, intent(INOUT) :: iseed

            real :: ran2

            ranu = a + ran2(iseed) * (b - a)

        end function ranu

end MODULE sourceph_mod
