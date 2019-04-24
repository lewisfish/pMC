module sourceph_mod

    use vector_class, only : vector, magnitude

    implicit none

    private
    public :: sourceph, ranu

    contains
        subroutine sourceph(xcell, ycell, zcell,raxi, dtoskin, iseed)

            use constants,   only : nxg, nyg, nzg, xmax, ymax, zmax, beam
            use photon_vars, only : xp, yp, zp

            implicit none

            integer,          intent(OUT)   :: xcell, ycell, zcell
            integer,          intent(INOUT) :: iseed
            real, intent(IN) :: raxi, dtoskin

            if(beam == "bessel")then
                call bessel(raxi, dtoskin, iseed)
            elseif(beam == "gaussian")then
                call gaussian_plano(iseed)
            end if

            xcell = int(nxg * (xp + xmax) / (2. * xmax)) + 1
            ycell = int(nyg * (yp + ymax) / (2. * ymax)) + 1
            zcell = int(nzg * (zp + zmax) / (2. * zmax)) + 1

        end subroutine sourceph


        subroutine bessel(Raxi, d, iseed)

            use constants,   only : nzg, xmax, zmax, ymax,pi, twopi, waist
            use opt_prop,    only : n, wavelength
            use photon_vars, only : xp, yp, zp, sint, cost, sinp, cosp, phi, phase, nxp, nyp, nzp, l

            implicit none

            integer, intent(INOUT) :: iseed
            real,    intent(IN)    :: raxi, d

            real :: r_pos, tana, x0, y0, z0, dist


            n = Sellmeier(wavelength * 1d6)

            phase = 0.d0
            tana = tan(getAlpha(n, waist, getfocal(4.6d-3, n)))!tan(5.d0 * pi / 180.d0)!

            !generate gaussian profile at top of axicon
            call rang(xp, yp, 0.d0, sqrt(2.d0)*waist/4.d0, iseed)
            r_pos = sqrt(xp**2 + yp**2)
            zp = (raxi - r_pos) * tana

            !distance through axicon
            phase = zp * n

            !sample mediums surface
            x0 = ranu(-xmax , xmax , iseed)
            y0 = ranu(-ymax , ymax , iseed)
            z0 = (r_pos* tana) + d + zp !11.1mm   => r_pos*tana is distane from lens surface to plane at tip of axicon

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

            if(l > 0)then
                phase = phase + (l*wavelength * modulo(phi, twopi / l)) / (twopi)
            end if

            !airy beam? -(((xp)**2 + (yp)**2)/(2.d0*f))-(10.)*a1*(((xp)**3+(yp)**3)/(r0**3))
         end subroutine bessel


        subroutine gaussian_plano(iseed)

            use constants,   only : xmax, ymax, zmax, nzg, waist
            use photon_vars, only : xp, yp, zp, phase, nxp, nyp, nzp, sint, cost, sinp, cosp, phi, l
            use opt_prop,    only : wavelength
            use inttau2,     only : reflect_refract
 
            use vector_class

            implicit none

            integer, intent(INOUT) :: iseed

            real         :: R, t, n1, n2, lensThickness, distLens, distAir, backfocal, zL, xi,yi,zi
            type(vector) :: orig, diri, centre, LensI, N, dirL, lensF
            logical      :: flag

            R = 4.6d-3
            n1 = 1.d0
            n2 = Sellmeier(wavelength * 1d6)
            phase = 0.d0
            lensThickness = 2.2d-3
            backfocal = getfocal(r, n2) - (lensThickness / n2) !8.5d-3

            !get x, y in gaussian dist
            call rang(xi, yi, 0.d0, sqrt(2.)*waist/4.d0, iseed)
            do

            zi = -10.d-3

            !init pos/direction of photon
            orig = vector(xi, yi, zi)
            dirI = vector(0.d0, 0.d0, 1.d0)

            !centre of lens defining sphere
            centre = vector(0.d0, 0.d0, R)

            !get intersection point as func of t
            if(intersect(orig, dirI, t, centre, R))exit
            end do

            !get real po on lens surface
            LensI = orig + t * dirI
            
            !get normal on lens surface
            N = (LensI - centre)
            N = N%magnitude()

            !do fresnel refraction
            flag = .false.
            call reflect_refract(dirI, N, n1, n2, iseed, flag)
            dirL = dirI

            !move photon to far side of lens
            zL = lensThickness
            distLens = (zL - LensI%z) / dirL%z
            if(distLens < 0.d0)then
                error stop
            end if
            LensF = LensI + distLens * dirL

            !sample on medium surface
            xp = ranu(-xmax, xmax, iseed)
            yp = ranu(-ymax, ymax, iseed)
            zp = (backfocal + lensThickness) - 1.5*zmax
            zi = zp

            distAir = sqrt((LensF%x - xp)**2 + (LensF%y - yp)**2 + (LensF%z - zp)**2)

            nxp = (xp - LensF%x) / distAir
            nyp = (yp - LensF%y) / distAir
            ! force -ive as thats how MCRT code is setup
            nzp = -abs((zp - LensF%z) / distAir)

            phi = atan2(nyp, nxp)
            cosp = cos(phi)
            sinp = sin(phi)

            cost = nzp
            sint = sqrt(1.d0 - cost**2)

            phase = distAir + LensI%z + distLens*n2
            zp = zmax - (1.e-5*(2.*zmax/nzg))


        end subroutine gaussian_plano


        subroutine gaussian_aspeheric(iseed)

            use constants,   only : xmax, ymax, zmax, nzg
            use photon_vars, only : xp, yp, zp, phase, nxp, nyp, nzp, sint, cost, sinp, cosp, phi
            use inttau2,     only : reflect_refract
 

            use vector_class

            implicit none

            integer, intent(INOUT) :: iseed

            real         :: xi, yi, zi, focalBack, lensThickness, R, n1, n2, distLens, distAir, radlens, as(5), k
            real         :: zadj, zL
            type(vector) :: orig, dirI, centre, lensI, N, dirL, lensF, final, dirF
            logical      :: flag

            !thor labs lens 354220-C
            as = [8.924167d-5, 4.38436d-7, 0.d0, 0.d0, 0.d0]
            k = -.73128d0
            lensThickness = 5.0d-3
            focalBack = 6.91d-3 + .9d-3
            R = 6.428132d-3
            n1 = 1.d0
            n2 = 1.584d0
            radlens = 5.5d-3 / 2.d0
            centre = vector(0.d0, 0.d0, (focalBack + lensThickness) - R)

            zi = focalBack + lensThickness + zmax
            do
                call rang(xi, yi, 0.d0, 1.d-3/4.d0, iseed)
                if(sqrt(xi**2 + yi**2) < radlens)exit
            end do
            orig = vector(xi, yi, zi)
            dirI = vector(0.d0, 0.d0, -1.d0)
            zi = focalBack + lensThickness ! @ plane from where sag is defined

            zadj = aspheric(sqrt(xi**2+yi**2), R, k, as)
            zi = zi - zadj

            lensI = vector(xi, yi, zi) !on front surface of lens
            !!get normal
            N = grad(lensI%x, lensI%y, R, k, as, 1d-8)!grad_aspher(xi, yi, R, k, as)

            dirL = dirI
            call reflect_refract(dirL, N, n1, n2, iseed, flag)

            zL = focalBack
            distLens = (zL - lensI%z) / dirL%z
 
            if(distLens < 0.d0)then
                error stop "neg distance"
            end if
            lensF = lensI + distLens*dirL !@ lend face back

            final = vector(ranu(-xmax, xmax, iseed), ranu(-ymax, ymax, iseed), zmax - (1.e-5*(2.*zmax/nzg))) !@ in medium
            distAir = sqrt((final%x - LensF%x)**2 + (final%y - LensF%y)**2 + (final%z - LensF%z)**2)

            dirF = (final - lensF) / distAir

            nxp = dirF%x
            nyp = dirF%y
            nzp = dirF%z
            if(nzp > 0.d0)error stop

            phi = atan2(nyp, nxp)
            cosp = cos(phi)
            sinp = sin(phi)

            cost = nzp
            sint = sqrt(1.d0 - cost**2)

            xp = final%x
            yp = final%y
            zp = final%z

            phase = distAir + distLens*n2 + (orig%z - lensI%z)

        end subroutine gaussian_aspeheric


        real function aspheric(rad, radcurve, k, as)

            implicit none

            real, intent(IN) :: as(:), radcurve, k, rad

            integer :: i, n

            n = size(as)

            !$\frac{Y^2}{R\left(1+\sqrt{1-(1+k)\frac{Y^2}{R^2}}\right)}$
            aspheric = (rad**2) / (radcurve * (1.d0 + sqrt(1.d0 - (1.d0 + k) * (rad / radcurve)**2)))

            do i = 1, n
                aspheric = aspheric + as(i)*rad**(2*n+2)
            end do

        end function aspheric


        type(vector) function grad(x, y, radcurve, k, as, step)

            implicit none

            real, intent(IN) :: x, y, radcurve, k, as(:), step

            real :: rdiffplus, r, dx, dy, dz, rdiffminus

            r = sqrt(x**2 + y**2)
            rdiffplus = sqrt((x + step)**2 + y**2)
            rdiffminus = sqrt((x - step)**2 + y**2)

           
            dx = (aspheric(rdiffplus, radcurve, k, as) - aspheric(rdiffminus, radcurve, k, as)) / (2.*step)
            r = sqrt(x**2 + y**2)
            rdiffplus = sqrt((y + step)**2 + x**2)
            rdiffminus = sqrt((y - step)**2 + x**2)
            dy = (aspheric(rdiffplus, radcurve, k, as) - aspheric(rdiffminus, radcurve, k, as)) / (2.*step)
            dz = -1.d0

            grad = vector(-dx, -dy, dz)
            grad = grad%magnitude()

        end function grad


        logical function intersect(orig, dir, t, centre, radius)
        ! calculates where a line, with origin:orig and direction:dir hits a sphere, centre:centre and radius:radius
        ! returns true if intersection exists
        ! returns t, the paramertised parameter of the line equation
        ! adapted from scratchapixel
            
            use vector_class, only : vector

            implicit none

            type(vector), intent(IN)  :: dir, orig, centre
            real,         intent(OUT) :: t
            real,         intent(IN)  :: radius

            type(vector) :: L
            real         :: t0, t1, a, b, c, tmp

            intersect = .false.

            L = orig - centre
            a = dir .dot. dir
            b = 2.d0 * (dir .dot. L)
            c = (l .dot. l) - radius**2

            if(.not. solveQuadratic(a, b, c, t0, t1))return
            if(t0 > t1)then
                tmp = t1
                t1 = t0
                t0 = tmp
            end if
            if(t0 < 0.d0)then
                t0 = t1
                if(t0 < 0.)return
            end if

            t = t0
            intersect = .true.
            return

        end function intersect


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


        real function Sellmeier(wave)
        ! Sellmeier equation for fused quatrz
        ! I. H. Malitson. Interspecimen comparison of the refractive index of fused silica
            implicit none

            real, intent(IN) :: wave
            real :: wave2, a, b ,c

            wave2 = wave**2

            a = (0.6961663d0*wave2)/(wave2 - 0.0684043d0**2)
            b = (0.4079426*wave2) / (wave2 - .1162414**2)
            c = (0.8974794*wave2) / (wave2 - 9.896161**2)

            Sellmeier = sqrt(1.d0 + (a + b + c))

        end function Sellmeier


        real function getAlpha(n, d, f)
        ! n is refractive index
        ! d is diameter waist 1/e^2
        ! f is focal length
            implicit none

            real, intent(IN) :: n, d, f

            getAlpha = (1.d0 / (n - 1.d0)) * asin((1.75d0 * d) / (4.d0 * f))

        end function getAlpha


        real function getfocal(r, n)

            implicit none

            real, intent(IN) :: r, n

            getfocal = r / (n - 1.d0)

        end function getfocal

end MODULE sourceph_mod