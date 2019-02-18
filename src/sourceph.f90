MODULE sourceph_mod

    use vector_class, only : vector, magnitude

    implicit none

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

            ! call bessel(raxi, dtoskin, iseed)
            ! do
            ! call gaussian_plano(iseed)
            call gaussian_aspeheric(iseed)
            xcell = int(nxg * (xp + xmax) / (2. * xmax)) + 1
            ycell = int(nyg * (yp + ymax) / (2. * ymax)) + 1
            zcell = int(nzg * (zp + zmax) / (2. * zmax)) + 1
            ! print*,xcell,ycell,zcell,nzg
            ! if(xcell <= nxg .and. ycell <= nyg .and. zcell <= nzg)exit
            ! end do
            ! stop
        end subroutine sourceph


        subroutine bessel(Raxi, d, iseed)

            use constants,   only : nzg, xmax, zmax, ymax,pi, twopi
            use opt_prop,    only : n, wavelength
            use photon_vars, only : xp, yp, zp, sint, cost, sinp, cosp, phi, phase, nxp, nyp, nzp, l

            implicit none

            integer, intent(INOUT) :: iseed
            real, intent(IN) :: raxi, d

            real :: r_pos, tana, x0, y0, z0, dist!,a1,r0,f

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
            z0 = (r_pos* tana) + d + zp!11.1mm

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

            ! A1 = wavelength
            ! R0 = 1d-3
            ! f = 10d-3

            phase = phase + l*phi*wavelength/((n - 1.d0)*twopi)  !higher order bessel shizz => helical axicon
            !airy beam? -(((xp)**2 + (yp)**2)/(2.d0*f))-(10.)*a1*(((xp)**3+(yp)**3)/(r0**3))
            ! phase = phase + l*modulo(phi, twopi/l)*wavelength/((n - 1.d0)*2.*twopi)  !higher order bessel shizz => helical axicon

           

         end subroutine bessel


        subroutine gaussian_plano(iseed)

            use constants,   only : xmax, ymax, zmax, nzg, pi
            use photon_vars, only : xp, yp, zp, phase, nxp, nyp, nzp, sint, cost, sinp, cosp, phi, initp
            use inttau2,     only : reflect_refract
 

            use vector_class

            implicit none

            integer, intent(INOUT) :: iseed

            real :: xi, yi, zi, focalBack, lensThickness, R, t, ZL, n1, n2, distLens, distAir, radlens, thetamin, thetamax, dist,c
            real :: ran2
            type(vector) :: orig, dirI, centre, lensI, N, dirL, lensF, final, dirF, dir
            logical :: flag

            lensThickness = 2.2d-3
            focalBack = 8.5d-3
            R = 4.6d-3
            n1 = 1.d0
            n2 = 1.4585d0
            centre = vector(0.d0, 0.d0, (focalBack + lensThickness) - R)
            radlens = 2.5d-3

            zi = focalBack + lensThickness + 2.*zmax
            ! do
            call rang(xi, yi, 0.d0, 1.d-3/4.d0, iseed)
                ! if(sqrt(xi**2 + yi**2) > radlens)then
                !     cycle
                ! else
                    ! exit
                ! end if
            ! end do
            orig = vector(xi, yi, zi)
            dirI = vector(0.d0, 0.d0, -1.d0)

            if(.not. intersect(orig, dirI, t, centre, R))error stop "No lens intersection!"
            lensI = orig + t*dirI

            N = (lensI - centre)
            N = N%magnitude()
            ! print*,N
            dirL = dirI

            call reflect_refract(dirL, N, n1, n2, iseed, flag)

            zL = focalBack
            distLens = (zL - lensI%z) / (dirL%z)
            ! print*,zl, lensI%z, dirL%z

            if(distLens < 0.d0)then
                print*,zl, lensI%z, dirL%z
                error stop "Negative distance travelled in lens!"
            end if
            lensF = lensI + distLens*dirL

            final = vector(ranu(-xmax, xmax, iseed), ranu(-ymax, ymax, iseed), zmax - (1.e-5*(2.*zmax/nzg)))

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

            initp = (1.+abs(cost))/2.

        end subroutine gaussian_plano


        subroutine gaussian_aspeheric(iseed)

            use constants,   only : xmax, ymax, zmax, nzg
            use photon_vars, only : xp, yp, zp, phase, nxp, nyp, nzp, sint, cost, sinp, cosp, phi
            use inttau2,     only : reflect_refract
 

            use vector_class

            implicit none

            integer, intent(INOUT) :: iseed

    real :: xi, yi, zi, focalBack, lensThickness, R, t, ZL, n1, n2, distLens, distAir, radlens, as(5),k,ti(3)
            real :: zadj, l0, m0, n0, g0,c,s0,z0,z0bar,e,E1,c1,L,m12,ml2,p2,f,fdash,g
            type(vector) :: orig, dirI, centre, lensI, N, dirL, lensF, final, dirF
            logical :: flag
            integer :: i

            ! do i = 1, 1000
            !thor labs lens 354220-C
            as = [4.789735d-4, 4.049692d-6, 3.128181d-8, -6.498699d-10, 0.d0]
            k = -0.925522
            lensThickness = 3.434d-3
            focalBack = 5.9d-3
            R = 4.638124d-3
            n1 = 1.d0
            n2 = 1.586d0

            radlens = 9.936d-3 / 2.d0

            call rang(xi, yi, 0.d0, 1.d-3/4.d0, iseed)
            zi = focalBack + lensThickness + 2.*zmax
            ! do i = 1, 1000
            ! zi = focalBack + lensThickness + .5d-3
            call rang(xi, yi, 0.d0, 1.d-3/4.d0, iseed)
            orig = vector(xi, yi, zi)
            dirI = vector(0.d0, 0.d0, -1.d0)

            centre = vector(0.d0, 0.d0, (focalBack + lensThickness) - R)
            ! print*,orig
            zi = focalBack + lensThickness ! @ plane from where sag is defined

            zadj = aspheric(sqrt(xi**2+yi**2), R, k, as)
            zi = zi - zadj

            lensI = vector(xi, yi, zi) !on front surface of lens
            ! print*,lensI
            !!get normal
            N = grad(lensI%x, lensI%y, R, k, as, .9d-3)!grad_aspher(xi, yi, R, k, as)


            dirL = dirI
            call reflect_refract(dirL, N, n1, n2, iseed, flag)

            zL = focalBack
            distLens = (zL - lensI%z) / dirL%z
 
            if(distLens < 0.d0)then
                error stop "neg distance"
            end if
            lensF = lensI + distLens*dirL !@ lend face back
            ! print*,lensF
                  
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

            phase = distAir + distLens*n2 + zadj !(orig%z - lensI%z)
            ! print*,xp,yp,zp
            ! print*,
            ! print*,
        ! end do
        ! stop
        end subroutine gaussian_aspeheric


        real function aspheric(rad, radcurve, k, as)

            implicit none

            real, intent(IN) :: as(:), radcurve, k, rad

            integer :: i, n

            n = size(as)

            !$\frac{Y^2}{R\left(1+\sqrt{1-(1+k)\frac{Y^2}{R^2}}\right)}$
            aspheric = (rad**2) / (radcurve * (1.d0 + sqrt(1.d0 - (1.d0 + k)) * (rad / radcurve)**2))

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

        type(vector) function grad_aspher(x, y, R, k, as)

            implicit none

            real, intent(IN) :: x, y, as(:), R, k
            real :: x1, x2, x3, x4, y1, y2, y3, y4
            real :: x3top, x3bot1, x3bot2, x4top, x4bot1, x4bot2
            real :: y3top, y3bot1, y3bot2, y4top, y4bot1, y4bot2

            grad_aspher = vector(0.d0, 0.d0, 0.d0)

            x1 = 4.d0*as(1)*x*(x**2+y**2)

            x2 = 6.d0*as(2)*x*(x**2+y**2)**2 + 8.d0*as(3)*x*(x**2+y**2)**3 + 10.d0*as(4)*x*(x**2+y**2)**4

            x3top = (1.d0+k)*x*sqrt(x**2+y**2)
            x3bot1 = R**3*sqrt(1.d0 - (((1.d0+k)*(x**2+y**2))/R**2))
            x3bot2 = (1.d0+sqrt(1.d0 - (((1.d0+k)*(x**2+y**2))/R**2)))**2

            x3 = x3top / (x3bot1 * x3bot2)

            x4top = x
            x4bot1 = R*sqrt(x**2+y**2) 
            x4bot2 = (1.d0+sqrt(1.d0 - (((1.d0+k)*(x**2+y**2))/R**2)))

            x4 = x4top / (x4bot1 * x4bot2)

            grad_aspher%x = (x1 + x2 - x3 - x4)


            y1 = 4.d0*as(1)*y*(x**2+y**2)
            y2 = 6.d0*as(2)*y*(x**2+y**2)**2 + 8.d0*as(3)*y*(x**2+y**2)**3 + 10.d0*as(4)*y*(x**2+y**2)**4

            y3top = (1.d0+k)*y*sqrt(x**2+y**2)
            y3bot1 = R**3*sqrt(1.d0 - (((1.d0+k)*(x**2+y**2))/R**2))
            y3bot2 = (1.d0+sqrt(1.d0 - (((1.d0+k)*(x**2+y**2))/R**2)))**2

            y3 = y3top / (y3bot1 * y3bot2)

            y4top = y
            y4bot1 = R*sqrt(x**2+y**2) 
            y4bot2 = (1.d0+sqrt(1.d0 - (((1.d0+k)*(x**2+y**2))/R**2)))

            y4 = y4top / (y4bot1 * y4bot2)

            grad_aspher%y = (y + y2 - y3 - y4)

            grad_aspher%z = -1.d0

            grad_aspher = grad_aspher%magnitude()

        end function grad_aspher

        logical function solveQuadratic(a, b, c, x0, x1)
        ! solves quadratic equation given coeffs a, b, and c
        ! returns true if real soln
        ! returns x0 and x1
        ! adapted from scratchapixel

            implicit none

            real, intent(IN)  :: a, b, c
            real, intent(OUT) :: x0, x1

            real :: discrim, q

            solveQuadratic = .false.

            discrim = b**2 - 4.d0 * a * c
            if(discrim < 0.d0)then
                return
            elseif(discrim == 0.d0)then
                x0 = -0.5*b/a
                x1 = x0
            else
                if(b > 0.d0)then
                    q = -0.5d0 * (b + sqrt(discrim))
                else
                    q = -0.5d0 * (b - sqrt(discrim))
                end if
                x0 = q / a
                x1 = c / q
            end if
            solveQuadratic = .true.
            return

        end function solveQuadratic


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

end MODULE sourceph_mod

            ! real :: xi, yi, zi, R, t, n1, n2, lensThickness, zl, xl, yl, distLens, distAir, backfocal
            ! type(vector) :: posi, diri, centre, lensS, N, dirL
            ! logical :: flag

            ! R = 8.7d-3!4.6d-3
            ! n1 = 1.d0
            ! n2 = 1.4338d0!1.458d0
            ! lensThickness = 4.3d-3!2.2d-3
            ! backfocal = 17.0d-3!8.5d-3

            ! !get x, y in gaussian dist
            ! call rang(xi, yi, 0.d0, 1.d-3/4.d0, iseed)
            ! zi = -10.d-3

            ! !init pos/direction of photon
            ! posi = vector(xi, yi, zi)
            ! diri = vector(0.d0, 0.d0, 1.d0)

            ! !centre of lens defining sphere
            ! centre = vector(0.d0, 0.d0, R)

            ! !get intersection point as func of t
            ! if(.not.intersect(posi, diri, t, centre, R))error stop

            ! !get real po on lens surface
            ! lensS = posi + t * diri
            
            ! !get normal on lens surface
            ! N = (lensS - centre)
            ! N = N%magnitude()

            ! !do fresnel refraction
            ! flag = .false.
            ! call reflect_refract(diri, N, n1, n2, iseed, flag)
            ! dirL = diri

            ! !move photon to far side of lens
            ! zL = lensThickness
            ! distLens = (zL - lensS%z) / dirL%z
            ! if(distLens < 0.d0)error stop
            ! xL = xi + distLens * dirL%x
            ! yL = yi + distLens * dirL%y

            ! !sample on medium surface
            ! xp = ranu(-xmax, xmax, iseed)
            ! yp = ranu(-ymax, ymax, iseed)
            ! zp = (backfocal + lensThickness) - zmax

            ! distAir = sqrt((xL - xp)**2 + (yL - yp)**2 + (zL - zp)**2)

            ! nxp = (xp - xL) / distAir
            ! nyp = (yp - yL) / distAir
            ! ! force -ive as thats how MCRT code is setup
            ! nzp = -abs((zp - zL) / distAir)

            ! phi = atan2(nyp, nxp)
            ! cosp = cos(phi)
            ! sinp = sin(phi)

            ! cost = nzp
            ! sint = sqrt(1.d0 - cost**2)

            ! phase = distAir + distLens*n2 + lensS%z
            ! zp = zmax - (1.e-5*(2.*zmax/nzg))