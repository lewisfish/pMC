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

            ! call bessel(raxi, dtoskin, iseed)
            call gaussian(iseed)

            xcell = int(nxg * (xp + xmax) / (2. * xmax)) + 1
            ycell = int(nyg * (yp + ymax) / (2. * ymax)) + 1
            zcell = int(nzg * (zp + zmax) / (2. * zmax)) + 1
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


        subroutine gaussian(iseed)

            use constants,   only : xmax, ymax, zmax, nzg
            use photon_vars, only : xp, yp, zp, phase, nxp, nyp, nzp, sint, cost, sinp, cosp, phi
            use inttau2,     only : reflect_refract

            use vector_class

            implicit none

            integer, intent(INOUT) :: iseed

            real         :: radius, n2, n1, LensMaxt, t, dist, x0, y0, z0, tmp, focald, Pplane
            logical      :: rflag
            type(vector) :: I, N, dir, orig, centre, phit, nhit, lensF

            ! thor labs conve planar lens:LA4249
            ! n from @500nm https://www.filmetrics.com/refractive-index-database/SiO2/Fused-Silica-Silicon-Dioxide-Thermal-Oxide-ThermalOxide

            radius   = 4.6d-3   ! radius of sphere that defines the convex part of lens
            n2       = 1.4585d0  ! refractive index of lens @ 587.6nm
            n1       = 1.d0     ! refractive index of ambient medium
            LensMaxt = 2.2d-3   ! lens thickness at thickest part
            focald = 1.d0/((n2-1.d0)*(1.d0/radius))
            Pplane = LensMaxt / n2

!                                      
!                    --|                    
!                   /  |                    
!                  /   |                    
!------------------    |                    
!                / \   |                    
!               /   \  |
!              |     \ |                 
!              |      \|--------------------
!               \      |                    
!                \     |                    
!                 \    |                    
!                  \   |
!                   \  |
!                    --|
!              0       LensMaxt   


            ! pull photon portion from gaussian dist
            call rang(xp, yp, 0.d0, 1.d-3/4.d0, iseed)
            ! init pos of photon
            orig = vector(xp, yp, -1.5d-3)
            ! init dir of photon. Towards lens
            dir = vector(0.d0, 0.d0, 1.d0)
            ! pos of len, based upon where centre of lens defining sphere is
            centre = vector(0.d0, 0.d0, radius)

            ! get intersection pos as function of t
            rflag = intersect(orig, dir, t, centre, radius)
            ! get real pos of intersection
            phit = orig + t*dir
            ! get normal at intersection for reflection/refraction
            nhit = (phit - centre)
            nhit = nhit%magnitude()

            I = dir
            I = I%magnitude()
            N = nhit
            rflag = .false.

            ! do fresnel calculation
            call reflect_refract(I, N, n1, n2, iseed, rflag)

            dir = I
            dist = (LensMaxt - phit%z) / I%z
            !pos on lens, planar side
            lensF = phit + dist * dir

            ! draw x, y randomly on surface of medium
            x0 = ranu(-xmax, xmax, iseed)
            y0 = ranu(-ymax, ymax, iseed)
            ! set z on surface of medium
            z0 = LensMaxt + (focald - Pplane) - zmax   ! lens thickness + mechanical focal length - half medium size, so that focal point falls in middle of medium

            !get distance from lens surface to surface of medium
            tmp = sqrt((x0 - lensF%x)**2 + (y0 - lensF%y)**2 + (z0 - lensF%z)**2)

            ! dist travelled in air before lens + distance travelled in lens + distance to surface of medium
            phase = phit%z + dist*n2 + tmp
            

            nxp = (x0 - lensF%x) / tmp
            nyp = (y0 - lensF%y) / tmp
            nzp = -abs((z0 - lensF%z) / tmp) !-ive due to way z pos is defined in MCRT code

            cost = nzp
            sint = sqrt(1.d0 - cost**2)

            phi = atan2(nyp, nxp)

            cosp = cos(phi)
            sinp = sin(phi)

            xp = x0
            yp = y0
            zp = zmax - (1.e-5*(2.*zmax/nzg))
        end subroutine gaussian


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
