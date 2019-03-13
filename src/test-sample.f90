module m

    use vector_class

    implicit none

    real, parameter :: PI=4.d0*atan(1.d0)

    type :: SphQuad
        type(vector) :: o, x, y, z !local ref system "R"
        real :: z0, z0sq                    !
        real :: x0, y0, y0sq                ! rectangle cords in "R"
        real :: x1, y1, y1sq                !
        real :: b0, b1, b0sq, k          ! misc precomputer consts
        real :: s                        ! solid angle of "Q"
    end type SphQuad
    
    contains
    
        real function clamp(a, x, y)

            implicit none

            real, intent(IN) :: x, y, a

            if(a < x)then
                clamp = x
                return
            elseif(a > y)then
                clamp = y
                return
            else
                clamp = a
                return
            end if

        end function clamp

        function SphQuadSample(squad, u, v)

            implicit none

            type(vector) :: SphQuadSample
            type(SphQuad), intent(IN) :: squad
            real,          intent(IN) :: u, v

            real :: au, fu, cu, xu, d, h0, h1, hv, hv2, yv, eps = 1.d-8, futmp

            au = u * squad%s + squad%k
            fu = (cos(au) * squad%b0 - squad%b1) / sin(au)
            if(fu > 0)then
                futmp = 1.d0
            else
                futmp = -1.d0
            end if
            cu = 1.d0 / sqrt(fu**2 + squad%b0sq) * (futmp)
            ! cu = clamp(cu, -1., 1.)

            xu = -(cu * squad%z0) / sqrt(1.d0 - cu**2)
            ! xu = clamp( xu, squad%x0, squad%x1)
            d = sqrt(xu**2 + squad%z0sq)
            h0 = squad%y0 / sqrt(d**2 + squad%y0sq)
            h1 = squad%y1 / sqrt(d**2 + squad%y1sq)
            hv = h0 + v * (h1 - h0)
            hv2 = hv**2
            if(hv2 < 1. - eps)then
                yv = (hv*d) / sqrt(1. -hv2)
            else
                yv = squad%y1
            end if 

            SphQuadSample = squad%o + xu*squad%x + yv*squad%y + squad%z0*squad%z

        end function SphQuadSample


        subroutine SphQuadInit(squad, s, ex, ey, o)

            implicit none

            type(SphQuad), intent(INOUT)   :: squad
            type(vector), intent(IN) :: s, ex, ey, o

            type(vector) :: d, v00,v10,v01,v11,n1,n0,n2,n3
            real :: ex1, ey1, g0,g1,g2,g3

            squad%o = o

            ex1 = ex%length()
            ey1 = ey%length()

            squad%x = ex / ex1
            squad%y = ey / ey1

            squad%z = squad%x .cross. squad%y

            d = s - o
            squad%z0 = d .dot. squad%z

            if(squad%z0 > 0.)then
                squad%z = squad%z * (-1.)
                squad%z0 = squad%z0 * (-1.)
            end if

            squad%z0sq = squad%z0**2
            squad%x0 = d .dot. squad%x
            squad%y0 = d .dot. squad%y
            squad%x1 = squad%x0 + ex1
            squad%y1 = squad%y0 + ey1
            squad%y0sq = squad%y0**2
            squad%y1sq = squad%y1**2

            v00 = vector(squad%x0, squad%y0, squad%z0)
            v01 = vector(squad%x0, squad%y1, squad%z0)
            v10 = vector(squad%x1, squad%y0, squad%z0)
            v11 = vector(squad%x1, squad%y1, squad%z0)

            n0 = (v00 .cross. v10)
            n0 = n0%magnitude()

            n1 = (v10 .cross. v11)
            n1 = n1%magnitude()

            n2 = (v11 .cross. v01)
            n2 = n2%magnitude()

            n3 = (v01 .cross. v00)
            n3 = n3%magnitude()


            g0 = acos(-(n0 .dot. n1))
            g1 = acos(-(n1 .dot. n2))
            g2 = acos(-(n2 .dot. n3))
            g3 = acos(-(n3 .dot. n0))

            squad%b0 = n0%z
            squad%b1 = n2%z
            squad%b0sq = squad%b0**2
            squad%k = 2.d0 * PI - g2 - g3
            squad%s = g0 + g1 - squad%k

        end subroutine SphQuadInit
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
end module m

program p

    use m

    implicit none
    
    type(SphQuad) :: squad
    type(vector)  :: ex, ey, s, o, ans
    real :: u, v, ran2,a,b,avg,sigma
    integer :: iseed, i

    iseed = -83432432

    ! location of screen vertex
    s = vector(0., 0., 10d-2)
    ! screen vectors
    ex = vector(.5d-3, 0., 0.)
    ey = vector(0., .5d-3, 0.)
    avg = 0.
    sigma = 1.d-3 / 4.
    !init squad struct
    do i = 1, 10000
    !origin of point
    call rang(a, b, avg, sigma, iseed)
    ! a = 0.d0
    ! b = a
    o = vector(a, b, 0.)
    call SphQuadInit(squad, s, ex, ey, o)

    !sample u, v
    u = ran2(iseed)
    v = ran2(iseed)
    ! get pos on screen
    ans = SphQuadSample(squad, u, v)
    print*,ans%x, ans%y
end do
end program p