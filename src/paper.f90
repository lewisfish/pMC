module m

    implicit none
    
    contains
   
        real function ranu(a, b, iseed)

            implicit none


            real, intent(IN)       :: a, b
            integer, intent(INOUT) :: iseed

            real :: ran2

            ranu = a + ran2(iseed) * (b - a)

        end function ranu

end module m

program testPaper

    use m, only: ranu

    implicit none
    
    real    :: F, a, wave, distS, xslit, yslit, binwid, x, y, z, xf, yf, zf, cost, phase, twopi,xmax,dist
    integer :: i, N, idx, idy, iseed, j, u, N2
    integer, parameter :: size=100

    complex :: S(-size:size, -size:size)
    real :: out(-size:size,-size:size)

    twopi = 2.d0*4.d0*atan(1.d0)
    iseed = -489732
    xmax = 3000d-6
    binwid = 2.*xmax / real(size)
    S = cmplx(0.d0,0.d0)
    N = 50000
    n2 = 50000
    a =  200.d-6
    xslit = a / 2.d0
    yslit = 100.d0 * xslit 
    wave = 351.d-9 
    F = .8d0
    distS = (((F / a)**2) / 2.)*wave
    distS = 1.d0 / distS
    F = a*sqrt(2.d0/(wave*distS))
    print*,F,distS

    do i = 1, N
        if(mod(i,1000)==0)print*,(real(i)/real(N))*100
        x = ranu(-xslit, xslit, iseed)
        y = ranu(-xslit, xslit, iseed)
        z = 0.d0

        do j = 1, N2
            xf = ranu(-xmax, xmax, iseed)
            yf = ranu(-xmax, xmax, iseed)
            zf = distS
    
            dist = sqrt((x - xf)**2 + (y - yf)**2 + (z - zf)**2)

            cost = -zf/dist

            idx = nint(xf/binwid)
            idy = nint(yf/binwid)
            phase = twopi * distS/(wave * cost)
            ! print*,idx,idy
            S(idx, idy) = S(idx, idy) + cmplx(cos(phase),sin(phase))
        end do
    end do

    open(newunit=u,file="int1.dat",status="replace")
    do i = -size,size
        write(u,*)(cabs(S(i,j))**2,j=-50,50)
    end do
    close(u)

    open(newunit=u,file="out.dat",access="stream",form="unformatted")
    write(u)cabs(S)**2
    close(u)

    out = cabs(S)**2 / maxval(cabs(S)**2)
    open(newunit=u,file="slice2.dat",status="replace")
    do i = -size,size
        write(u,*)binwid*i,out(i,0)
    end do
    close(u)

end program testPaper