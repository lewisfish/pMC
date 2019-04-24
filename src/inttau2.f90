module inttau2

   implicit none
   
   private
   public :: tauint1, find, reflect_refract

contains

    subroutine tauint1(xcell,ycell,zcell,tflag,iseed,delta)
    !optical depth integration subroutine
    !
    !
        use constants,   only : xmax, ymax, zmax, fact, energy
        use photon_vars, only : xp, yp, zp, phase
        use opt_prop, only : kappa
        use iarray, only : phasor

        implicit none

        real,    intent(IN)    :: delta
        integer, intent(INOUT) :: xcell, ycell, zcell, iseed
        logical, intent(INOUT) :: tflag

        real    :: tau, taurun, taucell, xcur, ycur, zcur, d, dcell, ran2
        integer :: celli, cellj, cellk
        logical :: dir(3)
        complex :: phasec

        xcur = xp + xmax
        ycur = yp + ymax
        zcur = zp + zmax

        celli = xcell
        cellj = ycell
        cellk = zcell

        taurun = 0.
        d = 0.
        dir = (/.FALSE., .FALSE., .FALSE./)

        tau = -log(ran2(iseed))
        do
  
            dir = (/.FALSE., .FALSE., .FALSE./)
            dcell = wall_dist(celli, cellj, cellk, xcur, ycur, zcur, dir)
            taucell = dcell * kappa!rhokap(celli, cellj, cellk)
            if(taurun + taucell < tau)then
                taurun = taurun + taucell
                d = d + dcell

                phase = phase + dcell

                phasec = energy * cmplx(cos(fact*phase), sin(fact*phase))
                phasor(celli,cellj,cellk) = phasor(celli,cellj,cellk) + phasec
                ! jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell

                call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta)
            else

                dcell = (tau - taurun) /kappa!rhokap(celli, cellj, cellk)
                d = d + dcell
                
                phase = phase + dcell
                phasec = energy * cmplx(cos(fact*phase), sin(fact*phase))
                phasor(celli,cellj,cellk) = phasor(celli,cellj,cellk) + phasec
                ! jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell
                call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .FALSE., dir, delta)
                exit
            end if
            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                tflag = .true.
                exit
            end if

        end do
   
        xp = xcur - xmax
        yp = ycur - ymax
        zp = zcur - zmax
        xcell = celli
        ycell = cellj
        zcell = cellk

    end subroutine tauint1
   

    real function wall_dist(celli, cellj, cellk, xcur, ycur, zcur, dir)
    !funtion that returns distant to nearest wall and which wall that is (x,y or z)
    !
    !
        use iarray,      only : xface, yface, zface
        use photon_vars, only : nxp, nyp, nzp

        implicit none

        real,    intent(INOUT) :: xcur, ycur, zcur
        logical, intent(INOUT) :: dir(:)
        integer, intent(INOUT) :: celli, cellj, cellk
        real                   :: dx, dy, dz


        if(nxp > 0.)then
            dx = (xface(celli+1) - xcur)/nxp
        elseif(nxp < 0.)then
            dx = (xface(celli) - xcur)/nxp
        elseif(nxp == 0.)then
            dx = 100000.
        end if

        if(nyp > 0.)then
            dy = (yface(cellj+1) - ycur)/nyp
        elseif(nyp < 0.)then
            dy = (yface(cellj) - ycur)/nyp
        elseif(nyp == 0.)then
            dy = 100000.
        end if

        if(nzp > 0.)then
            dz = (zface(cellk+1) - zcur)/nzp
        elseif(nzp < 0.)then
            dz = (zface(cellk) - zcur)/nzp
        elseif(nzp == 0.)then
            dz = 100000.
        end if

        wall_dist = min(dx, dy, dz)
        if(wall_dist < 0.)then
            print*,'dcell < 0.0 warning! ',wall_dist,dx,dy,dz,nxp,nyp,nzp
            ! error stop 1
        end if
        if(wall_dist == dx)dir=(/.TRUE., .FALSE., .FALSE./)
        if(wall_dist == dy)dir=(/.FALSE., .TRUE., .FALSE./)
        if(wall_dist == dz)dir=(/.FALSE., .FALSE., .TRUE./)
        if(.not.dir(1) .and. .not.dir(2) .and. .not.dir(3))print*,'Error in dir flag'
      
   end function wall_dist


    subroutine update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, wall_flag, dir, delta)
    !routine that upates postions of photon and calls fresnel routines if photon leaves current voxel
    !
    !
        use photon_vars, only : nxp, nyp, nzp
        use iarray,      only : xface, yface, zface
        use utils,       only : str, red, bold, colour

        implicit none
      
      real,    intent(INOUT) :: xcur, ycur, zcur
      real,    intent(IN)    :: dcell, delta
      integer, intent(INOUT) :: celli, cellj, cellk
      logical, intent(IN)    :: wall_flag, dir(:)
      character(len=32)      :: tmp  
      
      if(wall_flag)then
      
         if(dir(1))then
            if(nxp > 0.)then
               xcur = xface(celli+1) + delta
            elseif(nxp < 0.)then
               xcur = xface(celli) - delta
            else
               print*,'Error in x dir in update_pos', dir, nxp, nyp, nzp
            end if
            ycur = ycur + nyp*dcell 
            zcur = zcur + nzp*dcell
         elseif(dir(2))then
            xcur = xcur + nxp*dcell
            if(nyp > 0.)then
               ycur = yface(cellj+1) + delta
            elseif(nyp < 0.)then
               ycur = yface(cellj) - delta
            else
               print*,'Error in y dir in update_pos', dir, nxp, nyp, nzp
            end if
            zcur = zcur + nzp*dcell
         elseif(dir(3))then
            xcur = xcur + nxp*dcell
            ycur = ycur + nyp*dcell 
            if(nzp > 0.)then
               zcur = zface(cellk+1) + delta
            elseif(nzp < 0.)then
               zcur = zface(cellk) - delta
            else
               print*,'Error in z dir in update_pos', dir, nxp, nyp, nzp
            end if
         else
            tmp = colour('Error in update_pos... '//str(dir), red, bold)
            error stop 1
         end if
      else
      
         xcur = xcur + nxp*dcell
         ycur = ycur + nyp*dcell 
         zcur = zcur + nzp*dcell
      
      end if


      if(wall_flag)then
         call update_voxels(xcur, ycur, zcur, celli, cellj, cellk)
      end if
      
    end subroutine update_pos


    subroutine update_voxels(xcur, ycur, zcur, celli, cellj, cellk)
    !updates the current voxel based upon position
    !
    !
        ! use iarray, only : xface, yface, zface
        use constants, only : xmax, ymax, zmax, nxg, nyg, nzg

        implicit none

        real,    intent(IN)    :: xcur, ycur, zcur
        integer, intent(INOUT) :: celli, cellj, cellk

        celli = floor(nxg * (xcur) / (2. * xmax)) + 1!find(xcur, xface) 
        cellj = floor(nyg * (ycur) / (2. * ymax)) + 1!find(ycur, yface)
        cellk = floor(nzg * (zcur) / (2. * zmax)) + 1!find(zcur, zface) 

        if(celli > nxg .or. celli < 1)celli = -1
        if(cellj > nyg .or. cellj < 1)cellj = -1
        if(cellk > nzg .or. cellk < 1)cellk = -1


    end subroutine update_voxels


    integer function find(val, a)
    !searchs for bracketing indicies for a value val in an array a
    !
    !
        implicit none

        real, intent(IN) :: val, a(:)
        integer          :: n, lo, mid, hi

        n = size(a)
        lo = 0
        hi = n + 1

        if (val == a(1)) then
            find = 1
        else if (val == a(n)) then
            find = n-1
        else if((val > a(n)) .or. (val < a(1))) then
            find = -1
        else
            do
                if (hi-lo <= 1) exit
                mid = (hi+lo)/2
                if (val >= a(mid)) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            find = lo
        end if
    end function find


        subroutine reflect_refract(I, N, n1, n2, iseed, rflag)

            use vector_class

            implicit none

            type(vector), intent(INOUT) :: I
            type(vector), intent(INOUT) :: N
            real,         intent(IN)    :: n1, n2
            integer,      intent(INOUT) :: iseed
            logical,      intent(OUT)   :: rflag

            real :: ran2

            rflag = .FALSE.

            if(ran2(iseed) <= fresnel(I, N, n1, n2))then
                ! call reflect(I, N)
                rflag = .true.
            else
                call refract(I, N, n1/n2)
            end if

        end subroutine reflect_refract


        subroutine reflect(I, N)
        !   get vector of reflected photon
        !
        !
            use vector_class

            implicit none

            type(vector), intent(INOUT) :: I
            type(vector), intent(IN)    :: N

            type(vector) :: R

            R = I - 2. * (N .dot. I) * N
            I = R

        end subroutine reflect


        subroutine refract(I, N, eta)
        !   get vector of refracted photon
        !
        !
            use vector_class

            implicit none

            type(vector), intent(INOUT) :: I
            type(vector), intent(IN)    :: N
            real,         intent(IN)    :: eta

            type(vector) :: T, Ntmp

            real :: c1, c2

            Ntmp = N

            c1 = (Ntmp .dot. I)
            if(c1 < 0.)then
                c1 = -c1
            else
                Ntmp = (-1.) * N
            end if

            c2 = sqrt(1.d0 - (eta)**2 * (1.d0 - c1**2))

            T = eta*I + (eta * c1 - c2) * Ntmp 

            I = T

        end subroutine refract


        function fresnel(I, N, n1, n2) result (tir)
        !   calculates the fresnel coefficents
        !
        !
            use vector_class
            use ieee_arithmetic, only : ieee_is_nan

            implicit none

            real, intent(IN)         :: n1, n2
            type(vector), intent(IN) :: I, N

            real             ::  costt, sintt, sint2, cost2, tir, f1, f2

            costt = abs(I .dot. N)

            sintt = sqrt(1. - costt * costt)
            sint2 = n1/n2 * sintt
            if(sint2 > 1.)then
                tir = 1.0
                return
            elseif(costt == 1.)then
                tir = 0.
                return
            else
                sint2 = (n1/n2)*sintt
                cost2 = sqrt(1. - sint2 * sint2)
                f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
                f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

                tir = 0.5 * (f1 + f2)
            if(ieee_is_nan(tir) .or. tir > 1. .or. tir < 0.)print*,'TIR: ', tir, f1, f2, costt,sintt,cost2,sint2
                return
            end if
        end function fresnel
end module inttau2