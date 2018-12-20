module inttau2

   implicit none
   
   private
   public :: tauint1, find

contains

    subroutine tauint1(xcell,ycell,zcell,tflag,iseed,delta)
    !optical depth integration subroutine
    !
    !
        use constants,   only : xmax, ymax, zmax, fact
        use photon_vars, only : xp, yp, zp, phase
        use iarray,      only : phasor, rhokap
        ! use opt_prop,    only : kappa
        ! use taufind2

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
            taucell = dcell * rhokap(celli, cellj, cellk)
            if(taurun + taucell < tau)then
                taurun = taurun + taucell
                d = d + dcell

                phase = phase + dcell
                phasec = cmplx(cos(fact*phase), sin(fact*phase))
                phasor(celli,cellj,cellk) = phasor(celli,cellj,cellk) + phasec


                call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta)
            else

                dcell = (tau - taurun) /rhokap(celli, cellj, cellk)
                d = d + dcell
                
                phase = phase + dcell
                phasec = cmplx(cos(fact*phase), sin(fact*phase))
                phasor(celli,cellj,cellk) = phasor(celli,cellj,cellk) + phasec

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
            error stop 1
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


    ! subroutine taufind1(xcell,ycell,zcell,delta,taurun)
    ! !   routine to find tau from current position to edge of grid in a direction (nxp,nyp,nzp)
    ! !
    ! !
    !     use photon_vars, only : xp, yp, zp
    !     ! use iarray,      only : rhokap
    !     use opt_prop,    only : wavelength
    !     use constants,   only : xmax, ymax, zmax 
     
    !     implicit none

    !     real,    intent(IN)    :: delta
    !     integer, intent(INOUT) :: xcell, ycell, zcell

    !     real                   :: taurun, taucell, xcur, ycur, zcur, d, dcell, tau
    !     integer                :: celli, cellj, cellk, iseed
    !     logical                :: dir(3),tmp

    !     xcur = xp + xmax
    !     ycur = yp + ymax
    !     zcur = zp + zmax

    !     celli = xcell
    !     cellj = ycell
    !     cellk = zcell

    !     taurun = 0.
    !     taucell = 0.

    !     d = 0.
    !     dcell = 0.

    !     dir = (/.FALSE., .FALSE., .FALSE./)
    !     do
    !         dcell = wall_dist(celli, cellj, cellk, xcur, ycur, zcur, dir)
    !         ! taucell = dcell * rhokap(celli,cellj,cellk)

    !         taurun = taurun + taucell
    !         d = d + dcell
    !         call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta)

    !         if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
    !             exit
    !         end if
    !     end do
    ! end subroutine taufind1


    ! subroutine tauquick(xcell, ycell, zcell, zp, delta, taurun)

    !     use constants,   only : nzg, zmax
    !     use iarray,      only : zface!,rhokap
    !     use opt_prop,    only : wavelength
    !     use photon_vars, only : nzp

    !     implicit none


    !     integer, intent(IN)  :: xcell, ycell, zcell
    !     real,    intent(IN)  :: zp, delta
    !     real,    intent(OUT) :: taurun

    !     integer :: cellk
    !     real    :: tau, dcell, taucell, zcur

    !     taurun = 0.

    !     zcur = zp + zmax

    !     do cellk = zcell, 1 -1

    !         dcell = (zface(cellk) - zcur)/nzp
    !         ! taucell = dcell * rhokap(xcell, ycell, cellk)
    !         taurun = taurun + taucell
    !         zcur = zface(cellk) + delta
    !     end do
    ! end subroutine tauquick


    ! subroutine peeling(xcell,ycell,zcell,delta)
   
    !     ! use iarray,      only : image, imaget, imagethg, imagep
    !     use constants,   only : PI, xmax, ymax, zmax, v, costim, sintim, cospim, sinpim, nbins, twopi
    !     use photon_vars, only : xp, yp, zp, nxp, nyp, nzp, phase
    !     use opt_prop,    only : hgg, g2, wavelength
    !     use taufind2

    !     implicit none


    !     real,    intent(IN)    :: delta
    !     integer, intent(INOUT) :: xcell, ycell, zcell
    !     real                   :: cosa, tau1, prob, xim, yim, bin_wid,xpold,ypold
    !     real                   :: nxpold, nypold, nzpold, hgfact,zpold, tau3, dist
    !     integer                :: binx, biny, xcellold, ycellold, zcellold


    !     nxpold = nxp
    !     nypold = nyp
    !     nzpold = nzp

    !     xcellold = xcell
    !     ycellold = ycell
    !     zcellold = zcell

    !     xpold = xp
    !     ypold = yp
    !     zpold = zp
    !     cosa = nxp*v(1) + nyp*v(2) + nzp*v(3)!angle of peeled off photon

    !     nxp = v(1)
    !     nyp = v(2)
    !     nzp = v(3)
    !     xim = yp*cospim - xp*sinpim
    !     yim = zp*sintim - yp*costim*sinpim - xp*costim*cospim

    !     ! call taufind1(xcell, ycell, zcell, delta, tau3)
    !     call tau2(xcell,ycell,zcell,delta, tau3, dist)

    !     bin_wid = 4.*xmax/Nbins

    !     binx = floor(xim/bin_wid)
    !     biny = floor(yim/bin_wid)


    !     hgfact = (1.-g2) / ((4.*pi)*(1.+g2-2.*hgg*cosa)**(1.5))

    !     prob =  cos(twopi*dist/wavelength) !* tau3 !* hgfact
    !     ! image(binx, biny) = image(binx, biny) + prob
    !     ! imagep(binx, biny) = imagep(binx, biny) + prob + phase
    !     ! imaget(binx, biny) = imaget(binx, biny) + prob*tau3
    !     ! imagethg(binx, biny) = imagethg(binx, biny) + prob*tau3*hgfact

    !     nxp = nxpold
    !     nyp = nypold
    !     nzp = nzpold

    !     xcell = xcellold 
    !     ycell = ycellold 
    !     zcell = zcellold 

    ! end subroutine peeling
end module inttau2