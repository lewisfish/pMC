module gridset_mod

    implicit none

    contains
        subroutine gridset(id)

            use constants, only : nxg, nyg, nzg, xmax, ymax, zmax
            use iarray,    only : xface, yface, zface!, rhokap
            use ch_opt,    only : init_opt4
            ! use opt_prop,  only : kappa

            implicit none

            integer, intent(IN) :: id
            integer :: i, j, k
            real    :: x, y, z

            if(id == 0)then
                print*, ' '
                print *, 'Setting up density grid....'
            end if
            !**********  Linear Cartesian grid. Set up grid faces ****************
            do i = 1, nxg+1
                xface(i)=(i-1)*2.*xmax/nxg
            end do

            do i = 1, nyg+1
                yface(i)=(i-1)*2.*ymax /nyg
            end do

            do i = 1, nzg+1
                zface(i)=(i-1)*2.*zmax/nzg
            end do

            call init_opt4()
            !**************  Loop through x, y, and z to set up grid density and refractive index grid.  ****
            do i = 1, nxg
                do j = 1, nyg
                    do k = 1, nzg
                        x = xface(i)-xmax+xmax/nxg
                        y = yface(j)-ymax+ymax/nyg
                        z = zface(k)-zmax+zmax/nzg
                        !***********Call density setup subroutine
                        
                        ! if(i >= 45 .and. i<= 55)then
                        !    if(j >= 45 .and. j<= 55)then
                        !       if(k >= 150 .and. k<= 175)then
                        !          rhokap(i,j,k)=1.d20
                        !       else
                        !          rhokap(i,j,k)=kappa
                        !       end if
                        !    else
                        !       rhokap(i,j,k)=kappa
                        !    end if
                        ! else
                        !    rhokap(i,j,k)=kappa
                        ! end if
                    end do
                end do
            end do
        end subroutine gridset
end module gridset_mod
