MODULE gridset_mod

implicit none
save

CONTAINS
   subroutine gridset(id)

   use density_mod
   use constants, only : nxg, nyg, nzg, xmax, ymax, zmax
   use iarray,    only : rhokap, xface, yface, zface
   use opt_prop,  only : kappa
   use ch_opt

   implicit none

   integer, intent(IN) :: id
   integer :: i, j, k
   real    :: x, y, z, taueq1, taupole1, taueq2, taupole2

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
   call init_opt4
   !**************  Loop through x, y, and z to set up grid density and refractive index grid.  ****
   do i = 1, nxg
      do j = 1, nyg
         do k = 1, nzg
            x = xface(i)-xmax+xmax/nxg
            y = yface(j)-ymax+ymax/nyg
            z = zface(k)-zmax+zmax/nzg
!***********Call density setup subroutine
            
               ! if(i >= 95 .and. i<= 105)then
               !    if(j >= 95 .and. j<= 105)then
               !       if(k >= 150 .and. k<= 160)then
               !          rhokap(i,j,k)=100000.*kappa
               !       else
               !          rhokap(i,j,k)=kappa
               !       end if
               !    else
               !       rhokap(i,j,k)=kappa
               !    end if
               ! else
                  rhokap(i,j,k)=kappa
               ! end if
         end do
      end do
   end do

   !****************** Calculate equatorial and polar optical depths ****
   taueq1=0.
   taupole1=0.
   taueq2=0.
   taupole2=0.
   do i=1,nxg
      taueq1=taueq1+rhokap(i,nyg/2,nzg/2)
   end do
   do i=1,nzg
      taupole1=taupole1+rhokap(nxg/2,nyg/2,i)
   end do
   taueq1=taueq1*2.*xmax/nxg
   taupole1=taupole1*2.*zmax/nzg
   if(id == 0)then
      print'(A,F9.5,A,F9.5)',' taueq1 = ',taueq1,'  taupole1 = ',taupole1
   end if
   
   end subroutine gridset
end MODULE gridset_mod
