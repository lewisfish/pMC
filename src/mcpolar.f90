program mcpolar

use mpi

!shared data
use utils, only : green, blue, bold, colour, red, white_b, black, str
use constants
use photon_vars
use iarray
use opt_prop

!subroutines
use subs
use gridset_mod
use sourceph_mod
use inttau2
use ch_opt
use stokes_mod
use writer_mod

implicit none

integer             :: nphotons, iseed, j, xcell, ycell, zcell
logical             :: tflag
double precision    :: nscatt
real                :: ran, delta, start, finish, ran2

!mpi variables
integer :: id, error, numproc
real    :: nscattGLOBAL

!set directory paths
call directory

!allocate and set arrays to 0
call alloc_array
call zarray


call MPI_init(error)

call MPI_Comm_size(MPI_COMM_WORLD, numproc, error)

call MPI_Comm_rank(MPI_COMM_WORLD, id, error)

!**** Read in parameters from the file input.params
open(10,file=trim(resdir)//'input.params',status='old')
   read(10,*) nphotons
   read(10,*) xmax
   read(10,*) ymax
   read(10,*) zmax
   read(10,*) n1
   read(10,*) n2
   read(10,*) beam
   close(10)

! set seed for rnd generator. id to change seed for each process
iseed = -95648324 + id
iseed = -abs(iseed)  ! Random number seed must be negative for ran2

call init_opt4
wavelength = 1064.e-9
lambdainpx = wavelength/(2.*xmax/nxg)


if(id == 0)then
   print*, ''      
   print*,'# of photons to run',nphotons*numproc
end if

!***** Set up density grid *******************************************
call gridset(id)

!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
   delta = 1.e-8*(2.*zmax/nzg)
!   deltay = 1.e-5*(2.*ymax/nyg)          !!!!!!!!1! impliment!!!!!!!!!!!!!!!!!!!!!!!!!
!   deltaz = 1.e-5*(2.*zmax/nzg)
nscatt=0

call cpu_time(start)

!loop over photons 
call MPI_Barrier(MPI_COMM_WORLD, error)
print*,'Photons now running on core: ',colour(id, green)
do j=1,nphotons
  
   call init_opt4

   tflag=.FALSE.

   if(mod(j,1000) == 0)then
      print *, colour(j, blue, bold),' scattered photons completed on core: ',colour(id, str(30+mod(id,7)), bold)
   end if
    
!***** Release photon from point source *******************************
   call sourceph(xcell,ycell,zcell,iseed)

!****** Find scattering location

      call tauint1(xcell,ycell,zcell,tflag,iseed,delta)
            
!******** Photon scatters in grid until it exits (tflag=TRUE) 
   do while(tflag.eqv..FALSE.)
      ran = ran2(iseed)
      
      if(ran < albedo)then!interacts with tissue
            call stokes(iseed)
            nscatt = nscatt + 1        
         else
            tflag=.true.
            exit
      end if

!************ Find next scattering location

       call tauint1(xcell,ycell,zcell,tflag,iseed,delta)
       
   end do
end do      ! end loop over nph photons

deallocate(xface,yface,zface,rhokap)

call cpu_time(finish)
if(finish-start.ge.60.)then
 print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
else
      print*, 'time taken ~',colour(floor(finish-start/60.),red, bold),'s'
end if

allocate(jmeanGLOBAL(nxg,nyg,nzg))
jmeanGLOBAL = 0.
call MPI_REDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
deallocate(jmean)
if(id /= 0)deallocate(jmeanGLOBAL)
call MPI_BARRIER(MPI_COMM_WORLD, error)


allocate(intensityGLOBAL(nxg,nyg,nzg))
intensityGLOBAL = 0.
call MPI_REDUCE(intensity, intensityGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
deallocate(intensity)
if(id /= 0)deallocate(intensityGLOBAL)
call MPI_BARRIER(MPI_COMM_WORLD, error)


allocate(phasorGLOBAL(nxg,nyg,nzg))
phasorGLOBAL = 0.
call MPI_REDUCE(phasor, phasorGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
deallocate(phasor)
if(id /= 0)deallocate(phasorGLOBAL)
call MPI_BARRIER(MPI_COMM_WORLD, error)


call MPI_REDUCE(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)

call MPI_BARRIER(MPI_COMM_WORLD, error)

if(id == 0)then
   print*,'Average # of scatters per photon:',nscattGLOBAL/(nphotons*numproc)
   !write out files

   call writer(nphotons, numproc)
   print*,colour('write done',black,white_b,bold)
end if

call MPI_Finalize(error)
end program mcpolar
