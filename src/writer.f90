MODULE writer_mod

implicit none
save

CONTAINS
   subroutine writer(xmax,ymax,zmax,nphotons, numproc)

   use constants, only : nxg,nyg,nzg,fileplace
   use iarray,    only : jmeanGLOBAL, intensityGLOBAL

   implicit none

   integer           :: nphotons, i, u, numproc
   real              :: xmax, ymax, zmax
   
   jmeanGLOBAL =jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

   inquire(iolength=i)jmeanGLOBAL

   open(newunit=u,file=trim(fileplace)//'jmean/jmean.dat',access='direct',status='REPLACE',form='unformatted',&
   recl=i)
   write(u,rec=1) jmeanGLOBAL
   close(u)

  inquire(iolength=i)intensityGLOBAL

   open(newunit=u,file=trim(fileplace)//'jmean/int.dat',access='direct',status='REPLACE',form='unformatted',&
   recl=i)
   write(u,rec=1) intensityGLOBAL
   close(u)

   open(newunit=u,file='slice.dat',status='REPLACE')
   do i = 1,  nxg
    write(u,*)intensityGLOBAL(i,:,125)
end do
   end subroutine writer
end MODULE writer_mod