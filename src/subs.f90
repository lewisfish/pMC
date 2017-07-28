MODULE subs

implicit none
save

CONTAINS
   SUBROUTINE directory
!  subroutine defines vars to hold paths to various folders   
!   
!   
   use constants,only : cwd,homedir,fileplace,resdir
   
#ifdef intel
  use ifport
#endif
   
   implicit none

   integer :: io

   !get current working directory
#ifdef intel
  io = getcwd(cwd)
#else
   CALL getcwd(cwd)
#endif

   !get 'home' dir from cwd
   homedir=trim(cwd(1:len(trim(cwd))-3))
   !get data dir
   fileplace=trim(homedir)//'data/'
   
   !checks to see if data folder exists, if not creates it.
#ifdef intel
  io = chdir(fileplace)
#else
   call chdir(fileplace,io)
#endif
   if(io.ne.0)then
      print*,'data directory does not exist...'
      print*, 'creating directory...'
#ifdef intel
     io = system("mkdir "//fileplace)
     io = chdir(fileplace)
     io = system("mkdir jmean/")
     io = system("mkdir deposit/")
     io = system("mkdir im/")
#else
      call system("mkdir "//fileplace)
      call chdir(fileplace,io)
      call system("mkdir jmean/")
      call system("mkdir deposit/")
      call system("mkdir im/")
      print*, 'created directory ',trim(fileplace)
#endif
   end if
   
   !get res dir
   resdir=trim(homedir)//'res/'
   
   end SUBROUTINE directory
   
   SUBROUTINE zarray
   
   use iarray
   
   !sets all arrays to zero
   implicit none
   
   xface = 0.
   yface = 0.
   zface = 0.
   rhokap = 0.
   ! jmean = 0.
   ! intensity = 0.
   ! image = 0.
   ! imaget = 0.
   ! imagep = 0.
   ! imagethg = 0.
   phasor = 0.

   end SUBROUTINE zarray

   SUBROUTINE alloc_array(numproc, id)
!  subroutine allocates allocatable arrays
!   
!   
      use iso_fortran_env, only : int64
      use utils,           only : str, mem_free
      use constants,       only : nxg, nyg, nzg, nbins
      use iarray

      implicit none

      integer, intent(IN) :: numproc, id

      integer(int64) :: limit, cnt, i
      
      
      limit = mem_free()
      cnt = 0_int64

      allocate(xface(nxg+1))
      inquire(iolength=i)xface(:)
      call chck_mem(cnt, i, limit, 'xface', numproc)

      allocate(yface(nyg+1))
      inquire(iolength=i)yface
      call chck_mem(cnt, i, limit, 'yface', numproc)

      allocate(zface(nzg+1))
      inquire(iolength=i)zface
      call chck_mem(cnt, i, limit, 'zface', numproc)

      allocate(rhokap(nxg,nyg,nzg))
      inquire(iolength=i)rhokap
      call chck_mem(cnt, i, limit, 'rhokap', numproc)

      ! allocate(intensity(nxg,nyg,nzg))
      ! inquire(iolength=i)intensity
      ! call chck_mem(cnt, i, limit, 'intensity', numproc)
      
      ! allocate(image(-((Nbins-1)/2):((Nbins-1)/2), -((Nbins-1)/2):((Nbins-1)/2)))
      ! inquire(iolength=i)image
      ! call chck_mem(cnt, i, limit, 'image', numproc)

      ! allocate(imaget(-((Nbins-1)/2):((Nbins-1)/2), -((Nbins-1)/2):((Nbins-1)/2)))
      ! inquire(iolength=i)imaget
      ! call chck_mem(cnt, i, limit, 'imaget', numproc)

      ! allocate(imagethg(-((Nbins-1)/2):((Nbins-1)/2), -((Nbins-1)/2):((Nbins-1)/2)))
      ! inquire(iolength=i)imagethg
      ! call chck_mem(cnt, i, limit, 'imagethg', numproc)

      !       allocate(imagep(-((Nbins-1)/2):((Nbins-1)/2), -((Nbins-1)/2):((Nbins-1)/2)))
      ! inquire(iolength=i)imagep
      ! call chck_mem(cnt, i, limit, 'imagep', numproc)

      allocate(phasor(nxg,nyg,nzg))
      inquire(iolength=i)phasor
      call chck_mem(cnt, i, limit, 'phasor', numproc)

      ! allocate(jmean(nxg,nyg,nzg))
      ! inquire(iolength=i)jmean
      ! call chck_mem(cnt, i, limit, 'jmean', numproc)
      
      if(id == 0)print'(A,1X,F5.2,A)','allocated:',dble(cnt)/dble(limit)*100.d0,' % of total RAM'
   end subroutine alloc_array


   subroutine chck_mem(cur, new, limit, name, numproc)

      use iso_fortran_env, only : int64
      use utils,           only : str

      implicit none

      integer(int64), intent(IN)    :: new, limit
      integer(int64), intent(INOUT) :: cur 
      integer,        intent(IN)    :: numproc
      character(*),   intent(IN)    :: name

      integer :: error

      cur = cur + new * numproc
      if(cur > limit)then
         print*,'Need '//str(cur-limit)//' more memory to run. '//name
         call mpi_finalize(error)
         stop
      end if
   end subroutine chck_mem

end MODULE subs
