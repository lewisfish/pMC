MODULE iarray
!
!  Contains all array var names.
!
implicit none
save

real, allocatable :: xface(:), yface(:), zface(:)
real, allocatable :: rhokap(:,:,:)
! real, allocatable :: image(:,:), imageGLOBAL(:,:)
! real, allocatable :: imagep(:,:), imagepGLOBAL(:,:)

! real, allocatable :: imaget(:,:), imagetGLOBAL(:,:)
! real, allocatable :: imagethg(:,:), imagethgGLOBAL(:,:)

! real, allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:)
! real, allocatable :: intensity(:,:,:), intensityGLOBAL(:,:,:)
real, allocatable :: phasor(:,:,:), phasorGLOBAL(:,:,:)
end MODULE iarray
