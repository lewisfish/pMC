module iarray
!
!  Contains all array var names.
!
    implicit none

    real, allocatable :: xface(:), yface(:), zface(:)
    ! real, allocatable :: rhokap(:,:,:)
    ! real, allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:)
    
    ! complex, allocatable :: image(:,:), imageGLOBAL(:,:)
    complex, allocatable :: phasor(:,:,:), phasorGLOBAL(:,:,:)
end module iarray
