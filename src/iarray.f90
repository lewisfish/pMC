module iarray
!
!  Contains all array var names.
!
    implicit none

    real, allocatable :: xface(:), yface(:), zface(:)
    real, allocatable :: rhokap(:,:,:)
    real, allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:)

    ! complex, allocatable :: imaget(:,:), imagetGLOBAL(:,:)
    complex, allocatable :: imageb(:,:), imagebGLOBAL(:,:)
    complex, allocatable :: phasor(:,:,:), phasorGLOBAL(:,:,:)
end module iarray
