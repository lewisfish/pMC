Module precision

    use iso_fortran_env

    implicit none

    public

#ifdef r8p
    integer, parameter :: RP = real64
#else
    integer, parameter :: RP = real32
#endif

#ifdef i8p
    integer, parameter :: IP = INT64
#else
    integer, parameter :: IP = INT32
#endif
end module precision