module utils

    implicit none

    !foreground colours
    character(len=2), parameter :: black   = '30', &
                                   red     = '31', &
                                   green   = '32', &
                                   yellow  = '33', &
                                   blue    = '34', &
                                   magenta = '35', &
                                   cyan    = '36', &
                                   white   = '37'

    !background colours
    character(len=2), parameter :: black_b   = '40', &
                                   red_b     = '41', &
                                   green_b   = '42', &
                                   yellow_b  = '43', &
                                   blue_b    = '44', &
                                   magenta_b = '45', &
                                   cyan_b    = '46', &
                                   white_b   = '47'

    !styles
    character(len=2), parameter :: bold          = '01', &
                                   italic        = '03', &
                                   underline     = '04', &
                                   inverse       = '07', &
                                   strikethrough = '09'

    !ANSI control characters                               
    character(len=2), parameter :: start = achar(27)//'['
    character(len=3), parameter :: end = '[0m'


    !functions to add colour to output via ANSI colour codes
    interface colour
        module procedure colour_char
        module procedure colour_int32
        module procedure colour_int64
        ! module procedure colour_real4
        module procedure colour_real8
    end interface

    !functions to turn variables into strings
    interface str
        module procedure str_I32
        module procedure str_I64
        module procedure str_Iarray
        ! module procedure str_R4
        module procedure str_R8
        module procedure str_R8array
        module procedure str_logicalarray
    end interface str

    !subroutines to swap variables
    interface swap
        module procedure swap_I
        ! module procedure swap_R4
        module procedure swap_R8
    end interface swap

    private
    public :: str, swap, colour, mem_free
    public :: bold, italic, underline, strikethrough, black, red, green, yellow, blue, magenta, cyan, white
    public :: black_b, red_b, green_b, yellow_b, blue_b, magenta_b, cyan_b, white_b

    contains

        subroutine swap_I(a, b)

            implicit none

            integer, intent(INOUT) :: a, b
            integer :: tmp

            tmp = a
            a = b
            b = tmp
        end subroutine swap_I


        ! subroutine swap_R4(a, b)

        !     implicit none

        !     real, intent(INOUT) :: a, b
        !     real :: tmp

        !     tmp = a
        !     a = b
        !     b = tmp
        ! end subroutine swap_R4


        subroutine swap_R8(a, b)

            implicit none

            double precision, intent(INOUT) :: a, b
            double precision :: tmp

            tmp = a
            a = b
            b = tmp
        end subroutine swap_R8


        function str_I32(i, len)

            use iso_fortran_env, only : Int32

            implicit none

            integer(int32),    intent(IN)    :: i
            integer, optional, intent(IN) :: len

            character(len=:), allocatable :: str_I32
            character(len=100) :: string
            integer            :: lentmp, lenuse

            write(string,'(I100.1)') I

            if(present(len))then
                lentmp = len_trim(adjustl(string))
                lenuse = len

                if(len >= lentmp)then
                    str_I32 = repeat("0", lenuse - lentmp)//trim(adjustl(string))
                else
                    str_I32 = trim(adjustl(string))
                    str_I32 = trim(adjustl(str_I32(:len)))                        
                end if
            else
                str_I32 = trim(adjustl(string))
            end if
        end function str_I32


        function str_I64(i, len)

            use iso_fortran_env, only : Int64

            implicit none

            integer(int64),    intent(IN)    :: i
            integer, optional, intent(IN) :: len

            character(len=:), allocatable :: str_I64
            character(len=100) :: string
            integer            :: lentmp, lenuse

            write(string,'(I100.1)') I

            if(present(len))then
                lentmp = len_trim(adjustl(string))
                lenuse = len

                if(len >= lentmp)then
                    str_I64 = repeat("0", lenuse - lentmp)//trim(adjustl(string))
                else
                    str_I64 = trim(adjustl(string))
                    str_I64 = trim(adjustl(str_I64(:len)))                        
                end if
            else
                str_I64 = trim(adjustl(string))
            end if
        end function str_I64


        function str_iarray(i)

            implicit none

            integer, intent(IN) :: i(:)

            character(len=:), allocatable :: str_iarray
            character(len=100) :: string
            integer :: j

            do j = 1, size(i)
                write(string,'(I100.1)') I(j)
                str_iarray = str_iarray//' '//trim(adjustl(string))
            end do
            
        end function str_iarray


        ! function str_R4(i)

        !     implicit none

        !     real, intent(IN) :: i

        !     character(len=:), allocatable :: str_R4
        !     character(len=100) :: string

        !     write(string,'(f100.8)') I

        !     str_R4 = trim(adjustl(string))
        ! end function str_r4


        function str_R8(i, len)

            implicit none

            double precision,  intent(IN) :: i
            integer, optional, intent(IN) :: len

            character(len=:), allocatable :: str_R8
            character(len=100) :: string

            write(string,'(f100.16)') I

            if(present(len))then
                str_R8 = trim(adjustl(string))
                str_R8 = trim(adjustl(str_R8(:len)))
            else
                str_R8 = trim(adjustl(string))
            end if
        end function str_R8


        function str_R8array(a)

            implicit none

            double precision, intent(IN) :: a(:)

            character(len=:), allocatable :: str_R8array
            character(len=100) :: string
            integer :: i

            do i = 1, size(a)
                write(string,'(f100.16)') a(i)
                str_R8array = str_R8array//' '//trim(adjustl(string))
            end do

        end function str_R8array


        function str_logicalarray(a)

            implicit none

            logical, intent(IN) :: a(:)

            character(len=:), allocatable :: str_logicalarray
            character(len=100) :: string
            integer :: i

            do i = 1, size(a)
                write(string,'(L1)') a(i)
                str_logicalarray = str_logicalarray//' '//trim(adjustl(string))
            end do

        end function str_logicalarray


        function colour_char(string, fmt1, fmt2, fmt3, fmt4, fmt5) result(colourised)

            implicit none

            character(*),           intent(IN) :: string
            character(*), optional, intent(IN) :: fmt1, fmt2, fmt3, fmt4, fmt5
            character(len=:), allocatable      :: colourised

            colourised = string

            if(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4) .and. present(fmt5))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//';'//fmt5//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2))then
                colourised = start//fmt1//';'//fmt2//'m'//string//achar(27)//end
            elseif(present(fmt1))then
                colourised = start//fmt1//'m'//string//achar(27)//end
            end if
        end function colour_char


        function colour_int32(inte, fmt1, fmt2, fmt3, fmt4, fmt5) result(colourised)

            use iso_fortran_env, only : int32

            implicit none

            integer(kind=int32),    intent(IN) :: inte
            character(*), optional, intent(IN) :: fmt1, fmt2, fmt3, fmt4, fmt5

            character(len=:), allocatable :: colourised, string
            character(len=50)             :: tmp

            write(tmp,'(I50.1)') inte
            string = trim(adjustl(tmp))
            colourised = trim(adjustl(string))

            if(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4) .and. present(fmt5))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//';'//fmt5//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2))then
                colourised = start//fmt1//';'//fmt2//'m'//string//achar(27)//end
            elseif(present(fmt1))then
                colourised = start//fmt1//'m'//string//achar(27)//end
            end if
        end function colour_int32




        function colour_int64(inte, fmt1, fmt2, fmt3, fmt4, fmt5) result(colourised)

            use iso_fortran_env, only : int64

            implicit none

            integer(kind=int64),    intent(IN) :: inte
            character(*), optional, intent(IN) :: fmt1, fmt2, fmt3, fmt4, fmt5

            character(len=:), allocatable :: colourised, string
            character(len=50)             :: tmp

            write(tmp,'(I50.1)') inte
            string = trim(adjustl(tmp))
            colourised = trim(adjustl(string))

            if(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4) .and. present(fmt5))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//';'//fmt5//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2))then
                colourised = start//fmt1//';'//fmt2//'m'//string//achar(27)//end
            elseif(present(fmt1))then
                colourised = start//fmt1//'m'//string//achar(27)//end
            end if
        end function colour_int64


        ! function colour_real4(inte, fmt1, fmt2, fmt3, fmt4, fmt5) result(colourised)

        !     implicit none

        !     real,                   intent(IN) :: inte
        !     character(*), optional, intent(IN) :: fmt1, fmt2, fmt3, fmt4, fmt5

        !     character(len=:), allocatable :: colourised, string
        !     character(len=50)             :: tmp

        !     write(tmp,'(F50.8)') inte
        !     string = trim(adjustl(tmp))
        !     colourised = trim(adjustl(string))

        !     if(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4) .and. present(fmt5))then
        !         colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//';'//fmt5//'m'//string//achar(27)//end
        !     elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4))then
        !         colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//'m'//string//achar(27)//end
        !     elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3))then
        !         colourised = start//fmt1//';'//fmt2//';'//fmt3//'m'//string//achar(27)//end
        !     elseif(present(fmt1) .and. present(fmt2))then
        !         colourised = start//fmt1//';'//fmt2//'m'//string//achar(27)//end
        !     elseif(present(fmt1))then
        !         colourised = start//fmt1//'m'//string//achar(27)//end
        !     end if
        ! end function colour_real4


        function colour_real8(inte, fmt1, fmt2, fmt3, fmt4, fmt5) result(colourised)

            implicit none

            double precision,       intent(IN) :: inte
            character(*), optional, intent(IN) :: fmt1, fmt2, fmt3, fmt4, fmt5

            character(len=:), allocatable :: colourised, string
            character(len=50)             :: tmp

            write(tmp,'(F50.16)') inte
            string = trim(adjustl(tmp))
            colourised = trim(adjustl(string))

            if(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4) .and. present(fmt5))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//';'//fmt5//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3) .and. present(fmt4))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//';'//fmt4//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2) .and. present(fmt3))then
                colourised = start//fmt1//';'//fmt2//';'//fmt3//'m'//string//achar(27)//end
            elseif(present(fmt1) .and. present(fmt2))then
                colourised = start//fmt1//';'//fmt2//'m'//string//achar(27)//end
            elseif(present(fmt1))then
                colourised = start//fmt1//'m'//string//achar(27)//end
            end if
        end function colour_real8


        function mem_free()
        ! func that returns amount of free memory
        ! should be portable across most linux systems
        ! returns memory available in b

            use iso_fortran_env, only : int64 !as numbers are large

            implicit none

            integer(int64) :: mem_free, available, pagecache, active, inactive, sreclaimable, freeram

            integer(int64)    :: i, low
            character(len=15) :: tmp
            integer           :: u, io


            open(newunit=u,file='/proc/zoneinfo',status='old')
            low = 0
            sreclaimable = 0
            inactive = 0
            freeram = 0
            active = 0
            do 
                read(u,*, iostat=io)tmp, i
                if(IS_IOSTAT_END(io))exit
                if(verify(trim(tmp), 'low') == 0)then
                    low = low + i
                end if
            end do
            close(u)

            open(newunit=u,file='/proc/meminfo',status='old')
            mem_free = 0
            do 
                read(u,*,iostat=io)tmp, i
                if(IS_IOSTAT_END(io))exit

                if(verify(trim(tmp), 'MemAvailable:') == 0)then
                    mem_free = i * 1024_int64
                    return!return early if MemAvailable availble on OS. Else do calculation for MemAvailable
                end if
                if(verify(trim(tmp), 'MemFree:') == 0)then
                    freeram = i
                end if
                if(trim(tmp) == 'Active(file):' )then
                    active = i
                end if
                if(verify(tmp, 'Inactive(file):') == 0)then
                    inactive = i
                end if
                if(trim(tmp) == 'SReclaimable:')then
                    sreclaimable = i
                end if
            end do
            close(u)

            !algorithm from kernal source
        !https://git.kernel.org/pub/scm/linux/kernel/git/torvalds/linux.git/commit/?id=34e431b0ae398fc54ea69ff85ec700722c9da773
            available = freeram - low
            pagecache = active + inactive
            pagecache = pagecache - min(pagecache/2, low)
            available = available + pagecache
            available = available + (sreclaimable - min(sreclaimable/2, low))
            mem_free = available * 1024_int64 !convert from Kib to b
        end function mem_free
end module utils