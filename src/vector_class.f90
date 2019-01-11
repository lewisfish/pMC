Module vector_class

    type :: vector
        real :: x, y, z
        contains

        procedure :: magnitude       => magnitude_fn
        procedure :: print           => print_sub
        generic   :: operator(.dot.) => vec_dot
        generic   :: operator(/)     => vec_div_scal
        generic   :: operator(*)     => vec_mult_vec, vec_mult_scal, scal_mult_vec
        generic   :: operator(+)     => vec_add_vec, vec_add_scal, scal_add_vec
        generic   :: operator(-)     => vec_minus_vec

        procedure, pass(a), private :: vec_dot

        procedure, pass(a), private :: vec_div_scal

        procedure, pass(a), private :: vec_mult_vec
        procedure, pass(a), private :: vec_mult_scal
        procedure, pass(b), private :: scal_mult_vec

        procedure, pass(a), private :: vec_add_vec
        procedure, pass(a), private :: vec_add_scal
        procedure, pass(b), private :: scal_add_vec

        procedure, pass(a), private :: vec_minus_vec

    end type vector

    private
    public :: magnitude, vector, print

    contains

        type(vector) function vec_minus_vec(a, b)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b

            vec_minus_vec = vector(a%x - b%x, a%y - b%y, a%z - b%z)

        end function vec_minus_vec


        type(vector) function vec_add_scal(a, b)

            implicit none

            class(vector), intent(IN) :: a
            real,          intent(IN) :: b

            vec_add_scal = vector(a%x + b, a%y + b, a%z + b)

        end function vec_add_scal


        type(vector) function scal_add_vec(a, b)

            implicit none

            class(vector), intent(IN) :: b
            real,          intent(IN) :: a

            scal_add_vec = vector(b%x + a, b%y + a, b%z + a)

        end function scal_add_vec


        type(vector) function vec_add_vec(a, b)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b

            vec_add_vec = vector(a%x + b%x, a%y + b%y, a%z + b%z)

        end function vec_add_vec


        elemental function vec_dot(a, b) result (dot)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b
            real :: dot

            dot = (a%x * b%x) + (a%y * b%y) + (a%z * b%z)

        end function vec_dot


        type(vector) function vec_mult_vec(a, b)

            implicit none

            class(vector), intent(IN) :: a
            type(vector),  intent(IN) :: b

            vec_mult_vec = vector(a%x * b%x, a%y * b%y, a%z * b%z)

        end function vec_mult_vec


        type(vector) function vec_mult_scal(a, b)

            implicit none

            class(vector), intent(IN) :: a
            real,          intent(IN) :: b

            vec_mult_scal = vector(a%x * b, a%y * b, a%z * b)

        end function vec_mult_scal


        type(vector) function scal_mult_vec(a, b)

            implicit none

            class(vector), intent(IN) :: b
            real,          intent(IN) :: a

            scal_mult_vec = vector(a * b%x, a * b%y, a * b%z)

        end function scal_mult_vec


        type(vector) function vec_div_scal(a, b)

            implicit none

            class(vector), intent(IN) :: a
            real,         intent(IN) :: b

            vec_div_scal = vector(a%x / b, a%y / b, a%z / b)

        end function vec_div_scal


        type(vector) function magnitude_fn(this)

            implicit none

            class(vector) :: this

            real :: tmp

            tmp = sqrt(this%x**2 + this%y**2 + this%z**2)
            magnitude_fn = this / tmp

        end function magnitude_fn


        subroutine print_sub(this)

            implicit none

            class(vector) :: this

                print*,this%x, this%y, this%z

        end subroutine
end Module vector_class