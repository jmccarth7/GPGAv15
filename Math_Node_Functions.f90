!
! File:   Math_Node_Functions.f03
! Author: Dave Coulter [NASA Summer Intern under John R. Moisan]
!
! Created on August 2, 2013, 4:00 PM
!

module Math_Node_Functions

    interface
        function Compute (a,b)
            real(kind=8) :: Compute
            real(kind=8), intent(in) :: a,b
        end function
    end interface

    type function_pointer
        procedure(Compute), pointer, nopass :: f
    end type

    type(function_pointer), dimension(16):: math_funcs



    ! GP Functions: These functions will be used in the math_funcs array, which
    ! contains procedure pointers to procedures below to be used by the
    ! Tree_Math_Node type. Above each function is the array position that it will
    ! occupy.

    contains

    ! math_funcs(1)
    ! Addition: a + b

    real(kind=8) function f_Add(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        f_Add = a + b
    end function

    ! math_funcs(2)
    ! Subtraction: a - b
    real(kind=8) function f_Subtract(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        f_Subtract = a - b
    end function

    ! math_funcs(3)
    ! Multiply: a * b
    real(kind=8) function f_Multiply(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        f_Multiply = a * b

        !write(6,'(A,3(1x,E15.7))') 'f_Mul:  a, b, a*b ', a, b, a*b
    end function

    ! math_funcs(4)
    ! Protected Divide (only if b not equal to zero): a / b
    real(kind=8) function f_ProtectedDivide(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if (b .ne. 0.0D+0) then
            f_ProtectedDivide = a / b
        else
            f_ProtectedDivide = 0.0D+0
        endif
    end function

    ! math_funcs(5)
    ! Ivlev Grazing Function: (1 - e^-abs(a*b))

    real(kind=8) function f_IvlevGrazingFunction(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        real(kind=8) :: cff

        cff=abs(a*b)
        f_IvlevGrazingFunction = 1.0D+0 - exp(-1.0D+0*cff)
    end function

    ! math_funcs(6)
    !orig  Michaelis-Menton Term (modified for Forward-Backward): (1 / (abs(a) + abs(b)))
    ! Michaelis-Menton Term (modified for Forward-Backward): (abs(b) / (abs(a) + abs(b)))
    real(kind=8) function f_MichealisMenton(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        real(kind=8) :: cff

        cff=abs(a)+abs(b)
        if (cff .gt. 0.0D+0) then
            f_MichealisMenton=abs(b)/cff
        else
            f_MichealisMenton=0.0D+0
        endif
    end function

    ! math_funcs(7)
    ! Mayzaud-Poulet Grazing Function:  abs(a*b)*(1 - e^-abs(a*b))
    real(kind=8) function f_MayzaudPouletGrazingFunction(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        real(kind=8) :: cff

        cff=abs(a*b)
        f_MayzaudPouletGrazingFunction = cff*(1.0D+0 - exp(-1.0D+0*cff))
    end function

    ! math_funcs(8)
    ! Power: a ^ b
    real(kind=8) function f_Power(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        f_Power = a ** b
    end function

    ! math_funcs(9)
    ! Minimum: a or b, whichever is lower
    real(kind=8) function f_Minimize(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        f_Minimize = min(a,b)
    end function

    ! math_funcs(10)
    ! Maximum: a or b, whichever is greater
    real(kind=8) function f_Maximize(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        f_Maximize = max(a,b)
    end function

    ! math_funcs(11)
    ! EXP: exp(-*a*b)
    real(kind=8) function f_ExponentialDecay(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        real(kind=8) :: cff

        cff=abs(a*b)
        f_ExponentialDecay = exp(-1.0D+0*cff)
    end function

    ! math_funcs(12)
    ! IF a .ne. 0 THEN b ELSE 0
    real(kind=8) function f_IfThen(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if (a .ne. 0.D+0) then
            f_IfThen=b
        else
            f_IfThen=0.D+0
        endif
    end function

    ! math_funcs(13)
    ! IF a .GT. b THEN 1 ELSE 0
    real(kind=8) function f_IfGt(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if (a .GT. b) then
            f_IfGt=1.D+0
        else
            f_IfGt=0.D+0
        endif
    end function

    ! math_funcs(14)
    ! IF a .GE. b THEN 1 ELSE 0
    real(kind=8) function f_IfGte(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if (a .GE. b) then
            f_IfGte=1.D+0
        else
            f_IfGte=0.D+0
        endif
    end function

    ! math_funcs(15)
    ! IF a .LT. b THEN 1 ELSE 0
    real(kind=8) function f_IfLt(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if (a .LT. b) then
            f_IfLt=1.D+0
        else
            f_IfLt=0.D+0
        endif
    end function

    ! math_funcs(16)
    ! IF a .LE. b THEN 1 ELSE 0
    real(kind=8) function f_IfLte(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if (a .LE. b) then
            f_IfLte=1.D+0
        else
            f_IfLte=0.D+0
        endif
    end function


end module Math_Node_Functions
