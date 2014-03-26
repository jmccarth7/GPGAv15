!
! File:   Math_Node_Functions.f03
! Author: Dave Coulter [NASA Summer Intern under John R. Moisan]
!
! Created on August 2, 2013, 4:00 PM
!

module Math_Node_Functions

    interface
        function Compute(a,b)
            real(kind=8) :: Compute
            real(kind=8), intent(in) :: a,b
        end function
    end interface

    type function_pointer
        procedure(Compute), pointer, nopass :: f
    end type

    type(function_pointer), dimension(20):: math_funcs



    ! GP Functions: These functions will be used in the math_funcs array, which
    ! contains procedure pointers to procedures below to be used by the
    ! Tree_Math_Node type. Above each function is the array position that it will
    ! occupy.

    contains

    !-------------------------------------------------------

    ! math_funcs(1)
    ! Addition: a + b

    real(kind=8) function f_Add(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        f_Add = a + b
    end function


    !-------------------------------------------------------

    ! math_funcs(2)
    ! Subtraction: a - b
    real(kind=8) function f_Subtract(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        f_Subtract = a - b
    end function


    !-------------------------------------------------------

    ! math_funcs(3)
    ! Multiply: a * b
    real(kind=8) function f_Multiply(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        !write(*,*)'f3:  a, b ', a, b                                                                           
        if( isnan(a) .or. isnan(b) ) then                                                                           
            f_Multiply = 0.0D0  
            return                                                                                                  
        endif                           

        f_Multiply = a * b

        !write(6,'(A,3(1x,E24.16))') 'f_Mul:  a, b, a*b ', a, b, a*b
    end function


    !-------------------------------------------------------

    ! math_funcs(4)
    ! Protected Divide (only if b not equal to zero): a / b
    real(kind=8) function f_ProtectedDivide(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        !if( b .ne. 0.0D+0) then
        if( abs(b) >  0.0D+0) then
            f_ProtectedDivide = a / b
        else
            f_ProtectedDivide = 0.0D+0
        endif
    end function


    !-------------------------------------------------------

    ! math_funcs(5)
    ! Ivlev Grazing Function: (1 - e^-abs(a*b))

    real(kind=8) function f_IvlevGrazingFunction(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        real(kind=8) :: cff

        cff=abs(a*b)
        if( cff < 100.0d0 ) then
            f_IvlevGrazingFunction = 1.0D+0 - exp(-1.0D+0*cff)
        else
            f_IvlevGrazingFunction = 1.0D+0
        endif 
    end function


    !-------------------------------------------------------

    ! math_funcs(6)
    !orig  Michaelis-Menton Term (modified for Forward-Backward): (1 / (abs(a) + abs(b)))
    ! Michaelis-Menton Term (modified for Forward-Backward): (abs(b) / (abs(a) + abs(b)))
    real(kind=8) function f_MichealisMenton(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        real(kind=8) :: cff

        cff=abs(a)+abs(b)
        if( cff .gt. 0.0D+0) then
            f_MichealisMenton=abs(b)/cff
        else
            f_MichealisMenton=0.0D+0
        endif
    end function


    !-------------------------------------------------------

    ! math_funcs(7)
    ! Mayzaud-Poulet Grazing Function:  abs(a*b)*(1 - e^-abs(a*b))
    real(kind=8) function f_MayzaudPouletGrazingFunction(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        real(kind=8) :: cff

        cff=abs(a*b)
        if( cff < 100.0d0 ) then
            f_MayzaudPouletGrazingFunction = cff*( 1.0D+0 - exp(-1.0D+0*cff) )
        else
            f_MayzaudPouletGrazingFunction = cff
        endif 
    end function


    !-------------------------------------------------------

    ! math_funcs(8)
    ! Power: a ^ b
    real(kind=8) function f_Power(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if( isnan(a) .or. isnan(b) ) then
            !write(*,*)'f8:  a, b ', a, b
            f_Power = 0.0D0
            return
        endif 

        if( abs(a) <= 1.0D-99 ) then
            !write(*,*)'f8:  a, b ', a, b
            f_Power = 0.0d0
            return
        endif

        if( abs(b) <= 1.0D-99 ) then
            !write(*,*)'f8:  a, b ', a, b
            f_Power = 1.0d0
            return
        endif

        f_Power = abs(a)**b

        f_Power = min( f_Power, 1.0D+99 )
        f_Power = max( f_Power, 1.0D-99 )

        !write(*,*)'f8:  a, b, f_Power ', a, b, f_Power

    end function


    !-------------------------------------------------------

    ! math_funcs(9)
    ! EXP: exp(-abs(a*b))
    real(kind=8) function f_ExponentialDecay(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        real(kind=8) :: cff

        if( isnan(a) .or. isnan(b) ) then
            !write(*,*)'f9:  a, b ', a, b
            f_ExponentialDecay = 0.0D0
            return
        endif 

        cff=abs(a*b)
        if( cff < 100.0d0 ) then
            f_ExponentialDecay = exp(-1.0D+0*cff)
        else
            f_ExponentialDecay = 0.0D0
        endif 
        !write(*,*)'f9:  a, b, f_ExponentialDecay ', a, b, f_ExponentialDecay
    end function


    !-------------------------------------------------------

    ! math_funcs(10)
    ! Minimum: a or b, whichever is lower
    real(kind=8) function f_Minimize(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        f_Minimize = min(a,b)
    end function


    !-------------------------------------------------------

    ! math_funcs(11)
    ! Maximum: a or b, whichever is greater
    real(kind=8) function f_Maximize(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        f_Maximize = max(a,b)
    end function



    !-------------------------------------------------------

!orig     ! math_funcs(9)
!orig     ! Minimum: a or b, whichever is lower
!orig     real(kind=8) function f_Minimize(a, b)
!orig         implicit none
!orig         real(kind=8), intent(in) :: a,b
!orig 
!orig         f_Minimize = min(a,b)
!orig     end function

    !-------------------------------------------------------


!orig     ! math_funcs(10)
!orig     ! Maximum: a or b, whichever is greater
!orig     real(kind=8) function f_Maximize(a, b)
!orig         implicit none
!orig         real(kind=8), intent(in) :: a,b
!orig 
!orig         f_Maximize = max(a,b)
!orig     end function


    !-------------------------------------------------------

!orig     ! math_funcs(11)
!orig     ! EXP: exp(-*a*b)
!orig     real(kind=8) function f_ExponentialDecay(a, b)
!orig         implicit none
!orig         real(kind=8), intent(in) :: a,b
!orig         real(kind=8) :: cff
!orig 
!orig         cff=abs(a*b)
!orig         f_ExponentialDecay = exp(-1.0D+0*cff)
!orig     end function


    !-------------------------------------------------------

    ! math_funcs(12)
    ! IF a .ne. 0 THEN b ELSE 0
    real(kind=8) function f_IfThen(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if( a .ne. 0.D+0) then
            f_IfThen = b
        else
            f_IfThen = 0.D+0
        endif
    end function


    !-------------------------------------------------------

    ! math_funcs(13)
    ! IF a .GT. b THEN 1 ELSE 0
    real(kind=8) function f_IfGt(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if( a .GT. b) then
            f_IfGt = 1.D+0
        else
            f_IfGt = 0.D+0
        endif
    end function


    ! math_funcs(14)
    ! IF a .GE. b THEN 1 ELSE 0
    real(kind=8) function f_IfGte(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if( a .GE. b) then
            f_IfGte = 1.D+0
        else
            f_IfGte = 0.D+0
        endif
    end function


    !-------------------------------------------------------

    ! math_funcs(15)
    ! IF a .LT. b THEN 1 ELSE 0
    real(kind=8) function f_IfLt(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if( a .LT. b) then
            f_IfLt = 1.D+0
        else
            f_IfLt = 0.D+0
        endif
    end function


    !-------------------------------------------------------

    ! math_funcs(16)
    ! IF a .LE. b THEN 1 ELSE 0
    real(kind=8) function f_IfLte(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b

        if( a .LE. b) then
            f_IfLte = 1.D+0
        else
            f_IfLte = 0.D+0
        endif
    end function

    !-------------------------------------------------------

    ! math_funcs(17)
    ! EXP_LP: exp(a)
    real(kind=8) function f_ExponentialLeftPlus(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        !real(kind=8) :: cff

        !cff=abs(a)
        f_ExponentialLeftPlus = exp( a )

        f_ExponentialLeftPlus = min( f_ExponentialLeftPlus, 1.0D+99 )
        f_ExponentialLeftPlus = max( f_ExponentialLeftPlus, 1.0D-99 )
    end function


    !-------------------------------------------------------

    ! math_funcs(18)
    ! EXP_RP: exp(b)
    real(kind=8) function f_ExponentialRightPlus(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        !real(kind=8) :: cff

        !cff=abs(a)
        f_ExponentialRightPlus = exp( b )
        f_ExponentialRightPlus = min( f_ExponentialRightPlus, 1.0D+99 )
        f_ExponentialRightPlus = max( f_ExponentialRightPlus, 1.0D-99 )
    end function


    !-------------------------------------------------------


    ! math_funcs(19)
    ! EXP_LM: exp(-a)
    real(kind=8) function f_ExponentialLeftMinus(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        !real(kind=8) :: cff

        !cff=abs(a)
        f_ExponentialLeftMinus = exp( -1.0d0 * a )
        f_ExponentialLeftMinus = min( f_ExponentialLeftMinus, 1.0D+99 )
        f_ExponentialLeftMinus = max( f_ExponentialLeftMinus, 1.0D-99 )
    end function


    !-------------------------------------------------------

    ! math_funcs(20)
    ! EXP_RM: exp(b)
    real(kind=8) function f_ExponentialRightMinus(a, b)
        implicit none
        real(kind=8), intent(in) :: a,b
        !real(kind=8) :: cff

        !cff=abs(a)
        f_ExponentialRightMinus = exp( -1.0d0 * b )
        f_ExponentialRightMinus = min( f_ExponentialRightMinus, 1.0D+99 )
        f_ExponentialRightMinus = max( f_ExponentialRightMinus, 1.0D-99 )
    end function


    !-------------------------------------------------------

    ! math_funcs(21)

    function G_Lower() result(n1)
    
        use Tree_Node_Factory_module
        use Fasham_Variables_module
        use GP_variables_module
        
        implicit none
        type(Tree_Node), pointer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, &
                                    n12, n13, n14, n15, n16, n17, n18, n19, &
                                    n24, n25, n28, n29, n30, n31, n32, n33, &
                                    n34, n35, n56, n57, n60, n61
        
        write(6,'(A)') 'in G_Lower '

        n61 => GetParameterNode(2.D+0)
        n60 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS)
        n57 => GetParameterNode(2.D+0)
        n56 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA)
        n35 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA)
        n34 => GetParameterNode(p2)
        n33 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
        n32 => GetParameterNode(p1)
        n31 => GetParameterNode(p3)
        n30 => GetMathNode(Power, n60, n61)
        n29 => GetParameterNode(p2)
        n28 => GetMathNode(Power, n56, n57)
        n25 => GetParameterNode(2.D+0)
        n24 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
        n19 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS)
        n18 => GetParameterNode(p3)
        n17 => GetMathNode(Multiply, n34, n35)
        n16 => GetMathNode(Multiply, n32, n33)
        n15 => GetMathNode(Multiply, n30, n31)
        n14 => GetMathNode(Multiply, n28, n29)
        n13 => GetParameterNode(p1)
        n12 => GetMathNode(Power, n24, n25)
        n9 => GetMathNode(Multiply, n18, n19)
        n8 => GetMathNode(Add, n16, n17)
        n7 => GetMathNode(Add, n14, n15)
        n6 => GetMathNode(Multiply, n12, n13)
        n5 => GetParameterNode(aK3)
        n4 => GetMathNode(Add, n8, n9)
        n3 => GetMathNode(Add, n6, n7)
        n2 => GetMathNode(Multiply, n4, n5)
        n1 => GetMathNode(Add, n2, n3)
        
    end function G_Lower



    !-------------------------------------------------------------------------------------------
    
    ! math_funcs(22)
    function G1() result(n1)
    
        use Tree_Node_Factory_module
        use Fasham_Variables_module
        use GP_variables_module
        use Fasham_Tree_Interfaces, only : G_Lower
        
        implicit none
        type (Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9, n10, n11, n16, n17
        
        write(6,'(A)') 'in G1 '
        n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
        n16 => GetParameterNode(g)
        n11 => GetParameterNode(2.D+0)
        n10 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
        n9 => GetParameterNode(p1)
        n8 => GetMathNode(Multiply, n16, n17)
        n5 => GetMathNode(Power, n10, n11)
        n4 => GetMathNode(Multiply, n8, n9)
        n3 => G_Lower()
        n2 => GetMathNode(Multiply, n4, n5)
        n1 => GetMathNode(ProtectedDivide, n2, n3)
        
    end function G1
    
    !-------------------------------------------------------------------------------------------
    
    
    ! math_funcs(23)
    
    function G2() result(n1)
    
        use Tree_Node_Factory_module
        use Fasham_Variables_module
        use GP_variables_module
        use Fasham_Tree_Interfaces, only : G_Lower
        
        implicit none
        type (Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9, n10, n11, n16, n17
        
        write(6,'(A)') 'in G2'
        n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
        n16 => GetParameterNode(g)
        n11 => GetParameterNode(2.D+0)
        n10 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA)
        n9 => GetParameterNode(p2)
        n8 => GetMathNode(Multiply, n16, n17)
        n5 => GetMathNode(Power, n10, n11)
        n4 => GetMathNode(Multiply, n8, n9)
        n3 => G_Lower()
        n2 => GetMathNode(Multiply, n4, n5)
        n1 => GetMathNode(ProtectedDivide, n2, n3)
        
    end function G2
    
    
    !-------------------------------------------------------------------------------------------
    
    ! math_funcs(24)
    
    function G3() result(n1)
    
        use Tree_Node_Factory_module
        use Fasham_Variables_module
        use GP_variables_module
        use Fasham_Tree_Interfaces, only : G_Lower
        
        implicit none
        type (Tree_Node), pointer :: n1, n2, n3, n4, n5, n8, n9, n10, n11, n16, n17
        
        write(6,'(A)') 'in G3'
        n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
        n16 => GetParameterNode(g)
        n11 => GetParameterNode(2.D+0)
        n10 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS)
        n9 => GetParameterNode(p3)
        n8 => GetMathNode(Multiply, n16, n17)
        n5 => GetMathNode(Power, n10, n11)
        n4 => GetMathNode(Multiply, n8, n9)
        n3 => G_Lower()
        n2 => GetMathNode(Multiply, n4, n5)
        n1 => GetMathNode(ProtectedDivide, n2, n3)
        
    end function G3




end module Math_Node_Functions
