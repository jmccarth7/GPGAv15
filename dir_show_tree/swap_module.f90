module swap_module

implicit none

INTERFACE swap
    MODULE PROCEDURE swap_r, masked_swap_rs
END INTERFACE


contains


SUBROUTINE swap_r(a,b)
real(kind=8), INTENT(INOUT) :: a,b
real(kind=8) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_r

SUBROUTINE masked_swap_rs(a,b,mask)
real(kind=8), INTENT(INOUT) :: a,b
LOGICAL, INTENT(IN) :: mask
real(kind=8) :: swp
if( mask ) then
    swp=a
    a=b
    b=swp
endif
END SUBROUTINE masked_swap_rs

end module swap_module
