module tree_module

use class_Tree_Node

implicit none
public :: insert, print_tree

!type, public :: node
!    integer :: value
!    type(node), pointer ::  left  => null()
!    type(node), pointer ::  right => null()
!end type node



!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

recursive subroutine insert(t, counter, number)
!recursive subroutine insert(t, number)




type(Tree_node), pointer :: t  ! A tree
!type(Tree_node_pointer) :: t  ! A tree
integer, intent(in) :: number
integer, intent(in) :: counter 
integer, parameter  :: null_value = -9999


!--------------------------------------------

!if(  number <= null_value ) return 

! If (sub)tree is empty, put number at root
! assume first number is not null value

if( .not. associated(t)) then

    allocate (t)
    t%node_type = number

    print *, counter, t%node_type

! Otherwise, insert into correct subtree

!elseif( number < t%value) then
else  !if( number <= null_value ) then

    if( mod( counter,2) == 0 )then

        if( number > null_value )then
            call insert(t%left, counter, number)
            !call insert(t%left,  number)
        endif 
    else
        if( number > null_value )then
            call insert(t%right,counter,  number)
            !call insert(t%right,  number)
        endif 
    endif

endif ! .not. assoc



end subroutine insert


!---------------------------------------------------------------------------


recursive subroutine print_tree(t)

! Print tree in infix order
type(Tree_node), pointer :: t  ! A tree
!type(Tree_node_pointer) :: t  ! A tree

!----------------------------------------------------

if( associated(t)) then
    print *, t%node_type            ! new
    call print_tree(t % left)
    !orig print *, t%value
    call print_tree(t % right)
endif

end subroutine print_tree



end module tree_module
