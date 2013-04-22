subroutine GP_Check_Terminals(i_Error)

! This subroutine looks through a specific GP_Individual_Node_Type array
! for nodes that do not correctly set terminals.  If the terminals are not correctly set it returns i_Error=1
! else i_Error=0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters
use GA_Parameters
use GP_Variables
use GA_Variables

implicit none
integer(kind=4) :: i_Error
integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node
integer(kind=4) :: i_Level
integer(kind=4) :: i_Function
integer(kind=4) :: i_Node_Left
integer(kind=4) :: i_Node_Right

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

i_Error=0

do  i_Tree=1,n_Trees

    do  i_Level=1,n_Levels-1

        i_Function=(2**(i_Level-1))-1       ! calculated the function number at the right end of the upper level

        do  i_Node=2**i_Level,(2*(2**i_Level))-1,2   ! run through each function at the level

            i_Function=i_Function+1                  ! sets the 'function' node's index

            i_Node_Left=i_Node                       ! sets the 'left terminal' node's index; 
                                                     ! i_node_left=i_function*2 would also work

            i_Node_Right=i_Node+1                    ! sets the 'right terminal' node's index; 
                                                     ! i_node_right=(i_function*2)+1 would also work

            if( GP_Individual_Node_Type(i_Function,i_Tree) .gt. 0) then  ! It is a function node

                if( GP_Individual_Node_Type(i_Node_Left,i_Tree) .lt. -n_CODE_Equations ) then
                    write(*,*) 'L: ',i_Node_Left,GP_Individual_Node_Type(i_Function,i_Tree),&
                                                 GP_Individual_Node_Type(i_Node_Left,i_Tree)
                    i_Error=1
                endif ! GP_Individual_Node_Type(i_Node_Left,i_Tree) .lt. -n_CODE_Equations 

                if( GP_Individual_Node_Type(i_Node_Right,i_Tree) .lt. -n_CODE_Equations) then
                  write(*,*) 'R: ',i_Node_Right,GP_Individual_Node_Type(i_Function,i_Tree),&
                                                GP_Individual_Node_Type(i_Node_Right,i_Tree)
                  i_Error=1
                endif !   GP_Individual_Node_Type(i_Node_Right,i_Tree) .lt. -n_CODE_Equations
            endif !  GP_Individual_Node_Type(i_Function,i_Tree) .gt. 0

        enddo ! i_node

    enddo ! i_level 

enddo ! i_tree

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Check_Terminals
!23456789012345678901234567890123456789012345678901234567890123456789012
