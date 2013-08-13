subroutine GP_Check_Terminals(i_Error)

! This subroutine looks through a specific GP_Individual_Node_Type array
! for nodes that do not correctly set terminals.

! If the terminals are not correctly set it returns i_Error=1

! else i_Error=0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module
use GP_Parameters_module
use GP_model_parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

integer(kind=4) :: i_Error

integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node
integer(kind=4) :: i_level
integer(kind=4) :: i_function
integer(kind=4) :: i_Node_left
integer(kind=4) :: i_Node_right

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

i_Error=0

do  i_Tree=1,n_Trees

    do  i_Level=1,n_Levels-1

        i_Function = 2**(i_Level-1) -1               ! calculated the function number
                                                     ! at the right end of the upper level

        do  i_Node=2**i_Level, 2*(2**i_Level)-1 , 2  ! run through each function at the level

            i_Function=i_Function+1                  ! sets the 'function' node's index

            i_Node_Left=i_Node                       ! sets the 'left terminal' node's index;
                                                     ! i_node_left=i_function*2 would also work

            i_Node_Right=i_Node+1                    ! sets the 'right terminal' node's index;
                                                     ! i_node_right=(i_function*2)+1 would also work



            !write(GP_print_unit,'(/A,6(1x,I6))') &
            !      'gct: i_Tree, i_Level, i_Function, i_Node, i_Node_Left, i_Node_Right ', &
            !            i_Tree, i_Level, i_Function, i_Node, i_Node_Left, i_Node_Right
            !write(GP_print_unit,'(A,6(1x,I6))') 'gct: n_CODE_Equations ', n_CODE_Equations
            !write(GP_print_unit,'(A,3(1x,I6)/)') &
            !      'gct: i_function, i_tree, GP_Individual_Node_Type(i_Function, i_Tree)',&
            !            i_function, i_tree, GP_Individual_Node_Type(i_Function, i_Tree)

            if( GP_Individual_Node_Type(i_Function,i_Tree) .gt. 0) then

                ! It is a function node if > 0

                ! check Left node

                write(GP_print_unit,'(A,3(1x,I6))') &
                      'gct: i_Node_Left, i_tree, GP_Individual_Node_Type(i_Node_Left,i_Tree) ',&
                            i_Node_Left, i_tree, GP_Individual_Node_Type(i_Node_Left,i_Tree)

                if( GP_Individual_Node_Type(i_Node_Left,i_Tree) .lt. -n_CODE_Equations ) then

                    if( myid == 0 )then
                        !write(GP_print_unit,'(A,3(1x,I6))') &
                        !      'gct: Left ', &
                        !      i_Node_Left,GP_Individual_Node_Type(i_Function, i_Tree),&
                        !                  GP_Individual_Node_Type(i_Node_Left,i_Tree)
                        write(GP_print_unit,'(A,4(1x,I6))') &
                              'gct: i_Node_Left, i_function, i_tree, &
                              &GP_Individual_Node_Type(i_Function, i_Tree)',&
                                    i_Node_Left, i_function, i_tree, &
                               GP_Individual_Node_Type(i_Function, i_Tree)
                        write(GP_print_unit,'(A,3(1x,I6))') &
                              'gct: i_Node_Left, i_tree, &
                              &GP_Individual_Node_Type(i_Node_Left,i_Tree) ',&
                                    i_Node_Left, i_tree, &
                               GP_Individual_Node_Type(i_Node_Left,i_Tree)

                        call print_entire_tree( )

                    endif ! myid == 0

                    i_Error=1


                endif ! GP_Individual_Node_Type(i_Node_Left,i_Tree) .lt. -n_CODE_Equations


                ! check Right node

                write(GP_print_unit,'(A,3(1x,I6))') &
                      'gct: i_Node_Right, i_tree, GP_Individual_Node_Type(i_Node_Right,i_Tree) ',&
                            i_Node_Right, i_tree, GP_Individual_Node_Type(i_Node_Right,i_Tree)


                if( GP_Individual_Node_Type(i_Node_Right,i_Tree) .lt. -n_CODE_Equations) then

                    if( myid == 0 )then
                        write(GP_print_unit,'(A,3(1x,I6))') &
                          'gct:Right ',&
                          i_Node_Right,GP_Individual_Node_Type(i_Function,  i_Tree),&
                                       GP_Individual_Node_Type(i_Node_Right,i_Tree)
                        write(GP_print_unit,'(A,4(1x,I6))') &
                          'gct: i_Node_Right, i_function, i_tree, &
                          &GP_Individual_Node_Type(i_Function, i_Tree)',&
                                i_Node_Right, i_function, i_tree, &
                           GP_Individual_Node_Type(i_Function, i_Tree)
                        write(GP_print_unit,'(A,3(1x,I6))') &
                          'gct: i_Node_Right, i_tree, &
                           &GP_Individual_Node_Type(i_Node_Right,i_Tree) ',&
                                i_Node_Right, i_tree, &
                            GP_Individual_Node_Type(i_Node_Right,i_Tree)


                        call print_entire_tree( )

                    endif ! myid == 0

                    i_Error = 1

                endif !   GP_Individual_Node_Type(i_Node_Right,i_Tree) .lt. -n_CODE_Equations



            endif !  GP_Individual_Node_Type(i_Function,i_Tree) .gt. 0

        enddo ! i_node

    enddo ! i_level

enddo ! i_tree


return

end subroutine GP_Check_Terminals
