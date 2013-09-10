subroutine GP_Clean_Tree_Nodes
!
! This subroutine looks through the GP_Adult_Population_Node_Type array
! for nodes that have both terminals set to a parameter.  The routine replaces
! these nodes with a parameter setting and re-sets the terminals to that node as -9999
! This helps to maintain simplicity within the tree structures.


use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

integer(kind=4) :: i_GP_individual
integer(kind=4) :: i_Tree
integer(kind=4) :: i_level
integer(kind=4) :: i_function
integer(kind=4) :: i_Node
integer(kind=4) :: i_Node_Left
integer(kind=4) :: i_Node_Right

!---------------------------------------------------------------------------

do  i_GP_Individual=1,n_GP_Individuals

    !write(6,'(/A,1x,I6)') 'gpctn: i_GP_Individual ', i_GP_Individual 

    do  i_Tree=1,n_Trees
        !write(6,'(A,1x,I6)') 'gpctn: i_tree ', i_tree

        ! move up the tree structure from level "n_level-1" to level "1"

        do  i_Level = n_Levels-1, 1, -1

            !write(6,'(A,1x,I6)') 'gpctn: i_level ', i_level

            ! calculated the function number at the right end of the upper level

            i_Function = 2**(i_Level-1) - 1

            !write(6,'(A,1x,I6)') 'gpctn: i_function ', i_function

            ! run through each function at the level

            do  i_Node = 2**i_Level, 2*(2**i_Level)-1, 2

                i_Function=i_Function+1  ! sets the 'function' node's index

                !write(6,'(A,1x,I6)') 'gpctn: i_node     ', i_node     
                !write(6,'(A,1x,I6)') 'gpctn: i_function ', i_function

                i_Node_Left=i_Node       ! sets the 'left terminal' node's index;
                                         ! i_node_left=i_function*2 would also work
                i_Node_Right=i_Node+1    ! sets the 'right terminal' node's index;
                                         ! i_node_right=(i_function*2)+1 would also work

                !write(6,'(A,2(1x,I6))') &
                !      'gpctn:  i_Function, GP_Adult_Population_Node_Type(i_Function, i_Tree,i_GP_Individual)', &
                !               i_Function, GP_Adult_Population_Node_Type(i_Function, i_Tree,i_GP_Individual)
                !write(6,'(A,2(1x,I6))') &
                ! 'gpctn: i_Node_Left, GP_Adult_Population_Node_Type(i_Node_Left, i_Tree,i_GP_Individual)', &
                !         i_Node_Left, GP_Adult_Population_Node_Type(i_Node_Left, i_Tree,i_GP_Individual)
                !write(6,'(A,2(1x,I6))') &
                ! 'gpctn: i_Node_Right, GP_Adult_Population_Node_Type(i_Node_Right, i_Tree,i_GP_Individual)', &
                !         i_Node_Right, GP_Adult_Population_Node_Type(i_Node_Right, i_Tree,i_GP_Individual)


                if( GP_Adult_Population_Node_Type(i_Function,  i_Tree,i_GP_Individual) .gt. 0 .and. &
                    GP_Adult_Population_Node_Type(i_Node_Left, i_Tree,i_GP_Individual) .eq. 0 .and. &
                    GP_Adult_Population_Node_Type(i_Node_Right,i_Tree,i_GP_Individual) .eq. 0 ) then

                    GP_Adult_Population_Node_Type(i_Function,  i_Tree,i_GP_Individual) = 0
                    GP_Adult_Population_Node_Type(i_Node_Left, i_Tree,i_GP_Individual) = -9999
                    GP_Adult_Population_Node_Type(i_Node_Right,i_Tree,i_GP_Individual) = -9999

                endif ! GP_Adult_Population_Node_Type(i_Function,i_Tree,i_GP_Individual) .gt. 0 ...

            enddo ! i_node

        enddo ! i_level

    enddo ! i_tree

enddo !  i_GP_Individual


return

end subroutine GP_Clean_Tree_Nodes
