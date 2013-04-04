subroutine GP_Clean_Tree_Nodes
!
! This subroutine looks through the GP_Adult_Population_Node_Type array
! for nodes that have both terminals set to a parameter.  The routine replaces
! these nodes with a parameter setting and re-sets the terminals to that node as -9999
! This helps to maintain simplicity within the tree structures.
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none
!integer (kind=4) :: i_GP_Individual
!integer (kind=4) :: i_Tree
!integer (kind=4) :: i_Node
!integer (kind=4) :: i_Level
!integer (kind=4) :: i_Function
!integer (kind=4) :: i_Node_Left
!integer (kind=4) :: i_Node_Right

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

do i_GP_Individual=1,n_GP_Individuals

  do i_Tree=1,n_Trees

    do i_Level=n_Levels-1,1,-1                      ! move up the tree structure from level "n_level-1" to level "1"

      i_Function=(2**(i_Level-1))-1                 ! calculated the function number at the right end of the upper level

      do i_Node=2**i_Level,(2*(2**i_Level))-1,2     ! run through each function at the level

        i_Function=i_Function+1                     ! sets the 'function' node's index
        i_Node_Left=i_Node                          ! sets the 'left terminal' node's index; i_node_left=i_function*2 would also work
        i_Node_Right=i_Node+1                       ! sets the 'right terminal' node's index; i_node_right=(i_function*2)+1 would also work
          
        if( GP_Adult_Population_Node_Type(i_GP_Individual,i_Function,i_Tree)   .gt. 0 .and. &
            GP_Adult_Population_Node_Type(i_GP_Individual,i_Node_Left,i_Tree)  .eq. 0 .and. &
            GP_Adult_Population_Node_Type(i_GP_Individual,i_Node_Right,i_Tree) .eq. 0           ) then

            GP_Adult_Population_Node_Type(i_GP_Individual,i_Function,i_Tree)=0
            GP_Adult_Population_Node_Type(i_GP_Individual,i_Node_Left,i_Tree)=-9999
            GP_Adult_Population_Node_Type(i_GP_Individual,i_Node_Right,i_Tree)=-9999

        endif ! GP_Adult_Population_Node_Type(i_GP_Individual,i_Function,i_Tree) .gt. 0 ...

      enddo ! i_node
    enddo ! i_level
  enddo ! i_tree

enddo !  i_GP_Individual


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Clean_Tree_Nodes
!23456789012345678901234567890123456789012345678901234567890123456789012
