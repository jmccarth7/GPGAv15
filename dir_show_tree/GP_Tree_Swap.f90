subroutine GP_Tree_Swap
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=8) :: cff

integer(kind=4),dimension(n_nodes,2)  :: Child_Tree_Swap_Node_Type
integer(kind=4),dimension(n_nodes,2)  :: Node_Depth
integer(kind=4) :: i_Parent
integer(kind=4) :: i_Child
integer(kind=4) :: i_Parent_One
integer(kind=4) :: i_Parent_Two
!integer(kind=4) :: i_Level
integer(kind=4) :: i_Levels
integer(kind=4) :: i_Parent_Level
integer(kind=4) :: i_Child_Level
integer(kind=4) :: i_Node_at_Level
integer(kind=4) :: n_Nodes_at_Level
!integer(kind=4) :: i_Node
integer(kind=4) :: j_Node
integer(kind=4) :: k_Node
integer(kind=4) :: i
integer(kind=4) :: icnt
integer(kind=4) :: icff
integer(kind=4) :: icnt_parent_one_nodes
integer(kind=4) :: icnt_parent_two_nodes
integer(kind=4) :: i_parent_one_swap_node
integer(kind=4) :: i_parent_two_swap_node
integer(kind=4) :: i_parent_swap_node
integer(kind=4) :: i_child_swap_node
integer(kind=4) :: parent_max_swap_level
integer(kind=4) :: child_max_swap_level
integer(kind=4) :: i_parent_node_pt
integer(kind=4) :: i_child_node_pt
integer(kind=4) :: i_level_node
integer(kind=4) :: n_parent_one_swap_levels
integer(kind=4) :: n_parent_two_swap_levels
integer(kind=4) :: parent_one_max_swap_level
integer(kind=4) :: parent_two_max_swap_level
integer(kind=4),dimension(n_nodes) :: parent_two_swappable_nodes

logical node_not_found

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!test the code's ability to carry out a tree code swap
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Child_Tree_Swap_Node_Type=-9999

i_parent_one=1
i_parent_two=2

! determine from the bottom of the tree what the summed tree levels are for each node for both parents
do i=1,2
  if( i .eq. 1) i_parent=i_parent_one
  if( i .eq. 2) i_parent=i_parent_two
  do i_level=n_levels,1,-1
    n_nodes_at_level=int(2**(i_level-1))
    do i_node_at_level=1,n_nodes_at_level
      i_node=(n_nodes_at_level-1)+i_node_at_level
      node_depth(i_node,i_parent)=0
      if( parent_Tree_Swap_Node_Type(i_node,i_parent) .ge. -1) then
        if( i_level .eq. n_levels) then
          node_depth(i_node,i_parent)=1
        else
          j_node=i_node*2
          k_node=1+(i_node*2)
          if( node_depth(j_node,i_parent) .ge. 1 .and. node_depth(k_node,i_parent) .ge. 1) then
            if( node_depth(k_node,i_parent) .gt. node_depth(j_node,i_parent)) then
              j_node=k_node
            endif
            node_depth(i_node,i_parent)=node_depth(j_node,i_parent)+1
          else
            node_depth(i_node,i_parent)=1
          endif
        endif
      endif
    enddo
  enddo
enddo

! count the number of nodes on parent_one with code on it.
icnt_parent_one_nodes=0
do i_node=1,n_nodes
  if( Parent_Tree_Swap_Node_Type(i_node,i_parent_one) .ge. -1) then
    icnt_parent_one_nodes=icnt_parent_one_nodes+1
  endif
enddo

! randomly choose parent_one's node swap location
call random_number(cff) ! uniform random number generator
i_parent_one_swap_node=1+int(float(icnt_parent_one_nodes)*cff)
icnt=0
node_not_found=.true.
do i_node=1,n_nodes
  if( node_not_found) then
    if( parent_Tree_Swap_Node_Type(i_node,i_parent_one) .ge. -1) then
      icnt=icnt+1
      if( icnt .eq. i_parent_one_swap_node) then
        i_parent_one_swap_node=i_node
        node_not_found=.false.
      endif
    endif
  endif
enddo
!
n_parent_one_swap_levels=node_depth(i_parent_one_swap_node,i_parent_one)

! find parent_one's swap level
parent_one_max_swap_level=1
do i_level=2,n_levels
  icff=(2**(i_level-1))-1
  if( i_parent_one_swap_node .gt. icff) then
    parent_one_max_swap_level=i_level
  endif
enddo
parent_one_max_swap_level=n_levels-parent_one_max_swap_level+1

!.....................................................................
!     determine the range of levels that the swap can occur over
!.....................................................................
! NOTE: This is important because this code is limited in the number
!       of levels that can be created.  This was set up this way to eliminate the
!       likelihood of having the 'bloat' problem crop up and also to allow
!       for the GPCODE algorithm to be ported over to F90
!.....................................................................

! select the parent_two nodes that can be swapped with
i_node=0
icnt_parent_two_nodes=0
do i_level=1,n_levels-1
  n_nodes_at_level=int(2**(i_level-1))

  parent_two_max_swap_level=n_levels-i_level+1

  do i_level_node=1,n_nodes_at_level
    i_node=i_node+1
    if( parent_Tree_Swap_Node_Type(i_node,i_parent_two) .ge. -1 .and. &
        node_depth(i_node,i_parent_two) .le. parent_one_max_swap_level .and. &
        n_parent_one_swap_levels .le. parent_two_max_swap_level) then
      icnt_parent_two_nodes=icnt_parent_two_nodes+1
      parent_two_swappable_nodes(icnt_parent_two_nodes)=i_node
    endif
  enddo
enddo

! randomly choose parent two's node swap location
call random_number(cff) ! uniform random number generator
icff=1+int(float(icnt_parent_two_nodes)*cff)
i_parent_two_swap_node=parent_two_swappable_nodes(icff)

n_parent_two_swap_levels=node_depth(i_parent_two_swap_node,i_parent_two)

! find parent_two's swap level
parent_two_max_swap_level=1
do i_level=2,n_levels
  icff=(2**(i_level-1))-1
  if( i_parent_two_swap_node .gt. icff) then
    parent_two_max_swap_level=i_level
  endif
enddo
parent_two_max_swap_level=n_levels-parent_two_max_swap_level+1

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! do the binary tree swap to create two new child tree structures
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

do i_child=1,2

  if( i_child .eq. 1) then
    i_parent_swap_node=i_parent_two_swap_node
    i_child_swap_node=i_parent_one_swap_node
    parent_max_swap_level=parent_two_max_swap_level
    child_max_swap_level=parent_one_max_swap_level
  else
    i_parent_swap_node=i_parent_one_swap_node
    i_child_swap_node=i_parent_two_swap_node
    parent_max_swap_level=parent_one_max_swap_level
    child_max_swap_level=parent_two_max_swap_level
  endif

  if( i_child .eq. 1) then
    i_parent=i_parent_one
  else
    i_parent=i_parent_two
  endif
  child_Tree_Swap_Node_Type(1:n_Nodes,i_child)=parent_Tree_Swap_Node_Type(1:n_Nodes,i_parent)

! set each of the children to one of the parents initially
  if( i_child .eq. 1) then
    i_parent=i_parent_two
  else
    i_parent=i_parent_one
  endif
! swap the node from the start of the tree
  child_Tree_Swap_Node_Type(i_child_swap_node,i_child)=parent_Tree_Swap_Node_Type(i_parent_swap_node,i_parent)

  i_levels=0
! clean out the bottom of the tree branches on the child
  do i_child_level=n_levels-child_max_swap_level+2,n_levels
    i_levels=i_levels+1
    i_child_swap_node=i_child_swap_node*2
    do i_node=1,2**i_levels
      i_child_node_pt=i_child_swap_node-1+i_node
      child_Tree_Swap_Node_Type(i_child_node_pt,i_child)=-9999
    enddo
  enddo

  if( i_child .eq. 1) then
    i_parent_swap_node=i_parent_two_swap_node
    i_child_swap_node=i_parent_one_swap_node
    parent_max_swap_level=parent_two_max_swap_level
    child_max_swap_level=parent_one_max_swap_level
  else
    i_parent_swap_node=i_parent_one_swap_node
    i_child_swap_node=i_parent_two_swap_node
    parent_max_swap_level=parent_one_max_swap_level
    child_max_swap_level=parent_two_max_swap_level
  endif

  i_Levels=0
! add in the branches
  do i_Child_Level=n_Levels-Child_Max_Swap_Level+2,n_Levels
    i_Levels=i_Levels+1
    i_Parent_Level=n_Levels-Parent_Max_Swap_Level+1+i_Levels
    if( i_parent_level .le. n_levels) then
      i_child_swap_node=i_child_swap_node*2
      i_parent_swap_node=i_parent_swap_node*2
      do i_node=1,2**i_levels
        i_child_node_pt=i_child_swap_node-1+i_node
        i_parent_node_pt=i_parent_swap_node-1+i_node
        Child_Tree_Swap_Node_Type(i_Child_Node_Pt,i_Child)=Parent_Tree_Swap_Node_Type(i_Parent_Node_Pt,i_Parent)
      enddo
    endif

  enddo

enddo

Parent_Tree_Swap_Node_Type=Child_Tree_Swap_Node_Type

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Tree_Swap
!23456789012345678901234567890123456789012345678901234567890123456789012
