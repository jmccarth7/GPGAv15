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
integer(kind=4) :: i_Node_Count
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
integer(kind=4) :: i_parent_node_Point
integer(kind=4) :: i_child_node_Point
integer(kind=4) :: i_level_node
integer(kind=4) :: n_parent_one_swap_levels
integer(kind=4) :: n_parent_two_swap_levels
integer(kind=4) :: parent_one_max_swap_level
integer(kind=4) :: parent_two_max_swap_level
integer(kind=4),dimension(n_nodes) :: parent_two_swappable_nodes

logical :: MALE_CROSS
logical :: FEMALE_CROSS
logical :: NODE_NOT_FOUND

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!test the code's ability to carry out a tree code swap
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! Look to see if there is code on both Male and Female trees to CROSS

!--------------------------------------------------------------------------

MALE_CROSS=.true.

i_Node_Count=0

do  i_Node=1,n_Nodes
    if( Parent_Tree_Swap_Node_Type(i_Node,1) .ne. -9999) i_Node_Count=i_Node_Count+1
enddo

if( i_Node_Count .eq. 0) then
    MALE_CROSS=.false. ! there are no Trees on this GP_CODE
endif

!--------------------------------------------------------------------------

FEMALE_CROSS=.true.

i_Node_Count=0
do  i_Node=1,n_Nodes
    if( Parent_Tree_Swap_Node_Type(i_Node,2) .ne. -9999) i_Node_Count=i_Node_Count+1
enddo

if( i_Node_Count .eq. 0) then
    FEMALE_CROSS=.false. ! there are no Trees on this GP_CODE
endif

!--------------------------------------------------------------------------

Child_Tree_Swap_Node_Type=-9999

i_parent_one=1
i_parent_two=2

if( MALE_CROSS .and. FEMALE_CROSS) then

    ! determine from the bottom of the tree 
    ! what the summed tree levels are for each node for both parents

    do  i=1,2

        if( i .eq. 1) i_parent=i_parent_one
        if( i .eq. 2) i_parent=i_parent_two

        do  i_level=n_levels,1,-1
            n_nodes_at_level=int(2**(i_level-1))

            do  i_node_at_level=1,n_nodes_at_level
                i_node=(n_nodes_at_level-1)+i_node_at_level
                Node_depth(i_node,i_parent)=0

                !!if( parent_Tree_Swap_Node_Type(i_node,i_parent) .ge. -1) then
                if( Parent_Tree_Swap_Node_Type(i_Node,i_parent) .ne. -9999) then

                    if( i_level .eq. n_levels) then

                        Node_depth(i_node,i_parent)=1

                    else

                        j_node=i_node*2
                        k_node=1+(i_node*2)

                        if( node_depth(j_node,i_parent) .ge. 1 .and. &
                            Node_depth(k_node,i_parent) .ge. 1           ) then

                            if( node_depth(k_node,i_parent) .gt. &
                                Node_depth(j_node,i_parent)         ) then
                                j_node=k_node
                            endif !   node_depth(k_node,i_parent) .gt. node_depth...                          

                            Node_depth(i_node,i_parent)=node_depth(j_node,i_parent)+1

                        else

                            Node_depth(i_node,i_parent)=1

                        endif !   node_depth(j_node,i_parent) .ge. 1 .and. ...

                    endif ! i_level .eq. n_levels

                endif !   parent_Tree_Swap_Node_Type(i_node,i_parent) .ge. -1

            enddo ! i_node_at_level
        enddo ! i_level
    enddo ! i 

    ! count the number of nodes on parent_one with code on it.
    icnt_parent_one_nodes=0
    do  i_node=1,n_nodes
        !!if( Parent_Tree_Swap_Node_Type(i_node,i_parent_one) .ge. -1) then
        if( Parent_Tree_Swap_Node_Type(i_Node,i_Parent_One) .ne. -9999) then
            icnt_parent_one_nodes=icnt_parent_one_nodes+1
        endif !   Parent_Tree_Swap_Node_Type(i_Node,i_Parent_One) .ne. -9999
    enddo ! i_node
    
    ! randomly choose parent_one's node swap location
    call random_number(cff) ! uniform random number generator
    
    i_parent_one_swap_node=1+int(float(icnt_parent_one_nodes)*cff)

    if( i_Parent_One_Swap_Node .gt. icnt_parent_one_nodes) i_Parent_One_Swap_Node=icnt_Parent_One_Nodes

    icnt=0
    node_not_found=.true.
    do  i_node=1,n_nodes
        !if( node_not_found) then
        !!!if( parent_Tree_Swap_Node_Type(i_node,i_parent_one) .ge. -1) then
        if( parent_Tree_Swap_Node_Type(i_node,i_parent_one) .ne. -9999) then
            icnt=icnt+1
            if( icnt .eq. i_parent_one_swap_node) then
                i_parent_one_swap_node=i_node
                exit
                !node_not_found=.false.
            endif !   icnt .eq. i_parent_one_swap_node
        endif !   parent_Tree_Swap_Node_Type(i_node,i_parent_one) .ne. -9999
        !endif !   node_not_found
    enddo ! i_node
    
    n_parent_one_swap_levels=node_depth(i_parent_one_swap_node,i_parent_one)
    
    ! find parent_one's swap level
    parent_one_max_swap_level=1
    
    do  i_level=2,n_levels
        icff=(2**(i_level-1))-1
        if( i_parent_one_swap_node .gt. icff) then
            parent_one_max_swap_level=i_level
        endif !   i_parent_one_swap_node .gt. icff
    enddo ! i_level
    
    parent_one_max_swap_level=n_levels-parent_one_max_swap_level+1
    
    ! write(*,*) 'i_Parent_One_Swap_Node = ',i_Parent_One_Swap_Node
    ! write(*,*) 'Parent_One_Max_Swap_Level = ',Parent_One_Max_Swap_Level
    
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

    !!!do i_level=1,n_levels-1
    ! ??? !off do i_level=1,n_levels-1
    
    do  i_Level=1,n_Levels
    
        n_nodes_at_level=int(2**(i_level-1))
    
        Parent_two_max_swap_level=n_levels-i_level+1
    
        do  i_level_node=1,n_nodes_at_level
    
            i_Node=i_Node+1
    
            !!if( parent_Tree_Swap_Node_Type(i_node,i_parent_two) .ge. -1        .and. &

            if( Parent_Tree_Swap_Node_Type(i_Node,i_Parent_two) .ne. -9999     .and. &
                Node_depth(i_node,i_parent_two) .le. parent_one_max_swap_level .and. &
                n_parent_one_swap_levels        .le. parent_two_max_swap_level         ) then
    
                icnt_parent_two_nodes=icnt_parent_two_nodes+1
                Parent_two_swappable_nodes(icnt_parent_two_nodes)=i_node
    
            endif !  parent_Tree_Swap_Node_Type(i_node,i_parent_two) .ne. -9999...
    
        enddo ! i_level_node
    enddo ! i_level
    
    ! randomly choose parent two's node swap location
    call random_number(cff) ! uniform random number generator
    
    icff=1+int(float(icnt_parent_two_nodes)*cff)

    if( icff .gt. icnt_parent_two_nodes) icff=icnt_parent_two_nodes
    
    i_parent_two_swap_node=parent_two_swappable_nodes(icff)
    
    n_parent_two_swap_levels=node_depth(i_parent_two_swap_node,i_parent_two)
    
    ! find parent_two's swap level
    parent_two_max_swap_level=1
    
    do  i_level=2,n_levels
        icff=(2**(i_level-1))-1
        if( i_parent_two_swap_node .gt. icff) then
            parent_two_max_swap_level=i_level
        endif !   i_parent_two_swap_node .gt. icff
    enddo ! i_level
    
    Parent_two_max_swap_level=n_levels-parent_two_max_swap_level+1
    
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! do the binary tree swap to create two new child tree structures
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    do  i_Child=1,1  ! we only need to keep the cross from the first child;  2
    
        if( i_child .eq. 1) then
    
            i_parent_swap_node    = i_parent_two_swap_node
            i_child_swap_node     = i_parent_one_swap_node
            Parent_max_swap_level = Parent_two_max_swap_level
            Child_max_swap_level  = Parent_one_max_swap_level
    
        else
    
            i_parent_swap_node    = i_parent_one_swap_node
            i_child_swap_node     = i_parent_two_swap_node
            Parent_max_swap_level = Parent_one_max_swap_level
            Child_max_swap_level  = parent_two_max_swap_level
    
        endif !   i_child .eq. 1
    
        if( i_child .eq. 1) then
            i_parent=i_parent_one
        else
            i_parent=i_parent_two
        endif !   i_child .eq. 1
    
        Child_Tree_Swap_Node_Type(1:n_Nodes,i_child)=parent_Tree_Swap_Node_Type(1:n_Nodes,i_parent)
    
        ! set each of the children to one of the parents initially
    
        if( i_child .eq. 1) then
            i_parent=i_parent_two
        else
            i_parent=i_parent_one
        endif ! i_child .eq. 1
    
    
        ! swap the node from the start of the tree
    
        Child_Tree_Swap_Node_Type(i_child_swap_node,i_child) = &
           Parent_Tree_Swap_Node_Type(i_parent_swap_node,i_parent)
    
        i_levels=0
        ! clean out the bottom of the tree branches on the child
        !  that the branch will be added to        [NOTE: IS THIS NEEDED STILL?????]
    
        do  i_child_level=n_levels-child_max_swap_level+2,n_levels
    
            i_levels=i_levels+1
            i_child_swap_node=i_child_swap_node*2
    
            do  i_node=1,2**i_levels
                i_Child_Node_Point=i_Child_Swap_Node-1+i_Node
                Child_Tree_Swap_Node_Type(i_Child_Node_Point,i_Child)=-9999
            enddo ! i_node
    
        enddo !  i_child_level
    
        if( i_child .eq. 1) then
            i_parent_swap_node    = i_parent_two_swap_node
            i_child_swap_node     = i_parent_one_swap_node
            Parent_max_swap_level = Parent_two_max_swap_level
            Child_max_swap_level  = Parent_one_max_swap_level
        else
            i_parent_swap_node    = i_parent_one_swap_node
            i_child_swap_node     = i_parent_two_swap_node
            Parent_max_swap_level = Parent_one_max_swap_level
            Child_max_swap_level  = Parent_two_max_swap_level
        endif !   i_child .eq. 1
    
        ! add in the branches 

        i_Levels=0

        do  i_Child_Level=n_Levels-Child_Max_Swap_Level+2,n_Levels
            i_Levels=i_Levels+1
            i_Parent_Level=n_Levels-Parent_Max_Swap_Level+1+i_Levels
            if( i_parent_level .le. n_levels) then
                i_child_swap_node  = i_child_swap_node*2
                i_parent_swap_node = i_parent_swap_node*2
                do  i_node=1,2**i_levels
                    i_Child_Node_Point=i_Child_Swap_Node-1+i_Node
                    i_Parent_Node_Point=i_Parent_Swap_Node-1+i_Node
                    Child_Tree_Swap_Node_Type(i_Child_Node_Point,i_Child) = &
                    Parent_Tree_Swap_Node_Type(i_Parent_Node_Point,i_Parent)
                enddo ! i_node
            endif ! i_parent_level .le. n_levels
    
        enddo  ! i_Child_Level
    
    enddo ! i_child
 
elseif( MALE_CROSS .and. .not. FEMALE_CROSS) then  ! the Male tree is empty

    Child_Tree_Swap_Node_Type(1:n_Nodes,1)=-9999     ! Empty out the Male Tree

elseif( .not. MALE_CROSS .and. FEMALE_CROSS) then  ! the Male tree is empty

    ! pick a node from the Female Parent_Tree and set it into the Male Child

    ! determine from the bottom of the tree 
    ! what the summed tree levels are for each node for the female parent

    i_Parent=i_Parent_Two
    do  i_Level=n_Levels,1,-1

        n_Nodes_at_Level=int(2**(i_Level-1))

        do  i_Node_at_Level=1,n_Nodes_at_Level

            i_Node=(n_Nodes_at_Level-1)+i_Node_at_Level
            Node_Depth(i_Node,i_Parent)=0

            if( Parent_Tree_Swap_Node_Type(i_Node,i_Parent) .ne. -9999) then
                if( i_Level .eq. n_Levels) then
                    Node_Depth(i_Node,i_Parent)=1
                else
                    j_Node=i_Node*2
                    k_Node=1+(i_Node*2)
                    if( Node_Depth(j_Node,i_Parent) .ge. 1 .and. &
                        Node_Depth(k_Node,i_Parent) .ge. 1         ) then

                        if( Node_Depth(k_Node,i_Parent) .gt. Node_Depth(j_Node,i_Parent)) then
                            j_Node=k_Node
                        endif
                        Node_Depth(i_Node,i_Parent)=Node_Depth(j_Node,i_Parent)+1
                    else
                        Node_Depth(i_Node,i_Parent)=1
                    endif
                endif
            endif
        enddo
    enddo

    i_Parent_One_Swap_Node=1

    ! select the parent_two nodes that can be swapped with
    i_Node=0
    icnt_Parent_Two_Nodes=0
    do  i_Level=1,n_Levels
        n_Nodes_at_Level=int(2**(i_Level-1))
        do  i_level_Node=1,n_Nodes_at_Level
            i_Node=i_Node+1
            if( Parent_Tree_Swap_Node_Type(i_Node,i_Parent_Two) .ne. -9999) then
                icnt_Parent_Two_Nodes=icnt_Parent_Two_Nodes+1
                Parent_Two_Swappable_Nodes(icnt_Parent_Two_Nodes)=i_Node
            endif
        enddo
    enddo

    ! randomly choose parent two's node swap location

    call random_number(cff) ! uniform random number generator

    icff=1+int(cff*float(icnt_Parent_Two_Nodes))
  
    i_Parent_Two_Swap_Node=Parent_Two_Swappable_Nodes(icff)
    
    n_Parent_Two_Swap_Levels=Node_Depth(i_Parent_Two_Swap_Node,i_Parent_Two)
  
    ! find parent_two's swap level
  
    Parent_Two_Max_Swap_Level=1
    do  i_Level=2,n_Levels
        icff=(2**(i_Level-1))-1
        if( i_Parent_Two_Swap_Node .gt. icff) Parent_Two_Max_Swap_Level=i_Level
    enddo
  
    Parent_Two_Max_Swap_Level=n_Levels-Parent_Two_Max_Swap_Level+1
    
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! do the binary tree swap to create two new child tree structures
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
    i_Child=1
  
    i_Parent_Swap_Node    = i_Parent_Two_Swap_Node
    i_Child_Swap_Node     = i_Parent_One_Swap_Node
    Parent_Max_Swap_Level = Parent_Two_Max_Swap_Level
    Child_Max_Swap_Level  = Parent_One_Max_Swap_Level
  
    ! swap the node from the start of the tree
    Child_Tree_Swap_Node_Type(i_Child_Swap_Node,i_Child) =  &
    Parent_Tree_Swap_Node_Type(i_Parent_Swap_Node,i_Parent_Two)
  
    i_Levels=0
  
    ! add in the branches 
    do  i_Child_Level=1,n_Levels
        i_Levels=i_Levels+1
        i_Parent_Level=n_Levels-Parent_Max_Swap_Level+1+i_Levels

        if( i_Parent_Level .le. n_Levels) then

            i_Child_Swap_Node=i_Child_Swap_Node*2
            i_Parent_Swap_Node=i_Parent_Swap_Node*2
            do  i_Node=1,2**i_Levels
                i_Child_Node_Point=i_Child_Swap_Node-1+i_Node
                i_Parent_Node_Point=i_Parent_Swap_Node-1+i_Node
                Child_Tree_Swap_Node_Type(i_Child_Node_Point,i_Child) =  &
                Parent_Tree_Swap_Node_Type(i_Parent_Node_Point,i_Parent)
            enddo

        endif

    enddo

elseif( .not. MALE_CROSS .and. .not. FEMALE_CROSS) then  ! both Male and Female trees are empty

    Child_Tree_Swap_Node_Type(1:n_Nodes,1)=-9999    ! Empty out the Male Tree

endif !   MALE_CROSS .and. FEMALE_CROSS

Parent_Tree_Swap_Node_Type=Child_Tree_Swap_Node_Type

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Tree_Swap

