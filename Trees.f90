!     
! File:   Fasham_Trees.f90
! Author: Dave Coulter
!
! Created on June 24, 2013, 11:52 AM
!

subroutine Build_Trees( treeSlice )
    
use class_Tree_Node

use TreeNodeFactory 
use GP_variables_module
use Fasham_Variables
use Fasham_Tree_Interfaces

!use tree_module
    
type(Tree_Node_Pointer), dimension(n_Trees) :: treeSlice

!type(Tree_Node),  dimension(n_Trees) :: my_trees
!type(Tree_Node),  pointer :: p_trees
    
type(Tree_Node), pointer :: root

integer :: itree
integer :: inode

!-----------------------------------------------------------------------------------------------

!  create trees from the GP_Individual_Node_Type which was read in


write(6,'(/A/)') 'Trees: create trees from GP_Individual_Node_Type  '

!do  itree = 1, 1  !n_Trees
!
!    counter = 0 
!    do  inode = 1, n_nodes
!
!        counter =  counter + 1 
!
!        number = GP_Individual_Node_Type( inode, itree ) 
!
!        p_trees => my_trees(itree)
!
!        call insert( p_trees, counter, number )                 
!
!        write(6,'(/A/)') '-------------------------------------------'
! 
!        call print_tree( p_trees ) 
!
!    enddo ! inode
!
!enddo ! itree 

! should create trees from the GP_Individual_Node_Type and GP_Individual_Node_parameter arrays

write(6,'(/A/)')      'Trees: call Deserialize_Trees2 '
write(6,'(A,1x,I6)')  'Trees: n_Tracked_resources ', n_Tracked_resources
write(6,'(A,1x,I6/)') 'Trees: n_trees ', n_trees


call Deserialize_Trees2( treeSlice, n_Tracked_resources, n_trees    ) 

write(6,'(/A/)') 'Trees: aft call Deserialize_Trees2 '


!-----------------------------------------------------------------------------------------------

!    ! Column 1
!    treeSlice(8)%n => GetNonMotileDilution(SPECIES_NITRATE)
!    treeSlice(15)%n => GetNonMotileDilution(SPECIES_AMMONIUM)
!    treeSlice(22)%n => GetNonMotileDilution(SPECIES_DISSOLVED_ORGANIC_NITROGEN)
!    treeSlice(29)%n => GetNonMotileDetritusDilution()
!    treeSlice(36)%n => GetNonMotileDilution(SPECIES_BACTERIA)
!    treeSlice(43)%n => GetNonMotileDilution(SPECIES_PHYTOPLANKTON)
!    treeSlice(50)%n => GetMotileDilution() ! Zooplankton
!    
!    ! Column 2
!    treeSlice(1)%n => GetNitrateInjection() ! Initial Nitrate - [mmol N m-3]
!    
!    ! Column 3
!    treeSlice(38)%n => Bacterial_Mortality_To_NH4()
!    treeSlice(52)%n => Zooplankton_Sink_To_NH4()
!    
!    ! Column 4
!    treeSlice(32)%n => Detrital_Sink_To_DON()
!    treeSlice(46)%n => Phytoplankton_Exudation_To_DON()
!    treeSlice(53)%n => Zooplankton_Excretion_To_DON()
!    
!    ! Column 5
!    treeSlice(47)%n => Phytoplankton_Sink_To_DET()
!    treeSlice(54)%n => Zooplankton_Sink_To_Detritus()
!    
!    ! Column 6
!    treeSlice(19)%n => NH4_Sink_To_Bacteria()
!    treeSlice(26)%n => DON_Sink_To_Bacteria()
!    
!    ! Column 7
!    treeSlice(13)%n => Nitrate_Sink_To_Phytoplankton()
!    treeSlice(20)%n => Ammonium_Sink_To_Phytoplankton()
!    
!    ! Column 8
!    treeSlice(35)%n => G3()
!    treeSlice(42)%n => G2()
!    treeSlice(49)%n => G1()

! LV model

!--------------------------------------------------
!      X        !      1         |       2        |
!               !                |                |
!--------------------------------------------------
!      3        !      X         |       4        |
!               !                |                |
!--------------------------------------------------
!      5        !      6         |       X        |
!               !                |                |
!-------------------------------------------------


!write(6,'(/A)') 'Trees: enter Build_Trees' 
!write(6,'(A,1x,I6/)') 'Trees: n_Trees ', n_Trees



    
! Column 1



!write(6,'(/A/)') 'Trees: treeSlice(5)%n => make_tree_5()'
!treeSlice(5)%n => make_tree_5()
!write(6,'(/A/)') 'Trees: AFT treeSlice(5)%n => make_tree_5()'



! Column 2

!write(6,'(/A/)') 'Trees: treeSlice(1)%n => make_tree_1()'
!treeSlice(1)%n => make_tree_1() ! Initial Nitrate - [mmol N m-3]
!write(6,'(/A/)') 'Trees: AFT treeSlice(1)%n => make_tree_1()'


! Column 3

!write(6,'(/A/)') 'Trees:  treeSlice(4)%n => make_tree_4() '
!treeSlice(4)%n => make_tree_4()
!write(6,'(/A/)') 'Trees:  AFT treeSlice(4)%n => make_tree_4() '

! debug >>>>
write(6,'(A,2(1x,I6))') 'Trees: before tree loop '
do  i = 1, n_Trees

    write(6,'(A,2(1x,I6))') 'Trees: i, treeSlice(i)%n%node_type ', &
                                    i, treeSlice(i)%n%node_type
    !write(6,'(A,2(1x,I6))') 'Trees: i, treeSlice(i)%n%node_count ', &
    !                                i, treeSlice(i)%n%node_count

enddo  ! i 
write(6,'(A,2(1x,I6))') 'Trees: after tree loop '
! debug <<<<<


end subroutine Build_Trees
