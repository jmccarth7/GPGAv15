subroutine Build_Trees( treeSlice, buildtrees ) 

use mpi
use mpi_module

use class_Tree_Node

use Tree_Node_Factory_module
use GP_variables_module
use fasham_variables_module
use Fasham_Tree_Interfaces

!!!!!use fasham_tree_functions_module

type(Tree_Node_Pointer), dimension(n_Trees) :: treeSlice
    
logical, intent(in) :: buildtrees

!integer :: itree
!integer :: inode


!------------------------------------------------------------------------------------------------

write(6,'(A,1x,A /)') 'build_trees: model   ', trim(model) 

if( buildtrees )then

    !  create trees from the GP_Individual_Node_Type which was read in
    
    !if( myid == 1 )then
        write(6,'(/A/)')      'build_trees: create trees from GP_Individual_Node_Type  '
        write(6,'(/A/)')      'build_trees: call Deserialize_Trees2 '
        write(6,'(A,1x,I6)')  'build_trees: n_Tracked_resources ', n_Tracked_resources
        write(6,'(A,1x,I6/)') 'build_trees: n_trees ', n_trees
    !endif ! myid == 1
    
    
    ! Deserialize_Trees2 should create trees from
    ! the GP_Individual_Node_Type and GP_Individual_Node_parameter arrays
    
    call Deserialize_Trees2( treeSlice, n_Tracked_resources, n_trees    )
    
    
    !if( myid == 1 )then
        write(6,'(/A/)') 'build_trees: aft call Deserialize_Trees2 '
    !endif ! myid == 1
    
    
    
    ! debug >>>>
    
    !if( myid == 1 )then
    !
    !    write(6,'(A,2(1x,I6))') 'build_trees: before tree loop '
    !
    !    do  i = 1, n_Trees
    !        write(6,'(A,2(1x,I6))') 'build_trees: i, treeSlice(i)%n%node_type ', &
    !                                              i, treeSlice(i)%n%node_type
    !    enddo  ! i
    !
    !    write(6,'(A,2(1x,I6))') 'build_trees: after tree loop '
    !
    !endif ! myid == 1
    
    ! debug <<<<<

else 

    if( trim(model) == 'fasham' .or. &
        trim(model) == 'FASHAM'          )then
    
    
        !  Fasham specific trees
    
        !-----------------------------------------
        ! File:   Fasham_Trees.f90
        ! Author: Dave Coulter
        ! Created on June 24, 2013, 11:52 AM
        !-----------------------------------------
        
        ! Column 1
    
        write(6,'(A/)') 'build_trees:  set Fasham tree pointers '
    
        treeSlice( 8)%n => GetNonMotileDilution(SPECIES_NITRATE)
        treeSlice(15)%n => GetNonMotileDilution(SPECIES_AMMONIUM)
        treeSlice(22)%n => GetNonMotileDilution(SPECIES_DISSOLVED_ORGANIC_NITROGEN)
        treeSlice(29)%n => GetNonMotileDetritusDilution()
        treeSlice(36)%n => GetNonMotileDilution(SPECIES_BACTERIA)
        treeSlice(43)%n => GetNonMotileDilution(SPECIES_PHYTOPLANKTON)
        treeSlice(50)%n => GetMotileDilution() ! Zooplankton
        
        ! Column 2
        treeSlice( 1)%n => GetNitrateInjection() ! Initial Nitrate - [mmol N m-3]
        
        ! Column 3
        treeSlice(38)%n => Bacterial_Mortality_To_NH4()
        treeSlice(52)%n => Zooplankton_Sink_To_NH4()
        
        ! Column 4
        treeSlice(32)%n => Detrital_Sink_To_DON()
        treeSlice(46)%n => Phytoplankton_Exudation_To_DON()
        treeSlice(53)%n => Zooplankton_Excretion_To_DON()
        
        ! Column 5
        treeSlice(47)%n => Phytoplankton_Sink_To_DET()
        treeSlice(54)%n => Zooplankton_Sink_To_Detritus()
        
        ! Column 6
        treeSlice(19)%n => NH4_Sink_To_Bacteria()
        treeSlice(26)%n => DON_Sink_To_Bacteria()
        
        ! Column 7
        treeSlice(13)%n => Nitrate_Sink_To_Phytoplankton()
        treeSlice(20)%n => Ammonium_Sink_To_Phytoplankton()
        
        ! Column 8
        treeSlice(35)%n => f_G3()
        treeSlice(42)%n => f_G2()
        treeSlice(49)%n => f_G1()
        
    endif ! model == 'fasham'  

endif !  buildtrees 



return

end subroutine Build_Trees
