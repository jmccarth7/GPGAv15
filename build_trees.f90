subroutine Build_Trees( treeSlice ) 

use mpi
use mpi_module

use class_Tree_Node

use Tree_Node_Factory_module
use GP_variables_module
use fasham_variables_module
use Fasham_Tree_Interfaces


type(Tree_Node_Pointer), dimension(n_Trees) :: treeSlice

!type(Tree_Node), pointer :: root

!integer :: itree
!integer :: inode


!-----------------------------------------------------------------------------------------------

!  create trees from the GP_Individual_Node_Type which was read in

!if( myid == 1 )then
!    write(6,'(/A/)')      'build_trees: create trees from GP_Individual_Node_Type  '
!    write(6,'(/A/)')      'build_trees: call Deserialize_Trees2 '
!    write(6,'(A,1x,I6)')  'build_trees: n_Tracked_resources ', n_Tracked_resources
!    write(6,'(A,1x,I6/)') 'build_trees: n_trees ', n_trees
!endif ! myid == 1


! Deserialize_Trees2 should create trees from
! the GP_Individual_Node_Type and GP_Individual_Node_parameter arrays

call Deserialize_Trees2( treeSlice, n_Tracked_resources, n_trees    )


!if( myid == 1 )then
!    write(6,'(/A/)') 'build_trees: aft call Deserialize_Trees2 '
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


end subroutine Build_Trees
