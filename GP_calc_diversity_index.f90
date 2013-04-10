!subroutine GP_calc_diversity_index( i_indiv, n_indiv, indiv_node_type, i_diversity )
subroutine GP_calc_diversity_index( n_indiv, indiv_node_type, i_diversity )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! Mutations in this subroutine are targeted to the nodes only.  
! The terminals are optimized later on using GA_lmdif.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=8) :: cff

!integer(kind=4), intent(in)  :: i_indiv       
integer(kind=4), intent(in)  :: n_indiv       
integer(kind=4), intent(in), dimension( n_indiv,n_Nodes,n_Trees) :: indiv_node_type

integer(kind=4) :: i_indiv       
integer(kind=4) :: i_diversity   
integer(kind=4) :: icnt_Nodes
integer(kind=4) :: icnt_parms
integer(kind=4) :: icnt_vars 
integer(kind=4) :: icnt_ops  
integer(kind=4) :: icnt


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!write(6,'(A)') ' '
write(6,'(/A,1x,I6)') 'gcdi: i_GP_generation = ', i_GP_generation
write(6,'(A)') &
      'gcdi: i_indiv, i_diversity, icnt_Nodes, icnt_parms, icnt_vars, icnt_ops'

do  i_indiv = 1, n_GP_individuals


    i_diversity   = 0
        
    icnt_Nodes=0
    icnt_parms=0
    icnt_vars =0
    icnt_ops  =0
    
    ! count the number of nodes on the tree selected for a mutation.  
    ! Only nodes are mutated.
    
    do  i_tree=1,n_trees
        do  i_Node=1,n_Nodes
    
            if( indiv_node_type(i_indiv, i_Node, i_tree )  > -9999 ) then

                icnt_Nodes=icnt_Nodes+1
    
                if( indiv_node_type(i_indiv, i_Node, i_tree ) < 0 )then
                    icnt_parms = icnt_parms + 1
                endif ! indiv_node_type... < 0
    
                if( indiv_node_type(i_indiv, i_Node, i_tree ) == 0 )then
                    icnt_vars  = icnt_vars  + 1
                endif ! indiv_node_type... == 0
    
                if( indiv_node_type(i_indiv, i_Node, i_tree )  > 0 )then
                    icnt_ops   = icnt_ops   + 1
                endif ! indiv_node_type... > 0
    
            endif ! indiv_node_type...
    
        enddo ! i_node
    
    enddo !  i_tree
    
    
    i_diversity = icnt_Nodes + icnt_parms + icnt_vars + icnt_ops

    GP_diversity_index( i_indiv ) = i_diversity
    
    !write(6,'(A,6(1x,I6))') &
    !      'gcdi: i_indiv, i_diversity, icnt_Nodes, icnt_parms, icnt_vars, icnt_ops', &
    !             i_indiv, i_diversity, icnt_Nodes, icnt_parms, icnt_vars, icnt_ops
    write(6,'(6(1x,I6))') &
          i_indiv, i_diversity, icnt_Nodes, icnt_parms, icnt_vars, icnt_ops

enddo  ! i_indiv 

write(6,'(A)') ' '
write(6,'(A,1x,50(1x,I2))') &
         'gcdi:GP_diversity_index ',  GP_diversity_index(1:n_GP_individuals) 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_calc_diversity_index
!1234567890123456789012345678901234567890123456789012345678901234567890
