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


integer(kind=4), intent(in)  :: n_indiv
integer(kind=4), intent(in), dimension( n_indiv,n_Nodes,n_Trees) :: indiv_node_type

integer(kind=4) :: i_indiv
integer(kind=4) :: i_diversity
integer(kind=4) :: icnt_Nodes
integer(kind=4) :: icnt_parms
integer(kind=4) :: icnt_vars
integer(kind=4) :: icnt_ops
integer(kind=4) :: max_number_nodes
real(kind=8)    :: xmax_number_nodes

integer(kind=4) :: i_GP_individual                                                                           
integer(kind=4) :: i_GP_Best_Parent                                                                          
integer(kind=4) :: i_GP_Generation                                                                           
integer(kind=4) :: j_GP_Individual                                                                           
integer(kind=4) :: i_Tree                                                                                    
integer(kind=4) :: i_Node            

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

max_number_nodes = n_trees * n_nodes
xmax_number_nodes = real( max_number_nodes, kind=8 )


write(GP_print_unit,'(/A,1x,I6)') 'gcdi: i_GP_generation = ', i_GP_generation
write(GP_print_unit,'(A)') &
      'gcdi: i_indiv i_diversity  icnt_Nodes  icnt_parms  icnt_vars  icnt_ops'

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

    icnt_ops   = int( 100.0d0 * real( icnt_ops,   kind=8 ) / xmax_number_nodes )
    icnt_vars  = int( 100.0d0 * real( icnt_vars,  kind=8 ) / xmax_number_nodes )
    icnt_parms = int( 100.0d0 * real( icnt_parms, kind=8 ) / xmax_number_nodes )

    i_diversity = 100 * ( icnt_parms * 100  + icnt_vars ) + icnt_ops

    GP_diversity_index( i_indiv ) = i_diversity

    write(GP_print_unit,'(6(1x,I10))') &
          i_indiv, i_diversity, icnt_Nodes, icnt_parms, icnt_vars, icnt_ops

enddo  ! i_indiv

write(GP_print_unit,'(/A)') 'gcdi:GP_diversity_index '
write(GP_print_unit,'(10(1x,I6))') GP_diversity_index(1:n_GP_individuals)


return

end subroutine GP_calc_diversity_index

