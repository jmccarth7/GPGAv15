subroutine summary_GP_indiv( i_GP_generation, i_GP_indiv )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Adult_Population_Node_Type(    i_GP_indiv,1:n_Nodes,1:n_Trees  )

!  GP_Population_Initial_Conditions( i_GP_indiv,1:n_CODE_Equations   )

!  GP_Population_Node_Parameters(    i_GP_indiv,1:n_Nodes,1:n_Trees  )


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module


implicit none



integer :: i_code_eq


integer :: i_GP_indiv
integer(kind=4),intent(in)  :: i_GP_Generation                                                                           
integer(kind=4) :: i_Tree                                                                                    
integer(kind=4) :: i_Node            

!----------------------------------------------------------------------------------------

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Adult_Population_Node_Type(    i_GP_indiv,1:n_Nodes,1:n_Trees)
!  GP_Population_Initial_Conditions( i_GP_indiv,1:n_CODE_Equations)
!  GP_Population_Node_Parameters(    i_GP_indiv,1:n_Nodes,1:n_Trees)

!---------------------------------------------------
! assume this subroutine is called only by cpu 0
!---------------------------------------------------

!--------------------------------------------------------------------------------

! initial conditions


write(GP_print_unit,'(/A)')&
      'sgpi: i_GP_gen, i_GP_indiv, i_code_eq, &
            &GP_Pop_Init_Cond(i_GP_Indiv,i_code_eq) '

do  i_code_eq = 1, n_CODE_Equations

    write(GP_print_unit,'(3(1x,I10), 7x, E15.7)')&
          i_GP_generation, i_GP_indiv, i_code_eq, &
          GP_Population_Initial_Conditions( i_GP_indiv,i_code_eq )

    write(GP_summary_output_unit, '(3(1x,I10), 1x, E20.10,2x,A)')&
          i_GP_generation, i_GP_indiv, i_code_eq, &
          GP_Population_Initial_Conditions( i_GP_indiv,i_code_eq ), &
          'gen_indiv_eq'

enddo  ! i_code_eq


write(GP_summary_output_unit, '(A)') '> '



!--------------------------------------------------------------------------------

! print the node types if node /= -9999

write(GP_print_unit,'(A)')  &
      'sgpi: i_GP_gen  i_GP_indiv   i_node     i_tree    GP_Indiv_Node_Type'


do  i_Node=1,n_Nodes
    do  i_Tree=1,n_Trees

        if( GP_Individual_Node_Type(i_Node,i_Tree) .ne. -9999         ) then

            write(GP_print_unit,'(5(1x,I10))') &
                  i_GP_generation, i_GP_indiv,i_node, i_tree, &
                  GP_Individual_Node_Type(i_Node,i_Tree)

            write(GP_summary_output_unit,'(5(1x,I10))') &
                  i_GP_generation, i_GP_indiv,i_node, i_tree, &
                  GP_Individual_Node_Type(i_Node,i_Tree)

        endif !   GP_Node_Type_Answer(i_Node,i_Tree) .ne. -9999 ...


    enddo ! i_tree
enddo  ! i_node


write(GP_summary_output_unit, '(A)') '> '

!---------------------------------------------------------------------------------



! print the node parameters (if there are any)

!write(GP_print_unit,'(/A)')     'sgpi: GP_population_node_parameters(i_GP_indiv,:,:) '
!write(GP_print_unit,'(5(1x,E15.7))')   GP_population_node_parameters(i_GP_indiv,:,:)

if( any( abs( GP_population_node_parameters(i_GP_indiv,:,:) ) &
                                                            > 1.0d-20 ) )then

    write(GP_print_unit,'(/A/)') &
       'sgpi: i_GP_gen i_GP_indiv     node        tree  GP_population_node_parameters'

    do  i_tree=1,n_trees
        do  i_node=1,n_nodes

            ! print only non-zero parameters

            if( abs( GP_population_node_parameters( &
                       i_GP_indiv,i_node,i_tree) ) > 1.0d-20   )then

                write(GP_print_unit,'(4(1x,I10), 1x, E20.10)') &
                      i_GP_generation, i_GP_indiv,i_node, i_tree, &
                      GP_population_node_parameters(i_GP_indiv,i_node,i_tree)

            endif ! abs( GP_population_node_parameters...
        enddo ! i_node
    enddo  ! i_tree

endif ! any( abs( GP_population_node_parameters(i_GP_indiv,:,:) )> 1.0d-20 )

! write all parameters to output file

do  i_tree=1,n_trees
    do  i_node=1,n_nodes

        write(GP_summary_output_unit,'(4(1x,I10), 1x, E20.10)') &
              i_GP_generation, i_GP_indiv,i_node, i_tree, &
              GP_population_node_parameters(i_GP_indiv,i_node,i_tree)

    enddo ! i_node
enddo  ! i_tree



write(GP_summary_output_unit, '(A)') '>>'


!---------------------------------------------------------------------------------


return

end subroutine summary_GP_indiv
