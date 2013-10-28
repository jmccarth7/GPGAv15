subroutine summary_GP_indiv( i_GP_generation, i_GP_indiv )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Adult_Population_Node_Type(    1:n_Nodes,1:n_Trees,i_GP_indiv  )

!  GP_Population_Initial_Conditions( 1:n_CODE_Equations, i_GP_indiv )

!  GP_Population_Node_Parameters(    1:n_Nodes,1:n_Trees,i_GP_indiv  )


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


integer(kind=4),intent(in)  :: i_GP_Generation
integer(kind=4),intent(in)  :: i_GP_indiv

integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

logical :: Lprint

!----------------------------------------------------------------------------------------

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Adult_Population_Node_Type(    1:n_Nodes,1:n_Trees,i_GP_indiv)
!  GP_Population_Initial_Conditions( 1:n_CODE_Equations, i_GP_indiv)
!  GP_Population_Node_Parameters(    1:n_Nodes,1:n_Trees,i_GP_indiv)

!---------------------------------------------------
! assume this subroutine is called only by cpu 0
!---------------------------------------------------


!--------------------------------------------------------------------------------

! set Lprint so printing is done only under the conditions in the if-test

Lprint = .TRUE.

!if( i_GP_generation == 1                                  .or. &
!    mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
!    i_GP_generation == n_GP_generations                          )then
!    Lprint = .TRUE.
!endif ! i_GP_generation == 1 .or. ...

!--------------------------------------------------------------------------------


! write the summary file header for each individual
! which has n_GP_parameters >= n_code_equations

!if( Lprint )then
!    write(GP_print_unit, '(/A/6x,6(1x,I10))') &
!      'sgpi: n_code_equations  n_trees  n_nodes  n_levels  &
!           & i_GP_individual  i_GP_generation', &
!             n_code_equations, n_trees, n_nodes, n_levels, &
!             i_GP_indiv, i_GP_generation
!endif ! Lprint
!
!write(GP_summary_output_unit, '(6(1x,I6))') &
!             n_code_equations, n_trees, n_nodes, n_levels, &
!             i_GP_indiv, i_GP_generation

if( Lprint )then
    write(GP_print_unit, '(/A/7(1x,I10))') &
      'sgpi: i_GP_gen i_GP_indiv   n_code_eq  &
             &n_trees    n_nodes  n_levels    n_parms', &
             i_GP_generation, i_GP_indiv, &
             n_code_equations, n_trees, n_nodes, n_levels, &
             GP_Individual_N_GP_param(i_GP_indiv)
endif ! Lprint

write(GP_summary_output_unit, '(2x,6(1x,I6))') &
             i_GP_generation, i_GP_indiv, &
             n_code_equations, n_trees, n_nodes, n_levels


!--------------------------------------------------------------------------------

! initial conditions


if( Lprint )then
    write(GP_print_unit,'(/A)')&
      'sgpi: i_GP_gen  i_GP_indiv  i_code_eq  &
            &GP_Pop_Init_Cond(i_code_eq, i_GP_Indiv) '
endif ! Lprint

do  i_code_eq = 1, n_CODE_Equations

    if( Lprint )then
        write(GP_print_unit,'(3(1x,I10), 7x, E24.16)')&
              i_GP_generation, i_GP_indiv, i_code_eq, &
              GP_Population_Initial_Conditions( i_code_eq, i_GP_indiv )
    endif ! Lprint

    write(GP_summary_output_unit, '(2x,2(1x,I6),1x,I3, 1x, E24.16,2x,A)')&
          i_GP_generation, i_GP_indiv, i_code_eq, &
          GP_Population_Initial_Conditions( i_code_eq, i_GP_indiv ), &
          'gen_indiv_eq'

enddo  ! i_code_eq


write(GP_summary_output_unit, '(A,2(1x,I6))') '> ', i_GP_generation, i_GP_indiv



!!--------------------------------------------------------------------------------
!
!
!! print the node types if node /= -9999
!
!
!!if( Lprint )then
!!    write(GP_print_unit,'(A)')  &
!!          'sgpi: i_GP_gen  i_GP_indiv   i_tree     i_node    GP_Indiv_Node_Type'
!!endif ! Lprint
!
!
!do  i_Node=1,n_Nodes
!    do  i_Tree=1,n_Trees
!
!        if( GP_Individual_Node_Type(i_Node,i_Tree) .ne. -9999         ) then
!
!            !if( Lprint )then
!            !    write(GP_print_unit,'(2(1x,I10),3(1x,I10))') &
!            !          i_GP_generation, i_GP_indiv,i_tree, i_node, &
!            !          GP_Individual_Node_Type(i_Node,i_Tree)
!            !endif ! Lprint
!
!            write(GP_summary_output_unit, '(2x,2(1x,I6),3(1x,I3))') &
!                  i_GP_generation, i_GP_indiv,i_tree, i_node, &
!                  GP_Individual_Node_Type(i_Node,i_Tree)
!
!        endif !   GP_Node_Type_Answer(i_Node,i_Tree) .ne. -9999 ...
!
!
!    enddo ! i_tree
!enddo  ! i_node
!
!
!write(GP_summary_output_unit, '(A,2(1x,I6))') '> ', i_GP_generation, i_GP_indiv
!
!
!!---------------------------------------------------------------------------------
!
!
!! print the node parameters (if there are any)
!
!!write(GP_print_unit,'(/A)')     'sgpi: GP_population_node_parameters(:,:,i_GP_indiv) '
!!write(GP_print_unit,'(5(1x,E15.7))')   GP_population_node_parameters(:,:,i_GP_indiv)
!
!
!if( Lprint )then
!    write(GP_print_unit,'(/A/)') &
!       'sgpi: i_GP_gen i_GP_indiv     tree        node   GP_population_node_parameters'
!
!    do  i_tree=1,n_trees
!        do  i_node=1,n_nodes
!
!            ! print only non-zero parameters
!
!            if( GP_Individual_Node_Type(i_Node,i_Tree) == 0               ) then
!                write(GP_print_unit,'(2(1x,I10), 2(1x,I10),2x, E24.16)') &
!                      i_GP_generation, i_GP_indiv,i_tree, i_node, &
!                      GP_population_node_parameters(i_node,i_tree,i_GP_indiv)
!
!            endif  ! GP_Individual_Node_Type(i_Node,i_Tree) == 0
!
!        enddo ! i_node
!    enddo  ! i_tree
!
!endif ! Lprint
!
!
!!---------------------------------------------------------------------------------

!--------------------------------------------------------------------------------


! print the node types if node /= -9999


if( Lprint )then
    write(GP_print_unit,'(A)')  &
     'sgpi: i_GP_gen i_GP_indiv i_tree   i_node  Node_Type   &
     &GP_population_node_parameters'
endif ! Lprint


do  i_Node=1,n_Nodes
    do  i_Tree=1,n_Trees

        if( GP_Individual_Node_Type(i_Node,i_Tree) .ne. -9999         ) then


            if( GP_Individual_Node_Type(i_Node,i_Tree) == 0               ) then

                write(GP_print_unit,'(2(1x,I8), 3(1x,I8),8x, E24.16)') &
                      i_GP_generation, i_GP_indiv,i_tree, i_node, &
                      GP_Individual_Node_Type(i_Node,i_Tree), &
                      GP_population_node_parameters(i_node,i_tree,i_GP_indiv)

            else

                write(GP_print_unit,'(2(1x,I8),3(1x,I8))') &
                      i_GP_generation, i_GP_indiv,i_tree, i_node, &
                      GP_Individual_Node_Type(i_Node,i_Tree)

            endif  ! GP_Individual_Node_Type(i_Node,i_Tree) == 0



            write(GP_summary_output_unit, '(2x,2(1x,I6),3(1x,I3))') &
                  i_GP_generation, i_GP_indiv,i_tree, i_node, &
                  GP_Individual_Node_Type(i_Node,i_Tree)

        endif !   GP_Node_Type_Answer(i_Node,i_Tree) .ne. -9999 ...


    enddo ! i_tree
enddo  ! i_node


write(GP_summary_output_unit, '(A,2(1x,I6))') '> ', i_GP_generation, i_GP_indiv


!---------------------------------------------------------------------------------


! print the node parameters (if there are any)

!write(GP_print_unit,'(/A)')     'sgpi: GP_population_node_parameters(:,:,i_GP_indiv) '
!write(GP_print_unit,'(5(1x,E15.7))')   GP_population_node_parameters(:,:,i_GP_indiv)


!if( Lprint )then
!    write(GP_print_unit,'(/A/)') &
!       'sgpi: i_GP_gen i_GP_indiv     tree        node   GP_population_node_parameters'
!
!    do  i_tree=1,n_trees
!        do  i_node=1,n_nodes
!
!            ! print only non-zero parameters
!
!            if( GP_Individual_Node_Type(i_Node,i_Tree) == 0               ) then
!                write(GP_print_unit,'(2(1x,I10), 2(1x,I10),2x, E24.16)') &
!                      i_GP_generation, i_GP_indiv,i_tree, i_node, &
!                      GP_population_node_parameters(i_node,i_tree,i_GP_indiv)
!
!            endif  ! GP_Individual_Node_Type(i_Node,i_Tree) == 0
!
!        enddo ! i_node
!    enddo  ! i_tree
!
!endif ! Lprint


!---------------------------------------------------------------------------------

if( Lprint )then
    call print_trees( i_GP_generation, i_GP_indiv, i_GP_indiv, &
                      GP_Adult_Population_Node_Type, ' ' )
endif ! Lprint

!call print_trees( i_GP_generation, i_GP_indiv, i_GP_indiv, &
!                        GP_Individual_Node_Type, '    ' )

!---------------------------------------------------------------------------------


! write all non-zero parameters to output file


do  i_tree=1,n_trees
    do  i_node=1,n_nodes

        if( GP_Individual_Node_Type(i_Node,i_Tree) == 0               ) then

            write(GP_summary_output_unit,'(2x,2(1x,I6),2(1x,I3), 1x,E24.16)') &
                  i_GP_generation, i_GP_indiv,i_tree, i_node, &
                  GP_population_node_parameters( i_node,i_tree, i_GP_indiv)

        endif ! GP_Individual_Node_Type(i_Node,i_Tree) == 0               ) then

    enddo ! i_node
enddo  ! i_tree



write(GP_summary_output_unit, '(A,2(1x,I6))') '>>', i_GP_generation, i_GP_indiv


!! write for each indiv.  first write in 0*.f90
!write(GP_summary_output_unit, '(4(1x,I10))') n_code_equations, n_trees, n_nodes, n_levels


!---------------------------------------------------------------------------------


return

end subroutine summary_GP_indiv
