subroutine summary_GP_indiv( i_GP_indiv )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Adult_Population_Node_Type( i_GP_indiv,1:n_Nodes,1:n_Trees)

!  GP_Population_Initial_Conditions(i_GP_indiv,1:n_CODE_Equations)

!  GP_Population_Node_Parameters(i_GP_indiv,1:n_Nodes,1:n_Trees)


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


logical bad_node


integer :: i
integer :: i_code_eq 
integer :: nop
integer :: message_len

logical :: Lplot

integer :: i_GP_indiv

!real(kind=8), allocatable, dimension(:) :: answer
!real(kind=8), allocatable, dimension(:) :: output_array

real (kind=8) :: ssum, ssum2, totobs, dff

real (kind=8) :: xcount

!----------------------------------------------------------------------------------------

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Adult_Population_Node_Type(    i_GP_indiv,1:n_Nodes,1:n_Trees)
!  GP_Population_Initial_Conditions( i_GP_indiv,1:n_CODE_Equations)
!  GP_Population_Node_Parameters(    i_GP_indiv,1:n_Nodes,1:n_Trees)




! compute GP_Child_Individual_SSE(i_GP_indiv)

if( myid == 0 )then
    write(6,'(A)')  'sgpi: i_GP_generation, i_GP_indiv, i_node, i_tree, GP_Indiv_Node_Type '
endif ! myid == 0

dff = 0.0d0
xcount = 0.0d0
do  i_Node=1,n_Nodes
    do  i_Tree=1,n_Trees

        !if( myid == 0 )then
        !    write(6,'(A,4(1x,I6),1x,E15.7)') &
        !     'sgpi: i_node, i_tree, GP_Indiv_Node_Type, GP_Node_Type_Answer ',&
        !         i_node, i_tree, GP_Individual_Node_Type(i_Node,i_Tree), &
        !         GP_Node_Type_Answer(i_Node,i_Tree)
        !endif ! myid == 0


        if( GP_Individual_Node_Type(i_Node,i_Tree) .ne. -9999         ) then


            if( myid == 0 )then
                write(6,'(5x,5(1x,I6))') &
                      i_GP_generation, i_GP_indiv,i_node, i_tree, &
                      GP_Individual_Node_Type(i_Node,i_Tree)

            endif ! myid == 0


            !if( myid == 0 )then
            !    write(6,'(A,2(1x,I6))') &
            !          'sgpi: i_GP_generation, i_GP_indiv ', &
            !              i_GP_generation, i_GP_indiv
            !    write(6,'(A,4(1x,I6),1x,E15.7)') &
            !      'sgpi: i_node, i_tree, GP_Indiv_Node_Type'
            !          i_node, i_tree, &
            !          GP_Individual_Node_Type(i_Node,i_Tree), &
            !endif ! myid == 0

        endif !   GP_Node_Type_Answer(i_Node,i_Tree) .ne. -9999 ...


    enddo ! i_tree
enddo  ! i_node



!---------------------------------------------------------------------------------

if( myid == 0 )then


    ! print the node parameters (if there are any)

    !write(6,'(/A)')     'sgpi: GP_population_node_parameters(i_GP_indiv,:,:) '
    !write(6,'(5(1x,E15.7))')  GP_population_node_parameters(i_GP_indiv,:,:)

    if( any( abs( GP_population_node_parameters(i_GP_indiv,:,:) ) &
                                                                > 1.0d-20 ) )then

        write(6,'(/A/)') &
           'sgpi: i_GP_generation, i_GP_indiv, node  tree  GP_population_node_parameters'

        do  i_tree=1,n_trees
            do  i_node=1,n_nodes
                if( abs( GP_population_node_parameters( &
                           i_GP_indiv,i_node,i_tree) ) > 1.0d-20   )then

                    write(6,'(4(1x,I6), 1x, E20.10)') &
                          i_GP_generation, i_GP_indiv,i_node, i_tree, &
                          GP_population_node_parameters(i_GP_indiv,i_node,i_tree)
                endif
            enddo ! i_node
        enddo  ! i_tree

    endif ! any( abs( GP_population_node_parameters(i_GP_indiv,:,:) )> 1.0d-20 )

endif ! myid == 0



!---------------------------------------------------------------------------------

if( myid == 0 )then


    write(6,'(/A)')&
          'sgpi: i_GP_generation, i_GP_indiv, i_code_eq, &
                &GP_Pop_Init_Cond(i_GP_Indiv,i_code_eq) '
    do i_code_eq = 1, n_CODE_Equations

        write(6,'(3(1x,I6), 1x, E15.7)')&
                 i_GP_generation, i_GP_indiv, i_code_eq, &
                 GP_Population_Initial_Conditions( i_GP_indiv,i_code_eq )

    enddo  ! i_code_eq 

endif ! myid == 0



return

end subroutine summary_GP_indiv
