subroutine  print_gp_node_type_parm( )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
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


integer(kind=4) :: i_GP_individual
integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

!----------------------------------------------------------------------------------------


do  i_GP_individual=1,n_GP_individuals

    write(GP_print_unit,'(/A,1x,I6/)') &
          'pgntp: i_GP_individual = ', i_GP_individual
    write(GP_print_unit,'(A)') &
          'pgntp: i_node    i_tree  GP_Adult_Population_Node_Type(i_Node,i_Tree)'

    do  i_Tree=1,n_Trees
        do  i_Node=1,n_Nodes
            if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_individual) /= -9999 )then
                write(GP_print_unit,'(3(1x,I10))') &
                      i_node, i_tree, &
                      GP_Adult_Population_Node_Type( i_Node,i_Tree,i_GP_individual)
            endif ! GP_Adult_Pop_Node_Type(i_Node,i_Tree,i_GP_indiv) /= -9999
        enddo ! i_node
    enddo ! i_tree


    write(GP_print_unit,'(//A)') &
          'pgntp:  i_node  i_tree  GP_Population_Node_Parameters(i_Node,i_Tree)'
    do  i_Tree=1,n_Trees
        do  i_Node=1,n_Nodes

            if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_individual) == 0 )then
                write(GP_print_unit,'(2(1x,I10),2x,E15.7)') &
                      i_node, i_tree, &
                      GP_Population_Node_Parameters(i_Node,i_Tree,i_GP_individual)
            endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_individual)

        enddo ! i_node
    enddo ! i_tree
    write(GP_print_unit,'(/A)') ' '


enddo !  i_GP_individual


return

end subroutine  print_gp_node_type_parm
