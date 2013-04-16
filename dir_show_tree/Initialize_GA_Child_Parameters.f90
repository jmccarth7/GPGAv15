subroutine Initialize_GA_Child_Parameters(Child_Parameters)


!use mpi
!use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none


real(kind=8) :: Child_Parameters( n_GA_Individuals,n_trees, n_nodes )
real(kind=4) :: cff
real(kind=8) :: dff

integer(kind=4) :: i

!----------------------------------------------------------------------------

!Run_GA_lmdif=.true.


write(6,'(/A,2(1x, I6)/)')  'Init: n_trees, n_nodes ',  n_trees, n_nodes

!write(6,'(/A/)') &
!      'Init:  i_GA_individual  child parameters(i_GA_individual,:)  '


do  i_GA_Individual=1,n_GA_Individuals


    do  i_tree=1,n_trees
        do  i_node=1,n_nodes


        call random_real(cff) ! random real number generator
        dff = cff
        Child_Parameters(i_GA_Individual,i_tree, i_node) = dff


            write(6,'(A,1x, I6,1x,12(1x,E15.7))') &
                  'Init: ', i_GA_individual, &
                        child_parameters(i_GA_individual,i_tree,i_node)

        enddo ! i_node
    enddo ! i_tree


enddo ! i_ga_individual



return


end subroutine Initialize_GA_Child_Parameters
!234567890123456789012345678901234567890123456789012345678901234567890
