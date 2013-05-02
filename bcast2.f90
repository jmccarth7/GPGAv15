subroutine bcast2( )

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



integer :: message_len


!----------------------------------------------------------------------------------------




        ! broadcast



        message_len = n_GP_Individuals * n_Nodes * n_Trees
        call MPI_BCAST( GP_Child_Population_Node_Type, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

        !if( myid == 0 )then
        !    write(GP_print_unit,'(A,1x,I6)') &
        !      '0: aft broadcast  GP_Child_Pop_Node_Type  Generation = ',&
        !                                            i_GP_Generation
        !endif ! myid == 0

        ! GP_Adult_Population_Node_Type

        message_len = n_GP_Individuals * n_Nodes * n_Trees
        call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

        ! Parent_Tree_Swap_Node_Type

        message_len = 2 * n_Nodes
        call MPI_BCAST( Parent_Tree_Swap_Node_Type, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

        ! GP_Child_Individual_SSE

        call MPI_BCAST( GP_Child_Individual_SSE, n_GP_individuals,    &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

        ! GP_Integrated_Population_Ranked_Fitness

        call MPI_BCAST( GP_Integrated_Population_Ranked_Fitness, n_GP_individuals,    &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

        ! GP_Population_Ranked_Fitness

        call MPI_BCAST( GP_Population_Ranked_Fitness, n_GP_individuals,    &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

        ! Run_GP_Calculate_Fitness array

        call MPI_BCAST( Run_GP_Calculate_Fitness , n_GP_Individuals,    &
                        MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )




return

end subroutine bcast2
