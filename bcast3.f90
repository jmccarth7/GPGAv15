subroutine bcast3( )

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



integer :: buffer_length


!-------------------------------------------------------------------------------

! GP_Adult_Individual_SSE


call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          '0: broadcast GP_Adult_Individual_SSE = ',&
!                        GP_Adult_Individual_SSE
!endif ! myid == 0

call MPI_BCAST( GP_Adult_Individual_SSE, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!------------------------------------------------------------------------------

! GP_population_node_parameters


!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          '0: broadcast GP_population_node_parameters ',&
!                        GP_population_node_parameters
!endif ! myid == 0

!call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

buffer_length = n_nodes * n_trees * n_GP_individuals

call MPI_BCAST( GP_population_node_parameters,  buffer_length,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          '0: AFT broadcast GP_population_node_parameters ',&
!                            GP_population_node_parameters
!endif ! myid == 0

!------------------------------------------------------------------------------

! GP_Population_Ranked_Fitness


!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          '0: broadcast GP_Population_Ranked_Fitness ',&
!                        GP_Population_Ranked_Fitness
!endif ! myid == 0


!call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

call MPI_BCAST( GP_Population_Ranked_Fitness, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          '0: AFT broadcast GP_Population_Ranked_Fitness ',&
!                            GP_Population_Ranked_Fitness
!endif ! myid == 0

!------------------------------------------------------------------------------

! GP_Integrated_Population_Ranked_Fitness


!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          '0: broadcast GP_Integrated_Population_Ranked_Fitness ',&
!                        GP_Integrated_Population_Ranked_Fitness
!endif ! myid == 0

!call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

call MPI_BCAST( GP_Integrated_Population_Ranked_Fitness, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!------------------------------------------------------------------------------

call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?


return

end subroutine bcast3
