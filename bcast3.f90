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




integer :: i
integer :: i_diversity
integer :: message_len

!logical :: Lplot


real(kind=8), allocatable, dimension(:) :: answer
real(kind=8), allocatable, dimension(:) :: output_array


character(200) :: tree_descrip

!----------------------------------------------------------------------------------------



!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          '0: broadcast GP_Adult_Individual_SSE = ',&
!                        GP_Adult_Individual_SSE
!endif ! myid == 0

call MPI_BCAST( GP_Adult_Individual_SSE, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          '0: broadcast GP_Population_Ranked_Fitness ',&
!                        GP_Population_Ranked_Fitness
!endif ! myid == 0

call MPI_BCAST( GP_Population_Ranked_Fitness, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          '0: AFT broadcast GP_Population_Ranked_Fitness ',&
!                            GP_Population_Ranked_Fitness
!endif ! myid == 0


!if( myid == 0 )then
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
!          '0: broadcast GP_Integrated_Population_Ranked_Fitness ',&
!                        GP_Integrated_Population_Ranked_Fitness
!endif ! myid == 0

call MPI_BCAST( GP_Integrated_Population_Ranked_Fitness, n_GP_individuals,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )




return

end subroutine bcast3
