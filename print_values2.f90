subroutine print_values2( answer )

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
real(kind=8), dimension(n_parameters) :: answer


!----------------------------------------------------------------------------------------

write(GP_print_unit,'(/A,1x,I6)') 'Total Parameters for this run = ',n_parameters
write(GP_print_unit,'(A,1x,I6)')  '0: n_GA_individuals ', n_GA_individuals

!write(GP_print_unit,'(A,1x,I6)') 'number of initial conditions   = ',n_initial_conditions
write(GP_print_unit,'(A,1x,I6)')   '0: n_time_steps     ', n_time_steps
write(GP_print_unit,'(A,1x,I6)')   '0: n_GA_Generations ', n_GA_Generations
write(GP_print_unit,'(A,1x,I6)')   '0: n_levels         ', n_levels
write(GP_print_unit,'(A,1x,I6)')   '0: n_functions      ', n_functions
write(GP_print_unit,'(A,1x,I6)')   '0: n_NODE_functions ', n_NODE_functions
write(GP_print_unit,'(A,1x,I6)')   '0: n_CODE_equations ', n_CODE_equations
write(GP_print_unit,'(A,1x,I6)')   '0: n_parameters     ', n_parameters
write(GP_print_unit,'(A,1x,I6)')   '0: n_trees          ', n_trees
write(GP_print_unit,'(A,1x,I6)')   '0: n_nodes          ', n_nodes
write(GP_print_unit,'(A,1x,E15.7)')'0: dt (days)        ', dt
write(GP_print_unit,'(A,1x,E15.7)')'0: dt (min)         ', dt * 1440.d0
!write(GP_print_unit,*) individual_fitness
!write(GP_print_unit,*) data_array
write(GP_print_unit,'(A,1x,E15.7)') '0: GA_Crossover_Probability  ', &
                                        GA_Crossover_Probability
write(GP_print_unit,'(A,1x,E15.7)') '0: GA_Mutation_Probability   ', &
                                        GA_Mutation_Probability
write(GP_print_unit,'(A,1x,E15.7)') '0: GA_save_elites_Probability', &
                                        GA_save_elites_Probability

write(GP_print_unit,'(A)')  '0: code calls lmdif only for best individual on last generation'
write(GP_print_unit,'(A)')  '0: 2-range random_real initialization of child parameters '
write(GP_print_unit,'(A//)')'0: replacement of bad individuals used'


! calculate the generation interval for printing the list of children

GA_child_print_interval = n_GA_generations /  number_GA_child_prints
GP_child_print_interval = n_GP_generations /  number_GP_child_prints

write(GP_print_unit,'(A,1x,I6)')   '0: GP_child_print_interval ', &
                                       GP_child_print_interval
write(GP_print_unit,'(A,1x,I6)')   '0: GA_child_print_interval ', &
                                       GA_child_print_interval
!-----------------------------------------------------------------------------

write(GP_print_unit,'(/A)') ' '
do  i = 1, n_parameters
    write(GP_print_unit,'(A,1x,I6,2x,E24.16)') &
          '0: i, answer(i) ', i, answer(i)
enddo ! i
write(GP_print_unit,'(/A)') ' '


return

end subroutine print_values2