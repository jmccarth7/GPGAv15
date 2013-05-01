subroutine print_values1( )

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


!real(kind=8), allocatable, dimension(:) :: answer
!real(kind=8), allocatable, dimension(:) :: output_array


character(200) :: tree_descrip

integer(kind=4) :: i_GP_individual                                                                           
integer(kind=4) :: i_GP_Best_Parent                                                                          
integer(kind=4) :: i_GP_Generation                                                                           
integer(kind=4) :: j_GP_Individual                                                                           
integer(kind=4) :: i_Tree                                                                                    
integer(kind=4) :: i_Node            
!----------------------------------------------------------------------------------------


write(GP_print_unit,'(A,1x,I10)')    '0: n_GA_individuals           ', &
                                         n_GA_individuals
write(GP_print_unit,'(A,1x,I10)')    '0: n_time_steps               ', &
                                         n_time_steps
write(GP_print_unit,'(A,1x,I10)')    '0: n_GA_Generations           ', &
                                         n_GA_Generations
write(GP_print_unit,'(A,1x, E15.7)') '0: GA_Crossover_Probability   ', &
                                         GA_Crossover_Probability
write(GP_print_unit,'(A,1x, E15.7)') '0: GA_Mutation_Probability    ', &
                                         GA_Mutation_Probability
write(GP_print_unit,'(A,1x, E15.7)') '0: GA_save_elites_Probability ', &
                                         GA_save_elites_Probability

write(GP_print_unit,'(/A,1x, E15.7)')'0: GP_Tree_Probability        ', &
                                         GP_Tree_Probability

write(GP_print_unit,'(A,1x, E15.7)') '0: GP_Elitist_Probability     ', &
                                         GP_Elitist_Probability

write(GP_print_unit,'(A,1x, E15.7)') '0: GP_Crossover_Probability   ', &
                                         GP_Crossover_Probability

write(GP_print_unit,'(A,1x, E15.7)') '0: GP_Asexual_Reproduction_Probability ', &
                                         GP_Asexual_Reproduction_Probability

write(GP_print_unit,'(A,1x, E15.7)') '0: GP_Mutation_Probability    ', &
                                         GP_Mutation_Probability
write(GP_print_unit,'(/A,1x,I10)')   '0: n_gp_individuals           ', &
                                         n_gp_individuals
write(GP_print_unit,'(A,1x,I10)')    '0: n_gp_generations           ', &
                                         n_gp_generations



return

end subroutine print_values1
