module GA_parameters_module

! set the GA related parameters

implicit none

! set the integer parameters

integer(kind=4),parameter :: GA_print_unit = 60

! set the real parameters

!  NOTE: in the next 2 parameters:  GA_Crossover_Probability+GA_Mutation_Probability must be less than <= 1.

! Note: The next 4 parameters must add up to 1.0
!   GA_Crossover_Probability  ! probability of sexual crossing of parameter strings in GA_lmdif
!   GA_Mutation_Probability   ! probability of mutation in parameter string of GA_lmdif
!   GA_Elitist_Probability  !=0.1   ! Keeps the top n_GP_Elitists of the Best Fit Individuals from Generation to Generation
!   GA_Asexual_Reproduction_Probability=0.4   ! probability of asexual reproduction





real(kind=8) :: GA_Crossover_Probability  ! probability of sexual crossing of parameter strings in GA_lmdif
!real(kind=4), parameter :: GA_Crossover_Probability = 0.4  ! probability of sexual crossing of binary string 

real(kind=8) :: GA_Mutation_Probability   ! probability of mutation in parameter string of GA_lmdif
!real (kind=4), parameter :: GA_Mutation_Probability = 0.1   ! probability of mutation in binary string

real(kind=8) :: GA_save_elites_Probability ! probability of saving an individual as an elite individual

real(kind=8)  :: GA_Elitist_Probability  !=0.1   ! Keeps the top n_GP_Elitists of the Best Fit Individuals from Generation to Generation

real(kind=8)  :: GA_Asexual_Reproduction_Probability !=0.4   ! probability of asexual reproduction


! GP Probability of a Tree being assigned
!real (kind=4), parameter :: GP_Tree_Probability=0.5 ! Estimated from previous work by Joel Cohen

! probability of setting a terminal node to a parameter
real (kind=4), parameter :: GP_Set_Terminal_to_Parameter_Probability = 0.6  



! GA routine-specific variables that can be modified

integer(kind=4) :: n_GA_Generations !=  1000 !  5000 ! 1000

integer(kind=4) :: i_GA_Generation

integer(kind=4) :: n_GA_Individuals !=  500 ! 3  ! 100

integer(kind=4) :: i_GA_individual, j_GA_individual

real(kind=8) ::  min_sse

end module GA_parameters_module
