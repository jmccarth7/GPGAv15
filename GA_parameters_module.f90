module GA_parameters_module

! set the GA related parameters

implicit none

! set the integer parameters

integer(kind=4),parameter :: GA_print_unit = 6  ! 60


logical :: Lprint_lmdif
! set the real parameters

!  NOTE: in the next 2 parameters:
!   GA_Crossover_Probability + GA_Mutation_Probability must be less than <= 1.


! Note: The next 4 parameters must add up to 1.0

! GA_Crossover_Probability  ! probability of sexual crossing of parameter strings in GA_lmdif
! GA_Mutation_Probability   ! probability of mutation in parameter string of GA_lmdif
! GA_Elitist_Probability    ! Keeps the top n_GP_Elitists of the Best Fit Individuals
                            !from Generation to Generation
! GA_Asexual_Reproduction_Probability   ! probability of asexual reproduction





real(kind=8) :: GA_Crossover_Probability
real(kind=8) :: GA_Mutation_Probability
real(kind=8) :: GA_save_elites_Probability
real(kind=8) :: GA_Elitist_Probability
real(kind=8) :: GA_Asexual_Reproduction_Probability


! GP Probability of a Tree being assigned
! Estimated from previous work by Joel Cohen
!real (kind=4), parameter :: GP_Tree_Probability=0.5

! probability of setting a terminal node to a parameter
real (kind=4), parameter :: GP_Set_Terminal_to_Parameter_Probability = 0.6


! determines how the GA_Tournament subroutine works

! = 0  - swap unmodified segments of parents
! = 1  - swap segments of parents and randomly reset node at segment boundaries
! = 2  - swap segments of parents and reset node at segment boundaries using JM
!        formula involving the mean and std. dev

integer(kind=4) :: ga_tournament_style

real (kind=4), parameter :: PI = 3.141592653589793D0

! GA routine-specific variables that can be modified

integer(kind=4) :: n_GA_Generations
integer(kind=4) :: i_GA_Generation
integer(kind=4) :: n_GA_Individuals
integer(kind=4) :: i_GA_individual
integer(kind=4) :: j_GA_individual

real(kind=8) ::  min_sse

end module GA_parameters_module
