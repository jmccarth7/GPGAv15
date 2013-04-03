module GA_variables_module

use GA_parameters_module
implicit none

integer (kind=4) :: n_GA_Elitists
integer (kind=4) :: n_GA_Asexual_Repoductions
integer (kind=4) :: n_GA_Crossovers
integer (kind=4) :: n_GA_Mutations

integer (kind=4) :: n_GA_crossovers
integer (kind=4) :: n_GA_mutations
integer (kind=4) :: n_GA_save_elites

!integer(kind=4), dimension( n_GA_individuals ) :: individual_elites
integer(kind=4), allocatable, dimension(:) :: individual_elites
real (kind=4) :: GA_Individual_Lowest_SSE
real (kind=4) :: GA_Adult_Individual_SSE(n_GA_Individuals)       ! must be kept for re-evaluations of next generations
real (kind=4) :: GA_Child_Individual_SSE(n_GA_Individuals)       ! must be kept for re-evaluations of next generations
real (kind=4) :: GA_Integrated_SSE(n_GA_Individuals)             ! needed to support sexual and "tournament-style" reproduction
real (kind=4) :: GA_Individual_Ranked_Fitness(n_GA_Individuals)  ! must be kept for re-evaluations of next generations
real (kind=4) :: GA_Integrated_Ranked_Fitness(n_GA_Individuals)  ! must be kept for re-evaluations of next generations



!real (kind=8) :: individual_SSE(n_GA_individuals)   
!     must be kept for re-evaluations of next generations
real(kind=8),allocatable, dimension(:) :: individual_SSE  

!real (kind=8) :: integrated_SSE(n_GA_individuals)   
!     needed to support sexual and "tournament-style" reproduction
real(kind=8),allocatable, dimension(:) :: integrated_SSE 

!real (kind=8) :: individual_ranked_fitness(n_GA_individuals)  
!      must be kept for re-evaluations of next generations
real(kind=8),allocatable, dimension(:) :: individual_ranked_fitness 

!real (kind=8) :: integrated_ranked_fitness(n_GA_individuals)  
!      must be kept for re-evaluations of next generations
real(kind=8),allocatable, dimension(:) :: integrated_ranked_fitness


!real (kind=8) :: fitness_expectation_value(n_GA_individuals)
real(kind=8),allocatable, dimension(:) :: fitness_expectation_value


real(kind=8) :: sum_individual_SSE
real(kind=8) :: sum_individual_fit
real(kind=8) :: SSE0
real(kind=8) :: max_sse

!logical :: Run_GA_lmdif(n_GA_individuals)
logical,allocatable, dimension(:) :: Run_GA_lmdif
logical, dimension(n_GA_Individuals) :: Run_GA_lmdif

end module GA_variables_module
